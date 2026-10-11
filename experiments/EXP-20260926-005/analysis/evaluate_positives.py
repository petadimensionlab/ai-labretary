"""EXP-20260926-005: 新しい陽性対照セット（delete 型 / overexpress 型）を plan.md §6 の事前登録基準で評価する。

EXP-20260926-002 の `reanalysis_cosine.py`（コサイン候補指標・null 分位）と `vector_phase0.py`（ベクトル候補）
をそのまま import して使い、この実験のために次だけを足す（EXP-002 側のコードは変更しない）:
  - 陽性対照が config から来る（EXP-002 の POS 固定を差し替え）
  - C1'（陽性 N 中 K 以上が null 上位 10% に 3 seed 中 2 以上）、C2'（C1' 合格遺伝子 vs HK 9 組中 7 以上）
  - C6（C1' 合格遺伝子に ddet_pct ≤ 0.80 のものがあるか）、C7（run O の方向対照が null 下位 50% か）
  - 分類ヘッドの ΔP(PT_healthy)（dP_raw と、a2 と同じ長さビン補正の dP_a2）
盲検: 主要対象・参照（rc.BLIND）は読み込み直後に除外（assert）。
実行:
  .venv/bin/python experiments/EXP-20260926-005/analysis/evaluate_positives.py --config experiments/EXP-20260926-005/config-delete.yaml
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "experiments/EXP-20260926-002/analysis"))
import reanalysis_cosine as rc  # noqa: E402
import vector_phase0 as vp  # noqa: E402

E = ROOT / "experiments/EXP-20260926-005"
TOK = E / "work/tokenized.dataset"   # 自分の複製だけを読む（guardrails: 他実験の dataset を直接指さない）


def token_positions(cells: pd.DataFrame) -> pd.Series:
    """摂動した遺伝子の input_ids 内の位置（rc.token_positions と同じだが、この実験の複製を読む）。"""
    import pickle
    from datasets import load_from_disk
    from geneformer import TOKEN_DICTIONARY_FILE
    tokd = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
    ds = load_from_disk(str(TOK))
    uids = sorted(cells.cell_uid.unique())
    ids = dict(zip(uids, ds.select(uids)["input_ids"]))
    return pd.Series([ids[u].index(tokd[g]) for u, g in zip(cells.cell_uid, cells.gene)], index=cells.index)


def bin_correct(df: pd.DataFrame, base: str, out: str, nbins: int = 5) -> np.ndarray:
    """a2 と同じ補正（長さ nbins 分位ビンごとの null 平均を引く。null 遺伝子は leave-one-gene-out）。"""
    v = np.full(len(df), np.nan)
    null_genes = sorted(df.loc[df.role == "negative", "gene"].unique())
    groups = [(df.gene == g, (df.role == "negative") & (df.gene != g)) for g in null_genes]
    groups.append((df.role != "negative", df.role == "negative"))
    q = np.linspace(0, 1, nbins + 1)[1:-1]
    for tgt, fit in groups:
        F, T = df[fit], df[tgt]
        edges = np.quantile(F.length, q)
        bf, bt = np.digitize(F.length, edges), np.digitize(T.length, edges)
        mu = F.groupby(bf)[base].mean()
        v[tgt.to_numpy()] = T[base].to_numpy() - mu.reindex(bt).to_numpy()
    df[out] = v
    return v


def load_run(res: Path, run: str, want_logits: bool):
    cells = pd.read_csv(res / "cells.csv", dtype={"donor_id": str})
    cells = cells[(cells.run == run) & ~cells.gene.isin(rc.BLIND)].reset_index(drop=True)
    negs = set(json.loads((res / run / "negatives.json").read_text()))
    cells.loc[cells.gene.isin(negs) & (cells.role == "screen"), "role"] = "negative"
    assert not set(cells.gene) & rc.BLIND, "盲検の遺伝子が結果に含まれています"
    rd = res / run
    E0, E1, dP, order = [], [], [], []
    for gene, s in cells.groupby("gene", sort=False):
        o = np.load(rd / "genes" / gene / "emb" / "orig.npy").astype(np.float64)
        p = np.load(rd / "genes" / gene / "emb" / "pert.npy").astype(np.float64)
        per = pd.read_csv(rd / "genes" / gene / "cells.csv", dtype={"donor_id": str})
        assert len(o) == len(per) == len(s) and (per.cell_uid.to_numpy() == s.cell_uid.to_numpy()).all()
        E0.append(o); E1.append(p); order.append(s.index.to_numpy())
        if want_logits:
            lo = np.load(rd / "genes" / gene / "emb" / "logits_orig.npy").astype(np.float64)
            lp = np.load(rd / "genes" / gene / "emb" / "logits_pert.npy").astype(np.float64)
            sm = lambda x: np.exp(x[:, 1] - np.logaddexp(x[:, 0], x[:, 1]))  # noqa: E731  P(PT_healthy)
            dP.append(sm(lp) - sm(lo))
    cells = cells.loc[np.concatenate(order)].reset_index(drop=True)
    st = np.load(rd / "state_emb_vectors.npz")
    cen = {k: st[k].astype(np.float64).ravel() for k in st.files}
    return cells, np.vstack(E0), np.vstack(E1), (np.concatenate(dP) if want_logits else None), cen


def quantiles(g: pd.DataFrame, metric: str) -> dict:
    """(run, gene) → その run の null 20 遺伝子に対する分位。"""
    out = {}
    for run, s in g.groupby("run"):
        null = s[s.role == "negative"][metric].to_numpy()
        for _, r in s.iterrows():
            out[(run, r.gene)] = rc.null_pct(r[metric], null)
    return out


def assess(g: pd.DataFrame, metric: str, cfg: dict, names: dict, ann: dict) -> dict:
    crit = cfg["analysis"]["criteria"]
    pos = [x for x in cfg["controls"]["positive"]]
    hk = [x for x in cfg["controls"]["negative_fixed"]]
    roles = {grp["role"]: grp["genes"] for grp in cfg["controls"].get("extra_fixed", []) or []}
    panel, extra = roles.get("negctrl_panel", []), roles.get("prev_positive", []) + roles.get("reverse_direction", [])
    runs = sorted(g.run.unique())
    need = crit["C1_min_seeds"] if len(runs) >= 3 else len(runs)
    pct = quantiles(g, metric)
    v = {(r.run, r.gene): r[metric] for _, r in g.iterrows()}

    def n_hi(gene):
        return sum(pct.get((r, gene), 0) >= crit["C1_positive_null_quantile"] for r in runs)

    c1_pass = [x for x in pos if n_hi(x) >= need]
    c2 = {}
    for x in c1_pass:
        c2[names[x]] = sum(v[(r, x)] > v[(r, h)] for r in runs for h in hk)
    c3 = {names[h]: sum(crit["C3_hk_quantile_range"][0] <= pct.get((r, h), -1) <= crit["C3_hk_quantile_range"][1]
                        for r in runs) for h in hk}
    lo, hi = crit["C3_panel_quantile_range"]
    cp = {names[x]: sum(lo <= pct.get((r, x), -1) <= hi for r in runs) for x in panel}
    n_panel_ok = sum(n >= crit["C3_panel_min_seeds"] for n in cp.values())
    c6_genes = [names[x] for x in c1_pass if ann.get(names[x], {}).get("ddet_pct", 1.0) <= crit["C6_marker_pct_max"]]
    c7 = {names[x]: sum(pct.get((r, x), 1.0) <= crit["C7_reverse_quantile_max"] for r in runs)
          for x in roles.get("reverse_direction", [])}
    return {
        "metric": metric,
        "n_seeds": len(runs),
        "pos_seeds_ge_q90": {names[x]: n_hi(x) for x in pos},
        "pos_pct_mean": {names[x]: round(float(np.mean([pct[(r, x)] for r in runs])), 3) for x in pos},
        "C1p_pass_genes": [names[x] for x in c1_pass],
        "C1p": len(c1_pass) >= crit["C1p_min_genes"],
        "C2p_wins_of_9": c2,
        "C2p": all(w >= crit["C2p_min_wins"] for w in c2.values()) if c2 else False,
        "C3_seeds": c3, "C3_pct_mean": {names[h]: round(float(np.mean([pct[(r, h)] for r in runs])), 3) for h in hk},
        "C3": all(n >= crit["C3_min_seeds"] for n in c3.values()),
        "C3_panel_n_ok": f"{n_panel_ok}/{len(cp)}", "C3_panel": n_panel_ok >= crit["C3_panel_min_genes"],
        "C6_nonmarker_pass": c6_genes, "C6": len(c6_genes) > 0,
        "C7_seeds_low": c7,
        "extra_pct_mean": {names[x]: round(float(np.mean([pct[(r, x)] for r in runs])), 3) for x in extra},
        "null_frac_negative": round(float(np.mean([(g[(g.run == r) & (g.role == "negative")][metric] < 0).mean()
                                                   for r in runs])), 2),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--results", default=None, help="結果ディレクトリ（既定: config の output.results_dir）")
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    res = Path(a.results) if a.results else E / cfg["output"]["results_dir"]
    out = Path(a.out) if a.out else res / "eval"
    out.mkdir(parents=True, exist_ok=True)
    ann = {r["gene"]: r for r in json.loads((E / "work/gene_annotation.json").read_text())["genes"]}
    import pickle
    from geneformer import TOKEN_DICTIONARY_FILE
    with open(Path(TOKEN_DICTIONARY_FILE).with_name("gene_name_id_dict_gc104M.pkl"), "rb") as f:
        names = {v: k for k, v in pickle.load(f).items()}

    parts, descs = [], {}
    for run in sorted(pd.read_csv(res / "cells.csv", usecols=["run"]).run.unique()):
        cells, E0, E1, dP, cen = load_run(res, run, cfg["capture"].get("classifier_logits", False))
        cells, desc = vp.vector_metrics(cells, E0, E1, cen)
        alt = [x for x in cells.columns if x.startswith("shift_alt")][0]
        cells["alt"] = cells[alt]
        cells["pos"] = token_positions(cells)
        cells = rc.corrected(cells)
        if dP is not None:
            cells["dP_raw"] = dP
            bin_correct(cells, "dP_raw", "dP_a2", cfg["analysis"]["length_bins"])
        descs[run] = desc
        parts.append(cells)
    cells = pd.concat(parts, ignore_index=True)

    metrics = ["m_a2", "m_raw", "m_d", "m_d2", "m_c", "V2", "V2L"] + (["dP_a2", "dP_raw"] if "dP_raw" in cells else [])
    g = cells.groupby(["run", "gene", "role"])[metrics].mean().reset_index()
    g["gene_name"] = g.gene.map(names)
    rows = [assess(g, m, cfg, names, ann) for m in metrics]
    (out / "criteria.json").write_text(json.dumps(rows, ensure_ascii=False, indent=1, default=str) + "\n")
    df = pd.DataFrame(rows)
    for c in ["pos_pct_mean", "pos_seeds_ge_q90", "C2p_wins_of_9", "C3_pct_mean", "extra_pct_mean", "C7_seeds_low"]:
        df[c] = df[c].map(lambda d: " ".join(f"{k}:{v}" for k, v in d.items()))
    df["C1p_pass_genes"] = df.C1p_pass_genes.map(lambda x: ",".join(x) or "-")
    df["C6_nonmarker_pass"] = df.C6_nonmarker_pass.map(lambda x: ",".join(x) or "-")
    df.to_csv(out / "criteria.csv", index=False)
    (out / "drift.json").write_text(json.dumps(descs, indent=1) + "\n")
    g.to_csv(out / "gene_seed_metrics.csv", index=False)
    pd.set_option("display.width", 260); pd.set_option("display.max_columns", 40)
    print(f"config: {a.config}  runs: {sorted(cells.run.unique())}  cells: {len(cells)}  genes: {cells.gene.nunique()}")
    print(df[["metric", "pos_pct_mean", "C1p_pass_genes", "C1p", "C2p_wins_of_9", "C2p", "C3_pct_mean", "C3",
              "C3_panel_n_ok", "C3_panel", "C6_nonmarker_pass", "C6", "extra_pct_mean",
              "null_frac_negative"]].to_string(index=False))


if __name__ == "__main__":
    main()

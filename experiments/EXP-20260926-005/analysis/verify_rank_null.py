"""探索的（事前登録外）: 事前登録の判定（criteria.json）は変更しない。

元のトークン順位（rank, relpos = pos/length）交絡を差し引いた指標 m_a2_rank と
null 分布（rank マッチ null）で run O / run D を評価し直す検証スクリプト。
validation.md §4 の設計上の穴（overexpress で効きしろが効果量を決めるのに null を
そこにマッチさせていない）への事後検証で、事前登録の判定基準には影響しない。

方法:
  - 細胞表は evaluate_positives.py::main と同一に作る（ep.load_run → vp.vector_metrics
    → cells["alt"] → ep.token_positions → rc.corrected → relpos = pos/length）。
  - m_a2_rank = m_a2（長さビン補正済み）を relpos の 5 分位ビンで再度 null 平均を引く。
    ビン構造は ep.bin_correct と同一（各 null 遺伝子に leave-one-gene-out、最後に非 null
    全遺伝子 vs 全 null。境界は fit 側 null の relpos 分位、np.digitize）。
  - 分位（遺伝子単位平均、ep.quantiles 形式）:
      pct_full    … その run の全 null 遺伝子に対する分位（rc.null_pct）
      pct_matched … rank マッチ null（遺伝子単位 relpos に対する最近隣 k=10。target が
                     null 自身なら除外。マッチ数 <5 のとき pct_full を使い matched_n に記録）
       参考として補正前 m_a2 の pct_full も出力（既存 criteria.json との対応を見るため）。
   - 追加 variant（探索的方法開発。補正の作り方は null 側のみで決まる構造）:
       m_a2_rank_f … m_a2 を relpos の 10 分位ビン（nbins=10）で同じ LOO bin 補正。
       m_a2_dist   … m_a2 を pos（絶対トークン位置）の 5 分位ビン（境界は fit 側 null の pos 分位）で補正。
       マッチ距離の共変量は対応するビンの共変量（rank_f → relpos、dist → pos、k=10・自身除外・<5 退避）。
     診断も同様に pos 基準の gene-level Spearman・HK 分位・null 負割合集計を追加する。
  - 判定は config の analysis.criteria 閾値をそのまま使い ep.assess と同じロジック
    （C1'/C2'/C3/C3-panel/C6/C7）を pre-computed pct に対して実装（local 版 assess_pct）。
    盲検 rc.BLIND は ep.load_run が読み込み直後に除外するため追加実装不要（assert 残し）。

EXP-20260926-005 専用（ep.E を使う）。実行（repo root から）:
  .venv/bin/python experiments/EXP-20260926-005/analysis/verify_rank_null.py --config experiments/EXP-20260926-005/config-overexpress.yaml
  .venv/bin/python experiments/EXP-20260926-005/analysis/verify_rank_null.py --config experiments/EXP-20260926-005/config-delete.yaml
出力: <results>/eval/rank_null/{criteria_rank.json,criteria_rank.csv,gene_seed_metrics.csv,summary.json}
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import evaluate_positives as ep  # noqa: E402  (ep.E / ep.load_run / ep.token_positions / ep.vp / ep.rc を使う)

K_MATCH = 10          # rank マッチ null の最近隣数
MIN_MATCH = 5         # マッチ数が未満なら pct_full に退避


def null_pct_safe(x, null) -> float:
    """rc.null_pct と同定義（null のうち x 未満の割合）。NaN を除外。x が NaN なら NaN。"""
    null = np.asarray(null, float)
    null = null[~np.isnan(null)]
    if len(null) == 0 or (isinstance(x, float) and np.isnan(x)) or np.asarray(x).size != 1 or np.isnan(np.asarray(x, float)):
        return float("nan")
    return float((null < x).mean())


def mean_no_nan(vals) -> float:
    vals = np.asarray([v for v in vals], float)
    vals = vals[~np.isnan(vals)]
    return float(vals.mean()) if len(vals) else float("nan")


def bin_correct_rank(df: pd.DataFrame, base: str, out: str, cov: str, nbins: int = 5) -> np.ndarray:
    """ep.bin_correct と同じビン構造で、補正の共変量を relpos に変えた版（m_a2 → m_a2_rank）。"""
    v = np.full(len(df), np.nan)
    null_genes = sorted(df.loc[df.role == "negative", "gene"].unique())
    groups = [(df.gene == gg, (df.role == "negative") & (df.gene != gg)) for gg in null_genes]
    groups.append((df.role != "negative", df.role == "negative"))
    q = np.linspace(0, 1, nbins + 1)[1:-1]
    for tgt, fit in groups:
        F, T = df[fit], df[tgt]
        edges = np.quantile(F[cov], q)                       # fit 側 null の relpos 分位
        bf, bt = np.digitize(F[cov], edges), np.digitize(T[cov], edges)
        mu = F.groupby(bf)[base].mean()
        v[tgt.to_numpy()] = T[base].to_numpy() - mu.reindex(bt).to_numpy()
    df[out] = v
    return v


def build_pct(g: pd.DataFrame) -> pd.DataFrame:
    """遺伝子×run 表に分位列を付ける。
    既存列: pct_a2_full / pct_full / pct_matched / matched_n（m_a2_rank、relpos マッチ）。
    追加列（探索）: pct_a2rf_full / pct_a2rf_matched / matched_n_f（m_a2_rank_f、relpos マッチ）、
                    pct_a2d_full / pct_a2d_matched / matched_n_d（m_a2_dist、pos マッチ）。
    いずれもその run の全 null に対する分位（full）と最近隣マッチ null（k=10・自身除外）に対する分位（matched）。
    """
    def nearest(pool, cov, val):
        dist = (pool[cov] - val).abs().to_numpy()
        k = min(K_MATCH, len(pool))
        take = np.argpartition(dist, k)[:k] if k < len(pool) else np.arange(len(pool))
        return k, take
    rows = []
    for run, s in g.groupby("run"):
        nulls = s[s.role == "negative"]
        null_a2 = nulls["m_a2"].to_numpy()
        null_rank = nulls["m_a2_rank"].to_numpy()
        null_rank_f = nulls["m_a2_rank_f"].to_numpy()
        null_dist = nulls["m_a2_dist"].to_numpy()
        for _, r in s.iterrows():
            pool = nulls[nulls.gene != r.gene]               # target が null 自身のときは自身除外
            k, take = nearest(pool, "relpos", r["relpos"])
            if k < MIN_MATCH:                                # マッチ不足は pct_full 退避
                pm, pmf = r["pct_full"], null_pct_safe(r["m_a2_rank_f"], null_rank_f)
            else:
                pm = null_pct_safe(r["m_a2_rank"], pool["m_a2_rank"].to_numpy()[take])
                pmf = null_pct_safe(r["m_a2_rank_f"], pool["m_a2_rank_f"].to_numpy()[take])
            kd, take_d = nearest(pool, "pos", r["pos"])
            if kd < MIN_MATCH:                               # pos マッチも同様に退避
                pmd = null_pct_safe(r["m_a2_dist"], null_dist)
            else:
                pmd = null_pct_safe(r["m_a2_dist"], pool["m_a2_dist"].to_numpy()[take_d])
            rows.append({"run": run, "gene": r.gene,
                         "pct_a2_full": null_pct_safe(r["m_a2"], null_a2),
                         "pct_full": null_pct_safe(r["m_a2_rank"], null_rank),
                         "pct_matched": pm, "matched_n": k,
                         "pct_a2rf_full": null_pct_safe(r["m_a2_rank_f"], null_rank_f),
                         "pct_a2rf_matched": pmf, "matched_n_f": k,
                         "pct_a2d_full": null_pct_safe(r["m_a2_dist"], null_dist),
                         "pct_a2d_matched": pmd, "matched_n_d": kd})
    pct = pd.DataFrame(rows)
    return g.merge(pct, on=["run", "gene"], how="left")


def assess_pct(g: pd.DataFrame, cells: pd.DataFrame, value_col: str, pct_col: str, label: str,
               cfg: dict, names: dict, ann: dict) -> dict:
    """ep.assess と同じ判定ロジック（pre-computed pct を受け取る局所版）。閾値は config そのまま。"""
    crit = cfg["analysis"]["criteria"]
    pos = [x for x in cfg["controls"]["positive"]]
    hk = [x for x in cfg["controls"]["negative_fixed"]]
    roles = {grp["role"]: grp["genes"] for grp in cfg["controls"].get("extra_fixed", []) or []}
    panel, extra = roles.get("negctrl_panel", []), roles.get("prev_positive", []) + roles.get("reverse_direction", [])
    runs = sorted(g.run.unique())
    need = crit["C1_min_seeds"] if len(runs) >= 3 else len(runs)
    pct = {(r.run, r.gene): r[pct_col] for _, r in g.iterrows()}
    v = {(r.run, r.gene): r[value_col] for _, r in g.iterrows()}

    def n_hi(gene):
        return sum((not np.isnan(pct.get((r, gene), np.nan)) and pct[(r, gene)] >= crit["C1_positive_null_quantile"])
                   for r in runs)

    def pct_seq(keys):  # 全 NaN なら nan（np.mean の NaN 警告を避ける）
        return {names[x]: round(mean_no_nan([pct[(r, x)] for r in runs]), 3) for x in keys}

    c1_pass = [x for x in pos if n_hi(x) >= need]
    c2 = {names[x]: sum((not np.isnan(v[(r, x)])) and (not np.isnan(v[(r, h)])) and v[(r, x)] > v[(r, h)]
                        for r in runs for h in hk) for x in c1_pass}
    c3 = {names[h]: sum(crit["C3_hk_quantile_range"][0] <= pct.get((r, h), -1) <= crit["C3_hk_quantile_range"][1]
                        for r in runs) for h in hk}
    lo, hi = crit["C3_panel_quantile_range"]
    cp = {names[x]: sum(lo <= pct.get((r, x), -1) <= hi for r in runs) for x in panel}
    n_panel_ok = sum(n >= crit["C3_panel_min_seeds"] for n in cp.values())
    c6_genes = [names[x] for x in c1_pass if ann.get(names[x], {}).get("ddet_pct", 1.0) <= crit["C6_marker_pct_max"]]
    c7 = {names[x]: sum((not np.isnan(pct.get((r, x), np.nan))) and pct.get((r, x), 1.0) <= crit["C7_reverse_quantile_max"]
                        for r in runs) for x in roles.get("reverse_direction", [])}
    nfn = [float((cells[(cells.run == r) & (cells.role == "negative")][value_col] < 0).mean()) for r in runs]
    return {
        "metric": label, "value_col": value_col, "pct_col": pct_col,
        "n_seeds": len(runs),
        "pos_seeds_ge_q90": {names[x]: n_hi(x) for x in pos},
        "pos_pct_mean": pct_seq(pos),
        "C1p_pass_genes": [names[x] for x in c1_pass],
        "C1p": len(c1_pass) >= crit["C1p_min_genes"],
        "C2p_wins_of_9": c2,
        "C2p": all(w >= crit["C2p_min_wins"] for w in c2.values()) if c2 else False,
        "C3_seeds": c3, "C3_pct_mean": pct_seq(hk),
        "C3": all(n >= crit["C3_min_seeds"] for n in c3.values()),
        "C3_panel_n_ok": f"{n_panel_ok}/{len(cp)}", "C3_panel": n_panel_ok >= crit["C3_panel_min_genes"],
        "C3_panel_seeds_ok": {k: n for k, n in cp.items()},
        "C6_nonmarker_pass": c6_genes, "C6": len(c6_genes) > 0,
        "C7_seeds_low": c7,
        "extra_pct_mean": pct_seq(extra),
        "null_frac_negative": round(float(np.mean(nfn)), 2),
        "null_frac_negative_by_seed": {r: round(f, 3) for r, f in zip(runs, nfn)},
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    res = ep.E / cfg["output"]["results_dir"]
    out = res / "eval/rank_null"
    out.mkdir(parents=True, exist_ok=True)
    ann = {r["gene"]: r for r in json.loads((ep.E / "work/gene_annotation.json").read_text())["genes"]}
    from geneformer import TOKEN_DICTIONARY_FILE
    with open(Path(TOKEN_DICTIONARY_FILE).with_name("gene_name_id_dict_gc104M.pkl"), "rb") as f:
        names = {v: k for k, v in pickle.load(f).items()}
    assert not (set(cfg["controls"].get("positive", [])) | set(cfg["controls"].get("negative_fixed", []))) & ep.rc.BLIND

    parts = []
    for run in sorted(pd.read_csv(res / "cells.csv", usecols=["run"]).run.unique()):
        cells, E0, E1, dP, cen = ep.load_run(res, run, False)
        assert not set(cells.gene) & ep.rc.BLIND, "盲検の遺伝子が結果に含まれています"
        cells, _ = ep.vp.vector_metrics(cells, E0, E1, cen)
        cells["alt"] = cells[[x for x in cells.columns if x.startswith("shift_alt")][0]]
        cells["pos"] = ep.token_positions(cells)
        cells = ep.rc.corrected(cells)                       # m_a2 などを付ける
        cells["relpos"] = cells.pos / cells.length
        bin_correct_rank(cells, "m_a2", "m_a2_rank", "relpos", cfg["analysis"]["length_bins"])
        bin_correct_rank(cells, "m_a2", "m_a2_rank_f", "relpos", 10)                       # 探索: relpos 10 分位ビン
        bin_correct_rank(cells, "m_a2", "m_a2_dist", "pos", cfg["analysis"]["length_bins"])  # 探索: pos（絶対位置）5 分位ビン
        parts.append(cells)
    cells = pd.concat(parts, ignore_index=True)

    metrics = ["m_a2", "m_a2_rank", "m_a2_rank_f", "m_a2_dist", "relpos", "pos"]
    g = cells.groupby(["run", "gene", "role"])[metrics].mean().reset_index()
    g["gene_name"] = g.gene.map(names)
    g = build_pct(g)

    # --- 判定（variant ごとの seed 集約。閾値は事前登録 config そのまま・変更しない） ---
    variants = [("m_a2_pct_full（参考）", "m_a2", "pct_a2_full"),
                ("m_a2_rank_pct_full", "m_a2_rank", "pct_full"),
                ("m_a2_rank_pct_matched", "m_a2_rank", "pct_matched"),
                ("m_a2_rank_f_pct_full", "m_a2_rank_f", "pct_a2rf_full"),
                ("m_a2_rank_f_pct_matched", "m_a2_rank_f", "pct_a2rf_matched"),
                ("m_a2_dist_pct_full", "m_a2_dist", "pct_a2d_full"),
                ("m_a2_dist_pct_matched", "m_a2_dist", "pct_a2d_matched")]
    rows = [assess_pct(g, cells, vcol, pcol, lab, cfg, names, ann) for lab, vcol, pcol in variants]

    # --- 診断（guardrails 対応: 符号・向きの確認。relpos 基準と pos 基準の両方） ---
    diag = {}
    rho, rho_pos = {}, {}
    for run, s in g.groupby("run"):
        nul = s[s.role == "negative"]
        e_r, e_p = {}, {}
        for m in ("m_a2", "m_a2_rank", "m_a2_rank_f", "m_a2_dist"):
            for tag, sub in (("all", s), ("null20", nul)):
                rr = stats.spearmanr(sub.relpos, sub[m])
                e_r[f"{m}_{tag}"] = {"rho": round(float(rr.statistic), 3), "p": float(rr.pvalue)}
                rp = stats.spearmanr(sub.pos, sub[m])
                e_p[f"{m}_{tag}"] = {"rho": round(float(rp.statistic), 3), "p": float(rp.pvalue)}
        rho[run], rho_pos[run] = e_r, e_p
    nfn = {}
    for run, s in cells.groupby("run"):
        nul = s[s.role == "negative"]
        nfn[run] = {m: round(float((nul[m] < 0).mean()), 3)
                    for m in ("m_a2", "m_a2_rank", "m_a2_rank_f", "m_a2_dist")}
    hk = list(cfg["controls"]["negative_fixed"])

    def seed_pct(gene, col):
        return {r.run: (None if np.isnan(r[col]) else round(float(r[col]), 3))
                for _, r in g[g.gene == gene].iterrows()}
    hk_pct = {names[h]: {"pct_full": seed_pct(h, "pct_full"), "pct_matched": seed_pct(h, "pct_matched"),
                         "pct_rank_f_full": seed_pct(h, "pct_a2rf_full"), "pct_rank_f_matched": seed_pct(h, "pct_a2rf_matched"),
                         "pct_dist_full": seed_pct(h, "pct_a2d_full"), "pct_dist_matched": seed_pct(h, "pct_a2d_matched")}
              for h in hk}

    def mstats(col):
        a = g[col].to_numpy()
        return {"n_entries": int(len(a)), "min": int(np.min(a)), "median": float(np.median(a))}
    summary = {"config": str(a.config), "rho_relpos_gene": rho, "rho_pos_gene": rho_pos,
              "null_frac_negative_by_seed": nfn, "hk_pct": hk_pct,
              "matched_null_n": mstats("matched_n"), "matched_null_n_rank_f": mstats("matched_n_f"),
              "matched_null_n_dist": mstats("matched_n_d"),
              "nan_m_a2_rank_genes": sorted(g.loc[g.m_a2_rank.isna(), "gene_name"].unique().tolist()),
              "nan_m_a2_rank_f_genes": sorted(g.loc[g.m_a2_rank_f.isna(), "gene_name"].unique().tolist()),
              "nan_m_a2_dist_genes": sorted(g.loc[g.m_a2_dist.isna(), "gene_name"].unique().tolist()),
              "note": "探索的（事前登録外）。判定は参考表示のみで criteria.json（事前登録）は変更しない。"}

    # --- 出力 ---
    (out / "criteria_rank.json").write_text(json.dumps(rows, ensure_ascii=False, indent=1, default=str) + "\n")
    df = pd.DataFrame(rows)
    for c in ["pos_pct_mean", "pos_seeds_ge_q90", "C2p_wins_of_9", "C3_pct_mean", "extra_pct_mean",
              "C7_seeds_low", "C3_panel_seeds_ok", "null_frac_negative_by_seed"]:
        df[c] = df[c].map(lambda d: " ".join(f"{k}:{v}" for k, v in d.items()))
    df["C1p_pass_genes"] = df.C1p_pass_genes.map(lambda x: ",".join(x) or "-")
    df["C6_nonmarker_pass"] = df.C6_nonmarker_pass.map(lambda x: ",".join(x) or "-")
    df.to_csv(out / "criteria_rank.csv", index=False)
    g.to_csv(out / "gene_seed_metrics.csv", index=False)
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=1) + "\n")

    pd.set_option("display.width", 300); pd.set_option("display.max_columns", 40)
    print(f"config: {a.config}  runs: {sorted(cells.run.unique())}  cells: {len(cells)}  genes: {cells.gene.nunique()}")
    print(df[["metric", "pos_pct_mean", "C1p_pass_genes", "C1p", "C2p_wins_of_9", "C2p", "C3_pct_mean", "C3",
              "C3_panel_n_ok", "C3_panel", "C6_nonmarker_pass", "C6", "extra_pct_mean",
              "null_frac_negative", "null_frac_negative_by_seed"]].to_string(index=False))
    print(json.dumps(summary, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()

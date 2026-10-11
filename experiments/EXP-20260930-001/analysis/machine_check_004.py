"""EXP-20260930-001 — 2 台の計算機が同じ答えを出すことを、本番の前に確かめる（plan.md §8.2）。

構成:
  spark-cf5b  petadimensionlab@100.104.171.109  (seed 81 担当)
  spark-c083  coi-next1@100.94.57.17            (seed 82 → 83 担当)

48 遺伝子（T1 11 + T2 7 + 陽性対照 10 + 恒常発現 3 + 中立 7 + 対照 10）× 200 細胞 × seed 81
（`work/config-machinecheck.yaml`。**事前登録の対象外**）。

判定基準（plan.md §8.2。**結果を見る前に固定**）:
  H-1（主）  48 遺伝子の標準偏差単位 z について、機体間の差の
             平均（系統成分）の絶対値 <= 0.02、かつ 最大絶対差 <= 0.05
  H-2       細胞ごとの主指標 m_goal の 2 台間 Spearman 順位相関 rho >= 0.99
  H-3       状態重心ベクトル 4 本（FIB_MFIB / aStr / dFIB / dMFIB）の 2 台間コサイン >= 0.9999

不一致だったときの扱い（結果を見る前に決める）:
  (a) まず transformers を 4.46.3 にそろえて取り直す
  (b) そろえられないなら spark-c083 1 台で 3 seed すべてを回す（11.7 時間）
  (c) 「機体を要因として記録するだけ」は採らない（seed と機体が交絡していて分離できない）

実行: .venv/bin/python experiments/EXP-20260930-001/analysis/machine_check_004.py --compare DIR_A DIR_B
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260930-001"
SEED = 81
CRIT = {"H1_mean_abs_max": 0.02, "H1_max_abs_max": 0.05, "H2_rho_min": 0.99, "H3_cos_min": 0.9999}


def _load_evaluate():
    spec = importlib.util.spec_from_file_location("ev004", Path(__file__).resolve().parent / "evaluate_004.py")
    ev = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ev)
    return ev


def compare(a: Path, b: Path) -> None:
    ev = _load_evaluate()
    cfg = yaml.safe_load((E / "work/config-machinecheck.yaml").read_text())
    gj = json.loads((E / "work/genes_004.json").read_text())["genes"]
    sym = {g["ensembl"]: g["symbol"] for g in gj}

    out: dict = {"criteria": CRIT, "dirs": [str(a), str(b)]}
    cells, tabs = {}, {}
    for tag, d in (("A", a), ("B", b)):
        cfg2 = json.loads(json.dumps(cfg))
        cfg2["output"]["results_dir"] = str(d.relative_to(E))
        cs, drop = ev.load(cfg2)
        cols = ev.metric_columns(cs, list(cfg2["states"].get("alt_states", [])))
        cells[tag] = cs
        nulls_pool = json.loads((d / f"run-{SEED}-L-1/negatives.json").read_text())
        nl = {h: hd for h, hd in cs[cs.gene.isin(nulls_pool)].groupby("gene")}
        pm = cfg2["analysis"]["pos_matched_permutation"]
        el = ev.eligibility(cs, nl, pm["tolerance_log2"], pm["min_cover_frac"])
        from movement_correction import gene_stats
        t = gene_stats(cs, nl, el, cols["goal"], pm["min_eligible_null_genes"], use_dnorm=True)
        mu, sd = t[t.gene.isin(nulls_pool)].D_raw.mean(), t[t.gene.isin(nulls_pool)].D_raw.std(ddof=1)
        t["z"] = (t.D_raw - mu) / sd
        tabs[tag] = t.set_index("gene")
        out[f"drop_{tag}"] = drop

    g = sorted(set(tabs["A"].index) & set(tabs["B"].index))
    dz = (tabs["A"].loc[g, "z"] - tabs["B"].loc[g, "z"]).to_numpy()
    out["H1"] = {"n_genes": len(g), "mean_abs": float(abs(dz.mean())), "max_abs": float(np.abs(dz).max()),
                 "worst_genes": [sym.get(x, x) for x in np.array(g)[np.argsort(-np.abs(dz))[:5]]],
                 "pass": bool(abs(dz.mean()) <= CRIT["H1_mean_abs_max"]
                              and np.abs(dz).max() <= CRIT["H1_max_abs_max"])}

    ka = cells["A"].set_index(["gene", "cell_uid"]).m_goal
    kb = cells["B"].set_index(["gene", "cell_uid"]).m_goal
    idx = ka.index.intersection(kb.index)
    rho = float(stats.spearmanr(ka.loc[idx], kb.loc[idx]).statistic)
    out["H2"] = {"n_cells": int(len(idx)), "spearman_rho": rho, "pass": bool(rho >= CRIT["H2_rho_min"])}

    cs_names = {}
    for tag, d in (("A", a), ("B", b)):
        sz = np.load(d / f"run-{SEED}-L-1/state_emb_vectors.npz")
        cs_names[tag] = {k: sz[k].astype(np.float64) for k in sz.files}
    common = set(cs_names["A"]) & set(cs_names["B"])
    coss = {k: float(cs_names["A"][k] @ cs_names["B"][k]
                     / (np.linalg.norm(cs_names["A"][k]) * np.linalg.norm(cs_names["B"][k])))
            for k in common}
    out["H3"] = {"cosine": coss, "min": min(coss.values()), "pass": bool(min(coss.values()) >= CRIT["H3_cos_min"])}

    dn = {tag: float(cells[tag].dnorm.median()) for tag in ("A", "B")}
    out["H4_descriptive"] = {"移動量 ||Δ|| の中央値": dn, "比": round(dn["A"] / max(dn["B"], 1e-12), 4)}

    out["pass"] = bool(out["H1"]["pass"] and out["H2"]["pass"] and out["H3"]["pass"])
    out["不一致のときの扱い"] = ("(a) transformers を 4.46.3 にそろえて取り直す → "
                                 "(b) そろえられないなら spark-c083 1 台で 3 seed すべてを回す（11.7 時間）。"
                                 "(c)「機体を要因として記録するだけ」は採らない（seed と交絡して分離できない）")
    (E / "work/machine_check_004.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--compare", nargs=2, metavar=("DIR_A", "DIR_B"))
    args = ap.parse_args()
    if args.compare:
        compare(Path(args.compare[0]), Path(args.compare[1]))
    else:
        ap.error("--compare DIR_A DIR_B を指定してください")

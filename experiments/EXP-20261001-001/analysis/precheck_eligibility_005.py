#!/usr/bin/env python3
"""EXP-20261001-001 事前チェック — 判定する遺伝子ごとに「比べられる対照遺伝子」が
下限 20 個以上あるかを**摂動の前に**確かめる（guardrails 2026-09-27）。

適格性は**入力中の位置（lp）と入力の長さ（ll）だけ**で決まり、どちらも摂動の前に確定しているので、
ここで正確に計算できる。判定に使うのと**同じ関数**（evaluate_005.eligibility）を呼ぶ。

合否（**結果を見る前に固定**）:
  (1) **介入セット 119 個と前回の摂動対象 T1・T2 の、合わせて 135 遺伝子すべて**で、
      適格な対照が 20 個以上（3 seed とも）。
      満たさなければ着手しない（その遺伝子は構造的に判定不能になる）。
  (2) それ以外の役割（陽性対照など）で下限を割る遺伝子は、**列挙して記録する**。
      その遺伝子の単一遺伝子の順位 q は**欠測**として扱い、「不合格」とは数えない
      （guardrails 2026-09-27）。セットの検定はセット平均の z（D_raw から計算）で行うので
      セットには残す。ただし残差が外挿になるので、**その遺伝子を除いたセットでの再判定を
      入れ子の感度分析として事前登録する**。

使い方: .venv/bin/python analysis/precheck_eligibility_005.py [A_300|B_500]
出力:   work/precheck_eligibility_005-{版}.json
"""
from __future__ import annotations

import importlib.util
import json
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20261001-001"
sys.path.insert(0, str(ROOT / "scripts/geneformer"))
import run_isp as ri  # noqa: E402

spec = importlib.util.spec_from_file_location("ev5", E / "analysis/evaluate_005.py")
ev5 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ev5)


def main() -> None:
    from datasets import load_from_disk
    cfg = yaml.safe_load(open(E / "config.yaml"))
    pm = cfg["analysis"]["pos_matched_permutation"]
    seeds = cfg["robustness"]["seeds"]
    sel = cfg["cell_selection"]

    gj = json.loads((E / "work/genes_005.json").read_text())
    roles = {g["ensembl"]: g["role"] for g in gj["genes"]}
    npool = json.loads((E / "work/null_pool_005.json").read_text())["versions"]
    ver = sys.argv[1] if len(sys.argv) > 1 else "B_500"
    nulls = [x["ensembl"] for x in npool[ver]]
    judged = [e for e, r in roles.items() if r != "negative"]
    allg = judged + nulls

    tokd = ri.gene_token_dict()
    rcfg = {"_path": str(E / "config.yaml"), "states": cfg["states"], "data": cfg["data"],
            "cell_selection": sel}
    m, rows, donors = ri.start_matrix(E / "work/tokenized.dataset", rcfg)
    ds = load_from_disk(str(E / "work/tokenized.dataset"))
    ids = {int(r): np.asarray(x) for r, x in zip(rows, ds.select(rows.tolist())["input_ids"])}

    out = {"null_version": ver, "min_eligible_required": pm["min_eligible_null_genes"],
           "n_nulls": len(nulls), "n_judged": len(judged),
           "tolerance_log2": pm["tolerance_log2"], "min_cover_frac": pm["min_cover_frac"],
           "per_seed": {}}
    worst = []
    for s in seeds:
        recs = []
        for g in allg:
            t = tokd[g]
            idx = m[:, t].nonzero()[0].tolist()
            keep = ri.balanced_subset(idx, donors, sel["max_cells_per_gene"],
                                      sel["max_cells_per_donor"], random.Random(f"{s}-{g}"))
            for i in keep:
                a = ids[int(rows[i])]
                p = int(np.flatnonzero(a == t)[0])
                recs.append({"gene": g, "lp": float(np.log2(1 + p)), "ll": float(np.log2(len(a)))})
        cs = pd.DataFrame(recs)
        nl = {h: hd for h, hd in cs[cs.gene.isin(set(nulls))].groupby("gene")}
        el = ev5.eligibility(cs, nl, pm["tolerance_log2"], pm["min_cover_frac"])
        nj = {g: len(el[g]) for g in judged}
        out["per_seed"][str(s)] = {
            "min_eligible_judged": int(min(nj.values())),
            "median_eligible_judged": float(np.median(list(nj.values()))),
            "genes_below_min": sorted(k for k, v in nj.items()
                                      if v < pm["min_eligible_null_genes"]),
            "min_eligible_all": int(min(len(v) for v in el.values())),
            "per_gene_judged": nj}
        worst.append(min(nj.values()))
        print(f"seed {s}: 判定遺伝子の適格な対照 最小 {min(nj.values())} / "
              f"中央値 {np.median(list(nj.values())):.0f}、"
              f"下限割れ {out['per_seed'][str(s)]['genes_below_min']}", flush=True)
    targets = [e for e, r in roles.items() if r in ("intervention", "target_T1", "target_T2")]
    out["targets_below_min"] = sorted(
        {g for s in seeds for g in out["per_seed"][str(s)]["genes_below_min"] if g in set(targets)})
    out["others_below_min"] = sorted(
        {g for s in seeds for g in out["per_seed"][str(s)]["genes_below_min"]}
        - set(out["targets_below_min"]))
    out["min_eligible_targets"] = int(min(
        out["per_seed"][str(s)]["per_gene_judged"][g] for s in seeds for g in targets))
    out["pass"] = not out["targets_below_min"]          # **着手の条件は (1) だけ**
    out["worst_over_seeds"] = int(min(worst))
    (E / f"work/precheck_eligibility_005-{ver}.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "per_seed"}, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()

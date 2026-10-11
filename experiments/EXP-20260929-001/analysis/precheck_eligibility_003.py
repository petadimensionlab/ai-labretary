#!/usr/bin/env python3
"""EXP-20260929-001 事前チェック P-C: 対照 300 個で、判定する遺伝子ごとに
「比べられる対照遺伝子」が下限 20 個以上あるかを**摂動の前に**確かめる。

guardrails 2026-09-27:「過剰発現型の対象にできるかは発現細胞数だけでは決まらない。
元の位置の広がりが大きい遺伝子は、そろえられる相手が足りず構造的に判定不能になる」。
適格性は**入力中の位置（lp）と入力長（ll）だけ**で決まり、どちらも摂動の前に確定しているので、
ここで正確に計算できる。判定に使うのと**同じ関数**（evaluate_003.eligibility）を呼ぶ。

合否（結果を見る前に固定）: 判定する 67 遺伝子すべてで、適格な対照が **20 個以上**（3 seed とも）。
満たさない遺伝子があれば、対照を 500 に増やす（入れ子なので設計は変えずに済む）。

出力: work/precheck_eligibility_003.json
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
E = ROOT / "experiments/EXP-20260929-001"
sys.path.insert(0, str(ROOT / "scripts/geneformer"))
import run_isp as ri  # noqa: E402

spec = importlib.util.spec_from_file_location("ev3", E / "analysis/evaluate_003.py")
ev3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ev3)


def main() -> None:
    from datasets import load_from_disk
    cfg = yaml.safe_load(open(E / "config.yaml"))
    pm = cfg["analysis"]["pos_matched_permutation"]
    seeds = cfg["robustness"]["seeds"]
    sel = cfg["cell_selection"]

    gj = json.loads((E / "work/genes_003.json").read_text())
    ens = {g["symbol"]: g["ensembl"] for g in gj["genes"]}
    roles = {g["ensembl"]: g["role"] for g in gj["genes"]}
    npool = json.loads((E / "work/null_pool_003.json").read_text())["versions"]
    ver = sys.argv[1] if len(sys.argv) > 1 else "A_300"
    nulls300 = [x["ensembl"] for x in npool[ver]]
    judged = [e for e, r in roles.items() if r != "negative"]
    allg = judged + nulls300

    tokd = ri.gene_token_dict()
    rcfg = {"_path": str(E / "config.yaml"), "states": cfg["states"], "data": cfg["data"],
            "cell_selection": sel}
    m, rows, donors = ri.start_matrix(E / "work/tokenized.dataset", rcfg)
    ds = load_from_disk(str(E / "work/tokenized.dataset"))
    ids = {int(r): np.asarray(x) for r, x in zip(rows, ds.select(rows.tolist())["input_ids"])}

    out = {"null_version": ver, "min_eligible_required": pm["min_eligible_null_genes"], "n_nulls": len(nulls300),
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
        nl = {h: hd for h, hd in cs[cs.gene.isin(set(nulls300))].groupby("gene")}
        el = ev3.eligibility(cs, nl, pm["tolerance_log2"], pm["min_cover_frac"])
        nj = {g: len(el[g]) for g in judged}
        out["per_seed"][str(s)] = {
            "min_eligible_judged": int(min(nj.values())),
            "median_eligible_judged": float(np.median(list(nj.values()))),
            "genes_below_min": sorted(k for k, v in nj.items() if v < pm["min_eligible_null_genes"]),
            "min_eligible_all": int(min(len(v) for v in el.values())),
            "per_gene_judged": nj}
        worst.append(min(nj.values()))
    out["pass"] = all(not out["per_seed"][str(s)]["genes_below_min"] for s in seeds)
    out["worst_over_seeds"] = int(min(worst))
    (E / f"work/precheck_eligibility_003-{ver}.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "per_seed"}, ensure_ascii=False, indent=1))
    for s in seeds:
        d = out["per_seed"][str(s)]
        print(f"seed {s}: 判定遺伝子の適格な対照 最小 {d['min_eligible_judged']} / "
              f"中央値 {d['median_eligible_judged']:.0f}、下限割れ {d['genes_below_min']}")


if __name__ == "__main__":
    main()

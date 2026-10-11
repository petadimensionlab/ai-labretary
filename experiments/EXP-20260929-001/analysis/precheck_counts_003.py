#!/usr/bin/env python3
"""EXP-20260929-001 事前チェック P-A: **トークン化済み入力**での摂動対象の細胞数を数え直す。

HYP-20260929-01 §3 の数は「生カウントで 1 以上」から数えたもので、Geneformer の
rank 符号化と入力長 4,096 での切り捨てを考慮していない。**生カウントの検出は上限であって実測ではない。**
ここでは run_isp.py の `start_matrix()` と `balanced_subset()` をそのまま使い、
本番とまったく同じ抽出規則（提供者あたり 20 / 合計 200 / seed 61・62・63）で数える。

合否（**結果を見る前に固定した規則**。plan.md §5.1）:
  - 3 seed すべてで抽出細胞数 ≥ 120 → 合格
  - 分泌プログラム 25 のうち合格しなかった遺伝子は、その遺伝子をセットから外す
  - ただし **CXCL13 / CCL19 / CCL21 のいずれかが不合格なら、外さずに PI へ諮り直す**
    （この 3 遺伝子を摂動できることが本設計の目的そのものであるため。HYP §10 の危険）
  - 残った K が **20 未満**なら設計を諮り直す（検出下限が 1.1 標準偏差を超えるため）
  - 目標マーカー 7 個（陽性対照）で不合格が出たら、その遺伝子を外して K を減らす。
    **7 個すべてが不合格なら目標の軸の陽性対照が無くなるので諮り直す**

出力: work/precheck_counts_003.json
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20260929-001"
sys.path.insert(0, str(ROOT / "scripts/geneformer"))

import run_isp as ri  # noqa: E402

SEEDS = [61, 62, 63]
MIN_CELLS = 120
MAX_CELLS = 200
MAX_PER_DONOR = 20

CFG = {
    "_path": str(E / "config.yaml"),
    "states": {"start_state": "FIB_orgLow"},
    "data": {"state_key": "fib_org_state"},
    "cell_selection": {"donor_key": "donor_id", "max_cells_per_gene": MAX_CELLS,
                       "max_cells_per_donor": MAX_PER_DONOR, "min_cells_per_gene": MIN_CELLS},
}
TOK = E / "work/tokenized.dataset"


def main() -> None:
    import numpy as np

    groups = json.loads((E / "analysis/excluded_symbols_003.json").read_text())["groups"]
    sym2ens = json.loads((E / "work/symbol_to_ensembl_003.json").read_text())
    tokd = ri.gene_token_dict()
    m, rows, donors = ri.start_matrix(TOK, CFG)

    def count(sym: str) -> dict:
        ens = sym2ens.get(sym)
        tid = tokd.get(ens) if ens else None
        if tid is None or tid >= m.shape[1]:
            return {"symbol": sym, "ensembl": ens, "token": None, "in_dict": False,
                    "n_detected_tokenized": 0, "n_donors": 0, "per_seed": {}, "min_over_seeds": 0,
                    "pass": False}
        idx = m[:, tid].nonzero()[0].tolist()
        per = {}
        for s in SEEDS:
            keep = ri.balanced_subset(idx, donors, MAX_CELLS, MAX_PER_DONOR,
                                      random.Random(f"{s}-{ens}"))
            per[str(s)] = {"n_cells": len(keep),
                           "n_donors": len({donors[i] for i in keep})}
        mn = min(v["n_cells"] for v in per.values())
        return {"symbol": sym, "ensembl": ens, "token": int(tid), "in_dict": True,
                "n_detected_tokenized": len(idx),
                "n_donors": len({donors[i] for i in idx}),
                "per_seed": per, "min_over_seeds": mn, "pass": bool(mn >= MIN_CELLS)}

    res = {"start_state": "FIB_orgLow", "n_start_cells": int(m.shape[0]),
           "n_start_donors": len(set(donors)), "seeds": SEEDS,
           "rule": {"min_cells_per_gene": MIN_CELLS, "max_cells_per_gene": MAX_CELLS,
                    "max_cells_per_donor": MAX_PER_DONOR,
                    "require": "3 seed すべてで 120 細胞以上"},
           "sets": {}}
    for name in ["secretory25", "org_markers7", "myof_markers24", "housekeeping3",
                 "negctrl_panel7", "reference2"]:
        res["sets"][name] = [count(s) for s in groups[name]]

    sec = res["sets"]["secretory25"]
    chemo = {r["symbol"]: r for r in sec if r["symbol"] in ("CXCL13", "CCL19", "CCL21")}
    passed = [r["symbol"] for r in sec if r["pass"]]
    dropped = [r["symbol"] for r in sec if not r["pass"]]
    res["verdict"] = {
        "secretory_pass": passed, "secretory_dropped": dropped, "K_after": len(passed),
        "chemokine3": {k: v["min_over_seeds"] for k, v in chemo.items()},
        "chemokine3_all_pass": all(v["pass"] for v in chemo.values()),
        "org_markers7_pass": [r["symbol"] for r in res["sets"]["org_markers7"] if r["pass"]],
        "org_markers7_dropped": [r["symbol"] for r in res["sets"]["org_markers7"] if not r["pass"]],
        "myof24_pass": [r["symbol"] for r in res["sets"]["myof_markers24"] if r["pass"]],
        "myof24_dropped": [r["symbol"] for r in res["sets"]["myof_markers24"] if not r["pass"]],
    }
    v = res["verdict"]
    if not v["chemokine3_all_pass"]:
        v["action"] = "STOP: CXCL13/CCL19/CCL21 のいずれかが 120 細胞に届かない → PI に諮り直す"
    elif v["K_after"] < 20:
        v["action"] = f"STOP: K = {v['K_after']} < 20 → 設計を諮り直す"
    elif not v["org_markers7_pass"]:
        v["action"] = "STOP: 目標の軸の陽性対照が 1 つも立たない → 諮り直す"
    else:
        v["action"] = "GO"
    (E / "work/precheck_counts_003.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False) + "\n")

    print(f"出発点 {res['n_start_cells']} 細胞 / {res['n_start_donors']} 提供者\n")
    for name, rr in res["sets"].items():
        print(f"--- {name} ---")
        for r in rr:
            print(f"  {r['symbol']:10s} 生→トークン化 {r['n_detected_tokenized']:6d} 細胞 / "
                  f"{r['n_donors']:3d} 提供者  上限後 min(3 seed) = {r['min_over_seeds']:4d}  "
                  f"{'○' if r['pass'] else '×'}")
    print("\n判定:", json.dumps(v, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()

"""EXP-20260927-001 評価コードの検査（合成データ, ISP の結果は使わない）。
元の位置と長さに強く依存する見かけの効果（EXP-005 run O で見られた型の交絡）を入れた合成 cells.csv を作り、
evaluate_matched.py が (a) 交絡だけのときに陽性対照を合格させない、(b) 真の効果を 2 遺伝子に入れたときに検出する、
(c) 較正（null 遺伝子自身の分位）が保たれる、を確認する。
実行: .venv/bin/python experiments/EXP-20260927-001/analysis/synthetic_check.py
出力: work/synthetic/{null,effect}/cells.csv と eval/
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260927-001"
CFG = yaml.safe_load(open(E / "config.yaml"))


def make(effect: bool, seed: int = 7) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    ctl = CFG["controls"]
    genes = list(ctl["positive"]) + list(ctl["negative_fixed"]) + list(ctl["negative_sampler"]["genes"]) + \
        [x for grp in ctl["extra_fixed"] for x in grp["genes"]]
    donors = [f"D{i:02d}" for i in range(14)]
    donor_eff = dict(zip(donors, rng.normal(0, 2e-4, 14)))
    gene_eff = dict(zip(genes, rng.normal(0, 3e-4, len(genes))))
    gene_rel = dict(zip(genes, rng.uniform(0.15, 0.9, len(genes))))       # 遺伝子ごとの典型的な相対位置
    gene_len = dict(zip(genes, rng.uniform(10.5, 11.8, len(genes))))       # 遺伝子ごとの典型的な log2 長さ
    for hk, rel in zip(ctl["negative_fixed"], (0.95, 0.86, 0.55)):         # 実データの HK に近い相対位置
        gene_rel[hk] = rel
    true = set(ctl["positive"][:2]) if effect else set()
    rows, uid = [], 0
    for s in CFG["robustness"]["seeds"]:
        for g in genes:
            n = CFG["cell_selection"]["max_cells_per_gene"]
            L = np.clip(2 ** rng.normal(gene_len[g], 0.4, n), 300, 4096).astype(int)
            rel = np.clip(rng.beta(4 * gene_rel[g], 4 * (1 - gene_rel[g]), n), 0.001, 0.999)
            pos = np.maximum(1, (rel * L).astype(int))
            d = rng.choice(donors, n)
            # 交絡: 元の位置が低い（pos が大きい）ほど goal 方向へ大きく動く + 長さの効果
            y = 2.7e-3 + 2.5e-3 * (np.log2(1 + pos) - 10) / 2 + 1.0e-3 * (np.log2(L) - 11) \
                + gene_eff[g] + np.array([donor_eff[x] for x in d]) + rng.normal(0, 3e-3, n)
            if g in true:
                y = y + 1.5e-3
            for k in range(n):
                rows.append({"run": f"run-{s}-L-1", "gene": g, "cell_uid": uid, "donor_id": d[k], "length": L[k],
                             "shift_goal": y[k], "pos": pos[k]})
                uid += 1
    return pd.DataFrame(rows)


def main() -> None:
    for tag, eff in (("null", False), ("effect", True)):
        d = E / "work/synthetic" / tag
        d.mkdir(parents=True, exist_ok=True)
        make(eff).to_csv(d / "cells.csv", index=False)
        print(f"== {tag}", flush=True)
        subprocess.run([sys.executable, str(E / "analysis/evaluate_matched.py"), "--config", str(E / "config.yaml"),
                        "--cells", str(d / "cells.csv"), "--out", str(d / "eval"), "--no-annotation"], check=True)


if __name__ == "__main__":
    main()

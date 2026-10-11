#!/usr/bin/env python3
"""EXP-20261007-001 §5.3.1 — A2 で bfloat16 と float32 の食い違いを標準得点の単位で実測し、
判定を保留する範囲の幅を、結果を見る前に固定した規則（観測された最大の 2.00 倍）で決める。

**摂動は回さない。** cf5b で回した 2 つの結果（results-dtype-bf16 / results-dtype-fp32）を読むだけ。
比べる前に scripts/analysis/compare_conditions.py を通し、終了コード 1（記録なし）なら止める。

走らせ方:
    .venv/bin/python experiments/EXP-20261007-001/analysis/dtype_resolution_011.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261007-001"
sys.path.insert(0, str(ROOT / "experiments/EXP-20261004-001/analysis"))
sys.path.insert(0, str(ROOT / "experiments/EXP-20261001-001/analysis"))
sys.path.insert(0, str(ROOT / "scripts/analysis"))
import evaluate_009 as ev9  # noqa: E402

FACTOR = 2.00   # 幅 = 観測された最大の 2.00 倍（plan.md §5.3.1。結果を見る前に固定）


def d_raw(name: str) -> tuple[pd.DataFrame, dict]:
    cfg = yaml.safe_load((E / f"work/config-dtype-{name}.yaml").read_text())
    main_cfg = yaml.safe_load((E / "config.yaml").read_text())
    ev9.E = E
    ev9.BLOCKED = {e: e for e in main_cfg["blocked_genes"]}
    c, drop = ev9.load_cells(cfg, f"results-dtype-{name}")
    genes = list(cfg["controls"]["negative_sampler"]["genes"])
    frames = []
    for s in cfg["robustness"]["seeds"]:
        cs = c[c.seed == s].assign(m_goal=lambda d: d.m_proj)
        nl = {h: hd for h, hd in cs[cs.gene.isin(genes)].groupby("gene")}
        frames.append(ev9.gene_stats_stratified(cs, nl, "m_goal").assign(seed=s))
    return pd.concat(frames, ignore_index=True), drop


def main() -> None:
    out: dict = {"実験": "EXP-20261007-001", "節": "plan.md §5.3.1", "規則": f"幅 = 観測された最大の食い違い × {FACTOR}"}
    r = subprocess.run([sys.executable, str(ROOT / "scripts/analysis/compare_conditions.py"),
                        str(E / "results-dtype-bf16"), str(E / "results-dtype-fp32")], capture_output=True, text=True)
    out["実行条件の突き合わせ"] = {"終了コード": r.returncode, "出力": r.stdout[-3000:]}
    if r.returncode == 1:
        (E / "work/dtype_resolution_011.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
        raise SystemExit("**止める**: 実行条件の記録が無く、突き合わせられない（終了コード 1）")
    bf, dbf = d_raw("bf16")
    fp, dfp = d_raw("fp32")
    m = bf.merge(fp, on=["gene", "seed"], suffixes=("_bf16", "_fp32"))
    if len(m) != 36 or m.n_cells_bf16.ne(m.n_cells_fp32).any():
        raise SystemExit(f"**止める**: 突き合わせられた遺伝子×群が 36 でない、または細胞数が違う（{len(m)}）")
    sd = fp.groupby("seed").D_raw.std(ddof=1)
    m["食い違い"] = (m.D_raw_bf16 - m.D_raw_fp32).abs() / m.seed.map(sd)
    mx = float(m["食い違い"].max())
    sym = {g["ensembl"]: g["symbol"] for g in json.loads((E / "work/dtype_genes_011.json").read_text())["genes"]}
    out |= {
        "移動量 0 で外した細胞": {"bf16": dbf, "fp32": dfp},
        "群ごとの 12 個の D_raw(float32) の標準偏差": {str(k): float(v) for k, v in sd.items()},
        "食い違い（標準得点の単位）": {"個数": int(len(m)), "平均": round(float(m["食い違い"].mean()), 4),
                                    "標準偏差": round(float(m["食い違い"].std(ddof=1)), 4), "最大": round(mx, 4),
                                    "群ごとの最大": {str(k): round(float(v), 4)
                                                for k, v in m.groupby("seed")["食い違い"].max().items()}},
        "決めた幅": round(FACTOR * mx, 4),
        "参考: A1 の値（EXP-20261004-001）": {"平均": 0.194, "標準偏差": 0.243, "最大": 0.436, "幅": 0.87},
        "遺伝子ごと": [{"symbol": sym.get(g, g), "群": int(s), "食い違い": round(float(v), 4)}
                    for g, s, v in zip(m.gene, m.seed, m["食い違い"])]}
    (E / "work/dtype_resolution_011.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "遺伝子ごと"}, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()

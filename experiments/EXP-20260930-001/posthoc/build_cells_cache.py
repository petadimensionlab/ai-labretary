#!/usr/bin/env python3
"""EXP-20260930-001 の細胞ごとの表を作って parquet に貯める（事後診断・事前登録の外）。

**摂動は 1 回も回さない。** 既存の埋め込みを読み直すだけ。
事前登録済みの `analysis/evaluate_004.py` は **import するだけで一切書き換えない**。

出力: posthoc/cells_cache.parquet
"""
from __future__ import annotations

import importlib.util
import sys
import time
from pathlib import Path

import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
E = HERE.parent
ROOT = E.parents[1]
sys.path.insert(0, str(ROOT / "scripts/analysis"))


def load_evaluate_module():
    p = E / "analysis/evaluate_004.py"
    spec = importlib.util.spec_from_file_location("evaluate_004", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main() -> None:
    ev = load_evaluate_module()
    cfg = yaml.safe_load(open(E / "config.yaml"))
    t0 = time.time()
    c, drop = ev.load(cfg)
    print(f"load 完了 {time.time() - t0:.1f} 秒, {len(c):,} 行, 移動量0除外 {drop['n_zero_movement']}")
    out = HERE / "cells_cache.parquet"
    c.to_parquet(out, index=False)
    print(f"保存: {out}  列: {list(c.columns)}")


if __name__ == "__main__":
    main()

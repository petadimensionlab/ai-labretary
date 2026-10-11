#!/usr/bin/env python3
"""比べられる対照（eligibility）を seed ごとに 1 回だけ計算して貯める（事後診断）。

元の位置 lp と入力の長さ ll だけで決まり、**補正の版によらない**ので使い回せる。
事前登録済みの `analysis/evaluate_004.py` の `eligibility()` をそのまま import して使う。
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
E = HERE.parent
ROOT = E.parents[1]
sys.path.insert(0, str(ROOT / "scripts/analysis"))

spec = importlib.util.spec_from_file_location("evaluate_004", E / "analysis/evaluate_004.py")
ev = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ev)

cfg = yaml.safe_load(open(E / "config.yaml"))
pm = cfg["analysis"]["pos_matched_permutation"]
ensmap = {g["symbol"]: g["ensembl"] for g in json.loads((E / "work/genes_004.json").read_text())["genes"]}
nulls_all = set(ensmap[x["symbol"]] for x in
                json.loads((E / "work/null_pool_004.json").read_text())["versions"]["B_500"])

c = pd.read_parquet(HERE / "cells_cache.parquet")
for s in sorted(c.seed.unique()):
    out = HERE / f"eligibility-seed{s}.json"
    if out.exists():
        print(f"既にある: {out}")
        continue
    t0 = time.time()
    cs = c[c.seed == s]
    nl = {h: hd for h, hd in cs[cs.gene.isin(nulls_all)].groupby("gene")}
    el = ev.eligibility(cs, nl, pm["tolerance_log2"], pm["min_cover_frac"])
    out.write_text(json.dumps(el))
    print(f"seed{s}: {time.time() - t0:.1f} 秒, {len(el)} 遺伝子")

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
EXP-20261011-001 の緩和検算（**追記の感度分析。事前登録の判定を変えない**）.
PI 指示（2026-10-11）: 主判定は「判定不能（R-6）」のまま。
そのうえで「(ii) 退化30 の有意水準を緩めたら P-1 に合格するか」を検算し、合格する最小の有意水準を
「条件付き合格」として記録する. P-1 (i)(iii) と偽セット生成は evaluate_014.py と同一の式・同一の
乱数の鍵を使う。有意水準 q だけを変えて、最初に P-1 が True になる q を探す。
出力: results/relax_014.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
for p in ("experiments/EXP-20261004-001/analysis", "experiments/EXP-20261001-001/analysis",
          "scripts/analysis", "experiments/EXP-20261011-001/analysis"):
    sys.path.insert(0, str(ROOT / p))
import importlib.util

spec = importlib.util.spec_from_file_location(
    "ev014", ROOT / "experiments/EXP-20261011-001/analysis/evaluate_014.py")
ev = importlib.util.module_from_spec(spec)
if spec.loader is None:
    raise SystemExit("loader is None")
spec.loader.exec_module(ev)
from evaluate_005 import set_null  # noqa: E402

EXP = ev.EXP
SEEDS = ev.SEEDS
N_DRAW = ev.N_DRAW
PLAN = ev.PLAN
Key = PLAN["乱数の鍵"]

cfg = yaml.safe_load((EXP / "config.yaml").read_text())
G = ev.gene_sets_from_cfg_and_json(cfg)

rank = pd.read_csv(EXP / "results" / "gene_rank_014.csv").set_index("gene")
Zdf = pd.DataFrame(rank[["z_81", "z_82", "z_83"]].to_numpy(), index=rank.index, columns=SEEDS)

pool = [g for g in G["sets"]["プールA1"] if g in Zdf.index]
poolZ = Zdf.reindex(pool).to_numpy()
pc71 = [g for g in G["sets"]["活性化PC71"] if g in Zdf.index]
pc30 = [g for g in G["sets"]["退化PC30"] if g in Zdf.index]
ov_ens = {ens for ens, s in G["symbol"].items() if s in PLAN["重なり7"]}
pc71r = [g for g in pc71 if g not in ov_ens]
pc30r = [g for g in pc30 if g not in ov_ens]

obs71 = Zdf.reindex(pc71).to_numpy().mean(axis=0)
obs30 = Zdf.reindex(pc30).to_numpy().mean(axis=0)
obs71r = Zdf.reindex(pc71r).to_numpy().mean(axis=0)
obs30r = Zdf.reindex(pc30r).to_numpy().mean(axis=0)

null71 = set_null(poolZ, len(pc71), N_DRAW, np.random.default_rng(Key["P-1 活性化の偽セット"]), replace=True)
null30 = set_null(poolZ, len(pc30), N_DRAW, np.random.default_rng(Key["P-1 退化の偽セット"]), replace=True)
thr71 = np.quantile(null71, 0.99, axis=0)

null_diff = np.empty((3, N_DRAW))
for s in range(3):
    ia = np.random.default_rng(Key["P-1 (iii) の無作為 2 セット"] + s).integers(
        0, len(poolZ), size=(N_DRAW, len(pc71)))
    ib = np.random.default_rng(Key["P-1 (iii) の無作為 2 セット"] + 1000 + s).integers(
        0, len(poolZ), size=(N_DRAW, len(pc30)))
    null_diff[s] = poolZ[:, s][ia].mean(axis=1) - poolZ[:, s][ib].mean(axis=1)
thr_diff = np.quantile(null_diff, 0.99, axis=1)

c_i = bool(np.all(obs71 > 0) and np.all(obs71 > thr71))
c_iii = bool(np.all((obs71 - obs30) > thr_diff))

thr71r = np.quantile(set_null(poolZ, len(pc71r), N_DRAW,
                              np.random.default_rng(Key["P-1 活性化の偽セット"] + 2), replace=True),
                     0.99, axis=0)
c_ir = bool(np.all(obs71r > 0) and np.all(obs71r > thr71r))
c_iiir = bool(np.all((obs71r - obs30r) > thr_diff))

q_list = [0.01, 0.02, 0.025, 0.03, 0.04, 0.05, 0.075, 0.10]


def c_ii(q, obs, n, rng):
    thr = np.quantile(set_null(poolZ, n, N_DRAW, rng, replace=True), q, axis=0)
    return bool(np.all(obs < 0) and np.all(obs < thr)), thr


levels = []
first_pass = None
for q in q_list:
    c_ii_q, thr30_q = c_ii(q, obs30, len(pc30), np.random.default_rng(Key["P-1 退化の偽セット"]))
    p1_q = c_i and c_ii_q and c_iii
    if p1_q and first_pass is None:
        first_pass = q
    levels.append({"有意水準_q": q, "退化30_観測(3seed)": [round(x, 4) for x in obs30],
                   "退化30_偽セット点(3seed)": [round(x, 4) for x in thr30_q],
                   "(ii) 成り立つ": c_ii_q, "P-1": p1_q})

levels_r = []
first_pass_r = None
for q in q_list:
    c_iir_q, thr30r_q = c_ii(q, obs30r, len(pc30r), np.random.default_rng(Key["P-1 退化の偽セット"] + 2))
    p1r_q = c_ir and c_iir_q and c_iiir
    if p1r_q and first_pass_r is None:
        first_pass_r = q
    levels_r.append({"有意水準_q": q, "退化23_観測(3seed)": [round(x, 4) for x in obs30r],
                     "退化23_偽セット点(3seed)": [round(x, 4) for x in thr30r_q],
                     "(ii) 成り立つ": c_iir_q, "P-1": p1r_q})


def attained(obs):
    out = []
    for s in range(3):
        frac = float(np.mean(np.sort(null30[:, s]) < obs[s]))
        out.append(round(frac, 4))
    return out


out = {
    "実験": "EXP-20261011-001",
    "仮説": "HYP-20261011-01",
    "目的": "PI 指示（2026-10-11）の緩和検算。(ii) 退化30 の有意水準 q を緩めたとき、P-1 が合格する最小の q を「条件付き合格」として記録する。本登録の判定（判定不能 R-6）は変えない。",
    "先行": "本登録（results/criteria_014.json）の P-1=False： (i)=True・(ii)=False（1% 未満に達せず）・(iii)=True、最初に当たった反証条件 R-6、R-12 交絡 該当",
    "使ったもの": {
        "z": "results/gene_rank_014.csv（既存 626＋新規 30 の和集合の u_diff 標準化 z。重なり 10 個は A-4 で新規と max 差 0.0 のため既存側で可）",
        "プールA1": f"{len(pool)} 個（残欠 {len(G['sets']['プールA1']) - len(pool)}）",
        "活性化71": len(pc71), "退化30": len(pc30), "退化23（重なり7を除いた版）": len(pc30r),
        "偽セットの式・乱数の鍵": "evaluate_014.py と同一（P-1 (iii) の無作為 2 セット鍵ほか）"},
    "方法": "P-1 (i)(iii) と偽セット生成は本登録と同一・同一鍵。有意水準 q だけを 1%→10% に動かし、最初に P-1 が全て True になる q を探す",
    "本登録相当（q=0.01）": {"(i)": c_i, "(ii) 成り立つ": c_ii(0.01, obs30, len(pc30),
                                                 np.random.default_rng(Key["P-1 退化の偽セット"]))[0],
                          "(iii)": c_iii, "P-1": c_i and c_ii(0.01, obs30, len(pc30),
                                                 np.random.default_rng(Key["P-1 退化の偽セット"]))[0] and c_iii},
    "観測（退化30 z_diff 3seed）": [round(x, 4) for x in obs30],
    "観測（退化23 z_diff 3seed）": [round(x, 4) for x in obs30r],
    "有意水準を動かしたとき（退化30）": levels,
    "条件付き合格になる最小の有意水準（退化30）": first_pass,
    "退化30 の 3 seed 別に観測値が達した偽セット分位": attained(obs30),
    "有意水準を動かしたとき（退化23・重なり除去の版）": levels_r,
    "条件付き合格になる最小の有意水準（退化23）": first_pass_r,
    "limitations": "この結果は後の確認で採録した探索的検算。有意水準を緩めた分だけ I 型過誤の面で弱い。本登録（判定不能 R-6）と並記し、判断は求められた側（PI）が行う。",
}
with open(EXP / "results" / "relax_014.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("written:", EXP / "results" / "relax_014.json")
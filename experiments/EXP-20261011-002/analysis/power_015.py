#!/usr/bin/env python3
"""EXP-20261011-002 の検出力シミュレーション（結果を見る前に実行。PC の観測値は使わない）。

新規則（マッチ偽セット・主判定 P-1'）そのものをシミュレートし、最小検出可能効果量（MDE）を出す。
- 帰無分布: 対照プール 549 の公開済み標準得点（EXP-20261011-001/results/gene_rank_014.csv のプール行）と
  注釈テーブル（precheck_015 と同じ定義）から、マッチ偽セット 2,000 組を引く
- 検出力: 真の効果 μ（3 群同じ大きさ）を加えたセット平均 = マッチ偽セット 1 組 + μ とし、
  3 群すべてで規則を満たす割合を 200 回で見る
- 出力: (ii) 退化 23 の MDE・(i) 活性化 71 の MDE・(iii) の MDE・同時の MDE
  （マッチ規則と、参考の無作為規則の両方）
**PC のセット平均（観測）は一度も計算しない。**

走らせ方:
    .venv/bin/python experiments/EXP-20261011-002/analysis/power_015.py
出力: work/power_015.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
EXP = ROOT / "experiments/EXP-20261011-002"
A1_NEW = ROOT / "experiments/EXP-20261011-001"
sys.path.insert(0, str(EXP / "analysis"))
sys.path.insert(0, str(ROOT / "experiments/EXP-20261001-001/analysis"))
import evaluate_015 as ev015            # noqa: E402
from evaluate_005 import set_null      # noqa: E402
from precheck_015 import load_a1_annotations   # noqa: E402

SEEDS = [81, 82, 83]
N_NULL = 2000
N_REP = 200
POWER_TARGET = 0.80


def main() -> None:
    cfg = yaml.safe_load((EXP / "config.yaml").read_text())
    G = ev015.gene_sets(cfg)
    sets = G["sets"]
    ann = load_a1_annotations(cfg)
    rank = pd.read_csv(A1_NEW / "results/gene_rank_014.csv").set_index("gene")
    Zdf = rank[["z_81", "z_82", "z_83"]].rename(columns={"z_81": "81", "z_82": "82", "z_83": "83"})

    pool = [g for g in sets["プールA1"] if g in Zdf.index and g in ann.index]
    poolZ = Zdf.reindex(pool).to_numpy()

    pc71 = [g for g in sets["活性化PC71"] if g in ann.index and np.isfinite(ann.loc[g, "lp_med"])]
    pc23 = [g for g in sets["退化PC23"] if g in ann.index and np.isfinite(ann.loc[g, "lp_med"])]
    ov23 = {e for e in pc23 if e in set(pool)}
    ov71 = {e for e in pc71 if e in set(pool)}
    draw23 = [g for g in pool if g not in ov23]
    draw71 = [g for g in pool if g not in ov71]

    def _aligned(draw_pool):
        idx = [g for g in draw_pool
               if g in ann.index and g in Zdf.index
               and ann.loc[g].notna().all() and Zdf.loc[g].notna().all()]
        return ann.reindex(idx), Zdf.reindex(idx)

    ann_pool23, Zpool23_df = _aligned(draw23)
    ann_pool71, Zpool71_df = _aligned(draw71)
    Zpool23 = Zpool23_df.to_numpy()
    Zpool71 = Zpool71_df.to_numpy()
    ann_pc23 = ann.reindex(pc23).dropna()
    ann_pc71 = ann.reindex(pc71).dropna()

    # ---- 帰無分布（マッチ偽セット）
    rng23 = np.random.default_rng(ev015.PLAN["乱数の鍵"]["P-1 退化のマッチ偽セット"])
    rng71 = np.random.default_rng(ev015.PLAN["乱数の鍵"]["P-1 活性化のマッチ偽セット"])
    null23 = ev015.matched_sham_sets(Zpool23, ann_pool23, ann_pc23, N_NULL, rng23)
    null71 = ev015.matched_sham_sets(Zpool71, ann_pool71, ann_pc71, N_NULL, rng71)
    thr23 = np.quantile(null23, 0.01, axis=0)
    thr71 = np.quantile(null71, 0.99, axis=0)
    # 無作為（参考）
    rng_u = np.random.default_rng(ev015.PLAN["乱数の鍵"]["補正なし（無作為）の偽セット"])
    rnull23 = set_null(poolZ, len(pc23), N_NULL, rng_u, replace=True)
    rnull71 = set_null(poolZ, len(pc71), N_NULL, rng_u, replace=True)
    rthr23 = np.quantile(rnull23, 0.01, axis=0)
    rthr71 = np.quantile(rnull71, 0.99, axis=0)

    def power_ii(mu: float, null: np.ndarray, thr: np.ndarray) -> float:
        obs = null + mu
        return float((np.all(obs < 0, axis=1) & np.all(obs < thr, axis=1)).mean())

    def power_i(mu: float, null: np.ndarray, thr: np.ndarray) -> float:
        obs = null + mu
        return float((np.all(obs > 0, axis=1) & np.all(obs > thr, axis=1)).mean())

    def mde_grid(null: np.ndarray, thr: np.ndarray, side: str) -> float | None:
        mus = (np.arange(0.05, 2.01, 0.05) * (-1 if side == "deg" else 1))
        for mu in mus:
            p = power_ii(-abs(mu), null, thr) if side == "deg" else power_i(abs(mu), null, thr)
            if p >= POWER_TARGET:
                return round(abs(float(mu)), 3)
        return None

    # 同時 (i)∧(ii)∧(iii): 活性化 μ71=+m、退化 μ23=−m。差の帰無は null71−null23 の分布
    null_diff = null71[:N_REP] - null23[:N_REP]
    thr_diff = np.quantile(null71 - null23, 0.99, axis=0)

    def power_simul(m: float) -> float:
        n = N_REP
        obs71 = null71[:n] + m
        obs23 = null23[:n] - m
        ok_i = np.all(obs71 > 0, axis=1) & np.all(obs71 > thr71, axis=1)
        ok_ii = np.all(obs23 < 0, axis=1) & np.all(obs23 < thr23, axis=1)
        d = obs71 - obs23
        ok_iii = np.all(d > thr_diff, axis=1)
        return float((ok_i & ok_ii & ok_iii).mean())

    mus = np.arange(0.05, 2.01, 0.05)
    sim_powers = [power_simul(float(m)) for m in mus]
    mde_simul = next((round(float(m), 3) for m, p in zip(mus, sim_powers) if p >= POWER_TARGET), None)

    out = {
        "実験": "EXP-20261011-002", "仮説": "HYP-20261011-01",
        "対象": "主判定 P-1'（A1。マッチ偽セット規則）。PC の観測値（セット平均）は使っていない",
        "方法": "帰無分布はマッチ偽セット 2,000 組（プール 549 の公開済み標準得点＋注釈）。検出力は 1 組の帰無に効果 μ を足した観測で 3 群すべて規則を満たす割合（200 回）",
        "帰無分布の閾値（マッチ）": {
            "退化23 の 1% 点": [round(x, 4) for x in thr23],
            "活性化71 の 99% 点": [round(x, 4) for x in thr71]},
        "帰無分布の閾値（無作為・参考）": {
            "退化23 の 1% 点": [round(x, 4) for x in rthr23],
            "活性化71 の 99% 点": [round(x, 4) for x in rthr71]},
        "MDE（効果の大きさ。z 単位。大きいほど検出が難しい）": {
            "退化23（マッチ規則。条件 (ii)）": mde_grid(null23, thr23, "deg"),
            "活性化71（マッチ規則。条件 (i)）": mde_grid(null71, thr71, "ada"),
            "退化23（無作為規則・参考）": mde_grid(rnull23, rthr23, "deg"),
            "同時（マッチ規則。活性化 +m と退化 −m を同時。200 回）": mde_simul},
        "同時の検出力カーブ（m=0.05..2.00。マッチ規則）": [
            {"m": round(float(m), 2), "power": round(p, 3)} for m, p in zip(mus, sim_powers)],
        "注": ("この値は設計時のシミュレーションであり、判定規則そのものの MDE である。"
               "退化 23 は K が 30 から減るため、前回（K=30・無作為規則で ≈0.66）よりやや悪化する見込み。"
               "粗い式 4.5/√K は報告に使わない（guardrails 2026-09-30）。"),
    }
    (EXP / "work").mkdir(parents=True, exist_ok=True)
    (EXP / "work" / "power_015.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    print("written:", EXP / "work" / "power_015.json")


if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""EXP-20260929-001: 検出下限を、**事前登録する規則そのもの**で計算する（guardrails 2026-09-28）。
併せて、主判定の絶対値条件 `set_abs_min` を対照遺伝子の実測から較正する。

較正の条件（**これを必ず添えて報告する**）:
  規則 = 「3 seed すべてで (i) 対照遺伝子から**復元ありで**引いた K 個セットの 99 パーセンタイルを超え、
          かつ (ii) 絶対値 set_abs_min 以上」。偶然合格率 5% に較正した規則の値ではない。

使うデータ: **EXP-20260928-002 の対照遺伝子 980 個 × 3 seed の実測値**（目標の軸、対照 980 版で標準化した z）。
本実験の摂動は 1 回も実行していない。EXP-20260928-002 は出発点も目標も違うので、これは
**この系で得られる最良の事前推定**であり、本番後に同じ規則で実測し直す（plan.md §7）。
遺伝子ごとの効果のばらつき τ は EXP-20260928-001 の目標マーカー 24 個の遺伝子間標準偏差 2.051 を使う
（guardrails 2026-09-27: 合成データではなく実測のばらつきを使う）。

出力: work/power_003.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20260929-001"
SRC = ROOT / "experiments/EXP-20260928-002"

N_DRAW = 200_000
PCT = 0.99
TAU = 2.051                 # 遺伝子ごとの効果のばらつき（EXP-20260928-001 の目標マーカー 24 個の実測。保守的）
TAU_ALT = [0.0, 0.45]       # 0.0 = 効果が均一 / 0.45 = EXP-20260928-002 の経路セットの実測の超過ばらつき
ABS_GRID = np.round(np.arange(0.05, 1.21, 0.005), 3)
MAX_CHANCE = 0.02
# **事前登録する絶対値条件**: HYP-20260929-01 §7 が K = 25 で見込んだ 99 パーセンタイル。
# 結果を見る前に固定する。下の較正で「この値でも偶然合格率 0.02 以下」であることを確かめる。
SET_ABS_MIN = 0.46
POWER_TARGET = 0.80
SEED = 20260930


def znull() -> np.ndarray:
    """EXP-20260928-002 の対照遺伝子 980 個 × 3 seed の z（目標の軸）。"""
    g = pd.read_csv(SRC / "results/eval/gene_seed_002.csv")
    g = g[(g.col == "m_goal") & (g.role == "negative")]
    seeds = sorted(g.seed.unique())
    p = g.pivot_table(index="gene", columns="seed", values="D_raw")
    z = np.column_stack([(p[s] - p[s].mean()) / p[s].std(ddof=1) for s in seeds])
    return z[~np.isnan(z).any(axis=1)]


def main() -> None:
    z = znull()
    rng = np.random.default_rng(SEED)
    N, S = z.shape
    res = {"source": "EXP-20260928-002 の対照遺伝子（目標の軸, 980 版）",
           "n_null_genes": int(N), "n_seeds": int(S), "n_draw": N_DRAW,
           "sigma_between_genes": [round(float(z[:, s].std(ddof=1)), 4) for s in range(S)],
           "seed_correlation": round(float(np.corrcoef(z.T)[0, 1]), 4),
           "tau_between_genes_effect": TAU, "percentile": PCT,
           "calibration_note": "閾値は『対照から復元ありで引いた同じ大きさのセットの 99 パーセンタイル』。"
                               "偶然合格率 5% に較正した規則ではない（guardrails 2026-09-28）",
           "sets": {}}

    def draws(k: int) -> np.ndarray:
        idx = rng.integers(0, N, size=(N_DRAW, k))       # 復元あり。seed をまたいで同じ遺伝子を引く
        return z[idx].mean(axis=1)                        # (N_DRAW, S)

    def power(k: int, delta: float, thr: np.ndarray, abs_min: float, tau: float,
              n: int = 40_000) -> float:
        """真の効果 delta（遺伝子ごとに SD tau でばらつく）を対照の分布に足したときの合格率。"""
        idx = rng.integers(0, N, size=(n, k))
        eff = rng.normal(delta, tau, size=(n, k)) if tau > 0 else np.full((n, k), float(delta))
        obs = z[idx].mean(axis=1) + eff.mean(axis=1)[:, None]
        return float(((obs > thr[None, :]).all(axis=1) & (obs >= abs_min).all(axis=1)).mean())

    def floor_at(k: int, thr: np.ndarray, abs_min: float, tau: float) -> float:
        """事前登録する規則そのもので、検出力 80% になる効果の大きさ（二分法）。"""
        lo, hi = 0.0, 5.0
        for _ in range(28):
            mid = (lo + hi) / 2
            if power(k, mid, thr, abs_min, tau) >= POWER_TARGET:
                hi = mid
            else:
                lo = mid
        return (lo + hi) / 2

    ks = {"secretory_all": 24, "co_or_low": 19, "co_or_high": 5, "tls_chemokine": 4,
          "org_markers7": 7, "myof_markers24": 24, "housekeeping3": 3, "negctrl_panel7": 7}
    # --- K = 24（主判定）で set_abs_min を較正する
    d24 = draws(ks["secretory_all"])
    thr24 = np.quantile(d24, PCT, axis=0)
    pass_pct24 = (d24 > thr24[None, :]).all(axis=1)
    chance = {float(a): float((pass_pct24 & (d24 >= a).all(axis=1)).mean()) for a in ABS_GRID}
    res["set_abs_min"] = SET_ABS_MIN
    res["set_abs_min_source"] = ("HYP-20260929-01 §7 が K=25 で見込んだ 99 パーセンタイル 0.46。"
                                 "結果を見る前に固定した")
    res["chance_rate_at_abs_min"] = round(chance[SET_ABS_MIN], 5)
    res["chance_rate_percentile_only"] = round(float(pass_pct24.mean()), 5)
    res["abs_min_binding_in_prior_data"] = bool(SET_ABS_MIN > float(thr24.max()))
    res["smallest_abs_meeting_max_chance"] = float(
        min(a for a in ABS_GRID if chance[float(a)] <= MAX_CHANCE))
    res["abs_grid_examples"] = {str(a): round(chance[a], 5)
                                for a in [0.10, 0.20, 0.30, 0.40, 0.46, 0.50, 0.60, 0.80]}
    set_abs_min = SET_ABS_MIN

    for name, k in ks.items():
        dd = draws(k)
        thr = np.quantile(dd, PCT, axis=0)
        eff_thr = max(float(thr.max()), set_abs_min)
        cr = float(((dd > thr[None, :]).all(axis=1) & (dd >= set_abs_min).all(axis=1)).mean())
        floor_tau = floor_at(k, thr, set_abs_min, TAU)
        floors_alt = {str(t): round(floor_at(k, thr, set_abs_min, t), 3) for t in TAU_ALT}
        # 復元なしの閾値（対照 300 / 500 での差。guardrails 2026-09-28）
        nr = {}
        for npool in (300, 500):
            if k <= npool:
                sub = z[rng.integers(0, N, size=npool)]   # 300 / 500 個の有限母集団を模す
                i2 = np.argsort(rng.random((50_000, npool)), axis=1)[:, :k]
                nr[str(npool)] = round(float(np.quantile(sub[i2].mean(axis=1), PCT, axis=0).max()), 4)
        res["sets"][name] = {
            "K": k, "percentile_99": [round(float(t), 4) for t in thr],
            "effective_threshold": round(eff_thr, 4),
            "chance_pass_rate": round(cr, 5),
            "detection_floor_uniform_effect": floors_alt["0.0"],
            "detection_floor_by_tau": floors_alt,
            "detection_floor_with_gene_sd": round(floor_tau, 3),
            "percentile_99_without_replacement": nr}

    (E / "work/power_003.json").write_text(json.dumps(res, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "sets"}, indent=1, ensure_ascii=False))
    print(f"\n{'セット':16s} {'K':>3s} {'99%点':>8s} {'実効閾値':>9s} {'偶然合格率':>10s} "
          f"{'下限(均一)':>10s} {'下限(ばらつき込み)':>16s}")
    for k, v in res["sets"].items():
        print(f"{k:16s} {v['K']:3d} {max(v['percentile_99']):8.3f} {v['effective_threshold']:9.3f} "
              f"{v['chance_pass_rate']:10.5f} {v['detection_floor_uniform_effect']:10.3f} "
              f"{v['detection_floor_with_gene_sd']:16.3f}")


if __name__ == "__main__":
    main()

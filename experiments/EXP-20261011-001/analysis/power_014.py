#!/usr/bin/env python3
"""EXP-20261011-001 — 主判定 P-1 の最小検出可能効果量（MDE）を、規則そのものをシミュレートして事前に出す。

HYP-20261011-01 §5.2 は「最小検出可能効果量は、この判定規則そのものをシミュレートして事前に出す」と
書いており、粗い式（4.5/√K・6.63/√K）は計画段階の目安にしか使わない（guardrails 2026-09-30）。

判定規則（P-1。A1 の 3 群すべて）:
  (i)   活性化の陽性対照 71 個のセット平均 z_diff > 偽セット（71 個×1,000 組）の 99 パーセンタイル
  (ii)  退化の陽性対照 30 個のセット平均 z_diff < 偽セット（30 個×1,000 組）の 1 パーセンタイル
  (iii) (活性化 71 のセット平均 − 退化 30 のセット平均) が、大きさ 71 と 30 の無作為な 2 セットの差の
        99 パーセンタイル超

シミュレーションの仮定（結果を見る前に固定する）:
- 対照プールの z は、標準化済みなので平均 0・分散 1 の正規に近いと置き、**seed 間の相関 ρ** だけを
  既に公開されている A2 の対照プール（510 個）の z（活性化の向き `gene_group_011.csv`、退化の向き
  `gene_group_013.csv`）から実測して使う。z_diff の null 分布そのものは A2/A1 の実測がないため、
  正規モデルで模す（限界。実際の較正は評価時に A1 の対照プール 549 で実測し、判定には使わず記録する）。
- 偽セットは、対照プールから**復元あり**で引く（`evaluate_005.set_null` と同じ手続きを模す）。
- 群（seed）どうしは相関 ρ の多変量正規（同じ遺伝子の seed 間の値は相関する。guardrails 2026-09-27）。

出力: experiments/EXP-20261011-001/work/power_014.json（事前登録のハッシュに含める）
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
EXP = ROOT / "experiments/EXP-20261011-001"
OUT = EXP / "work/power_014.json"

N_POOL = 549            # A1 の対照プール
N_DRAW = 1000           # 偽セットの組数（評価と同じ）
N_REP = 400             # 検出力の繰り返し（各 μ について）
SEEDS = [81, 82, 83]
K71, K30 = 71, 30


def rho_from_csv(path: Path, pool: list[str], direction: str | None = None) -> float:
    d = pd.read_csv(path)
    if direction is not None:
        d = d[d["方向"] == direction]
    d = d[d.gene.isin(pool)]
    z = d.pivot(index="gene", columns="seed", values="z").reindex(columns=SEEDS)
    z = z.dropna()
    c = z.corr().to_numpy()
    off = c[np.triu_indices(3, 1)]
    return float(np.mean(off))


def per_seed_null_quantiles(rng, rho: float, k: int, q: float) -> np.ndarray:
    """対照プール N 個の z（平均 0・分散 1・seed 間相関 ρ）を作り、復元ありで N_DRAW 組（k 個）を引き、
    各 seed のセット平均の分位数 q を返す（3 群）。"""
    C = np.full((3, 3), rho); np.fill_diagonal(C, 1.0)
    zz = rng.multivariate_normal(np.zeros(3), C, size=N_POOL).T          # (3, N_POOL)
    means = np.empty((3, N_DRAW))
    for s in range(3):
        idx = rng.integers(0, N_POOL, size=(N_DRAW, k))
        means[s] = zz[s][idx].mean(axis=1)
    return np.quantile(means, q, axis=1)


def power_set(rng, rho: float, mu: float, k: int, thr: np.ndarray, direction: int) -> float:
    """方向 direction=+1（活性化）: all3(mean > thr)；-1（退化）: all3(mean < thr)。"""
    C = np.full((3, 3), rho); np.fill_diagonal(C, 1.0)
    hit = 0
    for _ in range(N_REP):
        zz = rng.multivariate_normal(np.full(3, mu), C, size=k).T        # (3, k)
        m = zz.mean(axis=1)
        ok = bool(np.all((m > thr) if direction > 0 else (m < thr)))
        hit += ok
    return hit / N_REP


def power_diff(rng, rho: float, mu: float, thr: np.ndarray) -> float:
    """(iii): 活性化 71 を +μ、退化 30 を −μ で動かしたとき、all3(diff > thr)。"""
    C = np.full((3, 3), rho); np.fill_diagonal(C, 1.0)
    hit = 0
    for _ in range(N_REP):
        za = rng.multivariate_normal(np.full(3, mu), C, size=K71).T.mean(axis=1)
        zd = rng.multivariate_normal(np.full(3, -mu), C, size=K30).T.mean(axis=1)
        hit += bool(np.all(za - zd > thr))
    return hit / N_REP


def power_joint(rng, rho: float, mu: float, t71: np.ndarray, t30: np.ndarray, tdiff: np.ndarray) -> float:
    C = np.full((3, 3), rho); np.fill_diagonal(C, 1.0)
    hit = 0
    for _ in range(N_REP):
        za = rng.multivariate_normal(np.full(3, mu), C, size=K71).T.mean(axis=1)
        zd = rng.multivariate_normal(np.full(3, -mu), C, size=K30).T.mean(axis=1)
        hit += bool(np.all(za > t71) and np.all(zd < t30) and np.all(za - zd > tdiff))
    return hit / N_REP


def mde_for(rng, rho: float, func, lo: float, hi: float) -> float | None:
    """単調（lo 側は低パワー、hi 側は高パワー）と仮定して二分探索で power=0.8 を探す。"""
    if func(lo) >= 0.8:
        return lo
    if func(hi) < 0.8:
        return None
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if func(mid) >= 0.8:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def main() -> None:
    # 既に公開されている A2 の対照プールの z から、seed 間相関 ρ を実測する（記録・検出力の入力）
    g11 = json.loads((ROOT / "experiments/EXP-20261007-001/work/genes_011.json").read_text())
    pool510 = [r["ensembl"] for r in g11["対照500"]] + list(g11["恒常発現3"].values()) + list(g11["中立7"].values())
    rho_ada = rho_from_csv(ROOT / "experiments/EXP-20261007-001/results/gene_group_011.csv", pool510)
    rho_deg = rho_from_csv(ROOT / "experiments/EXP-20261007-003/results/gene_group_013.csv", pool510,
                           direction="degenerate")

    rough71_45 = 4.5 / np.sqrt(K71)
    rough30_45 = 4.5 / np.sqrt(K30)
    rough71_663 = 6.63 / np.sqrt(K71)
    rough30_663 = 6.63 / np.sqrt(K30)

    out: dict = {
        "実験": "EXP-20261011-001", "仮説": "HYP-20261011-01", "対象": "主判定 P-1（A1）",
        "注": ("判定規則そのものをシミュレートした最小検出可能効果量（検出力 80%）。"
               "仮定: 対照プールの z は平均 0・分散 1 の正規で seed 間相関 ρ（A2 の対照プールの実測から）。"
               "z_diff の null 分布そのものは実測できないため正規モデルで模す。"
               "実際の較正（偽セット合格率）は評価時に A1 のプール 549 で実測し、判定には使わず記録する"),
        "seed 間相関 ρ の実測（公開済み A2 の対照プールの z）": {
            "活性化の向き（gene_group_011.csv）": round(rho_ada, 4),
            "退化の向き（gene_group_013.csv）": round(rho_deg, 4),
            "使う値（感度分析で ρ=0・実測・0.9 を並べる）": round(rho_ada, 4),
        },
        "粗い式（**目安。判定に使わない**）": {
            "4.5/√K（仮説 §5.2 に併記）": {"活性化71": round(rough71_45, 3), "退化30": round(rough30_45, 3)},
            "6.63/√K（仮説 §10 に併記）": {"活性化71": round(rough71_663, 3), "退化30": round(rough30_663, 3)},
        },
        "シミュレーション": {"対照プールの大きさ": N_POOL, "偽セットの組数": N_DRAW,
                              "検出力の繰り返し": N_REP, "群": SEEDS},
        "結果（MDE。z の単位。大きいほど検出が難しい）": {},
    }

    rng0 = np.random.default_rng(20261023)
    for rho, label in [(0.0, "ρ=0（独立 seed）"), (rho_ada, f"ρ=実測 {rho_ada:.3f}"), (0.9, "ρ=0.9（強い相関）")]:
        rng = np.random.default_rng(20261023 + int(round(rho * 1000)))
        t71 = per_seed_null_quantiles(rng, rho, K71, 0.99)
        t30 = per_seed_null_quantiles(rng, rho, K30, 0.01)
        tdiff = per_seed_null_quantiles(rng, rho, K71 + K30, 0.99)   # 差の尺度は同じ（和集合で近似しない）。下で再計算
        # (iii) の null は「71 個の平均 − 30 個の平均」の分布。独立に引いて差を取る
        C = np.full((3, 3), rho); np.fill_diagonal(C, 1.0)
        zz = rng.multivariate_normal(np.zeros(3), C, size=N_POOL).T
        diffs = np.empty((3, N_DRAW))
        for s in range(3):
            ia = rng.integers(0, N_POOL, size=(N_DRAW, K71))
            ib = rng.integers(0, N_POOL, size=(N_DRAW, K30))
            diffs[s] = (zz[s][ia].mean(axis=1) - zz[s][ib].mean(axis=1))
        tdiff = np.quantile(diffs, 0.99, axis=1)

        mde_i = mde_for(rng, rho, lambda m: power_set(rng, rho, m, K71, t71, +1), 0.05, 2.5)
        mde_ii = mde_for(rng, rho, lambda m: power_set(rng, rho, -m, K30, t30, -1), 0.05, 2.5)
        mde_iii = mde_for(rng, rho, lambda m: power_diff(rng, rho, m, tdiff), 0.05, 2.5)
        mde_joint = mde_for(rng, rho, lambda m: power_joint(rng, rho, m, t71, t30, tdiff), 0.05, 3.0)
        out["結果（MDE。z の単位。大きいほど検出が難しい）"][label] = {
            "MDE (i) 活性化 71 が正へ": None if mde_i is None else round(mde_i, 3),
            "MDE (ii) 退化 30 が負へ": None if mde_ii is None else round(mde_ii, 3),
            "MDE (iii) 2 セットの差（±μ 対称）": None if mde_iii is None else round(mde_iii, 3),
            "MDE 同時（(i)∧(ii)∧(iii) を 80%）": None if mde_joint is None else round(mde_joint, 3),
            "読み方": "同時の MDE が、この判定規則が 80% の確率で通すために必要なセット平均のずれ（z 単位）",
        }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
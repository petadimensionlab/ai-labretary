#!/usr/bin/env python3
"""EXP-20261001-001 — 検出下限を、**事前登録する規則そのもの**で計算する（guardrails 2026-09-28）。

**選定手順 §6 のとおり、EXP-20260930-001 の power_004.py と同じ規則・同じ対照分布・同じ τ = 2.051 を使う。**
並べ替えの seed も 20260930 のまま変えない（選定手順の出力 counts.json の検出下限を再現するため）。
再現しなければ停止する。

較正の条件（**これを必ず添えて報告する**）:
  規則 = 「**3 seed すべてで** (i) 対照遺伝子から**復元ありで**引いた K 個セットの
          99 パーセンタイルを超え、かつ (ii) 絶対値が `analysis.set_abs_min[セット]` 以上」。
  **偶然合格率 5% に較正した規則ではない。**

使うデータ: **EXP-20260929-001 の対照遺伝子 500 個 × 3 seed の実測値**（目標の軸, B_500 版で標準化）。
本実験の摂動は 1 回も実行していない。EXP-20260929-001 は**目標の定義も交絡の補正の版も違う**ので、
これは「この系で得られる最良の事前推定」であり、本番後に本実験の対照 500 個で実測し直す。
遺伝子ごとの効果のばらつき τ は EXP-20260928-001 の実測 2.051 を使う（保守的。guardrails 2026-09-27）。

出力: work/power_005.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20261001-001"
SRC = ROOT / "experiments/EXP-20260929-001"

N_DRAW = 200_000
PCT = 0.99
TAU = 2.051
TAU_ALT = [0.0, 0.45, 2.051]
POWER_TARGET = 0.80
SEED = 20260930


def znull() -> np.ndarray:
    g = pd.read_csv(SRC / "results/eval/gene_seed_003.csv")
    g = g[(g.col == "m_goal") & (g.role == "negative")]
    seeds = sorted(g.seed.unique())
    p = g.pivot_table(index="gene", columns="seed", values="D_raw")
    z = np.column_stack([(p[s] - p[s].mean()) / p[s].std(ddof=1) for s in seeds])
    return z[~np.isnan(z).any(axis=1)]


def main() -> None:
    cfg = yaml.safe_load(open(E / "config.yaml"))
    abs_min = cfg["analysis"]["set_abs_min"]
    sizes = json.loads((E / "work/subsets_005.json").read_text())["sizes"]
    z = znull()
    rng = np.random.default_rng(SEED)
    N, S = z.shape

    def draws(k: int, n: int = N_DRAW) -> np.ndarray:
        return z[rng.integers(0, N, size=(n, k))].mean(axis=1)

    def power(k: int, delta: float, thr: np.ndarray, am: float, tau: float, n: int = 40_000) -> float:
        idx = rng.integers(0, N, size=(n, k))
        eff = rng.normal(delta, tau, size=(n, k)) if tau > 0 else np.full((n, k), float(delta))
        obs = z[idx].mean(axis=1) + eff.mean(axis=1)[:, None]
        return float(((obs > thr[None, :]).all(axis=1) & (obs >= am).all(axis=1)).mean())

    def floor_at(k: int, thr: np.ndarray, am: float, tau: float) -> float:
        lo, hi = 0.0, 6.0
        for _ in range(28):
            mid = (lo + hi) / 2
            if power(k, mid, thr, am, tau) >= POWER_TARGET:
                hi = mid
            else:
                lo = mid
        return (lo + hi) / 2

    res = {"source": "EXP-20260929-001 の対照遺伝子 500 個（目標の軸, B_500 版）",
           "source_caveat": "目標の定義も交絡の補正の版も違う。**事前推定であり、本番後に実測し直す**",
           "n_null_genes": int(N), "n_seeds": int(S), "n_draw": N_DRAW, "percentile": PCT,
           "sigma_between_genes": [round(float(z[:, s].std(ddof=1)), 4) for s in range(S)],
           "seed_correlation": round(float(np.corrcoef(z.T)[0, 1]), 4),
           "tau_between_genes_effect": TAU,
           "calibration_note": "閾値は『対照から復元ありで引いた同じ大きさのセットの 99 パーセンタイル』"
                               "と『絶対値 set_abs_min 以上』を 3 seed すべてで満たす規則。"
                               "偶然合格率 5% に較正した規則ではない（guardrails 2026-09-28）",
           "sets": {}}

    for name, k in sizes.items():
        am = abs_min.get(name)
        dd = draws(k)
        thr = np.quantile(dd, PCT, axis=0)
        if am is None:                      # 中立性のセットには絶対値条件が無い
            res["sets"][name] = {"K": k, "percentile_99": [round(float(t), 4) for t in thr],
                                 "note": "中立性の検査。絶対値条件は無い"}
            continue
        cr = float(((dd > thr[None, :]).all(axis=1) & (dd >= am).all(axis=1)).mean())
        res["sets"][name] = {
            "K": k, "set_abs_min": am,
            "percentile_99": [round(float(t), 4) for t in thr],
            "effective_threshold": round(max(float(thr.max()), am), 4),
            "binding_condition": "絶対値" if am > float(thr.max()) else "99 パーセンタイル",
            "chance_pass_rate": round(cr, 6),
            "detection_floor_by_tau": {str(t): round(floor_at(k, thr, am, t), 3) for t in TAU_ALT},
            "detection_floor_with_gene_sd": round(floor_at(k, thr, am, TAU), 3)}

    # **選定手順の出力 counts.json の検出下限を再現することを確かめる**（1 つでも違えば停止）
    design = json.loads((ROOT / "outputs/design/2026-10-01-intervention-gene-set/counts.json").read_text())
    want = {"INT": design["候補"]["P1_経路セット120"]["検出下限"],
            "INT_S1": design["入れ子の感度分析"]["感度a_S1経路の構成因子"]["検出下限"],
            "INT_S2": design["入れ子の感度分析"]["感度b_S2経路の転写標的"]["検出下限"],
            "INT_S3": design["入れ子の感度分析"]["感度c_S3キュレーション"]["検出下限"],
            "INT_pathway_specific": design["入れ子の感度分析"]["感度d_経路特異的のみ"]["検出下限"],
            "INT_plus_T1T2": design["入れ子の感度分析"]["感度e_主判定とT1T2の和"]["検出下限"]}
    bad = {k: (res["sets"][k]["detection_floor_with_gene_sd"], v) for k, v in want.items()
           if abs(res["sets"][k]["detection_floor_with_gene_sd"] - v) > 0.02}
    if bad:
        raise SystemExit(f"選定手順 counts.json の検出下限を再現できない: {bad}")
    res["counts_json_agreement"] = {k: {"再計算": res["sets"][k]["detection_floor_with_gene_sd"],
                                        "counts.json": v} for k, v in want.items()}

    (E / "work/power_005.json").write_text(json.dumps(res, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "sets"}, indent=1, ensure_ascii=False))
    print(f"\n{'セット':22s} {'K':>3s} {'99%点':>8s} {'絶対値':>7s} {'実効閾値':>9s} "
          f"{'偶然合格率':>10s} {'下限(ばらつき込み)':>18s}")
    for k, v in res["sets"].items():
        if "set_abs_min" not in v:
            print(f"{k:22s} {v['K']:3d} {max(v['percentile_99']):8.3f}   （中立性の検査）")
            continue
        print(f"{k:22s} {v['K']:3d} {max(v['percentile_99']):8.3f} {v['set_abs_min']:7.2f} "
              f"{v['effective_threshold']:9.3f} {v['chance_pass_rate']:10.6f} "
              f"{v['detection_floor_with_gene_sd']:18.3f}")


if __name__ == "__main__":
    main()

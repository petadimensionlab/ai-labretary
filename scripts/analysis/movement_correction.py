#!/usr/bin/env python3
"""移動量の補正 — 主指標の交絡除去に「摂動でどれだけ動いたか」を加える。

## なぜあるか

過剰発現型の摂動では、**恒常発現遺伝子（ACTB・GAPDH・PPIA）が中立にならない**問題が
EXP-20260926-005 / EXP-20260927-003 / EXP-20260927-004 / EXP-20260929-001 の 4 実験で続いた。
恒常発現遺伝子は「細胞の種類によらず常に働いている」ので、目標の状態へ寄る理由がない。
それが寄って見えるなら、指標が何か別のものを測っている。

2026-09-30 に原因が特定された。**これらの遺伝子は摂動したときの移動量が飛び抜けて大きい。**

| 遺伝子 | 移動量 ‖Δ‖ の中央値 | 対照 500 個の中での分位 |
|---|---|---|
| ACTB | 0.2248 | 99.8% |
| GAPDH | 0.1875 | 98.6% |
| PPIA | 0.1770 | 97.8% |
| （対照遺伝子の中央値） | 0.1317 | 50% |

主指標は移動量で割って正規化しているが、**割り切れていない**
（対照 500 個で、目標の軸の指標と移動量の順位相関が +0.171 残る。EXP-20260929-001 validation.md §2.5(a)）。

## 何をするか

交絡を除く二次回帰の説明変数に、**移動量の対数 `ld = log2(‖Δ‖)`** を加える。

| 版 | 説明変数 | 列数 |
|---|---|---|
| 従来（EXP-20260927-004 以降） | 元の位置 `lp`、入力の長さ `ll` の二次 | 6 |
| **本モジュール** | `lp`・`ll`・**`ld`** の二次 | 10 |

`lp = log2(1 + 摂動した遺伝子の元の順位)`、`ll = log2(入力の長さ)`、`ld = log2(‖Δ‖)`。
どちらも同じ形（定数・1 次・2 乗・交差項をすべて入れた完全二次式）なので、
版の違いは「`ld` を入れたかどうか」だけになる。

## 使いどころと限界

- **この補正は「中立であるべきものが中立に戻る」ことを目的とする。**
  目標の軸の陽性対照が立たない問題は**これでは解決しない**（EXP-20260929-001 で 0.718 → 0.698、必要 0.86〜0.90）。
  そちらは目標の定義のやり直し（陽性対照 K ≥ 20）で対処する。
- **結果を見てから版を選んではいけない。** どちらの版を主判定に使うかは事前登録で固定する。
  本モジュールは両方を同時に計算して並べるので、**事前登録で決めた版だけを判定に使う**こと。

## 使い方

    from movement_correction import design, gene_stats

    f_old = gene_stats(cells, nulls, elig, "m_goal", min_elig, use_dnorm=False)
    f_new = gene_stats(cells, nulls, elig, "m_goal", min_elig, use_dnorm=True)

自己検査（合成データで、遅い素直な実装との一致と、補正の効きを確かめる）:

    .venv/bin/python scripts/analysis/movement_correction.py --self-test
"""
from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["design", "gene_stats", "gene_stats_reference", "self_test"]


# ---------------------------------------------------------------------------
# 設計行列
# ---------------------------------------------------------------------------
def design(lp, ll, ld=None) -> np.ndarray:
    """交絡を除く回帰の説明変数を作る。

    lp: 摂動した遺伝子の元の順位の対数 log2(1 + 順位)
    ll: 入力の長さの対数 log2(長さ)
    ld: 移動量の対数 log2(‖Δ‖)。None なら従来どおり lp・ll だけの 6 列を返す。

    どちらも完全二次式（定数・1 次・2 乗・交差項）。列の並びは版の間で前半が一致するので、
    係数を並べて比べられる。
    """
    lp = np.asarray(lp, dtype=float)
    ll = np.asarray(ll, dtype=float)
    base = [np.ones_like(lp), lp, ll, lp ** 2, ll ** 2, lp * ll]
    if ld is None:
        return np.column_stack(base)
    ld = np.asarray(ld, dtype=float)
    if ld.shape != lp.shape:
        raise ValueError(f"ld の長さが lp と違う: {ld.shape} vs {lp.shape}")
    if not np.all(np.isfinite(ld)):
        raise ValueError("ld に有限でない値がある（‖Δ‖ が 0 の細胞が混じっていないか確かめること）")
    return np.column_stack(base + [ld, ld ** 2, lp * ld, ll * ld])


def _cols(cs: pd.DataFrame, use_dnorm: bool):
    """設計行列に渡す列を取り出す。use_dnorm のときだけ移動量の対数を作る。"""
    if not use_dnorm:
        return cs.lp, cs.ll, None
    if "dnorm" not in cs.columns:
        raise KeyError("移動量の列 dnorm が無い。埋め込みから ‖pert − orig‖ を計算して入れること")
    d = np.asarray(cs.dnorm, dtype=float)
    if np.any(d <= 0):
        raise ValueError(f"移動量が 0 以下の細胞が {int((d <= 0).sum())} 個ある。対数が取れない")
    return cs.lp, cs.ll, np.log2(d)


# ---------------------------------------------------------------------------
# 遺伝子ごとの残差
# ---------------------------------------------------------------------------
def gene_stats(cs: pd.DataFrame, nulls: dict, elig: dict, col: str, min_elig: int,
               use_dnorm: bool = True) -> pd.DataFrame:
    """遺伝子ごとに、対照遺伝子（自分自身は除く）で当てた回帰からの残差の平均 D_raw と、
    対照遺伝子の中での順位 q を返す。

    EXP-20260928-001 / EXP-20260929-001 の実装を、設計行列だけ差し替えられるようにしたもの。
    速度のため正規方程式で解き、対照全体の X'X・X'y から自分の寄与を引く
    （素直な実装との一致は self_test() で確かめている）。
    """
    hs = list(nulls)
    Xn, yn = {}, {}
    for h, hd in nulls.items():
        a, b, c = _cols(hd, use_dnorm)
        Xn[h] = design(a, b, c)
        yn[h] = hd[col].to_numpy()
    G = {h: Xn[h].T @ Xn[h] for h in hs}
    B = {h: Xn[h].T @ yn[h] for h in hs}
    Xbar = np.vstack([Xn[h].mean(axis=0) for h in hs])
    ybar = np.array([yn[h].mean() for h in hs])
    hi = {h: i for i, h in enumerate(hs)}
    G_all = sum(G.values())
    B_all = sum(B.values())
    rows = []
    for gene, tg in cs.groupby("gene"):
        Gm, Bm = (G_all - G[gene], B_all - B[gene]) if gene in G else (G_all, B_all)
        beta = np.linalg.solve(Gm, Bm)
        a, b, c = _cols(tg, use_dnorm)
        Xg = design(a, b, c)
        Dg = float(tg[col].to_numpy().mean() - Xg.mean(axis=0) @ beta)
        e = elig[gene]
        if e:
            ii = np.array([hi[h] for h in e])
            Dh = ybar[ii] - Xbar[ii] @ beta
        else:
            Dh = np.zeros(0)
        q = float((Dh < Dg).mean()) if len(Dh) >= min_elig else np.nan
        rows.append({"gene": gene, "col": col, "D_raw": Dg, "q": q, "n_elig": len(Dh),
                     "n_cells": len(tg), "lp_med": float(tg.lp.median()),
                     "ll_med": float(tg.ll.median()),
                     # relpos_med は交絡の診断（M2）が使う。EXP-20260929-001 の実装にはあったが、
                     # 2026-09-30 に共有モジュールへ切り出したとき写し忘れ、評価が KeyError で止まった。
                     # **判定の閾値にも遺伝子セットにも関係しない、診断用の集計値である。**
                     "relpos_med": float(tg.relpos.median()) if "relpos" in tg else np.nan,
                     "dnorm_med": float(tg.dnorm.median()) if "dnorm" in tg else np.nan})
    return pd.DataFrame(rows)


#: 呼び出し側（evaluate_*.py）が必要とする列。切り出し元の EXP-20260929-001 の実装に合わせる。
#: **共有モジュールへ切り出すときは、元の実装の出力列と突き合わせること**（2026-09-30 の教訓）。
REQUIRED_COLUMNS = ("gene", "col", "D_raw", "q", "n_elig", "n_cells",
                    "lp_med", "ll_med", "relpos_med")


def gene_stats_reference(cs: pd.DataFrame, nulls: dict, elig: dict, col: str, min_elig: int,
                         use_dnorm: bool = True) -> pd.DataFrame:
    """素直な実装（遅い）。self_test() の照合にだけ使う。"""
    Xn, yn = {}, {}
    for h, hd in nulls.items():
        a, b, c = _cols(hd, use_dnorm)
        Xn[h] = design(a, b, c)
        yn[h] = hd[col].to_numpy()
    rows = []
    for gene, tg in cs.groupby("gene"):
        others = [h for h in nulls if h != gene]
        X = np.vstack([Xn[h] for h in others])
        y = np.concatenate([yn[h] for h in others])
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        a, b, c = _cols(tg, use_dnorm)
        Dg = float((tg[col].to_numpy() - design(a, b, c) @ beta).mean())
        Dh = np.array([float((yn[h] - Xn[h] @ beta).mean()) for h in elig[gene]])
        q = float((Dh < Dg).mean()) if len(Dh) >= min_elig else np.nan
        rows.append({"gene": gene, "col": col, "D_raw": Dg, "q": q, "n_elig": len(Dh),
                     "n_cells": len(tg)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 自己検査
# ---------------------------------------------------------------------------
def self_test() -> None:
    """(a) 速い実装と素直な実装が一致すること、(b) 移動量で動く交絡を入れると
    従来の補正では残り、移動量を入れた補正では消えること、を合成データで確かめる。"""
    rng = np.random.default_rng(0)

    # --- 合成データ: 指標が「移動量の対数」に線形に依存する（＝交絡）
    rows = []
    for gi in range(60):
        n = int(rng.integers(120, 201))
        lp = rng.normal(5, 1, n)
        ll = rng.normal(11, 0.3, n)
        # 遺伝子ごとに移動量の水準が違う（恒常発現遺伝子が極端に大きい状況を模す）
        lvl = rng.normal(0, 0.4)
        dnorm = np.exp2(rng.normal(-3.0 + lvl, 0.15, n))
        m = 0.3 * lp - 0.2 * ll + 0.8 * np.log2(dnorm) + rng.normal(0, 0.5, n)
        rows.append(pd.DataFrame({"gene": f"G{gi:03d}", "lp": lp, "ll": ll,
                                  "dnorm": dnorm, "m": m}))
    cs = pd.concat(rows, ignore_index=True)
    cs["relpos"] = cs.lp / (cs.ll + 1e-9)        # 呼び出し側は pos / length で作る。ここでは形だけ合わせる
    nulls = {g: d for g, d in cs.groupby("gene") if g < "G050"}     # 対照 50 個
    elig = {g: [h for h in nulls if h != g] for g in cs.gene.unique()}

    # (0) 呼び出し側が必要とする列がそろっているか
    # **これが無かったために、EXP-20260930-001 の評価が 8.6 時間の計算の後で止まった。**
    # モジュール単体の検査だけでは、呼び出し側の期待とのずれを捕まえられない。
    cols = set(gene_stats(cs, nulls, elig, "m", 20, use_dnorm=True).columns)
    missing = [c for c in REQUIRED_COLUMNS if c not in cols]
    assert not missing, f"呼び出し側が必要とする列が足りない: {missing}"
    print(f"(0) 必要な列がそろっている（{len(REQUIRED_COLUMNS)} 列）  OK")

    # (a) 実装の一致（両方の版で）
    for use in (False, True):
        fast = gene_stats(cs, nulls, elig, "m", 20, use_dnorm=use).sort_values("gene")
        ref = gene_stats_reference(cs, nulls, elig, "m", 20, use_dnorm=use).sort_values("gene")
        d = float(np.max(np.abs(fast.D_raw.to_numpy() - ref.D_raw.to_numpy())))
        assert d < 1e-9, f"use_dnorm={use} で速い実装と素直な実装がずれる: 最大差 {d:.3e}"
        print(f"(a) 実装の一致 use_dnorm={use}: 最大差 {d:.2e}  OK")

    # (b) 交絡が消えるか — 残差と移動量の順位相関で見る
    from scipy import stats
    out = {}
    for use in (False, True):
        f = gene_stats(cs, nulls, elig, "m", 20, use_dnorm=use)
        out[use] = abs(float(stats.spearmanr(f.D_raw, f.dnorm_med).statistic))
    print(f"(b) 残差と移動量の順位相関: 従来 {out[False]:.3f} → 移動量を入れる {out[True]:.3f}")
    assert out[True] < out[False] / 2, "移動量を入れても交絡が半分以下に減っていない"
    assert out[True] < 0.15, f"移動量を入れた後も相関が大きい: {out[True]:.3f}"
    print("(b) 交絡の除去  OK")

    # (c) 交絡が無いデータでは、余計な項を入れても結果がほとんど変わらないこと
    rows = []
    for gi in range(60):
        n = int(rng.integers(120, 201))
        lp = rng.normal(5, 1, n)
        ll = rng.normal(11, 0.3, n)
        dnorm = np.exp2(rng.normal(-3.0, 0.15, n))
        m = 0.3 * lp - 0.2 * ll + rng.normal(0, 0.5, n)          # 移動量に依存しない
        rows.append(pd.DataFrame({"gene": f"G{gi:03d}", "lp": lp, "ll": ll,
                                  "dnorm": dnorm, "m": m}))
    cs2 = pd.concat(rows, ignore_index=True)
    nulls2 = {g: d for g, d in cs2.groupby("gene") if g < "G050"}
    elig2 = {g: [h for h in nulls2 if h != g] for g in cs2.gene.unique()}
    a = gene_stats(cs2, nulls2, elig2, "m", 20, use_dnorm=False).sort_values("gene")
    b = gene_stats(cs2, nulls2, elig2, "m", 20, use_dnorm=True).sort_values("gene")
    r = float(stats.spearmanr(a.D_raw, b.D_raw).statistic)
    print(f"(c) 交絡が無いとき、2 版の残差の順位相関: {r:.3f}")
    assert r > 0.95, f"交絡が無いのに版で結果が変わりすぎる: {r:.3f}"
    print("(c) 余計な項の悪影響が小さいこと  OK")

    print("\n自己検査はすべて合格")


if __name__ == "__main__":
    import sys
    if "--self-test" in sys.argv:
        self_test()
    else:
        print(__doc__)

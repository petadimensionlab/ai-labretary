"""EXP-20260928-002 — 検出下限を「事前登録する規則そのもの」で計算する（guardrails 2026-09-28）。

guardrails 2026-09-28:
  「検出下限を書くときは『どの偶然合格率に較正した規則での値か』を必ず添える。
    事前登録に書く下限は、事前登録する規則そのもので計算する。」
  EXP-20260928-001 では 5% 較正の値（0.32）を書いてしまい、実際の規則（99 パーセンタイル）での
  下限は 0.59〜0.83 だった。同じ誤りを繰り返さないため、本スクリプトは
  **config.yaml の criteria をそのまま読んで**、その規則での下限を出す。

事前登録する規則（config.yaml の analysis.criteria）:
  3 seed すべてで、セット平均 z が
    (i) 対照遺伝子から **復元ありで** 引いた同じ大きさのセットの 99 パーセンタイルを超え、かつ
    (ii) set_abs_min（0.245）以上

使うデータ: EXP-20260928-001 の対照遺伝子 96 個 × 3 seed の主指標の実測値だけ
  （results/eval/gene_seed_001.csv の col=m_goal）。本実験の摂動は 1 回も実行していない。
  遺伝子ごとの効果のばらつき tau は、EXP-20260928-001 の目標マーカー 24 個の z の散らばりから推定する。

出力: work/power_002.json
実行: .venv/bin/python experiments/EXP-20260928-002/analysis/power_002.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260928-002"
SRC1 = ROOT / "experiments/EXP-20260928-001"
SRC1_SEEDS = [41, 42, 43]
N_DRAW = 100_000
N_POWER = 20_000
RNG = np.random.default_rng(20260929)


def zmat(col: str) -> pd.DataFrame:
    gs = pd.read_csv(SRC1 / "results/eval/gene_seed_001.csv")
    P = gs[gs.col == col]
    rows = []
    for s in SRC1_SEEDS:
        ps = P[P.seed == s]
        nul = ps[ps.role == "negative"]
        mu, sd = nul.D_raw.mean(), nul.D_raw.std(ddof=1)
        for _, r in ps.iterrows():
            rows.append({"symbol": r.symbol, "role": r.role, "seed": s, "z": (r.D_raw - mu) / sd})
    Z = pd.DataFrame(rows)
    return Z.pivot_table(index=["symbol", "role"], columns="seed", values="z").reset_index()


def sample_idx(N: int, k: int, n: int, replace: bool, pool_n: int | None,
               rng: np.random.Generator) -> np.ndarray:
    """対照遺伝子 N 個から k 個のセットを n 組つくり、A への行番号（n × k）を返す。

    replace=True  : 復元あり（事前登録する引き方。有限母集団の補正が入らない）
    replace=False : 復元なし。pool_n を指定したときは、実測の N 個から復元ありで pool_n 個の
                    「合成された母集団」を作り、そこから復元なしで k 個を引く。
                    対照がセットの何倍あるかによる閾値の甘さを近似するため（guardrails 2026-09-28）。
    """
    if replace:
        return rng.integers(0, N, size=(n, k))
    if pool_n is None:
        return np.argsort(rng.random((n, N)), axis=1)[:, :k]
    out = np.empty((n, k), dtype=np.int64)
    block = 2000
    for i in range(0, n, block):
        m = min(block, n - i)
        pool = rng.integers(0, N, size=(m, pool_n))
        take = np.argsort(rng.random((m, pool_n)), axis=1)[:, :k]
        out[i:i + m] = np.take_along_axis(pool, take, axis=1)
    return out


def main() -> None:
    cfg = yaml.safe_load((E / "config.yaml").read_text())
    cr = cfg["analysis"]["criteria"]
    pct, abs_min = cr["set_percentile"], cr["set_abs_min"]
    out: dict = {
        "規則": f"3 seed すべてで、セット平均 z が対照セットの {pct*100:.0f} パーセンタイル超 かつ {abs_min} 以上",
        "較正": f"閾値は『対照から復元ありで引いた同じ大きさのセットの {pct*100:.0f} パーセンタイル』であり、"
                f"偶然合格率 5% に較正した規則ではない（guardrails 2026-09-28）",
        "出典データ": "EXP-20260928-001 results/eval/gene_seed_001.csv（col=m_goal）の対照遺伝子 96 個 × 3 seed。"
                      "本実験の摂動は 1 回も実行していない。",
    }

    W = zmat("m_goal")
    A = W[W.role == "negative"][SRC1_SEEDS].to_numpy()
    mk = W[W.role == "marker_goal"][SRC1_SEEDS].to_numpy()
    N = A.shape[0]
    tau = float(mk.mean(axis=1).std(ddof=1))     # 遺伝子ごとの効果のばらつき（3 seed 平均の遺伝子間 SD）
    out["対照遺伝子"] = {"n": int(N),
                         "seed 内の標準偏差": [round(float(A[:, i].std(ddof=1)), 3) for i in range(3)],
                         "同じ遺伝子の seed 間の標準偏差（中央値）": round(float(np.median(A.std(1, ddof=1))), 3)}
    out["遺伝子ごとの効果のばらつき tau"] = {
        "値": round(tau, 3),
        "根拠": "EXP-20260928-001 の目標マーカー 24 個の（3 seed 平均の）z の遺伝子間標準偏差。"
                "validation.md §4.2(a) の 2.05 と同じ量。"}

    sets = json.loads((E / "work/subsets_002.json").read_text())["subsets"]
    ks = {"pathway_all (主判定)": len(sets["pathway_all"]),
          "nonmarker (副次)": len(sets["nonmarker"]),
          "drop_goalmarker (副次)": len(sets["drop_goalmarker"]),
          "pathway_specific (副次)": len(sets["pathway_specific"]),
          "tier_S1": len(sets["tier_S1"]), "tier_S2": len(sets["tier_S2"]), "tier_S3": len(sets["tier_S3"]),
          "myof_markers 24 (アッセイの妥当性)": 24}

    def limits(k: int, replace: bool, pool_n: int | None) -> dict:
        m = A[sample_idx(N, k, N_DRAW, replace, pool_n, RNG)].mean(axis=1)
        thr = np.quantile(m, pct, axis=0)
        eff = np.maximum(thr, abs_min)                      # 実効の閾値（絶対値条件が効くことが多い）
        chance = float(((m > thr) & (m >= abs_min)).all(axis=1).mean())

        def power(mu: float, t: float) -> float:
            idx = sample_idx(N, k, N_POWER, replace, pool_n, RNG)
            e = RNG.normal(mu, t, size=(N_POWER, k))        # 遺伝子ごとの効果（seed をまたいで固定）
            mm = (A[idx] + e[:, :, None]).mean(axis=1)
            return float(((mm > thr) & (mm >= abs_min)).all(axis=1).mean())

        res = {"閾値（seed ごとの 99 パーセンタイル）": [round(float(x), 3) for x in thr],
               "実効の閾値（絶対値条件 %.3f と合わせて）" % abs_min: [round(float(x), 3) for x in eff],
               "偶然合格率（対照セットに同じ規則）": round(chance, 5)}
        for label, t in (("効果が均一", 0.0), (f"遺伝子間 SD {tau:.2f} 込み", tau)):
            grid = np.arange(0.0, 3.0, 0.01)
            lo, hi = 0, len(grid) - 1
            while lo < hi:                                   # 検出力は mu について単調
                mid = (lo + hi) // 2
                if power(float(grid[mid]), t) >= 0.80:
                    hi = mid
                else:
                    lo = mid + 1
            res[f"検出下限（検出力 80%、{label}）"] = round(float(grid[lo]), 3)
        return res

    out["検出下限"] = {}
    for label, k in ks.items():
        out["検出下限"][f"K={k}  {label}"] = {
            "復元あり（事前登録する引き方）": limits(k, True, None),
            "復元なし・対照 980（主判定の版）": limits(k, False, 980),
            "復元なし・対照 600（案 B 相当）": limits(k, False, 600),
            "復元なし・対照 520（案 C 相当）": limits(k, False, 520),
        }

    # 閾値の差（復元あり vs なし）を一覧にする（PI 指示 2026-09-28）
    tbl = {}
    for label, k in ks.items():
        r = out["検出下限"][f"K={k}  {label}"]
        tbl[f"K={k}"] = {
            "復元あり": r["復元あり（事前登録する引き方）"]["閾値（seed ごとの 99 パーセンタイル）"],
            "復元なし 980": r["復元なし・対照 980（主判定の版）"]["閾値（seed ごとの 99 パーセンタイル）"],
            "復元なし 600": r["復元なし・対照 600（案 B 相当）"]["閾値（seed ごとの 99 パーセンタイル）"],
            "復元なし 520": r["復元なし・対照 520（案 C 相当）"]["閾値（seed ごとの 99 パーセンタイル）"]}
    out["閾値の比較（復元あり / なし × 対照の数）"] = tbl
    out["読み方"] = (
        "絶対値条件 0.245 が効くため、K を増やしても実効の閾値は 0.245 で頭打ちになる。"
        "復元なしにすると対照の数が少ないほど閾値が甘くなる（有限母集団の補正）。"
        "主判定は復元ありで引く（guardrails 2026-09-28）。"
        "**不合格は『効果がない』ではなく『この下限以上の効果は無い』としか読めない。**")

    (E / "work").mkdir(exist_ok=True)
    (E / "work/power_002.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps({"tau": out["遺伝子ごとの効果のばらつき tau"]["値"],
                      "閾値の比較": tbl,
                      "主判定 K=120": out["検出下限"]["K=120  pathway_all (主判定)"]},
                     ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()

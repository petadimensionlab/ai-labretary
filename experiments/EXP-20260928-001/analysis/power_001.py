"""EXP-20260928-001 — 検出力の見積もり（実データのばらつきだけを使う。合成データは使わない）。

問い: EXP-20260927-004 で測った「検出下限 2.73 標準偏差」の床は、遺伝子セット単位の検定（A-3）で
      本当に下げられるのか。下げられないなら、この方向は採らない（PI 指示 2026-09-28 の要件 2）。

使うデータ: EXP-20260927-004 の対照遺伝子 84 個 × 3 seed の主指標の実測値のみ
            （results/eval/gene_seed_004.csv の col=m_primary）。摂動は新たに実行しない。

出力: work/power_001.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "experiments/EXP-20260927-004"
OUT = ROOT / "experiments/EXP-20260928-001/work"
SEEDS = [31, 32, 33]
KS = [5, 8, 10, 12, 15, 20, 24, 30]
N_DRAW = 200_000
RNG = np.random.default_rng(20260928)


def zmat(col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """遺伝子 × seed の標準偏差単位の表を返す（対照遺伝子群の平均・標準偏差で標準化）。"""
    gs = pd.read_csv(SRC / "results/eval/gene_seed_004.csv")
    P = gs[gs.col == col]
    rows = []
    for s in SEEDS:
        ps = P[P.seed == s]
        nul = ps[ps.role == "negative"]
        mu, sd = nul.D_raw.mean(), nul.D_raw.std(ddof=1)
        for _, r in ps.iterrows():
            rows.append({"symbol": r.symbol, "role": r.role, "seed": s, "z": (r.D_raw - mu) / sd})
    Z = pd.DataFrame(rows)
    W = Z.pivot_table(index=["symbol", "role"], columns="seed", values="z").reset_index()
    return Z, W


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    out: dict = {"note": "EXP-20260927-004 の実測ばらつきのみ。合成データは使っていない（guardrails 2026-09-27）。"}

    Z, W = zmat("m_primary")
    nul = W[W.role == "negative"].reset_index(drop=True)
    N = len(nul)
    A = nul[SEEDS].to_numpy()                     # 対照遺伝子 84 × 3 seed の z
    out["null_pool"] = {"n_genes": int(N),
                        "seed 内の標準偏差": [round(float(A[:, i].std(ddof=1)), 3) for i in range(3)],
                        "3 seed 平均の標準偏差": round(float(A.mean(1).std(ddof=1)), 3),
                        "同じ遺伝子の seed 間の標準偏差（中央値）": round(float(np.median(A.std(1, ddof=1))), 3)}

    # ---- 単一遺伝子の床（EXP-004 の 2.73 標準偏差の再現）
    # 規則「3 seed すべてで z >= t」で、対照 84 個の偶然合格率が 0.05 以下になる t を求める
    def chance_rate(t: float, mat: np.ndarray) -> float:
        return float((mat >= t).all(axis=1).mean())

    t_grid = np.arange(0.5, 6.01, 0.01)
    t1 = next(float(t) for t in t_grid if chance_rate(t, A) <= 0.05)
    out["single_gene"] = {
        "偶然合格率 5% 以下にする閾値（3 seed すべて）": round(t1, 2),
        "EXP-004 の事前登録 2.5 標準偏差での偶然合格率": round(chance_rate(2.5, A), 4),
        "EXP-004 が報告した検出下限（検出力 80%）": 2.73,
        "注": "単一遺伝子の床は標準化された量なので、細胞数・seed 数では下がらない（EXP-004 validation.md §2.9）。"}

    # ---- セット単位の検定: 対照遺伝子から K 個を無作為に取った集合の「z の平均」の分布
    sets_tab = {}
    for K in KS:
        idx = RNG.integers(0, N, size=(N_DRAW, K))          # 復元なしに近い近似が必要なので下で入れ替え
        # 復元なしの K 部分集合
        idx = np.argsort(RNG.random((N_DRAW, N)), axis=1)[:, :K]
        m = A[idx].mean(axis=1)                              # N_DRAW × 3 seed
        thr95 = float(np.quantile(m.max(axis=1) * 0 + m.min(axis=1), 0.95))  # 3 seed すべて超える規則 → min を見る
        # 「3 seed すべてで閾値 t を超える」規則の偶然合格率が 0.05 になる t
        t_set = next(float(t) for t in np.arange(0.0, 3.001, 0.005)
                     if float((m >= t).all(axis=1).mean()) <= 0.05)
        # 検出力 80%: 真のセット平均 μ のとき、seed ごとのゆらぎを含めて 3 seed すべてが t を超える確率
        # セットの seed 間ゆらぎは、対照遺伝子セットの seed 間ゆらぎの実測で代用する
        seed_noise = float(np.median((A[idx] - A[idx].mean(2, keepdims=True)).std(2, ddof=1).mean(1)) / np.sqrt(K))
        mu_grid = np.arange(0.0, 3.0, 0.005)
        # 近似: セット平均は 3 seed それぞれ N(μ, seed_noise) → 3 seed すべて > t の確率
        from scipy.stats import norm
        pw = (1 - norm.cdf((t_set - mu_grid) / max(seed_noise, 1e-6))) ** 3
        mu80 = float(mu_grid[np.argmax(pw >= 0.80)])
        sets_tab[K] = {"偶然合格率 5% の閾値（セット平均 z、3 seed すべて）": round(t_set, 3),
                       "セット平均の seed 内標準偏差": round(float(m[:, 0].std(ddof=1)), 3),
                       "セットの seed 間ゆらぎ": round(seed_noise, 4),
                       "検出下限（検出力 80%、セット平均 z）": round(mu80, 3)}
    out["set_level"] = sets_tab

    # ---- 参考: EXP-20260927-004 の陽性対照候補 10 のセット平均（効果量の事前見積もりに使う）
    cand = W[W.role == "positive_goal"]
    out["prior_effect_estimate"] = {
        "EXP-004 陽性対照候補 10 のセット平均 z（3 seed）": [round(float(cand[s].mean()), 3) for s in SEEDS],
        "SYN3 を除く 9 個": [round(float(cand[cand.symbol != "SYN3"][s].mean()), 3) for s in SEEDS],
        "文献層 5 個": [round(float(cand[cand.symbol.isin(
            ["VCAM1", "ICAM1", "CXCL12", "IL7", "TNFSF13B"])][s].mean()), 3) for s in SEEDS],
        "データ層 5 個": [round(float(cand[cand.symbol.isin(
            ["FKBP5", "SLC2A3", "CCN1", "SYN3", "STEAP4"])][s].mean()), 3) for s in SEEDS],
        "注": "これは EXP-004 の結果を見た後の値であり、本実験の判定基準には使わない。"
              "検出力の見積もり（必要な K の決定）にだけ使う。"}

    # ---- 軸を替えても床が下がるか（A-2 の評価）: 筋線維芽細胞の軸での同じ計算
    Zm, Wm = zmat("m_myof")
    nm = Wm[Wm.role == "negative"][SEEDS].to_numpy()
    t1m = next(float(t) for t in t_grid if chance_rate(t, nm) <= 0.05)
    out["axis_dependence"] = {
        "目標の軸での単一遺伝子の床（偶然 5%）": round(t1, 2),
        "筋線維芽細胞の軸での単一遺伝子の床（偶然 5%）": round(t1m, 2),
        "筋線維芽細胞の軸での ACTA2（3 seed, 標準偏差単位）":
            [round(float(Wm[Wm.symbol == "ACTA2"][s].iloc[0]), 2) for s in SEEDS],
        "筋線維芽細胞の軸での TAGLN": [round(float(Wm[Wm.symbol == "TAGLN"][s].iloc[0]), 2) for s in SEEDS],
        "読み方": "床（標準偏差単位）は軸をどう選んでも ほぼ同じ。軸の良し悪しが変えるのは信号の側だけで、"
                  "その伸びは既存データからは見積もれない。"}

    # ---- A-1（複数遺伝子の同時摂動）の効果の上限を、既存の移動方向から予測する
    dd = pd.read_csv(OUT / "gene_direction_001.csv")
    # 方向の平均コサイン rho（confound_probe_001.py で測った全ペア中央値）を使った飽和則
    rho = json.loads((OUT / "confound_probe_001.json").read_text())["pairwise_direction"]["全ペアの cos 中央値"]
    gains = {k: round(float(np.sqrt(k) / np.sqrt(1 + (k - 1) * rho)), 3) for k in [2, 3, 5, 10, 20]}
    out["A1_combined_perturbation_prediction"] = {
        "遺伝子間の移動方向の平均コサイン rho": round(float(rho), 3),
        "k 個を同時に動かしたときの『軸との内積』の伸び（単一遺伝子の平均に対する倍率）": gains,
        "k → 無限大での上限": round(float(1 / np.sqrt(rho)), 3),
        "読み方": "方向が互いに無相関なら sqrt(k) で伸びるが、実測 rho = 0.46 なので倍率は "
                  "k=5 で 1.33、k=10 で 1.40、上限 1.48 で飽和する。"
                  "一方、無作為な k 遺伝子セットを対照にすると対照側も同じ倍率を受けるため、"
                  "標準偏差単位での利得は A-3（k 遺伝子を個別に動かしてセット平均を取る）と同じ sqrt(k) になる。"}

    (OUT / "power_001.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()

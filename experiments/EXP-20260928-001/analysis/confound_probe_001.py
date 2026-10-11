"""EXP-20260928-001 — 推奨 B の調査: SYN3 と RASA1 が同じ規則を満たした理由を、既存データから調べる。

読むもの（すべて EXP-20260927-004 の既存の結果。新しい摂動は一切しない）:
  results/eval/gene_seed_004.csv         指標ごとの D_raw / q / 元の位置・長さ
  results/run-{31,32,33}-L-1/genes/<ENSG>/emb/{orig,pert}.npy   摂動前後の CLS 埋め込み
  results/run-{31,32,33}-L-1/state_emb_vectors.npz              目標・出発点・対照の状態の重心
  work/extra_centroids-seed{31,32,33}.npz                       筋線維芽細胞の重心
  results/cells.csv                                             細胞 × 遺伝子の対応（提供者・長さ）

調べる観点:
  1. 主指標 m_primary は、一次近似で「移動の向き（単位ベクトル）と、固定の軸 w との内積」になる。
     w = P_o(goal_hat − (start_hat + alt_hat + myof_hat)/3) / ‖o‖
     → 遺伝子の効果量は「移動の向き」だけで決まる。向きが何に決まっているかを調べる。
  2. 遺伝子ごとの平均移動方向 d_hat の遺伝子間コサイン。SYN3 と RASA1 は同じ向きに動いているか。
  3. d_hat の seed 間再現性（遺伝子固有の向きがどれだけ安定か）。
  4. a_g = cos(d_hat_g, w_hat) と実際の z（標準偏差単位）の対応。近似が成り立つかの確認。
  5. 上位遺伝子の向きが 1 つの共通方向に集まっているか（第 1 主成分の寄与）。
  6. 元の位置・入力の長さ・細胞数・移動量との相関（M2 で見ている 3 つ以外の交絡の探索）。

出力: work/confound_probe_001.json, work/gene_direction_001.csv
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "experiments/EXP-20260927-004"
OUT = ROOT / "experiments/EXP-20260928-001/work"
SEEDS = [31, 32, 33]

# 盲検・主対象は本調査でも一切扱わない（EXP-004 は回していないので現れないはずだが、明示して確認する）
BLOCKED = {"ENSG00000111321", "ENSG00000104856", "ENSG00000077150", "ENSG00000124762",
           "ENSG00000141510", "ENSG00000118785", "ENSG00000148346"}


def unit(v: np.ndarray) -> np.ndarray:
    return v / (np.linalg.norm(v) + 1e-30)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gs = pd.read_csv(SRC / "results/eval/gene_seed_004.csv")
    assert not set(gs.gene) & BLOCKED, "回さないはずの遺伝子が入っている"
    P = gs[gs.col == "m_primary"].copy()
    sym = dict(zip(gs.gene, gs.symbol))
    role = dict(zip(gs.gene, gs.role))

    # 標準偏差単位 z（対照遺伝子群の平均・標準偏差で標準化。evaluate_004.py と同じ）
    zrows = []
    for s in SEEDS:
        ps = P[P.seed == s]
        nul = ps[ps.role == "negative"]
        mu, sd = nul.D_raw.mean(), nul.D_raw.std(ddof=1)
        for _, r in ps.iterrows():
            zrows.append({"gene": r.gene, "seed": s, "z": (r.D_raw - mu) / sd, "q": r.q,
                          "lp_med": r.lp_med, "ll_med": r.ll_med, "relpos_med": r.relpos_med,
                          "n_cells": r.n_cells})
    Z = pd.DataFrame(zrows)

    # ---- 遺伝子ごとの平均移動方向と、軸 w との整合
    rec: dict[tuple[str, int], dict] = {}
    dirs: dict[int, dict[str, np.ndarray]] = {s: {} for s in SEEDS}
    for s in SEEDS:
        run = SRC / f"results/run-{s}-L-1"
        sz = np.load(run / "state_emb_vectors.npz")
        xz = np.load(SRC / f"work/extra_centroids-seed{s}.npz")
        cen = {"goal": sz["FIB_org"].astype(np.float64), "start": sz["FIB_neg_TLShigh"].astype(np.float64),
               "alt": sz["FIB_neg_TLSlow"].astype(np.float64), "myof": xz["MYOF_neg"].astype(np.float64)}
        ch = {k: unit(v) for k, v in cen.items()}
        wraw = ch["goal"] - (ch["start"] + ch["alt"] + ch["myof"]) / 3.0
        for gd in sorted((run / "genes").glob("ENSG*")):
            g = gd.name
            if g in BLOCKED:
                raise SystemExit(f"回さないはずの遺伝子の出力がある: {g}")
            o = np.load(gd / "emb/orig.npy").astype(np.float64)
            p = np.load(gd / "emb/pert.npy").astype(np.float64)
            d = p - o
            # 細胞ごとの軸 w_i（o に直交する成分）。細胞平均を遺伝子の軸とする
            oh = o / (np.linalg.norm(o, axis=1, keepdims=True) + 1e-30)
            w = (wraw[None, :] - (oh @ wraw)[:, None] * oh) / (np.linalg.norm(o, axis=1, keepdims=True) + 1e-30)
            dn = np.linalg.norm(d, axis=1, keepdims=True) + 1e-30
            dhat_cells = d / dn
            a_cells = np.einsum("ij,ij->i", dhat_cells, w)           # 一次近似の m_primary
            dbar = unit(d.mean(axis=0))
            dirs[s][g] = dbar
            rec[(g, s)] = {
                "gene": g, "seed": s, "symbol": sym.get(g), "role": role.get(g),
                "a_lin_mean": float(a_cells.mean()),                  # 一次近似での主指標（細胞平均）
                "cos_dbar_w": float(np.dot(dbar, unit(w.mean(axis=0)))),
                "dnorm_mean": float(np.linalg.norm(d, axis=1).mean()),
                "coherence": float(np.linalg.norm(d.mean(axis=0)) / np.linalg.norm(d, axis=1).mean()),
                "n_cells_emb": int(len(d)),
            }

    R = pd.DataFrame(rec.values()).merge(Z, on=["gene", "seed"], how="left")
    R.to_csv(OUT / "gene_direction_001.csv", index=False)

    out: dict = {"note": "EXP-20260927-004 の既存データのみを使った事後の調査。新しい摂動は行っていない。"}

    # ---- 1. 一次近似が主指標を再現するか
    out["linear_approx_check"] = {
        str(s): {"spearman_a_vs_z": float(stats.spearmanr(g.a_lin_mean, g.z).statistic),
                 "pearson_a_vs_z": float(stats.pearsonr(g.a_lin_mean, g.z).statistic),
                 "n_genes": int(len(g))}
        for s, g in R.groupby("seed")}

    # ---- 2. 遺伝子ごとの移動方向の seed 間再現性
    common = sorted(set(dirs[31]) & set(dirs[32]) & set(dirs[33]))
    same = [float(np.dot(dirs[31][g], dirs[32][g])) for g in common]
    out["direction_reproducibility"] = {
        "cos(dbar_seed31, dbar_seed32) 中央値": float(np.median(same)),
        "同 最小": float(np.min(same)), "n_genes": len(common)}

    # ---- 3. 遺伝子間の向きの類似度（seed 31 の dbar で代表）
    M = np.stack([dirs[31][g] for g in common])
    C = M @ M.T
    iu = np.triu_indices(len(common), 1)
    pair_all = C[iu]
    idx = {g: i for i, g in enumerate(common)}
    s2e = {v: k for k, v in sym.items()}

    def pair(a: str, b: str) -> float:
        return float(C[idx[s2e[a]], idx[s2e[b]]])

    top = (R[R.seed == 31].sort_values("z", ascending=False).head(8).symbol.tolist())
    out["pairwise_direction"] = {
        "全ペアの cos 中央値": float(np.median(pair_all)),
        "全ペアの cos の 95 パーセンタイル": float(np.percentile(pair_all, 95)),
        "SYN3 vs RASA1": pair("SYN3", "RASA1"),
        "SYN3 vs IL7": pair("SYN3", "IL7"),
        "SYN3 vs ACTA2": pair("SYN3", "ACTA2"),
        "RASA1 vs IL7": pair("RASA1", "IL7"),
        "SYN3 vs GAPDH": pair("SYN3", "GAPDH"),
        "上位 8 遺伝子（seed31 の z 順）": top,
        "上位 8 の相互 cos 中央値": float(np.median([pair(a, b) for i, a in enumerate(top)
                                                      for b in top[i + 1:]])),
    }
    # 上位 8 の向きが 1 つの共通方向に集まっているか
    T = np.stack([dirs[31][s2e[a]] for a in top])
    sv = np.linalg.svd(T, compute_uv=False)
    out["pairwise_direction"]["上位 8 の第 1 主成分の寄与率"] = float((sv[0] ** 2) / (sv ** 2).sum())
    Mall = M - 0.0
    sva = np.linalg.svd(Mall, compute_uv=False)
    out["pairwise_direction"]["全 106 の第 1 主成分の寄与率"] = float((sva[0] ** 2) / (sva ** 2).sum())

    # ---- 4. 主指標の遺伝子順位と、向き以外の量との相関（M2 の 3 つ以外の交絡の探索）
    cors = {}
    for s in SEEDS:
        g = R[(R.seed == s) & (R.role == "negative")]
        cors[str(s)] = {k: float(stats.spearmanr(g[k], g.z).statistic)
                        for k in ("dnorm_mean", "coherence", "lp_med", "ll_med", "relpos_med", "n_cells")}
    out["null_gene_correlates_of_z"] = cors

    # ---- 5. SYN3 / RASA1 と、比較対象の遺伝子の素性
    picks = ["SYN3", "RASA1", "IL7", "CXCL12", "EPB41L3", "CCSER1", "CHCHD10", "PPA1",
             "ACTA2", "TAGLN", "GAPDH", "FKBP5", "CCN1", "SLC2A3"]
    tab = {}
    for a in picks:
        e = s2e.get(a)
        if e is None:
            continue
        r = R[R.gene == e]
        tab[a] = {"role": role.get(e),
                  "z_3seed": [round(float(x), 2) for x in r.sort_values("seed").z],
                  "cos_dbar_w_3seed": [round(float(x), 4) for x in r.sort_values("seed").cos_dbar_w],
                  "dnorm_mean": round(float(r.dnorm_mean.mean()), 5),
                  "coherence": round(float(r.coherence.mean()), 3),
                  "lp_med": round(float(r.lp_med.mean()), 2),
                  "ll_med": round(float(r.ll_med.mean()), 2),
                  "relpos_med": round(float(r.relpos_med.mean()), 3),
                  "n_cells": int(r.n_cells.mean())}
    out["gene_profiles"] = tab

    # ---- 6. 「向きの整合」だけで z を説明できるか（cos_dbar_w → z の単回帰の決定係数）
    fit = {}
    for s in SEEDS:
        g = R[R.seed == s]
        r = stats.pearsonr(g.cos_dbar_w, g.z)
        fit[str(s)] = {"pearson": float(r.statistic), "r2": float(r.statistic ** 2)}
    out["z_explained_by_direction"] = fit

    (OUT / "confound_probe_001.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()

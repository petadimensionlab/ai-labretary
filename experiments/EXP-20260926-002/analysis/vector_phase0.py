"""EXP-20260926-002 phase 0: 取り出した CLS 埋め込みでベクトル水準の削除ドリフト補正指標を評価する（対照のみ, seed 0）。
ベクトル演算はすべて float64（保存は fp32）。基準は plan.md §6.2（C1〜C3, C4(i), C5）。
実行: .venv/bin/python experiments/EXP-20260926-002/analysis/vector_phase0.py
出力: experiments/EXP-20260926-002/phase0/results-phase0/vector/{criteria.csv,criteria.json,cells_vector.csv,consistency.json}
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reanalysis_cosine as rc  # noqa: E402  同じ候補指標・同じ評価関数を使う

ROOT = rc.ROOT
E = ROOT / "experiments/EXP-20260926-002"
RES = E / "phase0/results-phase0"  # run_isp.py の exp_dir は config の置き場所（phase0/）
RUN = RES / "run-0-L-1"
OUT = RES / "vector"
SRC_CELLS = ROOT / "experiments/EXP-20260926-001/results/tls_high/cells.csv"
K_NN = 50        # 局所ドリフト: 元埋め込みが近い null 細胞 50 個の平均変位
N_PC = 3         # ドリフト部分空間: null 変位の上位 3 主成分


def unit(x):
    return x / np.linalg.norm(x, axis=-1, keepdims=True)


def load():
    cells = pd.read_csv(RES / "cells.csv", dtype={"donor_id": str})
    cells = cells[~cells.gene.isin(rc.BLIND)].reset_index(drop=True)
    assert not set(cells.gene) & rc.BLIND
    E0, E1 = [], []
    for (run, gene), s in cells.groupby(["run", "gene"], sort=False):
        o = np.load(RUN / "genes" / gene / "emb" / "orig.npy").astype(np.float64)
        p = np.load(RUN / "genes" / gene / "emb" / "pert.npy").astype(np.float64)
        per = pd.read_csv(RUN / "genes" / gene / "cells.csv", dtype={"donor_id": str})
        assert len(o) == len(per) == len(s) and (per.cell_uid.to_numpy() == s.cell_uid.to_numpy()).all()
        E0.append(o); E1.append(p)
    # cells の並びを groupby 順に揃える
    order = np.concatenate([s.index.to_numpy() for _, s in cells.groupby(["run", "gene"], sort=False)])
    cells = cells.loc[order].reset_index(drop=True)
    st = np.load(RUN / "state_emb_vectors.npz")
    return cells, np.vstack(E0), np.vstack(E1), {k: st[k].astype(np.float64).ravel() for k in st.files}


def consistency(cells, E0, E1, cen):
    """取り出したベクトルから計算したコサイン変化が、perturber の shift（cells.csv）と一致するか。"""
    g, s = unit(cen["PT_healthy"]), unit(cen["aPT_CKD_TLShigh"])
    e0, e1 = unit(E0), unit(E1)
    sg = e1 @ g - e0 @ g
    ss = e1 @ s - e0 @ s
    out = {"max_abs_diff_goal": float(np.abs(sg - cells.shift_goal).max()),
           "max_abs_diff_start": float(np.abs(ss - cells.shift_start).max()),
           "mean_abs_shift_goal": float(np.abs(cells.shift_goal).mean())}
    # 同じ細胞の元埋め込みは遺伝子によらず同じはず（バッチ構成の違いによる bf16 の揺れを確認）
    dev = []
    for uid, idx in cells.groupby("cell_uid").groups.items():
        if len(idx) > 1:
            x = unit(E0[np.asarray(idx)])
            dev.append(float(1 - (x @ x.T).min()))
    out["orig_emb_same_cell_max_1mcos"] = float(max(dev)) if dev else None
    out["n_cells_multi_gene"] = len(dev)
    # EXP-001 run 1 seed 0 と比べる（同じ細胞・同じ設定で再現するか）
    old = pd.read_csv(SRC_CELLS, dtype={"donor_id": str})
    old = old[(old.run == "run-0-L-1") & ~old.gene.isin(rc.BLIND)]
    m = cells.merge(old, on=["gene", "cell_uid"], suffixes=("", "_old"))
    out["exp001_matched_rows"] = f"{len(m)}/{len(cells)}"
    out["exp001_pearson_shift_goal"] = float(np.corrcoef(m.shift_goal, m.shift_goal_old)[0, 1])
    out["exp001_max_abs_diff_shift_goal"] = float(np.abs(m.shift_goal - m.shift_goal_old).max())
    return out


def vector_metrics(cells, E0, E1, cen):
    """V 系の指標（単位球上の変位 Δ = ê1 − ê0 を使う。コサインが尺度に依らないのと合わせる）。
    null から推定するドリフトは、null 遺伝子自身には leave-one-gene-out。"""
    e0, e1 = unit(E0), unit(E1)
    D = e1 - e0
    g, s = unit(cen["PT_healthy"]), unit(cen["aPT_CKD_TLShigh"])
    u = unit(g - s)                                   # 軸: start → goal
    cells = cells.copy()
    cells["V0"] = D @ u                               # 補正なしの射影
    isnull = (cells.role == "negative").to_numpy()
    genes = cells.gene.to_numpy()
    cols = {k: np.full(len(cells), np.nan) for k in ["V1", "V2", "V3", "V4"]}
    groups = [(genes == gg, isnull & (genes != gg)) for gg in sorted(set(genes[isnull]))]
    groups.append((~isnull, isnull))
    for tgt, fit in groups:
        Dn, Dt = D[fit], D[tgt]
        dbar = Dn.mean(0)
        dhat = unit(dbar)
        # V1: 全体ドリフト方向の成分を除いてから射影
        Dt1 = Dt - np.outer(Dt @ dhat, dhat)
        cols["V1"][tgt] = Dt1 @ u
        # V2: 全体ドリフトベクトル（平均変位）を引いてから射影
        cols["V2"][tgt] = (Dt - dbar) @ u
        # V3: 局所ドリフト（元埋め込みが近い null 細胞 K 個の平均変位）を引いてから射影
        sim = e0[tgt] @ e0[fit].T
        nn = np.argpartition(-sim, K_NN, axis=1)[:, :K_NN]
        cols["V3"][tgt] = np.einsum("ij,j->i", Dt - Dn[nn].mean(1), u)
        # V4: null 変位の上位 N_PC 主成分（ドリフト部分空間, 平均を含めず中心化前の二次モーメント）を除いて射影
        _, _, vt = np.linalg.svd(Dn, full_matrices=False)
        P = vt[:N_PC]
        Dt4 = Dt - (Dt @ P.T) @ P
        u4 = u - P.T @ (P @ u)
        cols["V4"][tgt] = Dt4 @ unit(u4)
    for k, v in cols.items():
        cells[k] = v
    # --- 事後（post hoc, V0〜V4 の結果を見た後に追加。選定には使わず、本番 run で確認する）---
    L = cells.length.to_numpy()
    v2l = np.full(len(cells), np.nan)
    for tgt, fit in groups:  # V2L: 長さ 5 分位ビンごとの null 平均変位を引いてから射影（a2 のベクトル版）
        edges = np.quantile(L[fit], [0.2, 0.4, 0.6, 0.8])
        bf, bt = np.digitize(L[fit], edges), np.digitize(L[tgt], edges)
        mu = np.vstack([D[fit][bf == b].mean(0) for b in range(5)])
        v2l[tgt] = (D[tgt] - mu[bt]) @ u
    cells["V2L"] = v2l
    Dt = D - (D * e0).sum(1, keepdims=True) * e0     # 接平面への射影（|Δ|² による一様な減少を除く）
    cells["T_u"] = Dt @ u
    cells["V5"] = (D @ u) / np.linalg.norm(D, axis=1)  # 向きだけ（大きさを除く）
    # 記述: ドリフトの向き
    dn = D[isnull]
    dhat = unit(dn.mean(0))
    desc = {"null_drift_norm_mean": float(np.linalg.norm(dn.mean(0))),
            "null_disp_norm_mean": float(np.linalg.norm(dn, axis=1).mean()),
            "frac_disp_along_global_drift_mean": float(np.mean(((dn @ dhat) ** 2) / (np.linalg.norm(dn, axis=1) ** 2))),
            "cos_drift_vs_u": float(dhat @ u),
            "cos_drift_vs_minus_mean_orig": float(dhat @ -unit(e0.mean(0))),
            "cos_drift_vs_goal": float(dhat @ g), "cos_drift_vs_start": float(dhat @ s)}
    rad = (D * e0).sum(1)[isnull]                   # ê0 方向成分 = cos(e0,e1) − 1 ≈ −|Δ|²/2
    desc["mean_shift_goal_null"] = float((D[isnull] @ g).mean())
    desc["radial_part_of_shift_goal_null"] = float((rad * (e0[isnull] @ g)).mean())
    sv = np.linalg.svd(dn - 0, compute_uv=False)
    desc["null_disp_top_pc_var_frac"] = [round(float(x), 3) for x in (sv[:5] ** 2 / (sv ** 2).sum())]
    return cells, desc


LABELS = {"V0": "V0: Δ·u（補正なし）", "V1": "V1: 全体ドリフト方向を除去", "V2": "V2: 平均ドリフトを減算",
          "V3": f"V3: 局所ドリフト（kNN{K_NN}）を減算", "V4": f"V4: ドリフト部分空間（PC{N_PC}）を除去",
          "V2L": "[事後] V2L: 長さビン別の平均ドリフトを減算", "T_u": "[事後] T_u: 接平面射影", "V5": "[事後] V5: cos(Δ,u)"}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cells, E0, E1, cen = load()
    cons = consistency(cells, E0, E1, cen)
    print("consistency:", json.dumps(cons, indent=1))
    cells, desc = vector_metrics(cells, E0, E1, cen)
    print("drift:", json.dumps(desc, indent=1))
    alt = [x for x in cells.columns if x.startswith("shift_alt")][0]
    cells["alt"] = cells[alt]
    cells["pos"] = rc.token_positions(cells)
    cells = rc.corrected(cells)
    det = json.loads((ROOT / "experiments/EXP-20260926-001/work/detection-aPT_CKD_TLShigh.json").read_text())
    rows = []
    for m, lab in list(rc.METRICS.items()) + list(LABELS.items()):
        res, _ = rc.evaluate(cells, det, m)
        rows.append({"label": lab, **res})
    df = pd.DataFrame(rows)
    for col in ["C1_seeds", "C1_pct_mean", "C3_seeds", "C3_pct_mean"]:
        df[col] = df[col].map(lambda d: " ".join(f"{k}:{v}" for k, v in d.items()))
    df.to_csv(OUT / "criteria.csv", index=False)
    (OUT / "criteria.json").write_text(json.dumps(rows, ensure_ascii=False, indent=1, default=str) + "\n")
    (OUT / "consistency.json").write_text(json.dumps({"consistency": cons, "drift": desc}, indent=1) + "\n")
    keep = ["run", "gene", "role", "cell_uid", "donor_id", "length", "pos"] + list(rc.METRICS) + list(LABELS)
    cells[keep].to_csv(OUT / "cells_vector.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print(df[["label", "C1_pct_mean", "C2_wins", "C2_margin_nullSD", "C3_pct_mean", "C4i_splithalf_rho",
              "C5_rho_len_gene", "C5_rho_det_gene", "C5_rho_len_cell_median", "null_frac_negative",
              "C1", "C2", "C3", "C4", "C5", "n_pass_all"]].to_string(index=False))


if __name__ == "__main__":
    main()

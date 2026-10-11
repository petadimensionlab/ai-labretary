#!/usr/bin/env python3
"""推す案（向きから共通成分を抜く）の裏取り。

(1) 仕組みの確認: 移動量が大きい遺伝子ほど「誰もが動く向き」に沿って動くか
(2) 共通の向きの作り方を変えても結果が変わらないか（偽薬の向きを含む）
(3) 主成分を 2〜10 本抜いたらどうなるか
(4) 提供者を 1 人ずつ除いたときの M5 と陽性対照
(5) 対照 500 個の引き直しに対して M5 がどれだけ揺れるか
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import diag  # noqa: E402
import run_variants as rv  # noqa: E402

HERE = Path(__file__).parent
NULLSET = set(diag.NULLS_B500)


def tabs_for(cells, elig, col, xcols=("lp", "ll"), form="full_quadratic", genes=None):
    tg, tu = {}, {}
    for s in diag.SEEDS:
        cs = cells[cells.seed == s]
        nl = {h: hd for h, hd in cs[cs.gene.isin(NULLSET)].groupby("gene")}
        tg[s] = diag.gene_stats(cs, nl, elig[s], col, list(xcols), form=form)
        tu[s] = diag.gene_stats(cs, nl, elig[s], "m_uniform", list(xcols), form=form)
    return tg, tu


def m5_p4(cells, elig, col):
    tg, tu = tabs_for(cells, elig, col)
    M = diag.evaluate(tg, tu, cells, n_perm=50000)
    return {"M5": [round(M["M5"][s], 4) for s in diag.SEEDS],
            "M2max": [round(M["M2_max_abs"][s], 3) for s in diag.SEEDS],
            "M4": [round(M["M4"][s], 3) for s in diag.SEEDS],
            "M3": M["M3_n_in_range"], "陽性対照71": M["positive_control71"]["z"],
            "合格": M["all_pass"]}


def main():
    cells = pd.read_parquet(diag.CACHE / "cells2.parquet")
    disp = diag.load_disp()
    elig = diag.get_elig(cells)
    states = json.loads((diag.CACHE / "states.json").read_text())
    out = {}

    # ---------------- (1) 仕組み
    mech = {}
    for s in diag.SEEDS:
        cs = cells[(cells.seed == s) & cells.gene.isin(NULLSET)]
        g = cs.groupby("gene").agg(dn=("dnorm", "median"), fc=("f_common", "mean"),
                                   vb=("v_base", "mean"), vu=("v_uperp", "mean"))
        r = lambda a, b: round(float(stats.spearmanr(a, b).statistic), 3)
        mech[s] = {"共通成分の割合の平均": round(float(g.fc.mean()), 4),
                   "rho(移動量, 共通成分の割合)": r(g.dn, g.fc),
                   "rho(移動量, 事前登録の指標)": r(g.dn, g.vb),
                   "rho(移動量, 共通成分を抜いた指標)": r(g.dn, g.vu),
                   "細胞ごとの相関 corr(移動量の対数, 共通成分の割合)":
                       round(float(np.corrcoef(np.log2(cs.dnorm), cs.f_common)[0, 1]), 3)}
    out["1_仕組み"] = mech
    out["1_û·ĉ"] = json.loads((HERE / "u_dot_c.json").read_text())

    # ---------------- (2)(3) 共通の向きの作り方
    info_unit = rv.build_common(cells, disp)
    U = {}
    for s in diag.SEEDS:
        st = {k: np.asarray(v) for k, v in states[f"run-{s}-L-1"].items()}
        vg, vs = st["goal"], st["start"]
        u = vg / np.linalg.norm(vg) - vs / np.linalg.norm(vs)
        U[s] = u / np.linalg.norm(u)

    # 変位そのものの平均、対照の半分、主成分、偽薬
    dirs = {}
    rng = np.random.default_rng(7)
    for s in diag.SEEDS:
        m = ((cells.seed == s) & cells.gene.isin(NULLSET)).to_numpy()
        idx = np.flatnonzero(m)
        gl = cells.gene.to_numpy()[idx]
        half = set(rng.choice(sorted(NULLSET), 250, replace=False))
        hm = np.isin(gl, list(half))
        craw = np.zeros(1152)
        chalf = np.zeros(1152)
        G = np.zeros((1152, 1152))
        for a in range(0, len(idx), 20000):
            b = min(a + 20000, len(idx))
            D = disp[idx[a:b]].astype(np.float64)
            craw += D.sum(axis=0)
            Un = D / np.linalg.norm(D, axis=1, keepdims=True)
            chalf += Un[hm[a:b]].sum(axis=0)
            G += Un.T @ Un
            del D, Un
        craw /= np.linalg.norm(craw)
        chalf /= np.linalg.norm(chalf)
        w, V = np.linalg.eigh(G)
        V = V[:, ::-1]
        rnd = rng.normal(size=1152)
        rnd /= np.linalg.norm(rnd)
        dirs[s] = {"raw": craw, "half": chalf, "pc": V, "rand": rnd,
                   "寄与率上位5": [round(float(x / w.sum()), 4) for x in w[::-1][:5]]}
        del G

    def proj_metric(name, getdir):
        v = np.full(len(cells), np.nan)
        for s in diag.SEEDS:
            u = U[s]
            B = getdir(s)                      # 列が抜く向き（正規直交）
            uo = u - B @ (B.T @ u)
            uo = uo / np.linalg.norm(uo)
            ii = np.flatnonzero((cells.seed == s).to_numpy())
            for a in range(0, len(ii), 20000):
                jj = ii[a:a + 20000]
                v[jj] = (disp[jj].astype(np.float64) @ uo) / cells.dnorm.to_numpy()[jj]
        cells[name] = v
        return m5_p4(cells, elig, name)

    alt = {}
    alt["共通の向き＝単位変位の平均（推す案。遺伝子を除く）"] = m5_p4(cells, elig, "v_uperp")
    alt["共通の向き＝変位そのものの平均"] = proj_metric(
        "v_raw", lambda s: dirs[s]["raw"][:, None])
    alt["共通の向き＝対照の半分だけで作る"] = proj_metric(
        "v_half", lambda s: dirs[s]["half"][:, None])
    for k in (2, 3, 5, 10):
        alt[f"主成分を {k} 本抜く"] = proj_metric(
            f"v_pc{k}", lambda s, k=k: dirs[s]["pc"][:, :k])
    alt["偽薬: 無作為な向きを 1 本抜く"] = proj_metric(
        "v_rand", lambda s: dirs[s]["rand"][:, None])
    out["2_共通の向きの作り方"] = alt
    out["2_主成分の寄与率"] = {s: dirs[s]["寄与率上位5"] for s in diag.SEEDS}

    # ---------------- (4) 提供者を 1 人ずつ除く
    loo = {}
    for col, label in (("v_base", "事前登録"), ("v_uperp", "推す案")):
        rows = []
        donors = sorted(cells.donor_id.unique())
        for d in donors:
            sub = cells[cells.donor_id != d]
            tg, _ = tabs_for(sub, elig, col)
            m5 = [float(stats.spearmanr(tg[s][tg[s].gene.isin(NULLSET)].dnorm_med,
                                        tg[s][tg[s].gene.isin(NULLSET)].D_raw).statistic)
                  for s in diag.SEEDS]
            rows.append(m5)
        A = np.array(rows)
        loo[label] = {"提供者数": len(donors),
                      "M5 の最小": [round(float(x), 4) for x in A.min(axis=0)],
                      "M5 の最大": [round(float(x), 4) for x in A.max(axis=0)],
                      "閾値 0.10 を外れた回数": [int((np.abs(A[:, i]) > 0.10).sum())
                                                 for i in range(3)]}
        print(label, loo[label], flush=True)
    out["4_提供者を1人ずつ除く"] = loo

    # ---------------- (5) 対照プールの引き直しに対する揺れ
    boot = {}
    for col, label in (("v_base", "事前登録"), ("v_uperp", "推す案")):
        tg, _ = tabs_for(cells, elig, col)
        vals = {s: tg[s][tg[s].gene.isin(NULLSET)][["gene", "dnorm_med", "D_raw"]]
                for s in diag.SEEDS}
        rg = np.random.default_rng(11)
        B = []
        for _ in range(2000):
            ii = rg.integers(0, 500, 500)
            B.append([float(stats.spearmanr(vals[s].dnorm_med.to_numpy()[ii],
                                            vals[s].D_raw.to_numpy()[ii]).statistic)
                      for s in diag.SEEDS])
        A = np.array(B)
        boot[label] = {"M5 の標準誤差（対照 500 個の引き直し）": [round(float(x), 4)
                                                                  for x in A.std(axis=0, ddof=1)],
                       "95% 区間": [[round(float(np.quantile(A[:, i], 0.025)), 4),
                                     round(float(np.quantile(A[:, i], 0.975)), 4)]
                                    for i in range(3)],
                       "|M5|>0.10 になる割合": [round(float((np.abs(A[:, i]) > 0.10).mean()), 3)
                                                for i in range(3)]}
    out["5_対照プールの揺れ"] = boot

    (HERE / "robust.json").write_text(json.dumps(out, ensure_ascii=False, indent=1,
                                                 default=float) + "\n")
    print(json.dumps(out, ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    main()

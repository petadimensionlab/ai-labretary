#!/usr/bin/env python3
"""tls_feasibility.py の補足: 仮説設計に直結する細胞数・関連を計算する（読み取り専用）。

出力（outputs/literature/2026-09-25-tls-feasibility/）
  special_cells.csv     推定サブセット（SA-T 様、ABC 様、CCL19/CCL21/CXCL13+ 線維芽細胞、TLS 群別 aPT）の細胞数
  donor_assoc_extra.csv donor 単位: aPT 比率・SA-T 様比率・ABC 様比率 と TLS スコア/年齢/疾患の関連
  start_state_genes.csv 主要 start state における候補遺伝子の検出率（donor 平均）を TLS-high vs TLS-low で比較
  cxcl13_source.csv     CXCL13 カウントの細胞種別内訳（どの細胞が CXCL13 を出しているか）
  tumour_cd8_cxcl13.csv 腫瘍サンプルごとの CXCL13+ CD8 T 比率と B+形質細胞比率（A3/A4）
  hyp_state_counts.csv  仮説で使う start/goal/alt 状態の細胞数（TLS 群 × 疾患/領域）と発現細胞数
usage: .venv/bin/python scripts/analysis/tls_feasibility_extra.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tls_feasibility as tf  # noqa: E402

OUT = tf.OUT
# 推定サブセットの定義（すべて「推定」ラベル。著者アノテーションではない）
T_LABELS = {"A1": ["T", "NKT"], "A2": ["T", "T-CYT", "NKT", "cycT"],
            "A6": ["Naïve Th Cell", "Effector Th Cell", "Effector Tc Cell", "Naïve Tc Cell", "Treg", "NKT Cell", "MAIT"],
            "A3": None, "A4": None, "A5": ["Tcell", "Tcell_CD8"]}
B_LABELS = {"A1": ["B"], "A2": ["B"], "A6": ["Memory B Cell", "Naïve B Cell"],
            "A3": ["Naïve_B.1", "Naïve_B.2", "Non-switched-memory_B", "Switched-memory_B", "Activated_B-AREG",
                   "Activated_B-RHOB", "Stress_B", "INF-response_B"], "A4": ["B cell"], "A5": ["Bcell"]}
FIB_LABELS = {"A1": ["FIB", "aFIB", "dFIB", "M-FIB", "dM-FIB", "MYOF", "cycMYOF"], "A2": ["FIB", "aFIB", "dFIB", "MyoF"],
              "A6": ["Fibroblast", "Myofibroblast"], "A3": None, "A4": ["CD45- Myofibroblast"], "A5": ["Peri", "vSMC"]}
APT = {"A1": (["aPT"], ["PT-S1", "PT-S2", "PT-S3"]), "A2": (["aPT"], ["PT-S1/S2", "PT-S3"]),
       "A6": (["PT Injured", "PT Inflammatory"], ["PT S1", "PT S2", "PT S3"])}
CAND = ["VCAM1", "HAVCR1", "LCN2", "SPP1", "CDKN1A", "CDKN2A", "CCL2", "CXCL10", "CXCL9", "LTB", "LTBR", "RELB",
        "NFKB2", "TNFSF13B", "TNFSF13", "CD40", "IL7R", "CCL19", "CCL21", "CXCL13", "CXCL12", "IGFBP7", "GDF15",
        "SERPINE1", "IL6", "CXCL8", "TP53", "GLB1", "LMNB1"]


def load(key):
    d = tf.DS[key]
    samp_cols = [d["sample"]] if isinstance(d["sample"], str) else list(d["sample"])
    need = list(dict.fromkeys(samp_cols + [d["state"], d["group"], d["age"], "cell_type", "disease", "donor_id"]
                              + (["broad_type"] if key == "A3" else [])))
    parts = [tf.extract(f, need) for f in d["files"]]
    obs = pd.concat([p[0] for p in parts], ignore_index=True)
    G = pd.concat([p[1] for p in parts], ignore_index=True).fillna(0)
    obs["sample"] = obs[samp_cols].astype(str).agg("|".join, axis=1)
    donor = pd.read_csv(OUT / f"donor_{key}.csv", index_col=0)
    obs["TLS_group"] = obs["sample"].map(donor["TLS_group"])
    obs["TLS_score"] = obs["sample"].map(donor["TLS_score"])
    return d, obs, G, donor


def main():
    sc_rows, assoc_rows, gene_rows, src_rows = [], [], [], []
    for key in ["A1", "A2", "A3", "A4", "A5", "A6"]:
        d, obs, G, donor = load(key)
        st = obs[d["state"]].astype(str)
        bt = obs["broad_type"].astype(str) if key == "A3" else None
        isT = st.isin(T_LABELS[key]) if T_LABELS[key] else bt.eq("T-cell") if key == "A3" else obs.cell_type.str.contains("T cell|T-cell", regex=True)
        isB = st.isin(B_LABELS[key])
        isF = st.isin(FIB_LABELS[key]) if FIB_LABELS[key] else bt.eq("Fibro")
        sat = isT & (G["TNFSF8"] > 0) & (G["PDCD1"] > 0)             # SA-T 様（CD153+PD-1+）
        abc = isB & (G["ITGAX"] > 0) & ((G["TBX21"] > 0) | (G["FCRL5"] > 0))  # ABC 様
        fl = isF & ((G["CCL19"] > 0) | (G["CCL21"] > 0) | (G["CXCL13"] > 0))  # TLS オーガナイザー様線維芽細胞
        subsets = {"T_all": isT, "SA-T-like(TNFSF8+PDCD1+ T)": sat, "B_all": isB, "ABC-like(ITGAX+ & TBX21+|FCRL5+ B)": abc,
                   "FIB_all": isF, "FIB_CCL19|CCL21|CXCL13+": fl, "FIB_CCL19+": isF & (G["CCL19"] > 0), "FIB_CCL21+": isF & (G["CCL21"] > 0)}
        if key in APT:
            subsets["startPT_injured(" + "+".join(APT[key][0]) + ")"] = st.isin(APT[key][0])
            subsets["goalPT_healthy(" + "+".join(APT[key][1]) + ")"] = st.isin(APT[key][1])
        for name, m in subsets.items():
            row = dict(dataset=key, subset=name, n_cells=int(m.sum()), n_samples=int(obs.loc[m, "sample"].nunique()),
                       n_TLS_high=int((m & obs.TLS_group.eq("high")).sum()), n_TLS_mid=int((m & obs.TLS_group.eq("mid")).sum()),
                       n_TLS_low=int((m & obs.TLS_group.eq("low")).sum()))
            for gl in obs[d["group"]].astype(str).unique():
                row[f"n_{gl}"] = int((m & obs[d["group"]].astype(str).eq(gl)).sum())
            sc_rows.append(row)
        # donor 単位の比率
        g = obs.assign(isT=isT, sat=sat, isB=isB, abc=abc, isF=isF, fl=fl).groupby("sample")
        dd = donor.copy()
        dd["frac_SAT_in_T"] = g["sat"].sum() / g["isT"].sum().replace(0, np.nan)
        dd["frac_ABC_in_B"] = g["abc"].sum() / g["isB"].sum().replace(0, np.nan)
        dd["frac_orgFIB_in_FIB"] = g["fl"].sum() / g["isF"].sum().replace(0, np.nan)
        if key in APT:
            inj, hea = st.isin(APT[key][0]), st.isin(APT[key][1])
            gg = obs.assign(inj=inj, hea=hea).groupby("sample")
            dd["frac_injPT_in_PT"] = gg["inj"].sum() / (gg["inj"].sum() + gg["hea"].sum()).replace(0, np.nan)
        ref = dd.in_score_reference.astype(bool)
        dr = dd[ref]
        for col, minT in [("frac_SAT_in_T", "isT"), ("frac_ABC_in_B", "isB"), ("frac_orgFIB_in_FIB", "isF"), ("frac_injPT_in_PT", None)]:
            if col not in dr:
                continue
            if minT:
                cnt = g[minT].sum()
                keep = dr.index[cnt.reindex(dr.index).fillna(0) >= 20]
            else:
                keep = dr.index
            x = dr.loc[keep]
            x = x[x[col].notna()]
            if len(x) < 6:
                assoc_rows.append(dict(dataset=key, variable=col, vs="TLS_score", n=len(x), rho=np.nan, p=np.nan, note="n<6 (>=20 parent cells)"))
                continue
            r = stats.spearmanr(x[col], x.TLS_score)
            assoc_rows.append(dict(dataset=key, variable=col, vs="TLS_score", n=len(x), rho=r.statistic, p=r.pvalue,
                                   note=f"median={x[col].median():.4f}"))
            if x.age.nunique() > 2:
                r = stats.spearmanr(x[col], x.age)
                assoc_rows.append(dict(dataset=key, variable=col, vs="age", n=len(x), rho=r.statistic, p=r.pvalue, note=""))
            lv = x.group.astype(str)
            if "normal" in set(lv) and lv.nunique() > 1:
                dis = x[lv != "normal"][col]
                nor = x[lv == "normal"][col]
                r = stats.mannwhitneyu(dis, nor)
                assoc_rows.append(dict(dataset=key, variable=col, vs="disease vs normal", n=len(x), rho=np.nan, p=r.pvalue,
                                       note=f"disease median={dis.median():.4f} (n={len(dis)}), normal median={nor.median():.4f} (n={len(nor)})"))
        # 候補遺伝子の検出率（start state 内、donor 平均; ≥20 細胞の donor のみ）
        if key in APT:
            m = st.isin(APT[key][0])
            sub = obs[m]
            pos = (G.loc[m, [c for c in CAND if c in G]] > 0)
            pos["sample"] = sub["sample"].values
            pos["grp"] = sub["TLS_group"].values
            cnt = pos.groupby("sample").size()
            per = pos.groupby("sample").mean(numeric_only=True)
            per = per[cnt.reindex(per.index) >= 20]
            grp = pos.groupby("sample")["grp"].first().reindex(per.index)
            for gname in per.columns:
                hi, lo = per.loc[grp == "high", gname], per.loc[grp == "low", gname]
                p = stats.mannwhitneyu(hi, lo).pvalue if len(hi) >= 3 and len(lo) >= 3 else np.nan
                gene_rows.append(dict(dataset=key, start_state="+".join(APT[key][0]), gene=gname,
                                      det_all=float((G.loc[m, gname] > 0).mean()),
                                      det_TLShigh_donor_mean=hi.mean(), n_high=len(hi), det_TLSlow_donor_mean=lo.mean(), n_low=len(lo), p_mwu=p))
        # CXCL13 の供給源
        tot = G["CXCL13"].sum()
        by = G["CXCL13"].groupby(st.values).sum().sort_values(ascending=False)
        for lab, v in by.head(6).items():
            src_rows.append(dict(dataset=key, state=lab, CXCL13_count_share=v / tot if tot else np.nan, total_CXCL13_counts=tot))
        print(f"[{key}] done", flush=True)
    pd.DataFrame(sc_rows).to_csv(OUT / "special_cells.csv", index=False)
    pd.DataFrame(assoc_rows).to_csv(OUT / "donor_assoc_extra.csv", index=False)
    pd.DataFrame(gene_rows).to_csv(OUT / "start_state_genes.csv", index=False)
    pd.DataFrame(src_rows).to_csv(OUT / "cxcl13_source.csv", index=False)


HYP_STATES = {
    # (dataset, 状態名, ラベル列の値, 追加条件) — 追加条件: None | "org" | "neg"（線維芽細胞の CCL19/CCL21/CXCL13 検出の有無）
    "A2": [("aPT", ["aPT"], None), ("PT_healthy", ["PT-S1/S2", "PT-S3"], None)],
    "A1": [("FIB_org", FIB_LABELS["A1"], "org"), ("FIB_neg", FIB_LABELS["A1"], "neg")],
    "A3": [("Activated_B", ["Activated_B-AREG", "Activated_B-RHOB"], None), ("Plasma-IgG", ["Plasma-IgG"], None),
           ("Plasma-IgA", ["Plasma-IgA"], None), ("memory_B", ["Non-switched-memory_B", "Switched-memory_B"], None)],
}
HYP_GENES = {"A2": ["TP53", "CDKN1A", "LTBR", "SPP1", "LCN2", "VCAM1", "HAVCR1"],
             "A1": ["LTBR", "RELB", "NFKB2", "CDKN1A", "IL6", "NGFR", "ALDH1A2"],
             "A3": ["TBX21", "ITGAX", "ZEB2", "FCRL5", "CDKN1A", "PRDM1", "XBP1"]}


def extra_tables():
    rows, trows = [], []
    for key, defs in HYP_STATES.items():
        d, obs, G, donor = load(key)
        st = obs[d["state"]].astype(str)
        org = (G["CCL19"] > 0) | (G["CCL21"] > 0) | (G["CXCL13"] > 0)
        grpv = obs[d["group"]].astype(str)
        for name, labs, cond in defs:
            m = st.isin(labs)
            if cond == "org":
                m &= org
            elif cond == "neg":
                m &= ~org
            for (tg, gl), sub in obs[m].groupby([obs.loc[m, "TLS_group"].fillna("NA"), grpv[m]]):
                r = dict(dataset=key, state=name, TLS_group=tg, group=gl, n_cells=len(sub), n_samples=sub["sample"].nunique())
                for g in HYP_GENES[key]:
                    if g in G:
                        r[f"n_expr_{g}"] = int((G.loc[sub.index, g] > 0).sum())
                rows.append(r)
    for key in ["A3", "A4"]:
        d, obs, G, donor = load(key)
        st = obs[d["state"]].astype(str)
        if key == "A3":
            bt = obs["broad_type"].astype(str)
            cd8, isB, isPL, ly = st.str.startswith("CD8+T"), bt.eq("B-cell"), bt.eq("Plasma"), bt.isin(["T-cell", "NK", "B-cell", "Plasma"])
            tum = obs[d["group"]].eq("Tumour")
        else:
            cd8, isB, isPL = st.str.startswith("CD8A+"), st.eq("B cell"), pd.Series(False, index=st.index)
            ly = obs.cell_type.str.contains("T cell|NK|B cell")
            tum = obs[d["group"]].eq("Tumor")
        x = pd.DataFrame({"cd8": cd8, "cx": cd8 & (G["CXCL13"] > 0), "bpl": isB | isPL, "pl": isPL, "ly": ly,
                          "s": obs["sample"], "donor": obs["donor_id"]})[tum]
        g = x.groupby("s")
        t = pd.DataFrame({"dataset": key, "donor": g.donor.first(), "n_CD8": g.cd8.sum(), "frac_CXCL13pos_in_CD8": g.cx.sum() / g.cd8.sum(),
                          "frac_BPL_in_lymph": g.bpl.sum() / g.ly.sum(), "frac_PL_in_lymph": g.pl.sum() / g.ly.sum()})
        trows.append(t[t.n_CD8 >= 50])
    pd.DataFrame(rows).to_csv(OUT / "hyp_state_counts.csv", index=False)
    pd.concat(trows).to_csv(OUT / "tumour_cd8_cxcl13.csv")


if __name__ == "__main__":
    if "--tables-only" not in sys.argv:
        main()
    extra_tables()

#!/usr/bin/env python3
"""腎 TLS × 老化 仮説のための feasibility 解析（読み取り専用）。

data/raw/kidney_tls/*.h5ad の raw.X（整数カウント）から、必要な遺伝子列と細胞ごとの総カウント・
発現遺伝子数だけを h5py でチャンク読みする（元ファイルは一切変更しない）。

出力: outputs/literature/2026-09-25-tls-feasibility/
  donor_<DS>.csv      サンプル（donor / region）単位の B+形質細胞比率、pseudo-bulk CPM、TLS スコア、TLS 群
  states_<DS>.csv     細胞状態ごとの細胞数（全体 / TLS-high / TLS-low donor）、老化・SASP スコア、陽性率、中央値遺伝子数
  assoc_summary.csv   TLS スコアと疾患・年齢・腫瘍/正常の関連（Spearman / Mann-Whitney / Kruskal）
  token_dict.csv      対象遺伝子の Geneformer トークン辞書への収載とデータセット var への収載
usage: .venv/bin/python scripts/analysis/tls_feasibility.py [DS ...]   (DS: A1 A2 A3 A4 A5 A6)
"""
import pickle, re, sys, warnings
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

warnings.filterwarnings("ignore")
sys.path = [p for p in sys.path if not p.endswith("/scripts") and not p.endswith("/scripts/analysis")]
ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw/kidney_tls"
OUT = ROOT / "outputs/literature/2026-09-25-tls-feasibility"
OUT.mkdir(parents=True, exist_ok=True)

GENESETS = {
    "tls_core": ["CXCL13", "CCL19", "CCL21", "LTB", "LTA", "CXCR5", "CCR7", "MS4A1"],   # TLS スコアに使う 8 遺伝子
    "tls_extra": ["CD79A", "CD79B", "CR2", "LAMP3", "CD19", "AICDA", "BCL6"],
    "chemokine12": ["CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CXCL9", "CXCL10", "CXCL11"],
    "sen_core": ["CDKN1A", "CDKN2A", "TP53", "LMNB1", "GLB1"],
    "sasp": ["IL6", "CXCL8", "IL1A", "IL1B", "SERPINE1", "MMP1", "MMP3", "IGFBP3", "IGFBP7", "GDF15"],
    "pt_injury": ["HAVCR1", "LCN2", "VCAM1", "SPP1"],
    "tlt_fib": ["NGFR", "PDPN", "ALDH1A2", "CXCL12"],
    "sa_t": ["TNFSF8", "TNFRSF8", "PDCD1"],
    "abc": ["ITGAX", "TBX21", "FCRL5", "ZEB2"],
    "plasma": ["JCHAIN", "MZB1", "XBP1", "PRDM1", "IGHG1", "IGHA1"],
    "lt_axis": ["LTBR", "RELB", "NFKB2", "TNFSF13B", "TNFSF13", "CD40", "CD40LG", "TNF", "IL7R"],
    "housekeeping": ["ACTB", "GAPDH", "RPLP0", "PPIA", "TBP"],
    "lineage": ["PTPRC", "CD3E"],
}
ALL_GENES = list(dict.fromkeys(g for v in GENESETS.values() for g in v))
TLS8 = GENESETS["tls_core"]
SASP = GENESETS["sasp"]
POS_FRAC_GENES = ["CXCL13", "CCL19", "CCL21", "LTB", "CDKN1A", "CDKN2A", "VCAM1", "HAVCR1", "TNFSF8",
                  "TNFRSF8", "ITGAX", "TBX21", "FCRL5", "NGFR", "PDPN", "LTBR", "TNFSF13B"]
TOK_PER_SEC = 7200.0  # MPS bf16 batch 8 の実測（PI 周知値）

# データセット定義: file(s), sample 列, state 列, 群（disease/tumour）列, 年齢の取り方
DS = {
    "A1": dict(files=["A1_KPMP_v1.5_snRNA_Lake2023"], sample="donor_id", state="subclass.l2", group="disease", age="Age_binned"),
    "A2": dict(files=["A2_KPMP_v1.5_scRNA_Lake2023"], sample="donor_id", state="author_cell_type", group="disease", age="Age_binned"),
    "A3": dict(files=["A3_RCC_LiR2022"], sample=("donor_id", "summaryDescription"), state="annotation", group="summaryDescription", age="development_stage",
               ref_groups=["Tumour", "Tumour-normal", "Normal kidney"]),
    "A4": dict(files=["A4_ccRCC_Krishna2021"], sample="author_sample", state="author_cluster_name", group="author_type", age="development_stage",
               ref_groups=["Tumor", "Normal"]),
    "A5": dict(files=["A5a_ccRCC_Zhang2021", "A5b_normal_Zhang2021"], sample=("donor_id", "disease"), state="author_cell_type", group="disease", age="development_stage"),
    "A6": dict(files=["A6_UUO_Reck2025"], sample="donor_id", state="celltype_l2", group="disease", age="development_stage"),
}
DECADE = {"second": 15, "third": 25, "fourth": 35, "fifth": 45, "sixth": 55, "seventh": 65, "eighth": 75, "nineth": 85, "ninth": 85}


def age_num(s: str) -> float:
    s = str(s)
    m = re.match(r"(\d+)-year-old", s)
    if m:
        return float(m.group(1))
    for k, v in DECADE.items():
        if s.startswith(k):
            return float(v)
    return np.nan


def read_obs(h: h5py.File, cols) -> pd.DataFrame:
    out = {}
    for c in cols:
        g = h["obs"][c]
        if isinstance(g, h5py.Group):  # categorical
            cats = np.array(g["categories"].asstr()[...]) if g["categories"].dtype.kind in "OS" else g["categories"][...].astype(str)
            codes = g["codes"][...]
            out[c] = np.where(codes >= 0, cats[np.clip(codes, 0, None)], "nan")
        else:
            v = g[...]
            out[c] = np.array(g.asstr()[...]) if v.dtype.kind in "OS" else v
    return pd.DataFrame(out)


def read_var_names(h: h5py.File):
    v = h["raw/var"]
    ens = np.array(v["_index"].asstr()[...])
    fn = v["feature_name"]
    if isinstance(fn, h5py.Group):
        cats = np.array(fn["categories"].asstr()[...])
        sym = cats[fn["codes"][...]]
    else:
        sym = np.array(fn.asstr()[...])
    return ens, sym


def extract(fname: str, obs_cols):
    """raw.X から対象遺伝子のカウント、総カウント、発現遺伝子数を取り出す。"""
    with h5py.File(RAW / f"{fname}.h5ad", "r") as h:
        obs = read_obs(h, obs_cols)
        ens, sym = read_var_names(h)
        sym2idx = {}
        for i, s in enumerate(sym):
            sym2idx.setdefault(s, i)
        present = [g for g in ALL_GENES if g in sym2idx]
        cols = np.array([sym2idx[g] for g in present])
        X = h["raw/X"]
        indptr = X["indptr"][...]
        n = len(indptr) - 1
        total = np.zeros(n, dtype=np.float64)
        ngenes = np.zeros(n, dtype=np.int32)
        blocks = []
        step = 20000
        for s in range(0, n, step):
            e = min(n, s + step)
            a, b = indptr[s], indptr[e]
            m = sp.csr_matrix((X["data"][a:b], X["indices"][a:b], indptr[s:e + 1] - a), shape=(e - s, len(ens)))
            total[s:e] = np.asarray(m.sum(1)).ravel()
            ngenes[s:e] = np.diff(m.indptr)
            blocks.append(m[:, cols].toarray().astype(np.float32))
        G = pd.DataFrame(np.vstack(blocks), columns=present)
        ens_map = {sym[sym2idx[g]]: ens[sym2idx[g]] for g in present}
    return obs, G, total, ngenes, ens_map


def is_B(ct): return bool(re.search(r"\bB cell|naive B|memory B", ct, re.I)) and "plasma" not in ct.lower()
def is_PL(ct): return "plasma" in ct.lower()
def is_lymph(ct): return bool(re.search(r"T cell|T-cell|NK|natural killer|B cell|plasma|lymphoid|MAIT|B-cell", ct, re.I))


def analyse(key: str):
    d = DS[key]
    samp_cols = [d["sample"]] if isinstance(d["sample"], str) else list(d["sample"])
    need = list(dict.fromkeys(samp_cols + [d["state"], d["group"], d["age"], "cell_type", "disease", "donor_id"]
                              + (["broad_type"] if key == "A3" else [])))
    parts = []
    for f in d["files"]:
        obs, G, total, ngenes, ens_map = extract(f, need)
        obs["file"] = f
        parts.append((obs, G, total, ngenes))
    obs = pd.concat([p[0] for p in parts], ignore_index=True)
    G = pd.concat([p[1] for p in parts], ignore_index=True).fillna(0)
    total = np.concatenate([p[2] for p in parts])
    ngenes = np.concatenate([p[3] for p in parts])
    obs["sample"] = obs[samp_cols].astype(str).agg("|".join, axis=1)
    obs["age_num"] = obs[d["age"]].map(age_num)
    # B / 形質細胞 / リンパ球の判定
    ct = obs["cell_type"].astype(str)
    if key == "A3":
        bt = obs["broad_type"].astype(str)
        obs["isB"], obs["isPL"] = bt.eq("B-cell"), bt.eq("Plasma")
        obs["isLy"] = bt.isin(["T-cell", "NK", "B-cell", "Plasma"])
    else:
        obs["isB"], obs["isPL"], obs["isLy"] = ct.map(is_B), ct.map(is_PL), ct.map(is_lymph)
    # 正規化（per-cell log1p CP10k）→ 老化・SASP スコア
    ln = np.log1p(G.values / np.maximum(total, 1)[:, None] * 1e4)
    LN = pd.DataFrame(ln, columns=G.columns)
    sen_g = [g for g in ["CDKN1A", "CDKN2A"] if g in LN]
    sasp_g = [g for g in SASP if g in LN]
    obs["sen_score"] = LN[sen_g].mean(1)
    obs["sasp_score"] = LN[sasp_g].mean(1)

    # ---------- サンプル単位 ----------
    grp = obs.groupby("sample", observed=True)
    df = pd.DataFrame({
        "group": grp[d["group"]].first(), "donor_id": grp["donor_id"].first(), "age": grp["age_num"].first(),
        "n_cells": grp.size(), "n_B": grp["isB"].sum(), "n_PL": grp["isPL"].sum(), "n_lymph": grp["isLy"].sum(),
    })
    df["frac_BPL_all"] = (df.n_B + df.n_PL) / df.n_cells
    df["frac_BPL_lymph"] = (df.n_B + df.n_PL) / df.n_lymph.replace(0, np.nan)
    tot_s = pd.Series(total).groupby(obs["sample"].values).sum()
    sums = G.groupby(obs["sample"].values).sum()
    cpm = sums.div(tot_s, axis=0) * 1e6
    for g in TLS8 + ["CDKN1A", "CDKN2A", "CD79A", "JCHAIN", "TNFSF8", "ITGAX", "TBX21"]:
        if g in cpm:
            df[f"CPM_{g}"] = cpm.loc[df.index, g]
    # z 化の基準: 腎組織サンプルのみ（A3 の血液・脂肪・副腎、A4 の PBMC・リンパ節は基準から外し、同じ変換を当てはめるだけ）
    ref = df.group.isin(d["ref_groups"]) if "ref_groups" in d else pd.Series(True, index=df.index)
    ref &= df.n_cells >= 300
    z = []
    for g in TLS8:
        if g in cpm:
            v = np.log1p(cpm.loc[df.index, g])
            z.append((v - v[ref].mean()) / (v[ref].std(ddof=0) or 1))
    df["TLS_score"] = pd.concat(z, axis=1).mean(1)
    df["TLS_genes_used"] = ",".join(g for g in TLS8 if g in cpm)
    df["in_score_reference"] = ref
    # TLS 群は基準サンプル数 ≥ 6 のときのみ（三分位）。細胞数 < 300 のサンプルと基準外サンプルは除外
    ok = ref
    df["TLS_group"] = "excluded(<300 cells or non-kidney)"
    if ok.sum() >= 6:
        q1, q2 = df.loc[ok, "TLS_score"].quantile([1 / 3, 2 / 3])
        df.loc[ok, "TLS_group"] = np.where(df.loc[ok, "TLS_score"] >= q2, "high", np.where(df.loc[ok, "TLS_score"] <= q1, "low", "mid"))
    df.sort_values("TLS_score", ascending=False).to_csv(OUT / f"donor_{key}.csv")
    obs["TLS_group"] = obs["sample"].map(df["TLS_group"])

    # ---------- 関連 ----------
    rows = []
    dd = df[ok]
    rows.append(dict(dataset=key, test="spearman TLS_score vs frac_BPL_all", n=len(dd),
                     stat=stats.spearmanr(dd.TLS_score, dd.frac_BPL_all).statistic, p=stats.spearmanr(dd.TLS_score, dd.frac_BPL_all).pvalue))
    if dd.age.notna().sum() >= 6 and dd.age.nunique() > 2:
        r = stats.spearmanr(dd.age, dd.TLS_score, nan_policy="omit")
        rows.append(dict(dataset=key, test="spearman TLS_score vs age", n=int(dd.age.notna().sum()), stat=r.statistic, p=r.pvalue))
        r = stats.spearmanr(dd.age, dd.frac_BPL_all, nan_policy="omit")
        rows.append(dict(dataset=key, test="spearman frac_BPL_all vs age", n=int(dd.age.notna().sum()), stat=r.statistic, p=r.pvalue))
    if key in ("A1", "A2", "A6") and dd.age.notna().sum() >= 10:
        import statsmodels.formula.api as smf
        m = smf.ols("TLS_score ~ age + C(group)", data=dd.assign(group=dd.group.astype(str))).fit()
        rows.append(dict(dataset=key, test="OLS TLS_score ~ age + disease: age coef (per year)", n=int(m.nobs), stat=m.params["age"], p=m.pvalues["age"]))
        nn = dd[dd.group.astype(str) == "normal"]
        if nn.age.nunique() > 2 and len(nn) >= 6:
            r = stats.spearmanr(nn.age, nn.TLS_score)
            rows.append(dict(dataset=key, test="spearman TLS_score vs age (normal only)", n=len(nn), stat=r.statistic, p=r.pvalue))
    levels = [g for g in dd.group.unique() if (dd.group == g).sum() >= 2]
    if len(levels) >= 2:
        arrs = [dd.loc[dd.group == g, "TLS_score"] for g in levels]
        r = stats.kruskal(*arrs)
        desc = "; ".join(f"{g}: n={len(a)}, median={a.median():.2f}" for g, a in zip(levels, arrs))
        rows.append(dict(dataset=key, test=f"kruskal TLS_score by {d['group']} [{desc}]", n=len(dd), stat=r.statistic, p=r.pvalue))
        ref = [g for g in levels if g.lower() in ("normal", "normal kidney")]
        for g in levels:
            if ref and g != ref[0]:
                r = stats.mannwhitneyu(dd.loc[dd.group == g, "TLS_score"], dd.loc[dd.group == ref[0], "TLS_score"])
                rows.append(dict(dataset=key, test=f"mannwhitney TLS_score {g} vs {ref[0]}", n=int((dd.group.isin([g, ref[0]])).sum()), stat=r.statistic, p=r.pvalue))
    # donor 内ペア（A3/A4/A5: 腫瘍 vs 正常の同一 donor 比較）
    if key in ("A3", "A5"):
        tum = "Tumour" if key == "A3" else "clear cell renal carcinoma"
        nor = "Normal kidney" if key == "A3" else "normal"
        pv = df.pivot_table(index="donor_id", columns="group", values="TLS_score")
        if tum in pv and nor in pv:
            pr = pv[[tum, nor]].dropna()
            if len(pr) >= 3:
                r = stats.wilcoxon(pr[tum], pr[nor])
                rows.append(dict(dataset=key, test=f"paired wilcoxon TLS_score {tum} vs {nor} (tumour>normal in {(pr[tum] > pr[nor]).sum()}/{len(pr)})", n=len(pr), stat=r.statistic, p=r.pvalue))

    # ---------- 細胞状態 ----------
    st = obs.groupby(d["state"], observed=True)
    sdf = pd.DataFrame({
        "n_cells": st.size(), "n_samples": st["sample"].nunique(),
        "n_TLS_high": st["TLS_group"].apply(lambda s: (s == "high").sum()),
        "n_TLS_low": st["TLS_group"].apply(lambda s: (s == "low").sum()),
        "n_samples_TLS_high_ge20cells": st.apply(lambda x: (x[x.TLS_group == "high"].groupby("sample").size() >= 20).sum()),
        "n_samples_TLS_low_ge20cells": st.apply(lambda x: (x[x.TLS_group == "low"].groupby("sample").size() >= 20).sum()),
        "sen_mean": st["sen_score"].mean(), "sen_sd": st["sen_score"].std(),
        "sasp_mean": st["sasp_score"].mean(), "sasp_sd": st["sasp_score"].std(),
        "median_genes": pd.Series(ngenes).groupby(obs[d["state"]].values).median(),
    })
    grp_levels = obs[d["group"]].astype(str).unique()
    for gl in grp_levels:
        sdf[f"n_{gl}"] = st[d["group"]].apply(lambda s, gl=gl: (s.astype(str) == gl).sum())
    for g in POS_FRAC_GENES:
        if g in G:
            sdf[f"pos_{g}"] = (G[g] > 0).groupby(obs[d["state"]].values).mean()
    # MPS 概算: 1 forward あたりの秒数（入力長 4096 で打ち切り）
    sdf["sec_per_forward_mps"] = np.minimum(sdf["median_genes"], 4096) / TOK_PER_SEC
    sdf.sort_values("n_cells", ascending=False).to_csv(OUT / f"states_{key}.csv")
    return pd.DataFrame(rows), ens_map


def token_table(ens_maps: dict):
    from geneformer import TOKEN_DICTIONARY_FILE
    with open(TOKEN_DICTIONARY_FILE, "rb") as f:
        tok = pickle.load(f)
    ens_all = {}
    for m in ens_maps.values():
        ens_all.update(m)
    rows = []
    for gset, genes in GENESETS.items():
        for g in genes:
            e = ens_all.get(g)
            rows.append(dict(gene=g, set=gset, ensembl=e, in_token_dict=bool(e and e in tok),
                             **{f"in_{k}": g in m for k, m in ens_maps.items()}))
    pd.DataFrame(rows).to_csv(OUT / "token_dict.csv", index=False)


if __name__ == "__main__":
    keys = sys.argv[1:] or list(DS)
    allrows, maps = [], {}
    for k in keys:
        print(f"[{k}] ...", flush=True)
        r, m = analyse(k)
        allrows.append(r)
        maps[k] = m
    summ = pd.concat(allrows, ignore_index=True)
    p = OUT / "assoc_summary.csv"
    if p.exists() and set(keys) != set(DS):
        old = pd.read_csv(p)
        summ = pd.concat([old[~old.dataset.isin(keys)], summ], ignore_index=True)
    summ.to_csv(p, index=False)
    if set(keys) == set(DS):
        token_table(maps)
    print(summ.to_string())

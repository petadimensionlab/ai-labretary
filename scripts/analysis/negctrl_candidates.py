#!/usr/bin/env python3
"""陰性対照（ネガティブコントロール）候補遺伝子の発現データによる評価（読み取り専用, GPU 不使用）。

PI 指示（2026-09-26）: HK（ACTB/GAPDH/PPIA）は合否判定に残しつつ、他の陰性対照候補（リボソームタンパク質
遺伝子など）を列挙し妥当性を評価する。**ISP / perturbation の結果は一切読まない**（指標検証の基準 C3 が
循環しないように、候補の選定は文献と発現データだけで行う）。

入力: experiments/EXP-20260926-001/work/input.h5ad（生カウント, isp_state ラベル付き。読み取りのみ）
      Geneformer V2 の token / gene median 辞書
出力: outputs/literature/2026-09-26-negative-controls/
  genome_stats.csv     全遺伝子（Geneformer 辞書収載）の状態別統計
  candidates.csv       事前固定の候補パネル（文献カテゴリ）＋データ駆動候補の統計と選定規則の判定
  panel.csv            事前登録用の選定規則を適用した推奨パネル（判定用 + 探索枠 RP 2 遺伝子）
  rule.json            規則の閾値と適用結果の要約
usage: .venv/bin/python scripts/analysis/negctrl_candidates.py
"""
import json, pickle, sys, warnings
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

warnings.filterwarnings("ignore")
sys.path = [p for p in sys.path if not p.endswith("/scripts") and not p.endswith("/scripts/analysis")]
from geneformer import TOKEN_DICTIONARY_FILE, GENE_MEDIAN_FILE  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
H5AD = ROOT / "experiments/EXP-20260926-001/work/input.h5ad"
OUT = ROOT / "outputs/literature/2026-09-26-negative-controls"
OUT.mkdir(parents=True, exist_ok=True)

START, GOAL, ALT = "aPT_CKD_TLShigh", "PT_healthy", "aPT_CKD_TLSlow"
STATES = {"start": START, "goal": GOAL, "alt": ALT}
MAX_GENES = 4094            # V2: model_input_size 4096 − <cls> − <eos>
MIN_CELLS_DONOR = 20        # pseudo-bulk に入れる donor の最小細胞数
EQ_MARGIN = 0.25            # TOST の同等性マージン（log2 CPM）

# ---- 事前固定の候補（文献のみで決定。カテゴリ → 遺伝子記号） ----------------------------
CANDIDATES = {
    "HK_current": ["ACTB", "GAPDH", "PPIA"],
    "HK_classic_qPCR": ["B2M", "HPRT1", "TBP", "YWHAZ", "SDHA", "HMBS", "GUSB", "PGK1", "TFRC",
                        "UBC", "PPIB", "POLR2A", "IPO8"],
    "HK_uniform_EisenbergLevanon2013": ["C1orf43", "CHMP2A", "EMC7", "GPI", "PSMB2", "PSMB4", "RAB7A",
                                        "REEP5", "SNRPD3", "VCP", "VPS29"],
    # Strauss 2021 PLoS One（ヒト腎生検 18 候補。YWHAZ のみ全データで安定, SLC4A1AP が次点）
    "HK_kidney_Strauss2021": ["SLC4A1AP", "ACTG1", "REL", "NOP10", "PGAM1", "RPL13", "HNRNPL",
                              "PCBP1", "RER1", "PLA2G4A"],
    # Herath 2020（ラット虚血・腎毒性: PABPN1×ACTB, HMBS）/ HRT Atlas の PBMC 例（EIF2S3, SF3B1）
    "HK_injury_or_HRT": ["PABPN1", "SF3B1", "EIF2S3"],
    # Lin 2019 scSEG 本文で例示された SEG（上位リストは未確認）と bulk HK だが SEG でない例
    "SEG_scSEG_Lin2019_examples": ["RPL26", "RPL36", "HINT1"],
    "RP_ribosomal": ["RPLP0", "RPL13A", "RPS18", "RPL19", "RPS27A", "RPL32", "RPS13", "RPL27",
                     "RPL30", "RPS3", "RPL10", "RPS6", "RPL41", "RPS29"],
}
# 既知の傷害・老化・TLS・IFN 応答・解離ストレス関連（陰性対照から除外。事前固定）
EXCLUDE_SYMBOLS = {
    # 本系の対照・対象・参照
    "VCAM1", "HAVCR1", "LCN2", "SPP1", "CDKN1A", "TP53", "LTBR",
    # 傷害 / adaptive PT マーカー（Lake 2023 ほか）
    "CDKN2A", "PROM1", "CD24", "VIM", "SOX9", "ITGB8", "DCDC2", "MMP7", "KRT7", "KRT8", "KRT18",
    "KRT19", "KRT23", "CCL2", "SERPINE1", "TGFB1", "IL6", "CXCL8", "CXCL1", "CXCL2", "CST3", "TPM1",
    "CP", "C3", "PDGFB", "FGB", "FGA", "FGG", "S100A6", "ANXA3", "SPP1", "TNC", "CLU", "LGALS3",
    # TLS / 免疫
    "CXCL13", "CCL19", "CCL21", "LTB", "LTA", "TNFSF8", "TNFRSF8", "RELB", "NFKB2", "NFKBIA",
    # 解離ストレス（immediate early）
    "FOS", "FOSB", "JUN", "JUNB", "EGR1", "ATF3", "IER2", "IER3", "HSPA1A", "HSPA1B", "DUSP1",
    # 性染色体特異（donor の性比の偏りを拾う）
    "XIST", "TSIX", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "ZFY",
}
EXCLUDE_SYMBOLS |= {
    # 細胞周期・増殖（傷害後の再生で変わる）
    "PCNA", "MKI67", "TOP2A", "CCNB1", "CCNB2", "CDK1", "CCNA2", "BIRC5", "TYMS", "MCM2", "MCM7",
    "STMN1", "HMGB2", "TUBA1B", "UBE2C", "CENPF",
    # ストレス応答（metallothionein・熱ショック・酸化ストレス）
    "MT2A", "HSP90AA1", "HSP90AB1", "HSPB1", "DNAJB1", "HMOX1", "TXNIP", "SOD2", "GADD45A", "GADD45B",
}
EXCLUDE_PREFIX = ("MT-", "MT1", "HSPA", "HSPH", "HLA-", "IGH", "IGK", "IGL", "H1-", "H2A", "H2B", "H3-",
                  "H4C", "HIST")
ROLL_WIN = 201              # 発現頻度（start の tokenized 検出細胞数）が近い 201 遺伝子の移動中央値で中心化
POS = ["VCAM1", "HAVCR1"]
REF = ["TP53", "CDKN1A", "LTBR", "SPP1", "LCN2"]

# ---- 事前登録する選定規則（本スクリプトの値がそのまま登録値） -----------------------------
RULE = {
    "R1_min_start_cells_tokenized": 1000,   # 200 細胞/遺伝子（donor ≤ 25）を確実に抽出できる
    "R1_min_start_donors_ge25": 8,          # 25 細胞以上発現している start donor 数
    "R2_abs_log2fc_goal_max": 0.25,         # pseudo-bulk（donor 平均 log2CPM）start vs goal
    "R2_abs_log2fc_alt_max": 0.50,          # start vs alt（alt は 8 donor のみ → 緩め）
    "R2_tost_p_max": 0.05,                  # start vs goal の同等性（±0.25 log2）を donor 水準で示す
    "R3_abs_det_diff_max": 0.05,            # 検出率の差 start − goal
    "R4_abs_relrank_diff_max": 0.03,        # Geneformer 入力内の相対順位（0 = 先頭）の差 start − goal
    "R5_donor_cv_max_quantile": 0.75,       # donor 間 CV（start, goal）が適格遺伝子の 75% 分位以下
    "R6_exclude_known": True,               # EXCLUDE_SYMBOLS / 接頭辞 / 陽性・対象・参照
}


def load():
    a = ad.read_h5ad(H5AD)        # 213 MB, メモリに載せる（ファイルは変更しない）
    X = a.X.tocsr().astype(np.float32)
    obs = a.obs[["isp_state", "donor_id", "n_counts"]].copy()
    var = a.var[["ensembl_id", "gene_name"]].copy()
    return X, obs, var


def main():
    X, obs, var = load()
    tok = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
    med = pickle.load(open(GENE_MEDIAN_FILE, "rb"))
    ens = var["ensembl_id"].astype(str).values
    sym = var["gene_name"].astype(str).values
    in_vocab = np.array([(e in tok) and (e in med) for e in ens])
    medv = np.array([med.get(e, np.nan) for e in ens], dtype=np.float64)
    nc = obs["n_counts"].values.astype(np.float64)
    rowsum = np.asarray(X.sum(1)).ravel()
    print(f"cells {X.shape[0]}, genes {X.shape[1]}, in vocab {in_vocab.sum()}, "
          f"n_counts==rowsum: {np.allclose(nc, rowsum)}")

    G = X.shape[1]
    res = pd.DataFrame({"ensembl_id": ens, "symbol": sym, "in_gf_vocab": in_vocab})
    state = obs["isp_state"].astype(str).values
    donor = obs["donor_id"].astype(str).values

    # --- Geneformer 入力の順位（tokenizer と同じ: count / n_counts * 1e4 / median, 降順, 4094 で切る）
    Xv = X[:, in_vocab].tocsr()
    medsub = medv[in_vocab]
    vidx = np.where(in_vocab)[0]
    Gv = len(vidx)
    relrank_sum = {k: np.zeros(Gv) for k in STATES}
    absrank_sum = {k: np.zeros(Gv) for k in STATES}
    tok_det = {k: np.zeros(Gv) for k in STATES}
    top10_cnt = {k: np.zeros(Gv) for k in STATES}
    ntok = np.zeros(X.shape[0])
    for i in range(X.shape[0]):
        s, e = Xv.indptr[i], Xv.indptr[i + 1]
        cols = Xv.indices[s:e]
        vals = Xv.data[s:e] / nc[i] * 1e4 / medsub[cols]
        order = cols[np.argsort(-vals, kind="stable")][:MAX_GENES]
        n = len(order)
        ntok[i] = n
        key = [k for k, v in STATES.items() if v == state[i]]
        if not key:
            continue
        k = key[0]
        r = np.arange(n)
        relrank_sum[k][order] += r / max(n - 1, 1)
        absrank_sum[k][order] += r
        tok_det[k][order] += 1
        top10_cnt[k][order[:10]] += 1

    # --- 状態ごとの検出率・CP10k・donor pseudo-bulk
    cp10k = sp.diags(1e4 / nc) @ X
    donor_tab = {}
    for k, st in STATES.items():
        m = state == st
        ncell = m.sum()
        res[f"det_{k}"] = np.asarray((X[m] > 0).sum(0)).ravel() / ncell
        res[f"cp10k_{k}"] = np.asarray(cp10k[m].mean(0)).ravel()
        tdet = np.full(G, np.nan); tdet[vidx] = tok_det[k]
        res[f"tok_cells_{k}"] = tdet
        rr = np.full(G, np.nan); rr[vidx] = np.where(tok_det[k] > 0, relrank_sum[k] / np.maximum(tok_det[k], 1), np.nan)
        res[f"relrank_{k}"] = rr
        ar = np.full(G, np.nan); ar[vidx] = np.where(tok_det[k] > 0, absrank_sum[k] / np.maximum(tok_det[k], 1), np.nan)
        res[f"absrank_{k}"] = ar
        t10 = np.full(G, np.nan); t10[vidx] = top10_cnt[k] / ncell
        res[f"frac_top10_{k}"] = t10
        # pseudo-bulk: donor ごとの合計カウント → CPM
        ds = [d for d in np.unique(donor[m]) if (m & (donor == d)).sum() >= MIN_CELLS_DONOR]
        pb = np.vstack([np.asarray(X[m & (donor == d)].sum(0)).ravel() for d in ds])
        cpm = pb / pb.sum(1, keepdims=True) * 1e6
        donor_tab[k] = np.log2(cpm + 1)
        res[f"pb_log2cpm_{k}"] = donor_tab[k].mean(0)
        mu = cpm.mean(0)
        res[f"donor_cv_{k}"] = np.where(mu > 0, cpm.std(0, ddof=1) / np.where(mu > 0, mu, 1), np.nan)
        # 25 細胞以上発現している donor 数
        res[f"donors_ge25_{k}"] = np.vstack([np.asarray((X[m & (donor == d)] > 0).sum(0)).ravel() >= 25
                                             for d in np.unique(donor[m])]).sum(0)
        res.attrs[f"n_donors_{k}"] = len(ds)
        print(f"{k}: {ncell} cells, {len(ds)} donors (≥{MIN_CELLS_DONOR} cells), "
              f"median tokens/cell {np.median(ntok[m]):.0f}")

    for k in ("goal", "alt"):
        a_, b_ = donor_tab["start"], donor_tab[k]
        res[f"log2fc_start_vs_{k}"] = a_.mean(0) - b_.mean(0)
        t, p = stats.ttest_ind(a_, b_, equal_var=False)
        res[f"welch_p_{k}"] = p
        res[f"det_diff_start_minus_{k}"] = res["det_start"] - res[f"det_{k}"]
        res[f"relrank_diff_start_minus_{k}"] = res["relrank_start"] - res[f"relrank_{k}"]
    # TOST（Welch, ±EQ_MARGIN）start vs goal
    a_, b_ = donor_tab["start"], donor_tab["goal"]
    d = a_.mean(0) - b_.mean(0)
    va, vb = a_.var(0, ddof=1) / len(a_), b_.var(0, ddof=1) / len(b_)
    se = np.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / (len(a_) - 1) + vb ** 2 / (len(b_) - 1) + 1e-300)
    with np.errstate(invalid="ignore", divide="ignore"):
        p_lo = stats.t.sf((d + EQ_MARGIN) / se, df)
        p_hi = stats.t.cdf((d - EQ_MARGIN) / se, df)
    res["tost_p_goal"] = np.maximum(p_lo, p_hi)
    # BH（検出率 ≥ 1% の遺伝子）
    ok = (res[["det_start", "det_goal"]].max(1) >= 0.01) & res["welch_p_goal"].notna()
    q = np.full(G, np.nan)
    q[ok.values] = stats.false_discovery_control(res.loc[ok, "welch_p_goal"].values)
    res["welch_q_goal"] = q

    # --- 除外フラグ
    excl = np.array([(s in EXCLUDE_SYMBOLS) or s.startswith(EXCLUDE_PREFIX) for s in sym])
    res["known_excluded"] = excl

    # --- 発現頻度マッチの期待値で中心化（null は発現頻度マッチで抽出されるため、陰性対照は
    #     「同じ頻度帯の遺伝子と比べて」状態中立であればよい。start は深く（UMI 中央値 3,405 vs 2,819）、
    #     発現遺伝子が多いので、生の log2FC・検出率差は全体に正へずれる）
    pool = res["in_gf_vocab"] & (res["tok_cells_start"] >= 300) & ~res["known_excluded"]
    # 期待値 = null 抽出プール（start で ≥ 300 細胞・除外なし）を検出細胞数で並べ、近傍 ROLL_WIN 遺伝子の
    # 移動中央値。プール外の遺伝子は検出細胞数で線形補間（端は端の値）
    pp = res.loc[pool].sort_values("tok_cells_start")
    xs = pp["tok_cells_start"].values.astype(float)
    exp_fn = {}
    for c in ("log2fc_start_vs_goal", "log2fc_start_vs_alt", "det_diff_start_minus_goal",
              "relrank_diff_start_minus_goal"):
        rm = pp[c].rolling(ROLL_WIN, center=True, min_periods=ROLL_WIN // 2).median().values
        # 同じ x の重複は平均して単調化
        ux, inv = np.unique(xs, return_inverse=True)
        uy = np.bincount(inv, weights=np.nan_to_num(rm)) / np.bincount(inv)
        expv = np.interp(res["tok_cells_start"].fillna(0).values, ux, uy)
        res[c + "_exp"] = expv
        res[c + "_c"] = res[c] - expv
        exp_fn[c] = (ux, uy)
    shift_goal = res["log2fc_start_vs_goal_exp"]
    with np.errstate(invalid="ignore", divide="ignore"):
        dc = d - shift_goal.values
        res["tost_p_goal_c"] = np.maximum(stats.t.sf((dc + EQ_MARGIN) / se, df), stats.t.cdf((dc - EQ_MARGIN) / se, df))

    # --- 規則の判定
    R = RULE
    res["R1"] = (res["tok_cells_start"] >= R["R1_min_start_cells_tokenized"]) & \
                (res["donors_ge25_start"] >= R["R1_min_start_donors_ge25"])
    res["R2"] = (res["log2fc_start_vs_goal_c"].abs() < R["R2_abs_log2fc_goal_max"]) & \
                (res["log2fc_start_vs_alt_c"].abs() < R["R2_abs_log2fc_alt_max"]) & \
                (res["tost_p_goal_c"] < R["R2_tost_p_max"])
    res["R3"] = res["det_diff_start_minus_goal_c"].abs() < R["R3_abs_det_diff_max"]
    res["R4"] = res["relrank_diff_start_minus_goal_c"].abs() < R["R4_abs_relrank_diff_max"]
    # 参考: 中心化しない（生の）判定
    res["R2_raw"] = (res["log2fc_start_vs_goal"].abs() < R["R2_abs_log2fc_goal_max"]) & \
                    (res["log2fc_start_vs_alt"].abs() < R["R2_abs_log2fc_alt_max"]) & (res["tost_p_goal"] < R["R2_tost_p_max"])
    res["R3_raw"] = res["det_diff_start_minus_goal"].abs() < R["R3_abs_det_diff_max"]
    res["R4_raw"] = res["relrank_diff_start_minus_goal"].abs() < R["R4_abs_relrank_diff_max"]
    elig = res["in_gf_vocab"] & res["R1"]
    cv_thr = {k: res.loc[elig, f"donor_cv_{k}"].quantile(R["R5_donor_cv_max_quantile"]) for k in ("start", "goal")}
    res["R5"] = (res["donor_cv_start"] <= cv_thr["start"]) & (res["donor_cv_goal"] <= cv_thr["goal"])
    res["R6"] = ~res["known_excluded"]
    rcols = ["R1", "R2", "R3", "R4", "R5", "R6"]
    res["n_rules_pass"] = res[rcols].sum(1)
    res["pass_all"] = res["in_gf_vocab"] & res[rcols].all(1)
    # 同等性の総合スコア（小さいほど状態中立。正規化した偏差の最大値）
    res["neutrality_score"] = np.nanmax(np.vstack([
        res["log2fc_start_vs_goal_c"].abs() / R["R2_abs_log2fc_goal_max"],
        res["log2fc_start_vs_alt_c"].abs() / R["R2_abs_log2fc_alt_max"],
        res["det_diff_start_minus_goal_c"].abs() / R["R3_abs_det_diff_max"],
        res["relrank_diff_start_minus_goal_c"].abs() / R["R4_abs_relrank_diff_max"]]), axis=0)
    res["pass_all_raw"] = res["in_gf_vocab"] & res[["R1", "R2_raw", "R3_raw", "R4_raw", "R5", "R6"]].all(1)

    keep = res["in_gf_vocab"] & (res[["det_start", "det_goal"]].max(1) >= 0.01)
    res.loc[keep].to_csv(OUT / "genome_stats.csv", index=False, float_format="%.5g")

    # --- 候補表（文献カテゴリ + 陽性・参照 + データ駆動）
    s2i = {}
    for i, s in enumerate(sym):
        s2i.setdefault(s, i)
    rows = []
    for cat, genes in CANDIDATES.items():
        for g in genes:
            if g not in s2i:
                rows.append({"category": cat, "symbol": g, "note": "not in var"}); continue
            rows.append({"category": cat, **res.iloc[s2i[g]].to_dict()})
    for cat, genes in (("positive_control", POS), ("ref_target", REF)):
        for g in genes:
            rows.append({"category": cat, **res.iloc[s2i[g]].to_dict()})
    # データ駆動: 規則をすべて満たす遺伝子を検出率の 3 帯（start の tokenized 検出率）で上位 5 ずつ
    pa = res[res["pass_all"]].copy()
    n_start = int((state == START).sum())
    pa["det_tier"] = pd.cut(pa["tok_cells_start"] / n_start, [0, 0.4, 0.7, 1.01], labels=["low", "mid", "high"])
    lit = {g for v in CANDIDATES.values() for g in v}
    dd = (pa[~pa["symbol"].isin(lit)].sort_values("neutrality_score")
          .groupby("det_tier", observed=True).head(5))
    for _, r in dd.iterrows():
        rows.append({"category": f"data_driven_{r['det_tier']}", **r.drop("det_tier").to_dict()})
    cand = pd.DataFrame(rows)
    cand.to_csv(OUT / "candidates.csv", index=False, float_format="%.5g")

    # --- 推奨パネル（事前登録の決定規則）
    # 1) 規則をすべて満たす文献候補（HK_classic / uniform / SEG / RP）を neutrality_score 順
    # 2) 検出率の 3 帯（low < 0.4 ≤ mid < 0.7 ≤ high）それぞれから、文献候補を優先して最大 3 遺伝子、
    #    帯に文献候補が足りなければデータ駆動で補う → 最大 9 遺伝子
    # 3) RP 遺伝子は各帯で 1 つまで（同じ機能群に偏らない）
    panel = []
    lit_pass = cand[cand["category"].isin([c for c in CANDIDATES if c != "HK_current"]) &
                    cand.get("pass_all", False).fillna(False).astype(bool)].copy()
    lit_pass = lit_pass.drop_duplicates("symbol")
    lit_pass["det_tier"] = pd.cut(lit_pass["tok_cells_start"] / n_start, [0, 0.4, 0.7, 1.01],
                                  labels=["low", "mid", "high"])
    dd2 = pa[~pa["symbol"].isin(lit)].sort_values("neutrality_score")
    for tier in ("low", "mid", "high"):
        chosen = []
        for _, r in lit_pass[lit_pass["det_tier"] == tier].sort_values("neutrality_score").iterrows():
            if r["category"] == "RP_ribosomal" and any(p["category"] == "RP_ribosomal" for p in chosen):
                continue
            chosen.append({**r.to_dict(), "source": "literature"})
            if len(chosen) == 3:
                break
        for _, r in dd2[dd2["det_tier"] == tier].iterrows():
            if len(chosen) >= 3:
                break
            chosen.append({"category": "data_driven", **r.to_dict(), "source": "data"})
        panel += chosen
    # 探索枠（合否には使わない）: RP 候補のうち TOST 以外の条件（R1, |log2FC_c|, R3, R4, R5, R6）を満たし
    # neutrality_score が小さい順に 2 遺伝子
    rp = cand[(cand["category"] == "RP_ribosomal") & cand["R1"].fillna(False).astype(bool)].drop_duplicates("symbol")
    rp_ok = rp[(rp["log2fc_start_vs_goal_c"].abs() < RULE["R2_abs_log2fc_goal_max"]) &
               (rp["log2fc_start_vs_alt_c"].abs() < RULE["R2_abs_log2fc_alt_max"]) &
               rp[["R3", "R4", "R5", "R6"]].astype(bool).all(1)].sort_values("neutrality_score").head(2)
    for _, r in rp_ok.iterrows():
        panel.append({**r.to_dict(), "source": "exploratory_RP"})
    pan = pd.DataFrame(panel)
    pan.to_csv(OUT / "panel.csv", index=False, float_format="%.5g")

    summ = {"rule": RULE, "eq_margin_log2": EQ_MARGIN, "cv_thresholds": cv_thr,
            "n_start_cells": n_start, "n_donors": {k: res.attrs[f"n_donors_{k}"] for k in STATES},
            "n_genes_vocab_det1pct": int(keep.sum()), "n_eligible_R1": int(elig.sum()),
            "n_pass_all": int(res["pass_all"].sum()),
            "n_pass_each": {c: int((elig & res[c]).sum()) for c in rcols},
            "panel": pan[["symbol", "ensembl_id", "category", "source"]].to_dict("records")}
    (OUT / "rule.json").write_text(json.dumps(summ, indent=1, ensure_ascii=False, default=float) + "\n")
    print(json.dumps(summ, indent=1, ensure_ascii=False, default=float))
    grid = [300, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000]
    summ["centering_expectation_at_tok_cells_start"] = {"grid": grid, **{c: np.interp(grid, *exp_fn[c]).round(4).tolist()
        for c in exp_fn}}
    summ["n_pass_all_raw"] = int(res["pass_all_raw"].sum())
    (OUT / "rule.json").write_text(json.dumps(summ, indent=1, ensure_ascii=False, default=float) + "\n")
    cols = ["category", "symbol", "tok_cells_start", "det_start", "det_goal", "log2fc_start_vs_goal",
            "log2fc_start_vs_goal_c", "welch_q_goal", "tost_p_goal_c", "log2fc_start_vs_alt_c", "det_diff_start_minus_goal_c",
            "relrank_start", "relrank_goal", "relrank_diff_start_minus_goal", "relrank_diff_start_minus_goal_c",
            "donor_cv_start", "n_rules_pass", "pass_all"]
    with pd.option_context("display.width", 250, "display.max_rows", 200, "display.max_columns", 30):
        print(cand[cols].round(3).to_string(index=False))


if __name__ == "__main__":
    main()

"""EXP-20260927-001 評価（事前登録。plan.md §6。ISP 実行前に sha256 を固定する）。

主検定: 移動距離をそろえた並べ替え検定（pos_matched_permutation）
  対象遺伝子 g・seed s の細胞 i（値 y_i = shift_goal, 位置 p_i, 長さ L_i）について、
  null 遺伝子 h（g 自身は除く）の細胞 j との距離 d(i,j) = max(|log2(1+p_i) − log2(1+p_j)|, |log2 L_i − log2 L_j|)。
  g の細胞の 80% 以上に d ≤ 0.5 の j がある h を「使える null 遺伝子」とする（20 未満なら NaN）。
  細胞ごとに偏り補正した値 y_i − f(p_i, L_i)（f は対象遺伝子を除く同じ seed の null 細胞で当てはめた
  log2(1+pos)・log2 len の 2 次式）の遺伝子平均 D_raw を、使える null 遺伝子の同じ値と比べる（遺伝子ラベルの並べ替え）。
  q = 使える null 遺伝子のうち D_raw が対象より小さい割合、補正後の差 D = D_raw(g) − 使える null 遺伝子の平均。
副検定: 区間内パーセンタイル検定（within_bin_percentile）。判定には使わない。
提供者単位の副指標・固定遺伝子集合の半分割再現性・EXP-005 と同じ a2（参考）も出す。判定は主検定だけで行う。

実行: .venv/bin/python experiments/EXP-20260927-001/analysis/evaluate_matched.py --config experiments/EXP-20260927-001/config.yaml
      [--cells <cells.csv>]（既定は <results_dir>/cells.csv。列 pos があれば tokenized を読まない = 合成データの検査用）
出力: <results_dir>/eval/{criteria.json, gene_seed.csv, donor_level.csv, split_half.json}
"""
from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
BLIND = {"ENSG00000141510", "ENSG00000124762", "ENSG00000111321", "ENSG00000118785", "ENSG00000148346"}
# TP53, CDKN1A, LTBR, SPP1, LCN2（盲検。読み込み直後に除外して assert）


# ---------------------------------------------------------------- 入力
def load_cells(cfg: dict, path: Path) -> pd.DataFrame:
    c = pd.read_csv(path)
    c = c[~c.gene.isin(BLIND)].reset_index(drop=True)
    assert not set(c.gene) & BLIND, "盲検の遺伝子が結果に含まれています"
    c["seed"] = c.run.astype(str).str.extract(r"run-(\d+)", expand=False).astype(int)
    if "pos" not in c.columns:
        from datasets import load_from_disk
        from geneformer import TOKEN_DICTIONARY_FILE
        tokd = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
        ds = load_from_disk(str(ROOT / cfg["data"]["tokenized_dataset"]))
        uids = sorted(c.cell_uid.unique())
        ids = dict(zip(uids, ds.select(uids)["input_ids"]))
        c["pos"] = [ids[u].index(tokd[g]) for u, g in zip(c.cell_uid, c.gene)]
    c["lp"] = np.log2(1 + c.pos.astype(float))
    c["ll"] = np.log2(c.length.astype(float))
    c["relpos"] = c.pos / c.length
    c["y"] = c[cfg["analysis"]["base_metric"].replace("m_raw", "shift_goal")]
    return c


def roles(cfg: dict) -> dict[str, list[str]]:
    ctl = cfg["controls"]
    r = {"positive": list(ctl["positive"]), "housekeeping": list(ctl["negative_fixed"]),
         "negative": list(ctl["negative_sampler"]["genes"])}
    for grp in ctl.get("extra_fixed", []) or []:
        r[grp["role"]] = list(grp["genes"])
    return r


def comparator_map(cfg: dict) -> dict[str, list[str]]:
    """陽性対照 → 比較遺伝子 5 つ（config のコメントではなく work/comparators.csv の ref 列で対応づける）。"""
    t = pd.read_csv(ROOT / "experiments/EXP-20260927-001/work/comparators.csv")
    name2id = dict(zip(pd.read_csv(ROOT / "outputs/literature/2026-09-26-negative-controls/genome_stats.csv").symbol,
                       pd.read_csv(ROOT / "outputs/literature/2026-09-26-negative-controls/genome_stats.csv").ensembl_id))
    out: dict[str, list[str]] = {}
    for ref, s in t.groupby("ref", sort=False):
        out[name2id[ref]] = s.ensembl_id.tolist()
    assert set(out) == set(cfg["controls"]["positive"]) and all(len(v) == 5 for v in out.values())
    return out


# ---------------------------------------------------------------- 主検定
def dist(a_lp, a_ll, b_lp, b_ll):
    return np.maximum(np.abs(a_lp[:, None] - b_lp[None, :]), np.abs(a_ll[:, None] - b_ll[None, :]))


def design(lp, ll):
    """偏り補正の回帰の説明変数（2 次まで）。"""
    lp, ll = np.asarray(lp, float), np.asarray(ll, float)
    return np.column_stack([np.ones_like(lp), lp, ll, lp ** 2, ll ** 2, lp * ll])


def fit_bias(ref: pd.DataFrame):
    """null 細胞（対象遺伝子を除く）で y ~ 2 次式(log2(1+pos), log2 len) を最小二乗で当てはめ、予測関数を返す。
    許容幅内に残る位置・長さの差の分だけ、引いた null 細胞の値をずらす（偏り補正つきマッチング, plan.md §6.1）。"""
    beta, *_ = np.linalg.lstsq(design(ref.lp, ref.ll), ref.y.to_numpy(), rcond=None)
    return lambda lp, ll: design(lp, ll) @ beta


def perm_test(tg: pd.DataFrame, nulls: dict[str, pd.DataFrame], exclude: str | None, p: dict, rng=None) -> dict:
    """1 遺伝子×1 seed の移動距離をそろえた並べ替え検定（遺伝子ラベルの並べ替え, 偏り補正つき）。
    D_x = mean_i [y_i − f(lp_i, ll_i)]（f は対象遺伝子を除く同じ seed の null 細胞で当てはめた 2 次式）。
    対象 g の細胞の min_cover_frac 以上に d ≤ tol の細胞を持つ null 遺伝子 h（= 移動距離と長さの分布が重なるもの）だけを
    並べ替えの相手にし、帰無仮説（g はそれらと交換可能）のもとで q = #{h: D_h < D_g} / #{h}。"""
    tol, frac = p["tolerance_log2"], p["min_cover_frac"]
    fb = fit_bias(pd.concat([hd for h, hd in nulls.items() if h != exclude]))
    Dg = float((tg.y.to_numpy() - fb(tg.lp, tg.ll)).mean())
    Dh = []
    for h, hd in nulls.items():
        if h == exclude:
            continue
        D = dist(tg.lp.to_numpy(), tg.ll.to_numpy(), hd.lp.to_numpy(), hd.ll.to_numpy())
        if (D <= tol).any(axis=1).mean() < frac:
            continue
        Dh.append(float((hd.y.to_numpy() - fb(hd.lp, hd.ll)).mean()))
    Dh = np.array(Dh)
    base = {"T": float(tg.y.mean()), "D_raw": Dg, "n_elig": int(len(Dh))}
    if len(Dh) < p["min_eligible_null_genes"]:
        return {**base, "q": np.nan, "D": np.nan, "null_mean": np.nan, "null_sd": np.nan}
    return {**base, "q": float((Dh < Dg).mean()), "D": float(Dg - Dh.mean()),
            "null_mean": float(Dh.mean()), "null_sd": float(Dh.std(ddof=1))}


def cell_residual(c: pd.DataFrame, k: int) -> np.ndarray:
    """r_i = y_i − （同じ seed の null 細胞のうち d が最も近い k 個の偏り補正済みの y の平均）。null 遺伝子自身の細胞は除く。"""
    r = np.full(len(c), np.nan)
    for s, cs in c.groupby("seed"):
        nl = cs[cs.role == "negative"]
        for g, tg in cs.groupby("gene"):
            ref = nl[nl.gene != g]
            fb = fit_bias(ref)
            D = dist(tg.lp.to_numpy(), tg.ll.to_numpy(), ref.lp.to_numpy(), ref.ll.to_numpy())
            nn = np.argpartition(D, k, axis=1)[:, :k]
            adj = (ref.y.to_numpy() - fb(ref.lp, ref.ll))[nn].mean(axis=1) + fb(tg.lp, tg.ll)
            r[tg.index.to_numpy()] = tg.y.to_numpy() - adj
    return r


# ---------------------------------------------------------------- 副検定・参考
def within_bin_score(tg: pd.DataFrame, ref: pd.DataFrame, nl: int, npb: int) -> float:
    le = np.quantile(ref.ll, np.linspace(0, 1, nl + 1)[1:-1])
    rb, tb = np.digitize(ref.ll, le), np.digitize(tg.ll, le)
    u = np.full(len(tg), np.nan)
    for b in np.unique(tb):
        R, m = ref[rb == b], tb == b
        if len(R) == 0:
            continue
        pe = np.quantile(R.lp, np.linspace(0, 1, npb + 1)[1:-1])
        rpb, tpb = np.digitize(R.lp, pe), np.digitize(tg.lp[m], pe)
        yt = tg.y.to_numpy()[m]
        uu = np.empty(m.sum())
        for k in range(len(yt)):
            ry = R.y.to_numpy()[rpb == tpb[k]]
            uu[k] = (ry < yt[k]).mean() if len(ry) else np.nan
        u[np.flatnonzero(m)] = uu
    return float(np.nanmean(u) - 0.5)


def a2_like(cs: pd.DataFrame, nbins: int = 5) -> dict[str, float]:
    """EXP-005 と同じ考え方の a2（長さ 5 分位ビンの null 平均を引く。null 遺伝子は自分を除く）。"""
    out = {}
    q = np.linspace(0, 1, nbins + 1)[1:-1]
    for g, tg in cs.groupby("gene"):
        F = cs[(cs.role == "negative") & (cs.gene != g)]
        e = np.quantile(F.length, q)
        mu = F.groupby(np.digitize(F.length, e)).y.mean()
        out[g] = float((tg.y.to_numpy() - mu.reindex(np.digitize(tg.length, e)).to_numpy()).mean())
    return out


def null_pct(x: float, null: np.ndarray) -> float:
    null = null[~np.isnan(null)]
    return float((null < x).mean()) if len(null) and not np.isnan(x) else np.nan


# ---------------------------------------------------------------- 判定
def assess(g: pd.DataFrame, qcol: str, dcol: str, R: dict, comp: dict, cr: dict, ann: dict, seeds: list[int]) -> dict:
    Q = g.pivot(index="gene", columns="seed", values=qcol)
    Dm = g.pivot(index="gene", columns="seed", values=dcol)
    need = cr["min_seeds"]

    def n_hi(x):
        return int((Q.loc[x] >= cr["q_hi"]).sum()) if x in Q.index else 0

    def in_rng(x, lo, hi):
        return int(((Q.loc[x] >= lo) & (Q.loc[x] <= hi)).sum()) if x in Q.index else 0  # NaN は範囲外

    per = {}
    for x in R["positive"]:
        c1 = n_hi(x) >= need
        wins = int(sum((Dm.loc[x, s] > Dm.loc[h, s]) for s in seeds for h in R["housekeeping"]
                       if not (np.isnan(Dm.loc[x, s]) or np.isnan(Dm.loc[h, s]))))
        beat = {s: int(sum(Dm.loc[x, s] > Dm.loc[cg, s] for cg in comp[x]
                           if not (np.isnan(Dm.loc[x, s]) or np.isnan(Dm.loc[cg, s])))) for s in seeds}
        c9 = sum(v >= cr["C9_min_beaten"] for v in beat.values()) >= cr["C9_min_seeds"]
        per[x] = {"seeds_q_ge_hi": n_hi(x), "C1_gene": c1, "C2_wins_of_9": wins, "C2_gene": wins >= cr["C2_min_wins"],
                  "C9_beaten_by_seed": beat, "C9_gene": c9, "pass": bool(c1 and wins >= cr["C2_min_wins"] and c9),
                  "q_by_seed": {int(s): (None if np.isnan(Q.loc[x, s]) else round(float(Q.loc[x, s]), 3)) for s in seeds},
                  "ddet_pct": ann.get(x)}
    lo, hi = cr["C3_range"]
    c3 = {x: in_rng(x, lo, hi) for x in R["housekeeping"]}
    cp = {x: in_rng(x, lo, hi) for x in R["negctrl_panel"]}
    n_pass = sum(v["pass"] for v in per.values())
    n_c1 = sum(v["C1_gene"] for v in per.values())
    # 較正（C8a）: null 遺伝子自身の分位
    qn = g[g.role == "negative"][qcol].to_numpy()
    qn = qn[~np.isnan(qn)]
    c8a = {"n": int(len(qn)), "frac_hi": float((qn >= 0.90).mean()) if len(qn) else np.nan,
           "frac_lo": float((qn < 0.10).mean()) if len(qn) else np.nan}
    c8a["pass"] = bool(len(qn) >= 60 and c8a["frac_hi"] <= cr["C8a_tail_frac_max"] and c8a["frac_lo"] <= cr["C8a_tail_frac_max"])
    # 交絡の除去（C8b）: null 遺伝子の D と元の位置
    rho = {}
    for s in seeds:
        n = g[(g.role == "negative") & (g.seed == s)].dropna(subset=[dcol])
        rho[int(s)] = {"lp": float(stats.spearmanr(n.lp_med, n[dcol]).statistic),
                       "relpos": float(stats.spearmanr(n.relpos_med, n[dcol]).statistic), "n": int(len(n))}
    ok = sum(abs(v["lp"]) <= cr["C8b_abs_rho_max"] and abs(v["relpos"]) <= cr["C8b_abs_rho_max"] for v in rho.values())
    c8b = {"rho_by_seed": rho, "pass": ok >= cr["C8b_min_seeds"]}
    rho_all = {int(s): float(stats.spearmanr(g[g.seed == s].dropna(subset=[dcol]).lp_med,
                                             g[g.seed == s].dropna(subset=[dcol])[dcol]).statistic) for s in seeds}
    c7 = {x: int((Q.loc[x] <= cr["C7_reverse_max"]).sum()) if x in Q.index else 0 for x in R.get("reverse_direction", [])}
    res = {
        "per_positive": per, "n_positive_pass": n_pass, "n_positive_C1_only": n_c1,
        "C1_C2_C9": n_pass >= cr["C1_min_genes"], "C1_null_only": n_c1 >= cr["C1_min_genes"],
        "C3_seeds_in_range": c3, "C3": all(v >= need for v in c3.values()),
        "C3_panel_seeds_in_range": cp, "C3_panel_n_ok": sum(v >= need for v in cp.values()),
        "C8a": c8a, "C8b": c8b, "rho_all_genes_lp_vs_D": rho_all,
        "C6_nonmarker_pass": [x for x, v in per.items() if v["pass"] and (ann.get(x) or 1.0) <= cr["C6_marker_pct_max"]],
        "C7_seeds_low": c7,
    }
    res["C3_panel"] = res["C3_panel_n_ok"] >= cr["C3_panel_min_genes"]
    res["C8"] = bool(c8a["pass"] and c8b["pass"])
    # 判定規則（plan.md §6.4）。提案であり、確定は PI（G3）
    if not res["C8"]:
        v = "判定不能（較正または交絡の除去が不合格: 移動距離をそろえても交絡が残る／検定が較正されていない）"
    elif not (res["C3"] and res["C3_panel"]):
        v = "判定不能（恒常発現遺伝子または中立対照パネルが中立でない）"
    elif res["C1_C2_C9"]:
        v = "支持（陽性対照 2 以上が null・恒常発現遺伝子・比較遺伝子群のすべてを上回る）"
    elif res["C1_null_only"]:
        v = "否定（null からは分離するが、恒常発現遺伝子または発現をそろえた比較遺伝子群と区別できない）"
    else:
        v = "否定（移動距離をそろえた null から陽性対照が分離しない）"
    res["verdict_proposal"] = v
    return res


# ---------------------------------------------------------------- 提供者単位・半分割
def donor_level(c: pd.DataFrame, R: dict, dl: dict) -> pd.DataFrame:
    x = c.groupby(["gene", "seed", "donor_id"]).r.agg(["mean", "size"]).reset_index()
    x = x[x["size"] >= dl["min_cells_per_donor"]]
    xb = x.groupby(["gene", "donor_id"])["mean"].mean().reset_index()
    rows = []
    for g, s in xb.groupby("gene"):
        v = s["mean"].to_numpy()
        w = stats.wilcoxon(v, alternative="greater") if len(v) >= 5 and np.any(v != 0) else None
        rows.append({"gene": g, "n_donors": len(v), "frac_donor_pos": float((v > 0).mean()),
                     "wilcoxon_p_greater": (float(w.pvalue) if w is not None else np.nan)})
    d = pd.DataFrame(rows)
    nulls = d[d.gene.isin(R["negative"])].frac_donor_pos.to_numpy()
    d["frac_donor_pos_null_pct"] = [null_pct(v, nulls[d[d.gene.isin(R["negative"])].gene.to_numpy() != g])
                                    for g, v in zip(d.gene, d.frac_donor_pos)]
    if dl.get("mixed_model"):
        import statsmodels.formula.api as smf
        mm = []
        for g, s in c.groupby("gene"):
            s = s.dropna(subset=["r"]).copy()
            s["seed_c"] = s.seed.astype(str)
            try:
                f = smf.mixedlm("r ~ C(seed_c, Sum)", s, groups=s["donor_id"]).fit(reml=True, method="lbfgs")
                vd, ve = float(f.cov_re.iloc[0, 0]), float(f.scale)
                mm.append({"gene": g, "mm_intercept": float(f.params["Intercept"]), "mm_se": float(f.bse["Intercept"]),
                           "mm_z": float(f.params["Intercept"] / f.bse["Intercept"]), "mm_var_donor": vd,
                           "mm_var_resid": ve, "mm_icc": vd / (vd + ve) if vd + ve > 0 else np.nan})
            except Exception as e:  # 収束しない遺伝子は記録して続ける
                mm.append({"gene": g, "mm_error": str(e)[:80]})
        d = d.merge(pd.DataFrame(mm), on="gene", how="left")
    return d


def split_half(c: pd.DataFrame, sh: dict, fixed18: list[str]) -> dict:
    rng = np.random.default_rng(sh["rng_seed"])
    out = {}
    for s, cs in c.groupby("seed"):
        donors = np.array(sorted(cs.donor_id.unique()))
        rs_all, rs_18 = [], []
        for _ in range(sh["n_splits"]):
            a = set(rng.choice(donors, size=len(donors) // 2, replace=False))
            m = cs.donor_id.isin(a)
            A = cs[m].groupby("gene").r.mean()
            Bm = cs[~m].groupby("gene").r.mean()
            j = A.index.intersection(Bm.index)
            rs_all.append(stats.spearmanr(A[j], Bm[j]).statistic)
            j18 = [x for x in j if x in fixed18]
            rs_18.append(stats.spearmanr(A[j18], Bm[j18]).statistic)
        out[int(s)] = {"all": {"mean": float(np.mean(rs_all)), "min": float(np.min(rs_all)), "max": float(np.max(rs_all)),
                               "spearman_brown": float(2 * np.mean(rs_all) / (1 + np.mean(rs_all))), "n_genes": int(len(j))},
                       "fixed18": {"mean": float(np.mean(rs_18)), "min": float(np.min(rs_18)), "max": float(np.max(rs_18)),
                                   "spearman_brown": float(2 * np.mean(rs_18) / (1 + np.mean(rs_18))), "n_genes": len(j18)},
                       "n_donors": int(len(donors)), "n_donors_per_half": [len(donors) // 2, len(donors) - len(donors) // 2],
                       "n_splits": sh["n_splits"]}
    return out


def lodo(c: pd.DataFrame, genes: list[str], p: dict, cr: dict) -> dict:
    out = {}
    for x in genes:
        changed = []
        for dnr in sorted(c.donor_id.unique()):
            cc = c[c.donor_id != dnr]
            n_hi = 0
            for s, cs in cc.groupby("seed"):
                nulls = {h: hd for h, hd in cs[cs.role == "negative"].groupby("gene")}
                tg = cs[cs.gene == x]
                if len(tg) == 0:
                    continue
                r = perm_test(tg, nulls, None, p)
                n_hi += (not np.isnan(r["q"])) and r["q"] >= cr["q_hi"]
            if int(n_hi) < cr["min_seeds"]:
                changed.append(dnr)
        out[x] = {"n_donors_flip_C1": len(changed), "donors": changed}
    return out


# ---------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--cells", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--no-annotation", action="store_true", help="合成データ検査用（ddet_pct を読まない）")
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    A = cfg["analysis"]
    E = ROOT / "experiments" / cfg["experiment_id"]
    res_dir = E / cfg["output"]["results_dir"]
    out = Path(a.out) if a.out else res_dir / "eval"
    out.mkdir(parents=True, exist_ok=True)
    c = load_cells(cfg, Path(a.cells) if a.cells else res_dir / "cells.csv")
    R = roles(cfg)
    comp = comparator_map(cfg)
    ann = {}
    if not a.no_annotation:
        name2id = dict(zip(pd.read_csv(ROOT / "outputs/literature/2026-09-26-negative-controls/genome_stats.csv").symbol,
                           pd.read_csv(ROOT / "outputs/literature/2026-09-26-negative-controls/genome_stats.csv").ensembl_id))
        for r in json.loads((ROOT / "experiments/EXP-20260926-005/work/gene_annotation.json").read_text())["genes"]:
            ann[name2id.get(r["gene"], r["ensembl"])] = r["ddet_pct"]
    role_of = {g: k for k, v in R.items() for g in v}
    c["role"] = c.gene.map(role_of).fillna("other")
    seeds = sorted(c.seed.unique().tolist())
    p = A["pos_matched_permutation"]


    rows = []
    for s in seeds:
        cs = c[c.seed == s]
        nulls = {h: hd for h, hd in cs[cs.role == "negative"].groupby("gene")}
        a2 = a2_like(cs)
        wb = {}
        for g, tg in cs.groupby("gene"):
            ref = cs[(cs.role == "negative") & (cs.gene != g)]
            wb[g] = within_bin_score(tg, ref, A["within_bin_percentile"]["length_bins"], A["within_bin_percentile"]["pos_bins"])
        for g, tg in cs.groupby("gene"):
            r = perm_test(tg, nulls, g, p)
            rows.append({"seed": s, "gene": g, "role": role_of.get(g, "other"), "n_cells": len(tg),
                         "lp_med": float(tg.lp.median()), "relpos_med": float(tg.relpos.median()),
                         **r, "W": wb[g], "a2": a2[g]})
    g = pd.DataFrame(rows)
    for col, nc in (("W", "q_W"), ("a2", "q_a2")):
        g[nc] = [null_pct(v, g[(g.seed == s) & (g.role == "negative") & (g.gene != gg)][col].to_numpy())
                 for s, gg, v in zip(g.seed, g.gene, g[col])]
    g["D_W"] = g["W"]
    g.to_csv(out / "gene_seed.csv", index=False)

    cr = A["criteria"]
    primary = assess(g, "q", "D_raw", R, comp, cr, ann, seeds)   # 遺伝子間の比較（C2・C9・C8b）は共通の補正式の D_raw で行う
    secondary = assess(g, "q_W", "D_W", R, comp, cr, ann, seeds)
    cont = assess(g, "q_a2", "a2", R, comp, cr, ann, seeds)

    c["r"] = cell_residual(c, A["cell_residual"]["match_k"])
    dl = donor_level(c, R, A["donor_level"])
    dl.to_csv(out / "donor_level.csv", index=False)
    fixed18 = R["positive"] + R["housekeeping"] + R["negctrl_panel"] + R.get("reverse_direction", [])
    sh = split_half(c, A["split_half"], fixed18)
    (out / "split_half.json").write_text(json.dumps(sh, ensure_ascii=False, indent=1) + "\n")
    lo = lodo(c, [x for x, v in primary["per_positive"].items() if v["pass"]], p, cr) if A["donor_level"].get("lodo") else {}

    crit = {"primary_pos_matched_permutation": primary, "secondary_within_bin_percentile": secondary,
            "reference_a2_length_only": cont, "lodo_primary_pass_genes": lo,
            "concordance_primary_vs_secondary": {
                "pass_primary": [x for x, v in primary["per_positive"].items() if v["pass"]],
                "pass_secondary": [x for x, v in secondary["per_positive"].items() if v["pass"]]},
            "n_cells_min": int(g.n_cells.min()), "n_eligible_min_fixed": int(g[g.role != "negative"].n_elig.min()),
            "seeds": seeds, "note": "判定は primary のみ（plan.md §6.4）。verdict_proposal は提案で、確定は PI（G3）"}
    (out / "criteria.json").write_text(json.dumps(crit, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps({k: primary[k] for k in ["n_positive_pass", "n_positive_C1_only", "C3", "C3_panel", "C8",
                                              "verdict_proposal"]}, ensure_ascii=False, default=float))
    print("C8a", primary["C8a"], "C8b", {s: (round(v["lp"], 2), round(v["relpos"], 2)) for s, v in primary["C8b"]["rho_by_seed"].items()})


if __name__ == "__main__":
    main()

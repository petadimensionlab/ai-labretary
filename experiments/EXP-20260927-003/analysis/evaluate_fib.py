"""EXP-20260927-003 評価（事前登録。plan.md §6。摂動の実行前に sha256 を固定する）。

主検定は EXP-20260927-001 §6.1 と同じ「移動距離をそろえた並べ替え検定」。
  対象遺伝子 g・seed s の細胞 i（値 y_i、元の位置 p_i、入力の長さ L_i）について、
  距離 d(i,j) = max(|log2(1+p_i) − log2(1+p_j)|, |log2 L_i − log2 L_j|)。
  g の細胞の 80% 以上に d ≤ 0.5 の null 細胞がある null 遺伝子 h を「使える null 遺伝子」とする（20 未満なら NaN）。
  偏り補正 f（対象遺伝子を除く同じ seed の null 細胞で当てはめた 2 次式）を引いた遺伝子平均 D_raw を、
  使える null 遺伝子の同じ値と比べる（遺伝子ラベルの並べ替え）。q = D_raw が対象より小さい null 遺伝子の割合。

本実験で加えたもの
  - 状態の対照: MYOF（筋線維芽細胞）重心への接近度 shift_myof についても同じ検定を行い、
    (d) D_raw(goal) > D_raw(MYOF) を合格条件に入れる（「線維化方向へ動いただけ」でないことの確認）。
  - 符号の条件: 合格には D_raw > 0 かつ 補正前の平均 shift_goal > 0 を求める（guardrails 2026-09-26）。
  - パイプラインの陽性対照: ACTA2 / TAGLN が MYOF 重心へ近づくか（C11）。

実行: .venv/bin/python experiments/EXP-20260927-003/analysis/evaluate_fib.py --config experiments/EXP-20260927-003/config.yaml
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
# 盲検（解除されていない遺伝子）。LTBR・CDKN1A は HYP-20260925-02 の対象で、PI の解除を plan.md §11 に記録する
BLIND_REMAINING = {"ENSG00000141510", "ENSG00000118785", "ENSG00000148346"}   # TP53, SPP1, LCN2


def rel(p: str) -> Path:
    q = Path(p)
    return q if q.is_absolute() else ROOT / q


# ---------------------------------------------------------------- 入力
def cos(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return (a @ b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b) + 1e-12)


def add_extra_shifts(cfg: dict, c: pd.DataFrame, res_dir: Path) -> pd.DataFrame:
    """保存した埋め込み（capture.embeddings）と追加の重心から、重心ごとの接近度を列に足す。
    主判定に使うのは state_control.centroid（MYOF 重心）だけで、他は参考（感度分析）。"""
    names = list(cfg["extra_centroids"]["sets"])
    vals = {n: np.full(len(c), np.nan) for n in names}
    for (run, gene), part in c.groupby(["run", "gene"], sort=False):
        seed = int(str(run).split("-")[1])
        npz = rel(cfg["extra_centroids"]["output"].format(seed=seed))
        if not npz.exists():
            raise SystemExit(f"{npz} がありません。analysis/extra_centroids.py を先に実行してください")
        z = np.load(npz)
        gd = res_dir / str(run) / "genes" / str(gene) / "emb"
        o, p = np.load(gd / "orig.npy").astype(np.float64), np.load(gd / "pert.npy").astype(np.float64)
        if len(o) != len(part):
            raise SystemExit(f"{gd}: 埋め込み {len(o)} 細胞 ≠ cells.csv {len(part)} 細胞")
        for n in names:
            if n in z.files:
                cen = z[n].astype(np.float64)
                vals[n][part.index.to_numpy()] = cos(p, cen) - cos(o, cen)
    for n in names:
        c[f"shift_c:{n}"] = vals[n]
    c["shift_myof"] = c[f"shift_c:{cfg['analysis']['state_control']['centroid']}"]
    return c


def load_cells(cfg: dict, path: Path, res_dir: Path) -> pd.DataFrame:
    c = pd.read_csv(path, dtype={"donor_id": str})
    assert not set(c.gene) & BLIND_REMAINING, "解除されていない盲検の遺伝子が結果に含まれています"
    c["seed"] = c.run.astype(str).str.extract(r"run-(\d+)", expand=False).astype(int)
    if "pos" not in c.columns:
        from datasets import load_from_disk
        from geneformer import TOKEN_DICTIONARY_FILE
        tokd = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
        ds = load_from_disk(str(rel(cfg["data"]["tokenized_dataset"])))
        uids = sorted(c.cell_uid.unique())
        ids = dict(zip(uids, ds.select(uids)["input_ids"]))
        c["pos"] = [list(ids[u]).index(tokd[g]) for u, g in zip(c.cell_uid, c.gene)]
    c["lp"] = np.log2(1 + c.pos.astype(float))
    c["ll"] = np.log2(c.length.astype(float))
    c["relpos"] = c.pos / c.length
    c = add_extra_shifts(cfg, c, res_dir)
    c["y"] = c[cfg["analysis"]["base_metric"]]
    return c


def roles(cfg: dict) -> dict[str, list[str]]:
    ctl = cfg["controls"]
    r = {"target": list(cfg["perturbation"]["genes_to_perturb"]),
         "positive": list(ctl["positive"]),
         "housekeeping": list(ctl["negative_fixed"]),
         "negative": list(ctl["negative_sampler"]["genes"])}
    for grp in ctl.get("extra_fixed", []) or []:
        r[grp["role"]] = list(grp["genes"])
    return r


def comparator_map(cfg: dict) -> dict[str, list[str]]:
    """対象遺伝子 → 比較遺伝子 5 つ（work/comparators.csv の ref 列で対応づける）。"""
    E = rel(cfg["output"]["experiment_dir"])
    t = pd.read_csv(E / "work/comparators.csv")
    g = pd.read_csv(E / "work/gene_stats_fib.csv")
    name2id = dict(zip(g.symbol, g.ensembl))
    out = {name2id[ref]: s.ensembl.tolist() for ref, s in t.groupby("ref", sort=False)}
    assert all(len(v) == 5 for v in out.values()), "比較遺伝子が 5 つでない対象があります"
    return out


# ---------------------------------------------------------------- 主検定
def dist(a_lp, a_ll, b_lp, b_ll):
    return np.maximum(np.abs(a_lp[:, None] - b_lp[None, :]), np.abs(a_ll[:, None] - b_ll[None, :]))


def design(lp, ll):
    lp, ll = np.asarray(lp, float), np.asarray(ll, float)
    return np.column_stack([np.ones_like(lp), lp, ll, lp ** 2, ll ** 2, lp * ll])


def fit_bias(ref: pd.DataFrame, col: str):
    beta, *_ = np.linalg.lstsq(design(ref.lp, ref.ll), ref[col].to_numpy(), rcond=None)
    return lambda lp, ll: design(lp, ll) @ beta


def perm_test(tg: pd.DataFrame, nulls: dict[str, pd.DataFrame], exclude: str | None, p: dict,
              col: str = "y") -> dict:
    tol, frac = p["tolerance_log2"], p["min_cover_frac"]
    fb = fit_bias(pd.concat([hd for h, hd in nulls.items() if h != exclude]), col)
    Dg = float((tg[col].to_numpy() - fb(tg.lp, tg.ll)).mean())
    Dh = []
    for h, hd in nulls.items():
        if h == exclude:
            continue
        D = dist(tg.lp.to_numpy(), tg.ll.to_numpy(), hd.lp.to_numpy(), hd.ll.to_numpy())
        if (D <= tol).any(axis=1).mean() < frac:
            continue
        Dh.append(float((hd[col].to_numpy() - fb(hd.lp, hd.ll)).mean()))
    Dh = np.array(Dh)
    base = {"T": float(tg[col].mean()), "D_raw": Dg, "n_elig": int(len(Dh))}
    if len(Dh) < p["min_eligible_null_genes"]:
        return {**base, "q": np.nan, "D": np.nan, "null_mean": np.nan, "null_sd": np.nan,
                "null_frac_negative": np.nan}
    return {**base, "q": float((Dh < Dg).mean()), "D": float(Dg - Dh.mean()),
            "null_mean": float(Dh.mean()), "null_sd": float(Dh.std(ddof=1)),
            "null_frac_negative": float((Dh < 0).mean())}


def cell_residual(c: pd.DataFrame, k: int) -> np.ndarray:
    r = np.full(len(c), np.nan)
    for s, cs in c.groupby("seed"):
        nl = cs[cs.role == "negative"]
        for g, tg in cs.groupby("gene"):
            ref = nl[nl.gene != g]
            fb = fit_bias(ref, "y")
            D = dist(tg.lp.to_numpy(), tg.ll.to_numpy(), ref.lp.to_numpy(), ref.ll.to_numpy())
            kk = min(k, D.shape[1] - 1)
            nn = np.argpartition(D, kk, axis=1)[:, :kk]
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


def null_pct(x: float, null: np.ndarray) -> float:
    null = null[~np.isnan(null)]
    return float((null < x).mean()) if len(null) and not np.isnan(x) else np.nan


# ---------------------------------------------------------------- 判定
def assess(g: pd.DataFrame, qcol: str, dcol: str, R: dict, comp: dict, cr: dict, seeds: list[int]) -> dict:
    Q = g.pivot(index="gene", columns="seed", values=qcol)
    Dm = g.pivot(index="gene", columns="seed", values=dcol)
    Dmy = g.pivot(index="gene", columns="seed", values="D_raw_myof")
    Tm = g.pivot(index="gene", columns="seed", values="T")
    need = cr["min_seeds"]
    sign_req = bool(cr.get("require_positive_sign", True))

    def n_hi(x):
        if x not in Q.index:
            return 0
        ok = (Q.loc[x] >= cr["q_hi"])
        if sign_req:
            ok &= (Dm.loc[x] > 0) & (Tm.loc[x] > 0)
        return int(ok.sum())

    def in_rng(x, lo, hi):
        return int(((Q.loc[x] >= lo) & (Q.loc[x] <= hi)).sum()) if x in Q.index else 0  # NaN は範囲外

    per = {}
    for x in R["target"] + R.get("target_ref", []):
        c1 = n_hi(x) >= need
        wins = int(sum((Dm.loc[x, s] > Dm.loc[h, s]) for s in seeds for h in R["housekeeping"]
                       if not (np.isnan(Dm.loc[x, s]) or np.isnan(Dm.loc[h, s]))))
        beat = {int(s): int(sum(Dm.loc[x, s] > Dm.loc[cg, s] for cg in comp.get(x, [])
                                if not (np.isnan(Dm.loc[x, s]) or np.isnan(Dm.loc[cg, s])))) for s in seeds}
        c9 = sum(v >= cr["C9_min_beaten"] for v in beat.values()) >= cr["C9_min_seeds"]
        spec = {int(s): (None if (np.isnan(Dm.loc[x, s]) or np.isnan(Dmy.loc[x, s]))
                         else bool(Dm.loc[x, s] > Dmy.loc[x, s])) for s in seeds}
        c10 = sum(bool(v) for v in spec.values()) >= cr["C10_min_seeds"]
        per[x] = {"seeds_pass_a": n_hi(x), "a_null_separation": c1,
                  "b_wins_vs_housekeeping_of_9": wins, "b_housekeeping": wins >= cr["C2_min_wins"],
                  "c_beaten_by_seed": beat, "c_comparator": c9,
                  "d_goal_gt_myof_by_seed": spec, "d_state_specificity": c10,
                  "pass": bool(c1 and wins >= cr["C2_min_wins"] and c9 and c10),
                  "q_by_seed": {int(s): (None if np.isnan(Q.loc[x, s]) else round(float(Q.loc[x, s]), 3))
                                for s in seeds},
                  "D_raw_by_seed": {int(s): (None if np.isnan(Dm.loc[x, s]) else float(Dm.loc[x, s]))
                                    for s in seeds},
                  "D_raw_myof_by_seed": {int(s): (None if np.isnan(Dmy.loc[x, s]) else float(Dmy.loc[x, s]))
                                         for s in seeds},
                  "mean_shift_goal_by_seed": {int(s): (None if np.isnan(Tm.loc[x, s]) else float(Tm.loc[x, s]))
                                              for s in seeds}}
    lo, hi = cr["C3_range"]
    c3 = {x: in_rng(x, lo, hi) for x in R["housekeeping"]}
    cp = {x: in_rng(x, lo, hi) for x in R["negctrl_panel"]}
    main = [x for x in R["target"]]
    n_pass = sum(per[x]["pass"] for x in main)
    n_a = sum(per[x]["a_null_separation"] for x in main)

    qn = g[g.role == "negative"][qcol].to_numpy()
    qn = qn[~np.isnan(qn)]
    c8a = {"n": int(len(qn)), "frac_hi": float((qn >= 0.90).mean()) if len(qn) else np.nan,
           "frac_lo": float((qn < 0.10).mean()) if len(qn) else np.nan}
    c8a["pass"] = bool(len(qn) >= cr["C8a_min_n"] and c8a["frac_hi"] <= cr["C8a_tail_frac_max"]
                       and c8a["frac_lo"] <= cr["C8a_tail_frac_max"])
    rho = {}
    for s in seeds:
        n = g[(g.role == "negative") & (g.seed == s)].dropna(subset=[dcol])
        rho[int(s)] = {"lp": float(stats.spearmanr(n.lp_med, n[dcol]).statistic),
                       "relpos": float(stats.spearmanr(n.relpos_med, n[dcol]).statistic),
                       "ll": float(stats.spearmanr(n.ll_med, n[dcol]).statistic), "n": int(len(n))}
    ok = sum(all(abs(v[k]) <= cr["C8b_abs_rho_max"] for k in ("lp", "relpos", "ll")) for v in rho.values())
    c8b = {"rho_by_seed": rho, "pass": ok >= cr["C8b_min_seeds"]}

    # C11 パイプラインの陽性対照: ACTA2 / TAGLN が MYOF 重心へ近づくか
    Qmy = g.pivot(index="gene", columns="seed", values="q_myof")
    c11 = {}
    for x in cr["C11_positive_control"]:
        ok_s = int(((Qmy.loc[x] >= cr["q_hi"]) & (Dmy.loc[x] > 0)).sum()) if x in Qmy.index else 0
        c11[x] = {"seeds_pass": ok_s, "pass": ok_s >= need,
                  "q_myof_by_seed": {int(s): (None if np.isnan(Qmy.loc[x, s]) else round(float(Qmy.loc[x, s]), 3))
                                     for s in seeds}}
    c11_pass = sum(v["pass"] for v in c11.values()) >= cr["C11_min_genes"]

    res = {
        "per_gene": per, "n_target_pass": n_pass, "n_target_a_only": n_a,
        "C3_seeds_in_range": c3, "C3": all(v >= need for v in c3.values()),
        "C3_panel_seeds_in_range": cp, "C3_panel_n_ok": sum(v >= need for v in cp.values()),
        "C8a_calibration": c8a, "C8b_confounding": c8b,
        "C11_pipeline_positive_control": c11, "C11": c11_pass,
        "null_frac_negative_by_seed": {int(s): float(g[(g.role == "negative") & (g.seed == s)]
                                                     .null_frac_negative.mean()) for s in seeds},
    }
    res["C3_panel"] = res["C3_panel_n_ok"] >= cr["C3_panel_min_genes"]
    res["C8"] = bool(c8a["pass"] and c8b["pass"])
    # 判定規則（plan.md §6.4）。提案であり、確定は PI（G3）
    if not res["C8"]:
        v = "判定不能（較正または交絡の除去が不合格: 移動距離をそろえても交絡が残る／検定が較正されていない）"
    elif not (res["C3"] and res["C3_panel"]):
        v = "判定不能（恒常発現遺伝子または中立対照パネルが中立でない）"
    elif not res["C11"]:
        v = "判定不能（パイプラインの陽性対照が働かない: 重心へ近づく効果を検出できることが示せていない）"
    elif n_pass >= cr["C1_min_genes"]:
        v = f"支持（主対象 3 中 {n_pass} が null・恒常発現遺伝子・比較遺伝子群・状態の対照のすべてを上回る）"
    elif n_pass == 1:
        v = "一部支持（主対象 1 遺伝子のみが全条件を満たす。仮説全体は支持としない）"
    elif n_a >= 1:
        v = "否定（null からは分離するが、恒常発現遺伝子・比較遺伝子群・状態の対照と区別できない）"
    else:
        v = "否定（移動距離をそろえた null から主対象が分離しない）"
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
    nl = d[d.gene.isin(R["negative"])]
    d["frac_donor_pos_null_pct"] = [null_pct(v, nl[nl.gene != g].frac_donor_pos.to_numpy())
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
            except Exception as e:
                mm.append({"gene": g, "mm_error": str(e)[:80]})
        d = d.merge(pd.DataFrame(mm), on="gene", how="left")
    return d


def split_half(c: pd.DataFrame, sh: dict, fixed: list[str]) -> dict:
    rng = np.random.default_rng(sh["rng_seed"])
    out = {}
    for s, cs in c.groupby("seed"):
        donors = np.array(sorted(cs.donor_id.unique()))
        rs_all, rs_f = [], []
        for _ in range(sh["n_splits"]):
            a = set(rng.choice(donors, size=len(donors) // 2, replace=False))
            m = cs.donor_id.isin(a)
            A, B = cs[m].groupby("gene").r.mean(), cs[~m].groupby("gene").r.mean()
            j = A.index.intersection(B.index)
            rs_all.append(stats.spearmanr(A[j], B[j]).statistic)
            jf = [x for x in j if x in fixed]
            rs_f.append(stats.spearmanr(A[jf], B[jf]).statistic if len(jf) >= 5 else np.nan)
        out[int(s)] = {"all": {"mean": float(np.mean(rs_all)), "min": float(np.min(rs_all)),
                               "max": float(np.max(rs_all)),
                               "spearman_brown": float(2 * np.mean(rs_all) / (1 + np.mean(rs_all))),
                               "n_genes": int(len(j))},
                       "fixed": {"mean": float(np.nanmean(rs_f)), "n_genes": len(jf)},
                       "n_donors": int(len(donors)),
                       "n_donors_per_half": [len(donors) // 2, len(donors) - len(donors) // 2],
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
        out[x] = {"n_donors_flip_a": len(changed), "donors": changed}
    return out


# ---------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--cells", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    A = cfg["analysis"]
    E = rel(cfg["output"]["experiment_dir"])
    res_dir = E / cfg["output"]["results_dir"]
    out = Path(a.out) if a.out else res_dir / "eval"
    out.mkdir(parents=True, exist_ok=True)
    c = load_cells(cfg, Path(a.cells) if a.cells else res_dir / "cells.csv", res_dir)
    R = roles(cfg)
    comp = comparator_map(cfg)
    role_of = {g: k for k, v in R.items() for g in v}
    c["role"] = c.gene.map(role_of).fillna("other")
    seeds = sorted(c.seed.unique().tolist())
    p = A["pos_matched_permutation"]

    rows = []
    for s in seeds:
        cs = c[c.seed == s]
        nulls = {h: hd for h, hd in cs[cs.role == "negative"].groupby("gene")}
        wb = {}
        for g, tg in cs.groupby("gene"):
            ref = cs[(cs.role == "negative") & (cs.gene != g)]
            wb[g] = within_bin_score(tg, ref, A["within_bin_percentile"]["length_bins"],
                                     A["within_bin_percentile"]["pos_bins"])
        for g, tg in cs.groupby("gene"):
            r = perm_test(tg, nulls, g, p)
            rmy = perm_test(tg, nulls, g, p, col="shift_myof")
            rows.append({"seed": s, "gene": g, "role": role_of.get(g, "other"), "n_cells": len(tg),
                         "n_donors": tg.donor_id.nunique(),
                         "lp_med": float(tg.lp.median()), "ll_med": float(tg.ll.median()),
                         "relpos_med": float(tg.relpos.median()), **r,
                         "q_myof": rmy["q"], "D_raw_myof": rmy["D_raw"], "T_myof": rmy["T"],
                         "W": wb[g]})
    g = pd.DataFrame(rows)
    g["q_W"] = [null_pct(v, g[(g.seed == s) & (g.role == "negative") & (g.gene != gg)]["W"].to_numpy())
                for s, gg, v in zip(g.seed, g.gene, g["W"])]
    g["sign_D_raw"] = np.sign(g.D_raw)
    g["sign_T"] = np.sign(g["T"])
    g.to_csv(out / "gene_seed.csv", index=False)

    cr = A["criteria"]
    primary = assess(g, "q", "D_raw", R, comp, cr, seeds)
    secondary = assess(g, "q_W", "W", R, comp, cr, seeds)

    c["r"] = cell_residual(c, A["cell_residual"]["match_k"])
    dl = donor_level(c, R, A["donor_level"])
    dl.to_csv(out / "donor_level.csv", index=False)
    fixed = R["target"] + R.get("target_ref", []) + R["housekeeping"] + R["negctrl_panel"] + R["positive"]
    sh = split_half(c, A["split_half"], fixed)
    (out / "split_half.json").write_text(json.dumps(sh, ensure_ascii=False, indent=1) + "\n")
    lo = lodo(c, [x for x in R["target"] if primary["per_gene"][x]["pass"]], p, cr) \
        if A["donor_level"].get("lodo") else {}

    crit = {"primary_pos_matched_permutation": primary,
            "secondary_within_bin_percentile": secondary,
            "lodo_primary_pass_genes": lo,
            "concordance_primary_vs_secondary": {
                "pass_primary": [x for x in R["target"] if primary["per_gene"][x]["pass"]],
                "pass_secondary": [x for x in R["target"] if secondary["per_gene"][x]["pass"]]},
            "n_cells_min": int(g.n_cells.min()),
            "n_eligible_min_non_null": int(g[g.role != "negative"].n_elig.min()),
            "seeds": seeds,
            "note": "判定は primary のみ（plan.md §6.4）。verdict_proposal は提案で、確定は PI（G3）。"
                    "D_raw・T は符号付き（正 = 目標状態へ近づく）"}
    (out / "criteria.json").write_text(json.dumps(crit, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps({k: primary[k] for k in ["n_target_pass", "n_target_a_only", "C3", "C3_panel", "C8",
                                              "C11", "verdict_proposal"]}, ensure_ascii=False, default=float))


if __name__ == "__main__":
    main()

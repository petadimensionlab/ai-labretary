#!/usr/bin/env python3
"""
summary.csv から定量的な検証指標を計算する（validator エージェントの材料）。
出力: <results>/validation_metrics.json と標準出力の Markdown 表

--config を渡し、config に analysis セクションがあれば、事前登録した仮説の判定指標も計算する
（EXP-20260926-001 で追加。hyp_metrics() 参照）:
  - seed ごとに、対象・陽性対照の shift を同 seed のランダム遺伝子（null_role）と比べた経験的 p と null 内の分位
  - 対象遺伝子（fdr_family）内の BH-FDR、条件を満たした seed 数
  - null の正規近似による p（経験的 p の下限 1/(n_null+1) を補う参考値。事前登録の主指標ではない）
  - analysis.metrics（例: [shift_goal, shift_goal_minus_start]）の指標ごとに判定する。
    shift_goal_minus_start = goal 方向シフト − start 方向シフト（削除ドリフトを差し引く, EXP-20260926-001 G2 で登録）
  - analysis.fdr_p: norm なら primary 判定の FDR に p_norm を使う（既定 emp）
  - --alt-results があれば shift(start) / shift(alt run) の比と donor ブートストラップ 95% CI（cells.csv が必要）
    context_dependence は「比 ≥ ratio_threshold かつ CI 下限 > 1」に加え **shift_A > 0（goal 方向）** を要求する。
    両方の shift が負（goal から離れる）でも比は 1 を超えうるため（例: −3.2e-4 / −1.6e-4 = 1.96）、
    符号を見ない比は仮説の支持にならない（EXP-20260926-001 validator で修正, guardrails 2026-09-26）
  - 交絡チェックの発現頻度: n_detections があればそれ、なければ work/detection-<start>.json（start 細胞のうち
    遺伝子が入力に含まれる細胞数, --config が必要）、どちらもなければ n_cells。cell_selection では n_cells が上限で
    頭打ち（全遺伝子 100 なら定数 → 相関が NaN）になるため
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def shift_col(df):
    return "shift_goal_mean" if "shift_goal_mean" in df.columns else "shift"


METRICS = {  # 指標名 → (summary.csv の列から作る式, cells.csv の列から作る式)
    "shift_goal": (lambda d: d[shift_col(d)], lambda c: c["shift_goal"]),
    "shift_goal_minus_start": (lambda d: d["shift_goal_mean"] - d["shift_start_mean"],
                               lambda c: c["shift_goal"] - c["shift_start"]),
}


def hyp_metrics(df: pd.DataFrame, an: dict, targets: list[str], positives: list[str], col: str | None = None) -> dict:
    """事前登録の判定（HYP-20260925-01 形式）を seed ごとに計算する。"""
    col = col or shift_col(df)
    null_role = an.get("null_role", "negative")
    q_prim, q_pos, q_fal = an.get("null_quantile_primary", 0.95), an.get("positive_null_quantile", 0.90), \
        an.get("falsify_null_quantile", 0.80)
    rows = []
    for run, sub in df.groupby("run"):
        null = sub.loc[sub.role == null_role, col].to_numpy()
        if len(null) < 5:
            continue
        mu, sd = null.mean(), null.std(ddof=1)
        for _, r in sub[sub.gene.isin(targets + positives)].iterrows():
            x = r[col]
            rows.append({"run": run, "gene": r.gene, "role": r.role, "shift": x,
                         "n_cells": r.get("n_cells"), "n_null": len(null),
                         "null_pct": float((null < x).mean()),
                         "p_emp": float((1 + (null >= x).sum()) / (len(null) + 1)),
                         "p_norm": float(stats.norm.sf((x - mu) / sd)) if sd > 0 else np.nan,
                         "above_q_primary": bool(x > np.quantile(null, q_prim)),
                         "above_q_positive": bool(x > np.quantile(null, q_pos)),
                         "le_q_falsify": bool(x <= np.quantile(null, q_fal)), "negative_direction": bool(x < 0),
                         "null_median": float(np.median(null)), "null_frac_negative": float((null < 0).mean())})
    t = pd.DataFrame(rows)
    if t.empty:
        return {"note": "null 遺伝子が 5 未満の run しかありません"}
    fam = targets if an.get("fdr_family", "targets") == "targets" else targets + positives
    for run, sub in t.groupby("run"):
        f = sub.gene.isin(fam)
        t.loc[sub.index[f], "fdr_emp"] = bh(sub.loc[f, "p_emp"])
        t.loc[sub.index[f], "fdr_norm"] = bh(sub.loc[f, "p_norm"])
    fdr_col = "fdr_norm" if an.get("fdr_p", "emp") == "norm" else "fdr_emp"
    t["primary_pass"] = t.above_q_primary & (t[fdr_col] < an.get("fdr_threshold", 0.05))
    per_gene = {}
    for g, sub in t.groupby("gene"):
        per_gene[g] = {"role": sub.role.iloc[0], "n_seeds": len(sub), "shift_mean": float(sub["shift"].mean()),
                       "shift_sd": float(sub["shift"].std(ddof=1)) if len(sub) > 1 else None,
                       "null_pct": [round(v, 3) for v in sub.null_pct],
                       "p_emp": [round(v, 4) for v in sub.p_emp], "fdr_emp": [round(v, 4) for v in sub.fdr_emp.fillna(np.nan)],
                       "p_norm": [round(v, 5) for v in sub.p_norm], "fdr_norm": [round(v, 5) for v in sub.fdr_norm.fillna(np.nan)],
                       "n_seeds_primary_pass": int(sub.primary_pass.sum()),
                       "n_seeds_above_q_positive": int(sub.above_q_positive.sum()),
                       "falsified_all_seeds_le_q80_or_negative": bool((sub.le_q_falsify | sub.negative_direction).all()),
                       # 反証条件の内訳（null 自体が負に偏る場合、向きの条件は情報を持たない → null_frac_negative を併記）
                       "falsify_all_seeds_le_q80": bool(sub.le_q_falsify.all()),
                       "falsify_all_seeds_negative": bool(sub.negative_direction.all()),
                       "null_frac_negative": [round(v, 3) for v in sub.null_frac_negative]}
    prim = an.get("primary_targets", targets)
    k = an.get("min_seeds_passing", 2)
    return {"per_gene": per_gene,
            "primary_supported": {g: per_gene[g]["n_seeds_primary_pass"] >= k for g in prim if g in per_gene},
            "positive_controls_pass": {g: per_gene[g]["n_seeds_above_q_positive"] >= k for g in positives if g in per_gene},
            "p_emp_floor": float(1 / (t.n_null.min() + 1)),
            "fdr_used": fdr_col,
            "note": "p_emp の下限は 1/(n_null+1)。BH は対象遺伝子内で seed ごと"}


def detection_counts(exp: Path, config: str | None, results: str) -> dict:
    """start 細胞のうち各遺伝子が入力に含まれる細胞数（run_isp.py start_detection のキャッシュ）。
    results が <dir>/<run 名> なら config の runs からその run の start_state を使う。"""
    if not config:
        return {}
    import yaml
    cfg = yaml.safe_load(open(config))
    start = cfg.get("states", {}).get("start_state")
    run = Path(results).name
    start = ((cfg.get("runs") or {}).get(run, {}).get("states", {}).get("start_state")) or start
    f = exp / "work" / f"detection-{start}.json"
    return json.loads(f.read_text()) if f.exists() else {}


def shift_ratio(cells_a: pd.DataFrame, cells_b: pd.DataFrame, genes: list[str], n_boot: int, seed: int,
                value: str = "shift_goal") -> dict:
    """shift(run A) / shift(run B)。細胞ごとの shift を seed 間で平均（同じ細胞が複数 seed に出るため）し、
    donor を run ごとに独立に復元抽出して比の 95% CI を出す。"""
    rng = np.random.default_rng(seed)
    out = {}
    for g in genes:
        per = []
        for c in (cells_a, cells_b):
            x = c[c.gene == g].groupby(["cell_uid", "donor_id"], as_index=False)[value].mean()
            per.append({d: v.to_numpy() for d, v in x.groupby("donor_id")[value]})
        if not per[0] or not per[1]:
            continue
        mean = lambda dd: np.concatenate(list(dd.values())).mean()
        est_a, est_b = mean(per[0]), mean(per[1])
        boots = []
        for _ in range(n_boot):
            bs = []
            for dd in per:
                keys = list(dd)
                pick = rng.choice(len(keys), size=len(keys), replace=True)
                bs.append(np.concatenate([dd[keys[i]] for i in pick]).mean())
            boots.append(bs[0] / bs[1] if bs[1] != 0 else np.nan)
        boots = np.array(boots)
        out[g] = {"shift_A": float(est_a), "shift_B": float(est_b),
                  "ratio": float(est_a / est_b) if est_b != 0 else None,
                  "ci95": [float(np.nanpercentile(boots, 2.5)), float(np.nanpercentile(boots, 97.5))],
                  "n_donors_A": len(per[0]), "n_donors_B": len(per[1]),
                  "denominator_nonpositive_frac": float(np.mean([b <= 0 for b in boots if np.isfinite(b)] or [np.nan])),
                  "sign_A": "+" if est_a > 0 else "-", "sign_B": "+" if est_b > 0 else "-",
                  "note": "比は shift_A > 0 かつ shift_B > 0 のときだけ「goal 方向シフトの倍率」と読める"}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", required=True, help="experiments/EXP-...")
    ap.add_argument("--results", default="results", help="results ディレクトリ名（results-mock 等）")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--config", default=None, help="analysis セクションがあれば仮説の判定指標も計算")
    ap.add_argument("--alt-results", default=None, help="比較する run（例: results/tls_low）")
    a = ap.parse_args()
    res = Path(a.exp) / a.results
    df = pd.read_csv(res / "summary.csv")
    mock = (res / "MOCK_DATA_DO_NOT_REPORT").exists()
    m = {"mock_data": mock, "n_runs": df["run"].nunique(), "n_genes": df["gene"].nunique()}

    # 1 遺伝子ずつ perturb した場合、stats CSV に pval/fdr/n_detections がない（Shift_to_goal_end のみ）
    spec = {"shift": (shift_col(df), "mean"), "sd": (shift_col(df), "std")}
    if "fdr" in df:
        spec["fdr"] = ("fdr", "max")
    spec["n_det"] = ("n_detections", "mean") if "n_detections" in df else ("n_cells", "mean") if "n_cells" in df else None
    spec = {k: v for k, v in spec.items() if v}
    agg = df.groupby(["gene", "role"], as_index=False).agg(**spec)
    for c in ("fdr", "n_det"):
        if c not in agg:
            agg[c] = np.nan
    n_det_src = spec["n_det"][0] if "n_det" in spec else None
    det = {} if "n_detections" in df else detection_counts(Path(a.exp), a.config, a.results)
    if det:  # cell_selection では n_cells は上限（100 など）で頭打ち → 発現頻度には start の検出細胞数を使う
        agg["n_det"] = agg["gene"].map(det)
        n_det_src = "work/detection-<start>.json"
    agg["pct_rank"] = agg["shift"].rank(pct=True, ascending=False)  # 0 に近いほど上位

    # 1. ポジティブコントロール
    pos = agg[agg.role == "positive"]
    fdr_avail = bool(agg.fdr.notna().any())  # 1 遺伝子ずつの ISP では fdr がない → FDR 系の判定は N/A
    m["positive"] = {"n": len(pos),
                     "median_pct_rank": float(pos.pct_rank.median()) if len(pos) else None,
                     "frac_fdr_sig": float((pos.fdr < a.fdr).mean()) if len(pos) else None,
                     "pass": bool(len(pos) and pos.pct_rank.median() < 0.2)}

    # 2. ネガティブコントロール / null
    neg, tgt = agg[agg.role == "negative"], agg[agg.role.isin(["target", "screen"])]
    if len(neg) >= 5 and len(tgt):
        u = stats.mannwhitneyu(tgt["shift"], neg["shift"], alternative="greater")
        thr = float(np.quantile(neg["shift"], 0.95))
        m["negative"] = {"n": len(neg), "null_q95": thr, "mwu_p_target_gt_null": float(u.pvalue),
                         "n_targets_above_null_q95": int((tgt["shift"] > thr).sum()),
                         "neg_frac_fdr_sig": float((neg.fdr < a.fdr).mean()) if fdr_avail else None,
                         "pass": bool((neg.fdr < a.fdr).mean() <= 0.1) if fdr_avail else None}
        if not fdr_avail:
            m["negative"]["note"] = "fdr 列がない（1 遺伝子ずつの ISP）ため偽陽性率は判定しない。hypothesis の null 分位を見る"
    else:
        m["negative"] = {"n": len(neg), "pass": False, "note": "ネガティブコントロールが 5 未満"}

    # 3. 頑健性: run 間 Spearman
    wide = df.pivot_table(index="gene", columns="run", values=shift_col(df))
    runs = list(wide.columns)
    rhos = [stats.spearmanr(wide[r1], wide[r2], nan_policy="omit").statistic
            for i, r1 in enumerate(runs) for r2 in runs[i + 1:]]
    m["robustness"] = {"pairwise_spearman": [round(float(r), 3) for r in rhos],
                       "min_rho": float(min(rhos)) if rhos else None,
                       "pass": bool(rhos and min(rhos) >= 0.7)}

    # 4. 交絡: 発現頻度との相関
    if agg["n_det"].notna().sum() > 10:
        r = stats.spearmanr(agg["n_det"], agg["shift"], nan_policy="omit").statistic
        m["confound_expression"] = {"spearman_ndet_vs_shift": float(r), "n_det_source": n_det_src,
                                    "n_genes": int(agg["n_det"].notna().sum()), "pass": bool(abs(r) < 0.3)}
    else:
        m["confound_expression"] = {"note": "発現頻度が得られない（または定数）", "pass": None}

    # 上位候補
    top = tgt.sort_values("shift", ascending=False).head(20)
    m["top_targets"] = top[["gene", "shift", "sd", "fdr"]].round(5).to_dict("records")

    if a.config:
        import yaml
        cfg = yaml.safe_load(open(a.config))
        an = cfg.get("analysis")
        if an:
            tg = list(cfg["perturbation"]["genes_to_perturb"])
            ps = list(cfg.get("controls", {}).get("positive", []))
            an = {**an, "fdr_threshold": a.fdr}
            metrics = an.get("metrics", ["shift_goal"])
            if a.alt_results:
                ca = pd.read_csv(res / "cells.csv", dtype={"donor_id": str})
                cb = pd.read_csv(Path(a.exp) / a.alt_results / "cells.csv", dtype={"donor_id": str})
            for name in metrics:
                f_sum, f_cell = METRICS[name]
                df["_metric"] = f_sum(df)
                h = hyp_metrics(df, an, tg, ps, col="_metric")
                if a.alt_results:
                    for c in (ca, cb):
                        c["_metric"] = f_cell(c)
                    h["shift_ratio"] = shift_ratio(ca, cb, tg + ps, an.get("bootstrap_n", 2000),
                                                   an.get("bootstrap_seed", 0), value="_metric")
                    thr = an.get("ratio_threshold", 1.2)
                    # 符号を見ない比は誤判定する（両方負でも比 > 1）。goal 方向（shift_A > 0）を必須にする
                    h["context_dependence"] = {g: bool(v["ratio"] is not None and v["shift_A"] > 0 and v["shift_B"] > 0
                                                       and v["ratio"] >= thr and v["ci95"][0] > 1)
                                               for g, v in h["shift_ratio"].items()}
                    h["context_dependence_rule"] = (f"ratio ≥ {thr} かつ CI95 下限 > 1 かつ shift_A > 0 かつ shift_B > 0 "
                                                    "（shift_B ≤ 0 で shift_A > 0 なら比は定義できないが方向は仮説に沿う → 別途記載）")
                m.setdefault("hypothesis", {})[name] = h
            df.drop(columns="_metric", inplace=True)

    (res / "validation_metrics.json").write_text(json.dumps(m, indent=2, ensure_ascii=False))
    print(f"# 定量チェック ({res})" + ("  ⚠️ MOCK DATA" if mock else ""))
    print("| 項目 | 判定 | 詳細 |\n|---|---|---|")
    for k in ("positive", "negative", "robustness", "confound_expression"):
        if k in m:
            d = {kk: vv for kk, vv in m[k].items() if kk != "pass"}
            verdict = "N/A" if m[k]["pass"] is None else "Pass" if m[k]["pass"] else "Fail"
            print(f"| {k} | {verdict} | {json.dumps(d, ensure_ascii=False)} |")
    print(f"\nsaved: {res / 'validation_metrics.json'}")


if __name__ == "__main__":
    main()

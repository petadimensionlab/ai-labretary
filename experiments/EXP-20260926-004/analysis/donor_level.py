#!/usr/bin/env python3
"""EXP-20260926-004 (c): donor 単位の指標（plan.md §C.1 に事前登録した定義）で対照の判定と信頼性を再評価する（CPU のみ）。
usage:
  .venv/bin/python experiments/EXP-20260926-004/analysis/donor_level.py \
      --results experiments/EXP-20260926-003/results-ft --out experiments/EXP-20260926-004/donor_level/ft [--logits] \
      [--check-gsm experiments/EXP-20260926-003/results-ft/eval/gene_seed_metrics.csv] [--donor-subset high]
出力: summary.json（判定・信頼性）, gene_seed_scores.csv（細胞単位・DL-med・DL-mean）, donor_means.csv, lodo.csv
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as cm  # noqa: E402

MIN_CELLS_DONOR = 3
N_SPLIT = 20


def donor_means(c: pd.DataFrame, m: str) -> pd.DataFrame:
    d = c.groupby(["run", "gene", "role", "donor_id"]).agg(x=(m, "mean"), n=(m, "size")).reset_index()
    return d[d.n >= MIN_CELLS_DONOR]


def scores(c: pd.DataFrame, m: str) -> pd.DataFrame:
    """遺伝子×seed の 3 種のスコア: cell（細胞平均, EXP-003 と同じ）, dl_med, dl_mean。"""
    cell = c.groupby(["run", "gene", "role"])[m].mean().rename("cell")
    d = donor_means(c, m)
    agg = d.groupby(["run", "gene", "role"]).x.agg(dl_med="median", dl_mean="mean", n_donors="size")
    return pd.concat([cell, agg], axis=1).reset_index()


def split_half(c: pd.DataFrame, m: str, donors_keep=None) -> dict:
    rng = np.random.default_rng(0)
    out = {"cell": [], "dl_med": [], "dl_mean": []}
    for run, s in c.groupby("run"):
        if donors_keep is not None:
            s = s[s.donor_id.isin(donors_keep)]
        donors = np.array(sorted(s.donor_id.unique()))
        for _ in range(N_SPLIT):
            half = set(rng.permutation(donors)[: len(donors) // 2])
            h = s.donor_id.isin(half)
            a, b = scores(s[h], m).set_index("gene"), scores(s[~h], m).set_index("gene")
            idx = a.index.intersection(b.index)
            for k in out:
                aa, bb = a.loc[idx, k], b.loc[idx, k]
                ok = aa.notna() & bb.notna()
                out[k].append(stats.spearmanr(aa[ok], bb[ok]).statistic)
    res = {}
    for k, v in out.items():
        r = float(np.nanmean(v))
        res[k] = {"rho": round(r, 3), "sd": round(float(np.nanstd(v)), 3), "spearman_brown": round(2 * r / (1 + r), 3)}
    return res


def sign_and_excess(c: pd.DataFrame, m: str) -> dict:
    d = donor_means(c, m)
    # null 遺伝子の同 seed・同 donor の平均
    nul = d[d.role == "negative"].groupby(["run", "donor_id"]).x.mean().rename("null_x")
    d = d.join(nul, on=["run", "donor_id"])
    xbar = d.groupby(["gene", "role", "donor_id"]).x.mean().reset_index()  # seed 平均
    frac = xbar.groupby(["gene", "role"]).x.apply(lambda v: float((v > 0).mean())).reset_index(name="frac_pos")
    nullfrac = frac[frac.role == "negative"].frac_pos.to_numpy()
    res = {}
    for gene, sym in {**cm.POS, **cm.HK}.items():
        if gene not in set(frac.gene):
            continue
        fp = float(frac[frac.gene == gene].frac_pos.iloc[0])
        dd = d[d.gene == gene].assign(y=lambda t: t.x - t.null_x).groupby("donor_id").y.mean()
        w = stats.wilcoxon(dd, alternative="greater") if len(dd) >= 5 else None
        res[sym] = {"frac_donors_pos": round(fp, 3), "null_pct_of_frac_pos": round(cm.null_pct(fp, nullfrac), 3),
                    "n_donors": int(len(dd)), "excess_frac_pos": round(float((dd > 0).mean()), 3),
                    "excess_median": float(dd.median()), "excess_wilcoxon_p_greater": None if w is None else float(w.pvalue)}
    res["null_frac_pos_median"] = round(float(np.median(nullfrac)), 3)
    return res


def c4c5(c: pd.DataFrame, g: pd.DataFrame, m: str) -> dict:
    """参考指標: C4ii（固定対照 5 遺伝子の seed 間順位一致）と C5（長さ・発現細胞数との相関）。EXP-002 と同じ定義。"""
    runs = sorted(g.run.unique())
    fixed = list(cm.POS) + list(cm.HK)
    piv = g[g.gene.isin(fixed)].pivot(index="gene", columns="run", values="cell")
    rr = [stats.spearmanr(piv[a], piv[b]).statistic for i, a in enumerate(runs) for b in runs[i + 1:]]
    L = c.groupby(["run", "gene"]).length.mean()
    v = c.groupby(["run", "gene"])[m].mean()
    rc = [stats.spearmanr(s.length, s[m]).statistic for _, s in c.groupby(["run", "gene"]) if s.length.nunique() > 5]
    return {"C4ii_seed_rho_fixed5": round(float(np.mean(rr)), 3) if rr else None,
            "C5_rho_len_gene": round(float(stats.spearmanr(L, v).statistic), 3),
            "C5_rho_len_cell_median": round(float(np.nanmedian(rc)), 3)}


def lodo(c: pd.DataFrame, m: str, full_g: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    rows = []
    fullc = cm.criteria(full_g, "dl_med")
    for dnr in sorted(c.donor_id.unique()):
        g = scores(c[c.donor_id != dnr], m)
        rhos = []
        for run in sorted(g.run.unique()):
            a = full_g[full_g.run == run].set_index("gene").dl_med
            b = g[g.run == run].set_index("gene").dl_med
            idx = a.index.intersection(b.index)
            rhos.append(stats.spearmanr(a[idx], b[idx]).statistic)
        cr = cm.criteria(g, "dl_med")
        c1s = {k: sum(x >= 0.9 for x in v) for k, v in cr["C1_pct"].items()}
        rows.append({"donor_left_out": dnr, "rho_min": min(rhos), "rho_mean": float(np.mean(rhos)),
                     "C1": cr["C1"], "C2": cr["C2"], "C3": cr["C3"], "C3_panel": cr["C3_panel"],
                     "VCAM1_C1_seeds": c1s["VCAM1"], "HAVCR1_C1_seeds": c1s["HAVCR1"], "C2_wins": cr["C2_wins"]})
    L = pd.DataFrame(rows)
    fc1 = {k: sum(x >= 0.9 for x in v) >= 2 for k, v in fullc["C1_pct"].items()}
    summ = {"rho_min": round(float(L.rho_min.min()), 3), "rho_median": round(float(L.rho_mean.median()), 3),
            "n_donors": int(len(L)),
            "n_flip": {k: int((L[k] != fullc[k]).sum()) for k in ["C1", "C2", "C3", "C3_panel"]},
            "n_flip_VCAM1_C1": int(((L.VCAM1_C1_seeds >= 2) != fc1["VCAM1"]).sum()),
            "n_flip_HAVCR1_C1": int(((L.HAVCR1_C1_seeds >= 2) != fc1["HAVCR1"]).sum())}
    return L, summ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--logits", action="store_true")
    ap.add_argument("--check-gsm", default=None, help="EXP-003 の gene_seed_metrics.csv（a2 の再実装の照合）")
    ap.add_argument("--extra-results", default=None, help="同じ細胞・重心の追加 run（DE マッチ）を統合する")
    ap.add_argument("--donors", default=None, help="donor を限定（カンマ区切りファイル or JSON リスト）")
    ap.add_argument("--exclude-role", default=None,
                    help="この role の遺伝子を除く（例: de_matched → EXP-003 と同じ 34 遺伝子で split-half を比べる。"
                         "a2 は null 細胞だけで補正するので他の遺伝子の値は変わらない）")
    a = ap.parse_args()
    res, out = cm.ROOT / a.results, cm.ROOT / a.out
    out.mkdir(parents=True, exist_ok=True)
    c = cm.load_cells(res, logits=a.logits)
    if a.donors:
        keep = set(json.loads(Path(a.donors).read_text()))
        c = c[c.donor_id.isin(keep)].reset_index(drop=True)
    if a.exclude_role:
        c = c[c.role != a.exclude_role].reset_index(drop=True)
    c = cm.add_metrics(c)
    summary = {"results": a.results, "n_cells_rows": int(len(c)), "n_donors": int(c.donor_id.nunique()),
               "min_cells_donor": MIN_CELLS_DONOR}
    if a.check_gsm:
        gsm = pd.read_csv(cm.ROOT / a.check_gsm)
        mine = c.groupby(["run", "gene"]).m_a2.mean()
        ref = gsm.set_index(["run", "gene"]).m_a2
        idx = mine.index.intersection(ref.index)
        mx = float((mine[idx] - ref[idx]).abs().max())
        summary["a2_reimpl_max_abs_diff_vs_gsm"] = mx
        assert len(idx) == len(ref) and mx <= 1e-12, f"a2 の再実装が一致しない: {mx}"
    metrics = ["m_a2", "m_raw"] + (["dP_goal_a2"] if a.logits else [])
    allg = []
    for m in metrics:
        g = scores(c, m)
        g["metric"] = m
        allg.append(g)
        summary[m] = {k: cm.criteria(g, k) for k in ["cell", "dl_med", "dl_mean"]}
        summary[m]["split_half"] = split_half(c, m)
        summary[m]["C4ii_C5"] = c4c5(c, g, m)
        summary[m]["sign_excess"] = sign_and_excess(c, m)
        if m == "m_a2":
            L, ls = lodo(c, m, g)
            L.to_csv(out / "lodo.csv", index=False)
            summary[m]["lodo"] = ls
            donor_means(c, m).assign(symbol=lambda t: t.gene.map(cm.SYM)).to_csv(out / "donor_means_a2.csv", index=False)
            # donor ごとの細胞数（遺伝子×seed）
            nd = donor_means(c, m).groupby(["run", "gene"]).size()
            summary["n_donors_per_gene_seed"] = {"min": int(nd.min()), "median": float(nd.median()), "max": int(nd.max())}
    if (c.role == "de_matched").any():  # 拡張 run: (a) と同じ A-P / A-S2 を同じ定義で再計算
        de = {}
        for m in metrics:
            g = scores(c, m)
            per = {}
            for run, s in g.groupby("run"):
                v = s.set_index("gene")
                for col in ["cell", "dl_med"]:
                    null = s[s.role == "negative"][col].to_numpy()
                    dev = v.loc[[k for k in cm.DE if k in v.index], col]
                    pv = v.loc["ENSG00000162692", col]
                    ph = v.loc["ENSG00000113249", col]
                    per.setdefault(col, {})[run] = {
                        "VCAM1_pct_in_DE": float((dev < pv).mean()), "HAVCR1_pct_in_DE": float((dev < ph).mean()),
                        "VCAM1_null_pct": cm.null_pct(pv, null), "HAVCR1_null_pct": cm.null_pct(ph, null),
                        "DE_null_pct": {cm.DE[k]: cm.null_pct(x, null) for k, x in dev.items()}}
            for col, d in per.items():
                runs = sorted(d)
                n_de = sum(sum(d[r]["DE_null_pct"][k] >= 0.90 for r in runs) >= 2 for k in cm.DE.values())
                de[f"{m}:{col}"] = {"VCAM1_pct_in_DE": [d[r]["VCAM1_pct_in_DE"] for r in runs],
                                    "HAVCR1_pct_in_DE": [d[r]["HAVCR1_pct_in_DE"] for r in runs],
                                    "n_DE_pass_C1": int(n_de),
                                    "DE_null_pct_mean": {k: round(float(np.mean([d[r]["DE_null_pct"][k] for r in runs])), 3)
                                                         for k in cm.DE.values()},
                                    "per_seed": d}
        summary["de_matched"] = de
        print("DE:", json.dumps({k: {x: v[x] for x in ("VCAM1_pct_in_DE", "HAVCR1_pct_in_DE", "n_DE_pass_C1")}
                                 for k, v in de.items()}, ensure_ascii=False))
    pd.concat(allg).to_csv(out / "gene_seed_scores.csv", index=False)
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=1, default=float) + "\n")
    for m in metrics:
        print(f"== {m}")
        for k in ["cell", "dl_med", "dl_mean"]:
            r = summary[m][k]
            print(f"  {k:8s} C1={r['C1']} {r['C1_pct']} C2={r['C2_wins']} C3={r['C3']} {r['C3_pct']} panel={r['C3_panel_n_ok']} nullneg={r['null_frac_negative']}")
        print("  split-half:", summary[m]["split_half"])
        print("  sign/excess:", {k: v for k, v in summary[m]["sign_excess"].items() if k in ("VCAM1", "HAVCR1", "null_frac_pos_median")})
        if "lodo" in summary[m]:
            print("  LODO:", summary[m]["lodo"])
    print({k: summary[k] for k in summary if k.startswith(("a2_", "n_"))})


if __name__ == "__main__":
    main()

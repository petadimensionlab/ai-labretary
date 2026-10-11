#!/usr/bin/env python3
"""EXP-20261003-001 — 新しい 6 個の z を、対照 1,000 個と同じ尺度で出し、V-1 / V-2 を当てる。

## 何をするか

既存の 2 つの結果（EXP-20261001-001 の 717 遺伝子、EXP-20261002-002 の対照 500 個）と、
この実験の 6 個を**細胞ごとの表のところから合算し**、対照 1,000 個で標準化する。
行を足すのでは足りない理由は EXP-20261002-002 plan.md §6.0 と同じ（対照プールは z の
標準化だけでなく、遺伝子ごとの残差 D_raw そのものに入っている）。

**ただし今回は対照プールを 1 個も変えない**ので、既存の 1,217 遺伝子の z は**1 つも動かないはず**
である（`gene_stats()` の回帰は対照だけで当てるため、対照でない遺伝子を足しても影響しない）。
これを **再現の検査**として機械的に確かめる。ずれたら**止める**。

## 手順

 0. 状態の重心が 3 つの結果で**同一**であることを確かめる（違えば**止める**）
 1. 細胞ごとの表を作る（摂動は回さない。保存済みの埋め込みを読むだけ）
 2. **再現の検査**: 既存 1,217 遺伝子の z と、セットの平均（P-1 / P-4 / P-5）が
    `criteria_007.json` を再現するか。しなければ**止める**
 3. 新しい 6 個の z を seed ごとに出し、**V-1**（3 seed すべてで z > 2.73）を当てる
 4. **V-2**（偶然の水準 0.9%、6 個での期待件数 0.054 件）と突き合わせる
 5. 事前登録した「解釈の規則」のどれに当たるかを出す

走らせ方:
    .venv/bin/python experiments/EXP-20261003-001/analysis/evaluate_008.py
  試し実行の出力で通すとき（**判定ではない**）:
    .venv/bin/python experiments/EXP-20261003-001/analysis/evaluate_008.py --smoke
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261003-001"
P01 = ROOT / "experiments/EXP-20261001-001"
P02 = ROOT / "experiments/EXP-20261002-002"
for p in (ROOT / "scripts/analysis", P01 / "analysis", P02 / "analysis"):
    sys.path.insert(0, str(p))
from movement_correction import gene_stats                       # noqa: E402
from evaluate_005 import eligibility                             # noqa: E402
from evaluate_007 import load_roots, check_centroids             # noqa: E402

BLOCKED = {"ENSG00000141510": "TP53", "ENSG00000118785": "SPP1", "ENSG00000148346": "LCN2"}
Z_PER_GENE_TOL = 0.002   # 既存 1,217 遺伝子の z の再現の許容差（事前登録。criteria と同じ値）


def main(smoke: bool = False) -> None:
    cfg = yaml.safe_load((E / "config.yaml").read_text())
    mg, dec = cfg["merge"], cfg["decision"]
    src_cfg = yaml.safe_load((P01 / "config.yaml").read_text())
    a = src_cfg["analysis"]
    pm = a["pos_matched_permutation"]
    mc = cfg["metric"]["confound_correction"]
    z_floor = float(dec["V-1"]["閾値"])

    # ---- 事前登録どおりの設定で動いているかの検査（違えば止める）
    assert mc["variables"] == ["lp", "ll"], f"補正の版が事前登録と違う（{mc['variables']}）"
    assert mc["drop_zero_movement"] is True, "移動量 0 の細胞を外す設定になっていない"
    assert cfg["metric"]["primary_form"].startswith("projection_"), "主指標が射影になっていない"
    assert z_floor == 2.73, "検出下限を動かしてはいけない"

    # ---- 遺伝子の名前と、対照 1,000 個
    gj = json.loads((ROOT / mg["genes"]).read_text())["genes"]
    ens = {g["symbol"]: g["ensembl"] for g in gj}
    subs = json.loads((ROOT / mg["subsets"]).read_text())["subsets"]
    sets_def = {k: [ens[x["symbol"]] for x in v] for k, v in subs.items()}
    null500 = [ens[x["symbol"]] for x in
               json.loads((ROOT / mg["null_pool"]).read_text())["versions"]["B_500"]]
    xtra = json.loads((ROOT / mg["extra_null"]).read_text())["追加する対照"]
    extra_null = [x["ensembl"] for x in xtra]
    ens |= {x["symbol"]: x["ensembl"] for x in xtra}
    nulls_all = null500 + extra_null
    assert len(set(nulls_all)) == 1000, f"対照が 1,000 個でない（{len(set(nulls_all))}）"
    targets = dict(dec["対象"])                       # 記号 → Ensembl
    t_sym = {v: k for k, v in targets.items()}
    ens |= targets
    sym = {v: k for k, v in ens.items()}

    new_res = E / ("results-smoke" if smoke else cfg["output"]["results_dir"])
    roots = [("既存(EXP-20261001-001)", ROOT / mg["existing_roots"][0]["results_dir"]),
             ("追加の対照(EXP-20261002-002)", ROOT / mg["existing_roots"][1]["results_dir"]),
             ("新しい 6 個(EXP-20261003-001)", new_res)]
    seeds = [cfg["robustness"]["seeds"][0]] if smoke else list(cfg["robustness"]["seeds"])

    out: dict = {
        "実験": "EXP-20261003-001",
        "問い": "目標の軸（参照状態の線維芽細胞 → 活性化線維芽細胞 + 筋線維芽細胞）は、"
                "本当に「筋線維芽細胞らしさ」を測っているか",
        "mode": "smoke" if smoke else "final",
        "seed": seeds,
        "判定の単位": dec["単位"],
    }
    if smoke:
        out["注（試し実行）"] = (
            "**試し実行の出力に対する評価であり、判定ではない。** 新しい 6 個は 20 細胞で、"
            "状態の重心も本番と違う条件（40 細胞）で作っている。目的は「評価コードが最後まで"
            "通ること」の確認だけである。")

    # ---------------------------------------------------------- 0. 状態の重心
    cen = check_centroids(roots, seeds, mg["centroid_max_abs_diff"])
    out["0. 状態の重心が既存と同一か"] = cen
    bad_cen = [k for k, v in cen.items() if not v.get("一致", False)]
    out["0. 合否"] = {"不一致": bad_cen, "pass": not bad_cen}
    if bad_cen and not smoke:
        (E / "results").mkdir(parents=True, exist_ok=True)
        (E / "results/criteria_008.json").write_text(
            json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
        raise SystemExit(f"**止める**: 状態の重心が既存と違う {bad_cen}。別の尺度で測っている")

    # ---------------------------------------------------------- 1. 細胞ごとの表
    c, drop = load_roots(cfg, roots, ROOT / cfg["data"]["tokenized_dataset"])
    c = c[c.seed.isin(seeds)].reset_index(drop=True)
    uni = cfg["metric"]["uniform_component_states"]
    c["m_goal"] = c.m_proj
    c["m_uniform"] = c[[f"s_{x}" for x in uni]].mean(axis=1)
    present = set(c.gene.unique())
    missing_t = [t_sym[g] for g in targets.values() if g not in present]
    out["1. 細胞ごとの表"] = {
        "遺伝子数": len(present), "移動量 0 で外した細胞": drop,
        "出どころごとの遺伝子数":
            {k: int(c[c.root_label == k].gene.nunique()) for k in c.root_label.unique()},
        "新しい 6 個が揃っているか": {"欠け": missing_t, "pass": not missing_t},
        "新しい 6 個の細胞数":
            {t_sym[g]: {str(s): int((c[(c.seed == s) & (c.gene == g)]).shape[0]) for s in seeds}
             for g in targets.values() if g in present}}
    if missing_t and not smoke:
        raise SystemExit(f"**止める**: 新しい遺伝子が結果に無い {missing_t}")

    # ---------------------------------------------------------- 2. 表を作る（対照 1,000 個）
    min_elig = 1 if smoke else pm["min_eligible_null_genes"]
    nl_set = set(nulls_all) & present
    assert len(nl_set) == 1000, f"対照が結果に 1,000 個そろっていない（{len(nl_set)}）"
    tabs = {}
    for s in seeds:
        cs = c[c.seed == s]
        nl = {h: hd for h, hd in cs[cs.gene.isin(nl_set)].groupby("gene")}
        el = eligibility(cs, nl, pm["tolerance_log2"], pm["min_cover_frac"])
        f = gene_stats(cs, nl, el, "m_goal", min_elig, use_dnorm=False)
        u = gene_stats(cs, nl, el, "m_uniform", min_elig, use_dnorm=False)
        n = f[f.gene.isin(nl_set)]
        tabs[s] = pd.DataFrame({
            "gene": f.gene.to_numpy(),
            "z": ((f.D_raw - n.D_raw.mean()) / n.D_raw.std(ddof=1)).to_numpy(),
            "q": f.q.to_numpy(), "n_cells": f.n_cells.to_numpy(),
            "dnorm_med": f.dnorm_med.to_numpy(), "lp_med": f.lp_med.to_numpy(),
            "ll_med": f.ll_med.to_numpy(), "relpos_med": f.relpos_med.to_numpy(),
            "uniform": u.set_index("gene").D_raw.reindex(f.gene).to_numpy()}).reset_index(drop=True)

    # ---------------------------------------------------------- 3. 再現の検査（**止まりうる**）
    rep: dict = {"意味": "対照プールを 1 個も変えていないので、既存 1,217 遺伝子の z は"
                         "1 つも動かないはず。動いたら合算の仕方が誤っている"}
    old = pd.concat([pd.read_csv(P02 / f"results/gene_seed_007-{s}.csv") for s in seeds])
    per_seed = {}
    for s in seeds:
        o = old[old.seed == s].set_index("gene").z
        n_ = tabs[s].set_index("gene").z.reindex(o.index)
        d = (n_ - o).abs()
        per_seed[str(s)] = {"照合した遺伝子": int(d.notna().sum()),
                            "z の最大の差": float(d.max()),
                            "z の中央の差": float(d.median())}
    rep["既存 1,217 遺伝子の z"] = per_seed
    rep["z の最大の差（全 seed）"] = max(v["z の最大の差"] for v in per_seed.values())
    rep["z の許容差"] = Z_PER_GENE_TOL

    ref = json.loads((ROOT / mg["criteria"]).read_text())["3. 対照を増やした評価（**これが判定に使う値**）"]
    pos_i = [list(cfg["robustness"]["seeds"]).index(s) for s in seeds]
    pick = lambda v: [v[i] for i in pos_i]  # noqa: E731
    setmean = lambda g: [round(float(np.mean([tabs[s].set_index("gene").z[x]  # noqa: E731
                                              for x in g if x in set(tabs[s].gene)])), 3)
                         for s in seeds]
    rep["セットの平均 z"] = {
        "主判定 P-1（介入 119）": {"今回": setmean(sets_def["INT"]),
                                   "前回": pick(ref["P-1 主判定"]["set_mean_z"])},
        "陽性対照 P-4（71）": {"今回": setmean(sets_def["positive_control71"]),
                               "前回": pick(ref["P-4 陽性対照"]["set_mean_z"])},
        "中立性 P-5（4）": {"今回": setmean(sets_def["neutral_adopted"]),
                            "前回": pick(ref["P-5 中立性"]["set_mean_z"])}}
    tol = mg["reproduction_tolerance"]
    bad = [k for k, v in rep["セットの平均 z"].items()
           if max(abs(np.array(v["今回"], float) - np.array(v["前回"], float))) > tol]
    rep |= {"セットの平均の許容差": tol, "一致しなかったセット": bad,
            "pass": bool(not bad and rep["z の最大の差（全 seed）"] <= Z_PER_GENE_TOL)}
    out["2. 再現の検査（既存 1,217 遺伝子が動いていないか）"] = rep
    if not rep["pass"] and not smoke:
        (E / "results").mkdir(parents=True, exist_ok=True)
        (E / "results/criteria_008.json").write_text(
            json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
        raise SystemExit(f"**止める**: 既存の値を再現できなかった（セット {bad} / "
                         f"z の最大の差 {rep['z の最大の差（全 seed）']}）")

    # ---------------------------------------------------------- 4. V-1（本判定）
    rows = []
    for s_name, g in targets.items():
        if g not in present:
            rows.append({"遺伝子": s_name, "ensembl": g, "状態": "結果が無い"})
            continue
        zs = [float(tabs[s].set_index("gene").z[g]) for s in seeds]
        rows.append({
            "遺伝子": s_name, "ensembl": g,
            "z": [round(x, 3) for x in zs],
            "z の最小": round(min(zs), 3), "z の平均": round(float(np.mean(zs)), 3),
            "V-1（3 seed すべてで z > 2.73）": bool(all(x > z_floor for x in zs)),
            "対照の中での順位 q": [round(float(tabs[s].set_index("gene").q[g]), 4) for s in seeds],
            "細胞数": [int(tabs[s].set_index("gene").n_cells[g]) for s in seeds],
            "元の順位の中央 lp_med": [round(float(tabs[s].set_index("gene").lp_med[g]), 3) for s in seeds],
            "移動量の中央 dnorm_med": [round(float(tabs[s].set_index("gene").dnorm_med[g]), 4)
                                       for s in seeds]})
    V1 = pd.DataFrame(rows)
    passed = [r["遺伝子"] for r in rows if r.get("V-1（3 seed すべてで z > 2.73）")]
    out["3. V-1 単一遺伝子の検出下限"] = {
        "条件": "3 seed すべてで z > 2.73（単一遺伝子の検出下限。EXP-20260928-002 の実測）",
        "遺伝子ごと": rows,
        "合格した遺伝子": passed, "合格数": len(passed)}

    # ---------------------------------------------------------- 5. V-2（偶然の水準）
    p = float(dec["V-2"]["偶然の水準 p"])
    k = len(passed)
    from scipy import stats as _st
    out["4. V-2 偶然の水準"] = {
        "偶然の水準 p（対照 1,000 個で実測）": p,
        "出どころ": dec["V-2"]["出どころ"],
        "6 個での期待件数": round(p * len(targets), 4),
        "実際の合格数 k": k,
        "P(k 件以上 | 偶然)": round(float(_st.binomtest(k, len(targets), p,
                                                        alternative="greater").pvalue), 5)
        if k > 0 else 1.0,
        "読む向き": dec["読む向き"]}

    # ---------------------------------------------------------- 6. 解釈の規則（事前登録）
    acta2 = next((r for r in rows if r["遺伝子"] == "ACTA2"), {})
    myh11 = next((r for r in rows if r["遺伝子"] == "MYH11"), {})
    a_ok = bool(acta2.get("V-1（3 seed すべてで z > 2.73）"))
    m_ok = bool(myh11.get("V-1（3 seed すべてで z > 2.73）"))
    if a_ok:
        rule = ("規則 1: **ACTA2 が V-1 を満たした → 目標の軸は筋線維芽細胞らしさを測っている。**")
    elif m_ok:
        rule = ("規則 2: **ACTA2 は満たさず MYH11 が満たした → 目標の軸が拾っているのは"
                "平滑筋らしさであって、筋線維芽細胞ではない。**")
    elif k == 0:
        rule = ("規則 3: **6 個とも V-1 を満たさなかった → 目標の軸は筋線維芽細胞らしさを"
                "測っていない。** これまでの結論の解釈を大きく変える。")
    else:
        rule = (f"規則 4（事前登録の 3 つに当てはまらない）: ACTA2 も MYH11 も満たさないが、"
                f"ほかの {k} 個（{passed}）が満たした。**事前登録の規則で読み切れないので、"
                f"判定は保留し、事実だけを報告する。**")
    out["5. 解釈の規則（結果を見る前に固定したもの）"] = {
        "ACTA2 の V-1": a_ok, "MYH11 の V-1": m_ok, "合格数": k, "当たる規則": rule}

    out["判定の提案"] = rule
    out["注"] = ("判定（関門 G3）は PI が行う。この道具は提案までにとどまる。"
                 "また、これはコンピュータ上の遺伝子摂動の結果であり、"
                 "実験的に確認された生物学的事実ではない。"
                 "**負の値は対照の負の裾が薄く校正されていないため解釈しない。**")

    (E / "results").mkdir(parents=True, exist_ok=True)
    name = "criteria_008-smoke.json" if smoke else "criteria_008.json"
    (E / "results" / name).write_text(
        json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    for s in seeds:
        tabs[s].assign(seed=s).to_csv(
            E / "results" / f"gene_{'smoke' if smoke else 'seed'}_008-{s}.csv", index=False)
    V1.to_csv(E / "results" / f"targets_{'smoke' if smoke else 'final'}_008.csv", index=False)

    print(json.dumps({k2: out[k2] for k2 in
                      ("mode", "seed", "0. 合否", "4. V-2 偶然の水準",
                       "5. 解釈の規則（結果を見る前に固定したもの）", "判定の提案")}
                     | {"2. 再現の検査": {k2: rep[k2] for k2 in
                                          ("z の最大の差（全 seed）", "セットの平均 z",
                                           "一致しなかったセット", "pass")}},
                     ensure_ascii=False, indent=1, default=float))
    print("\n--- V-1 遺伝子ごと ---")
    for r in rows:
        if "z" in r:
            print(f"{r['遺伝子']:8s} z = {r['z']}  最小 {r['z の最小']:+.3f}  "
                  f"細胞 {r['細胞数']}  V-1: {'合格' if r['V-1（3 seed すべてで z > 2.73）'] else '不合格'}")
        else:
            print(f"{r['遺伝子']:8s} {r['状態']}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="試し実行の出力で評価コードを通す（**判定ではない**）")
    main(**vars(ap.parse_args()))

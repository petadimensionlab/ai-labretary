#!/usr/bin/env python3
"""移動量の補正の検算（事後・事前登録の外）— EXP-20260929-001 の結果で確かめる。

**これは事前登録された判定ではない。** EXP-20260929-001 の判定（判定不能）はすでに確定しており、
本スクリプトはそれを変えない。目的は、**次の実験の主指標に移動量の補正を入れてよいか**を、
既存の摂動結果で確かめることである（PI 承認の段階 2、2026-09-30）。

## 確かめること

1. **従来の版が、確定した結果を再現するか**（再現しなければ、以降の比較は信用できない）
2. **移動量を入れた版で、中立であるべきセットが中立に戻るか**
3. **移動量を入れた版で、他のセットの判定が動かないか**（動くなら副作用が大きすぎる）

## 実行

    .venv/bin/python experiments/EXP-20260929-001/posthoc/movement_correction_check.py

摂動は 1 回も回さない。既存の埋め込みを読み直すだけ。出力は同じディレクトリの
`movement_correction_check.json`。細胞ごとの中間データは `cells_cache.parquet` に貯める
（2 回目以降はこれを使う。`--refresh` で作り直す）。
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
E = HERE.parent
ROOT = E.parents[1]
sys.path.insert(0, str(ROOT / "scripts/analysis"))

import movement_correction as mc  # noqa: E402


def load_evaluate_module():
    """事前登録済みの evaluate_003.py を**読むだけ**で取り込む（一切変更しない）。"""
    p = E / "analysis/evaluate_003.py"
    spec = importlib.util.spec_from_file_location("evaluate_003", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main() -> None:
    refresh = "--refresh" in sys.argv
    ev = load_evaluate_module()
    cfg = yaml.safe_load(open(E / "config.yaml"))
    a = cfg["analysis"]
    pm = a["pos_matched_permutation"]
    cr = a["criteria"]

    gj = json.loads((E / "work/genes_003.json").read_text())["genes"]
    ens = {g["symbol"]: g["ensembl"] for g in gj}
    sym = {v: k for k, v in ens.items()}
    subs = json.loads((E / "work/subsets_003.json").read_text())["subsets"]
    sets_def = {k: [ens[x["symbol"]] for x in v] for k, v in subs.items()}
    npool = json.loads((E / "work/null_pool_003.json").read_text())["versions"]
    null_versions = {k: [ens[x["symbol"]] for x in v] for k, v in npool.items()}
    NULL = a["primary_null_version"]

    # --- 細胞ごとのデータ（キャッシュあり）
    cache = HERE / "cells_cache.parquet"
    if cache.exists() and not refresh:
        print(f"細胞ごとのデータを読み込み: {cache}")
        c = pd.read_parquet(cache)
    else:
        print("埋め込みから細胞ごとのデータを作る（数分かかる）…")
        t0 = time.time()
        c = ev.load(cfg)
        print(f"  完了 {time.time() - t0:.1f} 秒, {len(c):,} 行")
        c.to_parquet(cache, index=False)
        print(f"  保存: {cache}")
    ev.metric_columns(c)                       # m_goal ほかを付ける
    seeds = sorted(c.seed.unique().tolist())
    assert seeds == sorted(cfg["robustness"]["seeds"]), f"seed が計画と違う: {seeds}"
    print(f"seed: {seeds}  遺伝子: {c.gene.nunique()}  細胞ごとの行: {len(c):,}")

    nulls_all = set(null_versions["B_500"])
    REPORT = ["housekeeping3", "negctrl_panel7", "secretory_all", "org_markers7", "myof_markers24"]
    LABEL = {"housekeeping3": "恒常発現遺伝子 3 個（中立であるべき）",
             "negctrl_panel7": "中立対照パネル 7 個（中立であるべき）",
             "secretory_all": "分泌プログラム 24 個（主判定）",
             "org_markers7": "目標マーカー 7 個（目標の軸の陽性対照）",
             "myof_markers24": "筋線維芽細胞マーカー 24 個（第 2 の軸の陽性対照）"}

    out: dict = {"目的": "移動量の補正を主指標に入れてよいかを、既存の摂動結果で確かめる",
                 "事前登録": "外（事後の感度分析）。EXP-20260929-001 の確定した判定は変えない",
                 "seeds": seeds, "null_version": NULL, "axis": "goal（目標に特異的な成分 / 移動量）"}

    # --- 摂動が起きなかった細胞を外す
    # 移動量がちょうど 0 の細胞は、埋め込みが摂動前と 1 ビットも変わっていない
    # ＝**摂動がそもそも起きていない**。原因は「元の順位がすでに最上位の遺伝子を
    # 過剰発現させても入力が変わらない」ことで、実測でも該当細胞は全て pos == 1 だった。
    # 対数が取れないだけでなく、指標に 0 を混ぜると効果を 0 のほうへ薄める。
    # **2 つの版で同じ細胞を外す**（外し方の違いが版の違いに化けないようにする）。
    zero = c.dnorm <= 0
    zero_info = {"n": int(zero.sum()), "全体に占める割合": round(float(zero.mean()), 6),
                 "元の順位 pos の内訳": {str(k): int(v) for k, v in
                                        c.loc[zero, "pos"].value_counts().items()},
                 "指標の最大絶対値（すべて 0 のはず）": float(
                     c.loc[zero, ["s_goal", "s_start", "s_alt", "s_myof"]].abs().to_numpy().max())
                 if zero.any() else 0.0}
    out["0_摂動が起きなかった細胞"] = zero_info
    print(f"摂動が起きなかった細胞: {zero_info['n']} / {len(c):,} "
          f"（{100 * zero.mean():.4f}%、元の順位はすべて {list(zero_info['元の順位 pos の内訳'])}）", flush=True)
    c_use = c[~zero]

    # 比べられる対照（eligibility）は元の位置と入力の長さだけで決まり、**版によらない**。
    # seed ごとに 1 回だけ計算して両方の版で使い回す（版ごとに計算すると倍かかる）。
    tabs: dict[tuple[bool, int], pd.DataFrame] = {}
    tabs_orig: dict[int, pd.DataFrame] = {}      # 再現の確認用（全細胞・従来の版）
    for s in seeds:
        t0 = time.time()
        el_cache = HERE / f"eligibility-seed{s}.json"
        cs_all = c[c.seed == s]
        nl_all = {h: hd for h, hd in cs_all[cs_all.gene.isin(nulls_all)].groupby("gene")}
        if el_cache.exists() and not refresh:
            el = json.loads(el_cache.read_text())
            print(f"  seed {s}: 比較対象の一覧をキャッシュから読み込み", flush=True)
        else:
            el = ev.eligibility(cs_all, nl_all, pm["tolerance_log2"], pm["min_cover_frac"])
            el_cache.write_text(json.dumps(el))
            print(f"  seed {s}: 比較対象の一覧 {time.time() - t0:.1f} 秒", flush=True)

        # (A) 再現の確認: 確定した判定とまったく同じ道筋（全細胞・従来の版）
        f = ev.gene_stats(cs_all, nl_all, el, "m_goal", pm["min_eligible_null_genes"])
        nul = f[f.gene.isin(set(null_versions[NULL]))]
        f["z"] = (f.D_raw - nul.D_raw.mean()) / nul.D_raw.std(ddof=1)
        tabs_orig[s] = f.set_index("gene")

        # (B)(C) 比較: 摂動が起きなかった細胞を両版から同じように外す
        cs = c_use[c_use.seed == s]
        nl = {h: hd for h, hd in cs[cs.gene.isin(nulls_all)].groupby("gene")}
        for use in (False, True):
            t1 = time.time()
            f = mc.gene_stats(cs, nl, el, "m_goal", pm["min_eligible_null_genes"], use_dnorm=use)
            nul = f[f.gene.isin(set(null_versions[NULL]))]
            f["z"] = (f.D_raw - nul.D_raw.mean()) / nul.D_raw.std(ddof=1)
            tabs[(use, s)] = f.set_index("gene")
            print(f"    {'移動量あり' if use else '従来    '}: {time.time() - t1:.1f} 秒", flush=True)

    def _mean(t: pd.DataFrame, name: str) -> float:
        g = [x for x in sets_def[name] if x in t.index]
        return float(t.loc[g, "z"].mean())

    def set_mean_z(use: bool, s: int, name: str) -> float:
        return _mean(tabs[(use, s)], name)

    # --- 1. 確定した判定と同じ道筋で、同じ値が出るか（出なければ以降の比較は信用できない）
    published = {"secretory_all": {61: 0.024, 62: -0.036, 63: 0.053}}
    rep = {}
    for name, exp in published.items():
        for s, v in exp.items():
            got = _mean(tabs_orig[s], name)
            rep[f"{name}@seed{s}"] = {"確定値": v, "再計算": round(got, 4),
                                      "差": round(got - v, 4), "一致": abs(got - v) < 0.01}
    out["1_再現の確認（全細胞・従来の版）"] = rep
    out["1_判定"] = all(v["一致"] for v in rep.values())

    # --- 1b. 摂動が起きなかった細胞を外しても結果が動かないこと（外し方が結論を作らない確認）
    ex = {}
    for name in ("secretory_all", "housekeeping3", "org_markers7"):
        ex[name] = {f"seed{s}": {"全細胞": round(_mean(tabs_orig[s], name), 4),
                                 "0を外す": round(set_mean_z(False, s, name), 4)}
                    for s in seeds}
    out["1b_摂動が起きなかった細胞を外した影響（従来の版）"] = ex

    # --- 2/3. 両版のセット平均 z
    table = {}
    for name in REPORT:
        row = {"説明": LABEL[name], "n": len(sets_def[name])}
        for s in seeds:
            row[f"seed{s}"] = {"従来": round(set_mean_z(False, s, name), 3),
                               "移動量あり": round(set_mean_z(True, s, name), 3)}
        row["3seed平均"] = {"従来": round(float(np.mean([set_mean_z(False, s, name) for s in seeds])), 3),
                            "移動量あり": round(float(np.mean([set_mean_z(True, s, name) for s in seeds])), 3)}
        table[name] = row
    out["2_セット平均z"] = table

    # --- 中立性（対照の中での分位 q が 0.05〜0.95 に入るか）
    lo, hi = cr["neutral_range"]
    neu = {}
    for name in ("housekeeping3", "negctrl_panel7"):
        d = {}
        for use in (False, True):
            n_in = 0
            per = {}
            for s in seeds:
                t = tabs[(use, s)]
                g = [x for x in sets_def[name] if x in t.index]
                qs = t.loc[g, "q"]
                ok = bool(qs.notna().all() and ((qs >= lo) & (qs <= hi)).all())
                per[f"seed{s}"] = {"q": [None if pd.isna(x) else round(float(x), 3) for x in qs],
                                   "範囲内": ok}
                n_in += int(ok)
            d["移動量あり" if use else "従来"] = {"seedごと": per, "範囲内のseed数": n_in,
                                                "合格": n_in >= cr["neutral_min_seeds"]}
        neu[name] = d
    out["3_中立性（範囲%.2f〜%.2f, 必要%dseed）" % (lo, hi, cr["neutral_min_seeds"])] = neu

    # --- 残差と移動量の順位相関（交絡が減ったか）
    from scipy import stats
    corr = {}
    for use in (False, True):
        corr["移動量あり" if use else "従来"] = {
            f"seed{s}": round(float(stats.spearmanr(
                tabs[(use, s)].loc[[g for g in null_versions[NULL] if g in tabs[(use, s)].index],
                                   "D_raw"],
                tabs[(use, s)].loc[[g for g in null_versions[NULL] if g in tabs[(use, s)].index],
                                   "dnorm_med"]).statistic), 3) for s in seeds}
    out["4_残差と移動量の順位相関（対照%d個）" % len(null_versions[NULL])] = corr

    # --- 恒常発現遺伝子の移動量の分位（原因の確認）
    dn = {}
    t = tabs[(False, seeds[0])]
    nullmed = t.loc[[g for g in null_versions[NULL] if g in t.index], "dnorm_med"]
    for g in sets_def["housekeeping3"]:
        if g in t.index:
            v = float(t.loc[g, "dnorm_med"])
            dn[sym[g]] = {"移動量の中央値": round(v, 4),
                          "対照の中での分位": round(float((nullmed < v).mean()), 3)}
    dn["（対照の中央値）"] = round(float(nullmed.median()), 4)
    out["5_恒常発現遺伝子の移動量（seed%d）" % seeds[0]] = dn

    out["注記"] = (
        "本スクリプトは事前登録の外にある事後の感度分析であり、EXP-20260929-001 の判定は変えない。"
        "どちらの版を主判定に使うかは、次の実験の事前登録で結果を見る前に固定すること。"
        "この補正は『中立であるべきものが中立に戻る』ことを目的としており、"
        "目標の軸の陽性対照が立たない問題は解決しない。")

    p = HERE / "movement_correction_check.json"
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    print(f"\n書き出し: {p}")

    # --- 画面表示
    print("\n=== 1. 確定した判定と同じ道筋で同じ値が出るか ===")
    for k, v in rep.items():
        print(f"  {k}: 確定 {v['確定値']:+.3f} / 再計算 {v['再計算']:+.3f} "
              f"（差 {v['差']:+.4f}）{'一致' if v['一致'] else '★不一致'}")
    print("\n=== 2. セット平均 z（目標の軸）===")
    print(f"  {'セット':<38} {'従来':>8} {'移動量あり':>11}")
    for name in REPORT:
        r = table[name]
        print(f"  {r['説明']:<36} {r['3seed平均']['従来']:>8.3f} {r['3seed平均']['移動量あり']:>11.3f}")
    print("\n=== 3. 残差と移動量の順位相関（小さいほど交絡が少ない）===")
    for k, v in corr.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()

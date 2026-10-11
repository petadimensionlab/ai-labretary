#!/usr/bin/env python3
"""EXP-20261004-001 — 共通細胞集合で測り直した 626 遺伝子に、事前登録した判定（plan.md §6）を当てる。

## この評価が旧（evaluate_005/007/008.py）と違う点

**交絡の補正を発現の有無で層別する**（plan.md §5.4、config.yaml の
`metric.confound_correction.form: stratified_by_expression`）。

common モード（全遺伝子を同じ 200 細胞で過剰発現）では、ある遺伝子がもともと入力に
載っていない細胞が大半を占める（candidate の多くは 200 細胞中 1〜3 細胞にしか載っていない）。
`lp = log2(1 + 元の順位)` はその細胞では定義できない（順位が無い）ので、

  - **載っている細胞**: `lp`・`ll` の完全二次式（6 項）＝ 旧の補正とまったく同じ式
  - **載っていない細胞**: `ll` だけの二次式（3 項）。`lp` は使わない

の 2 本を**対照遺伝子の細胞だけ**で別々に当て（自分自身は除く leave-one-out。
`scripts/analysis/movement_correction.gene_stats` と同じ正規方程式の仕立て）、
各遺伝子の D_raw は**両方の層の残差をプールした単純平均**にする
（movement_correction.gene_stats の平均ベースの定義をそのまま踏襲。層を分けても
「細胞ごとの残差の単純平均」という定義そのものは変えていない）。

**`eligibility()`・`q`（対照の中での順位）はこの実験では計算しない。** plan.md §6 の
判定規則（V-1〜V-4）はどれも z（対照 549 個の平均・標準偏差で標準化した値）だけで決まり、
対照どうしの位置・長さを突き合わせた q には依存しない。旧 evaluate_*.py にあった
`eligibility()`／`sham_sets()`（位置・長さを合わせた偽集合）も、ここでは単純な
「対照 549 個から復元ありで引く」偽集合（`set_null`。evaluate_005.py と同じ関数）に置き換えている
（plan.md §6.3 V-3a/b の定義どおり）。

## 対照・候補の遺伝子集合は config.yaml から機械的に再構成する

この評価コードを書く過程で、`work/genes_009.json` の「恒常発現・中立」フィールドが **6 個**
（housekeeping 3 と neutral_adopted の和集合）しかなく、config.yaml が実際に回す **11 個**
（negative_fixed 3 + extra_fixed の negctrl_panel 7 + housekeeping_extra 1）と一致しないことが分かった
（2026-10-04。genes_009.json の「測る遺伝子」が 621 個になり、config.yaml が実際に回す 626 個より
5 個少なかった）。**config.yaml 自体は自己無矛盾**であり（`run_isp.py:main()` と同じロジックで
再構成すると、対照 549（主対照 498・低発現 40・恒常発現と中立 11）+ 候補 8 + 陽性対照 71 −
重複 2（GSN・ELN）= 626 と一致する）、Labretary が `genes_009.json` を config.yaml から
作り直して訂正済み（同ファイルの「訂正（2026-10-04）」を参照。「測る遺伝子」は 626 個になった）。

→ **この評価コードは config.yaml から直接、本番と同じロジックで遺伝子集合を再構成する**
（`resolve_gene_sets()`）。genes_009.json が訂正された後の値とも一致することを確認済みだが、
二重の取り違えを防ぐため、`work/genes_009.json` は主対照 498・低発現 40 の内訳（D-2 用）の
参照にのみ使い、遺伝子集合そのものの定義には使わない。

## 手順（plan.md §6 のとおり）

 0. config の事前登録どおりの設定を検査（違えば止める）
 1. 細胞ごとの表を作る（摂動は回さない。保存済みの埋め込みを読むだけ）
 2. A-1〜A-8 の止める検査／記録する検査
 3. 遺伝子ごとに層別補正で z を出す（seed ごと）
 4. V-1 / V-1b / V-1c（候補 8 ＋ 記述として陽性対照・対照全部）
 5. V-2（この実験の対照 549 個で偶然の水準を実測）
 6. V-3a/b/c（陽性対照のものさし）・V-4（中立性）
 7. §6.6 負の側の目盛り・§6.7 旧の測り方との突き合わせ（記述）・§6.8 D-1/D-2（診断）
 8. §6.5 の解釈の規則（R-1〜R-7）のどれに当たるかを出す（**確定は PI**）
 9. §6.9 626 遺伝子すべての順位表

走らせ方:
    .venv/bin/python experiments/EXP-20261004-001/analysis/evaluate_009.py
  試し実行の出力で評価コードを通すとき（**判定ではない**。C-10）:
    .venv/bin/python experiments/EXP-20261004-001/analysis/evaluate_009.py --smoke
  層別補正の自己検査（データ不要）:
    .venv/bin/python experiments/EXP-20261004-001/analysis/evaluate_009.py --self-test
"""
from __future__ import annotations

import argparse
import json
import pickle
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261004-001"
P01 = ROOT / "experiments/EXP-20261001-001"
P02 = ROOT / "experiments/EXP-20261002-002"
sys.path.insert(0, str(ROOT / "scripts/analysis"))
sys.path.insert(0, str(P01 / "analysis"))
from movement_correction import design as design_full   # noqa: E402  6 項（lp, ll の完全二次）
from evaluate_005 import set_null                        # noqa: E402  復元ありで k 個引いて平均する

BLOCKED = {"ENSG00000141510": "TP53", "ENSG00000118785": "SPP1", "ENSG00000148346": "LCN2"}


# --------------------------------------------------------------------------- 設計行列
def design_reduced(ll) -> np.ndarray:
    """載っていない細胞用の補正（§5.4）: `ll` だけの二次式（3 項）。`lp` は定義しない。"""
    ll = np.asarray(ll, dtype=float)
    return np.column_stack([np.ones_like(ll), ll, ll ** 2])


def cos(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return (a @ b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b) + 1e-12)


# --------------------------------------------------------------------------- 遺伝子集合
def resolve_gene_sets(cfg: dict) -> dict:
    """config.yaml から、本番（run_isp.py:main）と同じロジックで遺伝子集合を再構成する。

    `work/genes_009.json` はいったん「恒常発現・中立」が 6 個（本来 11 個）のまま作られ、
    config.yaml の 626 個と一致しなかった期間があったが、Labretary が config.yaml から
    作り直して訂正済み（本ファイルの docstring 参照）。**ここでは config.yaml だけを真実の元にする**
    （genes_009.json を直接の定義には使わない）。
    """
    targets = list(cfg["perturbation"]["genes_to_perturb"])
    pos = list(cfg["controls"]["positive"])
    neg = list(cfg["controls"]["negative_fixed"])
    extra = [(x, grp["role"]) for grp in (cfg["controls"].get("extra_fixed") or []) for x in grp["genes"]]
    neg_all = set(neg) | {x for x, _ in extra}
    sampler = cfg["controls"]["negative_sampler"]
    n_rand = cfg["controls"]["negative_random_n"]
    if sampler.get("method") != "explicit":
        raise SystemExit(f"この実験の対照は explicit 方式のはず（{sampler.get('method')}）")
    rand = [x for x in sampler["genes"] if x not in set(targets) | set(pos) | neg_all][:n_rand or None]

    hk11 = sorted(neg_all)
    nulls_549 = sorted(set(rand) | neg_all)
    all_measured = sorted(set(targets) | set(pos) | set(rand) | neg_all)

    if len(set(targets)) != 8 or len(set(pos)) != 71 or len(hk11) != 11 or len(nulls_549) != 549:
        raise SystemExit(
            f"**止める**: 遺伝子集合の個数が事前登録と違う（候補 {len(set(targets))}/8、"
            f"陽性対照 {len(set(pos))}/71、恒常発現・中立 {len(hk11)}/11、対照 {len(nulls_549)}/549）")
    if len(all_measured) != 626:
        raise SystemExit(f"**止める**: 測る遺伝子の合計が 626 でない（{len(all_measured)}）")

    # D-2 用の内訳（498 主対照 / 40 低発現）。genes_009.json のこの 2 つのフィールドは
    # config.yaml の対照 538（= rand）と完全一致することを確かめて使う（本ファイル作成時に検算済み）。
    g009 = json.loads((E / "work/genes_009.json").read_text())
    main498, low40 = set(g009["主対照"]), set(g009["低発現の追加対照"])
    if (main498 | low40) != set(rand) or main498 & low40:
        raise SystemExit("**止める**: genes_009.json の主対照/低発現の内訳が config.yaml の対照 538 と一致しない")

    return {"候補": targets, "陽性対照71": pos, "恒常発現・中立11": hk11,
            "対照549": nulls_549, "主対照498": sorted(main498), "低発現40": sorted(low40),
            "測る遺伝子626": all_measured}


# --------------------------------------------------------------------------- 細胞ごとの表
def load_cells(cfg: dict, results_dir: str | None = None) -> tuple[pd.DataFrame, dict]:
    """細胞ごとの表を作る。**摂動は回さない**。保存済みの `cells.csv`・埋め込みを読むだけ。

    `lp`（元の順位の対数）は、その遺伝子がもともと入力に載っている細胞でしか定義できない。
    載っていない細胞は `expressed=False`・`lp=NaN` のまま残し、設計行列を作る側（§5.4）で
    層ごとに分ける。"""
    from datasets import load_from_disk
    from geneformer import TOKEN_DICTIONARY_FILE

    res = E / (results_dir or cfg["output"]["results_dir"])
    c = pd.read_csv(res / "cells.csv", dtype={"donor_id": str})
    n0 = len(c)
    c = c[~c.gene.isin(BLOCKED)].reset_index(drop=True)
    if n0 != len(c):
        raise SystemExit(f"**止める**(A-1): 盲検中の遺伝子が結果に {n0 - len(c)} 行ある")
    c["seed"] = c.run.astype(str).str.extract(r"run-(\d+)", expand=False).astype(int)

    tokd = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
    ds = load_from_disk(str(ROOT / cfg["data"]["tokenized_dataset"]))
    uids = sorted(c.cell_uid.unique())
    ids = {u: np.asarray(x) for u, x in zip(uids, ds.select(uids)["input_ids"])}
    pos = np.full(len(c), np.nan)
    expressed = np.zeros(len(c), dtype=bool)
    for gene, part in c.groupby("gene", sort=False):
        t = tokd.get(gene)
        for j, u in zip(part.index.to_numpy(), part.cell_uid):
            hit = np.flatnonzero(ids[u] == t) if t is not None else np.empty(0, dtype=int)
            if len(hit):
                pos[j] = hit[0]
                expressed[j] = True
    c["pos"] = pos
    c["expressed"] = expressed
    with np.errstate(divide="ignore", invalid="ignore"):
        c["lp"] = np.log2(1 + c["pos"])          # NaN のまま残す（載っていない細胞）
    c["ll"] = np.log2(c.length.astype(float))

    s_ = cfg["states"]
    names = {"goal": s_["goal_state"], "start": s_["start_state"]}
    names |= {a: a for a in s_.get("alt_states", [])}
    vals = {n: np.full(len(c), np.nan) for n in names}
    dnorm = np.full(len(c), np.nan)
    proj = np.full(len(c), np.nan)
    for (run, gene), part in c.groupby(["run", "gene"], sort=False):
        sz = np.load(res / run / "state_emb_vectors.npz")
        gd = res / run / "genes" / gene / "emb"
        o = np.load(gd / "orig.npy").astype(np.float64)
        p = np.load(gd / "pert.npy").astype(np.float64)
        if len(o) != len(part):
            raise SystemExit(f"{gd}: 埋め込み {len(o)} != cells.csv {len(part)}")
        idx = part.index.to_numpy()
        for n, cname in names.items():
            v = sz[cname].astype(np.float64)
            vals[n][idx] = cos(p, v) - cos(o, v)
        d = p - o
        dn = np.linalg.norm(d, axis=1)
        dnorm[idx] = dn
        vg, vs = sz[names["goal"]].astype(np.float64), sz[names["start"]].astype(np.float64)
        u = vg / np.linalg.norm(vg) - vs / np.linalg.norm(vs)
        u = u / np.linalg.norm(u)
        with np.errstate(divide="ignore", invalid="ignore"):
            proj[idx] = (d @ u) / dn
    for n in names:
        c[f"s_{n}"] = vals[n]
    c["dnorm"] = dnorm
    c["m_proj"] = proj

    zero = c.dnorm <= 0
    drop = {"n_cells_before": int(len(c)), "n_zero_movement": int(zero.sum()),
            "frac_zero_movement": round(float(zero.mean()), 6)}
    c = c[~zero].reset_index(drop=True)
    drop["n_cells_after"] = int(len(c))
    return c, drop


# --------------------------------------------------------------------------- 層別の遺伝子ごとの残差
def gene_stats_stratified(cs: pd.DataFrame, nulls: dict[str, pd.DataFrame], col: str) -> pd.DataFrame:
    """§5.4 の層別補正。対照遺伝子（自分自身を除く）の細胞だけで、発現の有無ごとに別の
    二次回帰を当て（載っている細胞: lp・ll の完全二次式 6 項／載っていない細胞: ll の二次式 3 項）、
    各遺伝子の D_raw は**両方の層の残差をプールした単純平均**にする
    （movement_correction.gene_stats と同じ「平均ベース」の定義。正規方程式による高速化も同じ仕立て）。

    全細胞が expressed のときは、この関数の出力は movement_correction.gene_stats(use_dnorm=False)
    と一致する（`--self-test` で確認）。
    """
    hs = list(nulls)
    Xf, yf, Xr, yr = {}, {}, {}, {}
    for h, hd in nulls.items():
        he, hu = hd[hd.expressed], hd[~hd.expressed]
        Xf[h] = design_full(he.lp.to_numpy(), he.ll.to_numpy()) if len(he) else np.zeros((0, 6))
        yf[h] = he[col].to_numpy()
        Xr[h] = design_reduced(hu.ll.to_numpy()) if len(hu) else np.zeros((0, 3))
        yr[h] = hu[col].to_numpy()
    Gf = {h: Xf[h].T @ Xf[h] for h in hs}
    Bf = {h: Xf[h].T @ yf[h] for h in hs}
    Gr = {h: Xr[h].T @ Xr[h] for h in hs}
    Br = {h: Xr[h].T @ yr[h] for h in hs}
    Gf_all, Bf_all = sum(Gf.values()), sum(Bf.values())
    Gr_all, Br_all = sum(Gr.values()), sum(Br.values())

    rows = []
    for gene, tg in cs.groupby("gene"):
        Gmf, Bmf = (Gf_all - Gf[gene], Bf_all - Bf[gene]) if gene in Gf else (Gf_all, Bf_all)
        Gmr, Bmr = (Gr_all - Gr[gene], Br_all - Br[gene]) if gene in Gr else (Gr_all, Br_all)
        beta_f = np.linalg.solve(Gmf, Bmf) if Gmf.any() else np.zeros(6)
        beta_r = np.linalg.solve(Gmr, Bmr) if Gmr.any() else np.zeros(3)
        te, tu = tg[tg.expressed], tg[~tg.expressed]
        resid = []
        if len(te):
            resid.append(te[col].to_numpy() - design_full(te.lp.to_numpy(), te.ll.to_numpy()) @ beta_f)
        if len(tu):
            resid.append(tu[col].to_numpy() - design_reduced(tu.ll.to_numpy()) @ beta_r)
        allr = np.concatenate(resid) if resid else np.array([np.nan])
        rows.append({
            "gene": gene, "col": col, "D_raw": float(np.mean(allr)),
            "n_cells": int(len(tg)), "n_expressed": int(len(te)), "n_unexpressed": int(len(tu)),
            "frac_expressed": float(len(te) / len(tg)) if len(tg) else np.nan,
            "lp_med": float(te.lp.median()) if len(te) else np.nan,
            "ll_med": float(tg.ll.median()),
            "dnorm_med": float(tg.dnorm.median()) if "dnorm" in tg else np.nan})
    return pd.DataFrame(rows)


def self_test() -> None:
    """全細胞が expressed のとき、movement_correction.gene_stats(use_dnorm=False) と一致するか。"""
    from movement_correction import gene_stats as gene_stats_old
    rng = np.random.default_rng(0)
    rows = []
    for gi in range(40):
        n = int(rng.integers(80, 160))
        lp = rng.normal(5, 1, n)
        ll = rng.normal(11, 0.3, n)
        m = 0.3 * lp - 0.2 * ll + rng.normal(0, 0.5, n)
        rows.append(pd.DataFrame({"gene": f"G{gi:03d}", "lp": lp, "ll": ll, "m": m, "expressed": True}))
    cs = pd.concat(rows, ignore_index=True)
    nulls = {g: d for g, d in cs.groupby("gene") if g < "G030"}
    elig = {g: [h for h in nulls if h != g] for g in cs.gene.unique()}
    old = gene_stats_old(cs, nulls, elig, "m", 1, use_dnorm=False).set_index("gene").sort_index()
    new = gene_stats_stratified(cs, nulls, "m").set_index("gene").sort_index()
    d = float(np.max(np.abs(old.D_raw.to_numpy() - new.loc[old.index].D_raw.to_numpy())))
    assert d < 1e-9, f"全細胞 expressed のとき旧実装と一致しない（最大差 {d:.3e}）"
    print(f"(a) 全細胞 expressed のとき movement_correction.gene_stats と一致: 最大差 {d:.2e}  OK")

    # 発現層を分けても、未発現細胞の残差が「その層の回帰」から出ていることを合成データで確認
    rows2 = []
    for gi in range(40):
        n_e = int(rng.integers(5, 40))
        n_u = int(rng.integers(100, 160))
        lp = rng.normal(5, 1, n_e)
        lle = rng.normal(11, 0.3, n_e)
        me = 0.3 * lp - 0.2 * lle + rng.normal(0, 0.5, n_e)
        llu = rng.normal(11, 0.3, n_u)
        mu = -0.1 * llu + rng.normal(0, 0.5, n_u) + 2.0   # 未発現層は水準が違う（挿入の効果を模す）
        rows2.append(pd.DataFrame({"gene": f"G{gi:03d}",
                                   "lp": list(lp) + [np.nan] * n_u, "ll": list(lle) + list(llu),
                                   "m": list(me) + list(mu),
                                   "expressed": [True] * n_e + [False] * n_u}))
    cs2 = pd.concat(rows2, ignore_index=True)
    nulls2 = {g: d for g, d in cs2.groupby("gene") if g < "G030"}
    out2 = gene_stats_stratified(cs2, nulls2, "m")
    assert out2.D_raw.abs().median() < 1.0, "層別補正が未発現層の水準シフトを吸収できていない"
    print("(b) 層別補正が未発現層の水準の違いを吸収する（D_raw の中央値の絶対値 "
          f"{out2.D_raw.abs().median():.3f} < 1.0）  OK")
    print("\n自己検査はすべて合格")


# --------------------------------------------------------------------------- 本体
def main(smoke: bool = False) -> None:
    cfg = yaml.safe_load((E / "config.yaml").read_text())
    dec, mc = cfg["decision"], cfg["metric"]["confound_correction"]

    # ---- 事前登録どおりの設定かの検査（違えば止める）
    assert mc["form"] == "stratified_by_expression", f"補正の版が事前登録と違う（{mc['form']}）"
    assert mc["variables"] == ["lp", "ll"], f"補正の説明変数が事前登録と違う（{mc['variables']}）"
    assert mc["drop_zero_movement"] is True, "移動量 0 の細胞を外す設定になっていない"
    assert cfg["metric"]["primary_form"].startswith("projection_"), "主指標が射影になっていない"
    assert cfg["cell_selection"]["mode"] == "common", "cell_selection.mode が common でない"
    z_floor = float(dec["V-1"]["閾値"])
    z_hold = float(dec["V-1b 判定保留の帯"]["閾値"])
    assert z_floor == 2.73, "検出下限を動かしてはいけない"
    assert z_hold == 3.6, "判定保留の帯の上端を動かしてはいけない"

    gs = resolve_gene_sets(cfg)
    g009 = json.loads((E / "work/genes_009.json").read_text())
    ens = {sym: ens_ for sym, ens_ in g009["候補"].items()}
    allg = (json.loads((P01 / "work/genes_005.json").read_text())["genes"]
            + json.loads((P02 / "work/extra_null_007.json").read_text())["追加する対照"])
    sym = {r["ensembl"]: r["symbol"] for r in allg} | {v: k for k, v in ens.items()}

    res = E / ("results-smoke" if smoke else cfg["output"]["results_dir"])
    seeds = [cfg["robustness"]["seeds"][0]] if smoke else list(cfg["robustness"]["seeds"])

    out: dict = {
        "実験": "EXP-20261004-001",
        "問い": "全遺伝子を同じ共通細胞集合で過剰発現したとき、探索で上位に出た候補 8 個は、"
                "較正された尺度でも最小検出可能効果量を超えるか",
        "mode": "smoke" if smoke else "final",
        "seed": seeds,
        "注（遺伝子集合）": ("この評価は config.yaml から直接 626 個を再構成している（docstring 参照）。"
                              "work/genes_009.json は 2026-10-04 に Labretary が config.yaml から"
                              "作り直して訂正済みだが、遺伝子集合そのものの定義にはここでも使っていない。")}
    if smoke:
        out["注（試し実行）"] = ("**試し実行の出力に対する評価であり、判定ではない。** "
                                 "遺伝子 40 個・共通細胞 30 個・seed 81 のみ。"
                                 "目的は『評価コードが最後まで通ること』の確認だけ。")

    # ---------------------------------------------------------- 1. 細胞ごとの表
    c, drop = load_cells(cfg, "results-smoke" if smoke else None)
    c = c[c.seed.isin(seeds)].reset_index(drop=True)
    c["m_goal"] = c.m_proj
    uni = cfg["metric"]["uniform_component_states"]
    c["m_uniform"] = c[[f"s_{x}" for x in uni]].mean(axis=1)
    present = set(c.gene.unique())
    out["1. 細胞ごとの表"] = {"遺伝子数": len(present), "移動量 0 で外した細胞": drop}

    # ---------------------------------------------------------- 2. A-1〜A-8
    checks: dict = {}
    # A-2: 共通細胞の同一性（check_common_cells_009.py を同じ設定で呼ぶ。本ファイルでは重複実装しない）
    cc_cfg = str(E / "work/config-smoke.yaml") if smoke else str(E / "config.yaml")
    cc_out = E / "work" / ("common_cells_check_009-smoke.json" if smoke else "common_cells_check_009.json")
    r = subprocess.run([sys.executable, str(E / "analysis/check_common_cells_009.py"),
                        "--results", str(res), "--out", str(cc_out), "--config", cc_cfg],
                       capture_output=True, text=True)
    checks["A-2 共通細胞の同一性"] = {"exit_code": r.returncode, "詳細": str(cc_out),
                                      "pass": r.returncode == 0}
    if r.returncode != 0 and not smoke:
        out["A チェック"] = checks
        (E / "results").mkdir(parents=True, exist_ok=True)
        (E / "results/criteria_009.json").write_text(
            json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
        raise SystemExit(f"**止める**(A-2): 共通細胞の同一性チェックが不合格。{r.stdout[-2000:]}")

    # A-3: state_rows の sha256 が事前登録と一致するか（機体間の再比較は C-11 で済ませている。
    # production では同じ seed を 2 台で重複計算していないため、ここでは「使ったファイルが
    # 事前登録したものと1ビットも違わないか」を機械的に確かめる = ファイルの入れ替わりを検出する）
    import hashlib
    pinned = cfg["merge"]["state_rows_sha256"]
    a3 = {}
    for s in seeds:
        p = ROOT / str(cfg["state_embs"]["rows_file"]).format(seed=s)
        got = hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None
        a3[str(s)] = {"path": str(p), "sha256": got, "pinned": pinned.get(str(s)),
                      "一致": got == pinned.get(str(s))}
    checks["A-3 state_rows の sha256 一致（machine agreement は C-11 で確認済み）"] = {
        "seed ごと": a3, "pass": all(v["一致"] for v in a3.values())}
    if not checks["A-3 state_rows の sha256 一致（machine agreement は C-11 で確認済み）"]["pass"] and not smoke:
        out["A チェック"] = checks
        (E / "results").mkdir(parents=True, exist_ok=True)
        (E / "results/criteria_009.json").write_text(
            json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
        raise SystemExit("**止める**(A-3): state_rows が事前登録の sha256 と一致しない")

    # A-5: 626/626（smoke では対象外）
    if not smoke:
        missing = sorted(set(gs["測る遺伝子626"]) - present)
        checks["A-5 全遺伝子が結果にあるか"] = {"欠け": missing, "pass": not missing}
        if missing:
            out["A チェック"] = checks
            (E / "results").mkdir(parents=True, exist_ok=True)
            (E / "results/criteria_009.json").write_text(
                json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
            raise SystemExit(f"**止める**(A-5): 結果に無い遺伝子がある: {missing[:10]}")

    # A-6（記録）
    checks["A-6 移動量 0 の割合"] = {"値": drop["frac_zero_movement"],
                                     "注意": "5% を超えたら記録", "超過": drop["frac_zero_movement"] > 0.05}
    # A-7（記録）: 3 seed の dFIB・dMFIB 重心を相互に比べる（機体ごとの一様成分の食い違いを見る）
    a7 = {}
    if len(seeds) > 1:
        cents = {}
        for s in seeds:
            z = np.load(res / f"run-{s}-L-1/state_emb_vectors.npz")
            cents[s] = {a: z[a].astype(float) for a in cfg["states"].get("alt_states", [])}
        for a in cfg["states"].get("alt_states", []):
            diffs = {f"{s1}-{s2}": float(np.max(np.abs(cents[s1][a] - cents[s2][a])))
                     for i, s1 in enumerate(seeds) for s2 in seeds[i + 1:]}
            a7[a] = diffs
    checks["A-7 dFIB・dMFIB 重心の seed 間の差（記録のみ。主指標には使っていない）"] = a7
    out["A チェック"] = checks

    # ---------------------------------------------------------- 3. 層別補正で z を出す
    nulls_549 = set(gs["対照549"]) & present
    frames = []
    for s in seeds:
        cs = c[c.seed == s]
        nl = {h: hd for h, hd in cs[cs.gene.isin(nulls_549)].groupby("gene")}
        f = gene_stats_stratified(cs, nl, "m_goal")
        n = f[f.gene.isin(nulls_549)]
        f = f.assign(seed=s,
                     z=(f.D_raw - n.D_raw.mean()) / n.D_raw.std(ddof=1),
                     symbol=f.gene.map(sym))
        frames.append(f)
    long = pd.concat(frames, ignore_index=True)
    gene_csv_path = (res / "gene_smoke_009.csv") if smoke else (E / "results/gene_seed_009.csv")
    gene_csv_path.parent.mkdir(parents=True, exist_ok=True)
    long.to_csv(gene_csv_path, index=False)

    def zmat(genes: list[str]) -> np.ndarray:
        gi = [g for g in genes if g in set(long.gene)]
        return np.column_stack([[float(long[(long.gene == g) & (long.seed == s)].z.iloc[0]) for g in gi]
                                for s in seeds]), gi

    # ---------------------------------------------------------- 4. V-1 / V-1b / V-1c（候補 8）
    rows = []
    for g in gs["候補"]:
        if g not in present:
            rows.append({"gene": g, "symbol": sym.get(g, g), "状態": "結果が無い"})
            continue
        zs = [float(long[(long.gene == g) & (long.seed == s)].z.iloc[0]) for s in seeds]
        all_sig = all(z > z_floor for z in zs)
        mean_z = float(np.mean(zs))
        verdict = ("確認された" if (all_sig and (smoke or mean_z >= z_hold)) else
                   "判定保留（bfloat16 の分解能の内側）" if (all_sig and not smoke) else
                   "確認されなかった")
        rows.append({"gene": g, "symbol": sym.get(g, g), "z": [round(z, 3) for z in zs],
                     "z の最小": round(min(zs), 3), "z の平均": round(mean_z, 3),
                     "3 seed すべて z > 2.73": bool(all_sig),
                     "3 seed 平均 z >= 3.6": bool(mean_z >= z_hold), "判定": verdict})
    out["4. V-1/V-1b/V-1c 候補 8 個"] = {"遺伝子ごと": rows,
                                         "確認された": [r["symbol"] for r in rows if r.get("判定") == "確認された"],
                                         "判定保留": [r["symbol"] for r in rows
                                                      if r.get("判定", "").startswith("判定保留")],
                                         "確認されなかった": [r["symbol"] for r in rows
                                                              if r.get("判定") == "確認されなかった"]}

    # ---------------------------------------------------------- 5. V-2 偶然の水準（対照 549 で実測）
    nl_tab = long[long.gene.isin(nulls_549)]
    def v1_pass(g: str) -> bool:
        zs = long[long.gene == g].set_index("seed").reindex(seeds).z.to_numpy()
        return bool(np.all(zs > z_floor) and (smoke or np.mean(zs) >= z_hold))
    null_genes = sorted(set(nl_tab.gene))
    n_v1_null = sum(v1_pass(g) for g in null_genes)
    p_hat = n_v1_null / len(null_genes) if null_genes else float("nan")
    se_p = (p_hat * (1 - p_hat) / len(null_genes)) ** 0.5 if null_genes else float("nan")
    out["5. V-2 偶然の水準"] = {
        "対照の個数": len(null_genes), "V-1 を満たした対照": n_v1_null,
        "p（実測の片側割合）": round(p_hat, 5), "p の標準誤差": round(se_p, 5),
        "候補8での期待件数": round(p_hat * 8, 4),
        "候補8での実際の合格数": len(out["4. V-1/V-1b/V-1c 候補 8 個"]["確認された"]),
        "P(1 件以上 | 偶然)": round(1 - (1 - p_hat) ** 8, 5) if p_hat == p_hat else None}

    # D-2: 主対照498 と 低発現40 を分けて偽陽性率を実測
    def rate(genes: list[str]) -> dict:
        gg = [g for g in genes if g in present]
        k = sum(v1_pass(g) for g in gg)
        return {"n": len(gg), "k": k, "rate": round(k / len(gg), 5) if gg else None}
    r498, r40 = rate(gs["主対照498"]), rate(gs["低発現40"])
    out["6.8 D-2 偽陽性率（主対照498 対 低発現40）"] = {
        "主対照498": r498, "低発現40": r40,
        "比（低発現 ÷ 主対照）": (round(r40["rate"] / r498["rate"], 3)
                                   if r498["rate"] else None) if r498.get("rate") else None,
        "注": "低発現40 は n=40 なので標準誤差は約 1.6%。この診断は粗い（止める条件にはしない）"}

    # ---------------------------------------------------------- 6. V-3a/b/c・V-4
    rng = np.random.default_rng(20261004)
    n_draw = 2000 if smoke else 1000
    zpos71, g71 = zmat(gs["陽性対照71"])
    znull, gnull = zmat(null_genes)
    def sham_pct(zset: np.ndarray, k: int) -> list[float]:
        nd = set_null(znull, k, n_draw, rng, replace=True)
        return [float(np.quantile(nd[:, i], 0.99)) for i in range(len(seeds))]
    v3a_obs = zpos71.mean(axis=0) if zpos71.size else np.array([np.nan] * len(seeds))
    v3a_thr = sham_pct(znull, len(g71)) if g71 else [np.nan] * len(seeds)
    v3a_pass = bool(zpos71.size and all(o > t for o, t in zip(v3a_obs, v3a_thr)))
    g69 = [g for g in g71 if sym.get(g) not in ("GSN", "ELN")]
    z69, _ = zmat(g69)
    v3b_obs = z69.mean(axis=0) if z69.size else np.array([np.nan] * len(seeds))
    v3b_thr = sham_pct(znull, len(g69)) if g69 else [np.nan] * len(seeds)
    v3b_pass = bool(z69.size and all(o > t for o, t in zip(v3b_obs, v3b_thr)))
    k71 = sum(v1_pass(g) for g in g71)
    pval71 = float(stats.binomtest(k71, len(g71), max(p_hat, 1e-9), alternative="greater").pvalue) \
        if g71 and p_hat == p_hat else None
    out["7. V-3 ものさし（陽性対照）"] = {
        "V-3a（71個）": {"set_mean_z": [round(float(x), 3) for x in v3a_obs],
                         "偽セット99%点": [round(float(x), 3) for x in v3a_thr], "pass": v3a_pass},
        "V-3b（69個。GSN・ELN除く）": {"set_mean_z": [round(float(x), 3) for x in v3b_obs],
                                       "偽セット99%点": [round(float(x), 3) for x in v3b_thr], "pass": v3b_pass},
        "V-3c（二項検定）": {"k": k71, "n": len(g71), "p（対照で実測）": round(p_hat, 5) if p_hat == p_hat else None,
                             "期待件数": round(len(g71) * p_hat, 3) if p_hat == p_hat else None,
                             "片側p値": pval71, "pass": bool(pval71 is not None and pval71 < 0.05)}}
    v3_overall_pass = bool(v3a_pass)  # V-3a が主。不合格なら候補を解釈しない（plan.md §6.3）

    hk3 = [g for g in cfg["controls"]["negative_fixed"] if g in present]
    out["8. V-4 中立性（恒常発現 ACTB/GAPDH/B2M）"] = {
        sym.get(g, g): {"z": [round(float(long[(long.gene == g) & (long.seed == s)].z.iloc[0]), 3)
                              for s in seeds], "V-1 を満たすか": v1_pass(g)} for g in hk3}

    # ---------------------------------------------------------- 9. §6.6 負の側の目盛り
    n_neg_null = sum(bool(np.all(long[long.gene == g].set_index("seed").reindex(seeds).z.to_numpy()
                                 < -z_floor)) for g in null_genes)
    p_neg = n_neg_null / len(null_genes) if null_genes else float("nan")
    out["9. 負の側の目盛り（対照549で実測）"] = {
        "正の側 p": round(p_hat, 5) if p_hat == p_hat else None,
        "負の側 p": round(p_neg, 5) if p_neg == p_neg else None,
        "同じ桁か": bool(p_neg == p_neg and p_hat == p_hat and p_hat > 0
                         and 0.1 <= (p_neg / p_hat if p_hat else 0) <= 10),
        "読み方": "桁が違う（10分の1未満）なら負の値は記録するだけで解釈しない（plan.md §6.6）"}

    # ---------------------------------------------------------- 10. §6.7 旧の測り方との突き合わせ（記述）
    old_csvs = [P02 / f"results/gene_seed_007-{s}.csv" for s in seeds]
    if all(p.exists() for p in old_csvs) and not smoke:
        old = pd.concat([pd.read_csv(p) for p in old_csvs], ignore_index=True)
        cmp_genes = sorted((set(gs["主対照498"]) | set(gs["陽性対照71"])) & set(old.gene) & present)
        rows_cmp = []
        for s in seeds:
            on = old[old.seed == s].set_index("gene").z.reindex(cmp_genes)
            nn = long[long.seed == s].set_index("gene").z.reindex(cmp_genes)
            d = (nn - on)
            rho = float(stats.spearmanr(on.dropna(), nn.reindex(on.dropna().index)).statistic)
            rows_cmp.append({"seed": s, "遺伝子数": int(d.notna().sum()),
                             "z の差の中央": round(float(d.median()), 4),
                             "z の差の標準偏差": round(float(d.std()), 4),
                             "新旧の順位相関(Spearman)": round(rho, 4)})
        out["10. 旧の測り方との突き合わせ（記述・探索。判定には使わない）"] = {
            "比べた遺伝子数": len(cmp_genes), "seed ごと": rows_cmp,
            "注": "CLAUDE.md §3.5 により探索であることを明記する。確認は別データ/別seedで行う"}

    # ---------------------------------------------------------- 11. 解釈の規則 R-1〜R-7
    axis_a = [sym[g] for g in gs["候補"] if sym.get(g) in ("MADCAM1", "CXCL13", "CCL21")]
    axis_b = [sym[g] for g in gs["候補"] if sym.get(g) not in ("MADCAM1", "CXCL13", "CCL21")]
    conf = set(out["4. V-1/V-1b/V-1c 候補 8 個"]["確認された"])
    hold = set(out["4. V-1/V-1b/V-1c 候補 8 個"]["判定保留"])
    a_hit, b_hit = bool(conf & set(axis_a)), bool(conf & set(axis_b))
    if not v3_overall_pass:
        rule, tag = "R-6: V-3a 不合格。判定不能。候補の値は報告するが解釈しない", "R-6"
    elif a_hit and b_hit:
        rule, tag = "R-1: 摂動に使う細胞を遺伝子によらずそろえても上位候補は残った（軸A・軸Bの両方で確認）", "R-1"
    elif a_hit:
        rule, tag = "R-2: 三次リンパ組織側（軸A）に限って確認された", "R-2"
    elif b_hit:
        rule, tag = "R-3: 軸Bのみ確認。三次リンパ組織の3個は共通細胞集合では最小検出可能効果量を超えなかった", "R-3"
    elif hold:
        rule, tag = "R-4: bfloat16 では決められない（判定保留のみ）。float32 での測り直しを PI に諮る", "R-4"
    elif not conf and not hold:
        rule, tag = "R-5: 候補はどれも確認/保留に至らず。『効果がない』とは書けない", "R-5"
    else:
        rule, tag = "R-7: 事前登録のどの型にも当たらない。判定を保留し事実のみ報告", "R-7"
    out["11. 解釈の規則（結果を見る前に固定したもの）"] = {
        "軸A(三次リンパ組織3個)": axis_a, "軸B(それ以外5個)": axis_b,
        "軸Aで確認": a_hit, "軸Bで確認": b_hit, "該当": tag, "読み方": rule}

    # ---------------------------------------------------------- 12. §6.9 626遺伝子の順位表（記述）
    rank = (long.groupby("gene").z.mean().rename("z_mean3seed").reset_index()
            .assign(symbol=lambda d: d.gene.map(sym)).sort_values("z_mean3seed", ascending=False))
    rank_path = (res / "gene_rank_009-smoke.csv") if smoke else (E / "results/gene_rank_009.csv")
    rank.to_csv(rank_path, index=False)
    out["12. 全遺伝子の順位表"] = {"行数": len(rank), "保存先": str(rank_path),
                                   "上位10": rank.head(10)[["symbol", "gene", "z_mean3seed"]].to_dict("records")}

    out["判定の提案"] = rule
    out["注"] = ("判定（ゲート G3）は PI が行う。この道具は提案までにとどまる。"
                 "G1（仮説採択）は立てない（PI 2026-10-04）。ゲート G1/G3 の対象ではないので、"
                 "hypotheses/tracker.md は変更しない。"
                 "これはコンピュータ上の遺伝子摂動の結果であり、実験的に確認された生物学的事実ではない。"
                 "盲検中の TP53・SPP1・LCN2 は摂動も解析もしていない。")

    (E / "results").mkdir(parents=True, exist_ok=True)
    name = "criteria_009-smoke.json" if smoke else "criteria_009.json"
    (E / "results" / name).write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps({k: out[k] for k in out if k not in ("12. 全遺伝子の順位表",)},
                     ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true", help="試し実行の出力で評価コードを通す（**判定ではない**）")
    ap.add_argument("--self-test", action="store_true", help="層別補正の自己検査（データ不要）")
    args = ap.parse_args()
    if args.self_test:
        self_test()
    else:
        main(smoke=args.smoke)

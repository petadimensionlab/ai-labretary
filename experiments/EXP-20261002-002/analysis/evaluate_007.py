#!/usr/bin/env python3
"""EXP-20261002-002 — 既存 717 遺伝子 + 追加の対照 500 個を合算し、対照 1,000 個で評価する。

## なぜ「表に 500 行足すだけ」では済まないか（**重要**）

対照プールは z の標準化だけでなく、**遺伝子ごとの残差 D_raw そのもの**に入っている。

  - `gene_stats()` は、対照遺伝子**全体**で当てた二次回帰（自分自身は除く）からの残差を D_raw にする
  - `eligibility()` は、元の順位と入力の長さで比べられる対照だけを選んで順位 q を作る
  - `z` は対照の平均と標準偏差で標準化する

したがって対照が 500 → 1,000 になると、**すべての遺伝子の D_raw・q・z が動く**。
`gene_seed_005.csv` に追加分の行を連結するのでは足りない。
**細胞ごとの表から作り直す**（摂動は回さない。保存済みの埋め込みを読むだけ）。

## 手順

 0. 状態の重心が EXP-20261001-001 と**同一**であることを確かめる（違えば**止める**）
 1. 既存 717 + 追加 500 の細胞ごとの表を作る（保存済みの埋め込みから）
 2. **再現の検査**: 対照を既存の 500 個だけに戻して評価し、`criteria_005.json` と
    一致するか見る。一致しなければ**止める**（合算の仕方が間違っている）
 3. 対照 1,000 個で S-1〜S-4 と主判定 P-1 を出す
 4. 主判定の値が 500 個のときからどれだけ動いたかを記録する

走らせ方:
    .venv/bin/python experiments/EXP-20261002-002/analysis/evaluate_007.py
  試し実行の出力で通すとき（**判定ではない**）:
    .venv/bin/python experiments/EXP-20261002-002/analysis/evaluate_007.py --smoke
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261002-002"
P01 = ROOT / "experiments/EXP-20261001-001"
sys.path.insert(0, str(ROOT / "scripts/analysis"))
sys.path.insert(0, str(P01 / "analysis"))
from movement_correction import gene_stats                      # noqa: E402
from bias_budget import bias_budget                             # noqa: E402
from evaluate_005 import cos, eligibility, set_null, sham_sets   # noqa: E402

BLOCKED = {"ENSG00000141510": "TP53", "ENSG00000118785": "SPP1", "ENSG00000148346": "LCN2"}


# --------------------------------------------------------------------------- 細胞ごとの表
def load_roots(cfg: dict, roots: list[tuple[str, Path]], tok_path: Path) -> tuple[pd.DataFrame, dict]:
    """複数の結果ディレクトリから細胞ごとの表を作る。

    evaluate_005.load() と**同じ計算**を、結果の置き場所が 2 か所ある場合に広げたもの。
    主指標は向きへの射影、移動量 0 の細胞は外す（どちらも EXP-20261001-001 と同一）。
    """
    import pickle
    from datasets import load_from_disk
    from geneformer import TOKEN_DICTIONARY_FILE

    parts = []
    for label, res in roots:
        c = pd.read_csv(res / "cells.csv", dtype={"donor_id": str})
        c["root_label"] = label
        c["root_dir"] = str(res)
        parts.append(c)
    c = pd.concat(parts, ignore_index=True)
    n0 = len(c)
    c = c[~c.gene.isin(BLOCKED)].reset_index(drop=True)
    if n0 != len(c):
        raise SystemExit(f"回さないはずの遺伝子が {n0 - len(c)} 行あった。実行設定を確認すること")
    c["seed"] = c.run.astype(str).str.extract(r"run-(\d+)", expand=False).astype(int)
    dup = c.groupby(["seed", "gene"]).root_label.nunique()
    if (dup > 1).any():
        raise SystemExit(f"同じ seed・遺伝子が 2 か所の結果にある: {dup[dup > 1].index.tolist()[:5]}")

    tokd = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
    ds = load_from_disk(str(tok_path))
    uids = sorted(c.cell_uid.unique())
    ids = {u: np.asarray(x) for u, x in zip(uids, ds.select(uids)["input_ids"])}
    pos = np.empty(len(c), dtype=np.int64)
    for gene, part in c.groupby("gene", sort=False):
        t = tokd[gene]
        pos[part.index.to_numpy()] = [int(np.flatnonzero(ids[u] == t)[0]) for u in part.cell_uid]
    c["pos"] = pos
    c["lp"] = np.log2(1 + c.pos.astype(float))
    c["ll"] = np.log2(c.length.astype(float))
    c["relpos"] = c.pos / c.length

    s_ = cfg["states"]
    names = {"goal": s_["goal_state"], "start": s_["start_state"]}
    names |= {a: a for a in s_.get("alt_states", [])}
    vals = {n: np.full(len(c), np.nan) for n in names}
    dnorm = np.full(len(c), np.nan)
    proj = np.full(len(c), np.nan)
    for (rd, run, gene), part in c.groupby(["root_dir", "run", "gene"], sort=False):
        res = Path(rd)
        sz = np.load(res / run / "state_emb_vectors.npz")
        gd = res / run / "genes" / gene / "emb"
        o = np.load(gd / "orig.npy").astype(np.float64)
        p = np.load(gd / "pert.npy").astype(np.float64)
        if len(o) != len(part):
            raise SystemExit(f"{gd}: 埋め込み {len(o)} ≠ cells.csv {len(part)}")
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


def check_centroids(roots: list[tuple[str, Path]], seeds: list[int], tol: float) -> dict:
    """**状態の重心が 2 か所で同一か。** 違えば合算してはいけない（別の尺度で測ったことになる）。"""
    out = {}
    base = dict(roots)[list(dict(roots))[0]]
    for s in seeds:
        run = f"run-{s}-L-1"
        ref = np.load(base / run / "state_emb_vectors.npz")
        for label, res in roots[1:]:
            p = res / run / "state_emb_vectors.npz"
            if not p.exists():
                out[f"seed{s}/{label}"] = {"状態": "ファイルが無い", "path": str(p)}
                continue
            z = np.load(p)
            diffs = {k: float(np.max(np.abs(ref[k].astype(float) - z[k].astype(float))))
                     for k in ref.files}
            out[f"seed{s}/{label}"] = {"重心ごとの最大差": diffs,
                                       "最大差": max(diffs.values()),
                                       "一致": bool(max(diffs.values()) <= tol)}
    return out


# --------------------------------------------------------------------------- 本体
def main(smoke: bool = False) -> None:
    cfg = yaml.safe_load(open(E / "config.yaml"))
    mg, bb = cfg["merge"], cfg["bias_budget"]
    src_cfg = yaml.safe_load(open(P01 / "config.yaml"))
    a = src_cfg["analysis"]
    cr, ABS = a["criteria"], a["set_abs_min"]
    mc = cfg["metric"]["confound_correction"]

    # ---- 事前登録どおりの設定で動いているかの検査（違えば止める）
    assert mc["variables"] == ["lp", "ll"], f"補正の版が事前登録と違う（{mc['variables']}）"
    assert mc["drop_zero_movement"] is True, "移動量 0 の細胞を外す設定になっていない"
    assert cfg["metric"]["primary_form"].startswith("projection_"), "主指標が射影になっていない"
    assert bb["budget"] == 0.10 and bb["decision_line"] == 0.41, \
        "歪みの予算が EXP-20261002-001 と違う（**動かしてはいけない**）"
    assert bb["covariates"] == ["dnorm_med", "lp_med", "ll_med", "relpos_med", "uniform"], \
        "交絡の候補が EXP-20261002-001 と違う"

    gj = json.loads((ROOT / mg["existing"]["genes"]).read_text())["genes"]
    ens = {g["symbol"]: g["ensembl"] for g in gj}
    bias_goal = {g["symbol"]: g["bias_goal"] for g in gj}
    subs = json.loads((ROOT / mg["existing"]["subsets"]).read_text())["subsets"]
    sets_def = {k: [ens[x["symbol"]] for x in v] for k, v in subs.items()}
    null500 = [ens[x["symbol"]] for x in
               json.loads((ROOT / mg["existing"]["null_pool"]).read_text())["versions"]["B_500"]]
    xtra = json.loads((ROOT / mg["extra_null"]).read_text())
    extra_null = [x["ensembl"] for x in xtra["追加する対照"]]
    ens |= {x["symbol"]: x["ensembl"] for x in xtra["追加する対照"]}
    bias_goal |= {x["symbol"]: x["bias_goal"] for x in xtra["追加する対照"]}
    sym = {v: k for k, v in ens.items()}
    assert not (set(null500) & set(extra_null)), "既存と追加の対照が重なっている"

    new_res = E / ("results-smoke" if smoke else cfg["output"]["results_dir"])
    roots = [("既存(EXP-20261001-001)", ROOT / mg["existing"]["results_dir"]),
             ("追加(EXP-20261002-002)", new_res)]
    seeds = [cfg["robustness"]["seeds"][0]] if smoke else list(cfg["robustness"]["seeds"])

    out: dict = {
        "実験": "EXP-20261002-002",
        "種類": "追加の対照を足して歪みの見積もりの区間を縮める",
        "mode": "smoke" if smoke else "final",
        "seed": seeds,
    }
    if smoke:
        out["注（試し実行）"] = (
            "**試し実行の出力に対する評価であり、判定ではない。** 追加の対照はごく少数で、"
            "状態の重心も本番と違う条件で作っている。目的は「評価コードが最後まで通ること」の"
            "確認だけである（plan.md §10）。")

    # ---------------------------------------------------------- 0. 状態の重心の一致
    cen = check_centroids(roots, seeds, mg["centroid_max_abs_diff"])
    out["0. 状態の重心が既存と同一か"] = cen
    bad_cen = [k for k, v in cen.items() if not v.get("一致", False)]
    if bad_cen and not smoke:
        (E / "results").mkdir(exist_ok=True)
        (E / "results/criteria_007.json").write_text(
            json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
        raise SystemExit(f"**止める**: 状態の重心が既存と違う {bad_cen}。別の尺度で測っている")

    # ---------------------------------------------------------- 1. 細胞ごとの表
    c, drop = load_roots(cfg, roots, ROOT / cfg["data"]["tokenized_dataset"])
    c = c[c.seed.isin(seeds)].reset_index(drop=True)
    uni = cfg["metric"]["uniform_component_states"]
    c["m_goal"] = c.m_proj
    c["m_uniform"] = c[[f"s_{x}" for x in uni]].mean(axis=1)
    present = set(c.gene.unique())
    out["1. 細胞ごとの表"] = {
        "遺伝子数": len(present), "移動量 0 で外した細胞": drop,
        "結果の出どころごとの遺伝子数":
            {k: int(c[c.root_label == k].gene.nunique()) for k in c.root_label.unique()}}

    pm = a["pos_matched_permutation"]
    min_elig = 1 if smoke else pm["min_eligible_null_genes"]

    def build(nulls: list[str]) -> dict[int, pd.DataFrame]:
        """対照プールを与えて、seed ごとの遺伝子ごとの表（z と交絡の候補）を作る。"""
        nl_set = set(nulls) & present
        tabs = {}
        for s in seeds:
            cs = c[c.seed == s]
            nl = {h: hd for h, hd in cs[cs.gene.isin(nl_set)].groupby("gene")}
            el = eligibility(cs, nl, pm["tolerance_log2"], pm["min_cover_frac"])
            f = gene_stats(cs, nl, el, "m_goal", min_elig, use_dnorm=False)
            u = gene_stats(cs, nl, el, "m_uniform", min_elig, use_dnorm=False)
            n = f[f.gene.isin(nl_set)]
            t = pd.DataFrame({
                "gene": f.gene.to_numpy(),
                "z": ((f.D_raw - n.D_raw.mean()) / n.D_raw.std(ddof=1)).to_numpy(),
                "q": f.q.to_numpy(), "n_cells": f.n_cells.to_numpy(),
                "dnorm_med": f.dnorm_med.to_numpy(), "lp_med": f.lp_med.to_numpy(),
                "ll_med": f.ll_med.to_numpy(), "relpos_med": f.relpos_med.to_numpy(),
                "uniform": u.set_index("gene").D_raw.reindex(f.gene).to_numpy()})
            tabs[s] = t.reset_index(drop=True)
        return tabs

    def evaluate(tabs: dict[int, pd.DataFrame], nulls: list[str], label: str) -> dict:
        """セットの規則・S-1〜S-4・主判定を、与えた対照プールで出す。"""
        nl = [x for x in nulls if x in present]
        Z = {s: tabs[s].set_index("gene").z for s in seeds}
        rng = np.random.default_rng(cr["null_draw_seed"])
        znull = np.column_stack([Z[s].reindex(sorted(nl)).to_numpy() for s in seeds])
        n_perm = 2000 if smoke else cr["n_permutation"]

        def set_test(genes, name):
            gs = [x for x in genes if x in Z[seeds[0]].index]
            obs = np.array([float(np.mean([Z[s][x] for x in gs])) for s in seeds])
            nd = set_null(znull, len(gs), n_perm, rng)
            thr = np.array([float(np.quantile(nd[:, i], cr["set_percentile"]))
                            for i in range(len(seeds))])
            am = ABS.get(name, -1e9)
            return {"k": len(gs), "set_mean_z": [round(float(x), 3) for x in obs],
                    "threshold": [round(float(x), 3) for x in thr], "abs_min": am,
                    "pass": bool(all(o > t for o, t in zip(obs, thr))
                                 and all(o >= am for o in obs))}

        r: dict = {"対照の個数": len(nl)}
        r["P-1 主判定"] = P1 = set_test(sets_def["INT"], "INT")
        r["P-4 陽性対照"] = P4 = set_test(sets_def["positive_control71"], "positive_control71")
        neu = [x for x in sets_def["neutral_adopted"] if x in Z[seeds[0]].index]
        obs_n = np.array([float(np.mean([Z[s][x] for x in neu])) for s in seeds])
        nd_n = set_null(znull, len(neu), n_perm, rng)
        pct_n = [float((nd_n[:, i] < obs_n[i]).mean()) for i in range(len(seeds))]
        r["P-5 中立性"] = P5 = {
            "k": len(neu), "set_mean_z": [round(float(x), 3) for x in obs_n],
            "null_percentile": [round(x, 4) for x in pct_n],
            "n_seeds_in_range": int(sum(cr["neutral_range"][0] <= x <= cr["neutral_range"][1]
                                        for x in pct_n))}
        P5["pass"] = bool(P5["n_seeds_in_range"] >= cr["neutral_min_seeds"])

        # ---- S-1 介入セットの歪みの予算
        S1 = {str(s): bias_budget(tabs[s], sets_def["INT"], nl, covariates=bb["covariates"],
                                  budget=bb["budget"], decision_line=bb["decision_line"])
              for s in seeds}
        r["S-1 歪みの予算"] = {
            "per_seed": S1, "B": [S1[str(s)]["B_joint"] for s in seeds],
            "se_B": [S1[str(s)]["se_B"] for s in seeds],
            "上側": [S1[str(s)]["upper_joint"] for s in seeds], "予算": bb["budget"],
            "pass": bool(all(S1[str(s)]["pass"] for s in seeds))}

        # ---- S-2 偽のセットの二項検定
        def tbl(genes):
            t = tabs[seeds[0]].set_index("gene")
            t = t.loc[[x for x in genes if x in t.index]].copy()
            t["symbol"] = [sym[x] for x in t.index]
            t["bias_goal"] = [bias_goal[s] for s in t.symbol]
            return t.reset_index()

        # 単一遺伝子の偶然合格率を、**この対照プールで測り直す**
        nd_c = set_null(znull, len(sets_def["INT"]), 2000 if smoke else cr["n_chance_draw"], rng)
        thr_c = [float(np.quantile(set_null(znull, len(sets_def["INT"]), n_perm, rng)[:, i],
                                   cr["set_percentile"])) for i in range(len(seeds))]
        hit = np.ones(len(nd_c), bool)
        for i in range(len(seeds)):
            hit &= (nd_c[:, i] > thr_c[i]) & (nd_c[:, i] >= ABS["INT"])
        p_chance = float(hit.mean())
        n_sham = 5 if smoke else cr["n_sham_sets"]
        shams = sham_sets(tbl(nl), tbl(sets_def["INT"]), len(sets_def["INT"]), n_sham,
                          np.random.default_rng(cr["null_draw_seed"] + 1), cr["sham_match_bias"])
        sz, k_pass = [], 0
        for ss in shams:
            rr = set_test([ens[x] for x in ss], "INT")
            sz.append(rr["set_mean_z"]); k_pass += int(rr["pass"])
        pval = float(stats.binomtest(k_pass, len(shams), max(p_chance, 1e-9),
                                     alternative="greater").pvalue)
        r["S-2 偽のセットの二項検定"] = {
            "N": len(shams), "p（単一の偶然合格率。この対照プールで実測）": round(p_chance, 6),
            "期待合格数": round(len(shams) * p_chance, 4), "k（実際の合格数）": k_pass,
            "片側 p 値": round(pval, 4),
            "偽セットの平均 z": np.round(np.array(sz, float).mean(axis=0), 4).tolist(),
            "pass": bool(pval >= bb["sham_pvalue_min"])}

        # ---- S-3 陽性対照（歪みを引いても立つか）
        S3 = {str(s): bias_budget(tabs[s], sets_def["positive_control71"], nl,
                                  covariates=bb["covariates"], budget=bb["budget"],
                                  decision_line=bb["decision_line"]) for s in seeds}
        adj = [round(abs(P4["set_mean_z"][i]) - S3[str(s)]["upper_joint"], 3)
               for i, s in enumerate(seeds)]
        r["S-3 陽性対照"] = {
            "セット平均 z": P4["set_mean_z"], "規則の合否": P4["pass"],
            "歪みの上側": [S3[str(s)]["upper_joint"] for s in seeds],
            "歪みを引いた値": adj, "必要な値": ABS["positive_control71"],
            "pass": bool(P4["pass"] and all(x >= ABS["positive_control71"] for x in adj))}
        r["S-4 中立性"] = P5

        # ---- 記述のみ（判定に使わない）
        r["記述のみ（判定に使わない）"] = {
            str(s): {
                "M1 較正（上裾 / 下裾）": [
                    round(float((tabs[s][tabs[s].gene.isin(set(nl))].q >= 0.90).mean()), 4),
                    round(float((tabs[s][tabs[s].gene.isin(set(nl))].q < 0.10).mean()), 4)],
                **{f"M{n} {jp}": round(float(stats.spearmanr(
                    tabs[s][tabs[s].gene.isin(set(nl))][col],
                    tabs[s][tabs[s].gene.isin(set(nl))].z).statistic), 4)
                   for n, jp, col in [(2, "元の順位", "lp_med"), (2, "入力の長さ", "ll_med"),
                                      (2, "相対位置", "relpos_med"), (4, "一様成分", "uniform"),
                                      (5, "移動量", "dnorm_med")]}}
            for s in seeds}

        health = all(r[k]["pass"] for k in ("S-1 歪みの予算", "S-2 偽のセットの二項検定",
                                            "S-3 陽性対照", "S-4 中立性"))
        r["健全性"] = {k.split()[0]: r[k]["pass"] for k in
                       ("S-1 歪みの予算", "S-2 偽のセットの二項検定",
                        "S-3 陽性対照", "S-4 中立性")} | {"すべて合格": health}
        floor = json.loads((ROOT / mg["existing"]["power"]).read_text())[
            "sets"]["INT"]["detection_floor_with_gene_sd"]
        r["検出下限"] = floor
        r["判定の提案"] = ("判定不能（健全性の検査が不合格）" if not health
                           else "支持" if P1["pass"]
                           else f"否定（検出下限 {floor} 標準偏差以上の効果は認められなかった）")
        return r

    # ---------------------------------------------------------- 2. 再現の検査（**止まりうる**）
    tabs500 = build(null500)
    base = evaluate(tabs500, null500, "既存の対照 500 個")
    ref = json.loads((ROOT / mg["existing"]["criteria"]).read_text())
    # 前回の値は seed 81/82/83 の順に並んでいる。**試し実行では seed を 1 つしか回さない**ので、
    # 回した seed の位置だけを取り出して比べる（位置をそろえずに比べると必ず不一致になる）
    pos = [list(cfg["robustness"]["seeds"]).index(s) for s in seeds]
    pick = lambda v: [v[i] for i in pos]  # noqa: E731
    rep = {"主判定 P-1": {"今回": base["P-1 主判定"]["set_mean_z"],
                          "前回": pick(ref["P1_primary"]["set_mean_z"])},
           "陽性対照 P-4": {"今回": base["P-4 陽性対照"]["set_mean_z"],
                            "前回": pick(ref["P4_positive_control"]["set_mean_z"])},
           "中立性 P-5": {"今回": base["P-5 中立性"]["set_mean_z"],
                          "前回": pick(ref["P5_neutral"]["neutral_adopted"]["set_mean_z"])}}
    tol = mg["reproduction_tolerance"]
    bad = [k for k, v in rep.items()
           if max(abs(np.array(v["今回"], float) - np.array(v["前回"], float))) > tol]
    rep |= {"許容差": tol, "一致": not bad,
            "意味": "追加の対照を混ぜた表から、対照を既存の 500 個だけに戻して評価すると、"
                    "EXP-20261001-001 の確定値が再現されるはず。再現できなければ合算の仕方が誤っている"}
    out["2. 再現の検査（対照を既存の 500 個に戻す）"] = rep
    out["2b. 対照 500 個での健全性（参考。EXP-20261002-001 と同じはず）"] = {
        "S-1 上側": base["S-1 歪みの予算"]["上側"], "S-1 B": base["S-1 歪みの予算"]["B"],
        "健全性": base["健全性"]}
    if bad and not smoke:
        (E / "results").mkdir(exist_ok=True)
        (E / "results/criteria_007.json").write_text(
            json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
        raise SystemExit(f"**止める**: 既存の値を再現できなかった: {bad}")

    # ---------------------------------------------------------- 3. 対照 1,000 個での評価（**本番**）
    nulls_all = null500 + extra_null
    tabs_all = build(nulls_all)
    main_r = evaluate(tabs_all, nulls_all, "対照 1,000 個")
    out["3. 対照を増やした評価（**これが判定に使う値**）"] = main_r

    # ---------------------------------------------------------- 4. 主判定がどれだけ動いたか
    def mv(key, sub="set_mean_z"):
        x = np.array(base[key][sub], float)
        y = np.array(main_r[key][sub], float)
        return {"対照 500": x.round(3).tolist(), "対照 1,000": y.round(3).tolist(),
                "動いた量": np.round(y - x, 4).tolist(),
                "動いた量の最大": round(float(np.max(np.abs(y - x))), 4)}

    out["4. 対照を増やしたことで動いた量"] = {
        "主判定 P-1": mv("P-1 主判定"), "陽性対照 P-4": mv("P-4 陽性対照"),
        "中立性 P-5": mv("P-5 中立性"),
        "S-1 歪みの点推定 B": mv("S-1 歪みの予算", "B"),
        "S-1 標準誤差 se_B": mv("S-1 歪みの予算", "se_B"),
        "S-1 上側": mv("S-1 歪みの予算", "上側"),
        "注": "**主判定の値は対照の変更で動く。** 動くことは計画 §2 に事前に書いてある。"
              "判定が覆るのは絶対値の条件 0.41 を超えたときだけで、"
              "その余裕は計画 §3.4 に見積もってある"}
    out["判定の提案"] = main_r["判定の提案"]
    out["注"] = ("判定は PI（関門 G3）が行う。この道具は提案までにとどまる。"
                 "また、これはコンピュータ上の遺伝子摂動の結果であり、"
                 "実験的に確認された生物学的事実ではない。")

    (E / "results").mkdir(exist_ok=True)
    name = "criteria_007-smoke.json" if smoke else "criteria_007.json"
    (E / "results" / name).write_text(
        json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    for s in seeds:
        tabs_all[s].assign(seed=s).to_csv(
            E / "results" / f"gene_{'smoke' if smoke else 'seed'}_007-{s}.csv", index=False)
    print(json.dumps({k: out[k] for k in out if k.startswith(("2.", "3.", "4.", "判定"))
                      and k != "3. 対照を増やした評価（**これが判定に使う値**）"}
                     | {"3. 健全性": main_r["健全性"],
                        "3. S-1 上側": main_r["S-1 歪みの予算"]["上側"],
                        "3. 主判定": main_r["P-1 主判定"]["set_mean_z"]},
                     ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="試し実行の出力で評価コードを通す（**判定ではない**）")
    main(**vars(ap.parse_args()))

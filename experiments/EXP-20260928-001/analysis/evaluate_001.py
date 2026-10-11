"""EXP-20260928-001 — 事前登録した判定の機械的な当てはめ（plan.md §6 をそのまま実装する）。

EXP-20260927-004 からの変更点:
  判定の単位を「1 遺伝子」から「遺伝子セット」へ変える。統計量はセットに属する遺伝子の
  標準偏差単位 z の平均。対照は「対照遺伝子から無作為に取った同じ大きさのセット」。
  さらに軸 × セットの交差（目標マーカーは目標の軸で、筋線維芽細胞マーカーは筋線維芽細胞の軸で
  大きくなること）を必須条件にする。

判定対象にしないもの（読み込み直後に除外し、1 行でも残っていれば停止する）:
  主対象 LTBR / RELB / NFKB2、参考対象 CDKN1A、盲検の TP53 / SPP1 / LCN2、
  目標ラベルを定義する CCL19 / CCL21 / CXCL13。

出力: results/eval/criteria_001.json, results/eval/gene_seed_001.csv, results/eval/set_null_001.csv
実行: .venv/bin/python experiments/EXP-20260928-001/analysis/evaluate_001.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260928-001"

BLOCKED = {"ENSG00000111321": "LTBR", "ENSG00000104856": "RELB", "ENSG00000077150": "NFKB2",
           "ENSG00000124762": "CDKN1A", "ENSG00000141510": "TP53",
           "ENSG00000118785": "SPP1", "ENSG00000148346": "LCN2",
           "ENSG00000172724": "CCL19", "ENSG00000137077": "CCL21", "ENSG00000156234": "CXCL13"}


# ---------------------------------------------------------------- 指標
def cos(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return (a @ b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b) + 1e-12)


def metric_columns(c: pd.DataFrame) -> dict[str, str]:
    """4 つの重心から、同じ関数形・同じ尺度で軸だけを入れ替えた指標を作る。
    A12 の形（目標に特異的な成分 / 移動量）。EXP-20260927-004 と同じ定義。"""
    g, st, al, my, dn = c.s_goal, c.s_start, c.s_alt, c.s_myof, c.dnorm + 1e-30

    def form(target, others):
        return (target - sum(others) / len(others)) / dn

    c["m_goal"] = form(g, [st, al, my])       # 主軸: 目標に特異的な成分
    c["m_myofaxis"] = form(my, [st, al, g])   # 第 2 軸: 筋線維芽細胞に特異的な成分
    c["m_placebo"] = form(al, [st, g, my])    # 置き換え（対照の状態）
    c["m_uniform"] = (g + st + al + my) / 4.0
    c["m_A11"] = g - (st + al + my) / 3.0
    return {"goal": "m_goal", "myofaxis": "m_myofaxis", "placebo": "m_placebo",
            "uniform": "m_uniform", "A11": "m_A11"}


SCALE_GROUP = {"m_goal": "specific", "m_myofaxis": "specific", "m_placebo": "specific",
               "m_uniform": "raw_shift", "m_A11": "specific_unnormalized"}


def compare_same_scale(a: str, b: str) -> None:
    ga, gb = SCALE_GROUP.get(a), SCALE_GROUP.get(b)
    assert ga is not None and gb is not None, f"尺度の登録がない列を比較している: {a}, {b}"
    assert ga == gb, f"尺度の違う量を比較しようとしている: {a}（{ga}）と {b}（{gb}）"


# ---------------------------------------------------------------- 主検定（EXP-004 と同じ実装）
def design(lp, ll):
    lp, ll = np.asarray(lp, float), np.asarray(ll, float)
    return np.column_stack([np.ones_like(lp), lp, ll, lp ** 2, ll ** 2, lp * ll])


def eligibility(cs: pd.DataFrame, nulls: dict, tol: float, cover: float) -> dict[str, list[str]]:
    out = {}
    for gene, tg in cs.groupby("gene"):
        lp, ll = tg.lp.to_numpy(), tg.ll.to_numpy()
        elig = []
        for h, hd in nulls.items():
            if h == gene:
                continue
            D = np.maximum(np.abs(lp[:, None] - hd.lp.to_numpy()[None, :]),
                           np.abs(ll[:, None] - hd.ll.to_numpy()[None, :]))
            if (D <= tol).any(axis=1).mean() >= cover:
                elig.append(h)
        out[gene] = elig
    return out


def gene_stats(cs: pd.DataFrame, nulls: dict, elig: dict, col: str, min_elig: int) -> pd.DataFrame:
    rows = []
    Xn = {h: design(hd.lp, hd.ll) for h, hd in nulls.items()}
    yn = {h: hd[col].to_numpy() for h, hd in nulls.items()}
    for gene, tg in cs.groupby("gene"):
        others = [h for h in nulls if h != gene]
        X = np.vstack([Xn[h] for h in others])
        y = np.concatenate([yn[h] for h in others])
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        Dg = float((tg[col].to_numpy() - design(tg.lp, tg.ll) @ beta).mean())
        Dh = np.array([float((yn[h] - Xn[h] @ beta).mean()) for h in elig[gene]])
        q = float((Dh < Dg).mean()) if len(Dh) >= min_elig else np.nan
        rows.append({"gene": gene, "col": col, "D_raw": Dg, "q": q, "n_elig": len(Dh),
                     "n_cells": len(tg), "lp_med": float(tg.lp.median()),
                     "ll_med": float(tg.ll.median()), "relpos_med": float(tg.relpos.median())})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- 読み込み
def load(cfg: dict) -> pd.DataFrame:
    import pickle
    res = E / cfg["output"]["results_dir"]
    c = pd.read_csv(res / "cells.csv", dtype={"donor_id": str})
    n0 = len(c)
    c = c[~c.gene.isin(BLOCKED)].reset_index(drop=True)
    if n0 != len(c):
        raise SystemExit(f"回さないはずの遺伝子が {n0 - len(c)} 行あった。実行設定を確認すること")
    assert not set(c.gene) & set(BLOCKED), "回さないはずの遺伝子が結果に入っている"
    c["seed"] = c.run.astype(str).str.extract(r"run-(\d+)", expand=False).astype(int)

    from datasets import load_from_disk
    from geneformer import TOKEN_DICTIONARY_FILE
    tokd = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
    ds = load_from_disk(str(E / "work/tokenized.dataset"))
    uids = sorted(c.cell_uid.unique())
    ids = dict(zip(uids, ds.select(uids)["input_ids"]))
    c["pos"] = [list(ids[u]).index(tokd[g]) for u, g in zip(c.cell_uid, c.gene)]
    c["lp"] = np.log2(1 + c.pos.astype(float))
    c["ll"] = np.log2(c.length.astype(float))
    c["relpos"] = c.pos / c.length

    names = ["goal", "start", "alt", "myof"]
    vals = {n: np.full(len(c), np.nan) for n in names}
    dnorm = np.full(len(c), np.nan)
    s_ = cfg["states"]
    for (run, gene), part in c.groupby(["run", "gene"], sort=False):
        seed = int(str(run).split("-")[1])
        sz = np.load(res / run / "state_emb_vectors.npz")
        xz = np.load(E / f"work/extra_centroids-seed{seed}.npz")
        cen = {"goal": sz[s_["goal_state"]], "start": sz[s_["start_state"]],
               "alt": sz[s_["alt_states"][0]], "myof": xz["MYOF_neg"]}
        gd = res / run / "genes" / gene / "emb"
        o = np.load(gd / "orig.npy").astype(np.float64)
        p = np.load(gd / "pert.npy").astype(np.float64)
        assert len(o) == len(part), f"{gd}: 埋め込み {len(o)} ≠ cells.csv {len(part)}"
        idx = part.index.to_numpy()
        for n in names:
            vals[n][idx] = cos(p, cen[n].astype(np.float64)) - cos(o, cen[n].astype(np.float64))
        dnorm[idx] = np.linalg.norm(p - o, axis=1)
    for n in names:
        c[f"s_{n}"] = vals[n]
    c["dnorm"] = dnorm
    return c


# ---------------------------------------------------------------- セット単位の検定
def set_null(z: np.ndarray, k: int, n_draw: int, rng: np.random.Generator) -> np.ndarray:
    """対照遺伝子（N × seed）から k 個を復元なしで取ったセットの z 平均（n_draw × seed）。"""
    idx = np.argsort(rng.random((n_draw, z.shape[0])), axis=1)[:, :k]
    return z[idx].mean(axis=1)


def sham_sets(nulls: pd.DataFrame, target: pd.DataFrame, k: int, n_sets: int,
              rng: np.random.Generator) -> list[list[str]]:
    """偽マーカーセット: 目標に偏っていない対照遺伝子から、元の位置・入力の長さ・細胞数を
    目標マーカーセットにそろえて k 個を選ぶ。合格したら『合格は偏りでなく位置・長さで説明される』ことになる。"""
    out, used = [], set()
    for _ in range(n_sets):
        pick = []
        for _, t in target.iterrows():
            cand = nulls[~nulls.symbol.isin(used | set(pick))]
            if cand.empty:
                break
            dist = (np.abs(cand.lp_med - t.lp_med) / 0.25
                    + np.abs(cand.ll_med - t.ll_med) / 0.25
                    + np.abs(np.log(cand.n_cells / max(t.n_cells, 1))) / np.log(1.5))
            order = cand.symbol.to_numpy()[np.argsort(dist.to_numpy(), kind="stable")]
            for s in order[:5]:
                if s not in pick:
                    pick.append(s)
                    break
        used |= set(pick)
        out.append(pick[:k])
    return out


# ---------------------------------------------------------------- 判定
def main() -> None:
    cfg = yaml.safe_load(open(E / "config.yaml"))
    a = cfg["analysis"]
    cr = a["criteria"]
    gj = json.loads((E / "work/genes_001.json").read_text())["genes"]
    roles = {g["symbol"]: g["role"] for g in gj}
    ens = {g["symbol"]: g["ensembl"] for g in gj}
    half = {g["symbol"]: g.get("half") for g in gj}
    sym = {v: k for k, v in ens.items()}
    nulls_set = {ens[s] for s, r in roles.items() if r == "negative"}
    sets_def = {"goal_markers": [ens[s] for s, r in roles.items() if r == "marker_goal"],
                "myof_markers": [ens[s] for s, r in roles.items() if r == "marker_myof"],
                "housekeeping": [ens[s] for s, r in roles.items() if r == "housekeeping"],
                "negctrl_panel": [ens[s] for s, r in roles.items() if r == "negctrl_panel"]}
    assert cr["require_all_seeds"] is True, "config の require_all_seeds が False（plan.md §6 と食い違う）"
    assert len(sets_def["goal_markers"]) == cr["K_set"], "目標マーカーセットの大きさが config と違う"

    c = load(cfg)
    cols = metric_columns(c)
    seeds = sorted(c.seed.unique().tolist())
    assert seeds == sorted(cfg["robustness"]["seeds"]), "seed が計画と違う"

    pm = a["pos_matched_permutation"]
    frames = {}
    for s in seeds:
        cs = c[c.seed == s]
        nl = {h: hd for h, hd in cs[cs.gene.isin(nulls_set)].groupby("gene")}
        el = eligibility(cs, nl, pm["tolerance_log2"], pm["min_cover_frac"])
        for name, col in cols.items():
            f = gene_stats(cs, nl, el, col, pm["min_eligible_null_genes"])
            f["seed"] = s
            frames.setdefault(name, []).append(f)
    tabs = {k: pd.concat(v, ignore_index=True) for k, v in frames.items()}
    long = pd.concat(tabs.values(), ignore_index=True)
    long["symbol"] = long.gene.map(sym)
    long["role"] = long.symbol.map(roles)
    long["half"] = long.symbol.map(half)
    (E / "results/eval").mkdir(parents=True, exist_ok=True)
    long.to_csv(E / "results/eval/gene_seed_001.csv", index=False)

    # --- 軸ごとに z（対照遺伝子群で標準化）
    Zt: dict[str, pd.DataFrame] = {}
    for axis in ("goal", "myofaxis", "placebo"):
        T = tabs[axis].copy()
        for s in seeds:
            m = T.seed == s
            nul = T[m & T.gene.isin(nulls_set)]
            T.loc[m, "z"] = (T.loc[m, "D_raw"] - nul.D_raw.mean()) / nul.D_raw.std(ddof=1)
        T["symbol"] = T.gene.map(sym)
        T["role"] = T.symbol.map(roles)
        Zt[axis] = T

    def zmat(axis: str, genes: list[str]) -> np.ndarray:
        T = Zt[axis]
        return np.column_stack([[float(T[(T.gene == g) & (T.seed == s)].z.iloc[0]) for g in genes]
                                for s in seeds])

    rng = np.random.default_rng(cr["null_draw_seed"])
    NUL = sorted(nulls_set)
    znull = {ax: zmat(ax, NUL) for ax in ("goal", "myofaxis", "placebo")}

    M: dict = {"seeds": seeds, "n_genes": len(gj), "n_nulls": len(NUL)}

    # --- M1〜M4（指標の健全性。EXP-20260927-004 と同一の閾値）
    P = tabs["goal"]
    NULP = P[P.gene.isin(nulls_set)]
    M["M1_calibration"] = {}
    for s in seeds:
        q = NULP[NULP.seed == s].q.dropna()
        M["M1_calibration"][s] = {"n_valid": int(len(q)), "frac_hi": float((q >= 0.90).mean()),
                                  "frac_lo": float((q < 0.10).mean())}
    M["M1_pass"] = all(v["n_valid"] >= cr["M1_min_n"] and v["frac_hi"] <= cr["M1_tail_frac_max"]
                       and v["frac_lo"] <= cr["M1_tail_frac_max"] for v in M["M1_calibration"].values())

    M["M2_confounding"] = {}
    for s in seeds:
        n = NULP[NULP.seed == s]
        M["M2_confounding"][s] = {k: float(stats.spearmanr(n[f"{k}_med"], n.D_raw).statistic)
                                  for k in ("lp", "ll", "relpos")}
    M["M2_pass"] = sum(all(abs(v) <= cr["M2_abs_rho_max"] for v in dd.values())
                       for dd in M["M2_confounding"].values()) >= cr["M2_min_seeds"]

    qv = lambda tab, g, s: (float(tab[(tab.gene == g) & (tab.seed == s)].q.iloc[0])
                            if len(tab[(tab.gene == g) & (tab.seed == s)]) else np.nan)
    lo, hi = cr["M3_range"]
    inr = lambda qs: sum(1 for x in qs if not np.isnan(x) and lo <= x <= hi) >= cr["M3_min_seeds"]
    M["M3_housekeeping"] = {sym[g]: [qv(P, g, s) for s in seeds] for g in sets_def["housekeeping"]}
    M["M3_panel"] = {sym[g]: [qv(P, g, s) for s in seeds] for g in sets_def["negctrl_panel"]}
    M["M3_pass"] = (all(inr(v) for v in M["M3_housekeeping"].values())
                    and sum(inr(v) for v in M["M3_panel"].values()) >= cr["M3_panel_min_genes"])

    M["M4_rho"] = {}
    for s in seeds:
        n = NULP[NULP.seed == s].set_index("gene")
        u = tabs["uniform"][tabs["uniform"].seed == s].set_index("gene").reindex(n.index)
        pl = tabs["placebo"][tabs["placebo"].seed == s].set_index("gene").reindex(n.index)
        M["M4_rho"][s] = {"uniform": float(stats.spearmanr(n.D_raw, u.D_raw).statistic),
                          "placebo": float(stats.spearmanr(n.D_raw, pl.D_raw).statistic)}
    M["M4_pass"] = all(abs(v["uniform"]) <= cr["M4_abs_rho_uniform_max"]
                       and abs(v["placebo"]) <= cr["M4_abs_rho_placebo_max"] for v in M["M4_rho"].values())
    M["metric_pass"] = bool(M["M1_pass"] and M["M2_pass"] and M["M3_pass"] and M["M4_pass"])

    # --- セット単位の検定
    compare_same_scale(cols["goal"], cols["myofaxis"])
    compare_same_scale(cols["goal"], cols["placebo"])

    # 無作為な対照セットの分布は (軸, K) ごとに 1 回だけ引き、すべての検定で同じものを使う
    # （偽マーカーセットや偶然合格率の実測が、目標マーカーセットと違う閾値で判定されないようにする）
    _nullcache: dict[tuple[str, int], np.ndarray] = {}

    def null_dist(axis: str, k: int) -> np.ndarray:
        key = (axis, k)
        if key not in _nullcache:
            _nullcache[key] = set_null(znull[axis], k, cr["n_permutation"], rng)
        return _nullcache[key]

    def set_test(axis: str, genes: list[str], label: str) -> dict:
        k = len(genes)
        obs = zmat(axis, genes).mean(axis=0)                       # seed ごとのセット平均 z
        nd = null_dist(axis, k)                                    # n_draw × seed
        pct = [float((nd[:, i] < obs[i]).mean()) for i in range(len(seeds))]
        thr = [float(np.quantile(nd[:, i], cr["set_percentile"])) for i in range(len(seeds))]
        passed = bool(all(o > t for o, t in zip(obs, thr))
                      and all(o >= cr["set_abs_min"] for o in obs))
        return {"label": label, "axis": axis, "k": k,
                "set_mean_z": [round(float(x), 3) for x in obs],
                "null_percentile": [round(x, 5) for x in pct],
                "null_threshold": [round(x, 3) for x in thr],
                "pass": passed}

    res = {}
    for name, genes in sets_def.items():
        for axis in ("goal", "myofaxis"):
            res[f"{name}@{axis}"] = set_test(axis, genes, name)
    M["set_tests"] = res

    # G-1 主判定 / G-3 アッセイの妥当性
    M["G1_goal_markers_on_goal_axis"] = res["goal_markers@goal"]["pass"]
    M["G3_myof_markers_on_myof_axis"] = res["myof_markers@myofaxis"]["pass"]

    # G-2 交差（特異性）
    gg = np.array(res["goal_markers@goal"]["set_mean_z"])
    gm = np.array(res["goal_markers@myofaxis"]["set_mean_z"])
    mm = np.array(res["myof_markers@myofaxis"]["set_mean_z"])
    mg = np.array(res["myof_markers@goal"]["set_mean_z"])
    M["G2_crossover"] = {"goal_set: 目標の軸 − 筋線維芽細胞の軸": [round(float(x), 3) for x in gg - gm],
                         "myof_set: 筋線維芽細胞の軸 − 目標の軸": [round(float(x), 3) for x in mm - mg],
                         "pass": bool(all(gg - gm > 0) and all(mm - mg > 0))}

    # G-4 中立性: 恒常発現遺伝子・中立対照パネルのセット平均が無作為セットの 5〜95 パーセンタイル内
    neu = {}
    for name in ("housekeeping", "negctrl_panel"):
        k = len(sets_def[name])
        obs = zmat("goal", sets_def[name]).mean(axis=0)
        nd = null_dist("goal", k)
        pc = [float((nd[:, i] < obs[i]).mean()) for i in range(len(seeds))]
        neu[name] = {"set_mean_z": [round(float(x), 3) for x in obs],
                     "null_percentile": [round(x, 4) for x in pc],
                     "n_seeds_in_range": int(sum(cr["neutral_range"][0] <= x <= cr["neutral_range"][1]
                                                 for x in pc))}
    M["G4_neutral"] = neu
    M["G4_pass"] = all(v["n_seeds_in_range"] >= cr["neutral_min_seeds"] for v in neu.values())

    # G-5 偶然合格率の実測（対照遺伝子から無作為に取った K 個セットに、G-1 と同じ規則を当てはめる）
    k = cr["K_set"]
    thr = [float(np.quantile(null_dist("goal", k)[:, i], cr["set_percentile"])) for i in range(len(seeds))]
    nd = set_null(znull["goal"], k, cr["n_chance_draw"], rng)   # 閾値を決めたのとは別に引き直す
    hit = np.ones(len(nd), bool)
    for i in range(len(seeds)):
        hit &= (nd[:, i] > thr[i]) & (nd[:, i] >= cr["set_abs_min"])
    M["G5_chance_pass_rate_random_sets"] = {"n_draw": int(len(nd)), "rate": float(hit.mean()),
                                            "threshold": [round(x, 3) for x in thr]}

    # G-6 偽マーカーセット（位置・長さ・細胞数をそろえた対照遺伝子のセット）
    tgt = P[(P.seed == seeds[0]) & (P.gene.isin(sets_def["goal_markers"]))][
        ["gene", "lp_med", "ll_med", "n_cells"]].copy()
    tgt["symbol"] = tgt.gene.map(sym)
    npool = P[(P.seed == seeds[0]) & (P.gene.isin(nulls_set))][
        ["gene", "lp_med", "ll_med", "n_cells"]].copy()
    npool["symbol"] = npool.gene.map(sym)
    shams = sham_sets(npool, tgt, k, cr["n_sham_sets"], np.random.default_rng(cr["null_draw_seed"] + 1))
    sh = []
    for i, names in enumerate(shams):
        genes = [ens[s] for s in names]
        r = set_test("goal", genes, f"sham_{i + 1}")
        r["genes"] = names
        sh.append(r)
    M["G6_sham_sets"] = sh
    M["G6_pass"] = bool(not any(r["pass"] for r in sh))
    pd.DataFrame([{"set": r["label"], "z": r["set_mean_z"], "pass": r["pass"]} for r in sh]).to_csv(
        E / "results/eval/set_null_001.csv", index=False)

    # --- 記述値（合否に入れない）
    zg = Zt["goal"]
    M["descriptive"] = {
        "単一遺伝子で 3 seed すべて 2.5 標準偏差以上のマーカー": {
            name: sorted(sym[g] for g in genes if all(
                float(zg[(zg.gene == g) & (zg.seed == s)].z.iloc[0]) >= 2.5 for s in seeds))
            for name, genes in sets_def.items()},
        "単一遺伝子の偶然合格率（対照 96 に同じ規則）": float(np.mean(
            (znull["goal"] >= 2.5).all(axis=1))),
        "目標マーカーセットの半分ごとのセット平均 z（上位 / 無作為）": {
            h: [round(float(x), 3) for x in zmat("goal", [ens[s] for s in roles
                                                          if roles[s] == "marker_goal" and half[s] == h]
                                                 ).mean(axis=0)] for h in ("top", "random")},
        "参照 ACTA2 / TAGLN（筋線維芽細胞の軸、標準偏差単位）": {
            s: [round(float(Zt["myofaxis"][(Zt["myofaxis"].symbol == s) & (Zt["myofaxis"].seed == sd)].z.iloc[0]), 2)
                for sd in seeds] for s in ("ACTA2", "TAGLN")},
        "偏り（bias）と z の順位相関（目標マーカー 24、目標の軸）": {
            str(s): float(stats.spearmanr(
                [next(g["bias_goal"] for g in gj if g["ensembl"] == x) for x in sets_def["goal_markers"]],
                [float(zg[(zg.gene == x) & (zg.seed == s)].z.iloc[0]) for x in sets_def["goal_markers"]]
            ).statistic) for s in seeds},
    }

    M["positive_control_established"] = bool(
        M["metric_pass"] and M["G1_goal_markers_on_goal_axis"] and M["G2_crossover"]["pass"]
        and M["G3_myof_markers_on_myof_axis"] and M["G4_pass"] and M["G6_pass"]
        and M["G5_chance_pass_rate_random_sets"]["rate"] <= cr["max_chance_rate"])
    M["conclusion_proposal"] = (
        "目標の軸での陽性対照が確立した（セット単位）" if M["positive_control_established"] else
        "アッセイが働いていない（筋線維芽細胞マーカーが自分の軸で合格しない）→ 判定不能"
        if not M["G3_myof_markers_on_myof_axis"] else
        "指標の合格条件を満たさない → 判定不能" if not M["metric_pass"] else
        "目標の軸では、マーカーセット単位でも対照と区別できない")
    M["note"] = ("判定の確定は PI（G3）。本スクリプトは事前登録の規則を当てはめて提案するだけ。"
                 "本実験は主対象 LTBR/RELB/NFKB2 を回していないので、HYP-20260925-02 の判定には使えない。"
                 "in silico の結果であり、実験的に確認された生物学的事実ではない。")

    (E / "results/eval/criteria_001.json").write_text(
        json.dumps(M, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps({k: v for k, v in M.items() if k.endswith("_pass") or k.startswith("G")
                      or k in ("positive_control_established", "conclusion_proposal")},
                     ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    main()

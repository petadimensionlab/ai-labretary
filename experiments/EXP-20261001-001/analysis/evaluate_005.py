#!/usr/bin/env python3
"""EXP-20261001-001 — 事前登録した判定の機械的な当てはめ（plan.md §6 をそのまま実装する）。

EXP-20260930-001 からの変更点（**すべて PI が 2026-10-01 に決め、結果を見る前に固定した**）:

A. **主指標を「向きへの射影」にする。**
   û = 正規化( 目標の重心/‖目標の重心‖ − 出発点の重心/‖出発点の重心‖ )、
   主指標 = (摂動後 − 摂動前) · û / ‖摂動後 − 摂動前‖。
   **2 つの重心への寄りの差 (s_goal − s_start) は使わない**（差し引く相手で符号が変わり、
   一様成分と算術的に相関するため。plan.md §6.1）。

A-2. **交絡を除く回帰に移動量 `ld` を入れない。** 説明変数は `lp`（元の順位）と
   `ll`（入力の長さ）の完全二次式だけ（`movement_correction.gene_stats(use_dnorm=False)`）。

A-3. **一様成分は、主指標の射影に使っていない状態から作る。**
   退化した線維芽細胞 `dFIB` と退化した髄質線維芽細胞 `dMFIB` への寄りの平均。
   主指標と一様成分を同じ 2 量の差と和で作ると算術的に相関する（理論値 0.28〜0.30 が
   実測と小数第 4 位まで一致した）。**閾値 0.40 は動かさない。変えるのは一様成分の定義である。**

B. **中立性 P-5 は、採用規則 H-1〜H-3 を通った中立対照だけで判定する**（plan.md §6.1.2）。
   恒常発現遺伝子は**記述にとどめる**（単一核では「どの細胞でも一定」の前提が成り立たない）。

**途中のチェックポイントでは主判定のセットの値を計算も保存もしない**（`--checkpoint`。plan.md §6.4）。

判定対象にしないもの: 盲検中の TP53 / SPP1 / LCN2（読み込み直後に除外し、1 行でも残れば停止）。
LTBR・CDKN1A は PI が 2026-09-27 に解除済みなので摂動も解析も行う。

出力: results/eval/criteria_005.json, results/eval/gene_seed_005.csv, results/eval/set_sham_005.csv
実行: .venv/bin/python experiments/EXP-20261001-001/analysis/evaluate_005.py [--checkpoint]
      試し実行の出力で評価コードを通すとき:
      .venv/bin/python .../evaluate_005.py --results-dir results-smoke --smoke
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
E = ROOT / "experiments/EXP-20261001-001"
sys.path.insert(0, str(ROOT / "scripts/analysis"))
from movement_correction import gene_stats  # noqa: E402

# 盲検中（摂動も解析もしない）。LTBR・CDKN1A は解除済みなので入れない。
BLOCKED = {"ENSG00000141510": "TP53", "ENSG00000118785": "SPP1", "ENSG00000148346": "LCN2"}

# チェックポイントで見てよいセット（plan.md §6.4）。**主判定とその内訳はここに入れない。**
CHECKPOINT_SETS = ("positive_control71", "positive_control65", "positive_control39",
                   "positive_control70", "positive_control38", "positive_control35",
                   "neutral_adopted", "neutral_adopted_cv055", "housekeeping", "negctrl_panel")

# 主判定とその入れ子（**チェックポイントで計算してはいけない**）
PRIMARY_FAMILY = ("INT", "INT_S1", "INT_S2", "INT_S3", "INT_pathway_specific", "INT_plus_T1T2",
                  "T1", "T2")


def cos(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return (a @ b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b) + 1e-12)


def eligibility(cs: pd.DataFrame, nulls: dict, tol: float, cover: float) -> dict[str, list[str]]:
    """遺伝子ごとに「元の位置と入力の長さで比べられる対照遺伝子」を数える。
    EXP-20260928-001 以降と同一の実装（**変えていない**）。"""
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


def load(cfg: dict, results_dir: str | None = None) -> tuple[pd.DataFrame, dict]:
    """細胞ごとの表を作る。

    **主指標は「向きへの射影」**（plan.md §6.1 の A）。
    **摂動が起きなかった細胞（移動量 0）は外す**（過剰発現で入力が 1 ビットも変わらなかった細胞）。
    """
    import pickle
    from datasets import load_from_disk
    from geneformer import TOKEN_DICTIONARY_FILE

    res = E / (results_dir or cfg["output"]["results_dir"])
    c = pd.read_csv(res / "cells.csv", dtype={"donor_id": str})
    n0 = len(c)
    c = c[~c.gene.isin(BLOCKED)].reset_index(drop=True)
    if n0 != len(c):
        raise SystemExit(f"回さないはずの遺伝子が {n0 - len(c)} 行あった。実行設定を確認すること")
    c["seed"] = c.run.astype(str).str.extract(r"run-(\d+)", expand=False).astype(int)

    tokd = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
    ds = load_from_disk(str(E / "work/tokenized.dataset"))
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
    for (run, gene), part in c.groupby(["run", "gene"], sort=False):
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
        # **主指標の向き û**（単位化した重心の差。コサインで測る「寄り」と同じ尺度にそろえる）
        vg = sz[names["goal"]].astype(np.float64)
        vs = sz[names["start"]].astype(np.float64)
        u = vg / np.linalg.norm(vg) - vs / np.linalg.norm(vs)
        u = u / np.linalg.norm(u)
        with np.errstate(divide="ignore", invalid="ignore"):
            proj[idx] = (d @ u) / dn
    for n in names:
        c[f"s_{n}"] = vals[n]
    c["dnorm"] = dnorm
    c["m_proj"] = proj

    # **移動量 0 の細胞を外す**（guardrails 2026-09-30。比べる版すべてで同じ細胞を外す）
    zero = c.dnorm <= 0
    drop = {"n_cells_before": int(len(c)), "n_zero_movement": int(zero.sum()),
            "frac_zero_movement": round(float(zero.mean()), 6),
            "genes_with_zero": sorted(c.loc[zero, "gene"].unique().tolist())[:20]}
    c = c[~zero].reset_index(drop=True)
    drop["n_cells_after"] = int(len(c))
    return c, drop


def metric_columns(c: pd.DataFrame, cfg: dict) -> dict[str, str]:
    """**主指標は向きへの射影 1 本だけ**（plan.md §6.1 A）。
    一様成分は**射影に使っていない状態**（dFIB・dMFIB）の寄りの平均（同 A-3）。"""
    uni = cfg["metric"]["uniform_component_states"]
    for a in uni:
        if f"s_{a}" not in c:
            raise SystemExit(f"一様成分に使う状態 {a} の寄りが無い（config の alt_states を確認）")
    c["m_goal"] = c.m_proj                                        # **主指標**
    c["m_uniform"] = c[[f"s_{a}" for a in uni]].mean(axis=1)      # **一様成分（健全性 M4）**
    cols = {"goal": "m_goal", "uniform": "m_uniform"}
    # ---- ここから下は **記述のみ。判定に一切使わない** ----
    c["m_desc_diff"] = (c.s_goal - c.s_start) / c.dnorm           # 前回（EXP-20260930-001）の主指標
    c["m_desc_uniform_old"] = (c.s_goal + c.s_start) / 2.0        # 前回の一様成分（算術的に相関する）
    cols |= {"desc_diff": "m_desc_diff", "desc_uniform_old": "m_desc_uniform_old"}
    for a in cfg["states"].get("alt_states", []):
        c[f"m_desc_{a}"] = (c[f"s_{a}"] - c.s_start) / c.dnorm
        cols[f"desc_{a}"] = f"m_desc_{a}"
    return cols


def set_null(z: np.ndarray, k: int, n_draw: int, rng, replace: bool = True) -> np.ndarray:
    N = z.shape[0]
    if replace:
        idx = rng.integers(0, N, size=(n_draw, k))
    else:
        if k > N:
            raise ValueError(f"復元なしでは対照 {N} 個から {k} 個を引けない")
        idx = np.argsort(rng.random((n_draw, N)), axis=1)[:, :k]
    return z[idx].mean(axis=1)


def sham_sets(nulls: pd.DataFrame, target: pd.DataFrame, k: int, n_sets: int, rng,
              match_bias: bool = True) -> list[list[str]]:
    """偽のセット: 対照遺伝子から、元の位置・入力の長さ・細胞数（と発現の偏り）を
    対象セットにそろえて k 個を選ぶ。組の間で遺伝子の重複は許し、組の中では許さない。"""
    out = []
    for _ in range(n_sets):
        pick: list[str] = []
        order_t = target.sample(frac=1.0, random_state=int(rng.integers(0, 2 ** 31 - 1)))
        for _, t in order_t.iterrows():
            cand = nulls[~nulls.symbol.isin(pick)]
            if cand.empty:
                break
            dist = (np.abs(cand.lp_med - t.lp_med) / 0.25
                    + np.abs(cand.ll_med - t.ll_med) / 0.25
                    + np.abs(np.log(cand.n_cells / max(t.n_cells, 1))) / np.log(1.5))
            if match_bias:
                dist = dist + np.abs(cand.bias_goal - t.bias_goal) / 0.25
            order = cand.symbol.to_numpy()[np.argsort(dist.to_numpy(), kind="stable")]
            take = order[:5]
            pick.append(str(take[int(rng.integers(0, len(take)))]))
        out.append(pick[:k])
    return out


def main(checkpoint: bool = False, results_dir: str | None = None, smoke: bool = False) -> None:
    cfg = yaml.safe_load(open(E / "config.yaml"))
    a = cfg["analysis"]
    cr = a["criteria"]
    mc = cfg["metric"]["confound_correction"]
    gj = json.loads((E / "work/genes_005.json").read_text())["genes"]
    roles = {g["symbol"]: g["role"] for g in gj}
    ens = {g["symbol"]: g["ensembl"] for g in gj}
    bias = {g["symbol"]: g["bias_goal"] for g in gj}
    sym = {v: k for k, v in ens.items()}

    npool = json.loads((E / "work/null_pool_005.json").read_text())["versions"]
    null_versions = {k: [ens[x["symbol"]] for x in v] for k, v in npool.items()}
    subs = json.loads((E / "work/subsets_005.json").read_text())["subsets"]
    sets_def = {k: [ens[x["symbol"]] for x in v] for k, v in subs.items()}
    PRIMARY, PRIMARY_NULL = a["primary_set"], a["primary_null_version"]
    PC, NEU = a["primary_positive_control"], a["primary_neutral_set"]
    ABS = a["set_abs_min"]                      # セットごとの絶対値条件（**結果を見る前に固定**）

    # ---- 事前登録どおりの設定で動いているかの検査（違えば止める） ----
    assert cr["require_all_seeds"] is True, "config の require_all_seeds が False"
    assert len(sets_def[PRIMARY]) == cr["K_set"], "主判定のセットの大きさが config と違う"
    assert PRIMARY == "INT" and PC == "positive_control71" and NEU == "neutral_adopted", \
        "主判定・陽性対照・中立性のセットが事前登録と違う"
    assert set(null_versions["A_300"]) < set(null_versions["B_500"]), "対照の部分集合が入れ子でない"
    for name, genes in sets_def.items():
        assert not (set(null_versions["B_500"]) & set(genes)), f"対照と {name} が重複している"
    # **A-2: 移動量 ld を入れない**
    assert mc["variables"] == ["lp", "ll"], \
        f"主判定の補正の版が事前登録と違う（{mc['variables']}。ld を入れてはいけない）"
    assert mc["drop_zero_movement"] is True, "移動量 0 の細胞を外す設定になっていない"
    # **A: 主指標は射影**
    assert cfg["metric"]["primary_form"].startswith("projection_"), "主指標が射影になっていない"
    # **A-3: 一様成分は射影に使っていない状態から作る**
    uni = cfg["metric"]["uniform_component_states"]
    proj_states = {cfg["metric"]["projection_states"]["from"], cfg["metric"]["projection_states"]["to"]}
    assert not ({"start", "goal"} & set(uni)) and proj_states == {"start", "goal"}, \
        "一様成分が主指標の射影に使った状態から作られている（算術的に相関する）"

    c, drop = load(cfg, results_dir)
    cols = metric_columns(c, cfg)
    seeds = sorted(c.seed.unique().tolist())
    if not (checkpoint or smoke):
        assert seeds == sorted(cfg["robustness"]["seeds"]), "seed が計画と違う"

    pm = a["pos_matched_permutation"]
    min_elig = 1 if smoke else pm["min_eligible_null_genes"]
    nulls_all = set(null_versions["B_500"]) & set(c.gene.unique())
    frames: dict[str, list] = {}
    frames_ld: dict[str, list] = {}
    for s in seeds:
        cs = c[c.seed == s]
        nl = {h: hd for h, hd in cs[cs.gene.isin(nulls_all)].groupby("gene")}
        el = eligibility(cs, nl, pm["tolerance_log2"], pm["min_cover_frac"])
        for name, col in cols.items():
            # **A-2: use_dnorm=False（lp・ll の完全二次式のみ）**
            f = gene_stats(cs, nl, el, col, min_elig, use_dnorm=False)
            f["seed"] = s
            frames.setdefault(name, []).append(f)
            if name == "goal":   # **記述のみ**: 移動量を入れた版（版の比較。判定に使わない）
                fo = gene_stats(cs, nl, el, col, min_elig, use_dnorm=True)
                fo["seed"] = s
                frames_ld.setdefault(name, []).append(fo)
    tabs = {k: pd.concat(v, ignore_index=True) for k, v in frames.items()}
    tab_ld = pd.concat(frames_ld["goal"], ignore_index=True)
    long = pd.concat(tabs.values(), ignore_index=True)
    long["symbol"] = long.gene.map(sym)
    long["role"] = long.symbol.map(roles)
    outdir = E / (results_dir or cfg["output"]["results_dir"]) / "eval"
    outdir.mkdir(parents=True, exist_ok=True)
    long.to_csv(outdir / "gene_seed_005.csv", index=False)

    def standardise(T: pd.DataFrame, nset: set) -> pd.DataFrame:
        T = T.copy()
        for s in seeds:
            m = T.seed == s
            nul = T[m & T.gene.isin(nset)]
            T.loc[m, "z"] = (T.loc[m, "D_raw"] - nul.D_raw.mean()) / nul.D_raw.std(ddof=1)
        T["symbol"] = T.gene.map(sym)
        T["role"] = T.symbol.map(roles)
        return T

    nv = {k: [g for g in v if g in set(c.gene.unique())] for k, v in null_versions.items()}
    Zt = {(ver, ax): standardise(tabs[ax], set(gs)) for ver, gs in nv.items() for ax in cols}
    Z_ld = standardise(tab_ld, set(nv[PRIMARY_NULL]))

    present = set(c.gene.unique())
    sets_run = {k: [g for g in v if g in present] for k, v in sets_def.items()}

    def zmat(ver: str, axis: str, genes: list[str]) -> np.ndarray:
        T = Zt[(ver, axis)]
        return np.column_stack([[float(T[(T.gene == g) & (T.seed == s)].z.iloc[0]) for g in genes]
                                for s in seeds])

    rng = np.random.default_rng(cr["null_draw_seed"])
    znull = {(v, ax): zmat(v, ax, sorted(nv[v])) for v in nv for ax in cols}
    _cache: dict[tuple, np.ndarray] = {}
    n_perm = 2000 if smoke else cr["n_permutation"]

    def null_dist(ver: str, axis: str, k: int, replace: bool) -> np.ndarray:
        key = (ver, axis, k, replace)
        if key not in _cache:
            _cache[key] = set_null(znull[(ver, axis)], k, n_perm, rng, replace)
        return _cache[key]

    def set_test(ver: str, axis: str, name: str, genes: list[str], abs_min: float) -> dict:
        k = len(genes)
        obs = zmat(ver, axis, genes).mean(axis=0)
        out = {"label": name, "axis": axis, "null_version": ver, "k": k, "abs_min": abs_min,
               "set_mean_z": [round(float(x), 3) for x in obs]}
        for rep in (True, False):
            if not rep and k > znull[(ver, axis)].shape[0]:
                continue
            nd = null_dist(ver, axis, k, rep)
            thr = [float(np.quantile(nd[:, i], cr["set_percentile"])) for i in range(len(seeds))]
            pct = [float((nd[:, i] < obs[i]).mean()) for i in range(len(seeds))]
            tag = "replace" if rep else "no_replace"
            out[f"null_threshold_{tag}"] = [round(x, 3) for x in thr]
            out[f"null_percentile_{tag}"] = [round(x, 5) for x in pct]
            out[f"pass_{tag}"] = bool(all(o > t for o, t in zip(obs, thr))
                                      and all(o >= abs_min for o in obs))
        out["pass"] = out["pass_replace"] if cr["set_draw_replace"] else out["pass_no_replace"]
        out["all_seeds_positive"] = bool(all(o > 0 for o in obs))
        return out

    def chance_rate(ver: str, axis: str, k: int, abs_min: float) -> dict:
        rep = cr["set_draw_replace"]
        thr = [float(np.quantile(null_dist(ver, axis, k, rep)[:, i], cr["set_percentile"]))
               for i in range(len(seeds))]
        nd = set_null(znull[(ver, axis)], k, 2000 if smoke else cr["n_chance_draw"], rng, rep)
        hit = np.ones(len(nd), bool)
        for i in range(len(seeds)):
            hit &= (nd[:, i] > thr[i]) & (nd[:, i] >= abs_min)
        return {"n_draw": int(len(nd)), "rate": float(hit.mean()),
                "threshold": [round(x, 3) for x in thr]}

    M: dict = {"mode": "checkpoint" if checkpoint else ("smoke" if smoke else "final"),
               "seeds": seeds, "n_genes": len(gj), "zero_movement": drop,
               "metric_version": {"primary_form": cfg["metric"]["primary_form"],
                                  "confound_variables": mc["variables"],
                                  "uniform_component_states": uni,
                                  "drop_zero_movement": mc["drop_zero_movement"]},
               "n_nulls": {k: len(v) for k, v in nv.items()},
               "primary_set": PRIMARY, "primary_null_version": PRIMARY_NULL,
               "primary_positive_control": PC, "primary_neutral_set": NEU}
    if smoke:
        M["smoke_note"] = ("**試し実行の出力に対する評価であり、判定ではない。** "
                           "遺伝子も細胞も本番より大幅に少なく、閾値にも意味は無い。"
                           "目的は『評価コードが最後まで通ること』の確認だけである（plan.md §10）。")

    # ---------------- 欠測（比べる相手が足りずに計算できなかった遺伝子）は「不合格」と数えない
    P = tabs["goal"]
    M["missing_q_by_seed"] = {int(k): int(v)
                              for k, v in P[P.q.isna()].groupby("seed").gene.count().to_dict().items()}

    # ---------------- 指標の健全性 M1〜M5
    NULP = P[P.gene.isin(set(nv[PRIMARY_NULL]))]
    M["M1_calibration"] = {}
    for s in seeds:
        q = NULP[NULP.seed == s].q.dropna()
        M["M1_calibration"][s] = {"n_valid": int(len(q)), "frac_hi": float((q >= 0.90).mean()),
                                  "frac_lo": float((q < 0.10).mean())}
    m1_min_n = 1 if smoke else cr["M1_min_n"]
    M["M1_pass"] = all(v["n_valid"] >= m1_min_n and v["frac_hi"] <= cr["M1_tail_frac_max"]
                       and v["frac_lo"] <= cr["M1_tail_frac_max"]
                       for v in M["M1_calibration"].values())
    M["M2_confounding"] = {s: {k: float(stats.spearmanr(NULP[NULP.seed == s][f"{k}_med"],
                                                        NULP[NULP.seed == s].D_raw).statistic)
                               for k in ("lp", "ll", "relpos")} for s in seeds}
    M["M2_pass"] = sum(all(abs(v) <= cr["M2_abs_rho_max"] for v in dd.values())
                       for dd in M["M2_confounding"].values()) >= cr["M2_min_seeds"]
    qv = lambda g, s: (float(P[(P.gene == g) & (P.seed == s)].q.iloc[0])
                       if len(P[(P.gene == g) & (P.seed == s)]) else np.nan)
    lo, hi = cr["M3_range"]
    inr = lambda qs: sum(1 for x in qs if not np.isnan(x) and lo <= x <= hi) >= cr["M3_min_seeds"]
    # **M3 は採用した中立対照で判定する。恒常発現・パネル全体は記述のみ**（plan.md §6.1.2）
    M["M3_adopted"] = {sym[g]: [qv(g, s) for s in seeds] for g in sets_run[NEU]}
    M["M3_descriptive"] = {
        "恒常発現 3（判定に使わない）": {sym[g]: [qv(g, s) for s in seeds]
                                        for g in sets_run["housekeeping"]},
        "中立対照パネル 7（判定に使わない）": {sym[g]: [qv(g, s) for s in seeds]
                                              for g in sets_run["negctrl_panel"]}}
    n_in = sum(inr(v) for v in M["M3_adopted"].values())
    M["M3_n_in_range"] = int(n_in)
    M["M3_pass"] = bool(n_in >= min(cr["M3_adopted_min_genes"], len(M["M3_adopted"])))
    M["M4_rho"] = {}
    for s in seeds:
        n = NULP[NULP.seed == s].set_index("gene")
        u = tabs["uniform"][tabs["uniform"].seed == s].set_index("gene").reindex(n.index)
        od = tabs["desc_uniform_old"][tabs["desc_uniform_old"].seed == s].set_index("gene").reindex(n.index)
        M["M4_rho"][s] = {
            "uniform": float(stats.spearmanr(n.D_raw, u.D_raw).statistic),
            "（記述）前回の一様成分": float(stats.spearmanr(n.D_raw, od.D_raw).statistic)}
    M["M4_pass"] = all(abs(v["uniform"]) <= cr["M4_abs_rho_uniform_max"] for v in M["M4_rho"].values())
    M["M5_residual_vs_movement"] = {}
    for s in seeds:
        n = NULP[NULP.seed == s]
        dmed = c[c.seed == s].groupby("gene").dnorm.median()
        M["M5_residual_vs_movement"][s] = float(
            stats.spearmanr(dmed.reindex(n.gene).to_numpy(), n.D_raw.to_numpy()).statistic)
    M["M5_pass"] = all(abs(v) <= cr["M5_abs_rho_movement_max"]
                       for v in M["M5_residual_vs_movement"].values())
    M["metric_pass"] = bool(M["M1_pass"] and M["M2_pass"] and M["M3_pass"] and M["M4_pass"]
                            and M["M5_pass"])

    # ---------------- P-4 健全性チェック（陽性対照）と P-5 中立性
    tests: dict[str, dict] = {}
    for name in CHECKPOINT_SETS:
        if sets_run[name]:
            tests[name] = set_test(PRIMARY_NULL, "goal", name, sets_run[name], ABS.get(name, -1e9))
    M["P4_positive_control"] = tests.get(PC)
    M["P4_pass"] = bool(tests[PC]["pass"]) if PC in tests else None
    M["P4_sensitivity"] = {k: tests[k] for k in
                           ("positive_control65", "positive_control39", "positive_control70",
                            "positive_control38", "positive_control35") if k in tests}
    neu = {}
    for name in (NEU, "neutral_adopted_cv055", "housekeeping", "negctrl_panel"):
        if not sets_run[name]:
            continue
        k = len(sets_run[name])
        obs = zmat(PRIMARY_NULL, "goal", sets_run[name]).mean(axis=0)
        nd = null_dist(PRIMARY_NULL, "goal", k, cr["set_draw_replace"])
        pc = [float((nd[:, i] < obs[i]).mean()) for i in range(len(seeds))]
        neu[name] = {"set_mean_z": [round(float(x), 3) for x in obs],
                     "null_percentile": [round(x, 4) for x in pc],
                     "n_seeds_in_range": int(sum(cr["neutral_range"][0] <= x <= cr["neutral_range"][1]
                                                 for x in pc)),
                     "判定に使うか": name == NEU}
    M["P5_neutral"] = neu
    M["P5_pass"] = bool(neu[NEU]["n_seeds_in_range"] >= min(cr["neutral_min_seeds"], len(seeds))) \
        if NEU in neu else None

    if checkpoint:
        # **主判定のセットの値は計算も保存もしない**（plan.md §6.4）
        M["note"] = ("チェックポイント。**主判定（介入セット 119 個）とその入れ子の値は"
                     "計算も保存もしていない。** 見てよいのは装置が働いているかの確認だけ"
                     "（陽性対照・中立性・指標の健全性）。")
        (outdir / "checkpoint_005.json").write_text(
            json.dumps(M, ensure_ascii=False, indent=1, default=float) + "\n")
        print(json.dumps({k: v for k, v in M.items() if k.endswith("_pass")},
                         ensure_ascii=False, indent=1))
        return

    # ---------------- 主判定 P-1 と入れ子の感度分析
    for name in PRIMARY_FAMILY:
        if sets_run[name]:
            tests[name] = set_test(PRIMARY_NULL, "goal", name, sets_run[name], ABS.get(name, -1e9))
    M["set_tests"] = tests
    M["P1_primary"] = tests[PRIMARY]
    M["P1_pass"] = tests[PRIMARY]["pass"]
    M["R3_direction_reversed"] = not tests[PRIMARY]["all_seeds_positive"]

    # ---------------- 偶然合格率（対照で実測）と偽のセットの二項検定（plan.md §6.3）
    M["chance_pass_rate"] = {name: chance_rate(PRIMARY_NULL, "goal", len(sets_run[name]),
                                               ABS.get(name, -1e9))
                             for name in (PRIMARY, PC) if sets_run[name]}
    M["chance_rate_ok"] = bool(M["chance_pass_rate"][PRIMARY]["rate"] <= cr["max_chance_rate"])

    def tbl(genes: set) -> pd.DataFrame:
        t = P[(P.seed == seeds[0]) & (P.gene.isin(genes))][["gene", "lp_med", "ll_med", "n_cells"]].copy()
        t["symbol"] = t.gene.map(sym)
        t["bias_goal"] = t.symbol.map(bias)
        return t

    npl = tbl(set(nv[PRIMARY_NULL]))
    k_primary = len(sets_run[PRIMARY])
    n_sham = 5 if smoke else cr["n_sham_sets"]
    shams = sham_sets(npl, tbl(set(sets_run[PRIMARY])), k_primary, n_sham,
                      np.random.default_rng(cr["null_draw_seed"] + 1), cr["sham_match_bias"])
    rows = []
    for i, nm in enumerate(shams):
        r = set_test(PRIMARY_NULL, "goal", f"sham_{i + 1}", [ens[x] for x in nm], ABS[PRIMARY])
        r["genes"] = nm
        rows.append(r)
    kk = int(sum(r["pass"] for r in rows))
    N = len(rows)
    p = float(M["chance_pass_rate"][PRIMARY]["rate"])
    pval = float(stats.binomtest(kk, N, max(p, 1e-12), alternative="greater").pvalue)
    sm = np.array([r["set_mean_z"] for r in rows], dtype=float)
    rnd = null_dist(PRIMARY_NULL, "goal", k_primary, cr["set_draw_replace"])
    M["P_sham_sets"] = {
        "N": N, "p_単一の偶然合格率（対照で実測）": round(p, 6), "期待合格数": round(N * p, 3),
        "k_実際の合格数": kk, "片側 p 値": round(pval, 6),
        "pass": bool(pval >= cr["sham_binomial_alpha"]),
        "注": "「N 組すべて不合格」は条件にしない。問題のない設計でも落ちるため（plan.md §6.3）",
        "較正": {"偽セットの平均": [round(float(x), 4) for x in sm.mean(axis=0)],
                 "無作為セットの平均": [round(float(x), 4) for x in rnd.mean(axis=0)],
                 "散らばりの比": [round(float(x / y), 3) for x, y in
                                  zip(sm.std(axis=0, ddof=1), rnd.std(axis=0, ddof=1))]}}
    M["P_sham_sets"]["較正"]["中心のずれが許容内"] = bool(
        all(abs(x - y) <= cr["sham_centre_tol"] for x, y in zip(sm.mean(axis=0), rnd.mean(axis=0))))
    M["P_sham_sets"]["較正"]["散らばりが許容内"] = bool(
        all(cr["sham_spread_range"][0] <= r <= cr["sham_spread_range"][1]
            for r in M["P_sham_sets"]["較正"]["散らばりの比"]))
    pd.DataFrame([{"set": r["label"], "z": r["set_mean_z"], "pass": r["pass"],
                   "genes": " ".join(r["genes"])} for r in rows]).to_csv(
        outdir / "set_sham_005.csv", index=False)

    # ---------------- 入れ子の感度分析（**主判定の置き換えには使わない**）
    power = json.loads((E / "work/power_005.json").read_text())["sets"]
    sens = {}
    for name in a["sensitivity_sets"]:
        if name not in tests:
            continue
        sens[name] = {"pass": tests[name]["pass"], "k": tests[name]["k"],
                      "set_mean_z": tests[name]["set_mean_z"],
                      "検出下限": power.get(name, {}).get("detection_floor_with_gene_sd"),
                      "検出力が足りない": name in a["underpowered_sets"]}
    for ver in a["secondary_null_versions"]:
        for name in (PRIMARY, PC):
            if sets_run[name]:
                sens[f"{name}（対照 {ver}）"] = {
                    "pass": set_test(ver, "goal", name, sets_run[name], ABS.get(name, -1e9))["pass"]}
    sens["中立対照（H-3 を 0.55 にした版）"] = neu.get("neutral_adopted_cv055")
    M["sensitivity"] = sens
    M["underpowered_note"] = (
        "感度 (b)（INT_S2、転写標的 13 個）は検出下限 1.84 で検出力が足りない。"
        "**不合格を「効果がない」の根拠に使わない**（plan.md §9）。")

    # ---------------- 記述のみ（合否に入れない）
    zg = Zt[(PRIMARY_NULL, "goal")]
    zv = lambda g, s: float(zg[(zg.gene == g) & (zg.seed == s)].z.iloc[0])
    desc = {
        "前回の主指標（寄りの差。**判定に使わない**）": {
            name: set_test(PRIMARY_NULL, "desc_diff", name, sets_run[name], -1e9)["set_mean_z"]
            for name in (PRIMARY, PC) if sets_run[name]},
        "退化した 2 状態の軸での寄り（**判定に使わない**。深さがそろっていない）": {
            f"desc_{aa}": set_test(PRIMARY_NULL, f"desc_{aa}", PRIMARY, sets_run[PRIMARY],
                                   -1e9)["set_mean_z"]
            for aa in cfg["states"].get("alt_states", []) if sets_run[PRIMARY]},
        "移動量 ld を入れた版（版の比較。**判定に使わない**）": {
            name: [round(float(np.mean([float(Z_ld[(Z_ld.gene == g) & (Z_ld.seed == s)].z.iloc[0])
                                        for g in sets_run[name]])), 3) for s in seeds]
            for name in (PRIMARY, PC, "housekeeping") if sets_run[name]},
        "落選した中立対照の単一遺伝子の順位（**記述のみ**）": {
            sym[g]: [qv(g, s) for s in seeds]
            for g in sets_run["housekeeping"] + sets_run["negctrl_panel"]
            if g not in set(sets_run[NEU])},
        "単一遺伝子で 3 seed すべて 2.5 標準偏差以上": {
            name: sorted(sym[g] for g in genes if all(zv(g, s) >= 2.5 for s in seeds))
            for name, genes in sets_run.items() if genes},
        "発現の偏り（bias_goal）と z の順位相関（介入セット）": {
            str(s): float(stats.spearmanr([bias[sym[g]] for g in sets_run[PRIMARY]],
                                          [zv(g, s) for g in sets_run[PRIMARY]]).statistic)
            for s in seeds} if len(sets_run[PRIMARY]) > 2 else None}
    M["descriptive"] = desc

    # ---------------- 判定の提案（**確定は PI の関門 G3**）
    ok = bool(M["metric_pass"] and M["P1_pass"] and M["P4_pass"] and M["P5_pass"]
              and M["P_sham_sets"]["pass"] and M["chance_rate_ok"])
    floor = power[PRIMARY]["detection_floor_with_gene_sd"]
    M["conclusion_proposal"] = (
        f"判定不能（目標の軸の陽性対照 P-4（{PC}）が不合格。目標の軸の目盛りが立っていないので、"
        "主判定の値を解釈しない）" if not M["P4_pass"] else
        "判定不能（中立のはずのセットが動いた → I-2）" if not M["P5_pass"] else
        "判定不能（指標の健全性が不合格）" if not M["metric_pass"] else
        "判定不能（偽のセットが偶然の範囲を超えて合格した）" if not M["P_sham_sets"]["pass"] else
        "棄却（R-3: 寄る向きが逆の seed がある）" if M["R3_direction_reversed"] else
        "主予測 P-1 は支持される（in silico の摂動として。実験的事実ではない）" if ok else
        f"主予測 P-1 は不合格。**「効果がない」とは書けない。**"
        f"「検出下限 {floor} 標準偏差以上の効果は認められなかった」としか読めない")
    M["note"] = (
        "判定の確定は PI（関門 G3）。本スクリプトは事前登録の規則を当てはめて提案するだけ。"
        "**この実験で答えられるのは「線維芽細胞が活性化した間質の状態へ動くか」であって、"
        "「三次リンパ組織が作られるか」ではない**（HYP-20261001-01 §0）。"
        "陽性対照 71 個は健全性チェックであって、仮説を支持する証拠ではない。"
        "検定の単位はセットであり、個々の遺伝子が効くことは示さない。"
        "盲検中の TP53・SPP1・LCN2 は摂動も解析もしていない。"
        "in silico の結果であり、実験的に確認された生物学的事実ではない。")

    (outdir / "criteria_005.json").write_text(
        json.dumps(M, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps({k: v for k, v in M.items()
                      if k.endswith("_pass") or k in ("P1_primary", "P_sham_sets",
                                                      "chance_pass_rate", "sensitivity",
                                                      "conclusion_proposal", "zero_movement",
                                                      "metric_version", "M3_n_in_range")},
                     ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoint", action="store_true")
    ap.add_argument("--results-dir", default=None,
                    help="試し実行の出力で評価コードを通すとき results-smoke を指定する")
    ap.add_argument("--smoke", action="store_true",
                    help="試し実行の出力に対する実行（下限や繰り返し回数を緩める。**判定ではない**）")
    args = ap.parse_args()
    main(checkpoint=args.checkpoint, results_dir=args.results_dir, smoke=args.smoke)

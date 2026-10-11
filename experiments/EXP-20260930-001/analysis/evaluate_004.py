#!/usr/bin/env python3
"""EXP-20260930-001 — 事前登録した判定の機械的な当てはめ（plan.md §6 をそのまま実装する）。

EXP-20260929-001 からの変更点（**すべて結果を見る前に決めた**）:

1. **主指標の差し引きに退化した状態（dFIB・dMFIB）の重心を入れない。**
   この 2 本は入力の長さの中央値が 6.71 / 6.50 で、目標・出発点の 7.42 と大きく違う
   （`work/depth_match_summary_004.json`）。「その軸では判定しない」と決めても、
   主指標の式で差し引けば結果に効く（guardrails 2026-09-30）。
   主指標は **(目標への接近 − 出発点への接近) / 移動量** の 1 本だけにする。
2. **移動量の補正を主指標の定義に入れる**（`scripts/analysis/movement_correction.py`）。
   説明変数は lp・ll・**ld = log2(移動量)** の完全二次式。
   **摂動が起きなかった細胞（移動量 0）は外す。比べる版すべてで同じ細胞を外す。**
3. **偽のセットの合格条件を二項検定にする。**「N 組すべて不合格」は、問題のない設計でも
   67% の確率で落ちる条件だった（plan.md §6.3）。
4. 対照遺伝子群を引き直した（陽性対照との重複を 0 にするため。plan.md §5.3）。
5. **目標の軸の陽性対照の主判定を、最小 42 個版（使える 39 個）から保守 76 個版（使える 71 個）に
   替えた**（**PI 判断 2026-09-30**。検出下限が 1.03 → 0.78 に下がるため。plan.md §9.4）。
   39 個版は 71 個版の部分集合なので、入れ子の感度分析として同じ 1 回の実行で評価できる。

**途中のチェックポイントでは主判定のセットの値を計算も保存もしない**（`--checkpoint`。plan.md §6.4）。

判定対象にしないもの: 盲検中の TP53 / SPP1 / LCN2（読み込み直後に除外し、1 行でも残っていれば停止）。
LTBR・CDKN1A は PI が 2026-09-27 に解除済みなので摂動も解析も行う。

出力: results/eval/criteria_004.json, results/eval/gene_seed_004.csv, results/eval/set_sham_004.csv
実行: .venv/bin/python experiments/EXP-20260930-001/analysis/evaluate_004.py [--checkpoint]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260930-001"
sys.path.insert(0, str(ROOT / "scripts/analysis"))
from movement_correction import gene_stats  # noqa: E402

# 盲検中（摂動も解析もしない）。LTBR・CDKN1A は解除済みなので入れない。
BLOCKED = {"ENSG00000141510": "TP53", "ENSG00000118785": "SPP1", "ENSG00000148346": "LCN2"}

# チェックポイントで見てよいセット（plan.md §6.4）。**主判定のセットはここに入れない。**
CHECKPOINT_SETS = ("positive_control71", "positive_control39", "positive_control65",
                   "positive_control70", "positive_control38", "positive_control35",
                   "housekeeping", "negctrl_panel")


def cos(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return (a @ b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b) + 1e-12)


def eligibility(cs: pd.DataFrame, nulls: dict, tol: float, cover: float) -> dict[str, list[str]]:
    """遺伝子ごとに「元の位置と入力の長さで比べられる対照遺伝子」を数える。
    EXP-20260928-001 / EXP-20260929-001 と同一の実装。"""
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


def load(cfg: dict) -> tuple[pd.DataFrame, dict]:
    """細胞ごとの表を作る。**摂動が起きなかった細胞（移動量 0）を外す**（plan.md §6.1）。"""
    import pickle
    from datasets import load_from_disk
    from geneformer import TOKEN_DICTIONARY_FILE

    res = E / cfg["output"]["results_dir"]
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
        dnorm[idx] = np.linalg.norm(p - o, axis=1)
    for n in names:
        c[f"s_{n}"] = vals[n]
    c["dnorm"] = dnorm

    # **移動量 0 の細胞を外す**（過剰発現で入力が 1 ビットも変わらなかった細胞。guardrails 2026-09-30）
    zero = c.dnorm <= 0
    drop = {"n_cells_before": int(len(c)), "n_zero_movement": int(zero.sum()),
            "frac_zero_movement": round(float(zero.mean()), 6),
            "genes_with_zero": sorted(c.loc[zero, "gene"].unique().tolist())[:20]}
    c = c[~zero].reset_index(drop=True)
    drop["n_cells_after"] = int(len(c))
    return c, drop


def metric_columns(c: pd.DataFrame, alt: list[str]) -> dict[str, str]:
    """**主指標は (目標への接近 − 出発点への接近) / 移動量 の 1 本だけ**（plan.md §6.1）。
    退化した状態の重心（深さがそろっていない）は差し引きに入れない。"""
    dn = c.dnorm
    c["m_goal"] = (c.s_goal - c.s_start) / dn                       # **主指標**
    c["m_uniform"] = (c.s_goal + c.s_start) / 2.0                   # 一様成分（健全性 M4）
    c["m_A11"] = c.s_goal - c.s_start                               # 移動量で割らない版（記述）
    cols = {"goal": "m_goal", "uniform": "m_uniform", "A11": "m_A11"}
    for a in alt:                                                   # **記述のみ。判定に使わない**
        c[f"m_desc_{a}"] = (c[f"s_{a}"] - c.s_start) / dn
        cols[f"desc_{a}"] = f"m_desc_{a}"
    if alt:
        c["m_desc_goal_minus_all"] = (c.s_goal - c[["s_start"] + [f"s_{a}" for a in alt]].mean(axis=1)) / dn
        cols["desc_goal_minus_all"] = "m_desc_goal_minus_all"
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


def main(checkpoint: bool = False) -> None:
    cfg = yaml.safe_load(open(E / "config.yaml"))
    a = cfg["analysis"]
    cr = a["criteria"]
    mc = cfg["metric"]["confound_correction"]
    gj = json.loads((E / "work/genes_004.json").read_text())["genes"]
    roles = {g["symbol"]: g["role"] for g in gj}
    ens = {g["symbol"]: g["ensembl"] for g in gj}
    bias = {g["symbol"]: g["bias_goal"] for g in gj}
    sym = {v: k for k, v in ens.items()}

    npool = json.loads((E / "work/null_pool_004.json").read_text())["versions"]
    null_versions = {k: [ens[x["symbol"]] for x in v] for k, v in npool.items()}
    subs = json.loads((E / "work/subsets_004.json").read_text())["subsets"]
    sets_def = {k: [ens[x["symbol"]] for x in v] for k, v in subs.items()}
    PRIMARY, PRIMARY_NULL = a["primary_set"], a["primary_null_version"]
    ABS = a["set_abs_min"]                      # セットごとの絶対値条件（**結果を見る前に固定**）

    assert cr["require_all_seeds"] is True, "config の require_all_seeds が False（plan.md §6 と食い違う）"
    assert len(sets_def[PRIMARY]) == cr["K_set"], "主判定のセットの大きさが config と違う"
    assert set(null_versions["A_300"]) < set(null_versions["B_500"]), "対照の部分集合が入れ子でない"
    for name, genes in sets_def.items():
        assert not (set(null_versions["B_500"]) & set(genes)), f"対照と {name} が重複している"
    assert mc["variables"] == ["lp", "ll", "ld"], "主判定の補正の版が事前登録と違う"
    assert mc["drop_zero_movement"] is True, "移動量 0 の細胞を外す設定になっていない"

    c, drop = load(cfg)
    alt = list(cfg["states"].get("alt_states", []))
    cols = metric_columns(c, alt)
    seeds = sorted(c.seed.unique().tolist())
    if not checkpoint:
        assert seeds == sorted(cfg["robustness"]["seeds"]), "seed が計画と違う"

    pm = a["pos_matched_permutation"]
    nulls_all = set(null_versions["B_500"])
    frames: dict[str, list] = {}
    frames_old: dict[str, list] = {}
    for s in seeds:
        cs = c[c.seed == s]
        nl = {h: hd for h, hd in cs[cs.gene.isin(nulls_all)].groupby("gene")}
        el = eligibility(cs, nl, pm["tolerance_log2"], pm["min_cover_frac"])
        for name, col in cols.items():
            f = gene_stats(cs, nl, el, col, pm["min_eligible_null_genes"], use_dnorm=True)
            f["seed"] = s
            frames.setdefault(name, []).append(f)
            if name == "goal":     # **記述のみ**: 移動量を入れない従来の版（版の比較用）
                fo = gene_stats(cs, nl, el, col, pm["min_eligible_null_genes"], use_dnorm=False)
                fo["seed"] = s
                frames_old.setdefault(name, []).append(fo)
    tabs = {k: pd.concat(v, ignore_index=True) for k, v in frames.items()}
    tab_old = pd.concat(frames_old["goal"], ignore_index=True)
    long = pd.concat(tabs.values(), ignore_index=True)
    long["symbol"] = long.gene.map(sym)
    long["role"] = long.symbol.map(roles)
    (E / "results/eval").mkdir(parents=True, exist_ok=True)
    long.to_csv(E / "results/eval/gene_seed_004.csv", index=False)

    def standardise(T: pd.DataFrame, nset: set) -> pd.DataFrame:
        T = T.copy()
        for s in seeds:
            m = T.seed == s
            nul = T[m & T.gene.isin(nset)]
            T.loc[m, "z"] = (T.loc[m, "D_raw"] - nul.D_raw.mean()) / nul.D_raw.std(ddof=1)
        T["symbol"] = T.gene.map(sym)
        T["role"] = T.symbol.map(roles)
        return T

    Zt = {(ver, ax): standardise(tabs[ax], set(gs))
          for ver, gs in null_versions.items() for ax in cols}
    Z_old = standardise(tab_old, set(null_versions[PRIMARY_NULL]))

    def zmat(ver: str, axis: str, genes: list[str]) -> np.ndarray:
        T = Zt[(ver, axis)]
        return np.column_stack([[float(T[(T.gene == g) & (T.seed == s)].z.iloc[0]) for g in genes]
                                for s in seeds])

    rng = np.random.default_rng(cr["null_draw_seed"])
    znull = {(v, ax): zmat(v, ax, sorted(null_versions[v])) for v in null_versions for ax in cols}
    _cache: dict[tuple, np.ndarray] = {}

    def null_dist(ver: str, axis: str, k: int, replace: bool) -> np.ndarray:
        key = (ver, axis, k, replace)
        if key not in _cache:
            _cache[key] = set_null(znull[(ver, axis)], k, cr["n_permutation"], rng, replace)
        return _cache[key]

    def set_test(ver: str, axis: str, name: str, genes: list[str], abs_min: float) -> dict:
        k = len(genes)
        obs = zmat(ver, axis, genes).mean(axis=0)
        out = {"label": name, "axis": axis, "null_version": ver, "k": k, "abs_min": abs_min,
               "set_mean_z": [round(float(x), 3) for x in obs]}
        for rep in (True, False):
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
        nd = set_null(znull[(ver, axis)], k, cr["n_chance_draw"], rng, rep)   # 閾値とは別に引き直す
        hit = np.ones(len(nd), bool)
        for i in range(len(seeds)):
            hit &= (nd[:, i] > thr[i]) & (nd[:, i] >= abs_min)
        return {"n_draw": int(len(nd)), "rate": float(hit.mean()),
                "threshold": [round(x, 3) for x in thr]}

    M: dict = {"mode": "checkpoint" if checkpoint else "final", "seeds": seeds,
               "n_genes": len(gj), "zero_movement": drop,
               "correction_version": {"variables": mc["variables"], "form": mc["form"],
                                      "drop_zero_movement": mc["drop_zero_movement"]},
               "n_nulls": {k: len(v) for k, v in null_versions.items()},
               "primary_set": PRIMARY, "primary_null_version": PRIMARY_NULL}

    # ---------------- 欠測（比べる相手が足りずに計算できなかった遺伝子）は「不合格」と数えない
    P = tabs["goal"]
    M["missing_q_by_seed"] = {int(k): int(v)
                              for k, v in P[P.q.isna()].groupby("seed").gene.count().to_dict().items()}

    # ---------------- 指標の健全性 M1〜M4（EXP-20260927-004 以降と同一の閾値）
    NULP = P[P.gene.isin(set(null_versions[PRIMARY_NULL]))]
    M["M1_calibration"] = {}
    for s in seeds:
        q = NULP[NULP.seed == s].q.dropna()
        M["M1_calibration"][s] = {"n_valid": int(len(q)), "frac_hi": float((q >= 0.90).mean()),
                                  "frac_lo": float((q < 0.10).mean())}
    M["M1_pass"] = all(v["n_valid"] >= cr["M1_min_n"] and v["frac_hi"] <= cr["M1_tail_frac_max"]
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
    M["M3_housekeeping"] = {sym[g]: [qv(g, s) for s in seeds] for g in sets_def["housekeeping"]}
    M["M3_panel"] = {sym[g]: [qv(g, s) for s in seeds] for g in sets_def["negctrl_panel"]}
    M["M3_pass"] = (all(inr(v) for v in M["M3_housekeeping"].values())
                    and sum(inr(v) for v in M["M3_panel"].values()) >= cr["M3_panel_min_genes"])
    M["M4_rho"] = {}
    for s in seeds:
        n = NULP[NULP.seed == s].set_index("gene")
        u = tabs["uniform"][tabs["uniform"].seed == s].set_index("gene").reindex(n.index)
        M["M4_rho"][s] = {"uniform": float(stats.spearmanr(n.D_raw, u.D_raw).statistic)}
    M["M4_pass"] = all(abs(v["uniform"]) <= cr["M4_abs_rho_uniform_max"] for v in M["M4_rho"].values())
    # 残差と移動量の順位相関（補正が効いているかの確認。plan.md §6.1）
    M["M5_residual_vs_movement"] = {}
    for s in seeds:
        n = NULP[NULP.seed == s]
        dmed = c[c.seed == s].groupby("gene").dnorm.median()
        M["M5_residual_vs_movement"][s] = float(
            stats.spearmanr(dmed.reindex(n.gene).to_numpy(), n.D_raw.to_numpy()).statistic)
    M["M5_pass"] = all(abs(v) <= cr["M5_abs_rho_movement_max"] for v in M["M5_residual_vs_movement"].values())
    M["metric_pass"] = bool(M["M1_pass"] and M["M2_pass"] and M["M3_pass"] and M["M4_pass"]
                            and M["M5_pass"])

    # ---------------- P-4 健全性チェック（陽性対照）と P-5 中立性
    tests: dict[str, dict] = {}
    for name in CHECKPOINT_SETS:
        tests[name] = set_test(PRIMARY_NULL, "goal", name, sets_def[name], ABS.get(name, -1e9))
    PC = a["primary_positive_control"]               # **PI 判断 2026-09-30: 保守 71 個版**
    assert PC == "positive_control71", "主判定の陽性対照が事前登録と違う"
    M["P4_positive_control"] = tests[PC]
    M["P4_pass"] = tests[PC]["pass"]
    M["P4_sensitivity"] = {k: tests[k] for k in
                           ("positive_control39", "positive_control65", "positive_control70",
                            "positive_control38", "positive_control35")}
    neu = {}
    for name in ("housekeeping", "negctrl_panel"):
        k = len(sets_def[name])
        obs = zmat(PRIMARY_NULL, "goal", sets_def[name]).mean(axis=0)
        nd = null_dist(PRIMARY_NULL, "goal", k, cr["set_draw_replace"])
        pc = [float((nd[:, i] < obs[i]).mean()) for i in range(len(seeds))]
        neu[name] = {"set_mean_z": [round(float(x), 3) for x in obs],
                     "null_percentile": [round(x, 4) for x in pc],
                     "n_seeds_in_range": int(sum(cr["neutral_range"][0] <= x <= cr["neutral_range"][1]
                                                 for x in pc))}
    M["P5_neutral"] = neu
    M["P5_pass"] = all(v["n_seeds_in_range"] >= cr["neutral_min_seeds"] for v in neu.values())

    if checkpoint:
        # **主判定のセットの値は計算も保存もしない**（plan.md §6.4）
        M["note"] = ("チェックポイント。**主判定のセット（T1・T2・T1∪T2）の値は計算も保存もしていない。**"
                     "見てよいのは装置が働いているかの確認だけ（陽性対照・中立性・指標の健全性）。")
        (E / "results/eval/checkpoint_004.json").write_text(
            json.dumps(M, ensure_ascii=False, indent=1, default=float) + "\n")
        print(json.dumps({k: v for k, v in M.items() if k.endswith("_pass")}, ensure_ascii=False, indent=1))
        return

    # ---------------- 主判定 P-1 と副次 P-2 / P-3
    for name in ("T1", "T2", "T1_T2", "T1_no_CXCL12"):
        tests[name] = set_test(PRIMARY_NULL, "goal", name, sets_def[name], ABS[name])
    M["set_tests"] = tests
    M["P1_primary"] = tests["T1"]
    M["P1_pass"] = tests["T1"]["pass"]
    M["P2_pass"] = tests["T2"]["pass"]
    M["P3_pass"] = tests["T1_T2"]["pass"]
    M["R3_direction_reversed"] = not tests["T1"]["all_seeds_positive"]

    # ---------------- 偶然合格率（対照で実測）と偽のセットの二項検定（plan.md §6.3）
    M["chance_pass_rate"] = {name: chance_rate(PRIMARY_NULL, "goal", len(sets_def[name]), ABS[name])
                             for name in ("T1", "T2", "T1_T2", PC)}
    M["chance_rate_ok"] = bool(M["chance_pass_rate"]["T1"]["rate"] <= cr["max_chance_rate"])

    def tbl(genes: set) -> pd.DataFrame:
        t = P[(P.seed == seeds[0]) & (P.gene.isin(genes))][["gene", "lp_med", "ll_med", "n_cells"]].copy()
        t["symbol"] = t.gene.map(sym)
        t["bias_goal"] = t.symbol.map(bias)
        return t

    npl = tbl(set(null_versions[PRIMARY_NULL]))
    shams = sham_sets(npl, tbl(set(sets_def["T1"])), len(sets_def["T1"]), cr["n_sham_sets"],
                      np.random.default_rng(cr["null_draw_seed"] + 1), cr["sham_match_bias"])
    rows = []
    for i, nm in enumerate(shams):
        r = set_test(PRIMARY_NULL, "goal", f"sham_{i + 1}", [ens[x] for x in nm], ABS["T1"])
        r["genes"] = nm
        rows.append(r)
    kk = int(sum(r["pass"] for r in rows))
    N = len(rows)
    p = float(M["chance_pass_rate"]["T1"]["rate"])
    pval = float(stats.binomtest(kk, N, max(p, 1e-12), alternative="greater").pvalue)
    sm = np.array([r["set_mean_z"] for r in rows], dtype=float)
    rnd = null_dist(PRIMARY_NULL, "goal", len(sets_def["T1"]), cr["set_draw_replace"])
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
        E / "results/eval/set_sham_004.csv", index=False)

    # ---------------- 入れ子の感度分析（主判定の置き換えには使わない）
    sens = {"T1 から CXCL12 を除いた 10 個": tests["T1_no_CXCL12"]["pass"],
            "陽性対照 39 個版（最小 42 版。PI 判断の前の主判定）": tests["positive_control39"]["pass"],
            "陽性対照 65 個版（71 から逆向きの 6 個を除く）": tests["positive_control65"]["pass"],
            "陽性対照 70 個版（71 から LTBP1 を除く）": tests["positive_control70"]["pass"],
            "陽性対照 35 個版（39 から逆向きの 4 個を除く）": tests["positive_control35"]["pass"],
            "陽性対照 38 個版（39 から LTBP1 を除く）": tests["positive_control38"]["pass"]}
    for ver in a["secondary_null_versions"]:
        for name in ("T1", "T2", "T1_T2", PC):
            sens[f"{name}（対照 {ver}）"] = set_test(ver, "goal", name, sets_def[name],
                                                    ABS[name])["pass"]
    M["sensitivity"] = sens

    # ---------------- 記述のみ（合否に入れない）
    zg = Zt[(PRIMARY_NULL, "goal")]
    zv = lambda g, s: float(zg[(zg.gene == g) & (zg.seed == s)].z.iloc[0])
    desc = {"退化した線維芽細胞の軸・髄質の退化した線維芽細胞の軸での寄り（**判定に使わない**。"
            "深さがそろっていないため）": {
                f"desc_{aa}": set_test(PRIMARY_NULL, f"desc_{aa}", "T1", sets_def["T1"], -1e9)["set_mean_z"]
                for aa in alt},
            "退化の 2 本も差し引いた版（**判定に使わない**）":
                set_test(PRIMARY_NULL, "desc_goal_minus_all", "T1", sets_def["T1"], -1e9)["set_mean_z"]
                if alt else None,
            "移動量を補正しない従来の版（版の比較）": {
                name: [round(float(np.mean([float(Z_old[(Z_old.gene == g) & (Z_old.seed == s)].z.iloc[0])
                                            for g in sets_def[name]])), 3) for s in seeds]
                for name in ("T1", "positive_control39", "housekeeping")},
            "単一遺伝子で 3 seed すべて 2.5 標準偏差以上": {
                name: sorted(sym[g] for g in genes if all(zv(g, s) >= 2.5 for s in seeds))
                for name, genes in sets_def.items()},
            "恒常発現の追加 1（PPIA。記述のみ）": {
                "z": [round(zv(ens["PPIA"], s), 3) for s in seeds]} if "PPIA" in ens else None,
            "発現の偏り（bias_goal）と z の順位相関（T1、目標の軸）": {
                str(s): float(stats.spearmanr([bias[sym[g]] for g in sets_def["T1"]],
                                              [zv(g, s) for g in sets_def["T1"]]).statistic)
                for s in seeds}}
    M["descriptive"] = desc

    # ---------------- 判定の提案（**確定は PI の関門 G3**）
    ok = bool(M["metric_pass"] and M["P1_pass"] and M["P4_pass"] and M["P5_pass"]
              and M["P_sham_sets"]["pass"] and M["chance_rate_ok"])
    floor = ABS["T1"]
    M["conclusion_proposal"] = (
        f"判定不能（目標の軸の陽性対照 P-4（{PC}）が不合格。目標の軸の目盛りが立っていないので、"
        "主判定の値を解釈しない）" if not M["P4_pass"] else
        "判定不能（中立のはずのセットが動いた → I-2）" if not M["P5_pass"] else
        "判定不能（指標の健全性が不合格）" if not M["metric_pass"] else
        "判定不能（偽のセットが偶然の範囲を超えて合格した）" if not M["P_sham_sets"]["pass"] else
        "棄却（R-3: 寄る向きが逆の seed がある）" if M["R3_direction_reversed"] else
        "主予測 P-1 は支持される（in silico の摂動として。実験的事実ではない）" if ok else
        f"主予測 P-1 は不合格。**「効果がない」とは書けない。**"
        f"「検出下限 {floor} 標準偏差以上の効果は認められなかった」としか読めない（HYP §10-2）")
    M["note"] = (
        "判定の確定は PI（関門 G3）。本スクリプトは事前登録の規則を当てはめて提案するだけ。"
        "**この実験で答えられるのは「線維芽細胞が活性化した間質の状態へ動くか」であって、"
        "「三次リンパ組織が作られるか」ではない**（HYP §0）。"
        "陽性対照 39 個は健全性チェックであって、仮説を支持する証拠ではない。"
        "検定の単位はセットであり、個々の遺伝子が効くことは示さない。"
        "盲検中の TP53・SPP1・LCN2 は摂動も解析もしていない。"
        "in silico の結果であり、実験的に確認された生物学的事実ではない。")

    (E / "results/eval/criteria_004.json").write_text(
        json.dumps(M, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps({k: v for k, v in M.items()
                      if k.endswith("_pass") or k in ("P1_primary", "P_sham_sets",
                                                      "chance_pass_rate", "sensitivity",
                                                      "conclusion_proposal", "zero_movement")},
                     ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    main(checkpoint="--checkpoint" in sys.argv)

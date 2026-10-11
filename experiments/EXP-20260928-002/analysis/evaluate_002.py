"""EXP-20260928-002 — 事前登録した判定の機械的な当てはめ（plan.md §6 をそのまま実装する）。

EXP-20260928-001 からの変更点:
  判定するのは「経路セット 120 個」の平均。対照遺伝子の上位集合は 980 個で、
  そこから入れ子の部分集合 600（案 B 相当）・520（案 C 相当）を事前に定義してある。
  **主判定は 980 版**、600 / 520 版は副次的な報告（3 つの版で判定が変わるかを見る）。
  無作為な対照セットは **復元ありで引く**（guardrails 2026-09-28）。復元なしの閾値も併記する。
  偽経路セットは 50 組以上つくり、元の位置・入力の長さ・細胞数に加えて bias_goal もそろえる。

判定対象にしないもの（読み込み直後に除外し、1 行でも残っていれば停止する）:
  元の主対象 LTBR / RELB / NFKB2、盲検の CDKN1A / TP53 / SPP1 / LCN2、
  目標ラベルを定義する CCL19 / CCL21 / CXCL13。

出力: results/eval/criteria_002.json, results/eval/gene_seed_002.csv, results/eval/set_null_002.csv
実行: .venv/bin/python experiments/EXP-20260928-002/analysis/evaluate_002.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260928-002"

BLOCKED = {"ENSG00000111321": "LTBR", "ENSG00000104856": "RELB", "ENSG00000077150": "NFKB2",
           "ENSG00000124762": "CDKN1A", "ENSG00000141510": "TP53",
           "ENSG00000118785": "SPP1", "ENSG00000148346": "LCN2",
           "ENSG00000172724": "CCL19", "ENSG00000137077": "CCL21", "ENSG00000156234": "CXCL13"}


# ---------------------------------------------------------------- 指標
def cos(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return (a @ b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b) + 1e-12)


def metric_columns(c: pd.DataFrame) -> dict[str, str]:
    """4 つの重心から、同じ関数形・同じ尺度で軸だけを入れ替えた指標を作る。
    A12 の形（目標に特異的な成分 / 移動量）。EXP-20260927-004 / EXP-20260928-001 と同じ定義。"""
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
    """遺伝子ごとに、対照遺伝子（自分自身は除く）で当てた二次回帰からの残差の平均 D_raw と、
    対照遺伝子の中での順位 q を返す。

    **EXP-20260928-001 の実装（遺伝子ごとに lstsq を解き直し、対照ごとに残差を全細胞ぶん計算する）と
    数学的に同一だが、対照が 980 個・遺伝子が 1,136 個になると桁違いに遅くなるので、
    (1) 最小二乗を正規方程式にして、対照全体の X'X・X'y から自分の寄与を引く、
    (2) 残差の平均は mean(y) − mean(X)·beta で計算する（平均だけが必要なので細胞ごとに引く必要がない）、
    の 2 つで書き換えてある。self_test() で元の実装との一致（最大差 < 1e-9）を確かめている。**"""
    hs = list(nulls)
    Xn = {h: design(hd.lp, hd.ll) for h, hd in nulls.items()}
    yn = {h: hd[col].to_numpy() for h, hd in nulls.items()}
    G = {h: Xn[h].T @ Xn[h] for h in hs}
    B = {h: Xn[h].T @ yn[h] for h in hs}
    Xbar = np.vstack([Xn[h].mean(axis=0) for h in hs])          # 対照ごとの design の平均（N × 6）
    ybar = np.array([yn[h].mean() for h in hs])                 # 対照ごとの y の平均（N,）
    hi = {h: i for i, h in enumerate(hs)}
    G_all = sum(G.values())
    B_all = sum(B.values())
    rows = []
    for gene, tg in cs.groupby("gene"):
        if gene in G:                                            # 自分が対照なら自分の寄与を引く
            Gm, Bm = G_all - G[gene], B_all - B[gene]
        else:
            Gm, Bm = G_all, B_all
        beta = np.linalg.solve(Gm, Bm)
        Xg = design(tg.lp, tg.ll)
        Dg = float(tg[col].to_numpy().mean() - Xg.mean(axis=0) @ beta)
        e = elig[gene]
        if e:
            ii = np.array([hi[h] for h in e])
            Dh = ybar[ii] - Xbar[ii] @ beta
        else:
            Dh = np.zeros(0)
        q = float((Dh < Dg).mean()) if len(Dh) >= min_elig else np.nan
        rows.append({"gene": gene, "col": col, "D_raw": Dg, "q": q, "n_elig": len(Dh),
                     "n_cells": len(tg), "lp_med": float(tg.lp.median()),
                     "ll_med": float(tg.ll.median()), "relpos_med": float(tg.relpos.median())})
    return pd.DataFrame(rows)


def gene_stats_reference(cs: pd.DataFrame, nulls: dict, elig: dict, col: str, min_elig: int) -> pd.DataFrame:
    """EXP-20260928-001 の実装そのまま（遅い）。self_test() の照合にだけ使う。"""
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


def self_test() -> None:
    """合成データで (a) 速い実装と元の実装の一致、(b) set_null の分布、(c) 偽セットの作り方 を確かめる。
    実行: .venv/bin/python experiments/EXP-20260928-002/analysis/evaluate_002.py --self-test"""
    rng = np.random.default_rng(0)
    rows = []
    for gi in range(40):
        n = int(rng.integers(120, 201))
        lp = rng.normal(5, 1, n); ll = rng.normal(11, 0.3, n)
        rows.append(pd.DataFrame({"gene": f"G{gi:03d}", "lp": lp, "ll": ll,
                                  "m": 0.3 * lp - 0.2 * ll + rng.normal(0, 1, n),
                                  "relpos": rng.random(n)}))
    cs = pd.concat(rows, ignore_index=True)
    nulls = {h: hd for h, hd in cs[cs.gene >= "G010"].groupby("gene")}
    el = eligibility(cs, nulls, 0.5, 0.8)
    a = gene_stats(cs, nulls, el, "m", 5).sort_values("gene").reset_index(drop=True)
    b = gene_stats_reference(cs, nulls, el, "m", 5).sort_values("gene").reset_index(drop=True)
    d = max(float(np.abs(a.D_raw - b.D_raw).max()), float(np.nanmax(np.abs(a.q - b.q))))
    print(f"(a) 速い実装と元の実装の最大差: {d:.3e}  → {'一致' if d < 1e-9 else '不一致'}")
    assert d < 1e-9
    z = rng.normal(size=(980, 3))
    for rep in (True, False):
        m = set_null(z, 120, 20000, np.random.default_rng(1), rep)
        print(f"(b) set_null replace={rep}: 標準偏差 {m[:, 0].std():.4f} "
              f"(理論 {1 / np.sqrt(120) * (1 if rep else np.sqrt(1 - 120 / 980)):.4f})")
    print("自己検査 OK")


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
    # list.index を 68 万回呼ぶと非常に遅いので、細胞ごとに numpy 配列にしてから探す（結果は同じ）
    ids = {u: np.asarray(x) for u, x in zip(uids, ds.select(uids)["input_ids"])}
    pos = np.empty(len(c), dtype=np.int64)
    for gene, part in c.groupby("gene", sort=False):
        t = tokd[gene]
        pos[part.index.to_numpy()] = [int(np.flatnonzero(ids[u] == t)[0]) for u in part.cell_uid]
    c["pos"] = pos
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
def set_null(z: np.ndarray, k: int, n_draw: int, rng: np.random.Generator,
             replace: bool = True) -> np.ndarray:
    """対照遺伝子（N × seed）から k 個を取ったセットの z 平均（n_draw × seed）。

    replace=True（事前登録する引き方。guardrails 2026-09-28）: 復元あり。
      有限母集団の補正が入らないので、対照の数がセットの大きさに対して小さくても閾値が甘くならない。
    replace=False: 復元なし。比較のために併せて計算し、閾値の差を報告する。"""
    N = z.shape[0]
    if replace:
        idx = rng.integers(0, N, size=(n_draw, k))
    else:
        if k > N:
            raise ValueError(f"復元なしでは対照 {N} 個から {k} 個を引けない")
        idx = np.argsort(rng.random((n_draw, N)), axis=1)[:, :k]
    return z[idx].mean(axis=1)


def sham_sets(nulls: pd.DataFrame, target: pd.DataFrame, k: int, n_sets: int,
              rng: np.random.Generator, match_bias: bool = True) -> list[list[str]]:
    """偽経路セット: 対照遺伝子から、元の位置・入力の長さ・細胞数（と bias_goal）を
    経路セットにそろえて k 個を選ぶ。1 組でも合格したら、経路セットの合格は
    『経路であること』ではなく、そろえた性質で説明されることになる（HYP §6.2, F-5）。

    対照は 980 個、セットは 120 個なので、**組ごとに使用済みを外していくと 8 組で尽きる**。
    そこで組ごとに独立に選び直し（組の間で遺伝子の重複を許す）、代わりに
    各組の中では重複させない。50 組以上つくる（guardrails 2026-09-28）。"""
    out = []
    for _ in range(n_sets):
        pick: list[str] = []
        order_t = target.sample(frac=1.0, random_state=int(rng.integers(0, 2**31 - 1)))
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
            # 同率の近傍から 1 つ無作為に取る（組ごとに違うセットになるようにする）
            take = order[:5]
            pick.append(str(take[int(rng.integers(0, len(take)))]))
        out.append(pick[:k])
    return out


# ---------------------------------------------------------------- 判定
def main() -> None:
    cfg = yaml.safe_load(open(E / "config.yaml"))
    a = cfg["analysis"]
    cr = a["criteria"]
    gj = json.loads((E / "work/genes_002.json").read_text())["genes"]
    roles = {g["symbol"]: g["role"] for g in gj}
    ens = {g["symbol"]: g["ensembl"] for g in gj}
    bias = {g["symbol"]: g.get("bias_goal") for g in gj}
    sym = {v: k for k, v in ens.items()}

    npool = json.loads((E / "work/null_pool_002.json").read_text())["versions"]
    null_versions = {k: [ens[x["symbol"]] for x in v] for k, v in npool.items()}
    subs = json.loads((E / "work/subsets_002.json").read_text())["subsets"]
    sets_def = {k: [ens[x["symbol"]] for x in v] for k, v in subs.items()}
    sets_def |= {"myof_markers": [ens[s] for s, r in roles.items() if r == "marker_myof"],
                 "housekeeping": [ens[s] for s, r in roles.items() if r == "housekeeping"],
                 "negctrl_panel": [ens[s] for s, r in roles.items() if r == "negctrl_panel"]}
    PRIMARY, PRIMARY_NULL = a["primary_set"], a["primary_null_version"]
    assert cr["require_all_seeds"] is True, "config の require_all_seeds が False（plan.md §6 と食い違う）"
    assert len(sets_def[PRIMARY]) == cr["K_set"], "経路セットの大きさが config と違う"
    assert set(null_versions["C_520"]) < set(null_versions["B_600"]) < set(null_versions["A_980"]), \
        "対照の部分集合が入れ子になっていない"
    assert not (set(null_versions["A_980"]) & set(sets_def[PRIMARY])), "対照と経路セットが重複している"

    c = load(cfg)
    cols = metric_columns(c)
    seeds = sorted(c.seed.unique().tolist())
    assert seeds == sorted(cfg["robustness"]["seeds"]), "seed が計画と違う"

    pm = a["pos_matched_permutation"]
    nulls_all = set(null_versions["A_980"])
    frames = {}
    for s in seeds:
        cs = c[c.seed == s]
        nl = {h: hd for h, hd in cs[cs.gene.isin(nulls_all)].groupby("gene")}
        el = eligibility(cs, nl, pm["tolerance_log2"], pm["min_cover_frac"])
        for name, col in cols.items():
            f = gene_stats(cs, nl, el, col, pm["min_eligible_null_genes"])
            f["seed"] = s
            frames.setdefault(name, []).append(f)
    tabs = {k: pd.concat(v, ignore_index=True) for k, v in frames.items()}
    long = pd.concat(tabs.values(), ignore_index=True)
    long["symbol"] = long.gene.map(sym)
    long["role"] = long.symbol.map(roles)
    (E / "results/eval").mkdir(parents=True, exist_ok=True)
    long.to_csv(E / "results/eval/gene_seed_002.csv", index=False)

    # --- 軸ごとに z（対照遺伝子群で標準化）。版ごとに標準化の母集団が変わる
    Zt: dict[tuple[str, str], pd.DataFrame] = {}
    for ver, nulls_set in null_versions.items():
        nset = set(nulls_set)
        for axis in ("goal", "myofaxis", "placebo"):
            T = tabs[axis].copy()
            for s in seeds:
                m = T.seed == s
                nul = T[m & T.gene.isin(nset)]
                T.loc[m, "z"] = (T.loc[m, "D_raw"] - nul.D_raw.mean()) / nul.D_raw.std(ddof=1)
            T["symbol"] = T.gene.map(sym)
            T["role"] = T.symbol.map(roles)
            Zt[(ver, axis)] = T

    def zmat(ver: str, axis: str, genes: list[str]) -> np.ndarray:
        T = Zt[(ver, axis)]
        return np.column_stack([[float(T[(T.gene == g) & (T.seed == s)].z.iloc[0]) for g in genes]
                                for s in seeds])

    rng = np.random.default_rng(cr["null_draw_seed"])
    znull = {(v, ax): zmat(v, ax, sorted(null_versions[v])) for v in null_versions
             for ax in ("goal", "myofaxis", "placebo")}

    M: dict = {"seeds": seeds, "n_genes": len(gj),
               "n_nulls": {k: len(v) for k, v in null_versions.items()},
               "primary_set": PRIMARY, "primary_null_version": PRIMARY_NULL,
               "set_draw_replace": cr["set_draw_replace"]}

    # --- 欠測（比べる相手が足りずに計算できなかった遺伝子）を「不合格」と数えない（guardrails 2026-09-27）
    P = tabs["goal"]
    miss = P[P.q.isna()].groupby("seed").gene.count().to_dict()
    M["missing_q_by_seed"] = {int(k): int(v) for k, v in miss.items()}
    M["missing_in_pathway_set"] = sorted(
        sym[g] for g in sets_def[PRIMARY]
        if any(bool(P[(P.gene == g) & (P.seed == s)].q.isna().iloc[0]) for s in seeds))

    # --- M1〜M4（指標の健全性。EXP-20260927-004 / EXP-20260928-001 と同一の閾値）
    NULP = P[P.gene.isin(set(null_versions[PRIMARY_NULL]))]
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

    # 無作為な対照セットの分布は (版, 軸, K) ごとに 1 回だけ引き、すべての検定で同じものを使う
    _nullcache: dict[tuple, np.ndarray] = {}

    def null_dist(ver: str, axis: str, k: int, replace: bool) -> np.ndarray:
        key = (ver, axis, k, replace)
        if key not in _nullcache:
            _nullcache[key] = set_null(znull[(ver, axis)], k, cr["n_permutation"], rng, replace)
        return _nullcache[key]

    def set_test(ver: str, axis: str, genes: list[str], label: str) -> dict:
        k = len(genes)
        obs = zmat(ver, axis, genes).mean(axis=0)
        out = {"label": label, "axis": axis, "null_version": ver, "k": k,
               "set_mean_z": [round(float(x), 3) for x in obs]}
        for rep in (True, False):
            nd = null_dist(ver, axis, k, rep)
            thr = [float(np.quantile(nd[:, i], cr["set_percentile"])) for i in range(len(seeds))]
            pct = [float((nd[:, i] < obs[i]).mean()) for i in range(len(seeds))]
            passed = bool(all(o > t for o, t in zip(obs, thr)) and all(o >= cr["set_abs_min"] for o in obs))
            tag = "replace" if rep else "no_replace"
            out[f"null_threshold_{tag}"] = [round(x, 3) for x in thr]
            out[f"null_percentile_{tag}"] = [round(x, 5) for x in pct]
            out[f"pass_{tag}"] = passed
        out["pass"] = out["pass_replace"] if cr["set_draw_replace"] else out["pass_no_replace"]
        return out

    # G-1 主判定: 経路セット 120 が目標の軸で、対照 980 版の規則に合格する
    res = {}
    for name in sets_def:
        for axis in ("goal", "myofaxis"):
            res[f"{name}@{axis}@{PRIMARY_NULL}"] = set_test(PRIMARY_NULL, axis, sets_def[name], name)
    M["set_tests_primary_null"] = res
    M["G1_pathway_on_goal_axis"] = res[f"{PRIMARY}@goal@{PRIMARY_NULL}"]["pass"]

    # 副次: 対照の部分集合（案 B 相当・案 C 相当）で判定が変わるか
    sec = {}
    for ver in a["secondary_null_versions"]:
        for name in [PRIMARY] + a["secondary_sets"] + ["myof_markers"]:
            sec[f"{name}@goal@{ver}"] = set_test(ver, "goal", sets_def[name], name)
        sec[f"myof_markers@myofaxis@{ver}"] = set_test(ver, "myofaxis", sets_def["myof_markers"], "myof_markers")
    M["set_tests_secondary_null"] = sec
    M["null_version_agreement"] = {
        "主判定（980 版）": M["G1_pathway_on_goal_axis"],
        "600 版": sec[f"{PRIMARY}@goal@B_600"]["pass"],
        "520 版": sec[f"{PRIMARY}@goal@C_520"]["pass"],
        "3 版で判定が一致": bool(M["G1_pathway_on_goal_axis"]
                                 == sec[f"{PRIMARY}@goal@B_600"]["pass"]
                                 == sec[f"{PRIMARY}@goal@C_520"]["pass"])}

    # 副次的な予測（循環でない証拠）: 非マーカー部分セットと、目標マーカー条件を抜いた版
    M["G_secondary_subsets"] = {name: res[f"{name}@goal@{PRIMARY_NULL}"]["pass"]
                                for name in a["secondary_sets"]}
    M["F3_nonmarker_pass"] = res[f"nonmarker@goal@{PRIMARY_NULL}"]["pass"]

    # G-2 交差（特異性）
    gg = np.array(res[f"{PRIMARY}@goal@{PRIMARY_NULL}"]["set_mean_z"])
    gm = np.array(res[f"{PRIMARY}@myofaxis@{PRIMARY_NULL}"]["set_mean_z"])
    mm = np.array(res[f"myof_markers@myofaxis@{PRIMARY_NULL}"]["set_mean_z"])
    mg = np.array(res[f"myof_markers@goal@{PRIMARY_NULL}"]["set_mean_z"])
    M["G2_crossover"] = {"経路セット: 目標の軸 − 筋線維芽細胞の軸": [round(float(x), 3) for x in gg - gm],
                         "筋線維芽細胞マーカー: 筋線維芽細胞の軸 − 目標の軸": [round(float(x), 3) for x in mm - mg],
                         "pass": bool(all(gg - gm > 0) and all(mm - mg > 0))}

    # G-3 アッセイの妥当性（F-6）: 筋線維芽細胞マーカー 24 が自分の軸で合格する
    M["G3_myof_markers_on_myof_axis"] = res[f"myof_markers@myofaxis@{PRIMARY_NULL}"]["pass"]

    # G-4 中立性（F-8）
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
    M["G4_neutral"] = neu
    M["G4_pass"] = all(v["n_seeds_in_range"] >= cr["neutral_min_seeds"] for v in neu.values())

    # G-5 偶然合格率の実測（guardrails 2026-09-28）
    k = cr["K_set"]
    rep = cr["set_draw_replace"]
    thr = [float(np.quantile(null_dist(PRIMARY_NULL, "goal", k, rep)[:, i], cr["set_percentile"]))
           for i in range(len(seeds))]
    nd = set_null(znull[(PRIMARY_NULL, "goal")], k, cr["n_chance_draw"], rng, rep)  # 閾値とは別に引き直す
    hit = np.ones(len(nd), bool)
    for i in range(len(seeds)):
        hit &= (nd[:, i] > thr[i]) & (nd[:, i] >= cr["set_abs_min"])
    M["G5_chance_pass_rate_random_sets"] = {"n_draw": int(len(nd)), "rate": float(hit.mean()),
                                            "threshold": [round(x, 3) for x in thr]}

    # G-6 偽経路セット（F-5）: 50 組以上。bias_goal もそろえる
    tgt = P[(P.seed == seeds[0]) & (P.gene.isin(sets_def[PRIMARY]))][
        ["gene", "lp_med", "ll_med", "n_cells"]].copy()
    tgt["symbol"] = tgt.gene.map(sym)
    tgt["bias_goal"] = tgt.symbol.map(bias)
    npl = P[(P.seed == seeds[0]) & (P.gene.isin(set(null_versions[PRIMARY_NULL])))][
        ["gene", "lp_med", "ll_med", "n_cells"]].copy()
    npl["symbol"] = npl.gene.map(sym)
    npl["bias_goal"] = npl.symbol.map(bias)
    def run_sham(target: pd.DataFrame, kk: int, tag: str, seed_off: int) -> tuple[dict, list]:
        shams = sham_sets(npl, target, kk, cr["n_sham_sets"],
                          np.random.default_rng(cr["null_draw_seed"] + seed_off),
                          cr.get("sham_match_bias", True))
        rows = []
        for i, names in enumerate(shams):
            r = set_test(PRIMARY_NULL, "goal", [ens[x] for x in names], f"{tag}_{i + 1}")
            r["genes"] = names
            rows.append(r)
        npass = sum(r["pass"] for r in rows)
        idx = npl.set_index("symbol")
        match = {c: round(float(np.mean([idx.loc[r["genes"], c].mean() for r in rows])
                                - target[c].mean()), 4)
                 for c in ("lp_med", "ll_med", "n_cells", "bias_goal")}
        return ({"n_sets": len(rows), "n_pass": int(npass), "pass_rate": round(npass / len(rows), 4),
                 "set_mean_z_range": [round(float(min(min(r["set_mean_z"]) for r in rows)), 3),
                                      round(float(max(max(r["set_mean_z"]) for r in rows)), 3)],
                 "そろい具合（50 組平均 − 対象セット）": match,
                 "注": "対照遺伝子は |bias_goal| < 0.5 に限られるので、経路セットの偏りの大きい側は"
                       "そろえきれない。bias_goal のずれが負に大きいほど、この検査は甘くなる。"},
                rows)

    M["G6_sham_sets"], sh = run_sham(tgt, k, "sham", 1)
    M["G6_pass"] = bool(M["G6_sham_sets"]["n_pass"] == 0)
    # 非マーカー部分セット（|bias_goal| < 0.5）は対照と同じ範囲にあるので、偏りまで完全にそろえられる
    tgt_nm = tgt[tgt.gene.isin(sets_def["nonmarker"])]
    M["G6_sham_sets_nonmarker"], sh_nm = run_sham(tgt_nm, len(sets_def["nonmarker"]), "sham_nonmarker", 2)
    M["G6_nonmarker_pass"] = bool(M["G6_sham_sets_nonmarker"]["n_pass"] == 0)
    pd.DataFrame([{"set": r["label"], "z": r["set_mean_z"], "pass": r["pass"],
                   "genes": " ".join(r["genes"])} for r in sh + sh_nm]).to_csv(
        E / "results/eval/set_null_002.csv", index=False)

    # --- 記述値（合否に入れない）
    zg = Zt[(PRIMARY_NULL, "goal")]
    zval = lambda g, s: float(zg[(zg.gene == g) & (zg.seed == s)].z.iloc[0])
    M["descriptive"] = {
        "単一遺伝子で 3 seed すべて 2.5 標準偏差以上": {
            name: sorted(sym[g] for g in genes if all(zval(g, s) >= 2.5 for s in seeds))
            for name, genes in sets_def.items()},
        "単一遺伝子の偶然合格率（対照 980 に同じ規則）": float(
            (znull[(PRIMARY_NULL, "goal")] >= 2.5).all(axis=1).mean()),
        "偏り（bias_goal）と z の順位相関（経路セット 120、目標の軸）": {
            str(s): float(stats.spearmanr([bias[sym[g]] for g in sets_def[PRIMARY]],
                                          [zval(g, s) for g in sets_def[PRIMARY]]).statistic)
            for s in seeds},
        "参照 ACTA2 / TAGLN（筋線維芽細胞の軸、標準偏差単位）": {
            t: [round(float(Zt[(PRIMARY_NULL, "myofaxis")][
                (Zt[(PRIMARY_NULL, "myofaxis")].symbol == t)
                & (Zt[(PRIMARY_NULL, "myofaxis")].seed == sd)].z.iloc[0]), 2) for sd in seeds]
            for t in ("ACTA2", "TAGLN")},
        "機体（seed ごと）": {str(s): "results/env-remote-seed%d.json を参照" % s for s in seeds},
    }

    M["primary_conclusion_proposal_inputs"] = {
        "metric_pass": M["metric_pass"], "G1": M["G1_pathway_on_goal_axis"],
        "G2": M["G2_crossover"]["pass"], "G3": M["G3_myof_markers_on_myof_axis"],
        "G4": M["G4_pass"], "G5_rate": M["G5_chance_pass_rate_random_sets"]["rate"],
        "G6": M["G6_pass"], "F3_nonmarker": M["F3_nonmarker_pass"]}
    ok = bool(M["metric_pass"] and M["G1_pathway_on_goal_axis"] and M["G3_myof_markers_on_myof_axis"]
              and M["G4_pass"] and M["G6_pass"]
              and M["G5_chance_pass_rate_random_sets"]["rate"] <= cr["max_chance_rate"])
    M["conclusion_proposal"] = (
        "判定不能（アッセイの妥当性が不合格。筋線維芽細胞マーカーが自分の軸で合格しない → F-6）"
        if not M["G3_myof_markers_on_myof_axis"] else
        "判定不能（指標の健全性が不合格 → F-7）" if not M["metric_pass"] else
        "判定不能（偽経路セットが合格した → F-5）" if not M["G6_pass"] else
        "判定不能（中立のはずのセットが動いた → F-8）" if not M["G4_pass"] else
        ("主予測は支持される。**かつ非マーカー部分セットも合格**するので、"
         "「目標の状態のマーカーだから動いた」という循環では説明できない"
         if ok and M["F3_nonmarker_pass"] else
         "主予測は支持されるが、**非マーカー部分セットは不合格**。"
         "循環の寄与を否定できない（F-3）。「経路が効く」とは書かない" if ok else
         "主予測は棄却。ただし「効果がない」ではなく「検出下限以上の効果は無い」としか読めない（F-1）"))
    M["note"] = (
        "判定の確定は PI（G3）。本スクリプトは事前登録の規則を当てはめて提案するだけ。"
        "本実験は LTBR / RELB / NFKB2 を回していないので、HYP-20260925-02 の判定には使えない。"
        "検定の単位はセットであり、個々の遺伝子が効くことは示さない。"
        "in silico の結果であり、実験的に確認された生物学的事実ではない。")

    (E / "results/eval/criteria_002.json").write_text(
        json.dumps(M, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps({k: v for k, v in M.items()
                      if k.endswith("_pass") or k.startswith("G") or k.startswith("F")
                      or k in ("null_version_agreement", "conclusion_proposal",
                               "missing_q_by_seed", "missing_in_pathway_set")},
                     ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    import sys
    if "--self-test" in sys.argv:
        self_test()
    else:
        main()

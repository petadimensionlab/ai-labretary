"""EXP-20260927-004 — 事前登録した判定の機械的な当てはめ。

EXP-20260927-003 の `analysis/evaluate_fib.py` から引き継がないこと（validation.md §11）:
  副検定の 状態の特異性 (d) が、区間内パーセンタイル（±0.05 の桁）を主検定の `D_raw_myof`
  （コサインの差、10⁻⁵ の桁）と比べていた。**本スクリプトは、比べる 2 つの量が同じ指標・同じ尺度から
  来ていることを `compare_same_scale()` で必ず検査してから比較する。**

判定対象にしないもの（読み込み直後に除外し、assert で確認する）:
  主対象 LTBR / RELB / NFKB2、参考対象 CDKN1A、盲検の TP53 / SPP1 / LCN2。
  本実験はこれらを摂動しない（plan.md §4）。

出力: results/eval/criteria_004.json, results/eval/gene_seed_004.csv
実行: .venv/bin/python experiments/EXP-20260927-004/analysis/evaluate_004.py
"""
from __future__ import annotations

import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260927-004"

BLOCKED = {"ENSG00000111321": "LTBR", "ENSG00000104856": "RELB", "ENSG00000077150": "NFKB2",
           "ENSG00000124762": "CDKN1A", "ENSG00000141510": "TP53",
           "ENSG00000118785": "SPP1", "ENSG00000148346": "LCN2"}


# ---------------------------------------------------------------- 指標
def cos(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return (a @ b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b) + 1e-12)


def metric_columns(c: pd.DataFrame, cfg: dict) -> dict[str, str]:
    """事前登録した指標の列を作る。返り値は 名前 → 列名。
    primary / placebo / myof_mirror は**同じ関数形・同じ尺度**で、重心の役割だけを入れ替えたもの。"""
    g, st, al, my, dn = c.s_goal, c.s_start, c.s_alt, c.s_myof, c.dnorm + 1e-30
    norm = cfg["analysis"]["primary_metric"] == "A12"
    def form(target, others):
        v = target - sum(others) / len(others)
        return v / dn if norm else v
    c["m_primary"] = form(g, [st, al, my])        # 目標に特異的な成分
    c["m_placebo"] = form(al, [st, g, my])        # 目標 ↔ 対照の状態 を入れ替えた鏡像
    c["m_myof"] = form(my, [st, al, g])           # 目標 ↔ 筋線維芽細胞 を入れ替えた鏡像
    c["m_uniform"] = (g + st + al + my) / 4.0     # どの重心へも一様に近づく成分（診断用）
    c["m_A11"] = g - (st + al + my) / 3.0         # 副指標（移動量で割らない形）
    return {"primary": "m_primary", "placebo": "m_placebo", "myof_mirror": "m_myof",
            "uniform": "m_uniform", "A11": "m_A11"}


SCALE_GROUP = {"m_primary": "specific", "m_placebo": "specific", "m_myof": "specific",
               "m_uniform": "raw_shift", "m_A11": "specific_unnormalized"}


def compare_same_scale(col_a: str, col_b: str) -> None:
    """**EXP-20260927-003 の不具合の再発防止**: 桁の違う 2 つの量を比べていないことを検査する。"""
    ga, gb = SCALE_GROUP.get(col_a), SCALE_GROUP.get(col_b)
    assert ga is not None and gb is not None, f"尺度の登録がない列を比較している: {col_a}, {col_b}"
    assert ga == gb, (f"尺度の違う量を比較しようとしている: {col_a}（{ga}）と {col_b}（{gb}）。"
                      f"EXP-20260927-003 validation.md §11 と同じ不具合。")


# ---------------------------------------------------------------- 主検定（EXP-003 と同じ実装）
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
    res = E / cfg["output"]["results_dir"]
    c = pd.read_csv(res / "cells.csv", dtype={"donor_id": str})
    n0 = len(c)
    c = c[~c.gene.isin(BLOCKED)].reset_index(drop=True)
    assert not set(c.gene) & set(BLOCKED), "回さないはずの遺伝子が結果に入っている"
    if n0 != len(c):
        raise SystemExit(f"回さないはずの遺伝子が {n0 - len(c)} 行あった。実行設定を確認すること")
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
    for (run, gene), part in c.groupby(["run", "gene"], sort=False):
        seed = int(str(run).split("-")[1])
        sz = np.load(res / run / "state_emb_vectors.npz")
        xz = np.load(E / f"work/extra_centroids-seed{seed}.npz")
        cen = {"goal": sz["FIB_org"], "start": sz["FIB_neg_TLShigh"],
               "alt": sz["FIB_neg_TLSlow"], "myof": xz["MYOF_neg"]}
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


# ---------------------------------------------------------------- 判定
def main() -> None:
    cfg = yaml.safe_load(open(E / "config.yaml"))
    a = cfg["analysis"]
    cr = a["criteria"]
    roles = {g["symbol"]: g["role"] for g in json.loads((E / "work/genes_004.json").read_text())["genes"]}
    ens = {g["symbol"]: g["ensembl"] for g in json.loads((E / "work/genes_004.json").read_text())["genes"]}
    sym = {v: k for k, v in ens.items()}
    nulls_set = {ens[s] for s, r in roles.items() if r == "negative"}
    pos_goal = [ens[s] for s, r in roles.items() if r == "positive_goal"]
    hk = [ens[s] for s, r in roles.items() if r == "housekeeping"]
    panel = [ens[s] for s, r in roles.items() if r == "negctrl_panel"]
    pos_myof = [ens[s] for s, r in roles.items() if r == "positive_myof"]

    # 本スクリプトは M5 を「3 seed すべて」で実装している。config がそれと食い違っていたら止める
    assert cr["require_all_seeds"] is True, \
        "config の require_all_seeds が False。M5 は 3 seed すべてで実装している（plan.md §6.4）"

    c = load(cfg)
    cols = metric_columns(c, cfg)
    seeds = sorted(c.seed.unique().tolist())
    assert len(seeds) == cfg["robustness"]["seeds"].__len__(), "seed の数が計画と違う"

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
    (E / "results/eval").mkdir(parents=True, exist_ok=True)
    long.to_csv(E / "results/eval/gene_seed_004.csv", index=False)

    P = tabs["primary"]
    NUL = P[P.gene.isin(nulls_set)]
    sd = {s: float(NUL[NUL.seed == s].D_raw.std(ddof=1)) for s in seeds}
    mu = {s: float(NUL[NUL.seed == s].D_raw.mean()) for s in seeds}

    def z(gene: str, s: int) -> float:
        r = P[(P.gene == gene) & (P.seed == s)]
        return float((r.D_raw.iloc[0] - mu[s]) / sd[s]) if len(r) else np.nan

    def qv(tab: pd.DataFrame, gene: str, s: int) -> float:
        r = tab[(tab.gene == gene) & (tab.seed == s)]
        return float(r.q.iloc[0]) if len(r) else np.nan

    def dv(tab: pd.DataFrame, gene: str, s: int) -> float:
        r = tab[(tab.gene == gene) & (tab.seed == s)]
        return float(r.D_raw.iloc[0]) if len(r) else np.nan

    # --- 指標の合格条件（M1〜M5）
    M = {}
    M["M1_calibration"] = {}
    for s in seeds:
        q = NUL[NUL.seed == s].q.dropna()
        M["M1_calibration"][s] = {"n_valid": int(len(q)), "frac_hi": float((q >= 0.90).mean()),
                                  "frac_lo": float((q < 0.10).mean())}
    M["M1_pass"] = all(v["n_valid"] >= cr["M1_min_n"] and v["frac_hi"] <= cr["M1_tail_frac_max"]
                       and v["frac_lo"] <= cr["M1_tail_frac_max"] for v in M["M1_calibration"].values())

    M["M2_confounding"] = {}
    for s in seeds:
        n = NUL[NUL.seed == s]
        M["M2_confounding"][s] = {k: float(stats.spearmanr(n[f"{k}_med"], n.D_raw).statistic)
                                  for k in ("lp", "ll", "relpos")}
    M["M2_pass"] = sum(all(abs(v) <= cr["M2_abs_rho_max"] for v in d.values())
                       for d in M["M2_confounding"].values()) >= cr["M2_min_seeds"]

    M["M3_housekeeping"] = {sym[g]: [qv(P, g, s) for s in seeds] for g in hk}
    M["M3_panel"] = {sym[g]: [qv(P, g, s) for s in seeds] for g in panel}
    lo, hi = cr["M3_range"]
    inrange = lambda qs: sum(1 for x in qs if not np.isnan(x) and lo <= x <= hi) >= cr["M3_min_seeds"]
    M["M3_pass"] = (all(inrange(v) for v in M["M3_housekeeping"].values())
                    and sum(inrange(v) for v in M["M3_panel"].values()) >= cr["M3_panel_min_genes"])

    # M4: 方向を測れていること（一様成分・置き換え指標との相関）
    M["M4_rho"] = {}
    for s in seeds:
        n = NUL[NUL.seed == s].set_index("gene")
        u = tabs["uniform"][tabs["uniform"].seed == s].set_index("gene").reindex(n.index)
        pl = tabs["placebo"][tabs["placebo"].seed == s].set_index("gene").reindex(n.index)
        M["M4_rho"][s] = {"uniform": float(stats.spearmanr(n.D_raw, u.D_raw).statistic),
                          "placebo": float(stats.spearmanr(n.D_raw, pl.D_raw).statistic)}
    M["M4_pass"] = all(abs(v["uniform"]) <= cr["M4_abs_rho_uniform_max"]
                       and abs(v["placebo"]) <= cr["M4_abs_rho_placebo_max"]
                       for v in M["M4_rho"].values())

    # M5: 目標の軸での陽性対照（最重要）
    compare_same_scale(cols["primary"], cols["placebo"])
    compare_same_scale(cols["primary"], cols["myof_mirror"])
    pc = {}
    for g in pos_goal:
        qs = [qv(P, g, s) for s in seeds]
        ds_ = [dv(P, g, s) for s in seeds]
        zs = [z(g, s) for s in seeds]
        spec = [bool(dv(P, g, s) > dv(tabs["placebo"], g, s)
                     and dv(P, g, s) > dv(tabs["myof_mirror"], g, s)) for s in seeds]
        passed = all((not np.isnan(q)) and q >= cr["q_hi"] and d > 0 and zz >= cr["min_sd_units"]
                     for q, d, zz in zip(qs, ds_, zs))
        pc[sym[g]] = {"q": qs, "D_raw": ds_, "sd_units": [round(v, 2) for v in zs],
                      "state_specific": spec, "pass": bool(passed)}
    M["M5_positive_goal"] = pc
    M["M5_n_pass"] = int(sum(v["pass"] for v in pc.values()))
    M["M5_pass"] = M["M5_n_pass"] >= cr["M5_min_genes"]

    # 参考: 筋線維芽細胞の軸での陽性対照（パイプラインが動いていることの確認。判定には使わない）
    ref = {}
    X = tabs["myof_mirror"]
    XN = X[X.gene.isin(nulls_set)]
    for g in pos_myof:
        zs = []
        for s in seeds:
            n = XN[XN.seed == s]
            zs.append(round(float((dv(X, g, s) - n.D_raw.mean()) / n.D_raw.std(ddof=1)), 2))
        ref[sym[g]] = {"q": [qv(X, g, s) for s in seeds], "sd_units": zs}
    M["reference_myof_axis"] = ref

    # 偶然合格率の実測（対照遺伝子に M5 と同じ規則を当てはめる。seed 間の独立を仮定しない）
    hits = 0
    usable = 0
    for g in sorted(nulls_set):
        qs = [qv(P, g, s) for s in seeds]
        if sum(not np.isnan(x) for x in qs) < len(seeds):
            continue
        usable += 1
        hits += int(all(q >= cr["q_hi"] and dv(P, g, s) > 0 and z(g, s) >= cr["min_sd_units"]
                        for q, s in zip(qs, seeds)))
    M["chance_pass_rate_measured"] = {"n_null_genes": usable, "n_pass": hits,
                                      "rate": (hits / usable) if usable else None}
    M["metric_pass"] = bool(M["M1_pass"] and M["M2_pass"] and M["M3_pass"] and M["M4_pass"])
    M["conclusion_proposal"] = (
        "目標の軸で摂動を検出できる（M5 合格）" if M["metric_pass"] and M["M5_pass"] else
        "指標は健全だが、目標の軸で動かせる遺伝子を 1 つも検出できない" if M["metric_pass"] else
        "指標の合格条件を満たさない（判定不能）")
    M["note"] = ("判定の確定は PI（G3）。本スクリプトは事前登録の規則を当てはめて提案するだけ。"
                 "本実験は主対象 LTBR/RELB/NFKB2 を回していないので、HYP-20260925-02 の判定には使えない。")

    (E / "results/eval/criteria_004.json").write_text(
        json.dumps(M, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps({k: v for k, v in M.items() if k.endswith("_pass") or k in
                      ("M5_n_pass", "chance_pass_rate_measured", "conclusion_proposal")},
                     ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()

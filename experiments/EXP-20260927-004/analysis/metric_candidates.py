"""EXP-20260927-004 — 指標候補の事後評価（EXP-20260927-003 の既存データを使う）。

重要な性質:
  - **主対象 LTBR / RELB / NFKB2 と参考対象 CDKN1A は読み込み直後に除外する**（読み込み後に assert）。
    盲検の TP53 / SPP1 / LCN2 も同様。指標の選定は対照遺伝子だけで行う。
  - これは EXP-20260927-003 の結果を見た後の解析である。ここで選ぶのは「候補の絞り込み」までで、
    合否の判定は EXP-20260927-004 の新しい seed のデータに対して行う（plan.md §6）。

出力: work/metric_candidates.json, work/metric_candidates.csv

実行:
  .venv/bin/python experiments/EXP-20260927-004/analysis/metric_candidates.py
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
SRC = ROOT / "experiments/EXP-20260927-003"
OUT = ROOT / "experiments/EXP-20260927-004/work"

# 除外（主対象・参考対象・未解除の盲検）。指標の選定に一切使わない
EXCLUDE = {
    "ENSG00000111321": "LTBR", "ENSG00000104856": "RELB", "ENSG00000077150": "NFKB2",
    "ENSG00000124762": "CDKN1A",
    "ENSG00000141510": "TP53", "ENSG00000118785": "SPP1", "ENSG00000148346": "LCN2",
}
HK = {"ENSG00000075624": "ACTB", "ENSG00000111640": "GAPDH", "ENSG00000196262": "PPIA"}
POS_MYOF = {"ENSG00000107796": "ACTA2", "ENSG00000149591": "TAGLN"}

TOL, COVER, MIN_ELIG = 0.5, 0.8, 20


def cos(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return (a @ b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b) + 1e-12)


def load() -> pd.DataFrame:
    cfg = yaml.safe_load(open(SRC / "config.yaml"))
    res = SRC / "results"
    c = pd.read_csv(res / "cells.csv", dtype={"donor_id": str})
    n0 = len(c)
    c = c[~c.gene.isin(EXCLUDE)].reset_index(drop=True)
    assert not set(c.gene) & set(EXCLUDE), "除外した遺伝子が残っている"
    print(f"除外前 {n0} 行 → 除外後 {len(c)} 行、遺伝子 {c.gene.nunique()}")
    c["seed"] = c.run.astype(str).str.extract(r"run-(\d+)", expand=False).astype(int)

    # 元の位置（入力内の順位）
    from datasets import load_from_disk
    from geneformer import TOKEN_DICTIONARY_FILE
    tokd = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
    ds = load_from_disk(str(SRC / "work/tokenized.dataset"))
    uids = sorted(c.cell_uid.unique())
    ids = dict(zip(uids, ds.select(uids)["input_ids"]))
    c["pos"] = [list(ids[u]).index(tokd[g]) for u, g in zip(c.cell_uid, c.gene)]
    c["lp"] = np.log2(1 + c.pos.astype(float))
    c["ll"] = np.log2(c.length.astype(float))
    c["relpos"] = c.pos / c.length

    # 5 つの重心への接近度と、移動量 ||Δ||
    names = ["goal", "start", "alt", "myof", "goal_unmatched"]
    vals = {n: np.full(len(c), np.nan) for n in names}
    dnorm = np.full(len(c), np.nan)
    for (run, gene), part in c.groupby(["run", "gene"], sort=False):
        seed = int(str(run).split("-")[1])
        sz = np.load(res / run / "state_emb_vectors.npz")
        xz = np.load(SRC / f"work/extra_centroids-seed{seed}.npz")
        cen = {"goal": sz["FIB_org"], "start": sz["FIB_neg_TLShigh"], "alt": sz["FIB_neg_TLSlow"],
               "myof": xz["MYOF_neg"], "goal_unmatched": xz["FIB_org_unmatched"]}
        gd = res / run / "genes" / gene / "emb"
        o = np.load(gd / "orig.npy").astype(np.float64)
        p = np.load(gd / "pert.npy").astype(np.float64)
        assert len(o) == len(part), f"{gd}: 埋め込み {len(o)} ≠ cells.csv {len(part)}"
        idx = part.index.to_numpy()
        for n in names:
            vals[n][idx] = cos(p, cen[n].astype(np.float64)) - cos(o, cen[n].astype(np.float64))
        dnorm[idx] = np.linalg.norm(p - o, axis=1)
        # 方向の射影: Δ の単位ベクトルと（重心A − 重心B）の単位ベクトルの内積
        for tag, (A, B) in {"proj_goal_start": ("goal", "start"),
                            "proj_goal_myof": ("goal", "myof"),
                            "proj_myof_start": ("myof", "start"),
                            "proj_alt_start": ("alt", "start"),
                            "proj_alt_myof": ("alt", "myof"),
                            "proj_goal_alt": ("goal", "alt"),
                            "proj_myof_alt": ("myof", "alt")}.items():
            ax = cen[A].astype(np.float64) / np.linalg.norm(cen[A]) - \
                 cen[B].astype(np.float64) / np.linalg.norm(cen[B])
            ax = ax / np.linalg.norm(ax)
            d = p - o
            v = (d @ ax) / (np.linalg.norm(d, axis=1) + 1e-30)
            c.loc[idx, tag] = v
    for n in names:
        c[f"s_{n}"] = vals[n]
    c["dnorm"] = dnorm
    return c


def add_candidates(c: pd.DataFrame) -> dict[str, str]:
    """候補の指標（細胞ごとの値）を列として足す。返り値は 候補名 → 説明。"""
    g, st, al, my = c.s_goal, c.s_start, c.s_alt, c.s_myof
    U = (g + st + al + my) / 4.0                      # 一様な移動量（4 重心の平均接近度）
    c["m_A0"] = g
    c["m_A1"] = g - st
    c["m_A2"] = g - my
    c["m_A3"] = g - al
    c["m_A5"] = g - (st + my) / 2.0
    c["m_A6"] = c.proj_goal_start
    c["m_A7"] = c.proj_goal_myof
    c["m_A9"] = (g - my) / (c.dnorm + 1e-30)
    c["m_A10"] = (g - st) / (c.dnorm + 1e-30)
    c["m_A11"] = g - (st + al + my) / 3.0
    c["m_A12"] = (g - (st + al + my) / 3.0) / (c.dnorm + 1e-30)
    c["m_A13"] = c.proj_goal_alt
    c["U"] = U
    # 置き換え検定用（goal を alt に差し替えた同じ形の指標）
    c["p_A0"] = al
    c["p_A1"] = al - st
    c["p_A2"] = al - my
    c["p_A3"] = g - al          # A3 の置き換えは定義できない（自己参照）→ NaN 扱い
    c["p_A5"] = al - (st + my) / 2.0
    c["p_A6"] = c.proj_alt_start
    c["p_A7"] = c.proj_alt_myof
    c["p_A9"] = (al - my) / (c.dnorm + 1e-30)
    c["p_A10"] = (al - st) / (c.dnorm + 1e-30)
    c["p_A11"] = al - (st + g + my) / 3.0
    c["p_A12"] = (al - (st + g + my) / 3.0) / (c.dnorm + 1e-30)
    c["p_A13"] = -c.proj_goal_alt
    # 筋線維芽細胞の軸に写した鏡像（ACTA2 / TAGLN の陽性対照の確認用。goal ↔ myof を入れ替える）
    c["x_A0"] = my
    c["x_A1"] = my - st
    c["x_A2"] = my - g
    c["x_A3"] = my - al
    c["x_A5"] = my - (st + g) / 2.0
    c["x_A6"] = c.proj_myof_start
    c["x_A7"] = -c.proj_goal_myof
    c["x_A9"] = (my - g) / (c.dnorm + 1e-30)
    c["x_A10"] = (my - st) / (c.dnorm + 1e-30)
    c["x_A11"] = my - (st + al + g) / 3.0
    c["x_A12"] = (my - (st + al + g) / 3.0) / (c.dnorm + 1e-30)
    c["x_A13"] = c.proj_myof_alt
    return {
        "A0": "目標への接近（EXP-003 の主指標）",
        "A1": "目標への接近 − 出発点への接近",
        "A2": "目標への接近 − 筋線維芽細胞への接近",
        "A3": "目標への接近 − 対照の状態（TLS 少）への接近",
        "A5": "目標への接近 − (出発点・筋線維芽細胞の平均)",
        "A6": "移動方向の射影（目標 − 出発点の軸, 単位ベクトル）",
        "A7": "移動方向の射影（目標 − 筋線維芽細胞の軸, 単位ベクトル）",
        "A9": "(目標 − 筋線維芽細胞) / 移動量 ||Δ||",
        "A10": "(目標 − 出発点) / 移動量 ||Δ||",
        "A11": "目標への接近 − (出発点・対照の状態・筋線維芽細胞の平均)",
        "A12": "A11 / 移動量 ||Δ||",
        "A13": "移動方向の射影（目標 − 対照の状態の軸, 単位ベクトル）",
    }


def design(lp, ll):
    lp, ll = np.asarray(lp, float), np.asarray(ll, float)
    return np.column_stack([np.ones_like(lp), lp, ll, lp ** 2, ll ** 2, lp * ll])


def eligibility(cs: pd.DataFrame, nulls: dict[str, pd.DataFrame]) -> dict[str, list[str]]:
    """遺伝子 → 使える対照遺伝子のリスト（y に依存しないので 1 度だけ計算する）。"""
    out = {}
    for gene, tg in cs.groupby("gene"):
        lp, ll = tg.lp.to_numpy(), tg.ll.to_numpy()
        elig = []
        for h, hd in nulls.items():
            if h == gene:
                continue
            D = np.maximum(np.abs(lp[:, None] - hd.lp.to_numpy()[None, :]),
                           np.abs(ll[:, None] - hd.ll.to_numpy()[None, :]))
            if (D <= TOL).any(axis=1).mean() >= COVER:
                elig.append(h)
        out[gene] = elig
    return out


def gene_stats(cs: pd.DataFrame, nulls: dict[str, pd.DataFrame], elig: dict, col: str) -> pd.DataFrame:
    """偏り補正 → 遺伝子平均 D_raw → 使える対照遺伝子の中での順位 q。"""
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
        q = float((Dh < Dg).mean()) if len(Dh) >= MIN_ELIG else np.nan
        rows.append({"gene": gene, "D_raw": Dg, "q": q, "n_elig": len(Dh),
                     "T": float(tg[col].mean()), "n_cells": len(tg),
                     "lp_med": float(tg.lp.median()), "ll_med": float(tg.ll.median()),
                     "relpos_med": float(tg.relpos.median())})
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / "cells_with_shifts.parquet"
    if cache.exists():
        c = pd.read_parquet(cache)
        assert not set(c.gene) & set(EXCLUDE), "除外した遺伝子が残っている"
        print(f"キャッシュを読み込み: {len(c)} 行、遺伝子 {c.gene.nunique()}")
    else:
        c = load()
        c.to_parquet(cache)
    desc = add_candidates(c)
    null_genes = set(yaml.safe_load(open(SRC / "config.yaml"))["controls"]["negative_sampler"]["genes"])
    panel = set(yaml.safe_load(open(SRC / "config.yaml"))["controls"]["extra_fixed"][1]["genes"])
    seeds = sorted(c.seed.unique().tolist())

    per_seed_elig, per_seed_nulls = {}, {}
    for s in seeds:
        cs = c[c.seed == s]
        nulls = {h: hd for h, hd in cs[cs.gene.isin(null_genes)].groupby("gene")}
        per_seed_nulls[s] = nulls
        per_seed_elig[s] = eligibility(cs, nulls)
        print(f"seed {s}: 対照遺伝子 {len(nulls)}、遺伝子 {cs.gene.nunique()}")

    frames = {}
    for tag in list(desc) + ["U"]:
        cols = {"m": f"m_{tag}" if tag != "U" else "U"}
        if tag != "U":
            cols["p"] = f"p_{tag}"
            cols["x"] = f"x_{tag}"
        for kind, col in cols.items():
            out = []
            for s in seeds:
                d = gene_stats(c[c.seed == s], per_seed_nulls[s], per_seed_elig[s], col)
                d["seed"] = s
                out.append(d)
            frames[(tag, kind)] = pd.concat(out, ignore_index=True)
        print(f"候補 {tag} 完了")

    # ---------------- 診断
    rep = {}
    rows = []
    for tag in desc:
        M, P, X = frames[(tag, "m")], frames[(tag, "p")], frames[(tag, "x")]
        Uf = frames[("U", "m")]
        rho_U, rho_P, cal_hi, cal_lo, sds = [], [], [], [], []
        conf = {"lp": [], "relpos": [], "ll": []}
        for s in seeds:
            m = M[(M.seed == s) & M.gene.isin(null_genes)].set_index("gene")
            u = Uf[(Uf.seed == s)].set_index("gene").reindex(m.index)
            p = P[(P.seed == s)].set_index("gene").reindex(m.index)
            rho_U.append(float(stats.spearmanr(m.D_raw, u.D_raw).statistic))
            if tag not in ("A3", "A13"):
                rho_P.append(float(stats.spearmanr(m.D_raw, p.D_raw).statistic))
            q = m.q.dropna()
            cal_hi.append(float((q >= 0.90).mean()))
            cal_lo.append(float((q < 0.10).mean()))
            sds.append(float(m.D_raw.std(ddof=1)))
            for k in conf:
                conf[k].append(float(stats.spearmanr(m[f"{k}_med"], m.D_raw).statistic))
        # 分散分解（対照遺伝子のみ）: 遺伝子固有で再現する成分 vs seed 間の誤差
        w = M[M.gene.isin(null_genes)].pivot(index="gene", columns="seed", values="D_raw").dropna()
        v_between = float(w.mean(axis=1).var(ddof=1))
        v_within = float(w.var(axis=1, ddof=1).mean())            # seed 間の誤差の分散
        v_gene = max(v_between - v_within / len(seeds), 0.0)
        # 恒常発現遺伝子・中立対照パネルの中立性
        neu = {}
        for gid, nm in HK.items():
            qq = M[M.gene == gid].set_index("seed").q
            neu[nm] = [None if np.isnan(v) else round(float(v), 3) for v in qq.reindex(seeds)]
        panel_ok = 0
        for gid in panel:
            qq = M[M.gene == gid].set_index("seed").q.reindex(seeds)
            panel_ok += int((( qq >= 0.10) & (qq <= 0.90)).sum() >= 2)
        # 筋線維芽細胞の軸での陽性対照（鏡像の指標で ACTA2 / TAGLN が上位に来るか）
        pos = {}
        for gid, nm in POS_MYOF.items():
            qq = X[X.gene == gid].set_index("seed").q.reindex(seeds)
            dd = X[X.gene == gid].set_index("seed").D_raw.reindex(seeds)
            pos[nm] = {"q": [None if np.isnan(v) else round(float(v), 3) for v in qq],
                       "n_seeds_pass": int(((qq >= 0.90) & (dd > 0)).sum()),
                       "sd_units": [round(float((dd[s] - X[(X.seed == s) & X.gene.isin(null_genes)].D_raw.mean())
                                                / X[(X.seed == s) & X.gene.isin(null_genes)].D_raw.std(ddof=1)), 2)
                                    for s in seeds]}
        sd_null = float(np.mean(sds))
        sd_pooled = float(np.sqrt(v_gene + v_within / len(seeds)))
        r = {"desc": desc[tag],
             "rho_uniform_by_seed": [round(v, 3) for v in rho_U],
             "rho_placebo_alt_by_seed": [round(v, 3) for v in rho_P] if rho_P else None,
             "calibration_frac_hi": [round(v, 3) for v in cal_hi],
             "calibration_frac_lo": [round(v, 3) for v in cal_lo],
             "confounding_rho": {k: [round(v, 3) for v in vv] for k, vv in conf.items()},
             "null_sd_per_seed": sd_null,
             "var_gene_specific_frac": round(v_gene / (v_gene + v_within) if v_gene + v_within > 0 else np.nan, 3),
             "sd_seed_error": float(np.sqrt(v_within)),
             "sd_pooled_3seed": sd_pooled,
             "housekeeping_q": neu,
             "panel_n_neutral_of_7": panel_ok,
             "myof_axis_positive_control": pos}
        rep[tag] = r
        rows.append({"cand": tag, "desc": desc[tag],
                     "rho_U_mean": round(float(np.mean(np.abs(rho_U))), 3),
                     "rho_placebo_mean": (round(float(np.mean(np.abs(rho_P))), 3) if rho_P else None),
                     "cal_hi": round(float(np.mean(cal_hi)), 3), "cal_lo": round(float(np.mean(cal_lo)), 3),
                     "conf_max": round(max(max(np.abs(v)) for v in conf.values()), 3),
                     "null_sd": sd_null, "gene_var_frac": r["var_gene_specific_frac"],
                     "sd_pooled": sd_pooled,
                     "GAPDH_q": neu["GAPDH"], "ACTB_q": neu["ACTB"], "PPIA_q": neu["PPIA"],
                     "panel_ok": panel_ok,
                     "ACTA2_myof_pass": pos["ACTA2"]["n_seeds_pass"],
                     "ACTA2_myof_sd": pos["ACTA2"]["sd_units"],
                     "TAGLN_myof_pass": pos["TAGLN"]["n_seeds_pass"]})

    pd.DataFrame(rows).to_csv(OUT / "metric_candidates.csv", index=False)
    (OUT / "metric_candidates.json").write_text(
        json.dumps({"note": "EXP-20260927-003 の既存データ（84 遺伝子 × 3 seed）から、"
                            "主対象 LTBR/RELB/NFKB2 と参考対象 CDKN1A を除いた 80 遺伝子で計算。"
                            "結果を見た後の解析であり、候補の絞り込みにのみ使う。",
                    "seeds": seeds, "n_genes": int(c.gene.nunique()),
                    "excluded": list(EXCLUDE.values()), "candidates": rep},
                   ensure_ascii=False, indent=1, default=float) + "\n")
    pd.set_option("display.width", 220, "display.max_columns", 40)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()

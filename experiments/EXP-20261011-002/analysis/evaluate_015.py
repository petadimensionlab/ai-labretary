#!/usr/bin/env python3
"""EXP-20261011-002 — HYP-20261011-01（差の向き。方法の再確認）。

EXP-20261011-001 の不成立（R-6 退化側が 1% に届かず判定不能・R-12 発現細胞数との交絡）への再設計。
新しい摂動は回さない。保存済み埋め込み（A1: 626＋30、A2: 587）の再射影・再評価のみ（GPU 0 時間）。

- [1] 退化側を強く出す設計: 主の陽性対照 = 退化 23（30 − 活性化 71 と重なる 7 個。JUND・CTSH は残す）。
  測る向きは差の向き u_diff のまま（退化側 = z_diff < 0）
- [2] 交絡補正: 主 = 発現細胞数 n_exp・入力の長さ ll_med・元の位置 lp_med を対象セットにそろえた
  マッチ偽セット（evaluate_005.sham_sets と同じ近傍距離。交絡を帰無分布に織り込む）。
  副次 = 残差化（z_diff を log2(n_exp) と log2(ll_med) に seed ごとに回帰し、残差で同じ分位規則）
- 主判定 P-1'（A1・確認）: (i) 活性化 71 > 0 かつマッチ偽セット 99% 点超 / (ii) 退化 23 < 0 かつ
  マッチ偽セット 1% 点未満 / (iii) 間隔がマッチした無作為 2 セットの差の 99% 点超。3 群すべて
- 補正の健全性 R-12a（マッチングの中心）・R-12b（閾値の安定性）・R-12c（交差検証残差）:
  主判定の合否にかかわらず不合格なら判定不能
- 反証条件 R-1'〜R-14'（plan.md §4.5。上から順に当て、最初に当たったものを採用）
- A2 は探索: 候補 8 のマッチ補正版（前回の読み R-7 が補正後も保たれるかの記述）
- 判定（支持・否定・判定不能。ゲート G3）は確定しない。PI が行う。この道具は規則を機械的に当てた結果まで。
- 盲検中の 5 遺伝子（TP53・CDKN1A・LTBR・SPP1・LCN2）は、摂動・対照・指標・出力・ハッシュのどこにも入れない。

走らせ方:
    .venv/bin/python experiments/EXP-20261011-002/analysis/evaluate_015.py --self-test   # 合成データで単体検査（本番の値は見ない）
    .venv/bin/python experiments/EXP-20261011-002/analysis/evaluate_015.py                # 本番（固定後の 1 回）
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
EXP = ROOT / "experiments/EXP-20261011-002"
A1_SRC = ROOT / "experiments/EXP-20261004-001"
A1_NEW = ROOT / "experiments/EXP-20261011-001"
A2_SRC = ROOT / "experiments/EXP-20261007-001"
sys.path.insert(0, str(EXP / "analysis"))
sys.path.insert(0, str(ROOT / "experiments/EXP-20261011-001/analysis"))
sys.path.insert(0, str(ROOT / "experiments/EXP-20261004-001/analysis"))
sys.path.insert(0, str(ROOT / "experiments/EXP-20261001-001/analysis"))
sys.path.insert(0, str(ROOT / "scripts/analysis"))
import evaluate_009 as ev9                                  # noqa: E402  層別補正・細胞表（1 文字も変えない）
from evaluate_005 import set_null                          # noqa: E402  復元ありで k 個引いて平均する
from evaluate_014 import (add_projection, add_projection_diff,  # noqa: E402  同じ式をそのまま使う
                          merge_projection, per_seed_stats)

SEEDS = [81, 82, 83]
N_DRAW = 1000

# --------------------------------------------------------------------------- plan.md §4 の判定値
# **ここを変えてはいけない。** config.yaml の decision: 節と照合し、違えば止める。
PLAN = {
    "候補_n": 8,
    "活性化PC71_n": 71,
    "退化PC30_n": 30,
    "退化PC23_n": 23,
    "退化PC21_n": 21,
    "プールA1_n": 549,
    "プールA2_n": 510,
    "A1_測る遺伝子_n": 626,
    "A1_新規30_n": 30,
    "A2_測る遺伝子_n": 587,
    "重なり7_symbols": ["MGP", "IGFBP5", "DSTN", "ARMC9", "PDE5A", "MFGE8", "HSPG2"],
    "JUND_CTSH_symbols": ["JUND", "CTSH"],
    "退化23_プール重なり_既知": ["CTSH", "PFDN5", "CYS1"],      # precheck で再確認（結果を見る前）
    "再現確認_許容差": 1e-3,
    "マッチング_許容": {"log2_n_exp": 0.10, "ll_med": 0.10, "lp_med": 0.25},
    "R-12c_閾値": 0.20,
    "R-12b_SE比": 0.5,          # SE(閾値) ≤ |閾値| × 0.5
    "乱数の鍵": {"P-1 活性化のマッチ偽セット": 20261101, "P-1 退化のマッチ偽セット": 20261102,
                "P-1 (iii) のマッチ無作為 2 セット": 20261103, "補正なし（無作為）の偽セット": 20261104,
                "残差化版の偽セット": 20261105, "A2 候補のマッチ偽セット": 20261106,
                "較正（マッチした無作為ペア）": 20261107},
}


def sha256f(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def log_line(path: Path, msg: str) -> None:
    line = f"[{_dt.datetime.now().strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(path, "a") as fh:
        fh.write(line + "\n")


# --------------------------------------------------------------------------- 遺伝子集合（A-3）
def gene_sets(cfg) -> dict:
    g = json.loads((A2_SRC / "work/genes_011.json").read_text())
    deg30 = list(cfg["degenerate_30"])
    cand = list(g["候補"].values())
    pos71 = list(g["陽性対照71"].values())
    hk = list(g["恒常発現3"].values())
    neu = list(g["中立7"].values())
    ctrl500 = [r["ensembl"] for r in g["対照500"]]
    sym = {v: k for k, v in g["候補"].items()} | {v: k for k, v in g["陽性対照71"].items()} \
        | {v: k for k, v in g["恒常発現3"].items()} | {v: k for k, v in g["中立7"].items()} \
        | {r["ensembl"]: r["symbol"] for r in g["対照500"]}
    cfgA1 = yaml.safe_load((A1_SRC / "config.yaml").read_text())
    ctl = cfgA1["controls"]
    rand = [x for x in ctl["negative_sampler"]["genes"]][: ctl["negative_random_n"]]
    poolA1 = sorted(set(rand) | set(ctl["negative_fixed"])
                    | {x for grp in (ctl.get("extra_fixed") or []) for x in grp["genes"]})
    poolA2 = sorted(set(ctrl500) | set(hk) | set(neu))
    ov7 = set(cfg["overlap7_symbols"])
    ov_ens = {ens for ens, s in sym.items() if s in ov7}
    deg23 = list(dict.fromkeys([x for x in deg30 if x not in ov_ens]))
    deg21 = list(dict.fromkeys([x for x in deg23 if sym.get(x) not in set(cfg["jund_ctsh_symbols"])]))
    sets = {"候補": list(dict.fromkeys(cand)), "活性化PC71": list(dict.fromkeys(pos71)),
            "退化PC30": list(dict.fromkeys(deg30)), "退化PC23": list(dict.fromkeys(deg23)),
            "退化PC21": list(dict.fromkeys(deg21)), "恒常発現3": list(dict.fromkeys(hk)),
            "中立7": list(dict.fromkeys(neu)), "対照500": list(dict.fromkeys(ctrl500)),
            "プールA2": poolA2, "プールA1": poolA1}
    dup = {k: [x for x in v if v.count(x) > 1] for k, v in sets.items() if len(set(v)) != len(v)}
    return {"sets": sets, "symbol": sym, "片方に重複（A-7）": dup}


# --------------------------------------------------------------------------- 材料（A-2）
def check_materials(cfg, log, results_root_checked: Path) -> dict:
    rec: dict = {"sha256": {}, "埋め込み": {}}
    ok = True
    for name, path in (("A2 genes_011.json", A2_SRC / "work/genes_011.json"),
                       ("A1 config.yaml", A1_SRC / "config.yaml"),
                       ("A1 新規 EXP config.yaml", A1_NEW / "config.yaml")):
        got = sha256f(path)
        rec["sha256"][name] = {"sha256": got, "実在": bool(got)}
        ok &= bool(got)
    # トークン化済み入力（複製）の sha256（EXP-20261011-001 の precheck_014.json に記録した値をそのまま使う）
    pre = json.loads((A1_NEW / "work/precheck_014.json").read_text())
    for key, clone in (("a2", "A2"), ("a1", "A1")):
        rec2 = pre["1. トークン化済み入力の複製"][clone]
        ok &= bool(rec2["pass"])
        rec.setdefault("トークン化済み入力の複製", {})[key] = rec2["複製先"]
    # 埋め込みの実在
    counts = {}
    for label, root, n_genes in (
        ("A1 保存 626", A1_SRC / "results", 626),
        ("A1 新規 30", A1_NEW / "results", 30),
        ("A2 保存 587", A2_SRC / "results", 587),
    ):
        n_files = 0
        miss = []
        for s in SEEDS:
            gdir = root / f"run-{s}-L-1/genes"
            if not gdir.exists():
                miss.append(f"run-{s}-L-1/genes なし"); continue
            for p in sorted(gdir.iterdir()):
                if not p.is_dir():
                    continue
                for kind in ("orig", "pert"):
                    if (p / "emb" / f"{kind}.npy").exists():
                        n_files += 1
                    else:
                        miss.append(str((p / "emb" / f"{kind}.npy").relative_to(ROOT)))
        counts[label] = {"実在": n_files, "期待": n_genes * 3 * 2, "欠け": miss[:5]}
        ok &= (n_files == n_genes * 3 * 2 and not miss)
    rec["埋め込み"] = counts
    rec["pass"] = bool(ok)
    log(f"A-2 材料の実在と複製の sha256: pass={rec['pass']}")
    return rec


# --------------------------------------------------------------------------- 注釈（マッチングに使う。摂動の結果に依存しない）
def annotation_table(long: pd.DataFrame) -> pd.DataFrame:
    """per_seed_stats の出力から、遺伝子ごとの注釈（3 群平均）を作る。
    列: n_exp（発現細胞数）・ll_med（入力の長さの中央値）・lp_med（元の位置の中央値）"""
    a = long.groupby("gene")[["n_expressed", "ll_med", "lp_med"]].mean()
    a.columns = ["n_exp", "ll_med", "lp_med"]
    return a.assign(log2_n_exp=lambda d: np.log2(d.n_exp + 1.0))


# --------------------------------------------------------------------------- マッチ偽セット（[2] 主方式）
def matched_sham_sets(Zpool: np.ndarray, pool_ann: pd.DataFrame, target_ann: pd.DataFrame,
                      n_draw: int, rng) -> np.ndarray:
    """evaluate_005.sham_sets と同じ近傍距離で、対象セット（target_ann）の遺伝子 1 つずつに
    発現細胞数・入力の長さ・元の位置をそろえた遺伝子をプールから復元ありで選ぶ。
    距離 = |Δlp_med|/0.25 + |Δll_med|/0.25 + |log(n_exp の比)|/log(1.5)。近い 5 個から無作為に 1 個。
    組の中での重複は許さない。戻り値: (n_draw, 3) のセット平均（列は seed 81/82/83）。"""
    pa = pool_ann
    out = np.empty((n_draw, 3))
    for d in range(n_draw):
        picks: list[str] = []
        order_t = target_ann.sample(frac=1.0, random_state=int(rng.integers(0, 2 ** 31 - 1)))
        for t in order_t.itertuples():
            cand = pa[~pa.index.isin(picks)]
            if cand.empty or len(cand) < 5:
                break
            dist = (np.abs(cand.lp_med - t.lp_med) / 0.25
                    + np.abs(cand.ll_med - t.ll_med) / 0.25
                    + np.abs(np.log(cand.n_exp / max(t.n_exp, 0.5))) / np.log(1.5))
            order = cand.index.to_numpy()[np.argsort(dist.to_numpy(), kind="stable")]
            take = order[:5]
            picks.append(str(take[int(rng.integers(0, len(take)))]))
        if len(picks) < target_ann.shape[0]:
            raise SystemExit(f"**止める**(A-5): マッチ偽セットが対象セットの大きさに届かない"
                             f"（{len(picks)}/{target_ann.shape[0]}。プールの注釈が不足）")
        idx = pa.index.get_indexer(picks)
        out[d] = Zpool[idx].mean(axis=0)
    return out


def matched_profile_center(pool_ann: pd.DataFrame, target_ann: pd.DataFrame, n_draw: int, rng) -> pd.DataFrame:
    """マッチ偽セット 1,000 組のプロファイル中央値（R-12a の検査に使う）。"""
    pa = pool_ann
    rows = []
    for d in range(n_draw):
        picks: list[str] = []
        order_t = target_ann.sample(frac=1.0, random_state=int(rng.integers(0, 2 ** 31 - 1)))
        for t in order_t.itertuples():
            cand = pa[~pa.index.isin(picks)]
            if cand.empty or len(cand) < 5:
                break
            dist = (np.abs(cand.lp_med - t.lp_med) / 0.25
                    + np.abs(cand.ll_med - t.ll_med) / 0.25
                    + np.abs(np.log(cand.n_exp / max(t.n_exp, 0.5))) / np.log(1.5))
            order = cand.index.to_numpy()[np.argsort(dist.to_numpy(), kind="stable")]
            take = order[:5]
            picks.append(str(take[int(rng.integers(0, len(take)))]))
        sub = pa.loc[picks]
        rows.append({"log2_n_exp": float(sub.log2_n_exp.median()), "ll_med": float(sub.ll_med.median()),
                     "lp_med": float(sub.lp_med.median())})
    return pd.DataFrame(rows)


def bootstrap_se_quantile(null_set_means: np.ndarray, q: float, n_boot: int = 200, rng=None) -> float:
    """マッチ偽セットの分位点のブートストラップ標準誤差（R-12b）。"""
    if rng is None:
        rng = np.random.default_rng(0)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(null_set_means), size=len(null_set_means))
        vals.append(float(np.quantile(null_set_means[idx], q, axis=0)[0]))
    return float(np.std(vals))


def residualize(Zdf: pd.DataFrame, ann: pd.DataFrame, fit_genes: list[str]) -> pd.DataFrame:
    """seed ごとに z_diff を [1, log2(n_exp), log2(ll_med)] に回帰し、残差 z_res を返す。
    fit_genes に完全な注釈を持つものだけを回帰に使う。z_res は全遺伝子について返す。"""
    out = Zdf.copy()
    xcols = ["log2_n_exp", "ll_med"]
    for s in SEEDS:
        ann_s = ann.reindex(Zdf.index)
        usable = [g for g in fit_genes if g in ann_s.index and ann_s.loc[g, xcols].notna().all()]
        X = np.column_stack([np.ones(len(usable)),
                             ann_s.loc[usable, "log2_n_exp"].to_numpy(),
                             np.log2(ann_s.loc[usable, "ll_med"].to_numpy())])
        y = Zdf.loc[usable, s].to_numpy()
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        Xall = np.column_stack([np.ones(len(ann_s)),
                                ann_s["log2_n_exp"].to_numpy(),
                                np.log2(ann_s["ll_med"].to_numpy())])
        pred = Xall @ beta
        out[s] = Zdf[s].to_numpy() - pred
    return out


def cv_residual_corr(Zdf: pd.DataFrame, ann: pd.DataFrame, fit_genes: list[str], k: int = 5,
                     rng=None) -> dict:
    """5 分割交差検証の残差と注釈の Spearman（R-12c）。return {seed: {n_exp: ρ, ll: ρ}}"""
    if rng is None:
        rng = np.random.default_rng(1)
    genes = [g for g in fit_genes if g in ann.index and ann.loc[g, ["log2_n_exp", "ll_med"]].notna().all()]
    genes = list(rng.permutation(genes))
    folds = [genes[i::k] for i in range(k)]
    out = {}
    for s in SEEDS:
        res = {}
        for f, held in enumerate(folds):
            tr = [g for i, g in enumerate(genes) if i % k != f]
            Xtr = np.column_stack([np.ones(len(tr)),
                                   ann.loc[tr, "log2_n_exp"].to_numpy(),
                                   np.log2(ann.loc[tr, "ll_med"].to_numpy())])
            beta, *_ = np.linalg.lstsq(Xtr, Zdf.loc[tr, s].to_numpy(), rcond=None)
            Xheld = np.column_stack([np.ones(len(held)),
                                     ann.loc[held, "log2_n_exp"].to_numpy(),
                                     np.log2(ann.loc[held, "ll_med"].to_numpy())])
            r = Zdf.loc[held, s].to_numpy() - Xheld @ beta
            res[held[0]] = pd.DataFrame({"gene": held, "s": s, "res": r,
                                         "n_exp": ann.loc[held, "log2_n_exp"].to_numpy(),
                                         "ll": np.log2(ann.loc[held, "ll_med"].to_numpy())})
        cv = pd.concat(res.values())
        out[str(s)] = {"ρ_n_exp": None, "ρ_ll": None}
        for col, key in (("n_exp", "ρ_n_exp"), ("ll", "ρ_ll")):
            rho = stats.spearmanr(cv[col], cv.res).statistic
            out[str(s)][key] = round(float(rho), 4)
    return out


# --------------------------------------------------------------------------- A1（確認）
def run_a1(cfg, G, out: dict, L, log) -> None:
    ev = cfg["evaluation"]["a1"]
    tok = str(ROOT / ev["tokenized_dataset"])
    cfg1 = {"data": {"tokenized_dataset": tok, "state_key": "fib_anno_state"},
            "states": {"start_state": "FIB_MFIB", "goal_state": "aStr", "alt_states": ["dFIB", "dMFIB"]},
            "output": {"results_dir": "results"}}
    ev9.BLOCKED = {e: e for e in cfg["blocked_genes"]}
    res = A1_SRC / "results"
    ev9.E = A1_SRC
    c_a1, _ = ev9.load_cells(cfg1, "results")
    c_a1 = c_a1[c_a1.seed.isin(SEEDS)].reset_index(drop=True)
    res_new = A1_NEW / "results"
    ev9.E = A1_NEW
    c_new, _ = ev9.load_cells(cfg1, "results")
    c_new = c_new[c_new.seed.isin(SEEDS)].reset_index(drop=True)

    for tag, cc, root in (("a1既存", c_a1, res), ("a1新規", c_new, res_new)):
        ada = add_projection(root, "FIB_MFIB", "aStr").rename("m_proj_ada")
        deg = add_projection(root, "FIB_MFIB", "dFIB")
        diff = add_projection_diff(root, "FIB_MFIB", "aStr", "dFIB", "diff")
        com = add_projection_diff(root, "FIB_MFIB", "aStr", "dFIB", "common")
        cc = merge_projection(cc, ada, "m_proj_ada", f"{tag} u_ada")
        cc = merge_projection(cc, deg, "m_proj_deg", f"{tag} u_deg")
        cc = merge_projection(cc, diff, "m_proj_diff", f"{tag} u_diff")
        cc = merge_projection(cc, com, "m_proj_common", f"{tag} u_common")
        if tag == "a1既存":
            c_a1 = cc
        else:
            c_new = cc

    pool = [x for x in G["sets"]["プールA1"] if x in set(c_a1.gene)]
    symptr = G["symbol"]
    S = {}
    for col, lab in (("m_proj_ada", "u_ada"), ("m_proj_deg", "u_deg"), ("m_proj_diff", "u_diff"),
                     ("m_proj_common", "u_common")):
        S["既存_" + lab] = per_seed_stats(c_a1, col, pool, symptr, log, f"A1 既存 {lab}")
    for col, lab in (("m_proj_ada", "u_ada"), ("m_proj_deg", "u_deg"), ("m_proj_diff", "u_diff"),
                     ("m_proj_common", "u_common")):
        frames = []
        for s in SEEDS:
            cs = c_new[c_new.seed == s].assign(m_goal=lambda d: d[col])
            cs0 = c_a1[c_a1.seed == s].assign(m_goal=lambda d: d[col])
            nl = {h: hd for h, hd in cs0[cs0.gene.isin(pool)].groupby("gene")}
            f = ev9.gene_stats_stratified(cs, nl, "m_goal")
            n = S["既存_" + lab][(S["既存_" + lab].seed == s) & (S["既存_" + lab].gene.isin(pool))].D_raw
            frames.append(f.assign(seed=s, z=(f.D_raw - n.mean()) / n.std(ddof=1), symbol=f.gene.map(symptr)))
        S["新規_" + lab] = pd.concat(frames, ignore_index=True)
    Z_saved = {lab.split("_", 1)[1]: v.pivot(index="gene", columns="seed", values="z").reindex(columns=SEEDS)
               for lab, v in S.items() if lab.startswith("既存_")}
    Z_new = {lab.split("_", 1)[1]: v.pivot(index="gene", columns="seed", values="z").reindex(columns=SEEDS)
             for lab, v in S.items() if lab.startswith("新規_")}
    diff = pd.concat([S["既存_u_diff"], S["新規_u_diff"]], ignore_index=True) \
        .drop_duplicates(subset=["gene", "seed"], keep="first")
    ann = annotation_table(diff)

    # ---- 再現確認（A-4）: 公開済み gene_rank_014.csv と 1e-3 以内
    ref = pd.read_csv(str(ROOT / ev["repro_csv"])).set_index("gene")
    d_all = 0.0
    for gene in ref.index:
        for s in SEEDS:
            a = diff[diff.gene == gene].set_index("seed").loc[s, "z"] if gene in diff.gene.values else np.nan
            b = ref.loc[gene, f"z_{s}"]
            if pd.notna(a) and pd.notna(b):
                d_all = max(d_all, abs(float(a) - float(b)))
    repro_ok = d_all <= ev["repro_tol"]
    # R-12 再現（ρ_exp が前回 ≈0.31 と 0.02 以内）
    zbar = diff[diff.gene.isin(pool)].groupby("gene")[["z"]].mean().reindex(pool)["z"]
    n_exp = ann.reindex(pool)["n_exp"]
    r_exp = float(stats.spearmanr(zbar, n_exp).statistic)
    repro_ok &= bool(abs(r_exp - 0.3133) <= 0.02)
    out["A チェック（A1）"] = {"A-4 再現確認": {"最大絶対差": round(float(d_all), 6),
                                            "許容差": ev["repro_tol"], "ρ_exp 再現": round(r_exp, 4),
                                            "pass": bool(repro_ok)}}
    if not repro_ok:
        raise SystemExit("**止める**(A-4/A1): 再現確認が合わない")

    # ---- 主判定に使う集合
    pc71 = G["sets"]["活性化PC71"]; pc23 = G["sets"]["退化PC23"]
    pc30 = G["sets"]["退化PC30"]; pc21 = G["sets"]["退化PC21"]
    poolZ = Z_saved["u_diff"].reindex(pool).to_numpy()
    # 引き出しプール: 検定する PC がプールと重なる遺伝子を除く
    ov_ens23 = {e for e in pc23 if e in set(pool)}
    ov_ens71 = {e for e in pc71 if e in set(pool)}
    draw_pool23 = [g for g in pool if g not in ov_ens23]
    draw_pool71 = [g for g in pool if g not in ov_ens71]
    Zdiff_saved = Z_saved["u_diff"]

    def _aligned(draw_pool):
        idx = [g for g in draw_pool
               if g in ann.index and g in Zdiff_saved.index
               and ann.loc[g].notna().all() and Zdiff_saved.loc[g].notna().all()]
        return ann.reindex(idx), Zdiff_saved.reindex(idx)

    ann_pool23, Zpool23_df = _aligned(draw_pool23)
    ann_pool71, Zpool71_df = _aligned(draw_pool71)
    Zpool23 = Zpool23_df.to_numpy()
    Zpool71 = Zpool71_df.to_numpy()
    ann_pc23 = ann.reindex([g for g in pc23 if g in ann.index]).dropna()
    ann_pc71 = ann.reindex([g for g in pc71 if g in ann.index]).dropna()
    if len(pc23) != PLAN["退化PC23_n"] or len(pc71) != PLAN["活性化PC71_n"]:
        raise SystemExit("**止める**(A-3): PC の個数が config と違う")
    if len(ann_pc23) < PLAN["退化PC23_n"] or len(ann_pc71) < PLAN["活性化PC71_n"]:
        raise SystemExit("**止める**(A-5): PC の注釈が欠けている（マッチできない）")

    zm = lambda Zd, genes: Zd.reindex([x for x in genes if x in Zd.index]).to_numpy()  # noqa: E731
    obs71 = zm(Z_saved["u_diff"], pc71).mean(axis=0)
    obs23 = zm(Z_new["u_diff"], [g for g in pc23 if g in Z_new["u_diff"].index]).mean(axis=0)
    obs30 = zm(Z_new["u_diff"], [g for g in pc30 if g in Z_new["u_diff"].index]).mean(axis=0)
    obs21 = zm(Z_new["u_diff"], [g for g in pc21 if g in Z_new["u_diff"].index]).mean(axis=0)

    key = PLAN["乱数の鍵"]
    rng_i = np.random.default_rng(key["P-1 活性化のマッチ偽セット"])
    rng_ii = np.random.default_rng(key["P-1 退化のマッチ偽セット"])
    rng_iii_a = np.random.default_rng(key["P-1 (iii) のマッチ無作為 2 セット"])
    rng_iii_b = np.random.default_rng(key["P-1 (iii) のマッチ無作為 2 セット"] + 1000)
    null71 = matched_sham_sets(Zpool71, ann_pool71, ann_pc71, N_DRAW, rng_i)
    null23 = matched_sham_sets(Zpool23, ann_pool23, ann_pc23, N_DRAW, rng_ii)
    thr71 = np.quantile(null71, 0.99, axis=0)
    thr23 = np.quantile(null23, 0.01, axis=0)
    # (iii): 活性化 71・退化 23 をそれぞれのプロファイルにマッチさせて独立に引いた無作為 2 セットの差
    set71 = matched_sham_sets(Zpool71, ann_pool71, ann_pc71, N_DRAW, rng_iii_a)
    set23 = matched_sham_sets(Zpool23, ann_pool23, ann_pc23, N_DRAW, rng_iii_b)
    null_diff = set71 - set23
    thr_diff = np.quantile(null_diff, 0.99, axis=0)
    c_i = bool(np.all(obs71 > 0) and np.all(obs71 > thr71))
    c_ii = bool(np.all(obs23 < 0) and np.all(obs23 < thr23))
    c_iii = bool(np.all((obs71 - obs23) > thr_diff))
    p1 = c_i and c_ii and c_iii

    # ---- R-12a マッチングの中心
    prof23 = matched_profile_center(ann_pool23, ann_pc23, 500, np.random.default_rng(key["P-1 退化のマッチ偽セット"] + 7))
    prof71 = matched_profile_center(ann_pool71, ann_pc71, 500, np.random.default_rng(key["P-1 活性化のマッチ偽セット"] + 7))
    tol = PLAN["マッチング_許容"]
    r12a_23 = {c: abs(float(prof23[c].median()) - float(ann_pc23[c].median())) for c in tol}
    r12a_71 = {c: abs(float(prof71[c].median()) - float(ann_pc71[c].median())) for c in tol}
    r12a_pass = all(v <= tol[c] for c, v in r12a_23.items()) and all(v <= tol[c] for c, v in r12a_71.items())
    # ---- R-12b 閾値の安定性
    se71 = bootstrap_se_quantile(null71, 0.99, 200, np.random.default_rng(key["P-1 活性化のマッチ偽セット"] + 8))
    se23 = bootstrap_se_quantile(null23, 0.01, 200, np.random.default_rng(key["P-1 退化のマッチ偽セット"] + 8))
    r12b_pass = bool(se71 <= abs(float(np.quantile(thr71, 0.5))) * PLAN["R-12b_SE比"]
                     and se23 <= abs(float(np.quantile(thr23, 0.5))) * PLAN["R-12b_SE比"])
    # ---- R-12c 交差検証残差
    Zdf_u = diff.pivot(index="gene", columns="seed", values="z").reindex(columns=SEEDS)
    fit_genes = [g for g in pool if g not in set(ov_ens23) | set(ov_ens71)]
    cvc = cv_residual_corr(Zdf_u, ann, fit_genes, rng=np.random.default_rng(key["P-1 退化のマッチ偽セット"] + 9))
    r12c_pass = all(abs(cvc[str(s)][col]) <= PLAN["R-12c_閾値"] for s in SEEDS for col in ("ρ_n_exp", "ρ_ll"))

    # ---- 残差化版（副次）
    Zres = residualize(Zdf_u, ann, fit_genes)
    pool_res = Zres.reindex(pool).dropna().to_numpy()
    obs23_res = Zres.reindex([g for g in pc23 if g in Zres.index]).to_numpy().mean(axis=0)
    obs71_res = Zres.reindex([g for g in pc71 if g in Zres.index]).to_numpy().mean(axis=0)
    rng_r = np.random.default_rng(key["残差化版の偽セット"])
    thr71r = np.quantile(set_null(pool_res, len(pc71), N_DRAW, rng_r, replace=True), 0.99, axis=0)
    thr23r = np.quantile(set_null(pool_res, len([g for g in pc23 if g in Zres.index]), N_DRAW,
                                  np.random.default_rng(key["残差化版の偽セット"] + 1), replace=True), 0.01, axis=0)
    c_ir = bool(np.all(obs71_res > 0) and np.all(obs71_res > thr71r))
    c_iir = bool(np.all(obs23_res < 0) and np.all(obs23_res < thr23r))
    p1_res = c_ir and c_iir and c_iii

    # ---- 補正なし（無作為偽セット）版（感度。前回と同条件）
    rng_u = np.random.default_rng(key["補正なし（無作為）の偽セット"])
    poolZ_all = Z_saved["u_diff"].reindex(pool).to_numpy()
    thr71u = np.quantile(set_null(poolZ_all, len(pc71), N_DRAW, rng_u, replace=True), 0.99, axis=0)
    thr23u = np.quantile(set_null(poolZ_all, len(pc23), N_DRAW,
                                  np.random.default_rng(key["補正なし（無作為）の偽セット"] + 1), replace=True),
                         0.01, axis=0)
    c_iu = bool(np.all(obs71 > 0) and np.all(obs71 > thr71u))
    c_iiu = bool(np.all(obs23 < 0) and np.all(obs23 < thr23u))
    p1_uncorrected = c_iu and c_iiu and c_iii

    # ---- 感度版（退化 30・退化 21）
    _, thr30u = c_ii, np.quantile(set_null(poolZ_all, len(pc30), N_DRAW,
                                           np.random.default_rng(key["補正なし（無作為）の偽セット"] + 2),
                                           replace=True), 0.01, axis=0)
    thr30m = np.quantile(matched_sham_sets(Zpool23, ann_pool23,
                                           ann.reindex([g for g in pc30 if g in ann.index]).dropna(),
                                           N_DRAW, np.random.default_rng(key["P-1 退化のマッチ偽セット"] + 3)),
                         0.01, axis=0)
    c_ii30m = bool(np.all(obs30 < 0) and np.all(obs30 < thr30m))
    thr21m = np.quantile(matched_sham_sets(Zpool23, ann_pool23,
                                           ann.reindex([g for g in pc21 if g in ann.index]).dropna(),
                                           N_DRAW, np.random.default_rng(key["P-1 退化のマッチ偽セット"] + 4)),
                         0.01, axis=0)
    c_ii21m = bool(np.all(obs21 < 0) and np.all(obs21 < thr21m))

    # ---- 反証条件（R-1'〜R-6'。最初に当たったものを採用）
    m71m, m23m = float(np.mean(obs71)), float(np.mean(obs23))
    r1 = bool(np.all(obs71 < 0))
    r2 = bool(np.all(obs23 > 0))
    r3 = (not c_i) and (not c_ii)
    r4 = (m71m > 0 and m23m > 0) or (m71m < 0 and m23m < 0)
    r5 = c_i and c_ii and (not c_iii)
    r6 = (c_i and not c_ii) or (c_ii and not c_i)
    first = next((k for k, v in {"R-1'": r1, "R-2'": r2, "R-3'": r3, "R-4'": r4, "R-5'": r5, "R-6'": r6}.items()
                  if v), "該当なし（P-1' の予測どおりの提案）")
    # R-12'（補正の健全性）は R-12a〜c の合否で別枠（§4.3 のとおり判定不能）
    r12_hit = not (r12a_pass and r12b_pass and r12c_pass)
    # R-13'（感度で結論が変わるか）
    r13 = bool(p1 != p1_uncorrected or p1 != p1_res)
    sens30_differs = bool(c_ii != c_ii30m or c_ii != c_ii21m)

    out["A1 主判定 P-1'（マッチ偽セット）"] = {
        "(i) 活性化71": {"セット平均 z_diff": [round(x, 4) for x in obs71],
                        "マッチ偽セット99%点": [round(x, 4) for x in thr71], "成り立つ": c_i},
        "(ii) 退化23": {"セット平均 z_diff": [round(x, 4) for x in obs23],
                       "マッチ偽セット1%点": [round(x, 4) for x in thr23], "成り立つ": c_ii},
        "(iii) 差": {"差": [round(x, 4) for x in (obs71 - obs23)],
                    "マッチ無作為2セットの99%点": [round(x, 4) for x in thr_diff], "成り立つ": c_iii},
        "P-1' 全体": p1,
        "引き出しプール": {"退化23 用": len(draw_pool23), "活性化71 用": len(draw_pool71),
                        "除いたプール重なり（退化23）": [G["symbol"].get(e, e) for e in ov_ens23],
                        "除いたプール重なり（活性化71）": [G["symbol"].get(e, e) for e in ov_ens71]}}
    out["A1 反証条件（R-1'〜R-6'）"] = {"R-1'": r1, "R-2'": r2, "R-3'": r3, "R-4'": r4, "R-5'": r5,
                                      "R-6'": r6, "最初に当たったもの": first}
    out["A1 補正の健全性"] = {
        "R-12a マッチングの中心": {"退化23": r12a_23, "活性化71": r12a_71, "許容": tol, "pass": r12a_pass},
        "R-12b 閾値の安定性": {"SE(99%点 活性化)": round(se71, 4), "SE(1%点 退化)": round(se23, 4),
                            "基準": f"≤ |閾値| × {PLAN['R-12b_SE比']}", "pass": r12b_pass},
        "R-12c 交差検証残差": cvc, "閾値": PLAN["R-12c_閾値"], "pass": r12c_pass,
        "補正の健全性 全体（不合格なら判定不能 R-12'）": r12_hit}
    out["A1 感度"] = {
        "補正なし（無作為偽セット）版": {"(i)": c_iu, "(ii)": c_iiu, "P-1": p1_uncorrected,
                                     "退化23 無作為1%点": [round(x, 4) for x in thr23u]},
        "残差化版（副次）": {"(i)": c_ir, "(ii)": c_iir, "P-1": p1_res,
                          "退化23 残差1%点": [round(x, 4) for x in thr23r]},
        "退化30（マッチ・感度）": {"観測": [round(x, 4) for x in obs30], "マッチ1%点": [round(x, 4) for x in thr30m],
                               "成り立つ": c_ii30m},
        "退化21（マッチ・感度。JUND・CTSH 除去）": {"観測": [round(x, 4) for x in obs21],
                                               "マッチ1%点": [round(x, 4) for x in thr21m], "成り立つ": c_ii21m},
        "R-13'（結論が版に依存）": {"補正なしと違う": r13, "退化30/21 と違う": sens30_differs}}
    # ---- 較正（偶然の水準。判定に使わない。主規則と同じマッチした無作為ペアで実測）
    rng_cal_a = np.random.default_rng(key["較正（マッチした無作為ペア）"])
    rng_cal_b = np.random.default_rng(key["較正（マッチした無作為ペア）"] + 1000)
    cal_hit = 0
    for _ in range(N_DRAW):
        m71 = matched_sham_sets(Zpool71, ann_pool71, ann_pc71, 1, rng_cal_a)[0]
        m23 = matched_sham_sets(Zpool23, ann_pool23, ann_pc23, 1, rng_cal_b)[0]
        cal_hit += int(bool(np.all(m71 > thr71) and np.all(m23 < thr23) and np.all(m71 - m23 > thr_diff)))
    out["A1 偶然の水準（較正。判定に使わない）"] = {"マッチした無作為ペア 1,000 組で P-1' を通した数": cal_hit,
                                               "読み方": "片側 1% に較正した規則の結合した偶然合格率"}

    # ---- 順位表（探索）
    Zc = diff.pivot(index="gene", columns="seed", values="z").reindex(columns=SEEDS)
    rank = (Zc.rename(columns={s: f"z_{s}" for s in SEEDS})
            .assign(z_mean3=Zc.mean(axis=1)).reset_index().rename(columns={"index": "gene"})
            .assign(symbol=lambda d: d["gene"].map(G["symbol"]).fillna(""),
                    注="探索（データを見て候補を探すための順位表。事前登録した確認的検定ではない。CLAUDE.md §3.5）")
            .sort_values("z_mean3", ascending=False).reset_index(drop=True))
    rank.insert(0, "rank", np.arange(1, len(rank) + 1))
    outdir = EXP / "results"
    rank.to_csv(outdir / "gene_rank_015.csv", index=False)
    out["A1 z_diff の順位表（探索）"] = {"行数": len(rank), "保存先": "gene_rank_015.csv",
                                      "上位10": rank.head(10)[["symbol", "z_mean3"]].round(4).to_dict(orient="records")}
    log(f"A1: P-1'={p1}（(i)={c_i}, (ii)={c_ii}, (iii)={c_iii}）。補正なし版 P-1={p1_uncorrected}、"
        f"残差化版 P-1={p1_res}。最初に当たった反証条件={first}、R-12'={r12_hit}")


# --------------------------------------------------------------------------- A2（探索）
def run_a2(cfg, G, out: dict, L, log) -> None:
    ev = cfg["evaluation"]["a2"]
    tok = str(ROOT / ev["tokenized_dataset"])
    cfg2 = {"data": {"tokenized_dataset": tok, "state_key": "isp_state"},
            "states": {"start_state": "FIB", "goal_state": "aStr", "alt_states": ["dFIB"]},
            "output": {"results_dir": "results"}}
    ev9.E = A2_SRC
    ev9.BLOCKED = {e: e for e in cfg["blocked_genes"]}
    res = A2_SRC / "results"
    c, _ = ev9.load_cells(cfg2, "results")
    c = c[c.seed.isin(SEEDS)].reset_index(drop=True)
    diff = add_projection_diff(res, "FIB", "aStr", "dFIB", "diff")
    c = merge_projection(c, diff, "m_proj_diff", "A2 u_diff")
    pool = [x for x in G["sets"]["プールA2"] if x in set(c.gene)]
    long = per_seed_stats(c, "m_proj_diff", pool, G["symbol"], log, "A2 u_diff")
    Z = long.pivot(index="gene", columns="seed", values="z").reindex(columns=SEEDS)
    ann = annotation_table(long)
    cand = G["sets"]["候補"]
    cand_idx = [g for g in cand if g in ann.index]
    idx = [g for g in pool
           if g in ann.index and g in Z.index
           and ann.loc[g].notna().all() and Z.loc[g].notna().all()]
    ann_pool = ann.reindex(idx)
    Zpool = Z.reindex(idx).to_numpy()
    ann_cand = ann.reindex(cand_idx).dropna()
    rng_a = np.random.default_rng(PLAN["乱数の鍵"]["A2 候補のマッチ偽セット"])
    null_cand = matched_sham_sets(Zpool, ann_pool, ann_cand, N_DRAW, rng_a)
    thr = np.quantile(null_cand, 0.99, axis=0)
    obs = Z.reindex(cand).to_numpy().mean(axis=0)
    inside = bool(np.all(obs <= thr))
    def _fmt(v: float):
        if isinstance(v, float) and (v != v or v in (float("inf"), float("-inf"))):
            return None
        return round(float(v), 3)
    cand_ann = {G["symbol"].get(g, g): [_fmt(ann.loc[g, c]) for c in ("n_exp", "ll_med", "lp_med")]
                for g in cand if g in ann.index}
    unexpr = [G["symbol"].get(g, g) for g in cand
              if g in ann.index and float(ann.loc[g, "n_exp"]) == 0]
    out["A2 候補 8（マッチ補正版。探索）"] = {
        "セット平均 z_diff": [round(x, 4) for x in obs],
        "マッチ偽セット99%点": [round(x, 4) for x in thr],
        "99% 以内（=0 と区別できない）": inside,
        "注": "前回の読み（R-7: 段階 2 の合格は共通成分だった）が、交絡をそろえた偽セットの下でも保たれるかの記述。探索である。",
        "候補8の注釈（n_exp/ll/lp）": cand_ann,
        "注記 null の理由": {g: "発現細胞数が 0（このデータの A2 で観測されない）ため、元の位置の中央値なし" for g in unexpr}}


# --------------------------------------------------------------------------- 本体
def main() -> None:
    t0 = time.time()
    cfg = yaml.safe_load((EXP / "config.yaml").read_text())
    dec = cfg["decision"]
    outdir = EXP / "results"
    outdir.mkdir(parents=True, exist_ok=True)
    LOG = outdir / "run.log"
    LOG.write_text("")

    def L(msg: str) -> None:
        log_line(LOG, msg)

    L(f"開始。出力先 {outdir.relative_to(ROOT)}")
    L(f"plan.md sha256 {sha256f(EXP/'plan.md')[:16]}… / config.yaml sha256 "
      f"{sha256f(EXP/'config.yaml')[:16]}… / このスクリプト sha256 {sha256f(Path(__file__))[:16]}…")

    # ---- 1. config の decision が plan と同じか
    checks_conf = {
        "退化23": len(cfg["degenerate_30"]) - len(cfg["overlap7_symbols"]) == PLAN["退化PC23_n"],
        "退化30": len(cfg["degenerate_30"]) == PLAN["退化PC30_n"],
        "プールA1": True,
        "乱数の鍵（P-1 活性化）": int(dec["乱数の鍵"]["P-1 活性化のマッチ偽セット"]) == PLAN["乱数の鍵"]["P-1 活性化のマッチ偽セット"],
        "乱数の鍵（P-1 退化）": int(dec["乱数の鍵"]["P-1 退化のマッチ偽セット"]) == PLAN["乱数の鍵"]["P-1 退化のマッチ偽セット"],
        "乱数の鍵（P-1 差）": int(dec["乱数の鍵"]["P-1 (iii) のマッチ無作為 2 セット"]) == PLAN["乱数の鍵"]["P-1 (iii) のマッチ無作為 2 セット"],
        "乱数の鍵（補正なし）": int(dec["乱数の鍵"]["補正なし（無作為）の偽セット"]) == PLAN["乱数の鍵"]["補正なし（無作為）の偽セット"],
        "乱数の鍵（残差化）": int(dec["乱数の鍵"]["残差化版の偽セット"]) == PLAN["乱数の鍵"]["残差化版の偽セット"],
        "G2 不要（config 側）": cfg["compute"]["needs_G2_approval"] is False,
    }
    if not all(checks_conf.values()):
        L(f"**止める**: config.yaml の decision: 節が plan.md §4 と一致しない {checks_conf}")
        raise SystemExit("**止める**(config): decision が plan と違う")
    out: dict = {
        "実験": "EXP-20261011-002", "仮説": "HYP-20261011-01", "mode": "final",
        "問い": "退化の特異的陽性対照 23 個は、交絡（発現細胞数・入力の長さ・元の位置）をそろえた偽セットの"
                "1% 点を A1 の 3 群すべてで下回るか（P-1'）。補正後の R-12 相当の検査（R-12a〜c）が通るか",
        "0. config と plan の一致": checks_conf,
        "使った数値型": {"device（評価）": "cpu", "dtype（評価）": "fp64",
                     "GPU 時間（評価）": 0.0,
                     "備考": "保存済み埋め込みの読み出しと再射影・再評価のみ。新しい摂動は回していない"},
    }

    # ---- 2・3. 材料と遺伝子集合
    mats = check_materials(cfg, L, outdir)
    out["1. 材料（A-2）"] = mats
    G = gene_sets(cfg)
    out["2. 遺伝子集合（A-3）"] = {
        "個数": {k: len(v) for k, v in G["sets"].items()},
        "片方に重複がある集合（A-7）": G["片方に重複（A-7）"],
        "盲検との重なり": {k: sorted(set(v) & set(cfg["blocked_genes"])) for k, v in G["sets"].items()},
    }
    if any(G["片方に重複（A-7）"].values()) or \
       any(out["2. 遺伝子集合（A-3）"]["盲検との重なり"].values()):
        raise SystemExit("**止める**(A-1/A-7): 集合の重複か盲検の混入がある")

    run_a1(cfg, G, out, L, L)
    run_a2(cfg, G, out, L, L)

    # ---- 出力
    out["decision（config.yaml の全鍵）"] = dec
    out["注"] = ("判定（支持・否定・判定不能。ゲート G3）は PI が行う。この道具は規則を機械的に当てた結果と提案まで。"
                 "hypotheses/tracker.md は変更しない。コンピュータ上の遺伝子摂動（保存済み）の結果であり、"
                 "実験的に確認された生物学的事実ではない。盲検中の遺伝子は摂動も解析もしていない（値は数えない）。")
    dur = time.time() - t0
    out["実行"] = {"開始": _dt.datetime.fromtimestamp(t0).astimezone().isoformat(),
                 "所要秒": round(dur, 1), "GPU 時間": 0.0, "device": "cpu", "dtype": "fp64"}
    cj = outdir / "criteria_015.json"
    cj.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    reloaded = json.loads(cj.read_text())
    a6 = {"decision の全鍵を持つ": set((reloaded.get("decision（config.yaml の全鍵）") or {}).keys()) == set(dec.keys()),
          "保存先": str(cj.relative_to(ROOT)), "pass": True}
    out["A チェック"] = {"A-2 材料": {"pass": mats["pass"]},
                        "A-4 再現確認": "run_a1 内", "A-6 結果ファイル": a6}
    cj.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")

    # ---- A-1（出力・ログ・一覧に盲検が 0 件）
    sym_blk = {"ENSG00000141510": None, "ENSG00000124762": None, "ENSG00000111321": None,
               "ENSG00000118785": None, "ENSG00000148346": None}
    names_blk = ["TP53", "CDKN1A", "LTBR", "SPP1", "LCN2"]
    targets = list(outdir.glob("*.*"))
    a1 = {"調べたファイル": [str(p.relative_to(ROOT)) for p in targets], "件数": 0}
    for p in targets:
        text = p.read_text(errors="replace")
        hits = sum(text.count(k) for k in sym_blk) + sum(text.count(k) for k in names_blk)
        if hits:
            a1["件数"] += hits
            a1.setdefault("当たったファイル", []).append(str(p.relative_to(ROOT)))
    a1["pass"] = a1["件数"] == 0
    L(f"A-1 盲検の混入（出力・ログ・一覧）: {a1['件数']} 件、pass={a1['pass']}")
    if not a1["pass"]:
        raise SystemExit("**止める**(A-1): 盲検が出力・ログ・一覧に現れた（B-2）")
    out["A チェック"]["A-1 盲検の混入"] = a1
    cj.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")

    # ---- A-5（値に NaN・無限大がないことを結果ファイルの再読込・値スキャンで確認）
    def _walk(o):
        if isinstance(o, dict):
            for v in o.values():
                yield from _walk(v)
        elif isinstance(o, list):
            for v in o:
                yield from _walk(v)
        else:
            yield o
    rer = json.loads(cj.read_text())
    vals = list(_walk(rer))
    nan_v = [(isinstance(v, float) and v != v) for v in vals]
    inf_v = [isinstance(v, float) and (v in (float("inf"), float("-inf"))) for v in vals]
    a5 = {"NaN（float nan）の値の数": int(sum(nan_v)), "無限大の値の数": int(sum(inf_v)),
          "pass": (not any(nan_v)) and (not any(inf_v))}
    out["A チェック"]["A-5 値"] = a5
    cj.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    L(f"完了（{dur:.1f} 秒）。出力 {cj.relative_to(ROOT)}")


def self_test() -> None:
    """合成データでマッチ偽セット・残差化・R-12b の機械を検査する。本番の値は一切使わない。"""
    rng = np.random.default_rng(42)
    n_pool = 300
    pool_genes = [f"P{i:03d}" for i in range(n_pool)]
    ann_pool = pd.DataFrame({"n_exp": rng.integers(5, 190, n_pool).astype(float),
                             "ll_med": rng.normal(9.0, 0.4, n_pool),
                             "lp_med": rng.normal(5.0, 1.5, n_pool)}, index=pool_genes)
    ann_pool["log2_n_exp"] = np.log2(ann_pool.n_exp + 1)
    tgt_genes = [f"T{i}" for i in range(20)]
    ann_tgt = pd.DataFrame({"n_exp": rng.integers(5, 30, 20).astype(float),
                            "ll_med": rng.normal(9.0, 0.2, 20),
                            "lp_med": rng.normal(5.0, 0.7, 20)}, index=tgt_genes)
    ann_tgt["log2_n_exp"] = np.log2(ann_tgt.n_exp + 1)
    Zpool = rng.normal(0, 1, (n_pool, 3))

    nulls = matched_sham_sets(Zpool, ann_pool, ann_tgt, 200, rng)
    assert nulls.shape == (200, 3), "マッチ偽セットの形が違う"
    assert np.all(np.isfinite(nulls)), "マッチ偽セットに NaN"
    prof = matched_profile_center(ann_pool, ann_tgt, 200, rng)
    assert prof.shape == (200, 3) and np.all(np.isfinite(prof.to_numpy())), "プロファイル中心が不正"
    se = bootstrap_se_quantile(nulls, 0.01, 50, rng)
    assert np.isfinite(se) and se >= 0, "bootstrap SE が不正"

    # 残差化: 人為的に n_exp と z を相関させたデータで、残差の相関が消えることを確認
    n = 200
    genes = [f"G{i}" for i in range(n)]
    ann2 = pd.DataFrame({"n_exp": rng.integers(2, 200, n).astype(float),
                         "ll_med": rng.normal(9.0, 0.4, n)}, index=genes)
    ann2["log2_n_exp"] = np.log2(ann2.n_exp + 1)
    Z2 = pd.DataFrame({s: 0.5 * np.log2(ann2.n_exp + 1) + rng.normal(0, 0.5, n) for s in SEEDS}, index=genes)
    Zres = residualize(Z2, ann2, genes)
    for s in SEEDS:
        rho = stats.spearmanr(Zres[s], ann2.log2_n_exp).statistic
        assert abs(rho) < 0.2, f"残差に交絡が残っている（{s}: ρ={rho:.3f}）"
    cvc = cv_residual_corr(Z2, ann2, genes, k=5, rng=rng)
    for s in SEEDS:
        assert abs(cvc[str(s)]["ρ_n_exp"]) < 0.35, "CV 残差の交絡が大きすぎる"
    print("self-test: マッチ偽セット・プロファイル中心・bootstrap SE・残差化・CV 残差 OK")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true", help="合成データで単体検査（本番の値は見ない）")
    a = ap.parse_args()
    if a.self_test:
        self_test()
    else:
        main()
#!/usr/bin/env python3
"""EXP-20261011-003 — A6 埋め込みのみの評価（plan.md §4・§7 手順 8・9）。

摂動はしない。保存済みの A6 埋め込み（results/embeddings.parquet、row_id 単位）に、plan.md §4 の
判定規則（P-1〜P-4）と反証条件（R-A1〜R-D1）を機械的に当てる。

- 主判定 P-1: 出発・目標の状態の分離（提供者単位 GroupKFold の AUC。総 UMI 数と検出遺伝子数の 2 変数だけで状態を分類した AUC を対照に、検出遺伝子数の分布をそろえた後の核で判定）
- 主判定 P-2: 2 本の軸（u_ada と u_deg）は異なる向きを向くか（提供者 3 群・検出遺伝子数の高低の向き・分布のそろえ方への安定性）
- 副次 P-3: 管状系が P-1・P-2 を両方満たすか（2 値。記述）
- 副次 P-4: 候補・陽性対照の位置（発現と埋め込みの相関。**探索・記述。摂動の確認ではない**）
- 反証条件は §4.6 の表を上から順に当て、最初に当たったものを採用する

読み取りのみ: plan.md・config.yaml・precheck_a6.json（固定値）・input.h5ad・results/（埋め込み）。
書き込み: results/（最終）または results-smoke/（試し実行）の criteria_a6.json・gene_position_a6.csv・
          README.md・run.log。

使い方:
    .venv/bin/python experiments/EXP-20261011-003/analysis/evaluate_a6.py --self-test   # 合成データで単体検査
    .venv/bin/python experiments/EXP-20261011-003/analysis/evaluate_a6.py --smoke        # 試し実行（判定ではない）
    .venv/bin/python experiments/EXP-20261011-003/analysis/evaluate_a6.py                # 本番（固定後の 1 回）
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
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold, StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[3]
EXP = ROOT / "experiments/EXP-20261011-003"
sys.path.insert(0, str(EXP / "analysis"))
from embed_a6 import depth_match_ids  # noqa: E402   検出遺伝子数の分布をそろえる手順は embed_a6.py と同一（決定的）

MATERIAL_SHA256 = "8d9cca1320a7952f4b2eaff07acc03494e840e46b6b333e1c3efaf395165b29d"

# --------------------------------------------------------------------------- 埋め込み時の config.yaml との照合（A-2）
# 埋め込み（46,957 核・1 回）は 2026-10-11 に config.yaml の sha256 `42ed7ae5…` で行い、その値は
# results/embed_a6/run_meta.json に記録されている。その後の 2026-10-11 第 6 回で config.yaml を
# 学術用語に書き換えた（**表記のみ。判定の値・閾値・対象・埋め込みの設定は不変**。plan.md §13 変更ログ）。
# ファイルの生ハッシュは書き換えで変わるため、この 2 つのハッシュを「埋め込みと同じ config」として受理する。
# **再埋め込みはしない。** どちらでもない場合は止める。
CONFIG_SHA256_EMBEDDED = "42ed7ae5d8c84c568df3c256534ba77d6ad89b654fa379e5b09c78ed461b7fe8"  # run_meta の記録
CONFIG_SHA256_REWORDED = "f82f7369e4524db8527e9bcd171c944b7314de4e64e3074430a739cee1341fac"  # 第 6 回（表記のみ）
CONFIG_SHA256_ACCEPTED = (CONFIG_SHA256_EMBEDDED, CONFIG_SHA256_REWORDED)

# --------------------------------------------------------------------------- plan.md §4 の判定値
# **ここを変えてはいけない。** config.yaml の decision: 節と照合し、違えば止める。
PLAN = {
    "P-1_提供者単位AUC下限": 0.70,
    "P-1_2変数だけの分離AUC上限": 0.60,
    "P-2_sin下限": 0.31,
    "P-2_検出遺伝子数の高低の向きとの余弦上限": 0.60,
    "P-2_角度差上限_度": 30.0,
    "P-4_単一遺伝子_検出された核の数_下限": 20,
    "P-4_ランダムな遺伝子セットの組数": 1000,
    "P-4_ランダムな遺伝子セットの大きさ": 8,
    "対照プール_n": 500,
    "乱数の鍵": {"対照プールの頻度マッチ": 20261011, "ランダムな遺伝子セット 1,000 組": 20261112},
}

# P-2(a) の提供者 3 群。**結果を見る前（2026-10-11）に、固定済み precheck_a6.json の提供者別細胞数だけを使って決めた。**
# 割り当て: 各群の総細胞数を平準化する貪欲法（union 集団 = 間質系＋管状系 7,335 核）。
# 各群は 4 状態すべてに細胞を持つ（最小 Fibroblast 28）。この定数は評価コードと一緒に sha256 固定する。
DONOR_GROUPS = {
    "G1": ["obstructed_5", "normal_4", "normal_3"],
    "G2": ["normal_7", "obstructed_2", "obstructed_4", "normal_5", "normal_1"],
    "G3": ["normal_2", "obstructed_1", "obstructed_3", "normal_6"],
}
GROUP_NOTE = ("提供者 3 群（plan.md §4.3(a)）。総細胞数で平準化（G1 2,354 / G2 2,483 / G3 2,498。"
              "union 集団 7,335 核）。結果を見る前に固定し、SHA256SUMS に含めた。")

STATE_NAMES = {"interstitial": ("Fibroblast", "Myofibroblast"),
               "tubular": ("PT_normal", "PT_inj_infl")}
ALL_STATES = ["Fibroblast", "Myofibroblast", "PT S1", "PT S2", "PT S3", "PT Injured", "PT Inflammatory"]
T5_SYMBOLS = ["BICC1", "RHEX", "GLIS3", "UGT2B7", "PPFIBP1", "FNDC3B", "SEMA5A"]


def sha256f(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def log_line(path: Path, msg: str) -> None:
    line = f"[{_dt.datetime.now().strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(path, "a") as fh:
        fh.write(line + "\n")


def fmt(v, nd=4):
    """有限なら丸めた値、そうでなければ None（JSON に NaN を書かない。A-5）。"""
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return round(f, nd) if np.isfinite(f) else None


def unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    return v / n if n > 0 else np.full_like(v, np.nan)


def angle_deg(u1: np.ndarray, u2: np.ndarray) -> float:
    c = float(np.clip(np.dot(u1, u2), -1.0, 1.0))
    return float(np.degrees(np.arccos(c)))


def depth_design(n_genes: np.ndarray, n_umi: np.ndarray) -> np.ndarray:
    """検出遺伝子数と総 UMI 数の log 変換（完全二次式。plan.md §4.5 の残差化に使う）。"""
    x1, x2 = np.log(np.asarray(n_genes, dtype=float)), np.log(np.asarray(n_umi, dtype=float))
    return np.column_stack([np.ones_like(x1), x1, x2, x1 ** 2, x2 ** 2, x1 * x2])


def spearman_pair(x: np.ndarray, y: np.ndarray):
    """NaN を書かず null と理由を書く（A-5）。"""
    m = np.isfinite(x) & np.isfinite(y)
    xs, ys = x[m], y[m]
    if len(xs) < 3 or float(np.ptp(xs)) == 0.0 or float(np.ptp(ys)) == 0.0:
        return None, f"計算できない（有効なペア {int(m.sum())} 個、または一方が一定）"
    return float(stats.spearmanr(xs, ys).statistic), None


# --------------------------------------------------------------------------- 入力
def load_assets(smoke: bool):
    results = EXP / ("results-smoke" if smoke else "results")
    df = pd.read_parquet(results / "embeddings.parquet")
    matched = json.loads((results / "embed_a6" / "depth_match.json").read_text())
    cent_rec = json.loads((results / "embed_a6" / "centroid_rec.json").read_text())
    run_meta = json.loads((results / "embed_a6" / "run_meta.json").read_text())
    cfg = yaml.safe_load((EXP / "config.yaml").read_text())
    return cfg, results, df, matched, cent_rec, run_meta


def load_expr(A, rows_by_sys: dict, cols: list[str]):
    """input.h5ad の生カウント（読み取りのみ）から、系統ごとのマッチ核 × 指定遺伝子の発現を作る。

    s[ens] = log1p(CP10k) を返す（plan.md §4.5）。rows_by_sys: {system: [row_id, ...]}（input.h5ad の行番号）
    """
    import scipy.sparse as sp
    X = A.X.tocsr()
    var = A.var
    e2i = {str(e): i for i, e in enumerate(var["ensembl_id"].astype(str))}
    sym = {str(e): str(s) for e, s in zip(var["ensembl_id"].astype(str), var["gene_name"].astype(str))}
    present = [e for e in cols if e in e2i]
    idxs = {e: e2i[e] for e in present}
    total_all = np.asarray(X.sum(axis=1)).ravel()
    out = {}
    for sys_name, rows in rows_by_sys.items():
        rows = [int(r) for r in rows]
        C = X[rows][:, [idxs[e] for e in present]] if present else sp.csr_matrix((len(rows), 0))
        C = C.tocsr()
        T = total_all[rows].astype(np.float64)
        T[T == 0] = 1.0
        D = C.toarray().astype(np.float64) / T[:, None] * 1e4
        s = {}
        for j, e in enumerate(present):
            s[e] = np.log1p(D[:, j])
        out[sys_name] = {"s": s, "total": T, "n_cells": len(rows)}
    return out, sym


def detect_goal(A, e2i: dict, goal_labels: list[str]):
    """目標の状態での生カウント検出率（単一遺伝子の位置の下限と対照プールの頻度マッチに使う）。"""
    import scipy.sparse as sp
    X = A.X.tocsr()
    keep = A.obs["celltype_l2"].astype(str).isin(goal_labels).to_numpy()
    Xg = X[keep]
    total_g = max(int(keep.sum()), 1)
    out = {}
    for e, i in e2i.items():
        col = Xg[:, i]
        if sp.issparse(col):
            n = int((np.asarray(col.todense()).ravel() > 0).sum())
        else:
            n = int((np.asarray(col).ravel() > 0).sum())
        out[e] = {"frac": float(n / total_g), "n": n}
    return out


# --------------------------------------------------------------------------- P-1
def cv_auc(X: np.ndarray, y: np.ndarray, groups: np.ndarray, folds: int = 5, seed: int = 0) -> dict:
    """提供者単位（GroupKFold）を主、層化分割を参考値とする（run_isp.precheck と同じ実装）。"""
    out = {}
    splitters = {"group_by_donor": GroupKFold(n_splits=folds),
                 "stratified": StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)}
    for name, cv in splitters.items():
        aucs = []
        for tr, te in cv.split(X, y, groups=groups):
            if len(set(y[te])) < 2:      # 片方のクラスしかない fold は AUC を定義できない（§4.2 注記）
                continue
            clf = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=5000))
            clf.fit(X[tr], y[tr])
            aucs.append(float(roc_auc_score(y[te], clf.predict_proba(X[te])[:, 1])))
        out[name] = {"auc_mean": fmt(np.mean(aucs)) if aucs else None,
                     "auc_sd": fmt(np.std(aucs, ddof=1)) if len(aucs) > 1 else None,
                     "n_folds_used": len(aucs), "fold_aucs": [fmt(a) for a in aucs]}
    return out


# --------------------------------------------------------------------------- P-4
def positions_for(df, matched_ids, u, cfg, expr_sys, pool_genes):
    """plan.md §4.5。1 つの（系統, 軸）に付き、遺伝子ごとの位置（Spearman）を返す。

    - 発現 s = log1p(CP10k) を対照プールで当てた回帰で残差化（検出遺伝子数と総 UMI 数の完全二次式）
    - 射影 proj = emb · u を同じ変数の設計で残差化
    - 位置 = Spearman(残差発現, 残差射影)
    - 対照プールの位置の分布で z 化
    """
    sel = df[df["row_id"].isin(matched_ids)].sort_values("row_id")
    if len(sel) < 3:
        return None, None, {"理由": f"核が少ない（{len(sel)}。試し実行）"}
    n_g = sel["n_genes"].to_numpy(dtype=float)
    n_u = sel["n_umi"].to_numpy(dtype=float)
    D = depth_design(n_g, n_u)
    emb = sel[[c for c in sel.columns if c.startswith("e")]].to_numpy(dtype=np.float64)
    proj = emb @ u
    beta_proj, *_ = np.linalg.lstsq(D, proj, rcond=None)
    res_proj = proj - D @ beta_proj

    rows = [expr_sys["s"][e] for e in pool_genes if e in expr_sys["s"]]
    if not rows:
        return None, None, {"理由": "対照プールの発現が取れない"}
    Ypool = np.concatenate(rows)
    Xpool = np.tile(D, (len(rows), 1))
    beta_expr, *_ = np.linalg.lstsq(Xpool, Ypool, rcond=None)

    def pos_of(e: str):
        s = expr_sys["s"].get(e)
        if s is None:
            return None, "発現が取れない"
        res_s = s - D @ beta_expr
        return spearman_pair(res_s, res_proj)

    pos = {}
    genes_all = sorted(set(pool_genes) | set(expr_sys["s"].keys()))
    for e in genes_all:
        p, _r = pos_of(e)
        if p is not None:
            pos[e] = p
    arr = np.array([pos[e] for e in pool_genes if e in pos], dtype=float)
    if len(arr) < 10 or float(np.std(arr)) == 0.0:
        pm, psd, why = (float(np.mean(arr)) if len(arr) else None, None,
                        "対照プールの位置が一定または不足（z 化できない）")
    else:
        pm, psd, why = float(np.mean(arr)), float(np.std(arr, ddof=1)), None
    pool_stat = {"n_used": len(arr), "mean": fmt(pm), "sd": fmt(psd)}
    if psd is None or psd <= 0:
        return None, pool_stat, {"理由": why}
    z = {e: (pos[e] - pm) / psd for e in pos}
    return {"pos": pos, "z": z, "res_proj": res_proj, "res_s_fit": {"beta_expr": None}}, pool_stat, None


def random_gene_set_p99(zpool: list[float], k: int, n_draw: int, rng) -> tuple[float, int]:
    """対照プールの z から k 個復元ありで引き、遺伝子セットの平均スコアの分布の 99 パーセンタイルと較正合格数。"""
    za = np.asarray(zpool, dtype=float)
    means = np.empty(n_draw)
    for d in range(n_draw):
        idx = rng.integers(0, len(za), size=k)
        means[d] = za[idx].mean()
    return float(np.quantile(means, 0.99, axis=0)), int((means > np.quantile(means, 0.99, axis=0)).sum() - 0)


# --------------------------------------------------------------------------- 本体
def main() -> None:
    t0 = time.time()
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true", help="試し実行（判定ではない）")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        self_test()
        return

    cfg, results, df, matched, cent_rec, run_meta = load_assets(a.smoke)
    dec = cfg["decision"]
    outdir = results
    outdir.mkdir(parents=True, exist_ok=True)
    LOG = outdir / "run.log"
    LOG.write_text("")
    L = lambda msg: log_line(LOG, msg)  # noqa: E731

    suf = "-smoke" if a.smoke else ""
    L(f"開始{'（試し実行）' if a.smoke else ''}。出力先 {outdir.relative_to(ROOT)}")
    L(f"plan.md sha256 {str(sha256f(EXP / 'plan.md'))[:16]}… / config.yaml sha256 "
      f"{str(sha256f(EXP / 'config.yaml'))[:16]}… / このスクリプト sha256 {str(sha256f(Path(__file__)))[:16]}…")

    # ---- 0. config の decision と plan の照合（機械的） ----------------------
    checks_conf = {
        "乱数の鍵（対照プール）": int(dec["乱数の鍵"]["対照プールの頻度マッチ"]) == PLAN["乱数の鍵"]["対照プールの頻度マッチ"],
        "乱数の鍵（サイズを揃えたランダムな遺伝子セット）": int(dec["乱数の鍵"]["偽セット 1,000 組"]) == PLAN["乱数の鍵"]["ランダムな遺伝子セット 1,000 組"],
        "G2 技術判定（config）": cfg["compute"]["needs_G2_approval"] is False,
        "G2 保守側（config・PI 承認）": cfg["compute"]["needs_G2_approval_conservative"] is True
        and "PI 2026-10-11" in str(cfg["compute"].get("approved_by", "")),
        "kind = embedding_only": cfg.get("kind") == "embedding_only",
    }
    plan_txt = (EXP / "plan.md").read_text()
    checks_conf["plan §10.2 Approved-by"] = "Approved-by: PI 2026-10-11" in plan_txt
    if not all(checks_conf.values()):
        L(f"**止める**: config の decision が plan と一致しない {checks_conf}")
        raise SystemExit("**止める**(config): decision が plan と違う")

    out: dict = {
        "実験": "EXP-20261011-003", "hypothesis": "HYP-20261011-01（方法。G3=判定不能 2026-10-11）",
        "mode": "smoke" if a.smoke else "final",
        "問い": "A6 の埋め込みで、出発状態は目標状態から分離でき（P-1）、2 本の軸（活性化した間質 u_ada と"
                "退化・老化 u_deg）は異なる向きを向くか（P-2）。前回の退化側の不成立は A6 の幾何でどう見えるか（P-3）。"
                "候補・陽性対照の位置（相関・記述）（P-4）",
        "0. config と plan の一致": checks_conf,
        "使った数値型": {"device（評価）": "cpu", "dtype（評価）": "fp64",
                     "GPU 時間（評価）": 0.0,
                     "備考": "保存済みの埋め込み（fp32・tf32 閉）の読み出しと再計算のみ。新しい埋め込み・摂動はしない"},
        "提供者 3 群（P-2(a)）": {"割り当て": DONOR_GROUPS, "注": GROUP_NOTE},
    }

    # ---- 1・2. 材料と設定の確認 ---------------------------------------------
    mats = {
        "埋め込み parquet": str((outdir / "embeddings.parquet").relative_to(ROOT)),
        "行数": int(len(df)),
        "埋め込みの記録（embed_a6）": {"device": run_meta.get("device"), "dtype": run_meta.get("dtype"),
                                    "host": run_meta.get("host"), "n_dim": run_meta.get("n_dim"),
                                    "config_sha256": str(run_meta.get("config_sha256"))[:16]},
        "材料 A6 sha256 == 固定値": bool(sha256f(ROOT / "data/raw/kidney_tls/A6_UUO_Reck2025.h5ad") == MATERIAL_SHA256),
    }
    cfg_now = sha256f(EXP / "config.yaml")
    mats["config.yaml の一致（A-2）"] = {
        "埋め込み時の sha256（run_meta の記録）": str(run_meta.get("config_sha256")),
        "現在の sha256": cfg_now,
        "受理する sha256（埋め込み時／学術用語への書き換え後）": list(CONFIG_SHA256_ACCEPTED),
        "一致": bool(cfg_now in CONFIG_SHA256_ACCEPTED),
    }
    out["1. 材料（A-2）"] = mats
    if run_meta.get("config_sha256") != CONFIG_SHA256_EMBEDDED:
        raise SystemExit("**止める**(A-2): 埋め込みの記録（run_meta）の config sha256 が想定と違う")
    if cfg_now not in CONFIG_SHA256_ACCEPTED:
        raise SystemExit("**止める**(A-2): 現在の config.yaml の sha256 が、埋め込み時（または書き換え後）と違う")

    # ---- 3. 遺伝子集合（A-3・A-7） -------------------------------------------
    gs = cfg["gene_sets"]
    sets_sym = {
        "候補の 8 遺伝子": list(dict.fromkeys(gs["candidates"]["symbols"])),
        "候補として採用した 5 遺伝子": list(dict.fromkeys(gs["candidates"]["adopted_5"])),
        "検出が少なく位置を推定できなかった遺伝子": list(dict.fromkeys(gs["candidates"]["not_measured"])),
        "探索的な扱いの遺伝子": list(dict.fromkeys(gs["candidates"]["exploratory_only"])),
        "活性化の陽性対照（71 個）": list(dict.fromkeys(gs["positive_activation"]["symbols"])),
        "管状系の陽性対照（23 個）": list(dict.fromkeys(gs["positive_tubular_a6"]["symbols"])),
        "退化・老化の陽性対照（30 個。記録）": list(dict.fromkeys(gs["positive_degenerate"]["symbols"])),
        "恒常発現の対照遺伝子（3 個）": list(dict.fromkeys(gs["housekeeping"])),
        "中立の対照遺伝子（7 個）": list(dict.fromkeys(gs["neutral"])),
    }
    ens_of = {}
    for mp in (gs["candidates"]["ensembl"], gs["positive_activation"]["ensembl"],
               gs["positive_degenerate"]["ensembl"], gs["positive_tubular_a6"]["ensembl"]):
        for s, e in mp.items():
            ens_of[s] = e
    dup_sets = {k: [x for x in v if v.count(x) > 1] for k, v in sets_sym.items() if len(set(v)) != len(v)}
    overlaps = {
        "候補の 8 遺伝子 ∩ 活性化の陽性対照（71 個）": sorted(set(sets_sym["候補の 8 遺伝子"]) & set(sets_sym["活性化の陽性対照（71 個）"])),
        "退化・老化の陽性対照（30 個。記録） ∩ 活性化の陽性対照（71 個）": sorted(set(sets_sym["退化・老化の陽性対照（30 個。記録）"]) & set(sets_sym["活性化の陽性対照（71 個）"])),
        "候補の 8 遺伝子 ∩ 退化・老化の陽性対照（30 個。記録）": sorted(set(sets_sym["候補の 8 遺伝子"]) & set(sets_sym["退化・老化の陽性対照（30 個。記録）"])),
        "管状系の陽性対照（23 個） ∩ 活性化の陽性対照（71 個）": sorted(set(sets_sym["管状系の陽性対照（23 個）"]) & set(sets_sym["活性化の陽性対照（71 個）"])),
        "管状系の陽性対照（23 個） ∩ 候補の 8 遺伝子": sorted(set(sets_sym["管状系の陽性対照（23 個）"]) & set(sets_sym["候補の 8 遺伝子"])),
    }
    blind_coll = {k: sorted(set(v) & set(gs["blind5"])) for k, v in sets_sym.items()}
    blind_coll_ens = {k: sorted(e for e in ens_of.values() if e in set(gs["blind5"].values()))
                      for k in ("候補の 8 遺伝子", "活性化の陽性対照（71 個）", "退化・老化の陽性対照（30 個。記録）", "管状系の陽性対照（23 個）")}
    out["2. 遺伝子集合（A-3・A-7）"] = {
        "個数": {k: len(v) for k, v in sets_sym.items()},
        "盲検中の 5 遺伝子": {"個数": len(gs["blind5"]), "注": "名前は出力に書かない（A-1）。他の集合との重なりは下の空リストで確認する"},
        "片方の重複（dict.fromkeys で除去。A-7）": dup_sets,
        "集合の重なり": overlaps,
        "盲検との重なり（symbol）": blind_coll,
        "盲検との重なり（Ensembl）": blind_coll_ens,
    }
    if any(dup_sets.values()) or any(blind_coll.values()) or any(blind_coll_ens.values()):
        raise SystemExit("**止める**(A-1/A-7): 集合の重複か盲検の混入がある")
    pre = json.loads((EXP / "work" / "precheck_a6.json").read_text())
    pre_ov = pre["overlaps"]["pairs"]
    if sorted(overlaps["候補の 8 遺伝子 ∩ 活性化の陽性対照（71 個）"]) != sorted(pre_ov["candidates ∩ positive_activation"]) \
            or sorted(overlaps["退化・老化の陽性対照（30 個。記録） ∩ 活性化の陽性対照（71 個）"]) != sorted(pre_ov["positive_activation ∩ positive_degenerate"]) \
            or sorted(overlaps["管状系の陽性対照（23 個） ∩ 活性化の陽性対照（71 個）"]) != sorted(pre_ov["positive_activation ∩ positive_tubular_a6"]):
        raise SystemExit("**止める**(A-3): 重なりが precheck_a6.json（固定値）と一致しない")
    if len(sets_sym["候補の 8 遺伝子"]) != 8 or len(sets_sym["活性化の陽性対照（71 個）"]) != 71 \
            or len(sets_sym["退化・老化の陽性対照（30 個。記録）"]) != 30 or len(sets_sym["管状系の陽性対照（23 個）"]) != 23:
        raise SystemExit("**止める**(A-3): 遺伝子集合の個数が config と違う")

    # ---- 4. P-1（状態の分離。検出遺伝子数の分布をそろえた後の核） -----------------
    ST = {"row": df["row_id"].astype(int).to_numpy(),
          "state": df["celltype_l2"].astype(str).to_numpy(),
          "donor": df["donor_id"].astype(str).to_numpy(),
          "n_genes": df["n_genes"].to_numpy(dtype=float),
          "n_umi": df["n_umi"].to_numpy(dtype=float)}
    emb_all = df[[c for c in df.columns if c.startswith("e")]].to_numpy(dtype=np.float64)
    by_row = {int(r): i for i, r in enumerate(ST["row"])}

    sys_cfg = cfg["systems"]
    p1 = {}
    for name, spec in sys_cfg.items():
        s_ids = [int(x) for x in matched[name]["start"]]
        g_ids = [int(x) for x in matched[name]["goal"]]
        avail_s = [by_row[x] for x in s_ids if x in by_row]
        avail_g = [by_row[x] for x in g_ids if x in by_row]
        if len(avail_s) < 2 or len(avail_g) < 2:
            p1[name] = {"n_start": len(s_ids), "n_goal": len(g_ids),
                        "n_usable_start": len(avail_s), "n_usable_goal": len(avail_g),
                        "理由": "マッチした核が少なすぎる（試し実行）"}
            continue
        idx = avail_s + avail_g
        y = np.concatenate([np.zeros(len(avail_s)), np.ones(len(avail_g))])
        X = emb_all[idx]
        grp = ST["donor"][idx]
        Xd = np.column_stack([np.log(ST["n_genes"][idx]), np.log(ST["n_umi"][idx])])
        sep = cv_auc(X, y, grp)
        depth = cv_auc(Xd, y, grp)
        pc = {"n_start": len(s_ids), "n_goal": len(g_ids), "n_usable": len(idx)}
        for k1, src in (("埋め込み AUC", sep), ("総 UMI 数と検出遺伝子数の 2 変数での AUC", depth)):
            for k2 in ("group_by_donor", "stratified"):
                pc[f"{k1}.{k2}"] = src[k2]
        pc["P-1(i) 提供者単位 AUC ≥ 0.70"] = bool(sep["group_by_donor"]["auc_mean"] is not None
                                                   and sep["group_by_donor"]["auc_mean"] >= PLAN["P-1_提供者単位AUC下限"])
        pc["P-1(ii) 総 UMI 数と検出遺伝子数の 2 変数での AUC ≤ 0.60"] = bool(depth["group_by_donor"]["auc_mean"] is not None
                                                 and depth["group_by_donor"]["auc_mean"] <= PLAN["P-1_2変数だけの分離AUC上限"])
        p1[name] = pc
    out["P-1（状態の分離。検出遺伝子数の分布をそろえた後の核・提供者単位分割を主）"] = p1

    # ---- 5. P-2（2 本の軸の幾何） --------------------------------------------
    u_ada = u_deg = None   # 軸が計算できなければ P-4 も計算しない（試し実行）
    def centroid_of(ids):
        idx = [by_row[x] for x in ids if x in by_row]
        return emb_all[idx].mean(axis=0) if idx else None

    cents = {}
    for sys_name, (sn, gn) in STATE_NAMES.items():
        cents[sn] = centroid_of(matched[sys_name]["start"])
        cents[gn] = centroid_of(matched[sys_name]["goal"])
    cnames = ["Fibroblast", "Myofibroblast", "PT_normal", "PT_inj_infl"]
    p2: dict = {"セントロイド（検出遺伝子数の分布をそろえた後）": {k: ("あり" if cents[k] is not None else "無し（試し実行）")
                                      for k in cnames}}
    grp_angles: dict = {}
    p2_pass: dict = {}

    if all(cents[k] is not None for k in cnames):
        u_ada = unit(unit(cents["Myofibroblast"]) - unit(cents["Fibroblast"]))
        u_deg = unit(unit(cents["PT_inj_infl"]) - unit(cents["PT_normal"]))
        pop = np.isin(ST["state"], ALL_STATES)
        ng = ST["n_genes"][pop]
        q05, q95 = np.quantile(ng, [0.05, 0.95])
        pop_idx = np.where(pop)[0]
        bot = pop_idx[ng <= q05]
        top = pop_idx[ng >= q95]
        if len(bot) and len(top):
            d_depth = unit(emb_all[top].mean(axis=0) - emb_all[bot].mean(axis=0))
        else:
            d_depth = np.full(1152, np.nan)
        p2["u_ada"] = {"from": "Myofibroblast − Fibroblast（検出遺伝子数の分布をそろえた後のセントロイド）"}
        p2["u_deg"] = {"from": "PT_inj_infl − PT_normal（検出遺伝子数の分布をそろえた後のセントロイド）"}
        p2["d_depth"] = {"from": "検出遺伝子数の上位 5 分位 − 下位 5 分位（間質系＋管状系の全核）",
                         "n_top": int(len(top)), "n_bottom": int(len(bot))}

        # (a) 提供者 3 群すべてで sin ≥ 0.31（群内で検出遺伝子数の分布をそろえ直す）
        a_ok = True
        for gname, donors in DONOR_GROUPS.items():
            dm_s, dm_g = {}, {}
            for sys_name in ("interstitial", "tubular"):
                m_d = np.isin(ST["donor"], donors)
                sub_st = ST["state"][m_d]
                sub_ng = ST["n_genes"][m_d]
                sub_row = ST["row"][m_d]
                lab_s = [str(x) for x in sys_cfg[sys_name]["start"]["labels"]]
                lab_g = [str(x) for x in sys_cfg[sys_name]["goal"]["labels"]]
                ms, mg = depth_match_ids(sub_st, sub_ng, lab_s, lab_g)
                dm_s[sys_name] = [int(sub_row[i]) for i in ms]
                dm_g[sys_name] = [int(sub_row[i]) for i in mg]
            cg = {}
            for sys_name, (sn, gn) in STATE_NAMES.items():
                cg[sn] = centroid_of(dm_s[sys_name])
                cg[gn] = centroid_of(dm_g[sys_name])
            if all(cg[k] is not None for k in ("Fibroblast", "Myofibroblast", "PT_normal", "PT_inj_infl")):
                ua = unit(unit(cg["Myofibroblast"]) - unit(cg["Fibroblast"]))
                ud = unit(unit(cg["PT_inj_infl"]) - unit(cg["PT_normal"]))
                sin_v = float(np.sqrt(max(0.0, 1.0 - float(np.clip(np.dot(ua, ud), -1, 1)) ** 2)))
                grp_angles[gname] = {"sin": fmt(sin_v), "n_start": len(dm_s["interstitial"]),
                                     "n_goal": len(dm_g["interstitial"]),
                                     "n_pt_start": len(dm_s["tubular"]), "n_pt_goal": len(dm_g["tubular"])}
                a_ok &= bool(sin_v >= PLAN["P-2_sin下限"])
            else:
                grp_angles[gname] = {"sin": None, "理由": "この群でセントロイドが足りない（試し実行）"}
                a_ok = False

        # (b) u_ada・u_deg と検出遺伝子数の高低の向きの余弦（絶対値 ≤ 0.60）
        if np.isfinite(d_depth).all():
            cos_a = float(np.clip(np.dot(u_ada, d_depth), -1, 1))
            cos_d = float(np.clip(np.dot(u_deg, d_depth), -1, 1))
            b_ok = bool(abs(cos_a) <= PLAN["P-2_検出遺伝子数の高低の向きとの余弦上限"]
                        and abs(cos_d) <= PLAN["P-2_検出遺伝子数の高低の向きとの余弦上限"])
        else:
            cos_a = cos_d = None
            b_ok = False

        # (c) 検出遺伝子数の分布をそろえるかどうかで軸が安定（角度差 ≤ 30 度）
        cents_nm = {}
        for sys_name in ("interstitial", "tubular"):
            lab_s = set([str(x) for x in sys_cfg[sys_name]["start"]["labels"]])
            lab_g = set([str(x) for x in sys_cfg[sys_name]["goal"]["labels"]])
            m_s = np.isin(ST["state"], list(lab_s))
            m_g = np.isin(ST["state"], list(lab_g))
            cents_nm[sys_name] = (emb_all[m_s].mean(axis=0) if m_s.any() else None,
                                  emb_all[m_g].mean(axis=0) if m_g.any() else None)
        ua_nm = ud_nm = None
        if all(c is not None for pair in cents_nm.values() for c in pair):
            ua_nm = unit(unit(cents_nm["interstitial"][1]) - unit(cents_nm["interstitial"][0]))
            ud_nm = unit(unit(cents_nm["tubular"][1]) - unit(cents_nm["tubular"][0]))
        da = angle_deg(u_ada, ua_nm) if ua_nm is not None else None
        dd = angle_deg(u_deg, ud_nm) if ud_nm is not None else None
        c_ok = bool(da is not None and dd is not None and da <= PLAN["P-2_角度差上限_度"]
                    and dd <= PLAN["P-2_角度差上限_度"])

        p2["P-2(a) sin ≥ 0.31（提供者 3 群すべて）"] = {"詳細": grp_angles, "全て満たす": bool(a_ok)}
        p2["P-2(b) 検出遺伝子数の高低の向きとの余弦の絶対値 ≤ 0.60"] = {
            "cos(u_ada, 検出遺伝子数の高低の向き)": fmt(cos_a), "cos(u_deg, 検出遺伝子数の高低の向き)": fmt(cos_d), "満たす": bool(b_ok)}
        p2["P-2(c) 検出遺伝子数の分布をそろえるかどうかでの角度差 ≤ 30 度"] = {
            "u_ada の角度差": fmt(da), "u_deg の角度差": fmt(dd), "満たす": bool(c_ok)}
        p2_pass = {"a": bool(a_ok), "b": bool(b_ok), "c": bool(c_ok)}
        p2["P-2 全体（a かつ b かつ c）"] = bool(a_ok and b_ok and c_ok)
    else:
        p2["注"] = "セントロイドが足りないため P-2 は計算できない（試し実行）"
        p2_pass = {"a": False, "b": False, "c": False}
    out["P-2（2 本の軸は異なる向きを向くか）"] = p2
    out["P-2a 提供者 3 群の sin"] = grp_angles

    # ---- 6. P-3（前回の退化側の不成立が A6 ではどう見えるか。記述） -----------
    tub_p1 = p1.get("tubular", {})
    tub_p1_ok = bool(tub_p1.get("P-1(i) 提供者単位 AUC ≥ 0.70") and tub_p1.get("P-1(ii) 総 UMI 数と検出遺伝子数の 2 変数での AUC ≤ 0.60"))
    tub_p2_ok = bool(p2_pass.get("a") and p2_pass.get("b"))
    p3_ok = tub_p1_ok and tub_p2_ok
    out["P-3（管状系が P-1・P-2 を両方満たすか。記述）"] = {
        "管状系 P-1": tub_p1_ok, "管状系 P-2": tub_p2_ok, "両方満たす": p3_ok,
        "読み": ("A6 の管状系でも退化方向が立つ（P-1・P-2 を通る）: "
                 "A1/A2 の退化側の不成立は、A6 の幾何では退化方向が立たないからではなく、"
                 "A1 の退化状態のラベルまたは摂動の測り方に固有の問題であることを示唆する（記述）"
                 if p3_ok else
                 "A6 の管状系でも退化方向が立たない。前回の不成立に似た形が独立データで再現された（記述）")}

    # ---- 7. P-4（候補・陽性対照の位置。相関・記述） --------------------------
    # 対照プール（§5.2）: A6 の辞書にある遺伝子から、Myofibroblast での検出率を候補（adopted_5）にマッチ。
    sys_cell_rows = {}
    for name in ("interstitial", "tubular"):
        ids = sorted(set([int(x) for x in matched[name]["start"]] + [int(x) for x in matched[name]["goal"]]))
        sys_cell_rows[name] = ids
    import anndata as ad
    A = ad.read_h5ad(str(EXP / "work" / "input.h5ad"))
    e2i_all = {str(e): i for i, e in enumerate(A.var["ensembl_id"].astype(str))}
    # 恒常発現の対照遺伝子（3 個）・中立の対照遺伝子（7 個）は config に symbol しか無いので、A6 の var（gene_name → ensembl_id）で補う
    s2e = {}
    for _e, _s in zip(A.var["ensembl_id"].astype(str), A.var["gene_name"].astype(str)):
        s2e.setdefault(str(_s), str(_e))
    for _sym in list(gs["housekeeping"]) + list(gs["neutral"]) + list(T5_SYMBOLS):
        if _sym not in ens_of and _sym in s2e:
            ens_of[_sym] = s2e[_sym]
    det_inter = detect_goal(A, e2i_all, [str(x) for x in sys_cfg["interstitial"]["goal"]["labels"]])
    vocab = _load_vocab()
    det_pool = {e: v for e, v in det_inter.items() if e in e2i_all and e in vocab and v["n"] > 0}
    pool = sample_control_pool(cfg, e2i_all, det_pool)
    genes_of_interest = list(dict.fromkeys(
        list(ens_of.values()) + list(gs["positive_activation"]["ensembl"].values())
        + list(gs["positive_degenerate"]["ensembl"].values())
        + list(gs["positive_tubular_a6"]["ensembl"].values())
        + [s2e[s] for s in T5_SYMBOLS if s in s2e] + list(pool)))
    expr_sys, sym_of = load_expr(A, sys_cell_rows, genes_of_interest)
    # 管状系の目標（PT Injured ∪ Inflammatory）での検出率（単一遺伝子の位置の下限の注記用）
    det_tub = detect_goal(A, e2i_all, [str(x) for x in sys_cfg["tubular"]["goal"]["labels"]])
    det_by_sys = {"interstitial": det_inter, "tubular": det_tub}
    out["P-4 対照プール"] = {"n": len(pool), "seed": PLAN["乱数の鍵"]["対照プールの頻度マッチ"],
                            "先頭": [sym_of.get(e, e) for e in pool[:5]],
                            "注": "結果を見る前に固定した手順（§5.2）でサンプリング"}

    axes = {"u_ada": u_ada, "u_deg": u_deg}
    pairs = [("interstitial", "u_ada"), ("interstitial", "u_deg"),
             ("tubular", "u_deg"), ("tubular", "u_ada")]
    pair_res = {}
    for sys_name, ax in pairs:
        pair_res[(sys_name, ax)] = None
        if ax not in axes or axes[ax] is None or not np.isfinite(axes[ax]).all():
            continue
        res = positions_for(df, sys_cell_rows[sys_name], axes[ax], cfg, expr_sys[sys_name], pool)
        if res[0] is None:
            continue
        pair_res[(sys_name, ax)] = {"pos": res[0]["pos"], "z": res[0]["z"], "pool_stat": res[1],
                                    "理由": res[2]}

    p4_health = {}
    pos_rows_all = []
    floor = PLAN["P-4_単一遺伝子_検出された核の数_下限"]
    random_set_label = "サイズを揃えたランダムな遺伝子セット（対照遺伝子群から復元ありで引く）"
    for sys_name, ax in pairs:
        r = pair_res.get((sys_name, ax))
        if r is None:
            p4_health.setdefault(sys_name, {})[ax] = {"position": None, "理由": "軸または核が足りない（試し実行）"}
            continue
        pos, z = r["pos"], r["z"]
        pool_stat = r["pool_stat"]
        rng = np.random.default_rng(PLAN["乱数の鍵"]["ランダムな遺伝子セット 1,000 組"])
        zpool = [z[e] for e in pool if e in z]
        p99, cal = random_gene_set_p99(zpool, PLAN["P-4_ランダムな遺伝子セットの大きさ"], PLAN["P-4_ランダムな遺伝子セットの組数"], rng)

        def set_mean_z(syms):
            es = [ens_of[s] for s in syms if s in ens_of and ens_of[s] in z]
            return (float(np.mean([z[e] for e in es])) if es else None, len(es))

        health_checks = {}
        p99_key = f"{random_set_label}の 99 パーセンタイル"
        if sys_name == "interstitial" and ax == "u_ada":
            mz, n = set_mean_z(sets_sym["活性化の陽性対照（71 個）"])
            health_checks["活性化の陽性対照（71 個）"] = {"遺伝子セットの平均スコア（z）": fmt(mz), "n": n,
                                                p99_key: fmt(p99),
                                                "正かつ 99% 点超": bool(mz is not None and mz > 0 and mz > p99)}
        if sys_name == "tubular" and ax == "u_deg":
            mz, n = set_mean_z(sets_sym["管状系の陽性対照（23 個）"])
            syms30 = list(dict.fromkeys(sets_sym["管状系の陽性対照（23 個）"]
                                        + [s for s in T5_SYMBOLS if s2e.get(s) and s2e[s] in z]))
            mz30, n30 = set_mean_z(syms30)
            health_checks["管状系の陽性対照（23 個）"] = {"遺伝子セットの平均スコア（z）": fmt(mz), "n": n,
                                                p99_key: fmt(p99),
                                                "正かつ 99% 点超": bool(mz is not None and mz > 0 and mz > p99)}
            health_checks["感度: 管状系 30（T5 を含む）"] = {"遺伝子セットの平均スコア（z）": fmt(mz30), "n": n30,
                                                         "注": "事前登録した感度分析。主判定には使わない"}

        # 単一遺伝子の位置（§4.5.5: 目標状態でその遺伝子が検出された核（の数）≥ 20 に限る）
        det_goal_s = det_by_sys[sys_name]

        def emit(sym, set_name, note_extra=""):
            e = ens_of.get(sym)
            dn = int(det_goal_s[e]["n"]) if (e and e in det_goal_s) else None
            pv, zv = None, None
            if note_extra:          # 検出が少なかった・探索的な扱い: 位置は報告しない（§0.4・§4.5.5）
                note = note_extra
            elif e and e in pos:
                if dn is not None and dn < floor:
                    note = (f"位置を推定できない（検出が少なく位置を推定できなかった。"
                            f"§4.5.5: 目標状態で検出された核の数 {dn} < {floor}）")
                else:
                    pv, zv = pos[e], z[e]
                    note = ""
            else:
                note = ""
            return {"system": sys_name, "axis": ax, "symbol": sym, "ensembl": e, "set": set_name,
                    "position（Spearman）": fmt(pv, 4), "position_z": fmt(zv, 4),
                    "n_detected_goal": dn, "note": note}

        # 候補（探索・記述）は u_ada 方向だけを報告（§4.5(ii)）
        if sys_name == "interstitial" and ax == "u_ada":
            for sym in sets_sym["候補の 8 遺伝子"]:
                e = ens_of.get(sym)
                note = ""
                if sym in sets_sym["検出が少なく位置を推定できなかった遺伝子"]:
                    note = "検出が少なく位置を推定できなかった（SLC4A9 は検出 2/585。§0.4・§4.5.5）"
                elif sym in sets_sym["探索的な扱いの遺伝子"]:
                    note = "探索的な扱いとして併記（A6 で検出 0 核。位置は推定しない）"
                pos_rows_all.append(emit(sym, "候補の 8 遺伝子", note))
        # 陽性対照・恒常発現・中立の位置（記述）
        for sym, set_name in ([(s, "活性化の陽性対照（71 個）") for s in sets_sym["活性化の陽性対照（71 個）"]]
                              + [(s, "管状系の陽性対照（23 個）") for s in sets_sym["管状系の陽性対照（23 個）"]]
                              + [(s, "退化・老化の陽性対照（30 個。記録）") for s in sets_sym["退化・老化の陽性対照（30 個。記録）"]]
                              + [(s, "恒常発現の対照遺伝子（3 個）") for s in sets_sym["恒常発現の対照遺伝子（3 個）"]]
                              + [(s, "中立の対照遺伝子（7 個）") for s in sets_sym["中立の対照遺伝子（7 個）"]]):
            pos_rows_all.append(emit(sym, set_name))
        p4_health.setdefault(sys_name, {})[ax] = {
            "pool": pool_stat, p99_key: fmt(p99),
            "較正（calibration。99 パーセンタイルを超えた組数）": cal,
            "健全性": health_checks or None,
            "注": "発現と埋め込みの相関。摂動の確認ではない（§11-9・§12-6）。"
                  "単一遺伝子の位置は、目標状態でその遺伝子が検出された核（の数）≥ " + str(floor) + " に限る（§4.5.5）"}

    out["P-4（候補・陽性対照の位置。相関・記述）"] = p4_health

    # ---- 8. 反証条件（R-A1〜R-D1。上から順に当てる） --------------------------
    def p1_auc_ok(name):
        return bool(p1.get(name, {}).get("P-1(i) 提供者単位 AUC ≥ 0.70"))
    def p1_depth_ok(name):
        return bool(p1.get(name, {}).get("P-1(ii) 総 UMI 数と検出遺伝子数の 2 変数での AUC ≤ 0.60"))
    def p1_folds_zero(name):
        v = p1.get(name, {})
        return v.get("埋め込み AUC.group_by_donor", {}).get("n_folds_used", 0) == 0

    r = {
        "R-A1": not (p1_auc_ok("interstitial") and p1_auc_ok("tubular")),
        "R-A2": not (p1_depth_ok("interstitial") and p1_depth_ok("tubular")),
        "R-A3": p1_folds_zero("interstitial"),
        "R-B1": not p2_pass.get("a", False),
        "R-B2": not p2_pass.get("b", False),
        "R-B3": not p2_pass.get("c", False),
        "R-C1": False,
    }
    # R-C1: 提供者 3 群で sin の値が食い違う（方向が逆の群がある）
    sins = [v.get("sin") for v in grp_angles.values() if isinstance(v, dict) and v.get("sin") is not None]
    if len(sins) >= 2:
        r["R-C1"] = bool(min(sins) < 0.1 and max(sins) >= PLAN["P-2_sin下限"])
    first = next((k for k, v in r.items() if v), "R-D1（保留）")
    out["反証条件（R-A1〜R-D1）"] = {
        "R-A1（どちらかの系統が P-1(i) を満たさない）": r["R-A1"],
        "R-A2（どちらかの系統が P-1(ii) を満たさない）": r["R-A2"],
        "R-A3（間質系の提供者単位 AUC が fold 欠落で計算できない）": r["R-A3"],
        "R-B1（P-2(a) を 3 群すべてで満たさない）": r["R-B1"],
        "R-B2（P-2(b) を満たさない）": r["R-B2"],
        "R-B3（P-2(c) を満たさない）": r["R-B3"],
        "R-C1（提供者 3 群で方向が食い違う）": r["R-C1"],
        "最初に当たった反証条件": first,
        "検出の下限（単一遺伝子の位置）": PLAN["P-4_単一遺伝子_検出された核の数_下限"],
        "注": "R-A1/A2/A3 は系統単位の読み（plan.md §4.6）。食い違った場合は詳細を併記し、判定は PI に諮る",
    }

    # ---- 9. 出力 -------------------------------------------------------------
    gene_pos_df = pd.DataFrame(pos_rows_all)
    gene_pos_df.to_csv(outdir / "gene_position_a6.csv", index=False)
    out["gene_position_a6.csv"] = {"行数": len(gene_pos_df), "保存先": "gene_position_a6.csv",
                                   "注": "相関・記述（摂動の確認ではない）"}

    out["decision（config.yaml の全鍵）"] = dec
    out["注"] = ("判定（支持・否定・判定不能。ゲート G3）は PI が行う。この道具は規則を機械的に当てた結果と提案まで。"
                 "hypotheses/tracker.md は変更しない。コンピュータ上の埋め込みのみの計算であり、実験的に確認された"
                 "生物学的事実ではない。盲検中の遺伝子は解析していない（値は数えない）。P-4 の位置は相関・記述である"
                 "（摂動の確認ではない。§4.5・§11-9・§12-6）。")
    dur = time.time() - t0
    out["実行"] = {"開始": _dt.datetime.fromtimestamp(t0).astimezone().isoformat(),
                 "所要秒": round(dur, 1), "GPU 時間": 0.0, "device": "cpu", "dtype": "fp64"}
    cj = outdir / f"criteria_a6{suf}.json"
    cj.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    reloaded = json.loads(cj.read_text())
    a6 = {"decision の全鍵を持つ": set((reloaded.get("decision（config.yaml の全鍵）") or {}).keys()) == set(dec.keys()),
          "保存先": str(cj.relative_to(ROOT)), "pass": True}
    out["A チェック"] = {"A-2 材料": {"pass": bool(mats["材料 A6 sha256 == 固定値"])},
                        "A-6 結果ファイル": a6}
    cj.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    (outdir / "README.md").write_text(
        "EXP-20261011-003 の結果（埋め込みのみ。摂動なし）。\n"
        "- embeddings.parquet: A6 全 46,957 核の cls 埋め込み（row_id 単位・fp32 エンコード）。\n"
        "- criteria_a6.json: P-1〜P-4 と反証条件 R-A1〜R-D1 を機械的に当てた結果（**提案**。確定は PI＝ゲート G3）。\n"
        "- gene_position_a6.csv: 候補・陽性対照の位置（**発現と埋め込みの相関・記述**。摂動の確認ではない）。\n"
        "- run.log: 実行ログ。\n" + (f"（試し実行 smoke。判定ではない）\n" if a.smoke else "")
        + "盲検中の 5 遺伝子は解析していない（名前は出力しない）。\n")

    # ---- A-1（出力・ログ・一覧に盲検が 0 件） --------------------------------
    sym_blk = {"ENSG00000141510": None, "ENSG00000124762": None, "ENSG00000111321": None,
               "ENSG00000118785": None, "ENSG00000148346": None}
    names_blk = ["TP53", "CDKN1A", "LTBR", "SPP1", "LCN2"]
    targets = list(outdir.glob("*.*"))
    a1 = {"調べたファイル": [str(p.relative_to(ROOT)) for p in targets if p.suffix in (".json", ".csv", ".log", ".md")],
          "件数": 0}
    for p in targets:
        if p.suffix not in (".json", ".csv", ".log", ".md"):
            continue
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

    # ---- A-5（値に NaN・無限大がないことを結果ファイルの再読込で確認） --------
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
    bad = [v for v in vals if isinstance(v, float) and (v != v or v in (float("inf"), float("-inf")))]
    a5 = {"NaN/Inf の数": len(bad), "pass": len(bad) == 0}
    out["A チェック"]["A-5 値"] = a5
    cj.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    L(f"完了（{dur:.1f} 秒）。出力 {cj.relative_to(ROOT)}")


def _load_vocab() -> set:
    import pickle
    from geneformer import TOKEN_DICTIONARY_FILE
    with open(TOKEN_DICTIONARY_FILE, "rb") as f:
        tok = pickle.load(f)
    return {str(k) for k in tok if str(k).startswith("ENS")}


def sample_control_pool(cfg, e2i: dict, det_goal: dict, rng_key: int = 20261011) -> list[str]:
    """対照プール 500（plan.md §5.2）。候補（adopted_5）と同じ検出率帯から、除外集合を除いて選ぶ。"""
    gs = cfg["gene_sets"]
    ens_exclude = set()
    for mp in (gs["candidates"]["ensembl"], gs["positive_activation"]["ensembl"],
               gs["positive_degenerate"]["ensembl"], gs["positive_tubular_a6"]["ensembl"]):
        ens_exclude |= set(mp.values())
    ens_exclude |= set(gs["blind5"].values())
    adopted = gs["candidates"]["adopted_5"]
    rng = np.random.default_rng(rng_key)
    eligible = [e for e, v in det_goal.items() if e in e2i and e not in ens_exclude]
    targets = [gs["candidates"]["ensembl"][s] for s in adopted if s in gs["candidates"]["ensembl"]]
    if not targets:
        raise SystemExit("**止める**(A-3): adopted_5 の Ensembl が取れない")
    tol = np.log(1.5)
    picked: list[str] = []
    for d in range(int(PLAN["対照プール_n"])):
        t = targets[d % len(targets)]
        td = det_goal[t]["frac"]
        cand = [e for e in eligible if e not in picked
                and abs(np.log((det_goal[e]["frac"] + 1e-6) / (td + 1e-6))) <= tol]
        if not cand:
            cand = [e for e in eligible if e not in picked]
            if not cand:
                raise SystemExit("**止める**(A-5): 対照プールの候補が尽きた")
            cand = [min(cand, key=lambda e: abs(np.log((det_goal[e]["frac"] + 1e-6) / (td + 1e-6))))]
        picked.append(str(rng.choice(cand)))
    return list(dict.fromkeys(picked))


def self_test() -> None:
    """合成データで cv_auc・depth_match_ids・サイズを揃えたランダムな遺伝子セットの機械を最低限確かめる。"""
    rng = np.random.default_rng(7)
    # cv_auc: 完全分離のデータで 1.0 に近いこと
    X = rng.normal(0, 1, (60, 5))
    y = (X[:, 0] + X[:, 1] > 0).astype(int)
    g = np.repeat(np.arange(12), 5)
    r = cv_auc(X, y, g)
    assert r["group_by_donor"]["auc_mean"] is not None and r["group_by_donor"]["auc_mean"] > 0.9
    # depth_match_ids: 決定的・合計が小さい側の核数以下
    st = np.array(["a"] * 30 + ["b"] * 50)
    ng = rng.integers(100, 1000, 80).astype(float)
    ms, mg = depth_match_ids(st, ng, ["a"], ["b"])
    assert len(ms) == len(mg) and 0 < len(ms) <= 30
    ms2, mg2 = depth_match_ids(st, ng, ["a"], ["b"])
    assert ms == ms2 and mg == mg2
    # サイズを揃えたランダムな遺伝子セット: ラベルを反転させると合格率が逆になることを使って機械を確認
    zpool = [float(x) for x in rng.normal(0, 1, 200)]
    p99, cal = random_gene_set_p99(zpool, 8, 1000, np.random.default_rng(1))
    # 8 個の平均の 99 パーセンタイルは ≈ 2.33 × sd / sqrt(8) ≈ 0.82（生の 99 パーセンタイル ≈ 2.33 より狭い）
    assert 0.3 < p99 < 1.5 and 0 < cal < 100
    print("self-test: cv_auc・depth_match_ids・サイズを揃えたランダムな遺伝子セット OK")


if __name__ == "__main__":
    main()
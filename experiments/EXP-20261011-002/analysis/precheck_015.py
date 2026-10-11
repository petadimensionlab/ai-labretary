#!/usr/bin/env python3
"""EXP-20261011-002 の事前確認（結果を見る前に実行。判定の値は一切計算しない）。

- 集合の重なり（総当たり）・盲検 0 件
- A1 の注釈テーブル（発現細胞数 n_exp・入力の長さ ll_med・元の位置 lp_med。摂動の結果に依存しない）
- マッチ偽セットが引けるか（PC 遺伝子ごとに完全な注釈を持つプール近傍 ≥ 5 個）
- 各 PC 版の検出率（n_exp / 200。判定はしない。記録のみ）
- 材料の sha256・トークン化済み入力の複製の照合

**z 値・セット平均・偽セットの閾値は計算しない。** 主判定の値は evaluate_015.py が固定後に 1 回計算する。

走らせ方:
    .venv/bin/python experiments/EXP-20261011-002/analysis/precheck_015.py
出力: work/precheck_015.json
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
EXP = ROOT / "experiments/EXP-20261011-002"
A1_SRC = ROOT / "experiments/EXP-20261004-001"
A1_NEW = ROOT / "experiments/EXP-20261011-001"
sys.path.insert(0, str(EXP / "analysis"))
sys.path.insert(0, str(ROOT / "experiments/EXP-20261004-001/analysis"))
sys.path.insert(0, str(ROOT / "experiments/EXP-20261001-001/analysis"))
sys.path.insert(0, str(ROOT / "scripts/analysis"))
import evaluate_009 as ev9                                   # noqa: E402
from evaluate_015 import gene_sets, sha256f, annotation_table  # noqa: E402

SEEDS = [81, 82, 83]
N_DRAW = 1000


def load_a1_annotations(cfg) -> pd.DataFrame:
    """A1（既存 626＋新規 30）の注釈だけを、摂動の結果に依存しない形で計算する。
    load_cells → 各遺伝子の n_exp（発現細胞数）・ll_med・lp_med の 3 群平均。z は計算しない。"""
    tok = str(ROOT / cfg["evaluation"]["a1"]["tokenized_dataset"])
    cfg1 = {"data": {"tokenized_dataset": tok, "state_key": "fib_anno_state"},
            "states": {"start_state": "FIB_MFIB", "goal_state": "aStr", "alt_states": ["dFIB", "dMFIB"]},
            "output": {"results_dir": "results"}}
    ev9.BLOCKED = {e: e for e in cfg["blocked_genes"]}
    frames = []
    for tag, E, res in (("既存626", A1_SRC, A1_SRC / "results"), ("新規30", A1_NEW, A1_NEW / "results")):
        ev9.E = E
        c, _ = ev9.load_cells(cfg1, "results")
        c = c[c.seed.isin(SEEDS)]
        # gene_stats_stratified と同じ定義（発現細胞数・入力の長さ・元の位置）で注釈を作る
        rows = []
        for gene, tg in c.groupby("gene"):
            te = tg[tg.expressed]
            rows.append({"gene": gene, "set": tag,
                         "n_cells": int(len(tg)), "n_exp": int(len(te)),
                         "ll_med": float(tg.ll.median()),
                         "lp_med": float(te.lp.median()) if len(te) else float("nan")})
        frames.append(pd.DataFrame(rows))
    ann = pd.concat(frames, ignore_index=True)
    # 3 群平均（seed をまたいだ平均。R-12 と同じ定義）
    g = ann.groupby("gene")[["n_exp", "n_cells", "ll_med", "lp_med"]].mean()
    g.columns = ["n_exp", "n_cells", "ll_med", "lp_med"]
    return g.assign(log2_n_exp=lambda d: np.log2(d.n_exp + 1.0))


def matching_feasibility(pool_ann: pd.DataFrame, pc_ann: pd.DataFrame, tol_near: float = 1.0) -> dict:
    """PC の遺伝子ごとに、マッチ距離 tol_near 以内のプール遺伝子が何個あるかを数える。"""
    res = {}
    for gene, t in pc_ann.iterrows():
        if not np.isfinite(t.lp_med) or not np.isfinite(t.ll_med) or not np.isfinite(t.n_exp):
            res[gene] = {"近傍数": 0, "理由": "注釈が欠けている"}
            continue
        cand = pool_ann.dropna()
        dist = (np.abs(cand.lp_med - t.lp_med) / 0.25 + np.abs(cand.ll_med - t.ll_med) / 0.25
                + np.abs(np.log(cand.n_exp / max(t.n_exp, 0.5))) / np.log(1.5))
        res[gene] = {"近傍数（距離 < 1.0）": int((dist < tol_near).sum()), "最小の 5 個の距離": [round(float(x), 3) for x in np.sort(dist.to_numpy())[:5]]}
    return res


def main() -> None:
    cfg = yaml.safe_load((EXP / "config.yaml").read_text())
    G = gene_sets(cfg)
    sets = G["sets"]
    sym = G["symbol"]
    out: dict = {
        "実験": "EXP-20261011-002", "仮説": "HYP-20261011-01",
        "目的": "結果を見る前の事前確認。主判定の値（z・セット平均・偽セットの閾値）は計算しない",
    }

    # ---- 1. 集合と重なり（総当たり）
    names = list(sets)
    ov = {}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            inter = sorted(set(sets[a]) & set(sets[b]))
            if inter:
                ov[f"{a} ∩ {b}"] = {"個数": len(inter), "記号": [sym.get(e, e) for e in inter]}
    out["1. 集合の個数"] = {k: len(v) for k, v in sets.items()}
    out["2. 重なり（結果を見る前に数えた。総当たり）"] = ov
    # 盲検
    blk = set(cfg["blocked_genes"])
    out["3. 盲検との重なり"] = {k: sorted(set(v) & blk) for k, v in sets.items()}

    # ---- 4. 注釈テーブル（A1）
    ann = load_a1_annotations(cfg).dropna(subset=["ll_med"])
    out["4. 注釈テーブル（A1。3 群平均）"] = {"遺伝子数": int(len(ann)),
                                          "n_exp の範囲": [int(ann.n_exp.min()), int(ann.n_exp.max())],
                                          "ll_med の範囲": [round(float(ann.ll_med.min()), 2), round(float(ann.ll_med.max()), 2)]}

    # ---- 5. 各 PC 版の検出率とマッチ近傍
    pool = [g for g in sets["プールA1"] if g in ann.index]
    pool_ann = ann.reindex(pool)
    pc_vers = {"退化PC23": sets["退化PC23"], "退化PC30": sets["退化PC30"], "退化PC21": sets["退化PC21"],
               "活性化PC71": sets["活性化PC71"]}
    feas = {}
    for name, genes in pc_vers.items():
        pc_ann = ann.reindex(genes)
        det = pc_ann.n_exp / pc_ann.n_cells
        feas[name] = {
            "検出率の中央値": round(float(det.median()), 4),
            "検出率 0.06 未満": [sym.get(g, g) for g in genes if det.get(g, 1.0) < 0.06],
            "注釈が欠けている": [sym.get(g, g) for g in genes if g not in pc_ann.index or not np.isfinite(pc_ann.loc[g, "lp_med"])],
            "プール重なり（引き出しから除く）": [sym.get(g, g) for g in genes if g in set(pool)],
        }
        if name in ("退化PC23", "活性化PC71"):
            feas[name]["マッチ近傍（PC 遺伝子ごと）"] = matching_feasibility(pool_ann, pc_ann)
    out["5. PC 版の検出率とマッチ近傍"] = feas

    # ---- 6. 材料の sha256
    sha = {}
    for name, p in (("A2 genes_011.json", A2_SRC_ROOT() / "work/genes_011.json"),
                    ("A1 config.yaml", A1_SRC / "config.yaml"),
                    ("A1 新規 EXP config.yaml", A1_NEW / "config.yaml"),
                    ("再現確認の参照 gene_rank_014.csv", A1_NEW / "results/gene_rank_014.csv")):
        sha[name] = {"sha256": sha256f(p), "実在": p.exists()}
    pre = json.loads((A1_NEW / "work/precheck_014.json").read_text())
    sha["トークン化済み入力の複製（A1）"] = pre["1. トークン化済み入力の複製"]["A1"]
    sha["トークン化済み入力の複製（A2）"] = pre["1. トークン化済み入力の複製"]["A2"]
    out["6. 材料の sha256"] = sha

    ok = (not any(out["3. 盲検との重なり"].values())
          and all(feas[k]["検出率の中央値"] >= 0.06 for k in ("退化PC23", "退化PC30", "退化PC21")))
    out["pass"] = bool(ok)
    out["注"] = ("結果を見る前の確認。z 値・セット平均・偽セットの閾値は計算していない。"
                 "不合格条件（盲検の重なり・検出率 0.06 未満の中央値）がなければ着手してよい。")
    (EXP / "work").mkdir(parents=True, exist_ok=True)
    (EXP / "work" / "precheck_015.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    print("written:", EXP / "work" / "precheck_015.json")


def A2_SRC_ROOT() -> Path:
    return ROOT / "experiments/EXP-20261007-001"


if __name__ == "__main__":
    main()
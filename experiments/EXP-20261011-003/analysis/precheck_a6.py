#!/usr/bin/env python3
"""EXP-20261011-003 の着手前チェック（埋め込みを 1 回も取らずに計算できる分）。

**この実験は埋め込みのみで摂動をしない。** したがってここで計算するのは
(1) 材料の固定、(2) 状態・提供者の数え上げ、(3) 盲検・遺伝子集合の重なり（循環）、
(4) 候補遺伝子の検出率、(5) 陽性対照の検出率とまとまり（K）である。
**状態の分離（AUC）と 2 本の軸のなす角は、埋め込みが無いと計算できない**（本スクリプトの対象外。
`work/embed_a6/` の後、`axis_feasibility.py` と `evaluate_a6.py` で扱う）。

読み取りのみ。書き込むのは `--out` の 1 ファイルだけ。
usage: .venv/bin/python analysis/precheck_a6.py --config config.yaml --out work/precheck_a6.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[3]


def rel(p: str) -> Path:
    q = Path(str(p))
    return q if q.is_absolute() else ROOT / q


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def geneformer_vocab() -> set[str]:
    # scripts/geneformer が本物の geneformer を隠すので外す
    sys.path = [p for p in sys.path if not p.endswith("/scripts")]
    from geneformer import TOKEN_DICTIONARY_FILE
    with open(TOKEN_DICTIONARY_FILE, "rb") as f:
        tok = pickle.load(f)
    return {k for k in tok if str(k).startswith("ENS")}


def zcols(C: np.ndarray) -> np.ndarray:
    L = np.log1p(C / np.maximum(C.sum(axis=1, keepdims=True), 1) * 1e4)
    sd = L.std(axis=0)
    return (L - L.mean(axis=0)) / np.where(sd == 0, 1.0, sd)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    import anndata as ad

    h5 = rel(cfg["data"]["input_h5ad"])
    out: dict = {"experiment_id": cfg["experiment_id"], "h5ad": str(h5.relative_to(ROOT)),
                 "h5ad_sha256": sha256(h5), "reads": "raw.X（整数カウント）、X は log 正規化"}

    A = ad.read_h5ad(str(h5))
    raw = A.raw
    X = raw.X.tocsr() if raw is not None else A.X.tocsr()
    var = raw.var if raw is not None else A.var
    ens = np.array(var.index.astype(str))
    sym = np.array(var["feature_name"].astype(str)) if "feature_name" in var.columns else ens
    e2i = {e: i for i, e in enumerate(ens)}
    vocab = geneformer_vocab()

    obs = A.obs
    ct = obs["celltype_l2"].astype(str).to_numpy()
    donor = obs["donor_id"].astype(str).to_numpy()
    disease = obs["disease"].astype(str).to_numpy()
    n_genes = obs["n_genes"].to_numpy(dtype=float)
    n_umi = obs["n_umi"].to_numpy(dtype=float)

    out["dataset"] = {"shape": list(A.shape), "n_vars": int(A.n_vars),
                      "assay": obs["assay"].astype(str).unique().tolist(),
                      "disease_counts": {k: int(v) for k, v in
                                         obs["disease"].astype(str).value_counts().items()}}
    out["geneformer_vocab"] = {"n_vars_in_vocab": int(sum(1 for e in ens if e in vocab)),
                               "n_vars": int(len(ens))}

    # ---- 状態の数え上げ -------------------------------------------------
    systems = cfg["systems"]
    state_tab = {}
    masks = {}
    for name, spec in systems.items():
        for role in ("start", "goal"):
            labels = spec[role]["labels"]
            m = np.isin(ct, labels)
            masks[(name, role)] = m
            state_tab[f"{name}.{role}"] = {
                "labels": labels, "n_cells": int(m.sum()),
                "n_donors": int(len(set(donor[m].tolist()))),
                "disease_counts": {k: int(v) for k, v in
                                   zip(*np.unique(disease[m], return_counts=True))},
                "n_genes_median": float(np.median(n_genes[m])),
                "n_umi_median": float(np.median(n_umi[m])),
                "donor_counts": {k: int(v) for k, v in
                                 zip(*np.unique(donor[m], return_counts=True))},
            }
    out["states"] = state_tab

    # ---- 深さをそろえたときに残る細胞数 --------------------------------
    def matched_counts(ma, mb, nbins=10):
        pool = np.concatenate([n_genes[ma], n_genes[mb]])
        edges = np.quantile(pool, np.linspace(0, 1, nbins + 1))
        edges[0] -= 1
        edges[-1] += 1
        ba, bb = np.digitize(n_genes[ma], edges), np.digitize(n_genes[mb], edges)
        return int(sum(min((ba == k).sum(), (bb == k).sum()) for k in range(1, nbins + 1)))

    depth = {}
    for name, spec in systems.items():
        ma, mb = masks[(name, "start")], masks[(name, "goal")]
        depth[name] = {"n_matched_per_state": matched_counts(ma, mb),
                       "start_n_genes_q": [float(x) for x in np.percentile(n_genes[ma], [5, 50, 95])],
                       "goal_n_genes_q": [float(x) for x in np.percentile(n_genes[mb], [5, 50, 95])]}
    out["depth_matchable"] = depth

    # ---- 遺伝子集合の重なり（循環） ------------------------------------
    gs = cfg["gene_sets"]
    pcsym = set(gs["positive_activation"]["symbols"])
    dsym = set(gs["positive_degenerate"]["symbols"])
    tubsym = set(gs["positive_tubular_a6"]["symbols"])
    candsym = set(gs["candidates"]["symbols"])
    bl = gs["blind5"]
    allsets = {"candidates": candsym, "positive_activation": pcsym,
               "positive_degenerate": dsym, "positive_tubular_a6": tubsym,
               "housekeeping": set(gs["housekeeping"]),
               "neutral": set(gs["neutral"])}
    overlaps = {}
    keys = list(allsets)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            ov = allsets[keys[i]] & allsets[keys[j]]
            overlaps[f"{keys[i]} ∩ {keys[j]}"] = sorted(ov)
    # 盲検（symbol と Ensembl の両方で照合）
    blind_ens = set(bl.values())
    blind_rows = []
    for nm, s in allsets.items():
        ens_set = set()
        for g in s:
            e = gs["candidates"]["ensembl"].get(g)
            if e:
                ens_set.add(e)
        if (s & set(bl)) or (ens_set & blind_ens):
            blind_rows.append(nm)
    out["overlaps"] = {"pairs": overlaps, "blind5_symbols": bl,
                       "blind5_collision_sets": blind_rows}

    # ---- 候補の検出率 ---------------------------------------------------
    floor = float(cfg["gene_sets"]["candidates"]["detection_floor_goal"])
    cand_det = {}
    cand_ens = gs["candidates"]["ensembl"]
    for g in gs["candidates"]["symbols"]:
        e = cand_ens.get(g)
        row = {"ensembl": e, "in_var": bool(e in e2i), "in_vocab": bool(e in vocab)}
        for (name, role), m in masks.items():
            if e not in e2i:
                continue
            col = np.asarray(X[m][:, e2i[e]].todense()).ravel()
            row[f"{name}.{role}"] = {"detect_frac": round(float((col > 0).mean()), 5),
                                     "cells": int((col > 0).sum()), "n": int(m.sum())}
        # 判定は「目標の状態」で（活性化の候補なので間質系の goal）
        d = row.get("interstitial.goal", {}).get("detect_frac", 0.0)
        row["a6_detected_by_floor"] = bool(d >= floor)
        cand_det[g] = row
    out["candidate_detection"] = {"floor": floor, "per_gene": cand_det,
                                  "n_detected": int(sum(v["a6_detected_by_floor"] for v in cand_det.values()))}

    # ---- 陽性対照の検出率・まとまり -------------------------------------
    def positive_metrics(symbols, ensmap, m):
        idx = [e2i[ensmap[s]] for s in symbols if ensmap.get(s) in e2i]
        missing = [s for s in symbols if ensmap.get(s) not in e2i]
        C = np.asarray(X[m][:, idx].todense(), dtype=np.float64)
        det = (C > 0).mean(axis=0)
        L = np.log1p(C / np.maximum(C.sum(axis=1, keepdims=True), 1) * 1e4)
        keep = det > 0
        if keep.sum() >= 2:
            R = np.corrcoef(L[:, keep].T)
            off = R[np.triu_indices(keep.sum(), 1)]
            corr_med = float(np.nanmedian(off))
        else:
            corr_med = None
        return {"n_found": len(idx), "n_missing": len(missing), "missing": missing,
                "detection_median": float(np.median(det)),
                "detection_min": float(det.min()), "detection_max": float(det.max()),
                "n_detected_ge_0.06": int((det >= 0.06).sum()),
                "corr_median": corr_med}

    out["positive_controls"] = {
        "activation_71_in_interstitial": positive_metrics(
            gs["positive_activation"]["symbols"], gs["positive_activation"]["ensembl"],
            masks[("interstitial", "start")] | masks[("interstitial", "goal")]),
        "degenerate_30_in_tubular": positive_metrics(
            gs["positive_degenerate"]["symbols"], gs["positive_degenerate"]["ensembl"],
            masks[("tubular", "start")] | masks[("tubular", "goal")]),
        # A6 用に組み直した管状系の陽性対照（PI 回答 2026-10-11）。目標の状態での検出率が床 0.10 を満たすかを見る
        "positive_tubular_a6_in_tubular": positive_metrics(
            gs["positive_tubular_a6"]["symbols"], gs["positive_tubular_a6"]["ensembl"],
            masks[("tubular", "goal")]),
        "positive_tubular_a6_floor": float(gs["positive_tubular_a6"]["detection_floor_goal"]),
        "K_activation": len(pcsym), "K_degenerate": len(dsym),
        "K_positive_tubular_a6": len(tubsym),
    }
    # 陽性対照の K は 20 以上か（着手条件）
    out["positive_controls"]["K_activation_ge_20"] = len(pcsym) >= 20
    out["positive_controls"]["K_degenerate_ge_20"] = len(dsym) >= 20
    out["positive_controls"]["K_positive_tubular_a6_ge_20"] = len(tubsym) >= 20

    out["note"] = ("状態の分離（AUC）・2 本の軸のなす角・遺伝子の位置は埋め込みが要る。"
                   "本ファイルは埋め込み前の固定の材料。")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str) + "\n")
    print(f"written: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

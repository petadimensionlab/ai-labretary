#!/usr/bin/env python3
"""EXP-20260927-002: 線維芽細胞に推定ラベル fib_tls_state を付ける（tokenize の前）。

config.yaml の `label:` セクションだけを参照する。
- tokenize の前に **生カウントであること** を再確認する（guardrails 2026-09-25）:
  X が整数、行和 ≈ obs['n_counts']。満たさなければ停止する。
- goal  : organizer_genes のいずれかが min_counts 以上
- start : organizer_genes がすべて 0、かつ where の条件
- alt   : organizer_genes がすべて 0、かつ where の条件
どれにも当たらない細胞は drop_unlabeled: true なら除く。

usage: label_fib_states.py --config experiments/EXP-20260927-002/config.yaml
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[3]


def rel(p: str) -> Path:
    p = Path(p)
    return p if p.is_absolute() else ROOT / p


def check_raw_counts(adata) -> dict:
    """生カウントの確認。整数でない・行和が n_counts と合わないなら停止する。"""
    X = adata.X
    data = X.data if hasattr(X, "data") else np.asarray(X).ravel()
    if data.size and not np.allclose(data, np.round(data)):
        sys.exit("X が整数カウントではありません（正規化値の可能性）。tokenize しないでください")
    rs = np.asarray(X.sum(axis=1)).ravel()
    out = {"counts_integer_like": True, "n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
           "row_sum_median": float(np.median(rs)), "min_row_sum": float(rs.min())}
    if "n_counts" in adata.obs.columns:
        nc = adata.obs["n_counts"].to_numpy(dtype=float)
        agree = float(np.mean(np.isclose(rs, nc, rtol=1e-5, atol=1e-3)))
        out["row_sum_vs_n_counts_frac_equal"] = round(agree, 6)
        if agree < 0.999:
            sys.exit(f"行和と n_counts が一致しません（一致率 {agree:.4f}）。入力を確認してください")
    else:
        sys.exit("obs に n_counts がありません（prepare_h5ad.py を先に実行してください）")
    return out


def matches(obs, where: dict | None) -> np.ndarray:
    m = np.ones(len(obs), dtype=bool)
    for col, vals in (where or {}).items():
        if col not in obs.columns:
            sys.exit(f"obs に列 {col} がありません")
        m &= obs[col].astype(str).isin([str(v) for v in vals]).to_numpy()
    return m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))["label"]
    import anndata as ad

    src = rel(cfg["input_h5ad"])
    adata = ad.read_h5ad(src)
    summary = {"input": str(src.relative_to(ROOT)), "raw_count_check": check_raw_counts(adata)}

    genes = list(cfg["organizer_genes"])
    missing = [g for g in genes if g not in set(adata.var_names.astype(str))]
    if missing:
        sys.exit(f"organizer_genes が var にありません: {missing}")
    cols = [adata.var_names.get_loc(g) for g in genes]
    sub = adata[:, cols].X
    vals = np.asarray(sub.todense() if hasattr(sub, "todense") else sub)
    thr = float(cfg.get("min_counts", 1))
    any_ge = (vals >= thr).any(axis=1)
    all_zero = (vals == 0).all(axis=1)
    summary["organizer_genes"] = {g: int((vals[:, i] >= thr).sum()) for i, g in enumerate(genes)}
    summary["min_counts"] = thr

    obs = adata.obs
    rules = {"goal": cfg["goal"], "start": cfg["start"], "alt": cfg["alt"]}
    labels = np.full(adata.n_obs, "unlabeled", dtype=object)
    for role in ("goal", "start", "alt"):          # 上から順に、最初に一致したものを採る
        r = rules[role]
        base = any_ge if r["rule"] == "any_ge_min" else all_zero
        m = (labels == "unlabeled") & base & matches(obs, r.get("where"))
        labels[m] = r["label"]
    key = cfg["state_key"]
    adata.obs[key] = labels.astype(str)

    for old, new in (cfg.get("rename_obs") or {}).items():
        if old in adata.obs.columns:
            adata.obs[new] = adata.obs[old].astype(str)
    if cfg.get("add_n_genes"):
        X = adata.X
        adata.obs[cfg["add_n_genes"]] = np.asarray((X > 0).sum(axis=1)).ravel().astype(int)

    if cfg.get("drop_unlabeled", True):
        adata = adata[labels != "unlabeled"].copy()

    dkey = "donor_id"
    summary["per_state"] = adata.obs[key].value_counts().to_dict()
    summary["donors_per_state"] = adata.obs.groupby(key, observed=True)[dkey].nunique().to_dict()
    summary["median_n_genes_per_state"] = adata.obs.groupby(key, observed=True)[cfg.get("add_n_genes", "n_genes")] \
        .median().to_dict()
    summary["disease_per_state"] = {s: adata.obs.loc[adata.obs[key] == s, "disease"].value_counts().to_dict()
                                    for s in summary["per_state"]}
    out = rel(cfg["output_h5ad"])
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    summary["output"] = str(out.relative_to(ROOT))
    summary["shape"] = list(adata.shape)
    (out.parent / "label_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

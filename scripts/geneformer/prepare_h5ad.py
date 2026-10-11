#!/usr/bin/env python3
"""
Labretary: 生データ h5ad を Geneformer tokenizer 用に整形する。

config.yaml の `prepare:` セクションだけを参照する:
  prepare:
    sources:                      # data/raw 配下（読み取り専用。書き込まない）
      - {path: data/raw/..., obs: {state_label: ccRCC}}   # obs: 全細胞に付与する定数列（任意）
    obs_keep: [cell_type, disease, donor_id]
    filter: {cell_type: [macrophage]}   # obs 値による絞り込み（任意）
    state_key: disease                  # max_cells_per_state の層別キー
    max_cells_per_state: 500            # 状態ごとの上限（null で全細胞）
    seed: 0
    output: experiments/EXP-.../work/input.h5ad
    # --- 以下は任意（EXP-20260926-001 で追加。書かなければ従来どおり） ---
    donor_maps:                         # 外部 CSV の列を obs に付ける（例: donor → 推定 TLS 群）
      - {path: outputs/.../donor_A2.csv, key: donor_id, csv_key: donor_id, column: TLS_group, as: TLS_group}
    state_rules:                        # 上から順に最初に一致したラベルを state_key 列に付ける
      - {label: aPT_CKD_TLShigh, where: {author_cell_type: [aPT], disease: [chronic kidney disease], TLS_group: [high]}}
    drop_unlabeled: true                # どの規則にも一致しない細胞を除く
    max_cells_per_donor: null           # (state_key, donor_key) ごとの上限（null で制限なし）
    donor_key: donor_id

処理: raw.X（なければ X）の整数カウントを X に置き、var に ensembl_id、obs に n_counts を付ける。
カウントが整数でなければ停止する（.claude/memory/guardrails.md）。
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]


def rel(p: str) -> Path:
    p = Path(p)
    return p if p.is_absolute() else ROOT / p


def load_counts(path: Path, obs_keep: list[str], const_obs: dict, filt: dict):
    import anndata as ad
    import scipy.sparse as sp
    a = ad.read_h5ad(path, backed="r")
    mask = np.ones(a.n_obs, dtype=bool)
    for col, vals in (filt or {}).items():
        mask &= a.obs[col].astype(str).isin([str(v) for v in vals]).to_numpy()
    idx = np.where(mask)[0]
    if a.raw is not None:
        X, var = a.raw.X[idx], a.raw.var.copy()
    else:
        X, var = a.X[idx], a.var.copy()
    X = sp.csr_matrix(X, dtype=np.float32)
    if X.nnz and not np.allclose(X.data[:200000], np.round(X.data[:200000])):
        sys.exit(f"{path}: 整数カウントではありません（正規化値の可能性）。tokenize しないでください")
    obs = a.obs.iloc[idx][[c for c in obs_keep if c in a.obs.columns]].copy()
    for c in obs.columns:
        obs[c] = obs[c].astype(str)
    for k, v in (const_obs or {}).items():
        obs[k] = str(v)
    obs["source"] = path.stem
    a.file.close()
    out = ad.AnnData(X=X, obs=obs, var=var[[]].copy())
    out.var["ensembl_id"] = out.var_names.astype(str)
    if "feature_name" in var.columns:
        out.var["gene_name"] = var["feature_name"].astype(str).to_numpy()
    return out


def add_donor_maps(adata, maps: list[dict]) -> None:
    """donor 単位の CSV（例: TLS 群の推定結果）を obs に写す。対応がない細胞は 'NA'。"""
    import pandas as pd
    for m in maps or []:
        tab = pd.read_csv(rel(m["path"]), dtype=str)
        lut = tab.set_index(m.get("csv_key", m["key"]))[m["column"]]
        if lut.index.duplicated().any():
            sys.exit(f"{m['path']}: {m.get('csv_key', m['key'])} が重複しています")
        adata.obs[m.get("as", m["column"])] = adata.obs[m["key"]].astype(str).map(lut).fillna("NA").to_numpy()


def apply_state_rules(adata, key: str, rules: list[dict], drop_unlabeled: bool):
    """state_rules を上から順に評価し、最初に一致した label を obs[key] に付ける。"""
    labels = np.full(adata.n_obs, "unlabeled", dtype=object)
    for r in rules:
        m = labels == "unlabeled"
        for col, vals in r["where"].items():
            m &= adata.obs[col].astype(str).isin([str(v) for v in vals]).to_numpy()
        labels[m] = r["label"]
    adata.obs[key] = labels.astype(str)
    return adata[labels != "unlabeled"].copy() if drop_unlabeled else adata


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))["prepare"]
    import anndata as ad

    parts = [load_counts(rel(s["path"]), cfg.get("obs_keep", []), s.get("obs"), cfg.get("filter"))
             for s in cfg["sources"]]
    adata = ad.concat(parts, join="inner", merge="first", index_unique="-")
    key, cap = cfg.get("state_key"), cfg.get("max_cells_per_state")
    add_donor_maps(adata, cfg.get("donor_maps"))
    if cfg.get("state_rules"):
        adata = apply_state_rules(adata, key, cfg["state_rules"], cfg.get("drop_unlabeled", True))
    dcap, dkey = cfg.get("max_cells_per_donor"), cfg.get("donor_key", "donor_id")
    if key and dcap:
        rng = np.random.default_rng(cfg.get("seed", 0))
        keep = []
        for _, pos in adata.obs.groupby([key, dkey], observed=True).indices.items():
            keep.extend(rng.choice(pos, size=min(dcap, len(pos)), replace=False))
        adata = adata[np.sort(keep)].copy()
    if key and cap:
        rng = np.random.default_rng(cfg.get("seed", 0))
        keep = []
        for _, pos in adata.obs.groupby(key, observed=True).indices.items():
            keep.extend(rng.choice(pos, size=min(cap, len(pos)), replace=False))
        adata = adata[np.sort(keep)].copy()
    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    out = rel(cfg["output"])
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    summary = {"output": str(out.relative_to(ROOT)), "shape": list(adata.shape),
               "per_state": adata.obs[key].value_counts().to_dict() if key else None,
               "per_source": adata.obs["source"].value_counts().to_dict()}
    dkey = cfg.get("donor_key", "donor_id")
    if key and dkey in adata.obs.columns:
        summary["donors_per_state"] = adata.obs.groupby(key, observed=True)[dkey].nunique().to_dict()
    out.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

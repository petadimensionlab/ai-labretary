#!/usr/bin/env python3
"""EXP-20260929-001: 追加の状態重心を Mac で作る（**摂動より先に済ませる**）。

作るのは 3 本。
  FIB_ecmHigh          対照の状態（第 3 の軸）。細胞外基質プログラム 7 個の得点が、
                       出発点の中で検出遺伝子数の十分位内 上位 5% の細胞
  MYOF_fib             筋線維芽細胞（第 2 の軸）。著者注釈 MYOF / cycMYOF
  FIB_orgHigh_unmatched 感度分析: 深さをそろえない目標の重心（**判定に使わない**）

細胞は `work/state_rows-seed{S}.json` の `extra_rows` / `unmatched_goal_rows` で指定済みで、
このスクリプトは行番号を読むだけ（選び方はここでは決めない）。**3 seed とも Mac 1 台で作る**ので、
第 2・第 3 の軸の定義が機体に依存しない（EXP-20260928-002 §8.4 と同じ考え方）。

出力: work/extra_centroids-seed{S}.npz, work/extra_centroids-seed{S}.json
実行: .venv/bin/python experiments/EXP-20260929-001/analysis/extra_centroids_003.py --config <cfg> --seed 61
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts/geneformer"))


def rel(p: str) -> Path:
    q = Path(p)
    return q if q.is_absolute() else ROOT / q


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--seed", type=int, required=True)
    a = ap.parse_args()
    import os
    import run_isp as ri
    from datasets import load_from_disk

    cfg = ri.load_config(Path(a.config))
    os.environ["GENEFORMER_DEVICE"] = ri.resolve_device(cfg)
    ri.apply_dtype(cfg)
    ri.single_process_datasets(cfg)

    E = rel(cfg["output"]["experiment_dir"])
    ds = load_from_disk(str(rel(cfg["data"]["tokenized_dataset"])))
    key, dkey = cfg["data"]["state_key"], cfg["cell_selection"]["donor_key"]
    donors = [str(x) for x in ds[dkey]]
    length = np.array(ds["length"], dtype=float)
    spec = json.loads((E / f"work/state_rows-seed{a.seed}.json").read_text())
    start_rows = set(spec["rows"][cfg["states"]["start_state"]])
    goal_rows = set(spec["rows"][cfg["states"]["goal_state"]])

    all_rows: dict[str, list[int]] = {k: [int(i) for i in v] for k, v in spec["extra_rows"].items()}
    all_rows["FIB_orgHigh_unmatched"] = [int(i) for i in spec["unmatched_goal_rows"]]

    out_info = {}
    for name, keep in all_rows.items():
        out_info[name] = {"n_cells": len(keep), "n_donors": len({donors[i] for i in keep}),
                          "median_length": float(np.median(length[keep])),
                          "overlap_with_start_centroid": len(set(keep) & start_rows),
                          "overlap_with_goal_centroid": len(set(keep) & goal_rows)}
    for name in ("FIB_ecmHigh", "MYOF_fib"):
        if out_info[name]["overlap_with_start_centroid"] or out_info[name]["overlap_with_goal_centroid"]:
            sys.exit(f"{name}: 目標・出発点の重心と細胞が重なっている（設計では 0 のはず）")

    rows = sorted({i for v in all_rows.values() for i in v})
    work = E / "work"
    dpath = work / f"extra-centroid-cells-seed{a.seed}.dataset"
    if not dpath.exists():
        sub = ds.select(rows)
        sub = sub.map(lambda ex, i: {"row_id": rows[i]}, with_indices=True)
        sub.save_to_disk(str(dpath))
    from geneformer import EmbExtractor
    m, p = cfg["model"], cfg["perturbation"]
    ee = ri.call_compat(EmbExtractor, {
        "model_type": m["type"], "num_classes": m.get("num_classes", 0), "filter_data": None,
        "max_ncells": None, "emb_mode": p.get("emb_mode", "cls").replace("_and_gene", ""),
        "emb_layer": m.get("emb_layer", -1), "emb_label": [key, dkey, "row_id"],
        "summary_stat": None, "forward_batch_size": p.get("forward_batch_size", 8),
        "nproc": p.get("nproc", 1), "model_version": m.get("version")}, "EmbExtractor")
    outdir = work / f"extra-centroid-embs-seed{a.seed}"
    outdir.mkdir(parents=True, exist_ok=True)
    embs = ee.extract_embs(str(rel(m["directory"])), str(dpath), str(outdir), "cell_embs")
    if "row_id" not in embs.columns:
        raise SystemExit("cell_embs に row_id 列がありません（ID で結合できません）")
    X = embs.drop(columns=[key, dkey, "row_id"]).to_numpy(dtype=np.float64)
    by_row = dict(zip(embs["row_id"].astype(int).to_numpy(), range(len(embs))))
    out_vec = {}
    for name, keep in all_rows.items():
        idx = [by_row[i] for i in keep]
        assert len(idx) == len(keep), f"{name}: row_id で結合できない細胞があります"
        out_vec[name] = X[idx].mean(axis=0).astype(np.float32)
    np.savez(str(rel(cfg["extra_centroids"]["output"].format(seed=a.seed))), **out_vec)

    names = list(out_vec)
    cos = {f"{i}×{j}": round(float(out_vec[i] @ out_vec[j]
                                   / (np.linalg.norm(out_vec[i]) * np.linalg.norm(out_vec[j]))), 5)
           for n, i in enumerate(names) for j in names[n + 1:]}
    (work / f"extra_centroids-seed{a.seed}.json").write_text(json.dumps(
        {"seed": a.seed, "sets": out_info, "dim": int(X.shape[1]),
         "centroid_cosines": cos}, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps({"sets": out_info, "cos": cos}, ensure_ascii=False))


if __name__ == "__main__":
    main()

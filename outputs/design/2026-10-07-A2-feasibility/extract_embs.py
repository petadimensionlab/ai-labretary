#!/usr/bin/env python3
"""A2 線維芽細胞系列の細胞埋め込みを float32 で抽出する（摂動は行わない）。

出力: work/cell_embs.npz（X: 細胞 × 次元、row_id: トークン化済みデータの行番号）
"""
import sys, time, json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts/geneformer"))
import os
import run_isp as ri
from datasets import load_from_disk
from geneformer import EmbExtractor

CFG = ROOT / "outputs/design/2026-10-07-A2-feasibility/run_config_extract.yaml"
WORK = ROOT / "outputs/design/2026-10-07-A2-feasibility/work"

def main(limit=None, out="cell_embs.npz"):
    cfg = ri.load_config(CFG)
    os.environ["GENEFORMER_DEVICE"] = ri.resolve_device(cfg)
    dtype = ri.apply_dtype(cfg)
    ri.single_process_datasets(cfg)
    ds = load_from_disk(str(WORK / "tokenized.dataset"))
    import json as _j
    rf = WORK / "emb_rows.json"
    if limit is None and rf.exists():
        rows = [int(i) for i in _j.loads(rf.read_text())]
    else:
        rows = list(range(len(ds))) if limit is None else list(range(limit))
    dpath = WORK / (f"embcells-{len(rows)}.dataset")
    if not dpath.exists():
        ds.select(rows).map(lambda ex, i: {"row_id": rows[i]}, with_indices=True).save_to_disk(str(dpath))
    m, p = cfg["model"], cfg["perturbation"]
    ee = ri.call_compat(EmbExtractor, {
        "model_type": m["type"], "num_classes": m.get("num_classes", 0), "filter_data": None,
        "max_ncells": None, "emb_mode": p.get("emb_mode", "cls").replace("_and_gene", ""),
        "emb_layer": m.get("emb_layer", -1), "emb_label": ["row_id"], "summary_stat": None,
        "forward_batch_size": p.get("forward_batch_size", 8), "nproc": 1,
        "model_version": m.get("version")}, "EmbExtractor")
    odir = WORK / f"embs-{len(rows)}"
    odir.mkdir(parents=True, exist_ok=True)
    t = time.time()
    embs = ee.extract_embs(str(ROOT / m["directory"]), str(dpath), str(odir), "cell_embs")
    el = time.time() - t
    rid = embs["row_id"].astype(int).to_numpy()
    X = embs.drop(columns=["row_id"]).to_numpy(dtype=np.float64)
    np.savez(WORK / out, X=X, row_id=rid)
    print(json.dumps({"n_cells": len(rid), "dim": X.shape[1], "seconds": round(el, 1),
                      "device": os.environ["GENEFORMER_DEVICE"], "dtype": dtype,
                      "out": str(WORK / out)}, ensure_ascii=False))

if __name__ == "__main__":
    lim = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1] != "all" else None
    main(lim, sys.argv[2] if len(sys.argv) > 2 else "cell_embs.npz")

#!/usr/bin/env python3
"""EXP-20260929-001: work/input.h5ad を Geneformer V2 の入力にトークン化する。

摂動は行わない。run_isp.py の tokenize() をそのまま呼ぶ（同じ実装を使うため）。
出力: work/tokenized.dataset
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20260929-001"
sys.path.insert(0, str(ROOT / "scripts/geneformer"))

import run_isp as ri  # noqa: E402

CFG = {
    "_path": str(E / "config.yaml"),
    "model": {"version": "V2"},
    "data": {
        "input_h5ad": "experiments/EXP-20260929-001/work/input.h5ad",
        "tokenized_dataset": "experiments/EXP-20260929-001/work/tokenized.dataset",
        "custom_attr_name_dict": {
            "fib_org_state": "fib_org_state", "ecm_high": "ecm_high", "donor_id": "donor_id",
            "subclass_l2": "subclass_l2", "disease": "disease", "region": "region",
            "n_genes": "n_genes", "n_counts": "n_counts", "a1_row": "a1_row"},
    },
    "perturbation": {"nproc": 8},
}

if __name__ == "__main__":
    out = ri.tokenize(CFG, E / "work")
    from datasets import load_from_disk
    ds = load_from_disk(str(out))
    print(ds)
    import collections
    print(collections.Counter(ds["fib_org_state"]))
    print("ecm_high:", collections.Counter(ds["ecm_high"]))

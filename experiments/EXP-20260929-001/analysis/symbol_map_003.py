#!/usr/bin/env python3
"""EXP-20260929-001: 遺伝子記号 → Ensembl ID の対応表を work/input.h5ad の var から作る。
同じ記号が複数行にある場合は最初の 1 つを採る（extract_a1_fib.py と同じ規則）。
Geneformer の辞書（token_dictionary）に入っているかも併せて記録する。
出力: work/symbol_to_ensembl_003.json, work/symbol_map_003_meta.json
"""
import json, pickle, sys
from pathlib import Path
import h5py, numpy as np
ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20260929-001"
sys.path.insert(0, str(ROOT / "scripts/geneformer"))
import run_isp as ri

import anndata as ad
a = ad.read_h5ad(E / "work/input.h5ad", backed="r")
ens = a.var["ensembl_id"].astype(str).to_numpy()
sym = a.var["gene_name"].astype(str).to_numpy()
a.file.close()
m = {}
for s, e in zip(sym, ens):
    m.setdefault(str(s), str(e))
tokd = ri.gene_token_dict()
(E / "work/symbol_to_ensembl_003.json").write_text(json.dumps(m, ensure_ascii=False) + "\n")
groups = json.loads((E / "analysis/excluded_symbols_003.json").read_text())["groups"]
meta = {"n_symbols": len(m), "token_dict_size": len(tokd),
        "groups": {k: {s: {"ensembl": m.get(s), "in_dict": bool(m.get(s) in tokd)} for s in v}
                   for k, v in groups.items()}}
(E / "work/symbol_map_003_meta.json").write_text(json.dumps(meta, indent=1, ensure_ascii=False) + "\n")
for k, v in meta["groups"].items():
    bad = [s for s, d in v.items() if not d["in_dict"]]
    print(f"{k}: {len(v)} 個、辞書に無い {len(bad)} {bad}")

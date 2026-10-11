"""EXP-20260927-003: 追加の状態重心（MYOF 重心・深さマッチなしの goal 重心）を作る。

摂動（ISP）は行わない。config.yaml の extra_centroids に書いた条件で細胞を選び、
同じモデル・同じ埋め込み種別（cls）・同じ精度で細胞埋め込みを取り、その平均（重心）を保存する。
細胞の選び方は select_genes_fib.py の深さマッチと同じ規則（log トークン長の 1:1 最近傍, 許容 0.10）。

出力: work/extra_centroids-seed{S}.npz（キー = 重心の名前。値 = 1 次元ベクトル）
      work/extra_centroids-seed{S}.json（使った細胞数・提供者数・トークン長の中央値・start との重なり）
実行: .venv/bin/python experiments/EXP-20260927-003/analysis/extra_centroids.py --config <config.yaml> --seed 21
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np
import yaml

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
    import run_isp as ri
    from datasets import load_from_disk

    cfg = ri.load_config(Path(a.config))
    ri.resolve_device(cfg)
    import os
    os.environ["GENEFORMER_DEVICE"] = ri.resolve_device(cfg)
    ri.apply_dtype(cfg)
    ri.single_process_datasets(cfg)

    sel = __import__("importlib.util", fromlist=["util"])
    spec = sel.spec_from_file_location("sgf", Path(__file__).with_name("select_genes_fib.py"))
    sgf = sel.module_from_spec(spec)
    spec.loader.exec_module(sgf)

    E = rel(cfg["output"]["experiment_dir"])
    tok = rel(cfg["data"]["tokenized_dataset"])
    ds = load_from_disk(str(tok))
    key, dkey = cfg["data"]["state_key"], cfg["cell_selection"]["donor_key"]
    states = np.array([str(x) for x in ds[key]])
    donors = [str(x) for x in ds[dkey]]
    subcl = np.array([str(x) for x in ds["subclass_l2"]])
    length = np.array(ds["length"], dtype=float)
    log_len = np.log(length)
    start_state = cfg["states"]["start_state"]

    start_rows = set(json.loads((E / f"work/state_rows-seed{a.seed}.json").read_text())["rows"][start_state])

    out_vec, out_info, all_rows = {}, {}, {}
    for name, s in cfg["extra_centroids"]["sets"].items():
        m = np.isin(states, s["states"])
        if s.get("subclass_l2"):
            m &= np.isin(subcl, s["subclass_l2"])
        rng = random.Random(f"{a.seed}-{name}")
        keep = sgf.balanced(np.flatnonzero(m).tolist(), donors, s["max_cells"], s["max_cells_per_donor"], rng)
        if s.get("depth_match_to"):
            ref = np.zeros(len(states), bool)
            ref[sorted(json.loads((E / f"work/state_rows-seed{a.seed}.json").read_text())
                       ["rows"][s["depth_match_to"]])] = True
            pos = np.zeros(len(states), bool)
            pos[keep] = True
            mm = sgf.depth_match(log_len, pos, ref, sgf.DEPTH_TOL)
            keep = sorted(np.flatnonzero(mm & pos).tolist())
        all_rows[name] = keep
        out_info[name] = {"n_cells": len(keep), "n_donors": len({donors[i] for i in keep}),
                          "median_length": float(np.median(length[keep])),
                          "overlap_with_start_centroid": len(set(keep) & start_rows),
                          "depth_match_to": s.get("depth_match_to")}

    # 埋め込みは 1 回の EmbExtractor 呼び出しでまとめて取り、名前ごとに平均する
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
        "model_type": m["type"], "num_classes": m.get("num_classes", 0), "filter_data": None, "max_ncells": None,
        "emb_mode": p.get("emb_mode", "cls").replace("_and_gene", ""), "emb_layer": cfg["model"].get("emb_layer", -1),
        "emb_label": [key, dkey, "row_id"], "summary_stat": None,
        "forward_batch_size": p.get("forward_batch_size", 8), "nproc": p.get("nproc", 1),
        "model_version": m.get("version")}, "EmbExtractor")
    outdir = work / f"extra-centroid-embs-seed{a.seed}"
    outdir.mkdir(parents=True, exist_ok=True)
    embs = ee.extract_embs(str(rel(m["directory"])), str(dpath), str(outdir), "cell_embs")
    if "row_id" not in embs.columns:
        raise SystemExit("cell_embs に row_id 列がありません（ID で結合できません）")
    X = embs.drop(columns=[key, dkey, "row_id"]).to_numpy(dtype=np.float64)
    by_row = dict(zip(embs["row_id"].astype(int).to_numpy(), range(len(embs))))
    for name, keep in all_rows.items():
        idx = [by_row[i] for i in keep]
        assert len(idx) == len(keep), f"{name}: row_id で結合できない細胞があります"
        out_vec[name] = X[idx].mean(axis=0).astype(np.float32)
    np.savez(str(rel(cfg["extra_centroids"]["output"].format(seed=a.seed))), **out_vec)
    (work / f"extra_centroids-seed{a.seed}.json").write_text(
        json.dumps({"seed": a.seed, "sets": out_info, "dim": int(X.shape[1])}, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps(out_info, ensure_ascii=False))


if __name__ == "__main__":
    main()

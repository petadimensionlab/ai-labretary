#!/usr/bin/env python3
"""EXP-20261011-003 — A6 全 46,957 核の埋め込み（plan.md §7 手順 3・4・7）。

摂動はしない。Geneformer（V2-316M・fp32・tf32 閉・cls・emb_label=["row_id"]・max_ncells=None）で
A6 全核を 1 回エンコードし、`results/embeddings.parquet` に row_id で保存する。

- 読み取りのみ: data/raw の A6 h5ad（sha256 固定値）、config.yaml、plan.md
- 書き込み: work/tokenized.dataset、work/tokenized-rowid.dataset、results/（または results-smoke/）
- --dry-run: 材料の固定値照合と、G2 の読み（needs_G2_approval）を表示するだけ（GPU も埋め込みもしない）
- --smoke: 埋め込みを Geneformer 内部のダウンサンプル（max_ncells=500）で試し、results-smoke/ に書く。
  判定ではない。evaluate_a6.py --smoke で最後まで通すために使う（ガードレール 2026-10-01）
- 本番は 1 回だけ。results/embeddings.parquet が既にあると止める（--force で上書き）

使い方:
    .venv/bin/python experiments/EXP-20261011-003/analysis/embed_a6.py --config config.yaml --dry-run
    LABRETARY_DEVICE=mps .venv/bin/python experiments/EXP-20261011-003/analysis/embed_a6.py --config config.yaml --smoke
    LABRETARY_DEVICE=cuda .venv/bin/python experiments/EXP-20261011-003/analysis/embed_a6.py --config config.yaml   # 本番（リモート GPU）
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
EXP = ROOT / "experiments/EXP-20261011-003"
sys.path.insert(0, str(ROOT / "scripts/geneformer"))
import run_isp as ri  # noqa: E402   call_compat・apply_dtype・resolve_device・tokenize をそのまま使う

MATERIAL_SHA256 = "8d9cca1320a7952f4b2eaff07acc03494e840e46b6b333e1c3efaf395165b29d"  # plan.md §2.1
STATE_COL = "celltype_l2"


def sha256f(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def rel(p) -> Path:
    q = Path(str(p))
    return q if q.is_absolute() else ROOT / q


def load_confirmed_hashes() -> dict:
    """SHA256SUMS-preregistered.txt の（ファイル, ハッシュ）表。固定値と突き合わせる（A-2）。"""
    out = {}
    f = EXP / "analysis" / "SHA256SUMS-preregistered.txt"
    if f.exists():
        for line in f.read_text().splitlines():
            h, name = line.split("  ", 1)
            out[name.strip()] = h
    return out


def depth_match_ids(state_labels: np.ndarray, n_genes: np.ndarray,
                    start_labels: list[str], goal_labels: list[str], nbins: int = 10):
    """n_genes の 10 分位マッチ（plan.md §6.1）。決定的: 各ビンで (n_genes, row_id) の昇順から
    min(na, nb) 個を採用する。返り値は dataset 上の行番号（row_id）。"""
    s = np.isin(state_labels, start_labels)
    g = np.isin(state_labels, goal_labels)
    ids = np.arange(len(state_labels))
    pool = np.concatenate([n_genes[s], n_genes[g]])
    edges = np.quantile(pool, np.linspace(0, 1, nbins + 1))
    edges[0] -= 1
    edges[-1] += 1
    bs = np.digitize(n_genes[s], edges)
    bg = np.digitize(n_genes[g], edges)
    s_ids, g_ids = ids[s], ids[g]
    ms, mg = [], []
    for k in range(1, nbins + 1):
        a = sorted(s_ids[bs == k], key=lambda i: (float(n_genes[i]), int(i)))
        b = sorted(g_ids[bg == k], key=lambda i: (float(n_genes[i]), int(i)))
        m = min(len(a), len(b))
        ms += [int(x) for x in a[:m]]
        mg += [int(x) for x in b[:m]]
    return ms, mg


def load_meta_from_dataset(ds) -> pd.DataFrame:
    """tokenized-rowid.dataset の行番号（row_id）ごとの注釈表。dtype をそろえて返す。"""
    cols = ["row_id", "celltype_l2", "celltype_l1", "donor_id", "disease", "assay", "n_genes", "n_umi", "length"]
    meta = pd.DataFrame({c: ds[c] for c in cols if c in ds.column_names})
    meta["row_id"] = meta["row_id"].astype(np.int64)
    for c in ("n_genes", "n_umi", "length"):
        if c in meta.columns:
            meta[c] = meta[c].astype(np.float64)
    return meta


def build_centroids(embs: np.ndarray, meta: pd.DataFrame, matched: dict) -> tuple[dict, dict]:
    """検出された遺伝子数の分布をそろえた（10 分位マッチ）セントロイドと、マッチの記録を作る。

    matched: {"interstitial": {"start": [...], "goal": [...]}, "tubular": {...}}（row_id の列）
    セントロイドの名前は axis-A6-*.yaml が参照する（plan.md §5.5）: Fibroblast / Myofibroblast / PT_normal / PT_inj_infl
    """
    by_id = {int(r): i for i, r in enumerate(meta["row_id"].to_numpy())}
    names = {"interstitial": ("Fibroblast", "Myofibroblast"),
             "tubular": ("PT_normal", "PT_inj_infl")}
    cents: dict[str, np.ndarray] = {}
    rec: dict = {}
    for sys_name, (start_nm, goal_nm) in names.items():
        sub = matched[sys_name]
        for nm, key in ((start_nm, "start"), (goal_nm, "goal")):
            idx = [by_id[int(x)] for x in sub[key] if int(x) in by_id]
            rec[f"{sys_name}.{key}"] = {"n_matched": len(sub[key]), "n_available": len(idx)}
            if len(idx):
                cents[nm] = embs[idx].mean(axis=0).astype(np.float64)
    return cents, rec


def write_env_json() -> None:
    """scripts/capture_env.sh をそのまま呼ぶ（config ハッシュ・dtype・seed を含む env.json を results/ に書く）。"""
    import subprocess
    subprocess.run(["bash", str(ROOT / "scripts" / "capture_env.sh"), str(EXP)],
                   check=False, capture_output=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--force", action="store_true", help="既存の embeddings.parquet を上書きする（本番は原則 1 回）")
    args = ap.parse_args()

    cfg_path = rel(args.config)
    cfg = yaml.safe_load(cfg_path.read_text())
    cfg.setdefault("perturbation", {"nproc": int(cfg.get("embedding", {}).get("nproc", 1))})
    # この config は device/dtype を compute: ではなく embedding: に持つ（plan.md §7・config.yaml）。
    # run_isp.resolve_device / apply_dtype は compute: を見るので、実行前に in-memory で写す。
    cfg.setdefault("compute", {}).setdefault("device", cfg.get("embedding", {}).get("device", "cuda"))
    cfg["compute"].setdefault("dtype", cfg.get("embedding", {}).get("dtype", "fp32"))
    fixed = load_confirmed_hashes()
    results = EXP / ("results-smoke" if args.smoke else "results")
    results.mkdir(parents=True, exist_ok=True)
    out_emb = results / "embed_a6"
    out_emb.mkdir(parents=True, exist_ok=True)

    # ---- A-2 材料の固定値 ------------------------------------------------
    h5_raw = rel(cfg["data"]["input_h5ad"])
    checks = {
        "材料 A6 h5ad の sha256 == 固定値": sha256f(h5_raw) == MATERIAL_SHA256,
    }
    want_cfg = fixed.get("experiments/EXP-20261011-003/config.yaml")
    if want_cfg:
        checks["config.yaml の sha256 == 事前登録（SHA256SUMS）"] = sha256f(cfg_path) == want_cfg
    want_plan = fixed.get("experiments/EXP-20261011-003/plan.md")
    if want_plan:
        checks["plan.md の sha256 == 事前登録（SHA256SUMS）"] = sha256f(EXP / "plan.md") == want_plan
    for k, ok in checks.items():
        print(("OK  " if ok else "NG  ") + k)
    if not all(checks.values()):
        raise SystemExit("**止める**(A-2): 材料または config/plan の sha256 が固定値と一致しない（B-1）")

    # 材料の sha256 確認は raw で行った（上）。tokenize は整形済みの work/input.h5ad を使う（plan.md §7）。
    # config.yaml の data.input_h5ad は raw を指す（読み取り専用のまま）ため、実行時のみ差し替える。
    cfg["data"]["input_h5ad"] = str(EXP / "work" / "input.h5ad")

    comp = cfg["compute"]
    print("G2（config の読み）:")
    print(f"  needs_G2_approval            = {comp.get('needs_G2_approval')}")
    print(f"  needs_G2_approval_conservative = {comp.get('needs_G2_approval_conservative')}")
    print(f"  approved_by                  = {comp.get('approved_by')}")
    print(f"  compute.device               = {comp.get('device')} / dtype = {comp.get('dtype')}")

    if args.dry_run:
        print("dry-run OK（埋め込みは実行しない）")
        return 0

    if not args.smoke and (results / "embeddings.parquet").exists() and not args.force:
        raise SystemExit(f"**止める**: 本番の埋め込みは既にあります {results / 'embeddings.parquet'}。"
                         f"1 回だけ回す計画です（--force で上書きしない限り再実行しない）")

    t_start = time.time()
    device = ri.resolve_device(cfg)
    os.environ["GENEFORMER_DEVICE"] = device
    dtype = ri.apply_dtype(cfg)  # fp32: torch の tf32 経路を明示的に閉じる（CLAUDE.md §6.1）
    ri.single_process_datasets(cfg)
    work = EXP / "work"
    work.mkdir(parents=True, exist_ok=True)

    # ---- tokenize ---------------------------------------------------------
    t0 = time.time()
    src = work / "h5ad_in"
    if src.exists():                     # 前回の失敗で残った symlink があると tokenizer が全ファイルを読む
        for p in src.iterdir():
            p.unlink()
    tok = ri.tokenize(cfg, work)   # work/tokenized.dataset（なければ作る）
    t_tok = time.time() - t0

    # ---- row_id 付き dataset（precheck と同じ流儀。追加の forward は無い）----
    tok_rowid = work / "tokenized-rowid.dataset"
    if not tok_rowid.exists():
        from datasets import load_from_disk
        ds0 = load_from_disk(str(tok))
        ds0 = ds0.map(lambda ex, idx: {"row_id": idx}, with_indices=True)
        ds0.save_to_disk(str(tok_rowid))
    from datasets import load_from_disk
    ds = load_from_disk(str(tok_rowid))

    # row_id == 入力 h5ad の行番号であることを確かめる（順序の照合。A-2）
    import anndata as ad
    A = ad.read_h5ad(str(rel(cfg["data"]["input_h5ad"])))
    seq_ds = [str(x) for x in ds["celltype_l2"]]
    seq_h5 = [str(x) for x in A.obs[STATE_COL].astype(str).to_numpy()]
    if seq_ds != seq_h5:
        raise SystemExit(f"**止める**: tokenized dataset と input.h5ad の細胞順が一致しない"
                         f"（row_id を対応付けられない。前処理を疑う。B-3）")
    del A

    # ---- 埋め込み（cls、fp32、tf32 閉） ----------------------------------
    from geneformer import EmbExtractor
    from geneformer import perturber_utils as pu
    max_ncells = 500 if args.smoke else None
    ee = ri.call_compat(EmbExtractor, {
        "model_type": cfg["model"]["type"], "num_classes": int(cfg["model"].get("num_classes", 0)),
        "filter_data": cfg["data"].get("filter_data"),
        "max_ncells": max_ncells,
        "emb_mode": cfg["embedding"]["emb_mode"],          # cls
        "emb_layer": int(cfg["model"].get("emb_layer", -1)),
        "emb_label": cfg["embedding"]["emb_label"],        # ["row_id"]
        "summary_stat": None,
        "forward_batch_size": int(cfg["embedding"].get("forward_batch_size", 100)),
        "nproc": int(cfg["embedding"].get("nproc", 1)),
        "model_version": cfg["model"]["version"],
    }, "EmbExtractor")
    t0 = time.time()
    embs_df = ee.extract_embs(str(rel(cfg["model"]["directory"])), str(tok_rowid), str(out_emb), "cell_embs")
    t_emb = time.time() - t0

    meta = load_meta_from_dataset(ds)
    n_dim = embs_df.shape[1] - 1   # row_id 以外
    emb_cols = [c for c in embs_df.columns if c != "row_id"]
    if len(emb_cols) != n_dim:
        raise SystemExit(f"**止める**: 埋め込み列が想定と違う（{len(emb_cols)} 列）")
    embs = embs_df[emb_cols].to_numpy(dtype=np.float64)
    df = pd.DataFrame({"row_id": embs_df["row_id"].astype(np.int64).to_numpy()})
    df = df.merge(meta, on="row_id", how="left")
    if df["row_id"].isna().any() or df["celltype_l2"].isna().any():
        raise SystemExit("**止める**: 埋め込みと注釈が row_id で結びつかない細胞がある")
    for j, c in enumerate(emb_cols):
        df[f"e{j:04d}"] = embs[:, j].astype(np.float32)
    df = df.sort_values("row_id").reset_index(drop=True)
    # 以後のセントロイドは row_id 昇順の df を使う（sort 前の embs は使わない。対応付けを誤らない）
    embs_sorted = df[[f"e{j:04d}" for j in range(n_dim)]].to_numpy(dtype=np.float64)
    df.to_parquet(results / "embeddings.parquet", index=False)

    # ---- セントロイド（検出された遺伝子数の分布をそろえた。axis_feasibility が参照） ---------------------
    st = np.array([str(x) for x in ds[STATE_COL]])
    ng = np.array([float(x) for x in ds["n_genes"]])
    sys_cfg = cfg["systems"]
    matched = {}
    for name, spec in sys_cfg.items():
        ms, mg = depth_match_ids(st, ng,
                                 [str(x) for x in spec["start"]["labels"]],
                                 [str(x) for x in spec["goal"]["labels"]])
        matched[name] = {"start": ms, "goal": mg}
    cents, cent_rec = build_centroids(embs_sorted, df, matched)
    np.savez(out_emb / "centroids.npz", **{k: v.astype(np.float32) for k, v in cents.items()})
    (out_emb / "depth_match.json").write_text(json.dumps(matched, indent=1) + "\n")
    (out_emb / "centroid_rec.json").write_text(json.dumps(cent_rec, indent=1) + "\n")

    # ---- attn_implementation の記録（モデルを 1 回読む。実測の数値条件） ----
    import torch
    import transformers
    m = pu.load_model(cfg["model"]["type"], int(cfg["model"].get("num_classes", 0)),
                      str(rel(cfg["model"]["directory"])), mode="eval")
    attn = getattr(getattr(m, "config", None), "_attn_implementation", "unknown")
    numerics = {
        "torch": torch.__version__, "transformers": transformers.__version__,
        "attn_implementation": attn,
        "cuda_matmul_allow_tf32": bool(torch.backends.cuda.matmul.allow_tf32),
        "cudnn_allow_tf32": bool(torch.backends.cudnn.allow_tf32),
        "float32_matmul_precision": torch.get_float32_matmul_precision(),
    }
    if torch.cuda.is_available():
        numerics["cuda_device"] = torch.cuda.get_device_name(0)
        numerics["cuda_version"] = torch.version.cuda
    meta_run = {
        "experiment_id": cfg["experiment_id"], "mode": "smoke" if args.smoke else "final",
        "device": device, "dtype": dtype, "config_sha256": sha256f(cfg_path),
        "seed": None, "note_seed": "摂動なし・単一エンコード。乱数は使わない（Geneformer のダウンサンプルは seed=42 固定）",
        "n_cells_embedded": int(len(df)), "n_dim": int(n_dim),
        "max_ncells": max_ncells,
        "embedding_seconds": round(t_emb, 1),
        "tokenize_seconds": round(t_tok, 1),
        "total_seconds": round(time.time() - t_start, 1),
        "host": os.uname().nodename,
        "numerics": numerics,
        "material_sha256": MATERIAL_SHA256,
        "n_donors": int(df["donor_id"].nunique()) if "donor_id" in df.columns else None,
    }
    (out_emb / "run_meta.json").write_text(json.dumps(meta_run, indent=1, ensure_ascii=False) + "\n")
    write_env_json()

    print(json.dumps({"rows": int(len(df)), "n_dim": n_dim, "device": device, "dtype": dtype,
                      "embedding_seconds": round(t_emb, 1),
                      "centroids": sorted(cents)}, ensure_ascii=False))
    print(f"書き込み: {results / 'embeddings.parquet'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit as e:
        if isinstance(e.code, int) or e.code is None:
            raise
        print(str(e.code), file=sys.stderr)
        raise
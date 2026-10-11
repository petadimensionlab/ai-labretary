"""EXP-20260926-003 phase 4: Pretrained（EXP-20260926-002 案 A）の対照 ISP を、FT 側と同じ goal 重心
（held-out 健常 9 donor だけ）で再計算する（plan.md §4.4）。
- goal 重心: config-pretrained-matched.yaml の state_embs（include_donors）で、Pretrained の EmbExtractor により seed 3/4/5 ごとに抽出（MPS）
- shift_goal: 案 A で保存した元/摂動後の CLS 埋め込み（genes/<ENSG>/emb/{orig,pert}.npy）から cos(pert, g) − cos(orig, g) を再計算（float64）
- 出力: results-pretrained-matched-recalc/run-*/（genes は案 A へのシンボリックリンク、state_emb_vectors.npz は goal だけ置き換え）と cells.csv
  → EXP-002 の analysis/evaluate_controls.py でそのまま評価できる
実行:
  .venv/bin/python experiments/EXP-20260926-003/analysis/phase4_pretrained_recalc.py --check-only   # CPU: 案 A の重心で案 A の値を再現するか
  .venv/bin/python experiments/EXP-20260926-003/analysis/phase4_pretrained_recalc.py               # MPS: 重心の再抽出 + 再計算
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
A = ROOT / "experiments/EXP-20260926-002/results"
E = ROOT / "experiments/EXP-20260926-003"
OUT = E / "results-pretrained-matched-recalc"
CFG = E / "config-pretrained-matched.yaml"
# scripts/geneformer/ が本物の geneformer パッケージを隠さないよう、run_isp は明示パスで読む
sys.path = [p for p in sys.path if not p.endswith("/scripts")]
sys.path.insert(0, str(ROOT / "scripts/geneformer"))


def unit(x):
    return x / np.linalg.norm(x, axis=-1, keepdims=True)


def recompute(run_dir: Path, goal: np.ndarray, cells: pd.DataFrame) -> pd.Series:
    g = unit(goal.astype(np.float64))
    out = pd.Series(np.nan, index=cells.index)
    for gene, s in cells.groupby("gene", sort=False):
        o = unit(np.load(run_dir / "genes" / gene / "emb" / "orig.npy").astype(np.float64))
        p = unit(np.load(run_dir / "genes" / gene / "emb" / "pert.npy").astype(np.float64))
        per = pd.read_csv(run_dir / "genes" / gene / "cells.csv", dtype={"donor_id": str})
        assert len(o) == len(s) and (per.cell_uid.to_numpy() == s.cell_uid.to_numpy()).all()
        out[s.index] = p @ g - o @ g
    return out


def extract_goal(cfg: dict, seed: int) -> np.ndarray:
    import run_isp as ri
    os.environ["GENEFORMER_DEVICE"] = ri.resolve_device(cfg)
    ri.apply_dtype(cfg)
    tok = ri.rel(cfg["data"]["tokenized_dataset"])
    work = E / "work" / "phase4"
    work.mkdir(parents=True, exist_ok=True)
    f = ri.state_emb_dataset(tok, cfg, seed, work / ri.emb_name(cfg, seed))
    embs = ri.get_state_embs(cfg, f, work / f"state_embs-seed{seed}", -1, max_ncells=None)
    v = embs[cfg["states"]["goal_state"]]
    return np.asarray(v.float().cpu() if hasattr(v, "float") else v).ravel()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check-only", action="store_true")
    a = ap.parse_args()
    cfg = yaml.safe_load(open(CFG))
    cells_all = pd.read_csv(A / "cells.csv", dtype={"donor_id": str})
    report = {"source": str(A.relative_to(ROOT)), "config": str(CFG.relative_to(ROOT)), "runs": {}}
    parts = []
    for run in sorted(cells_all.run.unique()):
        rd = A / run
        cells = cells_all[cells_all.run == run].copy()
        st = dict(np.load(rd / "state_emb_vectors.npz"))
        # 再現チェック: 案 A の goal 重心で案 A の shift_goal を再現できるか
        re = recompute(rd, st["PT_healthy"], cells)
        diff = float(np.abs(re - cells.shift_goal).max())
        rep = {"reproduce_max_abs_diff": diff}
        if not a.check_only:
            seed = int(run.split("-")[1])
            goal = extract_goal(cfg, seed)
            rep["cos_goal_A_vs_heldout"] = float(unit(goal.astype(np.float64)) @ unit(st["PT_healthy"].astype(np.float64)))
            new = cells.copy()
            new["shift_goal_orig_A"] = new.shift_goal
            new["shift_goal"] = recompute(rd, goal, cells).to_numpy()
            od = OUT / run
            od.mkdir(parents=True, exist_ok=True)
            if not (od / "genes").exists():
                (od / "genes").symlink_to(rd / "genes")
            for f in ("negatives.json",):
                if (rd / f).exists():
                    (od / f).write_text((rd / f).read_text())
            np.savez(od / "state_emb_vectors.npz", **{**st, "PT_healthy": goal.astype(np.float32)})
            parts.append(new)
        report["runs"][run] = rep
        print(run, json.dumps(rep))
    if not a.check_only:
        pd.concat(parts, ignore_index=True).to_csv(OUT / "cells.csv", index=False)
        (OUT / "phase4_report.json").write_text(json.dumps(report, indent=1) + "\n")
        print("saved", OUT)


if __name__ == "__main__":
    main()

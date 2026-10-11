#!/usr/bin/env python3
"""EXP-20261007-001 §5.4.1 — float32 で spark-cf5b と spark-c083 の数値が一致するかを確かめる。

同じ設定（work/config-dtype-fp32.yaml、対照 12 個 × 本番の 3 群の細胞）を 2 台で回した結果
results-dtype-fp32（cf5b）と results-dtype-fp32-c083（c083）を比べる。**摂動は回さない。**

合否の基準（結果を見る前に plan.md §5.4.1 に固定）:
  M-4: 遺伝子 × 群の D_raw の差 ÷ その群の 12 個の D_raw の標準偏差 の最大が 0.05 以下
  前提: compare_conditions.py が終了コード 1 でない、両機体とも float32・tf32 が閉じている、M-1（細胞の並び）が完全一致
M-2・M-3 は記録する。

走らせ方:
    .venv/bin/python experiments/EXP-20261007-001/analysis/machine_check_011.py
終了コード: 0 = 合格 / 2 = 不合格 / それ以外 = 実行エラー
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261007-001"
sys.path.insert(0, str(ROOT / "experiments/EXP-20261004-001/analysis"))
sys.path.insert(0, str(ROOT / "experiments/EXP-20261001-001/analysis"))
sys.path.insert(0, str(ROOT / "scripts/analysis"))
import evaluate_009 as ev9  # noqa: E402

LIMIT = 0.05
A, B = "results-dtype-fp32", "results-dtype-fp32-c083"


def d_raw(res: str) -> pd.DataFrame:
    cfg = yaml.safe_load((E / "work/config-dtype-fp32.yaml").read_text())
    ev9.E = E
    ev9.BLOCKED = {e: e for e in yaml.safe_load((E / "config.yaml").read_text())["blocked_genes"]}
    c, _ = ev9.load_cells(cfg, res)
    genes = list(cfg["controls"]["negative_sampler"]["genes"])
    fr = []
    for s in cfg["robustness"]["seeds"]:
        cs = c[c.seed == s].assign(m_goal=lambda d: d.m_proj)
        nl = {h: hd for h, hd in cs[cs.gene.isin(genes)].groupby("gene")}
        fr.append(ev9.gene_stats_stratified(cs, nl, "m_goal").assign(seed=s))
    return pd.concat(fr, ignore_index=True)


def main() -> int:
    out: dict = {"実験": "EXP-20261007-001", "節": "plan.md §5.4.1", "基準": f"M-4 の最大 ≤ {LIMIT}"}
    r = subprocess.run([sys.executable, str(ROOT / "scripts/analysis/compare_conditions.py"), str(E / A), str(E / B)],
                       capture_output=True, text=True)
    out["実行条件の突き合わせ"] = {"終了コード": r.returncode, "出力": r.stdout[-3000:]}
    dev = {}
    for res in (A, B):
        j = json.loads((E / res / "device.json").read_text())
        n = j.get("numerics") or {}
        dev[res] = {"host": j.get("host"), "dtype": j.get("dtype"),
                    "cuda_matmul_allow_tf32": n.get("cuda_matmul_allow_tf32"),
                    "cudnn_allow_tf32": n.get("cudnn_allow_tf32"),
                    "float32_matmul_precision": n.get("float32_matmul_precision")}
    out["device.json"] = dev
    pre_ok = r.returncode != 1 and all(v["dtype"] == "float32" and v["cuda_matmul_allow_tf32"] is False
                                       and v["float32_matmul_precision"] == "highest" for v in dev.values())
    out["前提"] = pre_ok
    # M-1・M-2・M-3
    m1, m2 = True, 0.0
    m3 = {}
    for run in sorted(p.name for p in (E / A).glob("run-*-L-1")):
        za, zb = np.load(E / A / run / "state_emb_vectors.npz"), np.load(E / B / run / "state_emb_vectors.npz")
        m3[run] = {k: float(np.max(np.abs(za[k].astype(float) - zb[k].astype(float)))) for k in ("FIB", "aStr", "dFIB")}
        for g in sorted(p.name for p in (E / A / run / "genes").iterdir() if p.is_dir()):
            ca = pd.read_csv(E / A / run / "genes" / g / "cells.csv")
            cb = pd.read_csv(E / B / run / "genes" / g / "cells.csv")
            m1 &= ca.cell_uid.tolist() == cb.cell_uid.tolist()
            for f in ("orig.npy", "pert.npy"):
                m2 = max(m2, float(np.max(np.abs(np.load(E / A / run / "genes" / g / "emb" / f).astype(float)
                                                 - np.load(E / B / run / "genes" / g / "emb" / f).astype(float)))))
    out["M-1 摂動に使った細胞の並びが完全一致"] = bool(m1)
    out["M-2 埋め込みの最大の絶対差（記録）"] = m2
    out["M-3 重心の最大の絶対差（記録）"] = m3
    a, b = d_raw(A), d_raw(B)
    m = a.merge(b, on=["gene", "seed"], suffixes=("_a", "_b"))
    sd = a.groupby("seed").D_raw.std(ddof=1)
    m["差"] = (m.D_raw_a - m.D_raw_b).abs() / m.seed.map(sd)
    out["M-4 主指標の差（標準得点の単位）"] = {"個数": int(len(m)), "平均": float(m["差"].mean()),
                                         "最大": float(m["差"].max())}
    ok = bool(pre_ok and m1 and len(m) == 36 and float(m["差"].max()) <= LIMIT)
    out["合格"] = ok
    (E / "work/machine_check_011.json").write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "実行条件の突き合わせ"}, ensure_ascii=False, indent=1, default=float))
    print(out["実行条件の突き合わせ"]["出力"][-1500:])
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())

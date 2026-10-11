#!/usr/bin/env python3
"""EXP-20261004-001 — C-11 機体の一致確認（**本番の前に通す。不合格なら本番を開始しない**）。

## なぜこれまでの機体確認（EXP-20260928-002 以降）と形が違うか

これまでは「新しい機体 1 台で回し直し、既存の保存済み結果と突き合わせる」形だった。
この実験は `cell_selection.mode: common`（全遺伝子を同じ細胞集合で過剰発現）を**本番で初めて使う**
ので、**突き合わせる既存の結果が無い**。

→ **cf5b と c083 で同じ `work/config-machinecheck.yaml`（対照 500 個の先頭 12 個 × 共通 200 細胞 ×
seed81）を独立に実行し、2 つの出力を直接突き合わせる。**

比べる前に `scripts/analysis/compare_conditions.py` を通し、**実行条件（dtype・注意機構の実装・
tf32・torch/transformers の版）が記録されているか**を確かめる（2026-10-04 追加。
「記録なし」どうしを一致と見なさない。Mac 照合が落ちた原因を bfloat16 と確かめず機体差と
断定した失敗の再発防止）。**条件が違うこと自体は止める理由にしない**
（2026-10-04、compare_conditions.py の設計を訂正。終了コード 1 ＝記録が無く突き合わせられない
ときだけ止める。終了コード 2 ＝条件が違うが数値を比べられるときは、H-1〜H-4 の数値だけで合否を決める）。

| # | 見るもの | 合格条件 |
|---|---|---|
| H-1 | 摂動前後の埋め込み（細胞ごと・成分ごと） | 最大の絶対差 0.0 |
| H-2 | 摂動に使った細胞の構成（cell_uid の集合と順序） | 完全一致 |
| H-3 | 状態の重心 2 本（FIB_MFIB・aStr。主指標が使う本） | 最大の絶対差 0.0 |
| H-4 | 主指標（向きへの射影）の遺伝子ごとの中央値 | 最大の絶対差 0.0 |
| （記録） | dFIB・dMFIB の重心 | 差を記録するが止めない（plan.md §5.5.4） |
| （記録） | 実行条件（torch・transformers の版） | 違いを記録する。止める条件にはしない |

使い方:
    .venv/bin/python experiments/EXP-20261004-001/analysis/machine_check_009.py \
        --a experiments/EXP-20261004-001/results-machinecheck-cf5b \
        --b experiments/EXP-20261004-001/results-machinecheck-c083
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261004-001"

TOL = 0.0  # 主指標が使う 2 本・埋め込み・主指標は最大差 0.0（bfloat16 でも同一カーネルなら一致する）


def main(a: str, b: str, seed: int = 81) -> None:
    pa, pb = ROOT / a, ROOT / b
    run = f"run-{seed}-L-1"
    out: dict = {"A": str(pa), "B": str(pb), "seed": seed,
                 "目的": "cell_selection.mode: common を 2 台で独立に回し、直接突き合わせる"}

    # ---- 条件がそろっているか（compare_conditions.py。先に通す）
    cc = subprocess.run([sys.executable, str(ROOT / "scripts/analysis/compare_conditions.py"),
                        str(pa), str(pb)], capture_output=True, text=True)
    # **終了コードの意味（2026-10-04 に compare_conditions.py の設計を直した）**:
    #   0 = 判定に効く条件がすべて一致 / 1 = 記録が無く突き合わせられない（止める）/
    #   2 = 条件が違う（止めない。H-1〜H-4 の数値そのものを合否の基準にし、条件の違いは記録する）。
    # 条件が違うこと自体は止める理由にしない。数値（H-1〜H-4）が一致するかどうかで判定する。
    out["0. 実行条件の突き合わせ（compare_conditions.py）"] = {
        "exit_code": cc.returncode, "stdout": cc.stdout,
        "記録が無く突き合わせられない（止める）": cc.returncode == 1,
        "条件が違う（記録するが止めない。数値で判定する）": cc.returncode == 2}
    if cc.returncode == 1:
        out["pass"] = False
        out["判定"] = "**不合格。実行条件が記録されておらず突き合わせられない。本番を開始しない**"
        (E / "work/machine_check_009.json").write_text(
            json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
        print(json.dumps(out, ensure_ascii=False, indent=1, default=float))
        raise SystemExit("**不合格（記録なし）。本番を開始しない。**")

    ca = pd.read_csv(pa / "cells.csv", dtype={"donor_id": str})
    cb = pd.read_csv(pb / "cells.csv", dtype={"donor_id": str})
    ca, cb = ca[ca.run == run], cb[cb.run == run]
    genes = sorted(set(ca.gene) & set(cb.gene))
    out["突き合わせた遺伝子"] = {"n": len(genes), "A のみ": sorted(set(ca.gene) - set(cb.gene)),
                                 "B のみ": sorted(set(cb.gene) - set(ca.gene))}

    # ---- H-3 状態の重心（主指標が使う 2 本だけを止める条件に。他は記録）
    za, zb = np.load(pa / run / "state_emb_vectors.npz"), np.load(pb / run / "state_emb_vectors.npz")
    d3 = {k: float(np.max(np.abs(za[k].astype(float) - zb[k].astype(float))))
          for k in za.files if k in zb.files}
    stop_states = ["FIB_MFIB", "aStr"]
    out["H-3 状態の重心"] = {
        "重心ごとの最大差": d3,
        "止める対象（主指標が使う2本）": {k: d3.get(k) for k in stop_states},
        "記録のみ（他の重心）": {k: v for k, v in d3.items() if k not in stop_states},
        "pass": bool(all(d3.get(k, np.inf) <= TOL for k in stop_states))}

    # ---- H-1 / H-2 / H-4
    h1, h2, h4 = [], [], []
    vg, vs = za["aStr"].astype(np.float64), za["FIB_MFIB"].astype(np.float64)
    u = vg / np.linalg.norm(vg) - vs / np.linalg.norm(vs)
    u = u / np.linalg.norm(u)
    for g in genes:
        pa_c = ca[ca.gene == g].reset_index(drop=True)
        pb_c = cb[cb.gene == g].reset_index(drop=True)
        h2.append({"gene": g, "n_a": len(pa_c), "n_b": len(pb_c),
                   "cell_uid 一致": bool(len(pa_c) == len(pb_c)
                                        and (pa_c.cell_uid.to_numpy() == pb_c.cell_uid.to_numpy()).all())})
        oa = np.load(pa / run / "genes" / g / "emb/orig.npy").astype(np.float64)
        qa = np.load(pa / run / "genes" / g / "emb/pert.npy").astype(np.float64)
        ob = np.load(pb / run / "genes" / g / "emb/orig.npy").astype(np.float64)
        qb = np.load(pb / run / "genes" / g / "emb/pert.npy").astype(np.float64)
        if oa.shape != ob.shape:
            h1.append({"gene": g, "形が違う": [list(oa.shape), list(ob.shape)]})
            continue
        h1.append({"gene": g, "orig 最大差": float(np.max(np.abs(oa - ob))),
                   "pert 最大差": float(np.max(np.abs(qa - qb)))})
        ma = (qa - oa) @ u / np.linalg.norm(qa - oa, axis=1)
        mb = (qb - ob) @ u / np.linalg.norm(qb - ob, axis=1)
        h4.append({"gene": g, "主指標の中央値の差": float(abs(np.median(ma) - np.median(mb)))})

    out["H-1 埋め込み"] = {
        "遺伝子ごと": h1,
        "最大差": max([max(x.get("orig 最大差", 0), x.get("pert 最大差", 0)) for x in h1] or [np.inf]),
        "pass": bool(all("形が違う" not in x and max(x["orig 最大差"], x["pert 最大差"]) <= TOL
                         for x in h1))}
    out["H-2 摂動に使った細胞の構成"] = {"遺伝子ごと": h2, "pass": bool(all(x["cell_uid 一致"] for x in h2))}
    out["H-4 主指標"] = {"遺伝子ごと": h4,
                         "最大差": max([x["主指標の中央値の差"] for x in h4] or [np.inf]),
                         "pass": bool(all(x["主指標の中央値の差"] <= TOL for x in h4))}
    # **合否は H-1〜H-4 の数値そのものだけで決める**（2026-10-04 に設計を訂正）。
    # 実行条件（torch・transformers の版）が違うことは、合否の条件ではなく記録する事実である。
    # この手続きが守ろうとしているのは「条件が違うせいで数値がずれていることに気づかないまま進むこと」を
    # 防ぐことであり、数値が一致しているのに条件が違うというだけで止めるのは、その判断から外れている。
    out["pass"] = bool(out["H-1 埋め込み"]["pass"] and out["H-2 摂動に使った細胞の構成"]["pass"]
                       and out["H-3 状態の重心"]["pass"] and out["H-4 主指標"]["pass"])
    cond_differs = out["0. 実行条件の突き合わせ（compare_conditions.py）"]["条件が違う（記録するが止めない。数値で判定する）"]
    out["判定"] = ("合格。cf5b と c083 は同じ尺度で測れる" if out["pass"]
                   else "**不合格。原因を重み・計算バックエンド・ライブラリ版・数値型の1つずつ切り分けてから報告する"
                        "（推測で『機体が違うから』と書かない）**")
    if out["pass"] and cond_differs:
        out["判定"] += ("。ただし実行条件（torch・transformers の版）は 2 台で違う。"
                        "この観測は 12 遺伝子・1 seed に限られ、違いが本番（626 遺伝子・3 seed）でも"
                        "効かないことを証明したものではない（plan.md §10 に限界として明記）")

    (E / "work/machine_check_009.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps({k: out[k] for k in out if k not in ("H-1 埋め込み",)}
                     | {"H-1 埋め込み": {"pass": out["H-1 埋め込み"]["pass"],
                                        "最大差": out["H-1 埋め込み"]["最大差"]}},
                     ensure_ascii=False, indent=1, default=float))
    if not out["pass"]:
        raise SystemExit("**不合格。本番を開始しない。Labretary へ報告して指示を待つ**")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True)
    ap.add_argument("--b", required=True)
    ap.add_argument("--seed", type=int, default=81)
    main(**{k: v for k, v in vars(ap.parse_args()).items()})

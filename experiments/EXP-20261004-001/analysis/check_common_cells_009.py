#!/usr/bin/env python3
"""EXP-20261004-001 — **摂動した細胞の顔ぶれが、遺伝子によらず同一か**を機械的に確かめる（検査 A-2）。

この実験が直そうとしている弱点そのものを見る検査なので、**守ろうとしている判断から直接導いている**。

| 性質 | 値 |
|---|---|
| きれいなデータでの合格率 | **1.000**（common モードは全遺伝子が同一ファイルを読むので、設計上必ず通る）|
| 鳴る条件 | 顔ぶれ・並び・提供者・長さのどれかが 1 つでも食い違ったとき（許容 0。閾値を持たない）|
| 閾値ちょうどのときの歪み | **閾値が無い**（完全一致か否かの 2 値）。判断への影響は「顔ぶれが違えば候補と対照を同じ尺度で比べられない」= 実験の前提そのもの |

**欠測は不一致として扱う**（guardrails 2026-10-04: 「記録なし」どうしを一致と見なさない）。
cells.csv が無い遺伝子が 1 つでもあれば不合格。

使い方:
    .venv/bin/python experiments/EXP-20261004-001/analysis/check_common_cells_009.py \
        [--results experiments/EXP-20261004-001/results] [--out work/common_cells_check_009.json] \
        [--config experiments/EXP-20261004-001/config.yaml]
    終了コード 0 = 合格 / 2 = 不合格。**不合格なら本番を止める。**

`--config` は既定で本番の config.yaml を読む。試し実行（C-9/C-10）の出力は細胞数が
本番と違う（30 個）ので、`work/config-smoke.yaml` を渡して n_want をそちらから読む
（2026-10-04、evaluate_009.py --smoke から呼ぶために追加。既定の挙動は変えていない）。
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261004-001"


def digest(df: pd.DataFrame) -> str:
    """cell_uid の並び・donor・長さをまとめて 1 つの指紋にする（並べ替えない）。"""
    s = "\n".join(f"{u}\t{d}\t{l}" for u, d, l in zip(df.cell_uid, df.donor_id, df.length))
    return hashlib.sha256(s.encode()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=str(E / "results"))
    ap.add_argument("--out", default=str(E / "work/common_cells_check_009.json"))
    ap.add_argument("--config", default=str(E / "config.yaml"),
                    help="n_want（max_cells_per_gene）をここから読む。既定は本番の config.yaml")
    a = ap.parse_args()
    cfg = yaml.safe_load(Path(a.config).read_text())
    if cfg["cell_selection"].get("mode") != "common":
        raise SystemExit("config の cell_selection.mode が common ではありません")
    n_want = int(cfg["cell_selection"]["max_cells_per_gene"])
    res = Path(a.results)

    out: dict = {"実験": "EXP-20261004-001",
                 "検査": "摂動した細胞の顔ぶれが遺伝子によらず同一か",
                 "見るもの": "各遺伝子の cells.csv の cell_uid・donor_id・length を、並びを変えずに突き合わせる",
                 "seed ごと": {}}
    ok = True
    for run in sorted(p for p in res.glob("run-*-L-1") if p.is_dir()):
        gdirs = sorted(p for p in (run / "genes").iterdir() if p.is_dir())
        missing = [p.name for p in gdirs if not (p / "cells.csv").exists()]
        digs: dict[str, list[str]] = {}
        sizes = set()
        for p in gdirs:
            if not (p / "cells.csv").exists():
                continue
            df = pd.read_csv(p / "cells.csv", dtype={"donor_id": str})
            digs.setdefault(digest(df), []).append(p.name)
            sizes.add(len(df))
        rec = {
            "遺伝子数": len(gdirs),
            "cells.csv が無い遺伝子": missing,
            "異なる顔ぶれの種類": len(digs),
            "細胞数の種類": sorted(sizes),
            "細胞数が設定どおりか": sizes == {n_want},
        }
        if len(digs) > 1:
            items = sorted(digs.items(), key=lambda kv: -len(kv[1]))
            rec["食い違った遺伝子（多数派以外）"] = {k[:12]: v[:10] for k, v in items[1:]}
        rec["pass"] = bool(not missing and len(digs) == 1 and sizes == {n_want} and gdirs)
        ok &= rec["pass"]
        out["seed ごと"][run.name] = rec

    # seed どうしは**違っていてよい**（seed ごとに引き直すのが設計）。違わなければ逆に異常
    ds = {k: v["異なる顔ぶれの種類"] for k, v in out["seed ごと"].items()}
    out["注"] = ("seed が違えば細胞集合も違う（引き直す設計）。同一であるべきなのは"
                 "**同じ seed の中で遺伝子をまたいだとき**だけである")
    out["pass"] = bool(ok and out["seed ごと"])
    Path(a.out).write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps(out, ensure_ascii=False, indent=1))
    sys.exit(0 if out["pass"] else 2)


if __name__ == "__main__":
    main()

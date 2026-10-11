#!/usr/bin/env python3
"""EXP-20261001-001 — トークン化済み入力を自分の work/ に用意する。

**摂動は 1 回も行わない。** 行うのは次の 2 つだけである。

1. EXP-20260930-001 のトークン化済み入力（線維芽細胞系列 21,769 細胞、著者注釈の状態ラベル
   `fib_anno_state` 付き）の**データ部分だけ**を自分の work/ へ複製する。
   `cache-*.arrow` は複製しない（複製すると元の実験のディレクトリが変更される。guardrails 2026-09-26）。
2. 遺伝子記号 → Ensembl ID の対応表を複製する。

**入力そのものは EXP-20260930-001 と同一である。** 目標・出発点・退化の 2 状態の定義も、
細胞数・提供者数も、1 つも変えていない（HYP-20261001-01 §2、選定手順 §1 の「変えないもの」）。
本実験が変えるのは **(a) 主指標の形、(b) 介入する遺伝子セット、(c) 対照遺伝子群** の 3 点だけである。

出力:
  work/tokenized.dataset          状態ラベル付きのトークン化済み入力（複製）
  work/input_states_005.json      状態ごとの細胞数・提供者数（仮説 §2 の実測値と照合して停止する）
  work/symbol_to_ensembl_005.json 遺伝子記号 → Ensembl ID
  work/source_files_005.json      複製元のファイルの sha256
"""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20261001-001"
SRC = ROOT / "experiments/EXP-20260930-001/work"          # **読み取りのみ**

# 仮説 HYP-20261001-01 §2 が実データで数えた値。**一致しなければ止める**
EXPECT = {"aStr": (5131, 77), "FIB_MFIB": (15477, 79), "dFIB": (534, 60),
          "dMFIB": (558, 14), "unused": (69, 23)}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main() -> None:
    from datasets import load_from_disk

    (E / "work").mkdir(parents=True, exist_ok=True)
    src_ds = SRC / "tokenized.dataset"
    out = E / "work/tokenized.dataset"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    # データ部分だけを APFS のクローンで複製する（cache-*.arrow は複製しない）
    names = ["data-00000-of-00001.arrow", "dataset_info.json", "state.json"]
    prov = {}
    for n in names:
        subprocess.run(["cp", "-c", str(src_ds / n), str(out / n)], check=True)
        prov[f"tokenized.dataset/{n}"] = {"from": str(src_ds / n), "sha256": sha256(out / n)}
        if sha256(src_ds / n) != prov[f"tokenized.dataset/{n}"]["sha256"]:
            sys.exit(f"複製の sha256 が元と違う: {n}")

    ds = load_from_disk(str(out))
    st = np.array([str(x) for x in ds["fib_anno_state"]])
    donors = np.array([str(x) for x in ds["donor_id"]])
    counts = {k: (int((st == k).sum()), int(len({d for d, s in zip(donors, st) if s == k})))
              for k in EXPECT}
    bad = {k: (counts[k], EXPECT[k]) for k in EXPECT if counts[k] != EXPECT[k]}
    if bad:
        sys.exit(f"仮説 §2 の実測値と一致しない: {bad}")

    # 遺伝子記号 → Ensembl ID（EXP-20260930-001 が作ったものをそのまま使う）
    m = json.loads((SRC / "symbol_to_ensembl_004.json").read_text())
    (E / "work/symbol_to_ensembl_005.json").write_text(json.dumps(m, ensure_ascii=False) + "\n")
    prov["symbol_to_ensembl"] = {"from": str(SRC / "symbol_to_ensembl_004.json"),
                                 "sha256": sha256(SRC / "symbol_to_ensembl_004.json")}

    (E / "work/source_files_005.json").write_text(
        json.dumps(prov, indent=1, ensure_ascii=False) + "\n")
    (E / "work/input_states_005.json").write_text(json.dumps(
        {"state_key": "fib_anno_state",
         "note": "入力は EXP-20260930-001 と同一（目標・出発点・退化の 2 状態の定義を変えていない）",
         "counts": {k: {"n_cells": v[0], "n_donors": v[1]} for k, v in counts.items()},
         "n_cells_total": int(len(st)), "n_symbols": len(m)},
        indent=1, ensure_ascii=False) + "\n")
    print("状態ラベル:", json.dumps(counts, ensure_ascii=False))
    print("記号 →", len(m), "個の Ensembl ID")


if __name__ == "__main__":
    main()

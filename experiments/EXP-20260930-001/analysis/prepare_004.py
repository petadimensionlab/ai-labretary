#!/usr/bin/env python3
"""EXP-20260930-001 — トークン化済み入力に「著者の注釈による状態ラベル」を付ける。

**摂動は 1 回も行わない。** 行うのは次の 3 つだけである。

1. EXP-20260929-001 のトークン化済み入力（21,769 細胞 × 線維芽細胞系列）の**データ部分だけ**を
   自分の work/ へ複製する（`cache-*.arrow` は複製しない。guardrails 2026-09-26）。
2. 複製した入力に、Lake ら 2023 の注釈（`subclass_l2`）から作った状態ラベルの列
   `fib_anno_state` を足す。**行の順序と `input_ids` は 1 つも変えない**ので、
   行番号で書かれた中間ファイル（`state_rows-seed*.json`）がそのまま使える。
3. 遺伝子記号 → Ensembl ID の対応表を作る（`work/input.h5ad` は**読み取りのみ**）。

状態ラベルの作り方（手順書 `outputs/design/2026-09-30-annotation-goal-state/procedure-v2.md` §3）:

    aStr      = aFIB + MYOF   目標（活性化した間質の状態。著者の「活性化・不適応・修復」の種別）
    FIB_MFIB  = FIB + M-FIB   出発点（著者の「参照」の種別）。**摂動するのはこの細胞だけ**
    dFIB      = dFIB          退化した線維芽細胞（記述と着手前チェックの分母。**主指標では使わない**）
    dMFIB     = dM-FIB        髄質の退化した線維芽細胞（同上）
    unused    = cycMYOF       使わない（69 細胞）

出力:
  work/tokenized.dataset          状態ラベルを足したトークン化済み入力
  work/input_states_004.json      状態ごとの細胞数・提供者数（仮説 §2 の実測値と照合して停止する）
  work/symbol_to_ensembl_004.json 遺伝子記号 → Ensembl ID
  work/source_files_004.json      複製元のファイルの sha256（どこから来たかを残す）
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
E = ROOT / "experiments/EXP-20260930-001"
SRC = ROOT / "experiments/EXP-20260929-001/work"          # **読み取りのみ**

# 状態ラベル（手順書 §3 で固定。**変えない**）
STATE_OF = {"aFIB": "aStr", "MYOF": "aStr",
            "FIB": "FIB_MFIB", "M-FIB": "FIB_MFIB",
            "dFIB": "dFIB", "dM-FIB": "dMFIB",
            "cycMYOF": "unused"}

# 仮説 HYP-20260930-01 §2 が実データで数えた値。**一致しなければ止める**
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
    tmp = E / "work/tokenized-src.dataset"
    out = E / "work/tokenized.dataset"
    if tmp.exists():
        shutil.rmtree(tmp)
    if out.exists():
        shutil.rmtree(out)
    tmp.mkdir(parents=True)

    # データ部分だけを APFS のクローンで複製する（cache-*.arrow は複製しない）
    names = ["data-00000-of-00001.arrow", "dataset_info.json", "state.json"]
    prov = {}
    for n in names:
        subprocess.run(["cp", "-c", str(src_ds / n), str(tmp / n)], check=True)
        prov[f"tokenized.dataset/{n}"] = {"from": str(src_ds / n), "sha256": sha256(tmp / n)}
        if sha256(src_ds / n) != prov[f"tokenized.dataset/{n}"]["sha256"]:
            sys.exit(f"複製の sha256 が元と違う: {n}")

    ds = load_from_disk(str(tmp))
    sub = np.array([str(x) for x in ds["subclass_l2"]])
    donors = np.array([str(x) for x in ds["donor_id"]])
    unknown = sorted(set(sub) - set(STATE_OF))
    if unknown:
        sys.exit(f"状態ラベルの対応が無い subclass_l2 がある: {unknown}")
    state = np.array([STATE_OF[s] for s in sub])

    counts = {k: (int((state == k).sum()), int(len({d for d, s in zip(donors, state) if s == k})))
              for k in EXPECT}
    bad = {k: (counts[k], EXPECT[k]) for k in EXPECT if counts[k] != EXPECT[k]}
    if bad:
        sys.exit(f"仮説 §2 の実測値と一致しない: {bad}")

    ds2 = ds.add_column("fib_anno_state", state.tolist())
    # 行の順序と input_ids が変わっていないことを検算する
    assert ds2["a1_row"] == ds["a1_row"], "行の順序が変わった"
    assert [len(x) for x in ds2["length"]] if False else ds2["length"] == ds["length"], "長さが変わった"
    ds2.save_to_disk(str(out))
    shutil.rmtree(tmp)

    # 遺伝子記号 → Ensembl ID（input.h5ad は読み取りのみ）
    import anndata as ad
    a = ad.read_h5ad(SRC / "input.h5ad", backed="r")
    ens = a.var["ensembl_id"].astype(str).to_numpy()
    sym = a.var["gene_name"].astype(str).to_numpy()
    a.file.close()
    m: dict[str, str] = {}
    for s, e in zip(sym, ens):
        m.setdefault(str(s), str(e))
    (E / "work/symbol_to_ensembl_004.json").write_text(json.dumps(m, ensure_ascii=False) + "\n")
    prov["input.h5ad"] = {"from": str(SRC / "input.h5ad"), "sha256": sha256(SRC / "input.h5ad")}

    (E / "work/source_files_004.json").write_text(
        json.dumps(prov, indent=1, ensure_ascii=False) + "\n")
    (E / "work/input_states_004.json").write_text(json.dumps(
        {"state_key": "fib_anno_state", "state_of_subclass_l2": STATE_OF,
         "counts": {k: {"n_cells": v[0], "n_donors": v[1]} for k, v in counts.items()},
         "n_cells_total": int(len(state)), "n_symbols": len(m)},
        indent=1, ensure_ascii=False) + "\n")
    print("状態ラベル:", json.dumps(counts, ensure_ascii=False))
    print("記号 →", len(m), "個の Ensembl ID")


if __name__ == "__main__":
    main()

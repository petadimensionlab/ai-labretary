#!/usr/bin/env python3
"""EXP-20260926-004 (c-2): 拡張 run の細胞単位の評価（plan.md §C.2「細胞単位は EXP-003 と同じ C1〜C5・C3-panel,
evaluate_controls.py の関数」）。EXP-002 の evaluate_controls.py をそのまま呼び、拡張 run 固有の 2 点だけ差し替える:
  1. C5 の発現細胞数（DET）を新しい start（aPT_CKD_TLShighmid）の detection に替える
  2. a4/d4 用の削除位置（rc.token_positions）を拡張 run の tokenized dataset で引く（EXP-002 の dataset は cell_uid が違う）
判定の定義・閾値・a2 の実装は変えない。--exclude-role de_matched で EXP-003 と同じ 34 遺伝子に揃える
（a2 の補正は null 細胞だけで推定するので、除いても他の遺伝子の値は変わらない。変わるのは C4i/C5 の遺伝子集合だけ）。
usage:
  .venv/bin/python experiments/EXP-20260926-004/analysis/eval_expanded_cells.py --out eval34 --exclude-role de_matched
  .venv/bin/python experiments/EXP-20260926-004/analysis/eval_expanded_cells.py --out eval44
"""
from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "experiments/EXP-20260926-002/analysis"))
import evaluate_controls as ec  # noqa: E402
import reanalysis_cosine as rc  # noqa: E402

EXPD = ROOT / "experiments/EXP-20260926-004/expanded"
TOK = EXPD / "work/tokenized.dataset"


def token_positions_expanded(cells: pd.DataFrame) -> pd.Series:
    from datasets import load_from_disk
    from geneformer import TOKEN_DICTIONARY_FILE
    tokd = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
    ds = load_from_disk(str(TOK))  # 本実験の work/（他の実験の dataset ではない）
    uids = sorted(cells.cell_uid.unique())
    ids = dict(zip(uids, ds.select(uids)["input_ids"]))
    pos = [ids[u].index(tokd[g]) for u, g in zip(cells.cell_uid, cells.gene)]
    return pd.Series(pos, index=cells.index)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="expanded/results-expanded/ 以下の出力ディレクトリ名")
    ap.add_argument("--exclude-role", default=None)
    a = ap.parse_args()
    ec.DET = EXPD / "work/detection-aPT_CKD_TLShighmid.json"
    assert ec.DET.exists(), ec.DET
    rc.token_positions = token_positions_expanded
    orig = ec.load_run

    def load_run(res, run):
        cells, E0, E1, cen = orig(res, run)
        # 3. ベクトル指標（vp.vector_metrics）は start 重心を 'aPT_CKD_TLShigh' の名前で引く → 拡張 run の start
        #    （aPT_CKD_TLShighmid, plan §C.2）の重心を同じ名前でも渡す（副次指標のみ。a2 は重心の名前に依らない）
        assert "aPT_CKD_TLShighmid" in cen and "aPT_CKD_TLShigh" not in cen
        cen = {**cen, "aPT_CKD_TLShigh": cen["aPT_CKD_TLShighmid"]}
        if a.exclude_role:
            keep = (cells.role != a.exclude_role).to_numpy()
            cells, E0, E1 = cells[keep].reset_index(drop=True), E0[keep], E1[keep]
        return cells, E0, E1, cen
    ec.load_run = load_run
    rel = "../EXP-20260926-004/expanded/results-expanded"  # evaluate_controls は EXP-002 からの相対パス
    sys.argv = [sys.argv[0], "--results", rel, "--out", f"{rel}/{a.out}",
                "--config", str(EXPD / "config-expanded.yaml")]
    ec.main()


if __name__ == "__main__":
    main()

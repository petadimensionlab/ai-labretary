"""EXP-20260927-001 事前確認（事前登録の一部）。ISP の結果は一切使わない（トークン化済み入力だけ）。

start 細胞（aPT_CKD_TLShigh）で、各遺伝子の「元の位置」（input_ids 内の index = 過剰発現型の摂動で
最上位へ動かす距離, 移動距離）の分布を遺伝子ごとに集計する。null 遺伝子プールの選定と
「移動距離をそろえた並べ替え検定」の実行可能性（マッチできる null の数）の確認に使う。
出力: work/pos_by_gene.parquet（遺伝子ごとの要約）
実行: .venv/bin/python experiments/EXP-20260927-001/analysis/precheck_pos.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts/geneformer"))
import run_isp as ri  # noqa: E402

E = ROOT / "experiments/EXP-20260927-001"
CFG = yaml.safe_load(open(E / "config.yaml"))
TOK = ROOT / CFG["data"]["tokenized_dataset"]


def main() -> None:
    from datasets import load_from_disk
    ds = load_from_disk(str(TOK))
    st = CFG["states"]["start_state"]
    state = ds[CFG["data"]["state_key"]]
    rows = [i for i, s in enumerate(state) if s == st]
    ids = ds.select(rows)["input_ids"]
    donors = [str(x) for x in ds.select(rows)["donor_id"]]
    tok, pos, length, don = [], [], [], []
    for k, x in enumerate(ids):
        x = np.asarray(x)
        tok.append(x); pos.append(np.arange(len(x))); length.append(np.full(len(x), len(x)))
        don.append(np.full(len(x), k))
    df = pd.DataFrame({"tok": np.concatenate(tok), "pos": np.concatenate(pos),
                       "length": np.concatenate(length), "cell": np.concatenate(don)})
    df["relpos"] = df.pos / df.length
    g = df.groupby("tok").agg(n_start=("cell", "size"), pos_med=("pos", "median"),
                              pos_q10=("pos", lambda v: np.quantile(v, 0.1)),
                              pos_q90=("pos", lambda v: np.quantile(v, 0.9)),
                              relpos_med=("relpos", "median"), len_med=("length", "median"))
    cd = pd.Series(donors).to_numpy()
    nd = df.assign(d=cd[df.cell]).groupby("tok").d.nunique()
    g["n_donor"] = nd
    inv = {v: k for k, v in ri.gene_token_dict().items()}
    g["ensembl"] = [inv.get(t) for t in g.index]
    g = g[g.ensembl.notna() & g.ensembl.astype(str).str.startswith("ENSG")]
    out = E / "work/pos_by_gene.parquet"
    g.reset_index().to_parquet(out)
    print(f"start {len(rows)} 細胞, 遺伝子 {len(g)}; 保存 {out}")


if __name__ == "__main__":
    main()

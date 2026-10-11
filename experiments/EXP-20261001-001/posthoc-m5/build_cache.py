#!/usr/bin/env python3
"""M5 の事後診断 — 細胞ごとの表と変位ベクトルの控えを作る（事前登録の外。摂動は回さない）。

**主判定（介入セット 119 個）とその入れ子の遺伝子は、読み込んだ直後にすべて捨てる。**
捨てたあとに残るのは、対照 500 個・目標の軸の陽性対照・中立対照・恒常発現だけなので、
この控えからは主判定の値を計算しようがない。

出力（すべて scratchpad。プロジェクトには置かない）:
  cells.parquet   細胞ごとの表（lp, ll, relpos, 移動量, 各状態への寄り, 射影 ほか）
  disp.f32        変位ベクトル d = pert - orig（float32, 行は cells.parquet と同じ順）
  states.json     run ごとの状態の重心ベクトル
"""
from __future__ import annotations

import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20261001-001"
OUT = Path(sys.argv[1] if len(sys.argv) > 1 else
           "/private/tmp/claude-503/-Users-petadimensionlab-workspace-research-ai-labretary/"
           "3a2e253b-7169-44cc-ba54-2207da1ae321/scratchpad/cache")
OUT.mkdir(parents=True, exist_ok=True)

BLOCKED = {"ENSG00000141510", "ENSG00000118785", "ENSG00000148346"}
# **主判定とその入れ子**（ここで捨てる）
PRIMARY_FAMILY = ("INT", "INT_S1", "INT_S2", "INT_S3", "INT_pathway_specific", "INT_plus_T1T2",
                  "T1", "T2")


def main() -> None:
    cfg = yaml.safe_load(open(E / "config.yaml"))
    res = E / cfg["output"]["results_dir"]
    gj = json.loads((E / "work/genes_005.json").read_text())["genes"]
    ens = {g["symbol"]: g["ensembl"] for g in gj}
    subs = json.loads((E / "work/subsets_005.json").read_text())["subsets"]
    forbidden = {ens[x["symbol"]] for k in PRIMARY_FAMILY for x in subs[k]}
    (OUT / "forbidden_genes.json").write_text(json.dumps(sorted(forbidden)))

    c = pd.read_csv(res / "cells.csv", dtype={"donor_id": str})
    assert not set(c.gene) & BLOCKED, "盲検中の遺伝子が入っている"
    n0 = len(c)
    c = c[~c.gene.isin(forbidden)].reset_index(drop=True)
    print(f"主判定の家族 {len(forbidden)} 遺伝子を除外: {n0} → {len(c)} 行、"
          f"遺伝子 {c.gene.nunique()} 個")
    c["seed"] = c.run.astype(str).str.extract(r"run-(\d+)", expand=False).astype(int)

    # --- 元の順位（evaluate_005.load と同じ手順）
    from datasets import load_from_disk
    from geneformer import TOKEN_DICTIONARY_FILE
    tokd = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
    ds = load_from_disk(str(E / "work/tokenized.dataset"))
    uids = sorted(c.cell_uid.unique())
    ids = {u: np.asarray(x) for u, x in zip(uids, ds.select(uids)["input_ids"])}
    pos = np.empty(len(c), dtype=np.int64)
    for gene, part in c.groupby("gene", sort=False):
        t = tokd[gene]
        pos[part.index.to_numpy()] = [int(np.flatnonzero(ids[u] == t)[0]) for u in part.cell_uid]
    c["pos"] = pos
    c["lp"] = np.log2(1 + c.pos.astype(float))
    c["ll"] = np.log2(c.length.astype(float))
    c["relpos"] = c.pos / c.length

    s_ = cfg["states"]
    names = {"goal": s_["goal_state"], "start": s_["start_state"]}
    names |= {a: a for a in s_.get("alt_states", [])}
    keys = list(names)

    N, D = len(c), 1152
    disp = np.lib.format.open_memmap(OUT / "disp.npy", mode="w+", dtype=np.float32, shape=(N, D))
    cols = {f"s_{k}": np.full(N, np.nan) for k in keys}          # cos(p,v) - cos(o,v)
    cols |= {f"dv_{k}": np.full(N, np.nan) for k in keys}        # d . v̂
    cols |= {f"ov_{k}": np.full(N, np.nan) for k in keys}        # o . v̂
    cols["dnorm"] = np.full(N, np.nan)
    cols["onorm"] = np.full(N, np.nan)
    cols["pnorm"] = np.full(N, np.nan)
    cols["m_proj"] = np.full(N, np.nan)
    cols["du"] = np.full(N, np.nan)                              # d . û
    cols["ou"] = np.full(N, np.nan)                              # o . û

    states = {}
    t0 = time.time()
    done = 0
    for (run, gene), part in c.groupby(["run", "gene"], sort=False):
        sz = np.load(res / run / "state_emb_vectors.npz")
        if run not in states:
            states[run] = {k: sz[v].astype(np.float64).tolist() for k, v in names.items()}
        gd = res / run / "genes" / gene / "emb"
        o = np.load(gd / "orig.npy").astype(np.float64)
        p = np.load(gd / "pert.npy").astype(np.float64)
        assert len(o) == len(part), f"{gd}: {len(o)} != {len(part)}"
        idx = part.index.to_numpy()
        d = p - o
        disp[idx] = d.astype(np.float32)
        on = np.linalg.norm(o, axis=1)
        pn = np.linalg.norm(p, axis=1)
        dn = np.linalg.norm(d, axis=1)
        cols["onorm"][idx], cols["pnorm"][idx], cols["dnorm"][idx] = on, pn, dn
        for k, cname in names.items():
            v = sz[cname].astype(np.float64)
            vh = v / np.linalg.norm(v)
            cols[f"s_{k}"][idx] = (p @ vh) / pn - (o @ vh) / on
            cols[f"dv_{k}"][idx] = d @ vh
            cols[f"ov_{k}"][idx] = o @ vh
        vg = sz[names["goal"]].astype(np.float64)
        vs = sz[names["start"]].astype(np.float64)
        u = vg / np.linalg.norm(vg) - vs / np.linalg.norm(vs)
        u = u / np.linalg.norm(u)
        cols["du"][idx] = d @ u
        cols["ou"][idx] = o @ u
        with np.errstate(divide="ignore", invalid="ignore"):
            cols["m_proj"][idx] = (d @ u) / dn
        done += 1
        if done % 300 == 0:
            print(f"  {done} 群  {time.time() - t0:.0f}s", flush=True)
    for k, v in cols.items():
        c[k] = v
    disp.flush()

    # 移動量 0 の細胞を外す（事前登録と同じ）
    zero = (c.dnorm <= 0).to_numpy()
    print(f"移動量 0 の細胞: {int(zero.sum())}")
    keep = np.flatnonzero(~zero)
    c2 = c.loc[keep].reset_index(drop=True)
    d2 = np.lib.format.open_memmap(OUT / "disp_kept.npy", mode="w+", dtype=np.float32,
                                   shape=(len(keep), D))
    step = 20000
    for i in range(0, len(keep), step):
        d2[i:i + step] = disp[keep[i:i + step]]
    d2.flush()
    del disp, d2
    (OUT / "disp.npy").unlink()
    c2.to_parquet(OUT / "cells.parquet")
    (OUT / "states.json").write_text(json.dumps(states))
    print(f"保存: {len(c2)} 行 / {c2.gene.nunique()} 遺伝子  {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()

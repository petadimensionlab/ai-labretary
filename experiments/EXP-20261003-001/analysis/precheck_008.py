#!/usr/bin/env python3
"""EXP-20261003-001 — 着手前チェック（摂動を 1 回も回さない）。

確かめること（**どれか 1 つでも落ちたら、その遺伝子を外すか、本番を開始しない**）:

  C-1  6 個それぞれが Geneformer の辞書（トークン辞書）にあるか
  C-2  出発点 FIB_MFIB で、提供者あたり 20 細胞の上限を当てた後に 120 細胞以上取れるか
       （seed 81 / 82 / 83 のすべてで。`run_isp.py: select_cells()` と同じ手順で数える）
  C-3  既存の判定セット（陽性対照 71・介入 119・対照 1,000・中立 4）と重ならないか。
       盲検中の 3 遺伝子（TP53 / SPP1 / LCN2）と一致しないか
  C-4  偶然の水準 V-2 を、**対照 1,000 個の確定済みの表から実測する**
       （z > 2.73 を 3 seed すべてで満たす対照の割合）

出力: work/precheck_008.json
実行: .venv/bin/python experiments/EXP-20261003-001/analysis/precheck_008.py
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261003-001"
P01 = ROOT / "experiments/EXP-20261001-001"
P02 = ROOT / "experiments/EXP-20261002-002"
sys.path.insert(0, str(ROOT / "scripts/geneformer"))

CANDIDATES = ["ACTA2", "TAGLN", "CNN1", "MYH11", "POSTN", "PDGFRB"]
SEEDS = [81, 82, 83]
MAX_CELLS = 200
MAX_PER_DONOR = 20
MIN_CELLS = 120
START_STATE = "FIB_MFIB"
Z_FLOOR = 2.73          # 単一遺伝子の検出下限（EXP-20260928-002 の実測。変えない）
BLOCKED = {"ENSG00000141510": "TP53", "ENSG00000118785": "SPP1", "ENSG00000148346": "LCN2"}


def main() -> None:
    import run_isp as ri
    from datasets import load_from_disk
    import scipy.sparse as sp

    sym2ens = json.loads((P01 / "work/symbol_to_ensembl_005.json").read_text())
    tokd = ri.gene_token_dict()

    out: dict = {
        "実験": "EXP-20261003-001",
        "注": "**摂動は 1 回も回していない。** 辞書・トークン化済み入力・確定済みの集計表を読むだけ",
        "候補": CANDIDATES,
        "条件": {"出発点": START_STATE, "細胞上限": MAX_CELLS, "提供者上限": MAX_PER_DONOR,
                 "細胞数の下限": MIN_CELLS, "seed": SEEDS, "検出下限 z": Z_FLOOR},
    }

    # ---------------------------------------------------------------- C-1 辞書
    c1 = {}
    for s in CANDIDATES:
        e = sym2ens.get(s)
        c1[s] = {"ensembl": e, "辞書にある": bool(e and e in tokd),
                 "token": int(tokd[e]) if e and e in tokd else None}
    out["C-1 Geneformer の辞書"] = c1
    ok1 = [s for s in CANDIDATES if c1[s]["辞書にある"]]

    # ---------------------------------------------------------------- C-3 重なり
    gj = json.loads((P01 / "work/genes_005.json").read_text())["genes"]
    ens = {g["symbol"]: g["ensembl"] for g in gj}
    subs = json.loads((P01 / "work/subsets_005.json").read_text())["subsets"]
    null500 = [ens[x["symbol"]] for x in
               json.loads((P01 / "work/null_pool_005.json").read_text())["versions"]["B_500"]]
    xtra = json.loads((P02 / "work/extra_null_007.json").read_text())["追加する対照"]
    extra_null = [x["ensembl"] for x in xtra]
    sets_named = {
        "陽性対照 71": {ens[x["symbol"]] for x in subs["positive_control71"]},
        "介入 119": {ens[x["symbol"]] for x in subs["INT"]},
        "中立 4": {ens[x["symbol"]] for x in subs["neutral_adopted"]},
        "対照 1,000": set(null500) | set(extra_null),
        "摂動済みの全 1,217": {g["ensembl"] for g in gj} | set(extra_null),
        "盲検中の 3": set(BLOCKED),
    }
    c3 = {}
    for s in CANDIDATES:
        e = c1[s]["ensembl"]
        hit = {k: (e in v) for k, v in sets_named.items()}
        c3[s] = {"重なったセット": [k for k, v in hit.items() if v], "重なり件数": sum(hit.values())}
    out["C-3 既存セットとの重なり"] = {
        "セットの大きさ": {k: len(v) for k, v in sets_named.items()}, "候補ごと": c3,
        "重なりの合計": int(sum(v["重なり件数"] for v in c3.values()))}
    ok3 = [s for s in ok1 if c3[s]["重なり件数"] == 0]

    # ---------------------------------------------------------------- C-2 細胞数
    ds = load_from_disk(str(P01 / "work/tokenized.dataset"))
    st = np.array([str(x) for x in ds["fib_anno_state"]])
    donors = np.array([str(x) for x in ds["donor_id"]])
    ids = [np.asarray(x) for x in ds["input_ids"]]
    indptr = np.cumsum([0] + [len(x) for x in ids])
    col = np.concatenate(ids)
    ntok = max(int(col.max()) + 1, max(tokd.values()) + 1)
    M = sp.csr_matrix((np.ones(len(col), bool), col, indptr),
                      shape=(len(ids), ntok)).tocsc()
    sidx = np.flatnonzero(st == START_STATE)
    Ms = M[sidx].tocsc()
    dn_s = donors[sidx]
    c2 = {}
    for s in ok3:
        e = c1[s]["ensembl"]
        t = tokd[e]
        idx = Ms.indices[Ms.indptr[t]:Ms.indptr[t + 1]].tolist()
        per = {}
        for sd in SEEDS:
            keep = ri.balanced_subset(idx, list(dn_s), MAX_CELLS, MAX_PER_DONOR,
                                      random.Random(f"{sd}-{e}"))
            per[str(sd)] = {"n_cells": len(keep), "n_donors": len({dn_s[i] for i in keep})}
        c2[s] = {"ensembl": e, "出発点で入力に載る細胞": len(idx),
                 "出発点の提供者数": len(set(dn_s[idx])) if idx else 0,
                 "per_seed": per, "最小細胞数": min(v["n_cells"] for v in per.values()),
                 "合格": min(v["n_cells"] for v in per.values()) >= MIN_CELLS}
    out["C-2 抽出できる細胞数（提供者上限 20 の後）"] = c2
    out["出発点の細胞数"] = int(len(sidx))
    adopted = [s for s in ok3 if c2[s]["合格"]]

    # ---------------------------------------------------------------- C-4 偶然の水準
    g = pd.concat([pd.read_csv(P02 / f"results/gene_seed_007-{sd}.csv") for sd in SEEDS])
    Z = g.pivot(index="gene", columns="seed", values="z")
    nl = sorted(set(null500) | set(extra_null))
    Zn = Z.reindex(nl)
    assert Zn.notna().all().all(), "対照 1,000 個の z に欠けがある"
    hit3 = (Zn > Z_FLOOR).all(axis=1)
    out["C-4 偶然の水準（対照 1,000 個で実測）"] = {
        "出どころ": "experiments/EXP-20261002-002/results/gene_seed_007-{81,82,83}.csv（確定値）",
        "対照の個数": int(len(Zn)),
        "3 seed すべてで z > 2.73 の対照": int(hit3.sum()),
        "割合 p": round(float(hit3.mean()), 5),
        "seed ごとに z > 2.73 の対照の割合":
            {str(c): round(float((Zn[c] > Z_FLOOR).mean()), 5) for c in Zn.columns},
        "候補 6 個での期待件数": round(float(hit3.mean()) * len(CANDIDATES), 4),
        "P(1 件以上 | 偶然)": round(float(1 - (1 - hit3.mean()) ** len(CANDIDATES)), 4),
        "参考: 負の側 z < -2.73 を 3 seed すべて":
            int((Zn < -Z_FLOOR).all(axis=1).sum()),
    }

    out["採用する遺伝子"] = adopted
    out["外した遺伝子"] = {s: ("辞書に無い" if s not in ok1 else
                               "既存セットと重なる" if s not in ok3 else
                               f"細胞数 {c2[s]['最小細胞数']} < {MIN_CELLS}")
                           for s in CANDIDATES if s not in adopted}
    out["着手してよいか"] = len(adopted) > 0 and "ACTA2" in adopted
    out["着手の条件"] = "**ACTA2 が採用されること。** ACTA2 が外れるならこの実験は成立しない"

    (E / "work").mkdir(parents=True, exist_ok=True)
    (E / "work/precheck_008.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "C-2 抽出できる細胞数（提供者上限 20 の後）"},
                     ensure_ascii=False, indent=1, default=float))
    print("\n--- C-2 ---")
    for s, v in c2.items():
        print(f"{s} ({v['ensembl']}): 出発点 {v['出発点で入力に載る細胞']} 細胞 / "
              f"提供者 {v['出発点の提供者数']} → 抽出 "
              f"{[v['per_seed'][str(sd)]['n_cells'] for sd in SEEDS]} "
              f"（提供者 {[v['per_seed'][str(sd)]['n_donors'] for sd in SEEDS]}）"
              f" {'合格' if v['合格'] else '**不合格**'}")
    print("\n採用:", adopted, "／外した:", out["外した遺伝子"])


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""複数の遺伝子を同時に摂動したとき、効果が大きくなるかを、既存のデータから見積もる。

## 問い（PI 2026-10-03）

> ここまでに構築した評価を用いた場合、もし複数遺伝子を同時に摂動した場合に、
> 条件をパスして効果が得られる期待はもてそうでしょうか？

## 考え方

主指標は「変位ベクトルが目標の向きとなす角の余弦」である。
遺伝子を同時に摂動したときの変位が、**個々の変位のおおよその和**になるとすると:

- **目標の向きに沿った成分**は、向きが揃っていれば **n 倍**に足し合わさる（coherent）
- **目標と無関係な成分**は、向きがばらばらなら **√n 倍**にしかならない（incoherent）

したがって **余弦は最大で √n 倍に増えうる。** ただしこれは「向きが揃っている」場合の話で、
揃っていなければ増えない。**揃っているかどうかは、手元の埋め込みから直接測れる。**

## 測るもの

1. **個々の余弦の平均** — いま測っているもの（1 遺伝子ずつ摂動した結果の平均）
2. **和の余弦** — 変位ベクトルを足し合わせてから余弦を取ったもの（＝同時摂動が
   おおよそ加法的だった場合の見積もり）
3. **増幅率** = 2 ÷ 1。これが √n に近ければ「向きが揃っている」、1 に近ければ「揃っていない」

**この見積もりは摂動を 1 回も回さずに計算できる。** ただし §限界を必ず読むこと。

## 限界（**結論に必ず添える**）

- **加法性は仮定であって、確かめていない。** 基盤モデルは非線形で、同時に多数の遺伝子を
  入力の上位へ動かせば飽和する。実際の同時摂動は、この見積もりより小さくなる可能性が高い
- 遺伝子ごとに摂動した細胞の顔ぶれが違う。本来は同じ細胞でそろえるべきである
- **これは探索的な見積もりであり、事前登録された判定ではない**

実行: `.venv/bin/python experiments/EXP-20261001-001/posthoc-combo/coherence.py`
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
E = HERE.parent
ROOT = E.parents[1]


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def main() -> None:
    cfg = yaml.safe_load(open(E / "config.yaml"))
    gj = json.loads((E / "work/genes_005.json").read_text())["genes"]
    ens = {g["symbol"]: g["ensembl"] for g in gj}
    sym = {v: k for k, v in ens.items()}
    role = {g["ensembl"]: g.get("role") for g in gj}
    subs = json.loads((E / "work/subsets_005.json").read_text())["subsets"]
    sets_def = {k: [ens[x["symbol"]] for x in v if x["symbol"] in ens] for k, v in subs.items()}

    out = {"問い": "複数遺伝子の同時摂動で効果が大きくなるか",
           "方法": "既存の埋め込みから、変位ベクトルの向きの揃い方を測る。摂動は回していない",
           "注意": "加法性は仮定。基盤モデルは非線形なので、実際はこの見積もりより小さくなりうる"}

    for seed in cfg["robustness"]["seeds"]:
        run = f"run-{seed}-L-1"
        res = E / "results" / run
        if not res.exists():
            continue
        sz = np.load(E / f"results/{run}/state_emb_vectors.npz")
        s_ = cfg["states"]
        c_start = sz[s_["start_state"]].astype(np.float64)
        c_goal = sz[s_["goal_state"]].astype(np.float64)
        # 主指標が使う向き: 単位化した重心の差
        u = unit(unit(c_goal) - unit(c_start))

        per_set = {}
        for name in ("INT", "INT_S1", "INT_S2", "INT_S3", "positive_control71", "T1", "neutral_adopted"):
            genes = sets_def.get(name)
            if not genes:
                continue
            deltas, cosines = [], []
            for g in genes:
                gd = res / "genes" / g / "emb"
                if not (gd / "orig.npy").exists():
                    continue
                o = np.load(gd / "orig.npy").astype(np.float64)
                p = np.load(gd / "pert.npy").astype(np.float64)
                d = (p - o).mean(axis=0)           # その遺伝子の平均の変位
                nd = np.linalg.norm(d)
                if nd <= 0:
                    continue
                deltas.append(d)
                cosines.append(float(d @ u) / nd)
            if len(deltas) < 2:
                continue
            D = np.vstack(deltas)
            n = len(deltas)
            mean_cos = float(np.mean(cosines))
            summed = D.sum(axis=0)
            sum_cos = float(summed @ u) / np.linalg.norm(summed)
            # 向きの揃い方そのもの（変位どうしの平均コサイン）
            U = D / np.linalg.norm(D, axis=1, keepdims=True)
            G = U @ U.T
            pair_cos = float((G.sum() - np.trace(G)) / (n * (n - 1)))
            per_set[name] = {
                "遺伝子数": n,
                "個々の余弦の平均": round(mean_cos, 4),
                "和の余弦（加法を仮定した見積もり）": round(sum_cos, 4),
                "増幅率": round(sum_cos / mean_cos, 2) if mean_cos != 0 else None,
                "√n（向きが完全に揃ったときの上限）": round(np.sqrt(n), 2),
                "変位どうしの平均コサイン（揃い方）": round(pair_cos, 4),
            }
        out[f"seed{seed}"] = per_set

    p = HERE / "coherence.json"
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    print(f"書き出し: {p}\n")

    for seed in cfg["robustness"]["seeds"]:
        k = f"seed{seed}"
        if k not in out:
            continue
        print(f"=== {k} ===")
        print("%-22s %6s %10s %10s %8s %8s %10s" %
              ("セット", "n", "個々平均", "和", "増幅率", "√n", "揃い方"))
        for name, v in out[k].items():
            print("%-22s %6d %10.4f %10.4f %8s %8.2f %10.4f" %
                  (name, v["遺伝子数"], v["個々の余弦の平均"], v["和の余弦（加法を仮定した見積もり）"],
                   v["増幅率"], v["√n（向きが完全に揃ったときの上限）"],
                   v["変位どうしの平均コサイン（揃い方）"]))
        print()


if __name__ == "__main__":
    main()

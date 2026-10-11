#!/usr/bin/env python3
"""EXP-20261004-001 — **bf16 で合否を論じられる z の下限**を、実測から決める。

## なぜ要るか

PI 指示（2026-10-04、CLAUDE.md §6.1）:「基準を決めるときは fp32、探索するときは bf16 で構わない。
ただし**検出下限の近く（z が 2.73〜3.6）の遺伝子について合否を論じるなら fp32 が要る**」。
この帯の下端 2.73 と上端 3.6 の差が、bf16 の誤差に対して十分かを数値で確かめる。

## 測り方（摂動は 1 回も回さない。保存済みの埋め込みを読むだけ）

`experiments/EXP-20261003-001` の機体・数値型の確認で、**同じ 12 遺伝子・同じ 200 細胞**を
4 条件で回した結果が残っている。そこから次の 3 つを同じ物差しで並べる。

1. **bf16 ↔ fp32**（同じ細胞）… 数値型だけの差
2. **bf16 ↔ bf16（別の機体）**（同じ細胞）… カーネルだけの差
3. **seed 違い**（同じ数値型・同じ機体）… 細胞の引き直しの差

3 の大きさは、対照 1,000 個の z の分散分解から **z 換算 0.283（単一 seed の SD）**と分かっている
（`experiments/EXP-20261003-001/results/gene_seed_008-{81,82,83}.csv`）。
1 と 3 の比を取れば、**bf16 の誤差を z の単位に換算できる**（単位の取り方に依存しない）。

出力: experiments/EXP-20261004-001/work/dtype_resolution_009.json
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
E3 = ROOT / "experiments/EXP-20261003-001"
P01 = ROOT / "experiments/EXP-20261001-001/results"
OUT = ROOT / "experiments/EXP-20261004-001/work/dtype_resolution_009.json"

CONDS = {
    "bf16 CUDA (cf5b)": E3 / "results-machinecheck-cf5b",
    "bf16 MPS (Mac)":   E3 / "results-machinecheck-mac",
    "fp32 MPS (Mac)":   E3 / "results-machinecheck-fp32",
}


def direction(root: Path, run: str) -> np.ndarray:
    z = np.load(root / run / "state_emb_vectors.npz")
    a = z["aStr"].astype(np.float64); s = z["FIB_MFIB"].astype(np.float64)
    u = a / np.linalg.norm(a) - s / np.linalg.norm(s)
    return u / np.linalg.norm(u)


def med(root: Path, run: str, gene: str, u: np.ndarray) -> float:
    o = np.load(root / run / "genes" / gene / "emb/orig.npy").astype(np.float64)
    p = np.load(root / run / "genes" / gene / "emb/pert.npy").astype(np.float64)
    d = p - o
    return float(np.median(d @ u / np.linalg.norm(d, axis=1)))


def main() -> None:
    run81 = "run-81-L-1"
    genes = sorted(p.name for p in (E3 / "results-machinecheck-fp32" / run81 / "genes").iterdir()
                   if (p / "emb/pert.npy").exists())
    genes = [g for g in genes if all((P01 / f"run-{s}-L-1" / "genes" / g / "emb/pert.npy").exists()
                                     for s in (81, 82, 83))]
    u81 = direction(P01, run81)
    m = {k: {g: med(v, run81, g, u81) for g in genes} for k, v in CONDS.items()}
    seedm = {s: {g: med(P01, f"run-{s}-L-1", g, direction(P01, f"run-{s}-L-1")) for g in genes}
             for s in (81, 82, 83)}

    d_seed = np.array([abs(seedm[a][g] - seedm[b][g]) for a, b in itertools.combinations((81, 82, 83), 2)
                       for g in genes])
    d_dtype = np.array([abs(m["bf16 MPS (Mac)"][g] - m["fp32 MPS (Mac)"][g]) for g in genes])
    d_kernel = np.array([abs(m["bf16 CUDA (cf5b)"][g] - m["bf16 MPS (Mac)"][g]) for g in genes])

    # ---- z 換算のための分散分解（対照 1,000 個、3 seed）
    zs = [pd.read_csv(E3 / f"results/gene_seed_008-{s}.csv") for s in (81, 82, 83)]
    t = zs[0][["gene"]]
    for s, df in zip((81, 82, 83), zs):
        t = t.merge(df[["gene", "z"]].rename(columns={"z": f"z{s}"}), on="gene")
    g005 = json.loads((ROOT / "experiments/EXP-20261001-001/work/genes_005.json").read_text())["genes"]
    gx = json.loads((ROOT / "experiments/EXP-20261002-002/work/extra_null_007.json").read_text())["追加する対照"]
    role = {r["ensembl"]: r["role"] for r in g005 + gx}
    t = t[t.gene.map(role).isin(["negative", "negative_extra"])]
    Z = t[["z81", "z82", "z83"]].to_numpy()
    v_within = float(Z.var(axis=1, ddof=1).mean())          # seed 内（細胞の引き直し＋その seed の bf16 の揺れ）
    v_between = float(Z.mean(axis=1).var(ddof=1))
    v_gene = v_between - v_within / 3                       # 遺伝子ごとに再現する成分
    sd_seed_z = v_within ** 0.5                             # 単一 seed の「引き直し」の SD（z 単位）

    # 平均 |差| → SD の換算。|X−Y|（独立・等分散）の期待値 = sqrt(2/pi)*sqrt(2)*sd = 1.1284*sd
    k = (2 / np.pi) ** 0.5 * (2 ** 0.5)
    sd_seed_raw = d_seed.mean() / k                         # 生の単位での「引き直し」1 回ぶんの SD
    scale = sd_seed_z / sd_seed_raw                         # 生の単位 → z
    sd_dtype_z = d_dtype.mean() / (2 / np.pi) ** 0.5 * scale   # fp32 を真値とみなした bf16 の 1 回ぶんの SD

    out = {
        "目的": "bf16 で合否を論じられる z の下限を、実測から決める",
        "注": "摂動は 1 回も回していない。EXP-20261003-001 の機体・数値型の確認で保存された埋め込みを読んだだけ",
        "突き合わせた遺伝子": len(genes),
        "生の単位での平均|差|（主指標の中央値）": {
            "bf16 ↔ fp32（同じ細胞・同じ Mac）": float(d_dtype.mean()),
            "bf16 ↔ bf16（CUDA ↔ MPS、同じ細胞）": float(d_kernel.mean()),
            "seed 違い（同じ数値型・同じ機体）": float(d_seed.mean()),
        },
        "比": {
            "数値型 ÷ seed 違い": float(d_dtype.mean() / d_seed.mean()),
            "カーネル ÷ seed 違い": float(d_kernel.mean() / d_seed.mean()),
        },
        "対照 1,000 個の z の分散分解": {
            "単一 seed の z の分散": v_gene + v_within,
            "seed 内（細胞の引き直し）": v_within,
            "遺伝子ごとに再現する成分": v_gene,
            "再現する割合": v_gene / (v_gene + v_within),
            "出どころ": "experiments/EXP-20261003-001/results/gene_seed_008-{81,82,83}.csv の対照 1,000 個",
        },
        "z 換算": {
            "細胞の引き直し 1 回ぶんの SD": sd_seed_z,
            "bf16 の遺伝子ごとの誤差 SD": float(sd_dtype_z),
            "bf16 の誤差 平均|差|": float(d_dtype.mean() * scale),
            "bf16 の誤差 観測最大|差|": float(d_dtype.max() * scale),
        },
        "判定保留の帯": {
            "下端": 2.73, "上端": 3.6, "幅": 0.87,
            "観測最大の bf16 誤差に対する比": 0.87 / float(d_dtype.max() * scale),
            "読み方": "上端 3.6 は、観測された bf16 の最大の食い違いより広い。"
                      "3 seed の平均 z がこの帯に入った遺伝子は bf16 では合否を論じない",
        },
        "seed で平均されるか（**確かめられていない**）": {
            "状況": "fp32 で回したのは seed81 の 12 遺伝子だけなので、bf16 の誤差が seed ごとに"
                    "独立か、遺伝子ごとに再現するかを**この実験の前に判別する材料がない**",
            "場合 A（seed ごとに独立）": "3 seed の平均では SD 0.243/√3 = 0.140 z まで縮む。"
                                         "ただし V-1 は**各 seed が 2.73 を超えること**を求めるので、"
                                         "1 seed ぶんの 0.243 z はそのまま効く",
            "場合 B（遺伝子ごとに再現）": "3 seed を重ねても 0.243 z のまま縮まない",
            "どちらでも": "2.73 のすぐ上では合否が bf16 の丸めで入れ替わりうる。"
                          "→ 判定保留の帯を置く",
            "やってはいけない推論": "seed 内の分散（0.0799 = SD 0.283）に bf16 の 0.243 が"
                                     "収まるからといって「独立だ」とは言えない。収まることは"
                                     "場合 A と整合するだけで、場合 B を否定しない",
        },
        "自己整合な bf16 の実行での効き方": "候補も対照も同じ bf16 で測るので、標準化そのものは"
                                            "内部で整合する。効くのは (1) 対照の散らばりに bf16 の"
                                            "成分が乗って z がわずかに薄まること、(2) 個々の遺伝子が"
                                            "±0.24 z 動くこと、の 2 つ",
        "計算量（実測）": {
            "GB10 の 1 回の forward（batch 8 × 2048 トークン）": {"bf16 秒": 0.220, "fp32 秒": 1.107,
                                                                 "fp32 ÷ bf16": 5.02,
                                                                 "測り方": "transformers の AutoModel を"
                                                                           "spark-cf5b で直接回した（2026-10-04）。"
                                                                           "fp32 は tf32 を閉じた"},
            "Mac MPS の ISP 1 遺伝子（200 細胞）": {"bf16 秒": 112.0, "fp32 秒": 303.3, "fp32 ÷ bf16": 2.71,
                                                   "出どころ": "experiments/EXP-20261003-001/"
                                                               "results-machinecheck-{mac,fp32}/timing.jsonl"},
        },
    }
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()

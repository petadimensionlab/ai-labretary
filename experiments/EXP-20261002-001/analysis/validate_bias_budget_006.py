#!/usr/bin/env python3
"""新しい健全性の検査（歪みの予算）を、EXP-20261001-001 の**実データ**で点検する。

**主判定（介入セット 119 個）の値は 1 回も計算しない。** 使うのは対照 500 個と、
対照だけから作ったセットと、目標の軸の陽性対照 71 個だけである。
介入セットの遺伝子は読み込んだ直後に落とす（下の `drop_primary`）。

走る検査:

  A 無作為な対照 119 個のセットは、**登録済みの指標**で合格するか（偶然落ちないか）
  B 対照だけで作った**偏った塊**で、見積もり B が**実際の歪み**と一致するか
    （塊の中身はすべて対照なので、観測されたセット平均 z がそのまま実際の歪みである）
  C 落とす検査（M5 など）が、同じ塊でどう振る舞うか（比較用）
  D （`--cache` を渡したときだけ）**きれいな指標 = 無作為な向きへの射影 60 本**での合格率

A〜C は `results/eval/gene_seed_005.csv`（事前登録済みの評価が残した遺伝子ごとの表）だけで
走る。D は事後診断の控え（`posthoc-m5/build_cache.py` が作る `cells2.parquet` と
`disp_kept.npy`）が要る。控えは一時領域にあるので、無ければ D は飛ばす。
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "experiments/EXP-20261001-001"
sys.path.insert(0, str(ROOT / "scripts/analysis"))
from bias_budget import (BUDGET, CANONICAL_COVARIATES, DECISION_LINE,  # noqa: E402
                         bias_budget, control_quantile)

SEEDS = (81, 82, 83)
K = 119
#: 主判定とその入れ子（**この点検では一切使わない**）
PRIMARY_FAMILY = ("INT", "INT_S1", "INT_S2", "INT_S3", "INT_pathway_specific",
                  "INT_plus_T1T2", "T1", "T2")
PS = (0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9)


# --------------------------------------------------------------------------- 読み込み
def load_tables() -> tuple[dict, list, list, set]:
    """seed -> 遺伝子ごとの表（z と交絡の候補）、対照 500 個、陽性対照 71 個を返す。"""
    gj = json.loads((SRC / "work/genes_005.json").read_text())["genes"]
    ens = {g["symbol"]: g["ensembl"] for g in gj}
    subs = json.loads((SRC / "work/subsets_005.json").read_text())["subsets"]
    nulls = [ens[x["symbol"]] for x in
             json.loads((SRC / "work/null_pool_005.json").read_text())["versions"]["B_500"]]
    pc71 = [ens[x["symbol"]] for x in subs["positive_control71"]]
    primary = {ens[x["symbol"]] for k in PRIMARY_FAMILY for x in subs[k]}

    g = pd.read_csv(SRC / "results/eval/gene_seed_005.csv")
    g = g[~g.gene.isin(primary)].reset_index(drop=True)      # **主判定は落とす**
    assert not (set(g.gene) & primary), "主判定の遺伝子が残っている"
    out = {}
    for s in SEEDS:
        a = g[(g.seed == s) & (g.col == "m_goal")].set_index("gene")
        u = g[(g.seed == s) & (g.col == "m_uniform")].set_index("gene").D_raw
        n = a.loc[[x for x in nulls if x in a.index]]
        t = pd.DataFrame({
            "gene": a.index, "D_raw": a.D_raw.to_numpy(),
            "z": ((a.D_raw - n.D_raw.mean()) / n.D_raw.std(ddof=1)).to_numpy(),
            "dnorm_med": a.dnorm_med.to_numpy(), "lp_med": a.lp_med.to_numpy(),
            "ll_med": a.ll_med.to_numpy(), "relpos_med": a.relpos_med.to_numpy(),
            "uniform": u.reindex(a.index).to_numpy()})
        out[s] = t.reset_index(drop=True)
    nulls = [x for x in nulls if x in set(out[SEEDS[0]].gene)]
    pc71 = [x for x in pc71 if x in set(out[SEEDS[0]].gene)]
    return out, nulls, pc71, primary


def skewed_block(tab: pd.DataFrame, nulls, cov: str, p: float, k: int = K) -> list:
    """対照 500 個を性質 `cov` の順で並べ、平均分位がおよそ p になる k 個の塊を取る。"""
    n = tab[tab.gene.isin(set(nulls))].sort_values(cov).reset_index(drop=True)
    st = max(0, min(len(n) - k, int(round(p * len(n) - k / 2))))
    return list(n.gene[st:st + k])


# --------------------------------------------------------------------------- A
def check_A(tabs, nulls, n_draw=2000, rng=None) -> dict:
    """無作為な対照 119 個のセットは、登録済みの指標で合格するか。"""
    rng = rng or np.random.default_rng(606)
    arr = np.array(sorted(nulls))
    rows, ups, dq = [], [], []
    for _ in range(n_draw):
        gs = list(rng.choice(arr, K, replace=False))
        r = [bias_budget(tabs[s], gs, nulls) for s in SEEDS]
        rows.append(all(x["pass"] for x in r))
        ups.append(max(x["upper_joint"] for x in r))
        dq.append(np.mean([abs(d["分位のずれ"]) for d in r[0]["covariates"].values()]))
    return {"引いたセットの数": n_draw,
            "3 seed すべて合格した割合": round(float(np.mean(rows)), 4),
            "上側の最大の中央値": round(float(np.median(ups)), 4),
            "上側の最大の 95 パーセンタイル": round(float(np.quantile(ups, 0.95)), 4),
            "分位のずれの平均絶対値の中央値": round(float(np.median(dq)), 4)}


# --------------------------------------------------------------------------- B, C
def check_BC(tabs, nulls) -> dict:
    rows = []
    for cov in CANONICAL_COVARIATES:
        for p in PS:
            gs = skewed_block(tabs[SEEDS[0]], nulls, cov, p)
            est, act, upp, m5 = [], [], [], []
            for s in SEEDS:
                r = bias_budget(tabs[s], gs, nulls)
                est.append(r["B_joint"])
                upp.append(r["upper_joint"])
                t = tabs[s].set_index("gene")
                act.append(float(t.z.reindex(gs).mean()))
                n = t.loc[[x for x in nulls]]
                m5.append(float(stats.spearmanr(n.dnorm_med, n.D_raw).statistic))
            rows.append({"交絡の候補": cov, "露出 p": p,
                         "実際の歪み（セット平均 z、3 seed 平均）": round(float(np.mean(act)), 4),
                         "見積もり B（3 seed 平均）": round(float(np.mean(est)), 4),
                         "食い違い": round(float(np.mean(est) - np.mean(act)), 4),
                         "上側の最大": round(float(np.max(upp)), 4),
                         "新しい検査の合否": bool(np.max(upp) <= BUDGET),
                         "参考: 落とす検査 M5（セットに依らない）": round(float(np.mean(m5)), 4)})
    a = np.array([r["実際の歪み（セット平均 z、3 seed 平均）"] for r in rows])
    e = np.array([r["見積もり B（3 seed 平均）"] for r in rows])
    sl, ic = np.polyfit(a, e, 1)
    return {"塊ごとの比較": rows,
            "較正（実際 → 見積もり の最小二乗）": {
                "傾き（1.0 が理想）": round(float(sl), 3),
                "切片": round(float(ic), 4),
                "相関": round(float(np.corrcoef(a, e)[0, 1]), 3),
                "食い違いの絶対値の最大": round(float(np.max(np.abs(e - a))), 4),
                "食い違いの標準偏差": round(float(np.std(e - a)), 4)}}


def check_E(tabs, nulls, n_split=20, n_set=160, k=100, rng=None) -> dict:
    """**取り置きでの較正**（この点検の要）。

    対照 500 個を半分に分け、**片方（250 個）だけで傾きを当て**、もう片方から作った
    セットの歪みを予測して、**そのセットで実際に観測されたセット平均 z** と比べる。
    セットは「無作為」と「ある性質に偏らせた塊」を半々で作る。

    実際の値には、どうやっても消えない雑音（対照 k 個のセット平均 z の揺れ）が乗る。
    そのため「実際 ≒ 傾き × 見積もり」の傾きは**見積もり自身の雑音で希釈される**。
    希釈の分（1 − 平均 se_B² / 見積もりの分散）で割り戻した値が 1.0 に近いかを見る。
    あわせて、報告している信頼区間が本当に当たるか（被覆率）も測る。
    """
    rng = rng or np.random.default_rng(611)
    arr = np.array(sorted(nulls))
    rows = []
    for isp in range(n_split):
        perm = rng.permutation(arr)
        A, Bh = list(perm[:len(arr) // 2]), list(perm[len(arr) // 2:])
        for _ in range(n_set):
            if rng.integers(0, 2) == 0:
                S = list(rng.choice(Bh, k, replace=False))
            else:
                cov = CANONICAL_COVARIATES[rng.integers(0, len(CANONICAL_COVARIATES))]
                p = float(rng.uniform(0.15, 0.85))
                t = tabs[SEEDS[0]].set_index("gene").loc[Bh].sort_values(cov)
                st = max(0, min(len(Bh) - k, int(round(p * len(Bh) - k / 2))))
                S = list(t.index[st:st + k])
            keep = set(A) | set(S)
            for s in SEEDS:
                r = bias_budget(tabs[s][tabs[s].gene.isin(keep)], S, A)
                rows.append((isp, r["B_joint"],
                             float(tabs[s].set_index("gene").z.reindex(S).mean()),
                             r["se_B"]))
    d = pd.DataFrame(rows, columns=["split", "est", "act", "se"])
    sl = float(np.polyfit(d.est, d.act, 1)[0])
    att = float(1 - (d.se ** 2).mean() / d.est.var())
    nf = float(np.mean([np.std([float(tabs[s].set_index("gene").z
                                     .reindex(rng.choice(arr, k, replace=False)).mean())
                                for _ in range(2000)]) for s in SEEDS]))
    tot = np.sqrt(d.se ** 2 + nf ** 2)
    # **分け方（split）ごとに傾きを出し、そのばらつきで区間を作る。**
    # 行を無作為に取り直す素朴なブートストラップは、同じ分け方の行どうしが
    # 相関しているため区間を狭く見せる（実際、分け方を変えると傾きが 0.50〜0.71 と動いた）。
    per = np.array([float(np.polyfit(g.est, g.act, 1)[0]) for _, g in d.groupby("split")])
    sd_sl = float(per.std(ddof=1) / np.sqrt(len(per)))
    return {"比べたセットの数": int(len(d) / len(SEEDS)), "1 セットの遺伝子数": k,
            "傾きを当てた対照の数": len(arr) // 2,
            "観測の傾き（実際 ← 見積もり）": round(sl, 3),
            "分け方ごとの傾き": [round(float(x), 3) for x in per],
            "観測の傾きの 95% 区間（分け方のばらつきから）":
                [round(sl - 1.96 * sd_sl, 3), round(sl + 1.96 * sd_sl, 3)],
            "見積もり自身の雑音による希釈": round(att, 3),
            "希釈を戻した傾き（1.0 が理想）": round(sl / att, 3),
            "希釈を戻した傾きの 95% 区間": [round((sl - 1.96 * sd_sl) / att, 3),
                                            round((sl + 1.96 * sd_sl) / att, 3)],
            "se_B の平均": round(float(d.se.mean()), 4),
            "見積もりの標準偏差": round(float(d.est.std()), 4),
            f"対照 {k} 個のセット平均 z の雑音": round(nf, 4),
            "信頼区間の被覆率（目標 0.95）":
                round(float(np.mean(np.abs(d.act - d.est) <= 1.959964 * tot)), 3),
            "系統的なずれ 平均(実際 − 見積もり)": round(float((d.act - d.est).mean()), 4)}


def noise_floor(tabs, nulls, n_draw=4000, rng=None) -> dict:
    """対照 119 個のセット平均 z の揺れ（これ以下は実際の歪みと区別できない）。"""
    rng = rng or np.random.default_rng(607)
    arr = np.array(sorted(nulls))
    sds = {}
    for s in SEEDS:
        z = tabs[s].set_index("gene").z
        v = np.array([float(z.reindex(rng.choice(arr, K, replace=False)).mean())
                      for _ in range(n_draw)])
        sds[s] = round(float(v.std(ddof=1)), 4)
    return {"対照 119 個のセット平均 z の標準偏差": sds}


# --------------------------------------------------------------------------- D
def check_D(tabs, nulls, pc71, cache: Path, n_dir=60, rng=None) -> dict:
    """きれいな指標 = 無作為な向きへの射影 60 本での合格率（事後診断の控えが要る）。"""
    rng = rng or np.random.default_rng(2026)
    cells = pd.read_parquet(cache / "cells2.parquet")
    disp = np.load(cache / "disp_kept.npy", mmap_mode="r")
    elig = {int(k): v for k, v in json.loads((cache / "elig.json").read_text()).items()}
    sys.path.insert(0, str(SRC / "posthoc-m5"))
    import diag  # noqa: E402

    D = disp.shape[1]
    W = rng.normal(size=(D, n_dir))
    W /= np.linalg.norm(W, axis=0, keepdims=True)
    t0 = time.time()
    Y = np.empty((len(cells), n_dir), dtype=np.float64)
    step = 40000
    for i in range(0, len(cells), step):
        Y[i:i + step] = (np.asarray(disp[i:i + step], dtype=np.float64) @ W)
    Y /= cells.dnorm.to_numpy()[:, None]
    print(f"  無作為な向きへの射影 {n_dir} 本を計算 {time.time() - t0:.0f}s", flush=True)

    # --- 遺伝子ごとの残差を 60 本まとめて解く（設計行列は向きに依らないので使い回せる）
    nullset = set(nulls)
    ztabs = {s: {} for s in SEEDS}
    for s in SEEDS:
        cs = cells[cells.seed == s]
        ix = cs.index.to_numpy()
        nl = {h: hd for h, hd in cs[cs.gene.isin(nullset)].groupby("gene")}
        Xn = {h: diag.design([hd.lp.to_numpy(), hd.ll.to_numpy()]) for h, hd in nl.items()}
        G = {h: Xn[h].T @ Xn[h] for h in nl}
        B = {h: Xn[h].T @ Y[hd.index.to_numpy()] for h, hd in nl.items()}
        G_all, B_all = sum(G.values()), sum(B.values())
        rows, genes = [], []
        for gene, tg in cs.groupby("gene"):
            Gm = G_all - G[gene] if gene in G else G_all
            Bm = B_all - B[gene] if gene in B else B_all
            beta = np.linalg.solve(Gm, Bm)                       # (15, n_dir)
            Xg = diag.design([tg.lp.to_numpy(), tg.ll.to_numpy()])
            rows.append(Y[tg.index.to_numpy()].mean(axis=0) - Xg.mean(axis=0) @ beta)
            genes.append(gene)
        Draw = np.vstack(rows)                                   # (n_genes, n_dir)
        gi = pd.Index(genes)
        med = cs.groupby("gene").agg(dnorm_med=("dnorm", "median"), lp_med=("lp", "median"),
                                     ll_med=("ll", "median"), relpos_med=("relpos", "median"),
                                     uniform=("m_uniform", "mean")).reindex(gi)
        mask = np.isin(np.asarray(genes), np.array(sorted(nullset)))
        for j in range(n_dir):
            d = Draw[:, j]
            ztabs[s][j] = pd.DataFrame({
                "gene": gi, "z": (d - d[mask].mean()) / d[mask].std(ddof=1),
                **{c: med[c].to_numpy() for c in
                   ("dnorm_med", "lp_med", "ll_med", "relpos_med", "uniform")}})
        print(f"  seed{s} 完了 {time.time() - t0:.0f}s", flush=True)

    rg = np.random.default_rng(608)
    arr = np.array(sorted(nulls))
    rnd_sets = [list(rg.choice(arr, K, replace=False)) for _ in range(20)]
    res = {"無作為な対照 119 個（20 組 × 60 本）": [], "目標の軸の陽性対照 71 個（60 本）": []}
    for j in range(n_dir):
        for gs in rnd_sets:
            res["無作為な対照 119 個（20 組 × 60 本）"].append(
                all(bias_budget(ztabs[s][j], gs, nulls)["pass"] for s in SEEDS))
        res["目標の軸の陽性対照 71 個（60 本）"].append(
            all(bias_budget(ztabs[s][j], pc71, nulls)["pass"] for s in SEEDS))
    return {k: {"件数": len(v), "3 seed すべて合格した割合": round(float(np.mean(v)), 4)}
            for k, v in res.items()}


# --------------------------------------------------------------------------- 主
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default=None, help="事後診断の控えの場所（D を走らせるとき）")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parents[1]
                                         / "work/bias_budget_validation_006.json"))
    a = ap.parse_args()

    tabs, nulls, pc71, primary = load_tables()
    out = {"データ": {"出所": "experiments/EXP-20261001-001/results/eval/gene_seed_005.csv",
                      "落とした主判定の遺伝子数": len(primary),
                      "残った遺伝子数": int(len(tabs[SEEDS[0]])),
                      "対照 500 個のうち使えたもの": len(nulls),
                      "目標の軸の陽性対照": len(pc71),
                      "seed": list(SEEDS)},
           "設定": {"交絡の候補": list(CANONICAL_COVARIATES), "予算": BUDGET,
                    "判断の線": DECISION_LINE, "セットの大きさ K": K}}

    out["雑音の下限"] = noise_floor(tabs, nulls)
    print("雑音の下限 " + json.dumps(out["雑音の下限"], ensure_ascii=False), flush=True)
    out["A 無作為な対照セットの合格率（登録済みの指標）"] = check_A(tabs, nulls)
    print("A " + json.dumps(out["A 無作為な対照セットの合格率（登録済みの指標）"],
                            ensure_ascii=False), flush=True)
    out["E 取り置きでの較正（この点検の要）"] = check_E(tabs, nulls)
    print("E " + json.dumps(out["E 取り置きでの較正（この点検の要）"], ensure_ascii=False),
          flush=True)
    out["BC 偏った塊での見積もりと実際（記述。1 塊の実際の値には雑音が乗る）"] = check_BC(tabs, nulls)

    # --- 実在のセットの偏り（記述。陽性対照は判定に使う前の参考）
    q = {}
    for cov in CANONICAL_COVARIATES:
        t = tabs[SEEDS[0]].set_index("gene")
        cv = control_quantile(t[cov].reindex(pc71).to_numpy(),
                              t[cov].reindex(nulls).to_numpy())
        q[cov] = round(float(cv.mean()), 4)
    out["目標の軸の陽性対照 71 個の平均分位"] = q
    out["目標の軸の陽性対照 71 個の歪み"] = {
        str(s): bias_budget(tabs[s], pc71, nulls) for s in SEEDS}

    if a.cache:
        out["D きれいな指標（無作為な向き 60 本）の合格率"] = check_D(
            tabs, nulls, pc71, Path(a.cache))
        print("D " + json.dumps(out["D きれいな指標（無作為な向き 60 本）の合格率"],
                                ensure_ascii=False), flush=True)
    else:
        out["D きれいな指標（無作為な向き 60 本）の合格率"] = "控えが無いので飛ばした"

    Path(a.out).write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    print(f"書き出し: {a.out}")


if __name__ == "__main__":
    main()

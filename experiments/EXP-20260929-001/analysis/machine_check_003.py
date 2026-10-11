"""EXP-20260929-001 — 2 台の DGX が同じ答えを出すことを、本番の前に確かめる（PI 指示 2026-09-28）。

2 台はパッケージの版が違う:
  spark-c083  torch 2.13 / transformers 4.46.3   （seed 62・53 を担当）
  spark-cf5b  torch 2.12 / transformers 4.57     （seed 61 を担当）

**seed と機体が 1 対 1 に対応しているので、機体差は「3 seed すべてで合格」という条件と
見分けられない。** したがって機体の一致は、本番の前に確かめなければならない。

手順:
  1) `--make` で共通の小さな設定 work/config-machinecheck.yaml を作る
     （seed 61 / 分泌 24 + 目標マーカー 7 + 筋線維芽細胞マーカー 24 + 恒常発現 3 + 中立 7 + 参照 2 + 対照 20 = 87 遺伝子 × 200 細胞 ≒ 0.31 時間）
  2) 両機体で同じ設定を回し、results-machinecheck-<host>/ に書く
  3) `--compare <dirA> <dirB>` で突き合わせ、work/machine_check_003.json に書く

判定基準（**結果を見る前に固定する**）:
  H-1（主）  87 遺伝子の標準偏差単位 z について、機体間の差の
             **平均（系統成分）の絶対値 <= 0.02**、かつ **最大絶対差 <= 0.05**
             根拠: セット平均に効くのは系統成分。0.02 は絶対値条件 0.245 の 8%、
             検出下限 0.44（K=120, 遺伝子間 SD 込み）の 4.5% にあたる。
             最大 0.05 でも 0.245 の 20% で、1 遺伝子の寄与は 0.05/120 = 0.0004 にすぎない。
  H-2       細胞ごとの主指標 m_goal の 2 台間 Spearman 順位相関 rho >= 0.99（約 17,400 細胞）
  H-3       状態重心ベクトル 4 本（FIB_orgHigh / FIB_orgLow / FIB_ecmHigh / MYOF_fib）の
             2 台間コサイン >= 0.9999
  H-4       記述値（合否に入れない）: 埋め込みの最大絶対差、移動量 ||Δ|| の中央値の比

不一致だったときの扱い（**これも結果を見る前に決める**）:
  (a) まず transformers を 4.46.3 にそろえて H-1〜H-3 を取り直す
  (b) そろえられないなら **spark-c083 1 台で 3 seed すべてを回す**（21.9 時間。上限 24 時間内で G2 も不要）
  (c) 「機体を要因として記録するだけ」は採らない。seed と機体が交絡していて分離できないため

実行: .venv/bin/python experiments/EXP-20260929-001/analysis/machine_check_003.py --make
      .venv/bin/python experiments/EXP-20260929-001/analysis/machine_check_003.py --compare A B
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260929-001"
SEED = 61      # 機体の一致確認は seed 61 の細胞で行う
CRIT = {"H1_mean_abs_max": 0.02, "H1_max_abs_max": 0.05, "H2_rho_min": 0.99, "H3_cos_min": 0.9999}


def make() -> None:
    """設定は analysis/make_configs_003.py が work/config-machinecheck.yaml に生成済み。
    ここでは中身を表示するだけにして、**設定の生成箇所を 1 つに保つ**。"""
    p = E / "work/config-machinecheck.yaml"
    cfg = yaml.safe_load(p.read_text())
    n = cfg["compute"]["n_gene_runs"]
    print(f"{p.relative_to(ROOT)}: {n} 遺伝子 × 200 細胞 × 1 seed"
          f" ≒ {cfg['compute']['gpu_hours_estimate']} 時間/機体")


def compare(a: Path, b: Path) -> None:
    import importlib.util
    spec = importlib.util.spec_from_file_location("ev", Path(__file__).resolve().parent / "evaluate_003.py")
    ev = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ev)
    cfg = yaml.safe_load((E / "work/config-machinecheck.yaml").read_text())
    gj = json.loads((E / "work/genes_003.json").read_text())["genes"]
    sym = {g["ensembl"]: g["symbol"] for g in gj}
    nulls = set(json.loads((E / "results-machinecheck-A.negatives.json").read_text())) \
        if (E / "results-machinecheck-A.negatives.json").exists() else None

    out: dict = {"criteria": CRIT, "dirs": [str(a), str(b)]}
    tabs, cells = {}, {}
    for tag, d in (("A", a), ("B", b)):
        cfg2 = json.loads(json.dumps(cfg))
        cfg2["output"]["results_dir"] = d.name
        c = ev.load(cfg2)
        ev.metric_columns(c)
        cells[tag] = c
        ns = nulls or set(json.loads((d / f"run-{SEED}-L-1/negatives.json").read_text()))
        nl = {h: hd for h, hd in c[c.gene.isin(ns)].groupby("gene")}
        el = ev.eligibility(c, nl, 0.5, 0.8)
        t = ev.gene_stats(c, nl, el, "m_goal", 1)
        mu, sd = t[t.gene.isin(ns)].D_raw.mean(), t[t.gene.isin(ns)].D_raw.std(ddof=1)
        t["z"] = (t.D_raw - mu) / sd
        tabs[tag] = t.set_index("gene")

    g = sorted(set(tabs["A"].index) & set(tabs["B"].index))
    dz = (tabs["A"].loc[g, "z"] - tabs["B"].loc[g, "z"]).to_numpy()
    out["H1"] = {"n_genes": len(g), "mean_abs": float(abs(dz.mean())), "max_abs": float(np.abs(dz).max()),
                 "worst_genes": [sym.get(x, x) for x in np.array(g)[np.argsort(-np.abs(dz))[:5]]],
                 "pass": bool(abs(dz.mean()) <= CRIT["H1_mean_abs_max"]
                              and np.abs(dz).max() <= CRIT["H1_max_abs_max"])}

    ka = cells["A"].set_index(["gene", "cell_uid"]).m_goal
    kb = cells["B"].set_index(["gene", "cell_uid"]).m_goal
    idx = ka.index.intersection(kb.index)
    rho = float(stats.spearmanr(ka.loc[idx], kb.loc[idx]).statistic)
    out["H2"] = {"n_cells": int(len(idx)), "spearman_rho": rho, "pass": bool(rho >= CRIT["H2_rho_min"])}

    cs = {}
    for tag, d in (("A", a), ("B", b)):
        sz = np.load(d / f"run-{SEED}-L-1/state_emb_vectors.npz")
        xz = np.load(E / f"work/extra_centroids-seed{SEED}.npz")
        cs[tag] = {k: sz[k].astype(np.float64) for k in sz.files} | {"MYOF_fib": xz["MYOF_fib"].astype(np.float64)}
    coss = {k: float(cs["A"][k] @ cs["B"][k] / (np.linalg.norm(cs["A"][k]) * np.linalg.norm(cs["B"][k])))
            for k in cs["A"] if k in cs["B"]}
    out["H3"] = {"cosine": coss, "min": min(coss.values()),
                 "pass": bool(min(coss.values()) >= CRIT["H3_cos_min"])}

    dn = {tag: float(cells[tag].dnorm.median()) for tag in ("A", "B")}
    out["H4_descriptive"] = {"移動量 ||Δ|| の中央値": dn, "比": round(dn["A"] / max(dn["B"], 1e-12), 4)}

    out["pass"] = bool(out["H1"]["pass"] and out["H2"]["pass"] and out["H3"]["pass"])
    out["不一致のときの扱い"] = ("(a) transformers を 4.46.3 にそろえて取り直す → "
                                 "(b) そろえられないなら spark-c083 1 台で 3 seed すべてを回す（11.0 時間）。"
                                 "(c)「機体を要因として記録するだけ」は採らない（seed と交絡して分離できない）")
    (E / "work/machine_check_003.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--make", action="store_true")
    ap.add_argument("--compare", nargs=2, metavar=("DIR_A", "DIR_B"))
    args = ap.parse_args()
    if args.make:
        make()
    elif args.compare:
        compare(Path(args.compare[0]), Path(args.compare[1]))
    else:
        ap.error("--make か --compare のどちらかを指定してください")

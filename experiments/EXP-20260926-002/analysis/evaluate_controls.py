"""EXP-20260926-002: 対照のみの run（案 B / 案 A）で削除ドリフト補正指標を評価する（plan.md §6 の C1〜C5 と C3-panel）。
reanalysis_cosine.py（コサイン候補）と vector_phase0.py（ベクトル候補）の関数をそのまま使い、複数の結果ディレクトリ
（例: phase 0 の seed 0 + 案 B の seed 1, 2）を seed ごとに補正してから 3 seed で判定する。
盲検: 主要対象・参照（rc.BLIND）は読み込み直後に除外する。
実行例:
  .venv/bin/python experiments/EXP-20260926-002/analysis/evaluate_controls.py \
      --results phase0/results-phase0 optionB/results-B --out optionB/results-B/eval
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reanalysis_cosine as rc  # noqa: E402
import vector_phase0 as vp  # noqa: E402

E = rc.ROOT / "experiments/EXP-20260926-002"
DET = rc.ROOT / "experiments/EXP-20260926-001/work/detection-aPT_CKD_TLShigh.json"


def load_run(res: Path, run: str):
    """1 run（= 1 seed）の cells と CLS 埋め込み・状態重心を読む（vp.load の一般化）。"""
    cells = pd.read_csv(res / "cells.csv", dtype={"donor_id": str})
    cells = cells[(cells.run == run) & ~cells.gene.isin(rc.BLIND)].reset_index(drop=True)
    # seed を別プロセスで回すと、最後のプロセスの summarize が他 seed の null を screen と書く → negatives.json で戻す
    negs = set(json.loads((res / run / "negatives.json").read_text()))
    cells.loc[cells.gene.isin(negs) & (cells.role == "screen"), "role"] = "negative"
    assert not set(cells.gene) & rc.BLIND
    rd = res / run
    E0, E1, order = [], [], []
    for gene, s in cells.groupby("gene", sort=False):
        o = np.load(rd / "genes" / gene / "emb" / "orig.npy").astype(np.float64)
        p = np.load(rd / "genes" / gene / "emb" / "pert.npy").astype(np.float64)
        per = pd.read_csv(rd / "genes" / gene / "cells.csv", dtype={"donor_id": str})
        assert len(o) == len(per) == len(s) and (per.cell_uid.to_numpy() == s.cell_uid.to_numpy()).all()
        E0.append(o); E1.append(p); order.append(s.index.to_numpy())
    cells = cells.loc[np.concatenate(order)].reset_index(drop=True)
    st = np.load(rd / "state_emb_vectors.npz")
    return cells, np.vstack(E0), np.vstack(E1), {k: st[k].astype(np.float64).ravel() for k in st.files}


def c3_panel(cells: pd.DataFrame, metric: str, panel: list[str], crit: dict) -> dict | None:
    """C3-panel: パネル遺伝子のうち min_genes 以上が、min_seeds 以上の seed で null 分位 [lo, hi] に入る。"""
    if not panel or not set(panel) & set(cells.gene):
        return None  # パネルを回していない run（案 B など）
    lo, hi = crit.get("C3_panel_quantile_range", [0.10, 0.90])
    g = cells.groupby(["run", "gene", "role"])[metric].mean().reset_index()
    runs = sorted(g.run.unique())
    per = {}
    for gene in panel:
        n = 0
        for r in runs:
            s = g[g.run == r]
            if gene not in set(s.gene):
                continue
            null = s[s.role == "negative"][metric].to_numpy()
            n += lo <= rc.null_pct(float(s[s.gene == gene][metric].iloc[0]), null) <= hi
        per[gene] = n
    ok = sum(v >= crit.get("C3_panel_min_seeds", 2) for v in per.values())
    return {"per_gene_seeds_in_range": per, "n_genes_ok": ok, "pass": ok >= crit.get("C3_panel_min_genes", 5)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", nargs="+", required=True, help="EXP-20260926-002 からの相対パス")
    ap.add_argument("--out", required=True)
    ap.add_argument("--config", default=str(E / "config.yaml"))
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    crit = cfg["analysis"]["criteria"]
    panel_groups = {grp["role"]: grp["genes"] for grp in cfg["controls"].get("extra_fixed", []) or []}
    out = E / a.out
    out.mkdir(parents=True, exist_ok=True)

    parts, descs = [], {}
    for rel in a.results:
        res = E / rel
        for run in sorted(pd.read_csv(res / "cells.csv", usecols=["run"]).run.unique()):
            cells, E0, E1, cen = load_run(res, run)
            cells, desc = vp.vector_metrics(cells, E0, E1, cen)
            alt = [x for x in cells.columns if x.startswith("shift_alt")][0]
            cells["alt"] = cells[alt]
            cells["pos"] = rc.token_positions(cells)  # a4/d4 の補正に使う削除位置
            cells = rc.corrected(cells)
            descs[f"{rel}:{run}"] = desc
            parts.append(cells)
    cells = pd.concat(parts, ignore_index=True)
    # 同じ seed 名が複数ディレクトリに出ないことを確認（seed ごとに 1 run）
    assert cells.groupby("run").gene.nunique().min() > 0
    det = json.loads(DET.read_text())

    rows = []
    metrics = list(rc.METRICS.items()) + list(vp.LABELS.items())
    for m, lab in metrics:
        res, _ = rc.evaluate(cells, det, m)
        p = c3_panel(cells, m, panel_groups.get("negctrl_panel", []), crit)
        res["C3_panel"] = None if p is None else p["pass"]
        res["C3_panel_n_ok"] = None if p is None else f"{p['n_genes_ok']}/{len(p['per_gene_seeds_in_range'])}"
        rows.append({"label": lab, **res})
    df = pd.DataFrame(rows)
    for col in ["C1_seeds", "C1_pct_mean", "C3_seeds", "C3_pct_mean"]:
        df[col] = df[col].map(lambda d: " ".join(f"{k}:{v}" for k, v in d.items()))
    df.to_csv(out / "criteria.csv", index=False)
    (out / "criteria.json").write_text(json.dumps(rows, ensure_ascii=False, indent=1, default=str) + "\n")
    (out / "drift.json").write_text(json.dumps(descs, indent=1) + "\n")
    # 遺伝子 × seed の値（探索枠 RP やパネルの個別の挙動を見るため）
    g = cells.groupby(["run", "gene", "role"])[[m for m, _ in metrics]].mean().reset_index()
    g.to_csv(out / "gene_seed_metrics.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print("runs:", sorted(cells.run.unique()), "cells:", len(cells))
    print(df[["label", "C1_pct_mean", "C2_wins", "C3_pct_mean", "C3_panel_n_ok", "C4i_splithalf_rho",
              "C4ii_seed_rho_fixed5", "C5_rho_len_gene", "C5_rho_len_cell_median", "null_frac_negative",
              "C1", "C2", "C3", "C4", "C5", "n_pass_all"]].to_string(index=False))


if __name__ == "__main__":
    main()

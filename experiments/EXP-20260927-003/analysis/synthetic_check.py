"""EXP-20260927-003: 評価コード（evaluate_fib.py）を合成データで検査する。事前登録の一部。

実際の摂動結果は使わない。config.yaml の遺伝子構成と work/precheck_coverage.json の細胞数だけを使い、
細胞ごとの値を人工的に作って evaluate_fib.py を走らせる。目的は 2 つ。
  (1) 交絡だけ（元の位置と入力の長さに強く依存する見かけの効果）のとき「否定」になること
  (2) 主対象 2 つに本当の効果を入れたとき「支持」になること
あわせて較正（null 遺伝子自身の分位の裾）と交絡の除去（残差相関）を確認する。

実行: .venv/bin/python experiments/EXP-20260927-003/analysis/synthetic_check.py
出力: work/synthetic/{A_confound_only,B_true_effect}/eval/criteria.json と標準出力の要約
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260927-003"
CFG = yaml.safe_load(open(E / "config.yaml"))
SEEDS = CFG["robustness"]["seeds"]
DIM = 64


def gene_roles() -> dict[str, str]:
    ctl = CFG["controls"]
    r = {g: "target" for g in CFG["perturbation"]["genes_to_perturb"]}
    r |= {g: "positive" for g in ctl["positive"]}
    r |= {g: "housekeeping" for g in ctl["negative_fixed"]}
    r |= {g: "negative" for g in ctl["negative_sampler"]["genes"]}
    for grp in ctl["extra_fixed"]:
        r |= {g: grp["role"] for g in grp["genes"]}
    return r


def build(scenario: str, base: Path) -> Path:
    rng = np.random.default_rng(7)
    cov = json.loads((E / "work/precheck_coverage.json").read_text())
    stats = pd.read_csv(E / "work/gene_stats_fib.csv").set_index("ensembl")
    roles = gene_roles()
    res = base / scenario
    shutil.rmtree(res, ignore_errors=True)
    (res).mkdir(parents=True)
    sym = dict(zip(stats.index, stats.symbol))
    true_effect = {CFG["perturbation"]["genes_to_perturb"][i]: 3e-4 for i in (0, 1)} if scenario.startswith("B") else {}
    myof_effect = {g: 3e-4 for g in CFG["controls"]["positive"]}   # ACTA2/TAGLN は MYOF 重心へ近づく
    cen = rng.normal(size=DIM)
    cen /= np.linalg.norm(cen)
    rows = []
    for s in SEEDS:
        for g, role in roles.items():
            n = cov["cells_by_gene_seed"][f"{sym[g]}@{s}"]
            pos = np.clip(rng.normal(stats.loc[g, "pos_med"], stats.loc[g, "pos_med"] * 0.35, n), 1, None)
            ln = np.clip(rng.normal(stats.loc[g, "len_med"], 250, n), 300, 4095)
            lp, ll = np.log2(1 + pos), np.log2(ln)
            # 交絡: 元の位置と長さに強く依存する見かけの効果（すべての遺伝子に共通で入る）
            y = 2.0e-4 - 2.2e-5 * lp + 1.1e-5 * ll + 1.4e-6 * lp ** 2 + rng.normal(0, 8e-5, n)
            y += rng.normal(0, 2e-5)                      # 遺伝子ごとのばらつき
            y += true_effect.get(g, 0.0)
            ymy = 1.0e-4 - 1.0e-5 * lp + rng.normal(0, 8e-5, n) + rng.normal(0, 2e-5) + myof_effect.get(g, 0.0)
            gd = res / f"run-{s}-L-1" / "genes" / g / "emb"
            gd.mkdir(parents=True, exist_ok=True)
            o = rng.normal(size=(n, DIM))
            o /= np.linalg.norm(o, axis=1, keepdims=True)
            co = o @ cen
            # cos(pert, cen) = cos(orig, cen) + ymy になるよう、重心方向に少しだけ足して正規化しなおす
            p = o + (ymy / (1 - co ** 2 + 1e-9))[:, None] * (cen[None, :] - co[:, None] * o)
            p /= np.linalg.norm(p, axis=1, keepdims=True)
            np.save(gd / "orig.npy", o.astype(np.float32))
            np.save(gd / "pert.npy", p.astype(np.float32))
            rows.append(pd.DataFrame({
                "run": f"run-{s}-L-1", "gene": g, "cell_uid": np.arange(n),
                "donor_id": [f"D{i % 14}" for i in range(n)], "length": ln.astype(int),
                "pos": pos.astype(int), "shift_goal": y, "shift_start": -y, "role": role}))
    pd.concat(rows, ignore_index=True).to_csv(res / "cells.csv", index=False)
    for s in SEEDS:
        np.savez(base / f"extra_centroids-seed{s}.npz", **{CFG["analysis"]["state_control"]["centroid"]: cen})
    return res


def main() -> None:
    base = E / "work/synthetic"
    base.mkdir(parents=True, exist_ok=True)
    out = {}
    for scenario in ("A_confound_only", "B_true_effect"):
        res = build(scenario, base)
        cfg = json.loads(json.dumps(CFG))
        cfg["output"] = {"experiment_dir": str(base.relative_to(ROOT)), "results_dir": scenario}
        cfg["extra_centroids"]["output"] = str((base / "extra_centroids-seed{seed}.npz").relative_to(ROOT))
        cfg["data"]["tokenized_dataset"] = CFG["data"]["tokenized_dataset"]
        (base / f"config-{scenario}.yaml").write_text(yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False))
        # comparators.csv / gene_stats_fib.csv は本体の work/ を参照させる
        (base / "work").mkdir(exist_ok=True)
        for f in ("comparators.csv", "gene_stats_fib.csv"):
            tgt = base / "work" / f
            if not tgt.exists():
                shutil.copy(E / "work" / f, tgt)
        r = subprocess.run([sys.executable, str(E / "analysis/evaluate_fib.py"),
                            "--config", str(base / f"config-{scenario}.yaml")],
                           capture_output=True, text=True)
        print(scenario, r.stdout.strip() or r.stderr[-2000:])
        cj = res / "eval" / "criteria.json"
        if cj.exists():
            p = json.loads(cj.read_text())["primary_pos_matched_permutation"]
            out[scenario] = {"verdict": p["verdict_proposal"], "n_target_pass": p["n_target_pass"],
                             "C8a": p["C8a_calibration"], "C8b_rho": p["C8b_confounding"]["rho_by_seed"],
                             "C11": p["C11"], "C3": p["C3"], "C3_panel": p["C3_panel"]}
    (base / "summary.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""EXP-20260930-001 — 機体の一致確認用の設定を作る（plan.md §8.2。**事前登録の対象外**）。

48 遺伝子 × 200 細胞 × seed 81 = T1 11（全部）+ T2 7（全部）+ 陽性対照 10（先頭10）
+ 恒常発現 3（全部）+ 中立パネル 7（全部）+ 対照 10（先頭10）。

config.yaml と同じ生成ロジック（analysis/make_configs_004.py）を土台にして、
遺伝子リストと results_dir・seeds だけを差し替える。判定基準・補正・チェックポイント設定は同一。

出力: work/config-machinecheck.yaml
実行: .venv/bin/python experiments/EXP-20260930-001/analysis/make_configs_machinecheck.py
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import yaml

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20260930-001"
ED = "experiments/EXP-20260930-001"
SEED = 81


def main() -> None:
    cfg = yaml.safe_load((E / "config.yaml").read_text())
    g = json.loads((E / "work/genes_004.json").read_text())
    ens = {x["symbol"]: x["ensembl"] for x in g["genes"]}
    by = g["by_role"]

    mc_genes = {
        "target_T1": by["target_T1"],
        "target_T2": by["target_T2"],
        "positive_control": by["positive_control"][:10],
        "housekeeping": by["housekeeping"],
        "negctrl_panel": by["negctrl_panel"],
        "negative": by["negative"][:10],
    }
    n_total = sum(len(v) for v in mc_genes.values())
    assert n_total == 48, f"48 遺伝子にならない: {n_total}"

    t1t2 = [ens[s] for s in mc_genes["target_T1"] + mc_genes["target_T2"]]
    null10 = [ens[s] for s in mc_genes["negative"]]

    c2 = copy.deepcopy(cfg)
    c2["robustness"]["seeds"] = [SEED]
    c2["state_embs"]["rows_file"] = f"{ED}/work/state_rows-seed{{seed}}.json"
    c2["perturbation"]["genes_to_perturb"] = t1t2
    c2["controls"] = {
        "positive": [ens[s] for s in mc_genes["positive_control"]],
        "negative_fixed": [ens[s] for s in mc_genes["housekeeping"]],
        "fixed_role": "housekeeping",
        "extra_fixed": [{"role": "negctrl_panel", "genes": [ens[s] for s in mc_genes["negctrl_panel"]]}],
        "negative_random_n": len(null10),
        "negative_sampler": {"method": "explicit", "genes": null10},
        "comparator": None,
    }
    c2["output"]["results_dir"] = "results-machinecheck"   # compare 時に host ごとへ差し替える
    c2["compute"]["n_gene_runs"] = n_total
    c2["compute"]["n_cell_gene"] = n_total * 200 * 1
    c2["compute"]["gpu_hours_estimate"] = round(
        c2["compute"]["n_cell_gene"] * 0.1154 / 3600
        + 1 * (950 * 2 + 270 + 135) * 0.1154 / 2 / 3600, 3)

    head = ("# EXP-20260930-001 — 機体の一致確認用（**事前登録の対象外**。plan.md §8.2）\n"
            "# 48 遺伝子（T1 11 + T2 7 + 陽性対照 10 + 恒常発現 3 + 中立 7 + 対照 10）× 200 細胞 × seed 81\n"
            "# 作成: analysis/make_configs_machinecheck.py（手で書き換えない）\n")
    (E / "work/config-machinecheck.yaml").write_text(head + yaml.safe_dump(c2, allow_unicode=True, sort_keys=False))

    (E / "work/machinecheck_genes.json").write_text(json.dumps(
        {"n_total": n_total, "seed": SEED, "by_role": mc_genes}, ensure_ascii=False, indent=1) + "\n")

    print("48 遺伝子確認:", n_total, "延べ 細胞×遺伝子:", c2["compute"]["n_cell_gene"],
          "見積もり:", c2["compute"]["gpu_hours_estimate"], "時間/機体")


if __name__ == "__main__":
    main()

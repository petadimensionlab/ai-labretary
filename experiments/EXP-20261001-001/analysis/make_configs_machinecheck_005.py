#!/usr/bin/env python3
"""EXP-20261001-001 — 機体の一致確認用の設定を作る（plan.md §8.2。**事前登録の対象外**）。

44 遺伝子 × 200 細胞 × seed 81 = 介入 20（先頭20）+ 陽性対照 10（先頭10）
+ 採用した中立対照 4（ACTB・CERS4・ANAPC16・EEF1D）+ 対照 10（先頭10）。

config.yaml と同じ値を土台にして、遺伝子リストと results_dir・seeds だけを差し替える。
判定基準・補正・チェックポイント設定は同一。

出力: work/config-machinecheck.yaml, work/machinecheck_genes_005.json
実行: .venv/bin/python experiments/EXP-20261001-001/analysis/make_configs_machinecheck_005.py
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import yaml

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20261001-001"
ED = "experiments/EXP-20261001-001"
SEED = 81
NEUTRAL_ADOPTED = ["ACTB", "CERS4", "ANAPC16", "EEF1D"]  # plan.md §6.1.2 で採用された4個


def main() -> None:
    cfg = yaml.safe_load((E / "config.yaml").read_text())
    g = json.loads((E / "work/genes_005.json").read_text())
    ens = {x["symbol"]: x["ensembl"] for x in g["genes"]}
    by = g["by_role"]

    mc_genes = {
        "intervention": by["intervention"][:20],
        "positive_control": by["positive_control"][:10],
        "neutral_adopted": NEUTRAL_ADOPTED,
        "negative": by["negative"][:10],
    }
    n_total = sum(len(v) for v in mc_genes.values())
    assert n_total == 44, f"44 遺伝子にならない: {n_total}"
    for role, syms in mc_genes.items():
        missing = [s for s in syms if s not in ens]
        assert not missing, f"{role} に ensembl が無い遺伝子: {missing}"

    perturb_genes = [ens[s] for s in mc_genes["intervention"] + mc_genes["positive_control"]
                      + mc_genes["neutral_adopted"]]
    null10 = [ens[s] for s in mc_genes["negative"]]

    c2 = copy.deepcopy(cfg)
    c2["robustness"]["seeds"] = [SEED]
    c2["state_embs"]["rows_file"] = f"{ED}/work/state_rows-seed{{seed}}.json"
    c2["perturbation"]["genes_to_perturb"] = perturb_genes
    c2["controls"] = {
        "positive": [ens[s] for s in mc_genes["positive_control"]],
        "negative_fixed": [],
        "fixed_role": "neutral_adopted",
        "extra_fixed": [{"role": "neutral_adopted", "genes": [ens[s] for s in mc_genes["neutral_adopted"]]}],
        "negative_random_n": len(null10),
        "negative_sampler": {"method": "explicit", "genes": null10},
        "comparator": None,
    }
    c2["output"]["results_dir"] = "results-machinecheck"  # compare 時に host ごとへ差し替える
    c2["compute"]["n_gene_runs"] = n_total
    c2["compute"]["n_cell_gene"] = n_total * 200 * 1
    c2["compute"]["gpu_hours_estimate"] = round(
        c2["compute"]["n_cell_gene"] * 0.1154 / 3600
        + 1 * (950 * 2 + 270 + 135) * 0.1154 / 2 / 3600, 3)

    head = ("# EXP-20261001-001 — 機体の一致確認用（**事前登録の対象外**。plan.md §8.2）\n"
            "# 44 遺伝子（介入 20 + 陽性対照 10 + 採用した中立対照 4 + 対照 10）× 200 細胞 × seed 81\n"
            "# 作成: analysis/make_configs_machinecheck_005.py（手で書き換えない）\n")
    (E / "work/config-machinecheck.yaml").write_text(head + yaml.safe_dump(c2, allow_unicode=True, sort_keys=False))

    (E / "work/machinecheck_genes_005.json").write_text(json.dumps(
        {"n_total": n_total, "seed": SEED, "by_role": mc_genes}, ensure_ascii=False, indent=1) + "\n")

    print("44 遺伝子確認:", n_total, "延べ 細胞×遺伝子:", c2["compute"]["n_cell_gene"],
          "見積もり:", c2["compute"]["gpu_hours_estimate"], "時間/機体")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""EXP-20260930-001 — 設定ファイルを 1 か所から機械的に作る（遺伝子リストを人が写さない）。

出力:
  config.yaml                    本体（3 seed 分。**事前登録の対象**）
  config-seed{81,82,83}.yaml     seed 単位で計算機を分けるとき（seeds と results_dir 以外は同一）
  work/config-smoke.yaml         試し実行専用（**事前登録の対象外**。min_cells_per_gene を下げる）

**遺伝子で分けてはいけない。分けるのは seed 単位のみ**（機体差が結論と同じ形で出るため）。
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import yaml

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20260930-001"
ED = "experiments/EXP-20260930-001"
SEEDS = [81, 82, 83]

# 盲検中（**解除されていない 3 個**。LTBR・CDKN1A は PI が 2026-09-27 に解除済み）
BLOCKED = {"TP53": "ENSG00000141510", "SPP1": "ENSG00000118785", "LCN2": "ENSG00000148346"}

# セットごとの絶対値条件（**結果を見る前に固定**）。
# T1 / T2 / T1∪T2 / 陽性対照は HYP-20260930-01 §6 が書いた値をそのまま使う
# （§5.1 が「76 に増やしても 0.52」と書いた保守版の値を、実際に使える 71 個で 0.53 とした）。
# 入れ子の感度分析のセットは、同じ式（約 4.5 / √K）で小数第 2 位まで。
# **2026-09-30 の PI 判断で、主判定の陽性対照を 39 個版から 71 個版に替えた。**
SET_ABS_MIN = {"T1": 1.36, "T2": 1.70, "T1_T2": 1.06, "T1_no_CXCL12": 1.42,
               "positive_control71": 0.53,    # **主判定（健全性チェック）**
               "positive_control39": 0.69, "positive_control65": 0.56,
               "positive_control70": 0.54, "positive_control38": 0.73,
               "positive_control35": 0.76}


def main() -> None:
    g = json.loads((E / "work/genes_004.json").read_text())
    ens = {x["symbol"]: x["ensembl"] for x in g["genes"]}
    by = g["by_role"]
    null500 = [ens[s] for s in by["negative"]]
    t1t2 = [ens[s] for s in by["target_T1"] + by["target_T2"]]

    cfg = {
        "experiment_id": "EXP-20260930-001",
        "hypothesis_id": "HYP-20260930-01",
        "model": {"version": "V2", "directory": "models/Geneformer/Geneformer-V2-316M",
                  "type": "Pretrained", "num_classes": 0, "emb_layer": -1},
        "data": {
            "input_h5ad": "experiments/EXP-20260929-001/work/input.h5ad",   # **読み取りのみ**
            "tokenized_dataset": f"{ED}/work/tokenized.dataset",
            "cell_type_col": "subclass_l2",
            "state_key": "fib_anno_state",
            "custom_attr_name_dict": {"fib_anno_state": "fib_anno_state", "donor_id": "donor_id",
                                      "subclass_l2": "subclass_l2", "disease": "disease",
                                      "region": "region", "n_genes": "n_genes",
                                      "n_counts": "n_counts", "a1_row": "a1_row"},
            "filter_data": None, "max_ncells": None},
        "states": {"start_state": "FIB_MFIB", "goal_state": "aStr", "alt_states": ["dFIB", "dMFIB"]},
        "state_embs": {"max_cells_per_state": 950, "max_cells_per_donor": 60,
                       "rows_file": f"{ED}/work/state_rows-seed{{seed}}.json"},
        "cell_selection": {"donor_key": "donor_id", "max_cells_per_gene": 200,
                           "max_cells_per_donor": 20, "min_cells_per_gene": 120},
        "perturbation": {"perturb_type": "overexpress", "perturb_rank_shift": None,
                         "controls_only": False, "genes_to_perturb": t1t2, "combos": 0,
                         "anchor_gene": None, "emb_mode": "cls", "forward_batch_size": 8,
                         "nproc": 1},
        "blocked_genes": list(BLOCKED.values()),
        "controls": {
            "positive": [ens[s] for s in by["positive_control"]],   # **71 個（主判定）**
            "negative_fixed": [ens[s] for s in by["housekeeping"]],
            "fixed_role": "housekeeping",
            "extra_fixed": [
                {"role": "negctrl_panel", "genes": [ens[s] for s in by["negctrl_panel"]]},
                {"role": "housekeeping_extra", "genes": [ens[s] for s in by["housekeeping_extra"]]}],
            "negative_random_n": len(null500),
            "negative_sampler": {"method": "explicit", "genes": null500},
            "comparator": None},
        "capture": {"embeddings": True, "classifier_logits": False},
        "metric": {"confound_correction": {
            "variables": ["lp", "ll", "ld"], "form": "full_quadratic",
            "drop_zero_movement": True}},
        "analysis": {
            "primary_metric": "goal_minus_start_over_movement",
            "primary_unit": "gene_set",
            "primary_set": "T1",
            "primary_positive_control": "positive_control71",   # **PI 判断 2026-09-30**
            "primary_null_version": "B_500",
            "primary_test": "set_mean_vs_random_control_sets",
            "pos_definition": "index_in_input_ids",
            "subtract_centroids": ["start"],
            "descriptive_only_centroids": ["dFIB", "dMFIB"],
            "set_abs_min": SET_ABS_MIN,
            "secondary_null_versions": ["A_300"],
            "pos_matched_permutation": {"tolerance_log2": 0.5, "min_cover_frac": 0.8,
                                        "min_eligible_null_genes": 20},
            "criteria": {
                "require_all_seeds": True, "K_set": len(by["target_T1"]),
                "n_permutation": 200000, "set_percentile": 0.99, "set_draw_replace": True,
                "n_chance_draw": 200000, "max_chance_rate": 0.02,
                "neutral_range": [0.05, 0.95], "neutral_min_seeds": 2,
                "n_sham_sets": 200, "sham_match_bias": True, "sham_centre_tol": 0.10,
                "sham_spread_range": [0.7, 1.3], "sham_binomial_alpha": 0.05,
                "null_draw_seed": 20260930,
                "M1_tail_frac_max": 0.20, "M1_min_n": 60,
                "M2_abs_rho_max": 0.25, "M2_min_seeds": 2,
                "M3_range": [0.10, 0.90], "M3_min_seeds": 2, "M3_panel_min_genes": 5,
                "M4_abs_rho_uniform_max": 0.40,
                "M5_abs_rho_movement_max": 0.10}},
        "axis_feasibility": {
            "experiment_id": "EXP-20260930-001-precheck-v3-annotation-goal",
            "expression": {"h5ad": "experiments/EXP-20260929-001/work/input.h5ad",
                           "gene_name_col": "gene_name", "state_col": "subclass_l2"},
            "states": {"start": {"labels": ["FIB", "M-FIB"]}, "goal": {"labels": ["aFIB", "MYOF"]},
                       "reference": {"labels": ["dFIB"]}},
            "population": ["start", "goal"],
            "markers": {
                "goal_defining": [],          # **目標は遺伝子で定義していない**（著者の注釈）
                "goal_positive_control": sorted(by["positive_control"]),   # **71 個（主判定）**
                "reference_positive_control": ["MGP", "BGN", "DCN", "COL1A2", "C7", "AEBP1"]},
            "centroids": {
                # 2026-09-30 の着手前チェック（第 3 版）が使ったのと**同じ埋め込み**から作った重心。
                # analysis/centroids_from_embs_004.py が作る。**摂動は 1 回も使っていない**。
                # 重心は陽性対照の選び方に依存しないので、セットを替えても作り直す必要がない
                # （軸の長さの比が第 3 版と小数第 9 位まで一致することを確認済み）。
                "mode": "precomputed", "seeds": SEEDS,
                "files": [f"{ED}/work/axisfeas_centroids-seed{{seed}}.npz"],
                "names": {"start": "FIB_MFIB", "goal": "aStr", "reference": "dFIB",
                          "others": ["dMFIB"]}},
            "perturbation": {"min_cells_per_gene": 120},
            "thresholds": {},
            "output": f"{ED}/work/axis_feasibility.json"},
        "stats": {"mode": "goal_state_shift", "fdr_threshold": 0.05},
        "robustness": {"seeds": SEEDS, "subsample_frac": 1.0, "emb_layers": [-1]},
        "output": {"results_dir": "results", "experiment_dir": ED},
        "compute": {"device": "mps", "dtype": "bf16", "cache_model": True,
                    "sec_per_cell_gene": 0.1154,
                    "n_gene_runs": len(t1t2) + len(by["positive_control"])
                    + len(by["housekeeping"])
                    + len(by["housekeeping_extra"]) + len(by["negctrl_panel"]) + len(null500),
                    "n_cell_gene": None, "gpu_hours_estimate": None}}
    n_gene_runs = cfg["compute"]["n_gene_runs"]
    cfg["compute"]["n_cell_gene"] = n_gene_runs * 200 * len(SEEDS)
    cfg["compute"]["gpu_hours_estimate"] = round(
        cfg["compute"]["n_cell_gene"] * 0.1154 / 3600
        + len(SEEDS) * (950 * 2 + 270 + 135) * 0.1154 / 2 / 3600, 2)

    head = ("# EXP-20260930-001 — すべてのパラメータをここに書く（plan.md §6 が事前登録の本体）\n"
            "# 作成: analysis/make_configs_004.py（手で書き換えない）\n"
            "# 主指標: (目標への接近 − 出発点への接近) / 移動量。**退化した状態の重心は差し引きに入れない**\n"
            "#         （深さがそろっていないため。guardrails 2026-09-30）\n"
            "# 交絡の補正: lp・ll・ld（= log2 移動量）の完全二次式。移動量 0 の細胞は外す\n")
    (E / "config.yaml").write_text(head + yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False))

    for s in SEEDS:
        c2 = copy.deepcopy(cfg)
        c2["robustness"]["seeds"] = [s]
        c2["output"]["results_dir"] = f"results-seed{s}"
        (E / f"config-seed{s}.yaml").write_text(
            f"# EXP-20260930-001 seed {s} 専用（seeds と results_dir 以外は config.yaml と同一）\n"
            + yaml.safe_dump(c2, allow_unicode=True, sort_keys=False))

    sm = copy.deepcopy(cfg)
    sm["robustness"]["seeds"] = [SEEDS[0]]
    sm["output"]["results_dir"] = "results-smoke"          # **本番の results/ を汚さない**
    sm["cell_selection"]["min_cells_per_gene"] = 20        # guardrails 2026-09-28
    sm["cell_selection"]["max_cells_per_gene"] = 50
    sm["perturbation"]["genes_to_perturb"] = t1t2[:2]
    sm["controls"] = {"positive": cfg["controls"]["positive"][:2],
                      "negative_fixed": cfg["controls"]["negative_fixed"],
                      "fixed_role": "housekeeping",
                      "extra_fixed": [{"role": "negctrl_panel",
                                       "genes": cfg["controls"]["extra_fixed"][0]["genes"][:1]}],
                      "negative_random_n": 2,
                      "negative_sampler": {"method": "explicit", "genes": null500[:2]},
                      "comparator": None}
    (E / "work/config-smoke.yaml").write_text(
        "# 試し実行専用（**事前登録の対象外**）。min_cells_per_gene を 20 に下げ、results-smoke/ に書く\n"
        + yaml.safe_dump(sm, allow_unicode=True, sort_keys=False))

    print("遺伝子の総数:", n_gene_runs, "／延べ 細胞×遺伝子:", cfg["compute"]["n_cell_gene"])
    print("実測モデルでの見積もり:", cfg["compute"]["gpu_hours_estimate"], "時間（1 台で 3 seed）")
    print("書いた:", [str(p.relative_to(ROOT)) for p in
                      [E / "config.yaml"] + [E / f"config-seed{s}.yaml" for s in SEEDS]
                      + [E / "work/config-smoke.yaml"]])


if __name__ == "__main__":
    main()

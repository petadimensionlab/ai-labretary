#!/usr/bin/env python3
"""EXP-20261001-001 — 設定ファイルを 1 か所から機械的に作る（遺伝子リストを人が写さない）。

出力:
  config.yaml                    本体（3 seed 分。**事前登録の対象**）
  config-seed{81,82,83}.yaml     seed 単位で計算機を分けるとき（seeds と results_dir 以外は同一）
  work/config-smoke.yaml         試し実行専用（**事前登録の対象外**。min_cells_per_gene を下げる）

**遺伝子で分けてはいけない。分けるのは seed 単位のみ**（機体差が結論と同じ形で出るため）。

EXP-20260930-001 の設定から変えるのは次の 3 点だけである（PI 2026-10-01。plan.md §6.1）。

 A  主指標を**向きへの射影**にする（2 つの寄りの差を使わない）
 A-2 交絡を除く回帰に**移動量 `ld` を入れない**（`variables: [lp, ll]`）
 A-3 一様成分を**主指標の射影に使っていない状態**（dFIB・dMFIB）から作る
"""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path

import yaml

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20261001-001"
ED = "experiments/EXP-20261001-001"
P30 = "experiments/EXP-20260930-001"
SEEDS = [81, 82, 83]

# 盲検中（**解除されていない 3 個**。LTBR・CDKN1A は PI が 2026-09-27 に解除済み）
BLOCKED = {"TP53": "ENSG00000141510", "SPP1": "ENSG00000118785", "LCN2": "ENSG00000148346"}

# 絶対値の条件は **4.5 / √K を小数第 2 位に丸める**（選定手順 §1。EXP-20260930-001 と同じ式）。
# **人が数字を書かない**（写し間違いを避ける）。中立性のセットには絶対値の条件を置かない。
NO_ABS_MIN = ("neutral_adopted", "neutral_adopted_cv055", "housekeeping", "negctrl_panel")


def abs_min(k: int) -> float:
    return round(4.5 / math.sqrt(k), 2)


def main() -> None:
    g = json.loads((E / "work/genes_005.json").read_text())
    ens = {x["symbol"]: x["ensembl"] for x in g["genes"]}
    by = g["by_role"]
    sizes = json.loads((E / "work/subsets_005.json").read_text())["sizes"]
    null500 = [ens[s] for s in by["negative"]]
    # **摂動対象としてスクリプトに渡すのは介入セットと前回の摂動対象 T1・T2 の 135 個**
    # （陽性対照・中立・対照は controls の側に渡す。run_isp.py の見積もりが二重に数えない）
    perturb = [ens[s] for s in by["intervention"] + by["target_T1"] + by["target_T2"]]

    SET_ABS_MIN = {k: abs_min(v) for k, v in sizes.items() if k not in NO_ABS_MIN}

    cfg = {
        "experiment_id": "EXP-20261001-001",
        "hypothesis_id": "HYP-20261001-01",
        "model": {"version": "V2", "directory": "models/Geneformer/Geneformer-V2-316M",
                  "type": "Pretrained", "num_classes": 0, "emb_layer": -1},
        "data": {
            # **読み取りのみ**。トークン化済み入力が work/ にあるので h5ad は読まれない
            "input_h5ad": "experiments/EXP-20260929-001/work/input.h5ad",
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
                         "controls_only": False, "genes_to_perturb": perturb, "combos": 0,
                         "anchor_gene": None, "emb_mode": "cls", "forward_batch_size": 8,
                         "nproc": 1},
        "blocked_genes": list(BLOCKED.values()),
        "controls": {
            "positive": [ens[s] for s in by["positive_control"]],      # **71 個**
            "negative_fixed": [ens[s] for s in by["housekeeping"]],
            "fixed_role": "housekeeping",
            "extra_fixed": [
                {"role": "negctrl_panel", "genes": [ens[s] for s in by["negctrl_panel"]]},
                {"role": "housekeeping_extra", "genes": [ens[s] for s in by["housekeeping_extra"]]}],
            "negative_random_n": len(null500),
            "negative_sampler": {"method": "explicit", "genes": null500},
            "comparator": None},
        "capture": {"embeddings": True, "classifier_logits": False},
        # ------------------------------------------------------------------ 主指標（**PI 2026-10-01**）
        "metric": {
            # **A-2: 移動量 `ld` を入れない**（movement_correction.py の use_dnorm=False）
            "confound_correction": {"variables": ["lp", "ll"], "form": "full_quadratic",
                                    "drop_zero_movement": True},
            # **A: 向きへの射影**。u = 目標の重心 / |目標の重心| − 出発点の重心 / |出発点の重心| を
            # 正規化した向き。主指標 = (摂動後 − 摂動前) · û / ‖摂動後 − 摂動前‖
            "primary_form": "projection_on_start_to_goal_direction_over_movement",
            "projection_states": {"from": "start", "to": "goal"},
            # **A-3: 一様成分は、主指標の射影に使っていない状態から作る**
            "uniform_component_states": ["dFIB", "dMFIB"]},
        "analysis": {
            "primary_metric": "projection_goal_direction_over_movement",
            "primary_unit": "gene_set",
            "primary_set": "INT",                      # **介入セット 119 個**
            "primary_positive_control": "positive_control71",
            "primary_neutral_set": "neutral_adopted",  # **規則 H-1〜H-3 を通ったものだけ**
            "primary_null_version": "B_500",
            "primary_test": "set_mean_vs_random_control_sets",
            "pos_definition": "index_in_input_ids",
            "projection_centroids": ["start", "goal"],
            "uniform_centroids": ["dFIB", "dMFIB"],
            "descriptive_only_centroids": ["dFIB", "dMFIB"],
            "set_abs_min": SET_ABS_MIN,
            "sensitivity_sets": ["INT_S1", "INT_S2", "INT_S3", "INT_pathway_specific",
                                 "INT_plus_T1T2", "positive_control65"],
            # 検出力が足りないセット（**不合格を「効果がない」の根拠に使わない**。plan.md §9）
            "underpowered_sets": ["INT_S2"],
            "secondary_null_versions": ["A_300"],
            "pos_matched_permutation": {"tolerance_log2": 0.5, "min_cover_frac": 0.8,
                                        "min_eligible_null_genes": 20},
            "criteria": {
                "require_all_seeds": True, "K_set": sizes["INT"],
                "n_permutation": 200000, "set_percentile": 0.99, "set_draw_replace": True,
                "n_chance_draw": 200000, "max_chance_rate": 0.02,
                "neutral_range": [0.05, 0.95], "neutral_min_seeds": 2,
                "n_sham_sets": 200, "sham_match_bias": True, "sham_centre_tol": 0.10,
                "sham_spread_range": [0.7, 1.3], "sham_binomial_alpha": 0.05,
                "null_draw_seed": 20261001,
                "M1_tail_frac_max": 0.20, "M1_min_n": 60,
                "M2_abs_rho_max": 0.25, "M2_min_seeds": 2,
                "M3_range": [0.10, 0.90], "M3_min_seeds": 2,
                # 採用した中立対照 4 個のうち 3 個以上（前回のパネル 5/7 = 71% と同じ割合）
                "M3_adopted_min_genes": 3,
                "M4_abs_rho_uniform_max": 0.40,        # **閾値は動かさない**（PI 2026-10-01）
                "M5_abs_rho_movement_max": 0.10}},
        "axis_feasibility": {
            # **選定手順 §7 で 2026-10-01T19:39 に実行済み**（判定: 判定保留、不合格の条件 0 件）。
            # 目標・出発点・陽性対照 71 個・seed を 1 つも変えていないので、同じ値が出る。
            "experiment_id": "EXP-20261001-001-precheck",
            "expression": {"h5ad": "experiments/EXP-20260929-001/work/input.h5ad",
                           "gene_name_col": "gene_name", "state_col": "subclass_l2"},
            "states": {"start": {"labels": ["FIB", "M-FIB"]}, "goal": {"labels": ["aFIB", "MYOF"]},
                       "reference": {"labels": ["dFIB"]}},
            "population": ["start", "goal"],
            "markers": {
                "goal_defining": [],          # **目標は遺伝子で定義していない**（著者の注釈）
                "goal_positive_control": sorted(by["positive_control"]),
                "reference_positive_control": ["MGP", "BGN", "DCN", "COL1A2", "C7", "AEBP1"]},
            "centroids": {
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
                    "n_gene_runs": g["n_total"],
                    "n_cell_gene": None, "gpu_hours_estimate": None}}

    n_gene_runs = cfg["compute"]["n_gene_runs"]
    cfg["compute"]["n_cell_gene"] = n_gene_runs * 200 * len(SEEDS)
    cfg["compute"]["gpu_hours_estimate"] = round(
        cfg["compute"]["n_cell_gene"] * 0.1154 / 3600
        + len(SEEDS) * (950 * 2 + 270 + 135) * 0.1154 / 2 / 3600, 2)

    head = ("# EXP-20261001-001 — すべてのパラメータをここに書く（plan.md §6 が事前登録の本体）\n"
            "# 作成: analysis/make_configs_005.py（手で書き換えない）\n"
            "# 主指標: 変位を『出発点の重心 → 目標の重心』の向きへ射影し、移動量で割ったもの\n"
            "#         （**2 つの寄りの差は使わない**。PI 2026-10-01）\n"
            "# 交絡の補正: lp・ll の完全二次式。**移動量 ld は入れない**（PI 2026-10-01）\n"
            "# 一様成分: 主指標の射影に使っていない状態（dFIB・dMFIB）の寄りの平均\n")
    (E / "config.yaml").write_text(head + yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False))

    for s in SEEDS:
        c2 = copy.deepcopy(cfg)
        c2["robustness"]["seeds"] = [s]
        c2["output"]["results_dir"] = f"results-seed{s}"
        (E / f"config-seed{s}.yaml").write_text(
            f"# EXP-20261001-001 seed {s} 専用（seeds と results_dir 以外は config.yaml と同一）\n"
            + yaml.safe_dump(c2, allow_unicode=True, sort_keys=False))

    sm = copy.deepcopy(cfg)
    sm["robustness"]["seeds"] = [SEEDS[0]]
    sm["output"]["results_dir"] = "results-smoke"          # **本番の results/ を汚さない**
    sm["cell_selection"]["min_cells_per_gene"] = 20        # guardrails 2026-09-28
    sm["cell_selection"]["max_cells_per_gene"] = 50
    # 試し実行は**評価コードを最後まで通すためのもの**なので、
    # 判定に使う役割を 1 つずつ必ず含める（介入・陽性対照・採用した中立対照・対照）。
    neu = json.loads((E / "work/neutral_rule_005.json").read_text())["採用"]
    sm["perturbation"]["genes_to_perturb"] = [ens[s] for s in by["intervention"][:3]]
    sm["controls"] = {"positive": cfg["controls"]["positive"][:3],
                      "negative_fixed": cfg["controls"]["negative_fixed"],
                      "fixed_role": "housekeeping",
                      "extra_fixed": [{"role": "negctrl_panel",
                                       "genes": [ens[x] for x in neu if x in by["negctrl_panel"]]},
                                      {"role": "housekeeping_extra",
                                       "genes": cfg["controls"]["extra_fixed"][1]["genes"]}],
                      "negative_random_n": 30,
                      "negative_sampler": {"method": "explicit", "genes": null500[:30]},
                      "comparator": None}
    (E / "work/config-smoke.yaml").write_text(
        "# 試し実行専用（**事前登録の対象外**）。min_cells_per_gene を 20 に下げ、results-smoke/ に書く\n"
        "# **評価コードを最後まで通すことが目的**なので、判定に使う役割をすべて 1 つ以上含める\n"
        + yaml.safe_dump(sm, allow_unicode=True, sort_keys=False))

    print("遺伝子の総数:", n_gene_runs, "／延べ 細胞×遺伝子:", cfg["compute"]["n_cell_gene"])
    print("実測モデルでの見積もり:", cfg["compute"]["gpu_hours_estimate"], "時間（1 台で 3 seed）")
    print("絶対値の条件 (4.5/√K):", json.dumps(SET_ABS_MIN, ensure_ascii=False))
    print("書いた:", [str(p.relative_to(ROOT)) for p in
                      [E / "config.yaml"] + [E / f"config-seed{s}.yaml" for s in SEEDS]
                      + [E / "work/config-smoke.yaml"]])


if __name__ == "__main__":
    main()

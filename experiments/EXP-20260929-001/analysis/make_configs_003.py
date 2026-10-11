#!/usr/bin/env python3
"""EXP-20260929-001: config.yaml と機体ごとの設定を 1 つの遺伝子リストから機械的に生成する。

**遺伝子で分割してはいけない**（機体差が「分泌プログラムは対照より動く」という結論と
まったく同じ形で現れ、見分けられなくなる。labretary-level 2026-09-28）。
分けるのは seed 単位だけで、生成する設定は `robustness.seeds` と `output.results_dir` 以外
完全に同一である。

生成物:
  config.yaml                      3 seed をまとめた基準の設定（dry-run と事前登録の対象）
  config-seed{61,62,63}.yaml       機体ごと（seed 単位の分割）
  work/config-smoke.yaml           試し実行専用（**事前登録の対象外**）
  work/config-machinecheck.yaml    機体の一致確認専用
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260929-001"
SEEDS = [61, 62, 63]
MACHINE = {61: "spark-cf5b", 62: "spark-c083", 63: "spark-c083"}
SEC_PER_CELL_GENE = 0.1154        # EXP-20260928-001 の実測（93,600 細胞×遺伝子で 3.00 時間）


def ylist(items, indent: int) -> str:
    pad = " " * indent
    return "\n".join(f"{pad}- {x}" for x in items)


def build(seeds, results_dir, genes, smoke=False, machinecheck=False) -> str:
    g = genes
    n_gene = (len(g["target"]) + len(g["marker_myof"]) + len(g["marker_org"])
              + len(g["housekeeping"]) + len(g["negctrl_panel"]) + len(g["reference_myof"])
              + len(g["negative"]))
    n_cell_gene = n_gene * 200 * len(seeds)
    hours = round(n_cell_gene * SEC_PER_CELL_GENE / 3600, 3)
    head = "# EXP-20260929-001 — すべてのパラメータをここに書く（plan.md §6 が事前登録の本体）\n"
    if len(seeds) == 1 and not smoke and not machinecheck:
        head = (f"# EXP-20260929-001 — seed {seeds[0]} を {MACHINE.get(seeds[0], '?')} で回すための設定"
                f"（analysis/make_configs_003.py が自動生成）\n"
                f"# config.yaml との違いは robustness.seeds と output.results_dir だけ。\n"
                f"# 遺伝子リスト（分泌 24 / 目標マーカー 7 / 筋線維芽細胞 24 / 恒常発現 3 /\n"
                f"# 中立パネル 7 / 参照 2 / 対照 500 = {n_gene}）は 3 つの設定で完全に同一である。\n"
                f"# **遺伝子で分割してはいけない**（機体差が結論と同じ形で現れる）。\n")
    return f"""{head}experiment_id: EXP-20260929-001
hypothesis_id: HYP-20260929-01

model:                              # EXP-20260927-003 / 004 / EXP-20260928-001 / 002 と同一（比較可能性のため変えない）
  version: V2
  directory: models/Geneformer/Geneformer-V2-316M
  type: Pretrained
  num_classes: 0
  emb_layer: -1

data:
  input_h5ad: experiments/EXP-20260929-001/work/input.h5ad     # analysis/build_input_003.py が作る。data/raw は読むだけ
  tokenized_dataset: experiments/EXP-20260929-001/work/tokenized.dataset
  cell_type_col: subclass_l2
  state_key: fib_org_state          # 目標・出発点・筋線維芽細胞の 3 値。対照の状態は obs の ecm_high で持つ
  custom_attr_name_dict: {{fib_org_state: fib_org_state, ecm_high: ecm_high, donor_id: donor_id,
                          subclass_l2: subclass_l2, disease: disease, region: region,
                          n_genes: n_genes, n_counts: n_counts, a1_row: a1_row}}
  filter_data: null
  max_ncells: null

states:
  start_state: FIB_orgLow           # 18,099 細胞 / 79 提供者
  goal_state: FIB_orgHigh           # 950 細胞 / 70 提供者
  alt_states: []                    # 対照の状態（FIB_ecmHigh）と筋線維芽細胞（MYOF_fib）は
                                    # extra_centroids_003.py が作る（出発点の細胞プールを削らないため）

state_embs:
  max_cells_per_state: 950
  max_cells_per_donor: 60
  rows_file: experiments/EXP-20260929-001/work/state_rows-seed{{seed}}.json   # 深さをそろえた行番号

cell_selection:
  donor_key: donor_id
  max_cells_per_gene: 200
  max_cells_per_donor: 20
  min_cells_per_gene: {5 if smoke else 120}

perturbation:
  perturb_type: overexpress         # 発現順位を最上位へ動かす。削除型は使わない
  perturb_rank_shift: null
  controls_only: false
  # 分泌プログラムのセット K = {len(g['target'])}（work/subsets_003.json の secretory_all と同一）。
  # **主判定はこのセットの「平均」である**（個々の遺伝子ではない）
  genes_to_perturb:
{ylist(g['target'], 4)}
  combos: 0                         # 遺伝子は 1 個ずつ。同時摂動は使わない
  anchor_gene: null
  emb_mode: cls                     # V2 は cls 必須（guardrails 2026-09-25）
  forward_batch_size: 8
  nproc: 1                          # GPU を使うプロセスで datasets の fork を避ける（guardrails 2026-09-26）

blocked_genes:                      # 盲検中。摂動しない。evaluate_003.py が読み込み直後に除外して停止条件にする
{ylist(g['blinded'], 2)}

controls:
  # 筋線維芽細胞マーカー 24（第 2 の軸の陽性対照。EXP-20260928-001 / 002 と同一の遺伝子）
  positive:
{ylist(g['marker_myof'], 4)}
  # 恒常発現遺伝子（中立性の検査）
  negative_fixed:
{ylist(g['housekeeping'], 4)}
  fixed_role: housekeeping
  extra_fixed:
    - role: marker_org              # **目標の軸の陽性対照**（目標の定義に使ったマーカー 7）。
      genes:                        # **健全性チェックであって、仮説を支持する証拠には使わない**
{ylist(g['marker_org'], 8)}
    - role: negctrl_panel           # 中立対照パネル 7（EXP-20260926-003 以降で共通）
      genes:
{ylist(g['negctrl_panel'], 8)}
    - role: reference_myof          # 参照（既測定。判定には使わない）
      genes:
{ylist(g['reference_myof'], 8)}
  # 対照遺伝子 500（**主判定**。work/null_pool_003.json の B_500 と同一）。
  # **先頭 300 が入れ子の部分集合**（A_300）で、感度分析に使う。
  # negative_random_n を 300 に変えるだけで段階実行できる（run_isp.py はリストの先頭から取る）。
  negative_random_n: {len(g['negative'])}
  negative_sampler:
    method: explicit
    genes:
{ylist(g['negative'], 6)}
  comparator: null

capture:
  embeddings: true                  # 元・摂動後の CLS 埋め込みを保存（3 本の軸を後から計算する）
  classifier_logits: false

extra_centroids:
  output: experiments/EXP-20260929-001/work/extra_centroids-seed{{seed}}.npz
  sets:                             # 細胞は work/state_rows-seed{{S}}.json の extra_rows で指定済み
    FIB_ecmHigh: {{note: 対照の状態（第 3 の軸）}}
    MYOF_fib: {{note: 筋線維芽細胞（第 2 の軸）}}
    FIB_orgHigh_unmatched: {{note: 感度分析（深さをそろえない目標の重心。判定に使わない）}}

analysis:                           # plan.md §6 の事前登録（**結果を見る前に確定。結果を見てから変えない**）
  primary_metric: A12               # (軸 − 他の 3 重心の平均) / 移動量。軸を入れ替えて 3 本つくる
  primary_unit: gene_set            # 判定の単位はセット
  primary_set: secretory_all        # 主判定の対象は分泌プログラム {len(g['target'])} 個
  primary_null_version: B_500       # 主判定は対照 500 版（事前チェック P-C の結果。plan.md §5.3）。
                                    # 300 版は入れ子の感度分析（500 の先頭 300 が 300 版）
  pos_definition: index_in_input_ids
  primary_test: set_mean_vs_random_control_sets
  axes: [goal, myofaxis, placebo]   # **3 本すべてで検定する**（EXP-20260928-002 validation.md §3.4(b)）
  pos_matched_permutation:
    bias_model: quadratic_lp_ll
    distance: max_abs_log2
    tolerance_log2: 0.5
    min_cover_frac: 0.8
    min_eligible_null_genes: 20
    statistic: D_raw
  criteria:
    require_all_seeds: true         # 3 seed すべて
    K_set: {len(g['target'])}
    n_permutation: 200000
    set_percentile: 0.99            # 無作為な対照セットの 99 パーセンタイルを超えること
    set_draw_replace: true          # **復元ありで引く**（guardrails 2026-09-28）。復元なしも併記する
    set_abs_min: 0.46               # 併せて満たす絶対値。HYP §7 が K=25 で見込んだ 99 パーセンタイル。
                                    # EXP-20260928-002 の対照 980 個で偶然合格率 0.0053（≤ 0.02）を確認済み
    n_chance_draw: 200000
    max_chance_rate: 0.02
    neutral_range: [0.05, 0.95]
    neutral_min_seeds: 2
    n_sham_sets: 200                # 偽のセット（guardrails 2026-09-28 の下限 50 を上回る）
    sham_match_bias: true
    sham_centre_tol: 0.10           # 偽セットの中心が無作為セットからずれてよい幅
    sham_spread_range: [0.7, 1.3]   # 偽セットの散らばり / 無作為セットの散らばり の許容範囲
    null_draw_seed: 20260930
    # --- 指標の健全性（EXP-20260927-004 / EXP-20260928-001 / 002 と同一の閾値。変えない）
    M1_tail_frac_max: 0.20
    M1_min_n: 60
    M2_abs_rho_max: 0.25
    M2_min_seeds: 2
    M3_range: [0.10, 0.90]
    M3_min_seeds: 2
    M3_panel_min_genes: 5
    M4_abs_rho_uniform_max: 0.40
    M4_abs_rho_placebo_max: 0.50
  secondary_sets:                   # 副次的な部分セット（work/subsets_003.json。結果を見る前に固定）
    - co_or_low                     # 共起オッズ比 2.0 未満（**循環でない証拠の要**。HYP P-6 / F-8）
    - co_or_high                    # 同 2.0 以上（5 個。記述値）
    - tls_chemokine                 # リンパ組織のケモカイン部分（4 個。HYP P-4。記述値）
    - non_tls                       # それ以外（20 個）
  secondary_null_versions: [A_300]  # 対照 300 版（入れ子の部分集合）での再判定
  descriptive_only:
    - single_gene_2p5sd
    - bias_vs_effect
    - reference_myof
    - donor_level
    - goal_unmatched_sensitivity
    - machine_per_seed

# run_isp.py が必須とするセクション。本実験の判定は analysis/evaluate_003.py が
# 埋め込みから直接計算するので、ここの出力は判定に使わない。
stats:
  mode: goal_state_shift
  fdr_threshold: 0.05

robustness:
  seeds: {json.dumps(seeds)}
  subsample_frac: 1.0
  emb_layers: [-1]

output:
  results_dir: {results_dir}
  experiment_dir: experiments/EXP-20260929-001

compute:
  device: mps                       # 本番はリモート計算機。isp_remote.sh が cuda で上書きする
  dtype: bf16                       # PI 標準（2026-09-26）。コサインは fp32（guardrails 2026-09-26）
  cache_model: true
  sec_per_cell_gene: {SEC_PER_CELL_GENE}
  n_gene_runs: {n_gene}
  n_cell_gene: {n_cell_gene}
  gpu_hours_estimate: {hours}
"""


def main() -> None:
    gj = json.loads((E / "work/genes_003.json").read_text())
    ens = {g["symbol"]: g["ensembl"] for g in gj["genes"]}
    by_role = gj["by_role"]
    npool = json.loads((E / "work/null_pool_003.json").read_text())["versions"]
    subs = json.loads((E / "work/subsets_003.json").read_text())["subsets"]
    groups = json.loads((E / "analysis/excluded_symbols_003.json").read_text())["groups"]
    sym2ens = json.loads((E / "work/symbol_to_ensembl_003.json").read_text())

    order = {
        "target": [x["ensembl"] for x in subs["secretory_all"]],
        "marker_org": [x["ensembl"] for x in subs["org_markers7"]],
        "marker_myof": [x["ensembl"] for x in subs["myof_markers24"]],
        "housekeeping": [x["ensembl"] for x in subs["housekeeping3"]],
        "negctrl_panel": [x["ensembl"] for x in subs["negctrl_panel7"]],
        "reference_myof": sorted(ens[s] for s in by_role["reference_myof"]),
        "negative": [x["ensembl"] for x in npool["B_500"]],
        "blinded": [sym2ens[s] for s in groups["blinded5"]],
    }
    allg = sum((order[k] for k in ("target", "marker_org", "marker_myof", "housekeeping",
                                   "negctrl_panel", "reference_myof", "negative")), [])
    assert len(allg) == len(set(allg)) == 567, f"遺伝子数が違う: {len(allg)}"
    assert not (set(allg) & set(order["blinded"])), "盲検中の遺伝子が入っている"

    (E / "config.yaml").write_text(build(SEEDS, "results", order))
    for s in SEEDS:
        (E / f"config-seed{s}.yaml").write_text(build([s], f"results-seed{s}", order))

    smoke = dict(order)
    smoke["target"] = order["target"][:2]
    smoke["marker_org"] = order["marker_org"][:1]
    smoke["marker_myof"] = order["marker_myof"][:1]
    smoke["negctrl_panel"] = order["negctrl_panel"][:1]
    smoke["reference_myof"] = order["reference_myof"][:1]
    smoke["negative"] = order["negative"][:2]
    (E / "work/config-smoke.yaml").write_text(
        build([61], "results-smoke", smoke, smoke=True))

    mc = dict(order)
    mc["negative"] = order["negative"][:20]
    (E / "work/config-machinecheck.yaml").write_text(
        build([61], "results-machinecheck", mc, machinecheck=True))

    print(f"生成: config.yaml, config-seed{{{','.join(map(str, SEEDS))}}}.yaml, "
          f"work/config-smoke.yaml, work/config-machinecheck.yaml")
    print(f"摂動する遺伝子 {len(allg)}（3 seed で {len(allg) * 200 * 3:,} 細胞×遺伝子 = "
          f"{len(allg) * 200 * 3 * SEC_PER_CELL_GENE / 3600:.2f} 時間）")


if __name__ == "__main__":
    main()

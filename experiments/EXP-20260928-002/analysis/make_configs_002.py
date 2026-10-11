"""EXP-20260928-002 — config.yaml と、機体・段階ごとの設定ファイルを生成する。

生成物はすべて事前登録の対象（SHA256SUMS-preregistered.txt に入れる）。
手で書き写すと 1,136 行の遺伝子リストで写し間違いが起きるので、work/genes_002.json から機械的に作る。

  config.yaml                     正本（3 seed すべて・対照 980）。dry-run と、統合後の評価はこれを使う
  config-seed{S}.yaml             機体ごと（seed 1 つ・対照 980・results-seed{S}/ に書く）
  config-seed{S}-stage{1,2}.yaml  段階実行（対照 520 / 600 まで。同じ results-seed{S}/ に書き、再開で続けられる）
  work/config-smoke.yaml          試し実行専用（事前登録の対象外。results-smoke/ に書く）

実行: .venv/bin/python experiments/EXP-20260928-002/analysis/make_configs_002.py
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260928-002"
EXP = "EXP-20260928-002"
SEEDS = [51, 52, 53]
STAGE_N = {1: 520, 2: 600, 3: 980}
MACHINE = {51: "spark-cf5b (petadimensionlab@100.104.171.109)",
           52: "spark-c083 (coi-next1@100.94.57.17)",
           53: "spark-c083 (coi-next1@100.94.57.17)"}

BLOCKED = [("ENSG00000111321", "LTBR   盲検（未解除）・元の主対象"),
           ("ENSG00000104856", "RELB   元の主対象。セットに入れない"),
           ("ENSG00000077150", "NFKB2  元の主対象。セットに入れない"),
           ("ENSG00000124762", "CDKN1A 盲検（未解除）"),
           ("ENSG00000141510", "TP53   盲検（未解除）"),
           ("ENSG00000118785", "SPP1   盲検（未解除）"),
           ("ENSG00000148346", "LCN2   盲検（未解除）"),
           ("ENSG00000172724", "CCL19  目標ラベルを定義する遺伝子（循環）"),
           ("ENSG00000137077", "CCL21  同上"),
           ("ENSG00000156234", "CXCL13 同上")]


def gene_lines(genes: list[dict], role: str, indent: str = "    ") -> str:
    return "\n".join(f"{indent}- {g['ensembl']}  # {g['symbol']}" for g in genes if g["role"] == role)


def build(genes: list[dict], seeds: list[int], n_neg: int, results_dir: str,
          header: str, sec_per_cell_gene: float) -> str:
    n_genes = 120 + 36 + n_neg
    n_cell_gene = n_genes * 200 * len(seeds)
    hours = round(n_cell_gene * sec_per_cell_gene / 3600, 2)
    return f"""{header}
experiment_id: {EXP}
hypothesis_id: HYP-20260928-01

model:                              # EXP-20260928-001 と同一（比較可能性のため変えない）
  version: V2
  directory: models/Geneformer/Geneformer-V2-316M
  type: Pretrained
  num_classes: 0
  emb_layer: -1

data:
  input_h5ad: experiments/EXP-20260927-002/work/input.h5ad      # 読み取りのみ。data/raw は触らない
  # EXP-20260928-001 の work/tokenized.dataset のデータ部分を cp -c で複製したもの（cache-*.arrow は複製しない。
  # guardrails 2026-09-26）。3 ファイルの sha256 が EXP-20260928-001 の事前登録と一致することを確認済み。
  # 筋線維芽細胞（MYOF / cycMYOF）の状態ラベルは `MYOF_<元の状態>` に付け替え済み（guardrails 2026-09-28）。
  tokenized_dataset: experiments/{EXP}/work/tokenized.dataset
  cell_type_col: subclass_l2
  state_key: fib_tls_state
  custom_attr_name_dict: {{fib_tls_state: fib_tls_state, donor_id: donor_id, subclass_l2: subclass_l2,
                          disease: disease, region: region, TLS_group: TLS_group,
                          n_genes: n_genes, n_counts: n_counts}}
  filter_data: null
  max_ncells: null

states:
  start_state: FIB_neg_TLShigh
  goal_state: FIB_org
  alt_states: [FIB_neg_TLSlow]

state_embs:
  max_cells_per_state: 800
  max_cells_per_donor: 60
  rows_file: experiments/{EXP}/work/state_rows-seed{{seed}}.json   # MYOF を除き深さをそろえた行番号

cell_selection:
  donor_key: donor_id
  max_cells_per_gene: 200
  max_cells_per_donor: 20
  min_cells_per_gene: 120

perturbation:
  perturb_type: overexpress
  perturb_rank_shift: null
  controls_only: false
  # 経路セット K = 120（work/pathway_set_002.json と同一。plan.md §4）
  # **本実験の主判定はこのセットの「平均」である**（個々の遺伝子ではない）
  genes_to_perturb:
{{PATHWAY}}
  combos: 0                         # 遺伝子は 1 個ずつ回す。同時摂動は使わない
  anchor_gene: null
  emb_mode: cls                     # V2 は cls 必須（guardrails 2026-09-25）
  forward_batch_size: 8
  nproc: 1                          # GPU 使用プロセスでは datasets の fork を避ける（guardrails 2026-09-26）

blocked_genes:                      # 摂動しない。evaluate_002.py が読み込み直後に除外して停止条件にする
{chr(10).join(f"  - {e}{' ' * max(1, 20 - len(e))}# {n}" for e, n in BLOCKED)}

controls:
  # 筋線維芽細胞マーカーセット 24（アッセイが働いていることの陽性対照。EXP-20260928-001 と同一）
  positive:
{{MYOF}}
  # 恒常発現遺伝子（中立性の検査）
  negative_fixed:
    - ENSG00000075624  # ACTB
    - ENSG00000111640  # GAPDH
    - ENSG00000196262  # PPIA
  fixed_role: housekeeping
  extra_fixed:
    - role: negctrl_panel          # 中立対照パネル 7（EXP-20260926-003/005・EXP-20260927-001/004・
      genes:                       # EXP-20260928-001 と同一）
        - ENSG00000090661  # CERS4
        - ENSG00000147687  # TATDN1
        - ENSG00000081870  # IFT25
        - ENSG00000101558  # VAPA
        - ENSG00000054148  # PHPT1
        - ENSG00000166295  # ANAPC16
        - ENSG00000104529  # EEF1D
    - role: reference_myof         # 参照（既測定。判定には使わない）
      genes:
        - ENSG00000107796  # ACTA2
        - ENSG00000149591  # TAGLN
  # 対照遺伝子。上位集合 980 を順番に並べてある（work/null_pool_002.json と同一）。
  #   先頭 520 = 案 C 相当 / 先頭 600 = 案 B 相当 / 全 980 = 主判定（案 A）。部分集合は入れ子。
  #   negative_random_n を 520 / 600 / 980 に変えるだけで段階実行できる（run_isp.py はリストの先頭から取る）。
  negative_random_n: {n_neg}
  negative_sampler:
    method: explicit
    genes:
{{NEG}}
  comparator: null

capture:
  embeddings: true                  # 元・摂動後の CLS 埋め込みを保存（軸ごとの接近度を後から計算する）
  classifier_logits: false

extra_centroids:
  output: experiments/{EXP}/work/extra_centroids-seed{{seed}}.npz
  sets:
    MYOF_neg:                       # 第 2 の軸: 筋線維芽細胞の重心
      states: [MYOF_FIB_neg_TLShigh, MYOF_FIB_neg_TLSlow]
      subclass_l2: [MYOF, cycMYOF]
      depth_match_to: FIB_neg_TLShigh
      max_cells: 800
      max_cells_per_donor: 60
    FIB_org_unmatched:              # 感度分析: 深さをそろえない目標重心（判定に使わない）
      states: [FIB_org]             # 状態ラベルを付け替えてあるので FIB_org に筋線維芽細胞は含まれない
      depth_match_to: null
      max_cells: 800
      max_cells_per_donor: 60

analysis:                           # plan.md §6 の事前登録（結果を見る前に確定。結果を見てから変えない）
  primary_metric: A12               # (軸 − 他の 3 重心の平均) / 移動量。軸を入れ替えて 3 本作る
  primary_unit: gene_set            # 判定の単位はセット
  primary_set: pathway_all          # 主判定の対象は経路セット 120
  primary_null_version: A_980       # 主判定は対照 980 版。600 / 520 版は副次的な報告
  pos_definition: index_in_input_ids
  primary_test: set_mean_vs_random_control_sets
  pos_matched_permutation:          # 遺伝子ごとの D_raw の作り方（EXP-20260927-003/004/EXP-20260928-001 と同じ実装）
    bias_model: quadratic_lp_ll
    distance: max_abs_log2
    tolerance_log2: 0.5
    min_cover_frac: 0.8
    min_eligible_null_genes: 20
    statistic: D_raw
  criteria:
    require_all_seeds: true         # 3 seed すべて
    K_set: 120                      # 経路セットの大きさ
    n_permutation: 200000           # 無作為な対照セットの引き直し回数（閾値の決定）
    set_percentile: 0.99            # セット平均 z が、無作為な対照セットの 99 パーセンタイルを超えること
    set_draw_replace: true          # **復元ありで引く**（guardrails 2026-09-28）。復元なしも併記する
    set_abs_min: 0.245              # 併せて満たす絶対値。EXP-20260927-004 の対照 84 個で K=24 のとき
                                    # 偶然合格率 5% になる値（別のデータから、結果を見る前に固定した）
    n_chance_draw: 200000           # 偶然合格率の実測に使う引き直し回数
    max_chance_rate: 0.02           # 偶然合格率の実測がこれ以下であること
    neutral_range: [0.05, 0.95]     # 恒常発現遺伝子・中立対照パネルのセットが入るべき分位
    neutral_min_seeds: 2
    n_sham_sets: 50                 # 偽経路セット（guardrails 2026-09-28。3 組では合格率を推定できない）
    sham_match_bias: true           # 偽経路セットは bias_goal もそろえる（HYP §6.2）
    null_draw_seed: 20260929
    # --- 指標の健全性（EXP-20260927-004 / EXP-20260928-001 と同一の閾値。変えない）
    M1_tail_frac_max: 0.20
    M1_min_n: 60
    M2_abs_rho_max: 0.25
    M2_min_seeds: 2
    M3_range: [0.10, 0.90]
    M3_min_seeds: 2
    M3_panel_min_genes: 5
    M4_abs_rho_uniform_max: 0.40
    M4_abs_rho_placebo_max: 0.50
  secondary_sets:                   # 副次的な部分セット（work/subsets_002.json。結果を見る前に固定）
    - nonmarker                     # |bias_goal|<0.5 かつ |bias_myof|<0.5 の経路遺伝子（循環でない証拠）
    - drop_goalmarker               # 目標マーカーの条件を満たす経路遺伝子を抜いた版
    - pathway_specific              # 一般機構（プロテアソーム・ユビキチン系）を抜いた版
    - general_mechanism             # 一般機構のみ（記述値。5 個しかないので判定には使わない）
    - tier_S1
    - tier_S2
    - tier_S3
  secondary_null_versions: [B_600, C_520]   # 対照の部分集合での再判定（案 B 相当・案 C 相当）
  descriptive_only:                 # 記録するが合否には入れない
    - single_gene_2p5sd
    - bias_vs_effect
    - reference_myof
    - donor_level
    - goal_unmatched_sensitivity
    - machine_per_seed

# run_isp.py が必須とするセクション（validate() と run_stats()）。
# 本実験の判定は analysis/evaluate_002.py が埋め込みから直接計算するので、ここの出力は判定に使わない。
stats:
  mode: goal_state_shift
  fdr_threshold: 0.05

robustness:
  seeds: {seeds}
  subsample_frac: 1.0
  emb_layers: [-1]

output:
  results_dir: {results_dir}
  experiment_dir: experiments/{EXP}

compute:
  device: mps                       # 実行はリモート計算機。isp_remote.sh が cuda で上書きする
  dtype: bf16                       # PI 標準（2026-09-26）。コサインは fp32（guardrails 2026-09-26）
  cache_model: true
  sec_per_cell_gene: {sec_per_cell_gene}          # EXP-20260928-001 の実測（93,600 細胞×遺伝子で 3.00 時間）
  sec_fixed_per_gene_run: 3.18
  sec_per_1e5_tokens: 4.49
  n_gene_runs: {n_genes * len(seeds)}
  n_cell_gene: {n_cell_gene}
  gpu_hours_estimate: {hours}
"""


def main() -> None:
    gj = json.loads((E / "work/genes_002.json").read_text())
    genes = gj["genes"]
    assert gj["seeds"] == SEEDS, gj["seeds"]
    path = gene_lines(genes, "pathway")
    myof = gene_lines(genes, "marker_myof")
    neg = gene_lines(genes, "negative", indent="      ")
    assert len([g for g in genes if g["role"] == "pathway"]) == 120
    assert len([g for g in genes if g["role"] == "marker_myof"]) == 24
    assert len([g for g in genes if g["role"] == "negative"]) == 980

    def write(p: Path, header: str, seeds: list[int], n_neg: int, res: str, spc: float = 0.1154) -> None:
        t = build(genes, seeds, n_neg, res, header, spc)
        t = t.replace("{PATHWAY}", path).replace("{MYOF}", myof).replace("{NEG}", neg)
        p.write_text(t)
        print(f"{p.relative_to(ROOT)}  seeds={seeds} 対照={n_neg} → {res}/")

    write(E / "config.yaml",
          "# EXP-20260928-002 — HYP-20260928-01 の検定（事前登録。plan.md を先に読むこと）\n"
          "#\n"
          "# これが正本。3 seed すべて・対照 980（案 A）。dry-run と、3 機体の結果を統合した後の評価に使う。\n"
          "# 実行は機体ごとに config-seed{S}.yaml を使う（seed 単位でのみ分割する。plan.md §8）。",
          SEEDS, 980, "results")

    for s in SEEDS:
        write(E / f"config-seed{s}.yaml",
              f"# EXP-20260928-002 — seed {s} を {MACHINE[s]} で回すための設定（自動生成）\n"
              f"# config.yaml との違いは robustness.seeds と output.results_dir だけ。\n"
              f"# 遺伝子リスト（経路 120 / 対照 980 / アッセイ対照 36）は 3 機体で完全に同一である。\n"
              f"# **遺伝子で分割してはいけない**（機体差が『経路セットは対照より動く』と同じ形で現れるため）。",
              [s], 980, f"results-seed{s}")
        for st in (1, 2):
            write(E / f"config-seed{s}-stage{st}.yaml",
                  f"# EXP-20260928-002 — seed {s} / 段階 {st}（対照 {STAGE_N[st]} まで）。自動生成\n"
                  f"# config-seed{s}.yaml との違いは controls.negative_random_n だけ。\n"
                  f"# 同じ results-seed{s}/ に書くので、次の段階は再開（済んだ遺伝子は飛ばす）で続けられる。",
                  [s], STAGE_N[st], f"results-seed{s}")

    # 試し実行専用（事前登録の対象外。本番とは別の results_dir に書く。guardrails 2026-09-28）
    (E / "work").mkdir(exist_ok=True)
    t = build(genes, [SEEDS[0]], 2, "results-smoke",
              "# EXP-20260928-002 — 試し実行（smoke）専用。**事前登録の対象外**。\n"
              "# 本番とは別の results_dir に書く（guardrails 2026-09-28）。\n"
              "# min_cells_per_gene は run_isp.py の --smoke が 5 に下げる（2026-09-28 の修正）。", 0.1154)
    t = t.replace("{PATHWAY}", path).replace("{MYOF}", myof).replace("{NEG}", neg)
    (E / "work/config-smoke.yaml").write_text(t)
    print("work/config-smoke.yaml  （事前登録の対象外）")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""EXP-20261005-001 の config.yaml を、EXP-20261004-001 の config.yaml から機械的に作る。

**遺伝子集合は 1 個も選び直さない。** EXP-20261004-001 の 626 個（候補 8・陽性対照 71・
対照 549、GSN と ELN が重複）をそのまま引き継ぐ。変えるのは次の 5 点だけ。

  1. seed を 84 / 85 / 86 にする（判定は新しい seed で行う。plan.md §4.2）
  2. 比べる相手の向きを `metric.comparator_directions` に事前登録する（plan.md §5）
  3. 状態の重心の行ファイルを この実験の work/ に向ける
  4. 出力先・実験 ID
  5. 判定規則 `decision` を特異性のものに差し替える

**比べる相手側の陽性対照 40 個**（plan.md §3.3）は `analysis/select_comparator_pc_010.py` が
あとから `controls.extra_fixed` に `role: comparator_positive_control` として足す。
このスクリプトは足さない（選び方を別のコードに分けて、事前登録の対象を明確にする）。

実行: .venv/bin/python experiments/EXP-20261005-001/analysis/make_config_010.py
"""
from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "experiments/EXP-20261004-001/config.yaml"
DST = ROOT / "experiments/EXP-20261005-001/config.yaml"
ED = "experiments/EXP-20261005-001"
SEEDS = [84, 85, 86]
SEC_PER_CELL_GENE = 0.1254   # spark-cf5b（NVIDIA GB10 / CUDA / bf16）の実測


def main() -> None:
    c = yaml.safe_load(SRC.read_text())

    c["experiment_id"] = "EXP-20261005-001"
    c["hypothesis_id"] = None
    c["kind"] = "specificity_of_goal_direction"
    c["note"] = (
        "EXP-20261004-001 で主判定を満たした候補 8 個について、**目標の向きへの動きが"
        "目標に特異的か、出発点から離れる動き全般か**を測る。摂動・細胞の選び方・遺伝子集合・"
        "主指標の形はすべて EXP-20261004-001 と同一で、変えるのは (1) seed、"
        "(2) 同じ形の指標を比べる相手の向きでも計算すること、の 2 点だけである。"
        "比べる相手の向きは、目標の向きと同じ作り方（深さ（入力の長さ）を 1 対 1 でそろえた"
        "細胞から作った重心の差）で作る。")

    c["output"] = {"experiment_dir": ED, "results_dir": "results"}
    c["data"]["tokenized_dataset"] = f"{ED}/work/tokenized.dataset"
    c["state_embs"]["rows_file"] = f"{ED}/work/state_rows-seed{{seed}}.json"
    c["robustness"]["seeds"] = SEEDS

    # ------------------------------------------------------------- 比べる相手の向き
    # すべて「深さ（入力の長さ）を 1 対 1 でそろえた細胞から作った重心の差」。
    # 重心の作り方の recipe は目標の向きとまったく同じで、ここが揃っていないと
    # 「入力の長さに由来する動き」と「生物としての類似」を分けられない（plan.md §5.1）。
    c["metric"]["comparator_directions"] = {
        "recipe": "unit(centroid_B) - unit(centroid_A) を単位長にする。"
                  "A・B は log(入力の長さ) の 1 対 1 最近傍マッチ（許容 0.10）で対応づけた細胞から作る",
        "goal": {"B": "aStr", "A": "start_for_aStr", "role": "主判定（再現）"},
        "comparators": {
            "dFIB":        {"B": "dFIB",          "A": "start_for_dFIB",  "role": "主な比べる相手（退化・皮質）"},
            "dMFIB":       {"B": "dMFIB",         "A": "start_for_dMFIB", "role": "記録（退化・髄質。提供者 14 人の偏りがある）"},
            "MFIB_vs_FIB": {"B": "MFIB_matched",  "A": "FIB_matched",     "role": "退化していない別の線維芽細胞の身元（皮質と髄質）"},
            "depth":       {"B": "start_shallow", "A": "start_deep",      "role": "入力の長さだけの軸"},
        },
        "record_only": {
            "aFIB": {"B": "aFIB_matched", "A": "start_for_aFIB", "role": "目標の片側（活性化線維芽細胞）"},
            "MYOF": {"B": "MYOF_matched", "A": "start_for_MYOF", "role": "目標の片側（筋線維芽細胞）"},
        },
        "orthogonalized": {
            "name": "spec",
            "against": ["dFIB", "dMFIB"],
            "note": "目標の向きから、退化した 2 状態の向きが張る面の成分を取り除いた残り。"
                    "副次。EXP-20261004-001 の重心での実測では目標の向きの長さの 91.9〜92.2% が残る",
        },
        "centroid_files": [
            f"{ED}/results/run-{{seed}}-L-1/state_emb_vectors.npz",
            f"{ED}/work/extra_centroids-seed{{seed}}.npz",
        ],
    }

    # ------------------------------------------------------------- 判定規則
    c["decision"] = {
        "単位": "個々の遺伝子（セット単位の検定は陽性対照にだけ当てる）",
        "規則の形": "3 seed すべてで標準得点 > 2.73 かつ 3 seed 平均 ≥ 3.6。"
                  "EXP-20261004-001 の V-1 とまったく同じ形で、閾値も同じ",
        "S-A 再現": {
            "統計量": "z_goal（目標の向きへの射影 ÷ 移動量 → 層別補正 → 対照 549 個で標準化）",
            "偶然の水準（EXP-20261004-001 の対照 549・seed 81/82/83 で実測）": 3 / 549,
        },
        "S-B 特異性（主判定）": {
            "統計量": "z_contrast(dFIB)。細胞ごとの (射影_goal − 射影_dFIB) ÷ 移動量 に"
                    "同じ層別補正を当て、対照 549 個で標準化したもの",
            "偶然の水準（同上で実測）": 1 / 549,
            "読み方": "S-A と S-B の両方を満たした遺伝子だけを「目標に特異的」と書く",
        },
        "S-C 記録": ["z_contrast(dMFIB)", "z_contrast(MFIB_vs_FIB)", "z_contrast(depth)",
                   "z_dFIB", "z_dMFIB", "z_depth", "z_spec", "z_aFIB", "z_MYOF"],
        "判定保留の範囲": {
            "条件": "3 seed すべてで > 2.73 だが 3 seed 平均 < 3.6",
            "根拠": "bf16 の丸めの誤差。z_goal で SD 0.243・観測最大 0.436（EXP-20261004-001 §5.6）。"
                  "差の統計量ではこれが約 0.28 になる（plan.md §6.2）",
        },
        "対象": {
            "MADCAM1": "ENSG00000099866", "CA12": "ENSG00000074410",
            "CXCL13": "ENSG00000156234", "CCL21": "ENSG00000137077",
            "ZFHX3": "ENSG00000140836", "SLC4A9": "ENSG00000113073",
            "GSN": "ENSG00000148180", "ELN": "ENSG00000049540",
        },
        "読む向き": "正の側を主判定にする。負の側は対照 549 個で目盛りが立っているかを"
                 "この実験で測り直してから、記録のみか解釈するかを決める（plan.md §6.6）",
    }

    # ------------------------------------------------------------- 計算量
    n_genes = 626   # 比べる相手側の陽性対照 40 個を足すと 666。select_comparator_pc_010.py が更新する
    n_cells = c["cell_selection"]["max_cells_per_gene"]
    c["compute"]["n_gene_runs"] = n_genes
    c["compute"]["n_cell_gene"] = n_genes * n_cells * len(SEEDS)
    c["compute"]["sec_per_cell_gene"] = SEC_PER_CELL_GENE

    c["merge"] = {
        "centroid_max_abs_diff": 0.0,
        "centroid_states_for_stop": ["FIB_MFIB", "aStr", "dFIB"],
        "注": "主判定 S-A・S-B が使うのは目標の向きと dFIB の向きである。"
             "止める条件はこの 3 本に当てる。dMFIB ほかは差を記録して続ける"
             "（EXP-20261004-001 §5.5.4 と同じ導き方）",
        "突き合わせ相手（記述のみ）": [
            "experiments/EXP-20261004-001/results/gene_seed_009.csv",
        ],
    }

    DST.write_text(yaml.safe_dump(c, allow_unicode=True, sort_keys=True))
    print(f"→ {DST}  遺伝子 {n_genes} 個 / seed {SEEDS} / "
          f"見込み {c['compute']['n_cell_gene'] * SEC_PER_CELL_GENE / 3600:.2f} GPU 時間（摂動のみ）")


if __name__ == "__main__":
    main()

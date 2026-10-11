#!/usr/bin/env python3
"""EXP-20261002-002 — 設定ファイルを 1 か所から機械的に作る（遺伝子リストを人が写さない）。

出力:
  config.yaml                      本体（3 seed 分。**事前登録の対象**）
  config-seed{81,82,83}.yaml       seed 単位で計算機を分けるとき（seeds と results_dir 以外は同一）
  work/config-smoke.yaml           試し実行専用（**事前登録の対象外**）
  work/config-machinecheck.yaml    機体の一致確認専用（**事前登録の対象外**）
  work/state_rows-seed{seed}.json  EXP-20261001-001 からの複製（sha256 で同一を確認）

## EXP-20261001-001 の設定から変えるのは次の 3 点だけである

| | 変えること | 変えない理由 / 変える理由 |
|---|---|---|
| 1 | **摂動する遺伝子を、追加の対照 500 個だけにする** | 既存の 717 個は摂動済み。回し直さない |
| 2 | 出力先を EXP-20261002-002 にする | 既存の結果を上書きしない |
| 3 | 実測の秒数・遺伝子数（見積もり用） | 見積もりの入力であって、計算の中身ではない |

**変えないもの（1 文字も）**: モデル・トークン化済み入力・状態の重心の行（rows_file）・
seed 81/82/83・細胞の選び方（200 細胞・提供者 20・下限 120）・摂動の型（過剰発現）・
emb_mode・主指標の形・交絡の補正（lp・ll の完全二次式）・一様成分の状態。
**既存の結果と突き合わせるので、ここが 1 つでも違うと合算できない。**
"""
from __future__ import annotations

import copy
import hashlib
import json
import shutil
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261002-002"
ED = "experiments/EXP-20261002-002"
P01D = "experiments/EXP-20261001-001"
SEEDS = [81, 82, 83]

# 盲検中（**解除されていない 3 個**）
BLOCKED = {"TP53": "ENSG00000141510", "SPP1": "ENSG00000118785", "LCN2": "ENSG00000148346"}

# 実測（EXP-20261001-001 の 717 遺伝子 × 200 細胞 × 3 seed から）。見積もりに使う
SEC_PER_CELL_GENE = 0.1154


def main() -> None:
    src = yaml.safe_load(open(ROOT / P01D / "config.yaml"))
    extra = json.loads((E / "work/extra_null_007.json").read_text())
    genes = [x["ensembl"] for x in extra["追加する対照"]]
    assert len(genes) == len(set(genes)) == extra["n_extra"], "追加の対照に重複がある"
    assert not (set(genes) & set(BLOCKED.values())), "盲検中の遺伝子が入っている"

    # ---- 状態の重心の行は EXP-20261001-001 と**同一のファイル**を使う（複製して sha256 照合）
    rows_hash = {}
    for s in SEEDS:
        a = ROOT / P01D / f"work/state_rows-seed{s}.json"
        b = E / f"work/state_rows-seed{s}.json"
        shutil.copyfile(a, b)
        ha = hashlib.sha256(a.read_bytes()).hexdigest()
        hb = hashlib.sha256(b.read_bytes()).hexdigest()
        assert ha == hb, f"state_rows-seed{s}.json の複製が一致しない"
        rows_hash[s] = ha

    cfg = {
        "experiment_id": "EXP-20261002-002",
        "hypothesis_id": "HYP-20261001-01",     # Status は変えない（関門は PI が判断する）
        "kind": "extra_controls_for_bias_budget",
        "note": ("追加の対照 500 個だけを摂動する。既存の 717 個（EXP-20261001-001）は摂動済みで"
                 "回し直さない。**合算して対照 1,000 個で健全性の検査 S-1〜S-4 をやり直す**"),
        # ---- 以下 model / data / states / state_embs / cell_selection は EXP-20261001-001 と同一
        "model": copy.deepcopy(src["model"]),
        "data": copy.deepcopy(src["data"]),
        "states": copy.deepcopy(src["states"]),
        "state_embs": {"max_cells_per_state": src["state_embs"]["max_cells_per_state"],
                       "max_cells_per_donor": src["state_embs"]["max_cells_per_donor"],
                       "rows_file": f"{ED}/work/state_rows-seed{{seed}}.json"},
        "cell_selection": copy.deepcopy(src["cell_selection"]),
        "perturbation": copy.deepcopy(src["perturbation"]),
        "blocked_genes": sorted(BLOCKED.values()),
        # ---- 対照は 1 つも渡さない。**摂動するのは追加の対照 500 個だけ**
        #      （run_isp.py の見積もりが二重に数えないよう、controls 側は空にする）
        "controls": {"positive": [], "negative_fixed": [], "negative_random_n": 0,
                     "extra_fixed": [], "comparator": None},
        "capture": copy.deepcopy(src["capture"]),
        "metric": copy.deepcopy(src["metric"]),
        "stats": copy.deepcopy(src["stats"]),
        "robustness": copy.deepcopy(src["robustness"]),
        "output": {"results_dir": "results", "experiment_dir": ED},
        "compute": {"device": "mps", "dtype": "bf16", "cache_model": True,
                    "sec_per_cell_gene": SEC_PER_CELL_GENE,
                    "n_gene_runs": len(genes),
                    "n_cell_gene": len(genes) * src["cell_selection"]["max_cells_per_gene"]
                    * len(SEEDS)},
    }
    # **トークン化済み入力は EXP-20261001-001 のものをそのまま使う**（読むだけ）。
    # 作り直すと行の並びが変わる危険があり、既存の結果と行番号で突き合わせられなくなる。
    cfg["data"]["tokenized_dataset"] = f"{P01D}/work/tokenized.dataset"
    cfg["perturbation"]["genes_to_perturb"] = genes

    # ---- 合算と検査のための設定（評価スクリプトが読む）
    cfg["merge"] = {
        "existing": {"experiment": "EXP-20261001-001",
                     "results_dir": f"{P01D}/results",
                     "gene_table": f"{P01D}/results/eval/gene_seed_005.csv",
                     "genes": f"{P01D}/work/genes_005.json",
                     "subsets": f"{P01D}/work/subsets_005.json",
                     "null_pool": f"{P01D}/work/null_pool_005.json",
                     "power": f"{P01D}/work/power_005.json",
                     "criteria": f"{P01D}/results/eval/criteria_005.json"},
        "extra_null": f"{ED}/work/extra_null_007.json",
        # **状態の重心が既存と同一であることを確かめる**（違えば合算してはいけない）
        "centroid_max_abs_diff": 0.0,
        # 既存の 717 遺伝子の値を、対照 500 個で測り直したときに再現できるか
        "reproduction_tolerance": 0.002,
        "state_rows_sha256": {str(k): v for k, v in rows_hash.items()},
    }
    # ---- 健全性の検査。**EXP-20261002-001 で固定した S-1〜S-4 をそのまま使う**
    bb = yaml.safe_load(open(ROOT / "experiments/EXP-20261002-001/config.yaml"))["bias_budget"]
    cfg["bias_budget"] = copy.deepcopy(bb)      # **予算 0.10 を動かさない**

    (E / "config.yaml").write_text(
        "# EXP-20261002-002 — すべてのパラメータをここに書く（plan.md §6 が事前登録の本体）\n"
        "# 作成: analysis/make_configs_007.py（手で書き換えない）\n"
        "# **摂動するのは追加の対照 500 個だけ。** 既存の 717 個は EXP-20261001-001 で摂動済み。\n"
        "# 健全性の検査 S-1〜S-4 と予算 0.10 は EXP-20261002-001 で固定したものをそのまま使う。\n"
        + yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True))

    # ---- seed ごと（**分けるのは seed 単位のみ。遺伝子で分けてはいけない**）
    for s in SEEDS:
        c = copy.deepcopy(cfg)
        c["robustness"]["seeds"] = [s]
        c["output"]["results_dir"] = f"results-seed{s}"
        c["compute"]["n_cell_gene"] = len(genes) * cfg["cell_selection"]["max_cells_per_gene"]
        (E / f"config-seed{s}.yaml").write_text(
            f"# EXP-20261002-002 seed {s} 専用（config.yaml から機械的に作る。手で書き換えない）\n"
            "# **seeds と results_dir 以外は config.yaml と同一**\n"
            + yaml.safe_dump(c, sort_keys=False, allow_unicode=True))

    # ---- 試し実行（**事前登録の対象外**。数値は意味を持たない。コードが通るかだけを見る）
    sm = copy.deepcopy(cfg)
    sm["perturbation"]["genes_to_perturb"] = genes[:6]
    sm["robustness"]["seeds"] = [SEEDS[0]]
    sm["cell_selection"] = dict(sm["cell_selection"], max_cells_per_gene=20,
                                max_cells_per_donor=5, min_cells_per_gene=5)
    sm["state_embs"] = {"max_cells_per_state": 40, "max_cells_per_donor": 5}  # rows_file を使わない
    sm["output"]["results_dir"] = "results-smoke"
    sm["compute"] = dict(sm["compute"], n_gene_runs=6, n_cell_gene=6 * 20)
    (E / "work/config-smoke.yaml").write_text(
        "# **試し実行専用。事前登録の対象外。** 状態の重心を 40 細胞で作るので\n"
        "# 数値は本番と比べられない。**評価コードが最後まで通るか**だけを見る。\n"
        + yaml.safe_dump(sm, sort_keys=False, allow_unicode=True))

    # ---- 機体の一致確認（**事前登録の対象外**）
    #   既存の対照 500 個から 12 個を選び、**EXP-20261001-001 と同じ設定で回し直して、
    #   保存済みの出力と突き合わせる**。機体どうしを比べるだけでなく、
    #   「**この機体が、既存の結果を作った機体と同じ値を出すか**」を直接見る。
    #   これが合わなければ、新しい対照を足しても既存の対照と比べられない。
    null500 = [x["ensembl"] for x in
               json.loads((ROOT / P01D / "work/null_pool_005.json").read_text())["versions"]["B_500"]]
    mcg = null500[:12]
    mc = copy.deepcopy(cfg)
    mc["perturbation"]["genes_to_perturb"] = mcg
    mc["robustness"]["seeds"] = [SEEDS[0]]
    mc["output"]["results_dir"] = "results-machinecheck"
    mc["compute"] = dict(mc["compute"], n_gene_runs=len(mcg),
                         n_cell_gene=len(mcg) * cfg["cell_selection"]["max_cells_per_gene"])
    mc["machine_check"] = {
        "対象": "既存の対照 500 個の先頭 12 個（**EXP-20261001-001 で摂動済み**）",
        "突き合わせる相手": f"{P01D}/results/run-81-L-1/genes/<gene>/emb/{{orig,pert}}.npy",
        "注": "機体どうしの比較ではなく、**既存の結果との比較**である。"
              "追加の対照を既存の対照と同じ尺度で測れることが要点",
    }
    (E / "work/config-machinecheck.yaml").write_text(
        "# **機体の一致確認専用。事前登録の対象外。**\n"
        "# 既存の対照 12 個を回し直し、EXP-20261001-001 の保存済みの出力と突き合わせる。\n"
        + yaml.safe_dump(mc, sort_keys=False, allow_unicode=True))

    print(f"摂動する遺伝子: {len(genes)}（追加の対照のみ）")
    print(f"細胞・遺伝子の延べ: {cfg['compute']['n_cell_gene']:,}"
          f" = {len(genes)} × {cfg['cell_selection']['max_cells_per_gene']} × {len(SEEDS)} seed")
    print("state_rows の sha256（EXP-20261001-001 と同一）:",
          {k: v[:12] for k, v in rows_hash.items()})
    print("書いた設定:", [p.name for p in sorted(E.glob("config*.yaml"))],
          [p.name for p in sorted((E / "work").glob("config-*.yaml"))])


if __name__ == "__main__":
    main()

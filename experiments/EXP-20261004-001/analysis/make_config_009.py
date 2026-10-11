#!/usr/bin/env python3
"""EXP-20261004-001 の config.yaml を作る（**手で書き換えない**）。

土台は `experiments/EXP-20261003-001/config.yaml`。変えた箇所だけをこのスクリプトが一覧で出す。
遺伝子の一覧は `analysis/precheck_009.py` が書いた `work/genes_009.json` から読む
（guardrails 2026-09-27: 設定に書く遺伝子リストは実行スクリプトが読める形で設定に直接書く）。

使い方: `.venv/bin/python experiments/EXP-20261004-001/analysis/make_config_009.py`
"""
from __future__ import annotations

import json
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261004-001"
BASE = ROOT / "experiments/EXP-20261003-001/config.yaml"
SRC = ROOT / "experiments/EXP-20261001-001/config.yaml"

# **bf16 を使う理由は plan.md §5.6 に書く。** ここでは値だけを置く
DTYPE = "bf16"
SEC_PER_CELL_GENE = 0.1254   # cf5b（NVIDIA GB10 / CUDA / bf16）の実測。
# 出どころ: experiments/EXP-20261003-001/results-machinecheck-cf5b/timing.jsonl
# 12 遺伝子 × 200 細胞 = 2,400 細胞遺伝子を 301 秒 → 0.1254 s。
# EXP-20261001-001 が使っていた 0.1154 より**大きいほうを採る**（guardrails 2026-09-27）。


def main() -> None:
    cfg = yaml.safe_load(BASE.read_text())
    src = yaml.safe_load(SRC.read_text())
    g = json.loads((E / "work/genes_009.json").read_text())
    changes: list[str] = []

    def setv(path: str, value):
        *parents, key = path.split(".")
        cur = cfg
        for k in parents:
            cur = cur.setdefault(k, {})
        old = cur.get(key)
        cur[key] = value
        shown = (f"{len(value)} 個" if isinstance(value, list) else value)
        shown_old = (f"{len(old)} 個" if isinstance(old, list) else old)
        changes.append(f"{path}: {shown_old} → {shown}")

    cand = list(g["候補"].values())
    pos = g["陽性対照71"]
    nulls = g["主対照"] + g["低発現の追加対照"]

    setv("experiment_id", "EXP-20261004-001")
    setv("hypothesis_id", None)   # 方法の実験。候補の確認は新しい仮説（G1）が要る。plan.md §0
    setv("kind", "common_cell_set_confirmation")
    setv("output.experiment_dir", "experiments/EXP-20261004-001")
    setv("data.tokenized_dataset", "experiments/EXP-20261004-001/work/tokenized.dataset")
    setv("state_embs.rows_file", "experiments/EXP-20261004-001/work/state_rows-seed{seed}.json")

    # ---- この実験の中心: 共通の細胞集合
    setv("cell_selection.mode", "common")
    setv("cell_selection.min_cells_per_gene", 200)   # 共通集合なので必ず 200。届かなければ止める

    # ---- 測る遺伝子
    setv("perturbation.genes_to_perturb", cand)
    setv("controls.positive", pos)
    setv("controls.negative_fixed", src["controls"]["negative_fixed"])
    setv("controls.fixed_role", src["controls"]["fixed_role"])
    setv("controls.extra_fixed", src["controls"]["extra_fixed"])
    setv("controls.negative_random_n", len(nulls))
    setv("controls.negative_sampler", {"method": "explicit", "genes": nulls})

    # ---- 計算量
    setv("compute.device", "cuda")
    setv("compute.dtype", DTYPE)
    setv("compute.sec_per_cell_gene", SEC_PER_CELL_GENE)
    n_genes = len(set(cand) | set(pos) | set(src["controls"]["negative_fixed"]) |
                  {x for grp in src["controls"]["extra_fixed"] for x in grp["genes"]} | set(nulls))
    setv("compute.n_gene_runs", n_genes)
    setv("compute.n_cell_gene", n_genes * 200 * 3)

    # ---- 既存の結果との突き合わせ（**合算はしない。別の尺度で測るため**）
    cfg["merge"] = {
        "注": "この実験は共通の細胞集合で測り直すので、**既存の 1,217 遺伝子とは合算しない**。"
              "既存は遺伝子ごとに細胞の顔ぶれが違う測り方であり、同じ尺度ではない",
        "突き合わせ相手（記述のみ）": {
            "既存の z（遺伝子ごと・3 seed）": [
                "experiments/EXP-20261002-002/results/gene_seed_007-{81,82,83}.csv",
                "experiments/EXP-20261003-001/results/gene_seed_008-{81,82,83}.csv",
            ],
        },
        "state_rows_sha256": cfg["merge"]["state_rows_sha256"],
        "centroid_states_for_stop": ["FIB_MFIB", "aStr"],
        "centroid_max_abs_diff": 0.0,
    }

    # ---- 判定（plan.md §6 が本体。ここは機械が読む形）
    cfg["decision"] = {
        "単位": "個々の遺伝子（セット単位の検定は副次）",
        "V-1": {"条件": "3 seed すべてで z > 2.73", "閾値": 2.73,
                "出どころ": "単一遺伝子の検出下限（検出力 80%）。EXP-20260928-002 の実測"
                            "（hypotheses/HYP-20260928-01.md, experiments/EXP-20260928-001/work/power_001.json）"},
        "V-1b 判定保留の帯": {
            "条件": "3 seed すべてで z > 2.73 だが 3 seed の平均 z < 3.6",
            "閾値": 3.6,
            "読み方": "**合否を論じない（判定保留）。** bf16 の遺伝子ごとの誤差（平均 0.30 z・"
                      "観測最大 0.64 z、work/dtype_resolution_009.json）が 2.73 の線をまたぐ帯であり、"
                      "seed を重ねても消えない（bf16 の誤差は遺伝子ごとに再現する）",
            "出どころ": "PI 指示 2026-10-04。実測の裏付けは plan.md §5.6",
        },
        "V-2 偶然の水準": {
            "測り方": "**この実験の対照（共通細胞集合で測り直したもの）**に V-1 を当てて実測する。"
                      "既存の 0.9% は別の測り方の値なので使わない",
            "報告": "期待件数 = 実測の片側割合 × 候補の個数 と、実際の件数を必ず並べる",
        },
        "読む向き": "**正の側だけ。** 負の側に目盛りが立っているかを、この実験の対照で測り直してから決める"
                    "（既存の測り方では 1,000 個中 1 個しか z < −2.73 が無かった）",
        "対象": g["候補"],
    }
    # ---- 着手前チェック（CLAUDE.md §6）。目標・出発点・陽性対照は EXP-20261001-001 と同じなので
    #      その節をそのまま持ってくる（experiment_id だけ変える）
    cfg["axis_feasibility"] = dict(src["axis_feasibility"])
    cfg["axis_feasibility"]["experiment_id"] = "EXP-20261004-001-precheck"
    changes.append("axis_feasibility: EXP-20261001-001 の節をそのまま（experiment_id だけ変更）")

    cfg["note"] = (
        "**この実験の主目的は、遺伝子ごとに摂動する細胞の顔ぶれが違うという弱点を直すこと。** "
        "cell_selection.mode: common により、出発点 FIB_MFIB から seed ごとに 1 つだけ作った "
        "200 細胞の集合で、候補・陽性対照・対照・恒常発現のすべてを摂動する。"
        "Geneformer 本体は未発現の遺伝子の過剰発現を支援しており（perturber_utils.overexpress_tokens、"
        "in_silico_perturber.py:525 で overexpress は細胞を絞らない）、絞っていたのは run_isp.py 側の選択だった")

    E.joinpath("config.yaml").write_text(
        yaml.safe_dump(cfg, allow_unicode=True, sort_keys=True, width=100))
    print(f"土台: {BASE.relative_to(ROOT)}")
    print(f"遺伝子の一覧: {(E / 'work/genes_009.json').relative_to(ROOT)}")
    print("変えた箇所:")
    for c in changes:
        print("  -", c)
    print(f"  - merge / decision / note を書き換え")
    print(f"測る遺伝子 {n_genes} 個 × 200 細胞 × 3 seed = {n_genes * 600:,} 細胞遺伝子")
    print(f"見積もり（摂動のみ）: {n_genes * 600 * SEC_PER_CELL_GENE / 3600:.2f} GPU 時間")


if __name__ == "__main__":
    main()

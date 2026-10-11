#!/usr/bin/env python3
"""EXP-20261003-001 — 設定を作る（手で書き換えない）。

EXP-20261002-002 の `config.yaml` を土台にし、**変える箇所だけ**を変える。
変えた箇所は最後に一覧で出力し、plan.md §5 の表と一致することを人が確かめる。

作るもの:
  config.yaml               本番（6 遺伝子 × 200 細胞 × seed 81/82/83）
  work/config-machinecheck.yaml  機体の一致確認（既存の対照 12 個を Mac で回し直す）
  work/config-smoke.yaml    試し実行（評価コードを最後まで通すため。判定ではない）
  work/state_rows-seed{81,82,83}.json  EXP-20261001-001 から複製（sha256 で照合）

実行: .venv/bin/python experiments/EXP-20261003-001/analysis/make_configs_008.py
"""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261003-001"
P01 = ROOT / "experiments/EXP-20261001-001"
P02 = ROOT / "experiments/EXP-20261002-002"

SEEDS = [81, 82, 83]
# 着手前チェック（work/precheck_008.json）に 6 個とも合格した。順は plan.md §2 の表と同じ
TARGETS = [
    ("ACTA2", "ENSG00000107796"),
    ("TAGLN", "ENSG00000149591"),
    ("CNN1", "ENSG00000130176"),
    ("MYH11", "ENSG00000133392"),
    ("POSTN", "ENSG00000133110"),
    ("PDGFRB", "ENSG00000113721"),
]
# 機体の一致確認に使う既存の対照 12 個（EXP-20261002-002 の確認と同じ顔ぶれ）
MACHINE_CHECK_GENES = [
    "ENSG00000149792", "ENSG00000104412", "ENSG00000104325", "ENSG00000148730",
    "ENSG00000182117", "ENSG00000205045", "ENSG00000071909", "ENSG00000075292",
    "ENSG00000148835", "ENSG00000037749", "ENSG00000112584", "ENSG00000176261",
]
NOTE = ("**摂動するのは新しい 6 個だけ。** 既存の 1,217 遺伝子（介入 119・陽性対照 71・対照 1,000 ほか）は "
        "EXP-20261001-001 / EXP-20261002-002 で摂動済みで、回し直さない。"
        "目的は「目標の軸が筋線維芽細胞らしさを測っているか」を、"
        "筋線維芽細胞を定義する遺伝子そのものを過剰発現して確かめること")


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> None:
    base = yaml.safe_load((P02 / "config.yaml").read_text())
    changes: list[tuple[str, str, str]] = []

    def setv(cfg: dict, path: str, val, why: str) -> None:
        keys = path.split(".")
        d = cfg
        for k in keys[:-1]:
            d = d[k]
        old = d.get(keys[-1], "（無し）")
        olds = (f"{len(old)} 個のリスト" if isinstance(old, list) else str(old))[:70]
        news = (f"{len(val)} 個のリスト" if isinstance(val, list) else str(val))[:70]
        d[keys[-1]] = val
        changes.append((path, f"{olds} → {news}", why))

    # ---------------------------------------------------------------- rows ファイルを複製
    (E / "work").mkdir(parents=True, exist_ok=True)
    rows_sha = {}
    for s in SEEDS:
        src = P01 / f"work/state_rows-seed{s}.json"
        dst = E / f"work/state_rows-seed{s}.json"
        shutil.copy2(src, dst)
        assert sha256(src) == sha256(dst)
        rows_sha[str(s)] = sha256(dst)
        assert rows_sha[str(s)] == base["merge"]["state_rows_sha256"][str(s)], \
            f"seed {s} の rows ファイルが EXP-20261002-002 の記録と違う"

    # ---------------------------------------------------------------- 本番の設定
    cfg = yaml.safe_load(yaml.safe_dump(base))
    setv(cfg, "experiment_id", "EXP-20261003-001", "実験 ID")
    setv(cfg, "kind", "myofibroblast_marker_axis_validation", "実験の種類")
    setv(cfg, "note", NOTE, "何を回すか")
    setv(cfg, "perturbation.genes_to_perturb", [e for _, e in TARGETS],
         "**ここが本体**: 筋線維芽細胞・平滑筋・活性化線維芽細胞・間質のマーカー 6 個")
    setv(cfg, "state_embs.rows_file",
         "experiments/EXP-20261003-001/work/state_rows-seed{seed}.json",
         "複製した rows ファイル（中身は EXP-20261001-001 と同一。sha256 照合済み）")
    setv(cfg, "output.experiment_dir", "experiments/EXP-20261003-001", "出力先")
    setv(cfg, "compute.n_gene_runs", len(TARGETS), "見積もり用（遺伝子数）")
    setv(cfg, "compute.n_cell_gene", len(TARGETS) * 200 * len(SEEDS), "見積もり用（細胞 × 遺伝子）")
    cfg["merge"] = {
        "existing_roots": [
            {"label": "既存(EXP-20261001-001)", "results_dir": "experiments/EXP-20261001-001/results",
             "n_genes": 717},
            {"label": "追加の対照(EXP-20261002-002)", "results_dir": "experiments/EXP-20261002-002/results",
             "n_genes": 500},
        ],
        "genes": "experiments/EXP-20261001-001/work/genes_005.json",
        "subsets": "experiments/EXP-20261001-001/work/subsets_005.json",
        "null_pool": "experiments/EXP-20261001-001/work/null_pool_005.json",
        "extra_null": "experiments/EXP-20261002-002/work/extra_null_007.json",
        "criteria": "experiments/EXP-20261002-002/results/criteria_007.json",
        "centroid_max_abs_diff": 0.0,
        "reproduction_tolerance": 0.002,
        "state_rows_sha256": rows_sha,
    }
    changes.append(("merge", "EXP-20261001-001 のみ → 2 つの結果（717 + 500）を読む",
                    "対照 1,000 個を同じ尺度で再現するため"))
    cfg["decision"] = {
        "単位": "個々の遺伝子（**セット単位の検定はしない**）",
        "V-1": {"条件": "3 seed すべてで z > 2.73", "閾値": 2.73,
                "出どころ": "単一遺伝子の検出下限。EXP-20260928-002 の実測（hypotheses/HYP-20260928-01.md）"},
        "V-2": {"偶然の水準 p": 0.009, "6 個での期待件数": 0.054,
                "P(1 件以上|偶然)": 0.0528,
                "出どころ": "experiments/EXP-20261003-001/work/precheck_008.json の C-4"
                            "（対照 1,000 個の確定済みの z から実測）"},
        "読む向き": "**正の側だけ。** 対照の z 分布は負の裾が薄く（3 seed すべてで z < -2.73 の対照は "
                    "1,000 個中 1 個）校正されていないため、負の値は解釈しない",
        "対象": {s: e for s, e in TARGETS},
    }
    changes.append(("decision", "（無し） → V-1 / V-2", "この実験の判定基準（結果を見る前に固定）"))
    cfg["_note_header"] = None
    del cfg["_note_header"]
    header = (f"# EXP-20261003-001 — すべてのパラメータをここに書く（plan.md §6 が事前登録の本体）\n"
              f"# 作成: analysis/make_configs_008.py（手で書き換えない）\n"
              f"# 土台: experiments/EXP-20261002-002/config.yaml。変えた箇所は plan.md §5 の表のとおり\n")
    (E / "config.yaml").write_text(header + yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False))

    # ---------------------------------------------------------------- 機体の一致確認
    mc = yaml.safe_load(yaml.safe_dump(cfg))
    mc["note"] = ("**機体の一致確認専用。事前登録の対象外。** 既存の対照 12 個をこの機体で回し直し、"
                  "EXP-20261001-001 の保存済みの出力と突き合わせる")
    mc["perturbation"]["genes_to_perturb"] = MACHINE_CHECK_GENES
    mc["robustness"]["seeds"] = [81]
    mc["output"]["results_dir"] = "results-machinecheck"
    mc["compute"]["n_gene_runs"] = len(MACHINE_CHECK_GENES)
    mc["compute"]["n_cell_gene"] = len(MACHINE_CHECK_GENES) * 200
    mc.pop("decision", None)
    (E / "work/config-machinecheck.yaml").write_text(
        "# 機体の一致確認専用。事前登録の対象外。\n"
        + yaml.safe_dump(mc, allow_unicode=True, sort_keys=False))

    # ---------------------------------------------------------------- 試し実行
    sm = yaml.safe_load(yaml.safe_dump(cfg))
    sm["note"] = ("**試し実行専用。事前登録の対象外。判定ではない。** "
                  "目的は評価コードが最後まで通ることの確認だけ")
    sm["robustness"]["seeds"] = [81]
    sm["cell_selection"]["max_cells_per_gene"] = 20
    sm["cell_selection"]["min_cells_per_gene"] = 10
    sm["state_embs"]["rows_file"] = None
    sm["state_embs"]["max_cells_per_state"] = 40
    sm["state_embs"]["max_cells_per_donor"] = 5
    sm["output"]["results_dir"] = "results-smoke"
    sm["compute"]["n_gene_runs"] = len(TARGETS)
    sm["compute"]["n_cell_gene"] = len(TARGETS) * 20
    sm.pop("decision", None)
    (E / "work/config-smoke.yaml").write_text(
        "# 試し実行専用。事前登録の対象外。**判定ではない。**\n"
        + yaml.safe_dump(sm, allow_unicode=True, sort_keys=False))

    print("変えた箇所（土台: experiments/EXP-20261002-002/config.yaml）")
    for p, d, why in changes:
        print(f"  - {p}: {d}\n      理由: {why}")
    print("\n摂動する 6 個:")
    for s, e in TARGETS:
        print(f"  {s:8s} {e}")
    print("\nrows ファイルの sha256:", json.dumps(rows_sha, indent=1))
    print("書いたもの:", [str(p.relative_to(ROOT)) for p in
                          (E / "config.yaml", E / "work/config-machinecheck.yaml",
                           E / "work/config-smoke.yaml")])


if __name__ == "__main__":
    main()

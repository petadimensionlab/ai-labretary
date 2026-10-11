#!/usr/bin/env python3
"""EXP-20261003-001 — seed 83 を既存と同じ機体（spark-c083）で回し直すための設定を作る。

## なぜ作るか（2026-10-03。**判定基準は 1 文字も変えない**）

評価が plan.md §6.1 の「状態の重心の検査」で止まった（中断条件 A-2）。
seed 83 の重心 dMFIB だけが EXP-20261001-001 と一致しなかった
（最大差 0.0009765625 = bf16 の 1 刻み、1,152 成分中 72 個、コサイン 0.99999996）。
原因は機体の違いで、**既存の seed 83 は spark-c083、今回は spark-cf5b** で回したためである。
dMFIB は 137 細胞と最も小さい状態で、そこでだけ丸めの差が出た。

→ **Labretary の決定（2026-10-03、案 A）: seed 83 を spark-c083 で回し直す。**
  事前登録を 1 文字も変えずに通る見込みである。
  **seed 81・82（spark-cf5b、照合済み）は回し直さない。**

**この設定は事前登録の対象外である**（`work/config-machinecheck.yaml` と同じ扱い）。
固定済みの `config.yaml` から**機械的に**作り、変えた箇所は 2 つだけであることを出力で示す。

実行: .venv/bin/python experiments/EXP-20261003-001/analysis/make_config_seed83_008.py
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261003-001"

# 固定済みの config.yaml の sha256（work/SHA256SUMS-preregistered.txt と一致すること）
FROZEN = "1357a0eede61629081d0e7c5a04b4c6ca6297a369de7b48babfc0c8f0b3a7a21"


def main() -> None:
    src = E / "config.yaml"
    got = hashlib.sha256(src.read_bytes()).hexdigest()
    if got != FROZEN:
        raise SystemExit(f"**止める**: config.yaml が固定したものと違う\n  固定 {FROZEN}\n  いま {got}")
    print(f"土台の config.yaml は固定したものと一致: {got}")

    cfg = yaml.safe_load(src.read_text())
    before = {"robustness.seeds": list(cfg["robustness"]["seeds"]),
              "output.results_dir": cfg["output"]["results_dir"]}
    cfg["robustness"]["seeds"] = [83]
    cfg["output"]["results_dir"] = "results-seed83-c083"
    cfg["compute"]["n_cell_gene"] = len(cfg["perturbation"]["genes_to_perturb"]) * 200
    after = {"robustness.seeds": cfg["robustness"]["seeds"],
             "output.results_dir": cfg["output"]["results_dir"]}

    # 変わっていないことを機械的に確かめる（**判定に関わるものは 1 つも触らない**）
    base = yaml.safe_load(src.read_text())
    must_match = ["experiment_id", "model", "data", "states", "state_embs", "cell_selection",
                  "perturbation", "blocked_genes", "controls", "capture", "metric", "stats",
                  "merge", "decision"]
    diff = [k for k in must_match if cfg.get(k) != base.get(k)]
    if diff:
        raise SystemExit(f"**止める**: 変えてはいけない節が変わっている: {diff}")

    out = E / "work/config-seed83-c083.yaml"
    out.write_text(
        "# **seed 83 を既存と同じ機体 spark-c083 で回し直すための設定。事前登録の対象外。**\n"
        "# 作成: analysis/make_config_seed83_008.py（手で書き換えない）\n"
        "# 土台: 固定済みの config.yaml（sha256 照合済み）。変えたのは seeds と results_dir だけ。\n"
        "# 判定基準（decision）・遺伝子・細胞の選び方・指標・対照は 1 文字も変えていない。\n"
        + yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False))

    print("\n変えた箇所（この 2 つ + 見積もり用の細胞数だけ）")
    for k in before:
        print(f"  - {k}: {before[k]} → {after[k]}")
    print(f"  - compute.n_cell_gene: {base['compute']['n_cell_gene']} → {cfg['compute']['n_cell_gene']}（見積もり用）")
    print(f"\n判定に関わる節（{len(must_match)} 個）はすべて土台と同一であることを確認した")
    print("書いたもの:", out.relative_to(ROOT))


if __name__ == "__main__":
    main()

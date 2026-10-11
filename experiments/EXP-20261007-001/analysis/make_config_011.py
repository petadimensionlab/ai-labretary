#!/usr/bin/env python3
"""EXP-20261007-001 の config.yaml（本番）と work/config-smoke.yaml（試し実行）を、
work/genes_011.json・work/donor_groups_011.json から機械的に作る。

遺伝子集合・群・判定の数値は**このファイルと plan.md にだけ書く**。手で config.yaml を直さない。

走らせ方:
    .venv/bin/python experiments/EXP-20261007-001/analysis/make_config_011.py
"""
from __future__ import annotations

import json
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261007-001"
W = E / "work"
REL = "experiments/EXP-20261007-001"


def main() -> None:
    g = json.loads((W / "genes_011.json").read_text())
    dg = json.loads((W / "donor_groups_011.json").read_text())["groups"]
    af = yaml.safe_load((ROOT / "outputs/design/2026-10-07-A2-feasibility/axis_feasibility_A2.yaml").read_text())
    cand = g["候補"]
    pos = list(g["陽性対照71"].values())
    hk = list(g["恒常発現3"].values())
    neu = list(g["中立7"].values())
    ctrl = [r["ensembl"] for r in g["対照500"]]
    n_genes = len(set(cand.values()) | set(pos) | set(hk) | set(neu) | set(ctrl))
    # 本番は float32 で、判定を保留する範囲を置かない（PI 判断 2026-10-08。plan.md §5.3）
    assert n_genes == 587, n_genes

    cfg = {
        "experiment_id": "EXP-20261007-001",
        "hypothesis_id": "HYP-20261007-01",
        "kind": "replication_in_single_cell_rnaseq",
        "note": ("HYP-20261007-01 の確認的検定。A2（KPMP 単一細胞 RNA-seq）の出発点 FIB を、提供者で分けた"
                 "互いに重ならない 3 群に分け（cell_selection.donor_groups。seed 81/82/83 は群の名前）、"
                 "群ごとに遺伝子によらず同一の細胞集合で候補 8・陽性対照 71・対照 500・恒常発現 3・中立 7 を過剰発現する。"
                 "状態の重心は 3 群で同一（state_rows の内容が同じ）。"),
        "data": {
            "input_h5ad": f"{REL}/work/input.h5ad",
            "tokenized_dataset": f"{REL}/work/tokenized.dataset",
            "cell_type_col": "author_cell_type",
            "state_key": "isp_state",
            "custom_attr_name_dict": {"isp_state": "isp_state", "donor_id": "donor_id",
                                      "author_cell_type": "author_cell_type", "disease": "disease",
                                      "n_counts": "n_counts"},
            "filter_data": None, "max_ncells": None},
        "model": {"version": "V2", "directory": "models/Geneformer/Geneformer-V2-316M",
                  "type": "Pretrained", "num_classes": 0, "emb_layer": -1},
        "states": {"start_state": "FIB", "goal_state": "aStr", "alt_states": ["dFIB"]},
        "perturbation": {"genes_to_perturb": list(cand.values()), "perturb_type": "overexpress",
                         "perturb_rank_shift": None, "combos": 0, "anchor_gene": None,
                         "emb_mode": "cls", "forward_batch_size": 8, "nproc": 1, "controls_only": False},
        "controls": {"positive": pos, "negative_fixed": hk, "fixed_role": "housekeeping",
                     "extra_fixed": [{"role": "negctrl_panel", "genes": neu}],
                     "negative_random_n": len(ctrl),
                     "negative_sampler": {"method": "explicit", "genes": ctrl},
                     "comparator": None},
        "blocked_genes": g["盲検（除外にだけ使う）"],
        "cell_selection": {"mode": "common", "donor_key": "donor_id", "max_cells_per_donor": 20,
                           "max_cells_per_gene": 60, "min_cells_per_gene": 50,
                           "donor_groups": {k: v for k, v in dg.items()}},
        "state_embs": {"rows_file": f"{REL}/work/state_rows-seed{{seed}}.json",
                       "max_cells_per_state": 460, "max_cells_per_donor": 20},
        "robustness": {"seeds": [81, 82, 83], "emb_layers": [-1], "subsample_frac": 1.0},
        "stats": {"mode": "goal_state_shift", "fdr_threshold": 0.05},
        "compute": {"device": "cuda", "dtype": "fp32", "cache_model": True,
                    # float32 の単価: bfloat16 の 0.1254 秒に、A2 で実測した所要の比 2.55（14.8 分 ÷ 5.8 分、plan.md §5.3.1）を掛けた
                    "sec_per_cell_gene": 0.32, "n_gene_runs": n_genes,
                    "n_cell_gene": n_genes * 173},
        "capture": {"embeddings": True, "classifier_logits": False},
        "metric": {"primary_form": "projection_on_start_to_goal_direction_over_movement",
                   "projection_states": {"from": "start", "to": "goal"},
                   "confound_correction": {"form": "stratified_by_expression", "variables": ["lp", "ll"],
                                           "drop_zero_movement": True},
                   "uniform_component_states": ["dFIB"]},
        "decision": {
            "単位": "群（提供者で分けた互いに重ならない 3 群。seed 81/82/83 は群の名前）",
            "標準化": "群ごとに、対照 510 個（対照 500 + 恒常発現 3 + 中立 7）の平均と標準偏差で D_raw を標準得点にする",
            "P-1 主判定（セット）": {
                "対象": "候補 8 個のセット平均の標準得点", "偽セット": "対照 510 個から復元ありで 8 個 × 1,000 組",
                "条件(i)": "3 群すべてで偽セットの 99 パーセンタイルを超える",
                "条件(ii)": "3 群すべてで 1.59 以上（A2 の対照で実測した最小検出可能効果量が 2.34 を超えたら、その実測値以上）",
                "閾値(ii)": 1.59, "置き換えの基準": 2.34},
            "P-2 副次（遺伝子ごと）": {
                "条件": "3 群すべてで標準得点 > 2.73（判定を保留する範囲は置かない）",
                "閾値": 2.73, "置き換えの基準": 2.73, "予測": "8 個のうち 3 個以上"},
            "判定保留の範囲": {"置かない": True,
                             "理由": "本番を float32 で回すため（PI 判断 2026-10-08。plan.md §5.3）。"
                                     "bfloat16 の誤差の実測は work/dtype_resolution_011.json に記録として残す"},
            "P-3 副次（三次リンパ組織の 3 個）": {"対象": ["MADCAM1", "CXCL13", "CCL21"],
                                             "条件": "1 個以上が P-2 の合格"},
            "R-6 陽性対照": "陽性対照 71 個のセット平均が 3 群すべてで偽セット（71 個 × 1,000 組）の 99 パーセンタイルを超えなければ判定不能",
            "R-7 偽セットの較正": "新しく引いた偽セット 1,000 組（8 個）に P-1 を当て、合格数が 1% の二項分布を片側 p < 0.05 で上回れば判定不能",
            "乱数の鍵": {"P-1 の偽セット": 20261007, "R-6 の偽セット": 20261008, "R-7 の偽セット": 20261009,
                       "最小検出可能効果量の計算": 20261010}},
        "output": {"experiment_dir": REL, "results_dir": "results"},
        "axis_feasibility": {**af["axis_feasibility"],
                             "注": "着手前チェックは outputs/design/2026-10-07-A2-feasibility/ で実施済み（判定保留、PI が着手を判断 2026-10-07）"},
    }
    (E / "config.yaml").write_text(
        "# EXP-20261007-001 — make_config_011.py が生成。手で直さない。\n"
        + yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False, width=120))

    # ---- 試し実行（判定ではない）: 役割を全部 1 つ以上含め、対照 35 個。出力先は results-smoke/
    smoke = json.loads(json.dumps(cfg))
    smoke["note"] = "試し実行（判定ではない）。make_config_011.py が生成。出力は results-smoke/ のみ"
    smoke["perturbation"]["genes_to_perturb"] = [cand["MADCAM1"], cand["CA12"]]
    smoke["controls"]["positive"] = pos[:2]
    smoke["controls"]["negative_fixed"] = hk[:1]
    smoke["controls"]["extra_fixed"] = [{"role": "negctrl_panel", "genes": neu[:1]}]
    zero = [r["ensembl"] for r in g["対照500"] if r["層"] == "出発点で0"][:7]
    expr = [r["ensembl"] for r in g["対照500"] if r["層"] == "発現あり"][:28]
    smoke["controls"]["negative_sampler"]["genes"] = expr + zero
    smoke["controls"]["negative_random_n"] = 35
    smoke["cell_selection"]["max_cells_per_gene"] = 10
    smoke["cell_selection"]["min_cells_per_gene"] = 10
    # 重心は各状態 40 細胞に絞る（本番の行の先頭 40。形式は本番と同じ）
    spec = json.loads((W / "state_rows-seed81.json").read_text())
    spec_s = {"rows": {k: v[:40] for k, v in spec["rows"].items()}, "注": "試し実行用（本番の行の先頭 40）"}
    for s_ in (81, 82, 83):
        (W / f"state_rows-smoke-seed{s_}.json").write_text(json.dumps(spec_s, ensure_ascii=False) + "\n")
    smoke["state_embs"]["rows_file"] = f"{REL}/work/state_rows-smoke-seed{{seed}}.json"
    smoke["robustness"]["seeds"] = [81, 82, 83]
    smoke["compute"]["device"] = "mps"
    smoke["compute"]["dtype"] = "fp32"
    smoke["output"]["results_dir"] = "results-smoke"
    (W / "config-smoke.yaml").write_text(
        "# 試し実行。make_config_011.py が生成。手で直さない。\n"
        + yaml.safe_dump(smoke, allow_unicode=True, sort_keys=False, width=120))

    # ---- bfloat16 と float32 の食い違いの実測（plan.md §5.3.1）。対照 12 個だけ、本番の細胞集合で
    import random as _r
    rr = _r.Random(20261011)
    l1 = [r["ensembl"] for r in g["対照500"] if r["層"] == "発現あり"]
    l2 = [r["ensembl"] for r in g["対照500"] if r["層"] == "出発点で0"]
    dt_genes = rr.sample(l1, 10) + rr.sample(l2, 2)
    for name, dtype in (("bf16", "bf16"), ("fp32", "fp32")):
        d = json.loads(json.dumps(cfg))
        d.pop("axis_feasibility", None)
        d["note"] = f"bfloat16 と float32 の食い違いの実測（plan.md §5.3.1）。数値型 {dtype}。make_config_011.py が生成"
        d["perturbation"]["genes_to_perturb"] = []
        d["perturbation"]["controls_only"] = True
        d["controls"] = {"positive": [], "negative_fixed": [], "fixed_role": "housekeeping", "extra_fixed": [],
                         "negative_random_n": len(dt_genes),
                         "negative_sampler": {"method": "explicit", "genes": dt_genes}, "comparator": None}
        d["compute"]["dtype"] = dtype
        d["output"]["results_dir"] = f"results-dtype-{name}"
        (W / f"config-dtype-{name}.yaml").write_text(
            "# bfloat16 と float32 の食い違いの実測。make_config_011.py が生成。手で直さない。\n"
            + yaml.safe_dump(d, allow_unicode=True, sort_keys=False, width=120))
    # ---- 2 台に群単位で分ける（plan.md §5.4）。config.yaml の機械的な派生
    for host, seeds_ in (("cf5b", [81]), ("c083", [82, 83])):
        d = json.loads(json.dumps(cfg))
        d["_derived_from"] = "config.yaml（make_config_011.py）。robustness.seeds と output.results_dir だけを変えた"
        d["robustness"]["seeds"] = seeds_
        d["output"]["results_dir"] = f"results-{host}"
        (W / f"config-{host}.yaml").write_text(
            f"# 本番を {host} で回す設定。make_config_011.py が生成。手で直さない。\n"
            + yaml.safe_dump(d, allow_unicode=True, sort_keys=False, width=120))
    (W / "dtype_genes_011.json").write_text(json.dumps(
        {"規則": "対照 500 の層 1 から 10 個・層 2 から 2 個。random.Random(20261011).sample",
         "genes": [{"ensembl": e, "symbol": next(r["symbol"] for r in g["対照500"] if r["ensembl"] == e)} for e in dt_genes]},
        ensure_ascii=False, indent=1) + "\n")
    print("config.yaml・work/config-smoke.yaml・work/config-dtype-{bf16,fp32}.yaml を書いた。測る遺伝子", n_genes)


if __name__ == "__main__":
    main()

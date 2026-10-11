"""EXP-20260929-001 — 3 機体（seed ごと）の結果を results/ に統合し、統合後に照合する。

seed 単位でのみ分割しているので、results-seed{S}/run-{S}-L-1/ はぶつからない。
統合は「移動ではなく複製」にして、元の results-seed{S}/ はそのまま残す（やり直せるようにする）。

統合後の照合（PI 指示 2026-09-28）:
  - seed ごとの遺伝子数と役割ごとの内訳
  - 欠測（skip された遺伝子）の一覧
  - 埋め込み emb/{orig,pert}.npy が全遺伝子 × 全 seed で実在すること
  - 摂動してはいけない遺伝子が 1 つも入っていないこと
  - どの seed をどの機体で計算したか（results/env-remote-seed{S}.json と results/machines.json）

実行: .venv/bin/python experiments/EXP-20260929-001/analysis/merge_results_003.py [--link]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260929-001"
SEEDS = [61, 62, 63]
MACHINE = {61: "spark-cf5b", 62: "spark-c083", 63: "spark-c083"}
# 盲検中の 5 遺伝子のみ。CXCL13 / CCL19 / CCL21 は**本実験では摂動の対象**。
BLOCKED = {"ENSG00000111321", "ENSG00000124762", "ENSG00000141510",
           "ENSG00000118785", "ENSG00000148346"}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--link", action="store_true", help="複製せずに APFS clone（cp -c）で置く")
    ap.parse_args()

    out = E / "results"
    out.mkdir(parents=True, exist_ok=True)
    cells, summ, timing = [], [], []
    report: dict = {"seeds": SEEDS, "machines": MACHINE, "per_seed": {}}

    for s in SEEDS:
        src = E / f"results-seed{s}"
        if not src.exists():
            raise SystemExit(f"{src} がありません。seed {s} の計算が回収されていない")
        rd = src / f"run-{s}-L-1"
        if not rd.exists():
            raise SystemExit(f"{rd} がありません")
        dst = out / rd.name
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(rd, dst, copy_function=shutil.copy2)
        for name, box in (("cells.csv", cells), ("summary.csv", summ)):
            f = src / name
            if f.exists():
                df = pd.read_csv(f, dtype={"donor_id": str})
                df["seed"] = s
                df["machine"] = MACHINE[s]
                box.append(df)
        if (src / "timing.jsonl").exists():
            for line in (src / "timing.jsonl").read_text().splitlines():
                if line.strip():
                    d = json.loads(line)
                    d["machine"] = MACHINE[s]
                    timing.append(d)
        for name in ("env.json", "env-remote.json", "device.json", "config.sha256"):
            f = src / name
            if f.exists():
                shutil.copy2(f, out / f"{Path(name).stem}-seed{s}{Path(name).suffix}")

        genes = sorted(p.name for p in (dst / "genes").iterdir() if p.is_dir())
        missing_emb = [g for g in genes
                       if not ((dst / "genes" / g / "emb/orig.npy").exists()
                               and (dst / "genes" / g / "emb/pert.npy").exists())]
        report["per_seed"][s] = {"machine": MACHINE[s], "n_genes": len(genes),
                                 "blocked_present": sorted(set(genes) & BLOCKED),
                                 "n_missing_emb": len(missing_emb),
                                 "missing_emb": missing_emb[:20]}

    if cells:
        pd.concat(cells, ignore_index=True).to_csv(out / "cells.csv", index=False)
    if summ:
        S = pd.concat(summ, ignore_index=True)
        S.to_csv(out / "summary.csv", index=False)
        report["roles"] = {str(k): int(v) for k, v in
                           S.drop_duplicates(["seed", "gene"]).groupby(["seed", "role"]).size().items()}
    if timing:
        (out / "timing.jsonl").write_text("\n".join(json.dumps(d, ensure_ascii=False) for d in timing) + "\n")

    cfg = yaml.safe_load((E / "config.yaml").read_text())
    expect = (len(cfg["perturbation"]["genes_to_perturb"]) + len(cfg["controls"]["positive"])
              + len(cfg["controls"]["negative_fixed"]) + cfg["controls"]["negative_random_n"]
              + sum(len(x["genes"]) for x in cfg["controls"]["extra_fixed"]))
    report["expected_genes_per_seed"] = expect
    report["all_seeds_complete"] = all(v["n_genes"] == expect for v in report["per_seed"].values())
    report["no_blocked_genes"] = all(not v["blocked_present"] for v in report["per_seed"].values())
    report["no_missing_embeddings"] = all(v["n_missing_emb"] == 0 for v in report["per_seed"].values())
    report["merged_files_sha256"] = {n: sha256(out / n) for n in ("cells.csv", "summary.csv")
                                     if (out / n).exists()}
    (out / "machines.json").write_text(json.dumps(
        {"assignment": MACHINE,
         "note": "seed 単位でのみ分割している。1 機体が 1 seed の全 1,136 遺伝子を担当する。"
                 "遺伝子で分割すると、機体差が『経路セットは対照より動く』と同じ形で現れて見分けられない。"},
        ensure_ascii=False, indent=1) + "\n")
    (out / "merge_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "per_seed"}, ensure_ascii=False, indent=1))
    for s, v in report["per_seed"].items():
        print(f"seed {s} ({v['machine']}): 遺伝子 {v['n_genes']} / 期待 {expect}, "
              f"埋め込み欠け {v['n_missing_emb']}, 遮断遺伝子 {v['blocked_present']}")


if __name__ == "__main__":
    main()

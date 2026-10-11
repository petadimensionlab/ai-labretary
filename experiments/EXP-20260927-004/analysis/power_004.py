"""EXP-20260927-004 — 検出力の見積もり（実データのばらつきから。合成データは使わない）。

guardrails 2026-09-27「合成データで検出力を事前に確かめるときは、遺伝子間のばらつきを実データから
推定して入れる」に従い、EXP-20260927-003 の対照遺伝子 80 × 3 seed の実測のばらつきだけを使う。

出すもの（候補の指標 A0 / A11 / A12 それぞれについて）:
  1. 対照遺伝子群のばらつきの分解（遺伝子固有で再現する成分 / seed 間の誤差）
  2. 偶然合格率の実測（対照遺伝子に主判定と同じ規則を当てはめる。seed 間の独立を仮定しない）
  3. 検出下限（1 遺伝子 100 細胞＝EXP-003 の実測、および 200 細胞＝本実験の設計）
  4. 参考: 筋線維芽細胞の軸での陽性対照 ACTA2 / TAGLN の効果の大きさ（標準偏差単位）

**主対象 LTBR / RELB / NFKB2 と CDKN1A、盲検の 3 遺伝子は metric_candidates.py が読み込み直後に
除外しており、本スクリプトもその除外後のデータしか触らない。**

出力: work/power_004.json
実行: .venv/bin/python experiments/EXP-20260927-004/analysis/power_004.py
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260927-004"
SRC = ROOT / "experiments/EXP-20260927-003"

_spec = importlib.util.spec_from_file_location("mc", E / "analysis/metric_candidates.py")
mc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mc)

CANDIDATES = ["A0", "A11", "A12"]
Q_HI, MIN_SEEDS, N_SEEDS = 0.90, 2, 3
N_CELLS_EXP003, N_CELLS_EXP004 = 100, 200


def power_at(delta: float, sd_null: float, sd_seed: float, thr_z: float = 1.2816) -> float:
    """真の効果 delta のとき、3 seed 中 2 以上で「対照遺伝子群の上位 10%」を超える確率。"""
    p = float(stats.norm.sf((thr_z * sd_null - delta) / sd_seed))
    return 3 * p ** 2 * (1 - p) + p ** 3


def solve_delta(target: float, sd_null: float, sd_seed: float) -> float:
    f = lambda x: power_at(x, sd_null, sd_seed) - target
    return float(optimize.brentq(f, 0.0, 20 * sd_null))


def main() -> None:
    cache = E / "work/cells_with_shifts.parquet"
    assert cache.exists(), f"{cache} が無い。先に metric_candidates.py を実行すること"
    c = pd.read_parquet(cache)
    assert not set(c.gene) & set(mc.EXCLUDE), "除外した遺伝子が残っている"
    mc.add_candidates(c)
    cfg = yaml.safe_load(open(SRC / "config.yaml"))
    null_genes = set(cfg["controls"]["negative_sampler"]["genes"])
    seeds = sorted(c.seed.unique().tolist())

    per_seed_nulls, per_seed_elig = {}, {}
    for s in seeds:
        cs = c[c.seed == s]
        nulls = {h: hd for h, hd in cs[cs.gene.isin(null_genes)].groupby("gene")}
        per_seed_nulls[s] = nulls
        per_seed_elig[s] = mc.eligibility(cs, nulls)

    out = {
        "note": "EXP-20260927-003 の実測（対照遺伝子 80 × 3 seed、1 遺伝子 100 細胞）から見積もった。"
                "主対象 LTBR/RELB/NFKB2・CDKN1A・盲検の 3 遺伝子は含まない。"
                "EXP-20260927-004 の 200 細胞への外挿は、seed 間の誤差の分散が細胞数に反比例すると仮定した。",
        "n_cells_exp003": N_CELLS_EXP003, "n_cells_exp004": N_CELLS_EXP004,
        "criteria": {"q_hi": Q_HI, "min_seeds": MIN_SEEDS, "n_seeds": N_SEEDS},
        "candidates": {},
    }
    for tag in CANDIDATES:
        frames = []
        for s in seeds:
            f = mc.gene_stats(c[c.seed == s], per_seed_nulls[s], per_seed_elig[s], f"m_{tag}")
            f["seed"] = s
            frames.append(f)
        M = pd.concat(frames, ignore_index=True)
        N = M[M.gene.isin(null_genes)]

        # 1. ばらつきの分解
        w = N.pivot(index="gene", columns="seed", values="D_raw").dropna()
        v_between = float(w.mean(axis=1).var(ddof=1))
        v_within = float(w.var(axis=1, ddof=1).mean())
        v_gene = max(v_between - v_within / len(seeds), 0.0)
        sd_null_100 = float(np.mean([N[N.seed == s].D_raw.std(ddof=1) for s in seeds]))
        sd_seed_100 = float(np.sqrt(v_within))

        # 2. 偶然合格率の実測（seed 間の相関をそのまま含む）
        q = N.pivot(index="gene", columns="seed", values="q")
        dpos = N.pivot(index="gene", columns="seed", values="D_raw") > 0
        passes = ((q >= Q_HI) & dpos).sum(axis=1)
        n_valid = q.notna().sum(axis=1)
        chance = float((passes[n_valid >= MIN_SEEDS] >= MIN_SEEDS).mean())
        thr_emp = [float(np.percentile(N[N.seed == s].D_raw, 90)) for s in seeds]

        # 3. 検出下限
        det = {}
        for n in (N_CELLS_EXP003, N_CELLS_EXP004):
            k = N_CELLS_EXP003 / n
            sd_null = float(np.sqrt(v_gene + v_within * k))
            sd_seed = float(np.sqrt(v_within * k))
            det[str(n)] = {
                "sd_null_per_seed": sd_null, "sd_seed_error": sd_seed,
                "delta_power50": solve_delta(0.50, sd_null, sd_seed),
                "delta_power80": solve_delta(0.80, sd_null, sd_seed),
                "delta_power80_in_sd": solve_delta(0.80, sd_null, sd_seed) / sd_null,
            }

        # 3b. 事前登録する判定規則（q ≥ 0.90 かつ D > 0 かつ D ≥ K_SD × 対照の標準偏差 を 3 seed すべて）
        # のもとでの検出下限。順位の条件より大きさの条件のほうが効くので、こちらが実効的な下限になる。
        K_SD = 2.5
        for n in (N_CELLS_EXP003, N_CELLS_EXP004):
            e = det[str(n)]
            thr = max(1.2816, K_SD) * e["sd_null_per_seed"]
            # 3 seed すべてで超える確率が 0.80 になる真の効果
            z = float(stats.norm.isf(1 - 0.80 ** (1 / 3)))
            d80 = thr + z * e["sd_seed_error"]
            e["rule_k_sd"] = K_SD
            e["delta_power80_under_rule"] = d80
            e["delta_power80_under_rule_in_sd"] = d80 / e["sd_null_per_seed"]
        print(f"      事前登録の規則（{K_SD} 標準偏差・3 seed すべて）での検出下限: "
              f"{det[str(N_CELLS_EXP004)]['delta_power80_under_rule']:.3e}"
              f"（{det[str(N_CELLS_EXP004)]['delta_power80_under_rule_in_sd']:.2f} 標準偏差, 200 細胞）")

        # 4. 筋線維芽細胞の軸での陽性対照（鏡像の指標。参照値）
        xf = []
        for s in seeds:
            f = mc.gene_stats(c[c.seed == s], per_seed_nulls[s], per_seed_elig[s], f"x_{tag}")
            f["seed"] = s
            xf.append(f)
        X = pd.concat(xf, ignore_index=True)
        ref = {}
        for gid, nm in mc.POS_MYOF.items():
            rows = X[X.gene == gid].set_index("seed")
            base_m = {s: X[(X.seed == s) & X.gene.isin(null_genes)].D_raw.mean() for s in seeds}
            base_s = {s: X[(X.seed == s) & X.gene.isin(null_genes)].D_raw.std(ddof=1) for s in seeds}
            ref[nm] = {"D_raw": [float(rows.D_raw[s]) for s in seeds],
                       "sd_units": [round(float((rows.D_raw[s] - base_m[s]) / base_s[s]), 2) for s in seeds]}

        # 2b. 判定規則の候補ごとの偶然合格率と、陽性対照候補 10 個での「1 つ以上合格」の確率
        # （seed 間の相関は実測のまま。10 個ぶんは対照遺伝子から 10 個を復元抽出して数える）
        rng = np.random.default_rng(0)
        rules = {}
        for qh in (0.90, 0.95, 0.98):
            for ms in (2, 3):
                ps = ((q >= qh) & dpos).sum(axis=1)
                ok = (ps >= ms) & (n_valid >= ms)
                per_gene = float(ok[n_valid >= ms].mean())
                arr = ok[n_valid >= ms].to_numpy()
                fam = float(np.mean([arr[rng.integers(0, len(arr), 10)].any() for _ in range(20000)]))
                rules[f"q>={qh}_seeds>={ms}"] = {"per_gene": round(per_gene, 4),
                                                 "family_of_10": round(fam, 4)}
        # 効果の大きさ（対照遺伝子群の標準偏差単位）を足した規則。
        # 対照遺伝子群のばらつきの 83% が「遺伝子固有で再現する成分」なので、順位だけの規則では
        # seed を重ねても偶然合格率が下がらない。絶対的な大きさの条件を足して分離する。
        Dn = N.pivot(index="gene", columns="seed", values="D_raw")
        Z = Dn.copy()
        for s in seeds:
            Z[s] = (Dn[s] - Dn[s].mean()) / Dn[s].std(ddof=1)
        for kk in (2.0, 2.5, 3.0):
            ok = ((q >= 0.90) & dpos & (Z >= kk)).sum(axis=1) >= 3
            ok = ok[n_valid >= 3]
            arr = ok.to_numpy()
            fam = float(np.mean([arr[rng.integers(0, len(arr), 10)].any() for _ in range(20000)]))
            rules[f"q>=0.9_seeds=3_and_{kk}SD"] = {"per_gene": round(float(ok.mean()), 4),
                                                   "family_of_10": round(fam, 4)}
        print(f"      判定規則ごとの偶然合格率（1 遺伝子 / 10 個中 1 つ以上）: " +
              ", ".join(f"{k} {v['per_gene']:.3f}/{v['family_of_10']:.3f}" for k, v in rules.items()))

        out["candidates"][tag] = {
            "rule_variants": rules,
            "desc": {"A0": "目標への接近（EXP-003 の主指標）",
                     "A11": "目標への接近 − (出発点・対照の状態・筋線維芽細胞の平均)",
                     "A12": "A11 / 移動量 ||Δ||"}[tag],
            "sd_null_per_seed_100cells": sd_null_100,
            "sd_seed_error_100cells": sd_seed_100,
            "sd_gene_specific": float(np.sqrt(v_gene)),
            "var_gene_specific_frac": round(v_gene / (v_gene + v_within), 3) if v_gene + v_within > 0 else None,
            "empirical_pass_line_p90": thr_emp,
            "chance_pass_rate": chance,
            "n_null_genes_used": int((n_valid >= MIN_SEEDS).sum()),
            "detection_limit": det,
            "myof_axis_reference": ref,
        }
        print(f"[{tag}] 対照 SD(100 細胞) {sd_null_100:.3e}、seed 誤差 {sd_seed_100:.3e}、"
              f"遺伝子固有の割合 {out['candidates'][tag]['var_gene_specific_frac']}")
        print(f"      偶然合格率 {chance:.3f}（対照 {out['candidates'][tag]['n_null_genes_used']} 遺伝子）")
        for n in (N_CELLS_EXP003, N_CELLS_EXP004):
            e = det[str(n)]
            print(f"      {n} 細胞: 検出下限 検出力50% {e['delta_power50']:.3e} / "
                  f"80% {e['delta_power80']:.3e}（{e['delta_power80_in_sd']:.2f} 標準偏差）")
        print(f"      参考 ACTA2（筋線維芽細胞の軸）: {ref['ACTA2']['sd_units']} 標準偏差")

    (E / "work/power_004.json").write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    print(f"\n書き出し: {E / 'work/power_004.json'}")


if __name__ == "__main__":
    main()

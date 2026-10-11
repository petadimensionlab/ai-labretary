#!/usr/bin/env python3
"""EXP-20261007-002 段階 1 — 著者ラベルの退化した線維芽細胞（dFIB）は、参照（FIB）より SenMayo 得点が高いか。

**摂動を使わない。Geneformer も使わない。** 生カウント（`raw.X`）だけを読む。GPU 0 時間。

この道具は plan.md §4 の判定規則をそのまま実装し、`config.yaml` の `decision:` 節が plan と同じ値かを検査する。
**判定（支持／否定／判定不能）は確定しない。** ゲート G3 は PI が行う。ここが出すのは値と提案まで。

走らせ方:
    .venv/bin/python experiments/EXP-20261007-002/analysis/score_senescence_012.py            # 本実行（A2 + A1）
    .venv/bin/python experiments/EXP-20261007-002/analysis/score_senescence_012.py --smoke    # A2 の小さな試し実行（判定ではない）

主な入力: plan.md（§2〜§7・§12）、config.yaml、senmayo_125.csv、lake2023_table12.csv、
          A2/A1 の h5ad（読み取り専用）、Geneformer の遺伝子記号の辞書。
主な出力: results/criteria_012.json、results/scores_012.csv、results/README.md、results/run.log
          （試し実行は results-smoke/ に分ける）。
"""
from __future__ import annotations

import argparse
import csv
import datetime as _dt
import hashlib
import json
import pickle
import random
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
EXP = ROOT / "experiments" / "EXP-20261007-002"

# ---------------------------------------------------------------------------
# plan.md §4 の判定値（**ここを変えてはいけない。config.yaml と照合する**）
PLAN = {
    "S1": {
        "test": "welch_t_two_sample_one_sided",
        "effect_size": "cohens_d_pooled",
        "d_threshold": 0.5,
        "p_threshold": 0.05,
        "sided": "one_sided",
        "min_detectable_effect_d": 0.38,
        "power_at_d_0_5": 0.95,
    },
    "S2": {
        "n_random_sets": 1000,
        "percentile": 99,
        "percentile_method": "linear_numpy_default",
        "random_key": 20261012,
        "shortage_rule": "nearest_stratum_with_surplus",
        "strata": "deciles_of_senmayo_frequency",
        "freq_definition": "fraction_of_interval_cells_with_raw_count_gt_zero",
        "draw_without_replacement": True,
    },
    "S3": {"ratio_threshold": 0.3333333333},
    "S4": {"random_key": 20261013},
    "T5_min_gene_count": 25,
    "T6_min_cells": 30,
    "evaluation_order": ["T5", "T6", "T2", "T1", "T3", "T4"],
    "metric": {
        "score": "fraction_of_list_genes_detected",
        "detection_rule": "raw_count_gt_zero",
        "denominator": 122,
        "standardization": "reference_only_within_interval",
        "sd_denominator": "unbiased_n_minus_1",
    },
    "expected_reference_cells": 81,
    "expected_degenerate_cells": 93,
}

_A2_INTERVAL = (1000, 2000)        # plan §2.1
_A1_INTERVAL = (500, 1100)         # plan §4.4
_A2_SENSITIVITY = (800, 2500)      # plan §6.2（記述のみ）
_A1_SENSITIVITY = (600, 1200)      # plan §6.2（記述のみ）
N_STRATA = 10

T_KEY = {"T1": "T1_effect_below_0.5", "T2": "T2_direction_reversed", "T3": "T3_not_specific",
         "T4": "T4_depth_explains", "T5": "T5_indicator_not_viable", "T6": "T6_too_few_cells"}


# ---------------------------------------------------------------------------
class Log:
    def __init__(self) -> None:
        self.lines: list[str] = []

    def __call__(self, msg: str) -> None:
        stamp = _dt.datetime.now().strftime("%H:%M:%S")
        line = f"[{stamp}] {msg}"
        print(line, flush=True)
        self.lines.append(line)

    def text(self) -> str:
        return "\n".join(self.lines) + "\n"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def read_column(g):
    """h5ad の列（categorical の Group か dataset）を配列に直す。文字列は str にする。"""
    if isinstance(g, h5py.Group):
        codes = g["codes"][:]
        cats = g["categories"][:]
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in cats], dtype=object)
        return cats[codes]
    arr = g[:]
    if arr.dtype.kind in ("S", "O", "U"):
        return np.asarray([a.decode() if isinstance(a, bytes) else str(a) for a in arr], dtype=object)
    return arr


def read_raw_rows(h5path: Path, rows) -> dict:
    """生カウント（raw/X, CSR）の指定した行だけを読む。{row: (列番号, 値)}"""
    rows = sorted({int(r) for r in rows})
    out: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    with h5py.File(h5path, "r") as f:
        indptr = f["raw/X/indptr"][:]
        data = f["raw/X/data"]
        indices = f["raw/X/indices"]
        for r in rows:
            s, e = int(indptr[r]), int(indptr[r + 1])
            if e > s:
                out[r] = (indices[s:e], data[s:e])
            else:
                out[r] = (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.float32))
    return out


def col_counts(rows_dict: dict, row_order) -> dict:
    """区間の細胞で、列（遺伝子）ごとに raw カウント > 0 の細胞数を数える。"""
    cnt: dict[int, int] = {}
    for r in row_order:
        ids, vals = rows_dict[int(r)]
        if len(ids) == 0:
            continue
        for c in np.unique(ids[vals > 0]):
            c = int(c)
            cnt[c] = cnt.get(c, 0) + 1
    return cnt


def detection_matrix(rows_dict: dict, row_order, cols_sorted: np.ndarray) -> np.ndarray:
    """指定した列集合について、細胞 × 遺伝子の検出（raw > 0）の真偽行列を作る。"""
    col_pos = {int(c): i for i, c in enumerate(cols_sorted)}
    D = np.zeros((len(row_order), len(cols_sorted)), dtype=bool)
    for ri, r in enumerate(row_order):
        ids, vals = rows_dict[int(r)]
        if len(ids) == 0:
            continue
        m = vals > 0
        for c in ids[m]:
            p = col_pos.get(int(c))
            if p is not None:
                D[ri, p] = True
    return D


def cohens_d_pooled(deg: np.ndarray, ref: np.ndarray) -> float:
    """プールした標準偏差（分母 n−1）で割った Cohen's d。plan §4.1。"""
    n1, n2 = len(deg), len(ref)
    if n1 < 2 or n2 < 2:
        return float("nan")
    v1 = float(np.var(deg, ddof=1))
    v2 = float(np.var(ref, ddof=1))
    sp2 = ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)
    if sp2 <= 0:
        return float("nan")
    return float((float(np.mean(deg)) - float(np.mean(ref))) / np.sqrt(sp2))


def welch_one_sided(deg: np.ndarray, ref: np.ndarray) -> float:
    """Welch の 2 標本片側 t 検定（対立仮説: 退化 > 参照）。"""
    if len(deg) < 2 or len(ref) < 2:
        return float("nan")
    return float(stats.ttest_ind(deg, ref, equal_var=False, alternative="greater").pvalue)


# ---------------------------------------------------------------------------
def stratum_index(freq: float, qs: np.ndarray):
    """層 j（1〜10）。plan §4.2: 層 j は q_{j-1} ≤ freq < q_j、最上位は freq ≤ q_10。"""
    if freq < qs[0] or freq > qs[10]:
        return None
    p = int(np.searchsorted(qs, freq, side="right"))
    return max(1, min(p, N_STRATA))


def draw_random_sets(pool_freq: dict, senmayo_freq, key, n_sets, denom, log):
    """頻度をそろえた無作為 122 遺伝子集合を n_sets 組つくる（plan §4.2）。"""
    qs = np.percentile(np.asarray(senmayo_freq, dtype=float), np.arange(0, 101, 10))
    k = np.zeros(N_STRATA + 1, dtype=int)
    for f_ in senmayo_freq:
        j = stratum_index(float(f_), qs)
        if j is not None:
            k[j] += 1
    if int(k[1:].sum()) != denom:
        raise SystemExit(f"SenMayo の層分けの合計が {int(k[1:].sum())}（{denom} でない）。中止する")

    strata: dict[int, list[str]] = {j: [] for j in range(1, N_STRATA + 1)}
    outside = 0
    for name, f_ in pool_freq.items():
        j = stratum_index(float(f_), qs)
        if j is None:
            outside += 1
        else:
            strata[j].append(name)
    for j in strata:
        strata[j].sort()
    pool_median = {j: (float(np.median([pool_freq[g] for g in strata[j]])) if strata[j] else float("nan"))
                   for j in range(1, N_STRATA + 1)}

    avail = {j: list(strata[j]) for j in range(1, N_STRATA + 1)}
    shortage_records = []
    for j in range(1, N_STRATA + 1):
        need = int(k[j]) - len(avail[j])
        if need <= 0:
            continue
        donors = [d for d in range(1, N_STRATA + 1)
                  if d != j and (len(avail[d]) - int(k[d])) > 0 and not np.isnan(pool_median[d])]
        donors.sort(key=lambda d: abs(pool_median[d] - pool_median[j]))
        for d in donors:
            if need <= 0:
                break
            surplus = len(avail[d]) - int(k[d])
            take = min(need, surplus)
            moved = [avail[d].pop() for _ in range(take)]   # 決定的に末尾から
            avail[j].extend(moved)
            shortage_records.append({"deficit_stratum": j, "donor_stratum": d, "moved": int(take)})
            need -= take
        if need > 0:
            raise SystemExit(f"層 {j} のプールが足りず、補っても {need} 個不足。中止して報告する（plan §4.2 手順 6）")

    rng = random.Random(int(key))
    sets = []
    for _ in range(int(n_sets)):
        chosen: list[str] = []
        for j in range(1, N_STRATA + 1):     # 層の番号の昇順に sample を呼ぶ
            kj = int(k[j])
            if kj > 0:
                chosen.extend(rng.sample(avail[j], kj))
        sets.append(chosen)
    if not all(len(s) == denom for s in sets):
        raise SystemExit("無作為集合の個数が 122 でない。中止する")
    info = {
        "boundaries": [float(x) for x in qs],
        "senmayo_stratum_counts": {str(j): int(k[j]) for j in range(1, N_STRATA + 1)},
        "pool_stratum_counts": {str(j): len(strata[j]) for j in range(1, N_STRATA + 1)},
        "pool_outside_strata": int(outside),
        "shortage_records": shortage_records,
        "n_sets": int(n_sets),
        "random_key": int(key),
        "note": "各層のプールは遺伝子記号の辞書順に並べ、層の番号の昇順に random.Random(key).sample を呼ぶ（再現のため）",
    }
    return sets, info


# ---------------------------------------------------------------------------
def check_config(cfg: dict) -> list[str]:
    errs: list[str] = []

    def eq(a, b, label):
        if a != b:
            errs.append(f"{label}: config={a!r} / plan={b!r}")

    d = cfg["decision"]
    s1 = d["S-1 主判定"]
    eq(s1.get("test"), PLAN["S1"]["test"], "decision.S-1.test")
    eq(s1.get("effect_size"), PLAN["S1"]["effect_size"], "decision.S-1.effect_size")
    eq(s1.get("effect_size_threshold"), PLAN["S1"]["d_threshold"], "decision.S-1.effect_size_threshold")
    eq(s1.get("p_threshold"), PLAN["S1"]["p_threshold"], "decision.S-1.p_threshold")
    eq(s1.get("p_sided"), PLAN["S1"]["sided"], "decision.S-1.p_sided")
    eq(s1.get("min_detectable_effect_d"), PLAN["S1"]["min_detectable_effect_d"], "decision.S-1.min_detectable_effect_d")
    eq(s1.get("power_at_d_0_5"), PLAN["S1"]["power_at_d_0_5"], "decision.S-1.power_at_d_0_5")

    s2 = d["S-2 特異性"]
    eq(s2.get("n_random_sets"), PLAN["S2"]["n_random_sets"], "decision.S-2.n_random_sets")
    eq(s2.get("percentile"), PLAN["S2"]["percentile"], "decision.S-2.percentile")
    eq(s2.get("percentile_method"), PLAN["S2"]["percentile_method"], "decision.S-2.percentile_method")
    eq(s2.get("random_key"), PLAN["S2"]["random_key"], "decision.S-2.random_key")
    fm = s2.get("frequency_matching", {})
    eq(fm.get("strata"), PLAN["S2"]["strata"], "decision.S-2.strata")
    eq(fm.get("freq_definition"), PLAN["S2"]["freq_definition"], "decision.S-2.freq_definition")
    eq(bool(fm.get("draw_without_replacement")), PLAN["S2"]["draw_without_replacement"], "decision.S-2.draw_without_replacement")
    eq(fm.get("shortage_rule"), PLAN["S2"]["shortage_rule"], "decision.S-2.shortage_rule")

    eq(d["S-3 深さ依存"].get("ratio_threshold"), PLAN["S3"]["ratio_threshold"], "decision.S-3.ratio_threshold")
    eq(d["T-5"]["condition"], "usable_gene_count_lt_25", "decision.T-5.condition")
    eq(d["T-6"]["condition"], "cells_lt_30_in_either_state", "decision.T-6.condition")
    eq(list(d["evaluation_order"]), PLAN["evaluation_order"], "decision.evaluation_order")

    m = cfg["metric"]
    eq(m.get("score"), PLAN["metric"]["score"], "metric.score")
    eq(m.get("detection_rule"), PLAN["metric"]["detection_rule"], "metric.detection_rule")
    eq(m.get("denominator"), PLAN["metric"]["denominator"], "metric.denominator")
    eq(m.get("standardization"), PLAN["metric"]["standardization"], "metric.standardization")
    eq(m.get("sd_denominator"), PLAN["metric"]["sd_denominator"], "metric.sd_denominator")

    eq(cfg["genes"].get("expected_gene_count"), PLAN["metric"]["denominator"], "genes.expected_gene_count")
    eq(cfg["genes"].get("min_gene_count_T5"), PLAN["T5_min_gene_count"], "genes.min_gene_count_T5")
    eq(cfg["data"]["A2"]["interval"]["min"], _A2_INTERVAL[0], "data.A2.interval.min")
    eq(cfg["data"]["A2"]["interval"]["max"], _A2_INTERVAL[1], "data.A2.interval.max")
    eq(cfg["data"]["A1"]["interval"]["min"], _A1_INTERVAL[0], "data.A1.interval.min")
    eq(cfg["data"]["A1"]["interval"]["max"], _A1_INTERVAL[1], "data.A1.interval.max")
    eq(cfg["data"]["A2"]["min_cells_per_state"], PLAN["T6_min_cells"], "data.A2.min_cells_per_state")
    eq(cfg["data"]["A2"]["depth_col"], "nFeature_RNA", "data.A2.depth_col")
    eq(cfg["data"]["A1"]["depth_col"], "nFeature_RNA", "data.A1.depth_col")
    return errs


def build_122(cfg: dict, log: Log) -> list[str]:
    """plan §3.1 の規則で senmayo_125.csv から 122 個を機械的に作る。config と同じかを検査する。"""
    rows = list(csv.DictReader((ROOT / cfg["genes"]["senmayo_file"]).open()))
    col = cfg["genes"]["senmayo_gene_column"]
    in_a2_col = cfg["genes"]["senmayo_in_A2_column"]
    in_a2 = [r[col] for r in rows if r[in_a2_col] == "True"]
    excl_not_in_a2 = set(cfg["genes"]["exclude_not_in_A2"])
    excl_blind = set(cfg["genes"]["exclude_blind_from_list"])
    got = sorted(g for g in in_a2 if g not in excl_not_in_a2 and g not in excl_blind)
    expected = list(cfg["genes"]["senmayo_122"])
    if got != expected:
        only_mech = sorted(set(got) - set(expected))
        only_cfg = sorted(set(expected) - set(got))
        raise SystemExit(f"122 個が plan §3.1 の機械的な作り方と違う（機械のみ {only_mech} / config のみ {only_cfg}）。中止する")
    if len(got) != cfg["genes"]["expected_gene_count"]:
        raise SystemExit(f"122 個でない（{len(got)}）。中止する")
    log(f"122 遺伝子リスト: senmayo_125.csv から機械的に作成し、config と一致（{len(got)} 個）")
    return got


# ---------------------------------------------------------------------------
def analyze_dataset(tag, dcfg, senmayo_names, blocked, geneformer_symbols, log,
                    smoke=False, limit_per_state=None, n_random_override=None):
    """1 データ（A2 または A1）の SenMayo 得点と S-1〜S-3・LODO・記述を計算する。"""
    h5path = ROOT / dcfg["input_h5ad"]
    lo, hi = dcfg["interval"]["min"], dcfg["interval"]["max"]
    ref_labels = list(dcfg["reference_labels"])
    deg_labels = list(dcfg["degenerate_labels"])
    obs = {}
    with h5py.File(h5path, "r") as f:
        for name, target in (("state", dcfg["state_col"]), ("nfeature", dcfg["depth_col"]),
                             ("donor", dcfg["donor_col"])):
            obs[name] = read_column(f["obs"][target])
        obs["cell_id"] = read_column(f["obs"]["_index"]) if "_index" in f["obs"] else np.asarray(
            [f"cell{i}" for i in range(len(obs["state"]))], dtype=object)
    n = len(obs["state"])
    state, nfeat, donor, cell_id = obs["state"], obs["nfeature"].astype(float), obs["donor"], obs["cell_id"]
    log(f"[{tag}] 細胞 {n} 個。参照ラベル {ref_labels} / 比較対象ラベル {deg_labels}")

    labels_present = set(np.unique(state).tolist())
    for lab in ref_labels + deg_labels:
        if lab not in labels_present:
            raise SystemExit(f"[{tag}] A-2: 状態ラベル {lab!r} がデータに無い。中止する")

    is_ref = np.isin(state, ref_labels)
    is_deg = np.isin(state, deg_labels)
    in_int = (nfeat >= lo) & (nfeat <= hi)
    ref_int = np.where(is_ref & in_int)[0]
    deg_int = np.where(is_deg & in_int)[0]
    if limit_per_state:
        ref_int = ref_int[:limit_per_state]
        deg_int = deg_int[:limit_per_state]
    interval_cells = np.concatenate([ref_int, deg_int])
    n_ref, n_deg = len(ref_int), len(deg_int)
    donors_ref = sorted(set(donor[ref_int].tolist()))
    donors_deg = sorted(set(donor[deg_int].tolist()))
    log(f"[{tag}] 区間 {lo}〜{hi}: 参照 {n_ref} 細胞 / {len(donors_ref)} 提供者、比較対象 {n_deg} 細胞 / {len(donors_deg)} 提供者")

    expected = dcfg.get("expected_cells_in_interval") or {}
    cells_report = {
        "interval": [lo, hi],
        "reference_cells": int(n_ref), "degenerate_cells": int(n_deg),
        "reference_donors": len(donors_ref), "degenerate_donors": len(donors_deg),
        "expected_reference": expected.get("reference"), "expected_degenerate": expected.get("degenerate"),
        "matches_expected": bool(expected.get("reference") == n_ref and expected.get("degenerate") == n_deg),
        "all_reference_cells": int(is_ref.sum()), "all_degenerate_cells": int(is_deg.sum()),
        "limited_for_smoke": bool(limit_per_state),
    }
    if expected and not cells_report["matches_expected"]:
        log(f"[{tag}] 注意: 区間の細胞数が仮説 §2.2 の実測（参照 {expected.get('reference')} / "
            f"比較対象 {expected.get('degenerate')}）と違う。記録し、plan の規則どおり進める（区間は変えない）")

    # T-6（着手しない条件。試し実行では止めない）
    t6 = (n_ref < PLAN["T6_min_cells"]) or (n_deg < PLAN["T6_min_cells"])
    if t6 and not smoke:
        log(f"[{tag}] T-6 該当: 区間のどちらかの状態が {PLAN['T6_min_cells']} 細胞未満（参照 {n_ref} / 比較対象 {n_deg}）。実行を始めない")
        raise SystemExit("T-6: 区間の取り方を先に決め直して PI に諮る（中止）")

    # 生カウントを読む（区間＋そろえる前＋感度区間）
    if tag == "A2":
        read_mask = is_ref | is_deg
        sens_lo, sens_hi, sens_label = _A2_SENSITIVITY[0], _A2_SENSITIVITY[1], "A2_800_2500"
    else:
        sens_lo, sens_hi, sens_label = _A1_SENSITIVITY[0], _A1_SENSITIVITY[1], "A1_600_1200"
        read_mask = (is_ref | is_deg) & (nfeat >= min(lo, sens_lo)) & (nfeat <= max(hi, sens_hi))
    read_rows = np.where(read_mask)[0]
    log(f"[{tag}] 生カウントを読む細胞: {len(read_rows)} 個")
    rows = read_raw_rows(h5path, read_rows)

    # A-3 生カウントの検査（整数様・行和と nCount_RNA の整合）
    sample_vals = np.concatenate([rows[r][1] for r in read_rows[:200]]) if len(read_rows) else np.array([0.0])
    integer_like = bool(np.all(sample_vals == np.round(sample_vals)))
    rowsum = np.array([float(rows[r][1].sum()) for r in read_rows])
    with h5py.File(h5path, "r") as f:
        ncount = read_column(f["obs"]["nCount_RNA"]).astype(float)[read_rows]
    ratio = rowsum / np.maximum(ncount, 1.0)
    med_ratio = float(np.median(ratio)) if len(ratio) else float("nan")
    a3 = {"counts_integer_like": integer_like, "median_rowsum_over_nCount": med_ratio, "n_checked": int(len(read_rows))}
    if not integer_like:
        raise SystemExit(f"[{tag}] A-3: raw.X が整数様でない。中止する（正規化値を数えている可能性）")
    if not (0.95 <= med_ratio <= 1.05):
        raise SystemExit(f"[{tag}] A-3: 行和と nCount_RNA が整合しない（中央値比 {med_ratio:.4f}）。中止する")
    log(f"[{tag}] A-3 合格: 整数様={integer_like}、行和/nCount_RNA の中央値比={med_ratio:.4f}")

    with h5py.File(h5path, "r") as f:
        feat = read_column(f["var"]["feature_name"])
    name_to_cols: dict[str, list[int]] = {}
    for i, g in enumerate(feat):
        name_to_cols.setdefault(str(g), []).append(i)
    n_dup = sum(1 for v in name_to_cols.values() if len(v) > 1)

    # A-4 / T-5
    usable = [g for g in senmayo_names if g in name_to_cols]
    usable_in_dict = [g for g in usable if g in geneformer_symbols]
    n_usable = len(usable_in_dict)
    a4 = {"senmayo_122": len(senmayo_names), "usable_in_var_and_dict": n_usable,
          "n_missing_in_var": len(senmayo_names) - len(usable),
          "duplicate_feature_names_in_var": int(n_dup)}
    if n_usable < PLAN["T5_min_gene_count"]:
        log(f"[{tag}] T-5 該当: 使える遺伝子が {n_usable} 個（<25）。指標が成り立たない")
        raise SystemExit("T-5: 指標が成り立たない（中止）")
    if n_usable < len(senmayo_names):
        raise SystemExit(f"[{tag}] A-4: 使える SenMayo 遺伝子が {n_usable}/{len(senmayo_names)}。25 個以上だが不足。中止して報告する")
    log(f"[{tag}] A-4 合格: 122 個すべてが var と辞書にある（重複名 {n_dup} 件）")

    gene122_cols = [name_to_cols[g][0] for g in senmayo_names]

    interval_rows = list(interval_cells)
    cnt = col_counts(rows, interval_rows)
    n_int = len(interval_rows)
    senmayo_freq = np.array([cnt.get(int(c), 0) for c in gene122_cols], dtype=float) / float(n_int)
    log(f"[{tag}] 区間の細胞で頻度を計算（分母 {n_int}）")

    row_index = {int(r): i for i, r in enumerate(read_rows)}
    ref_rows_i = np.array([row_index[int(r)] for r in ref_int])
    deg_rows_i = np.array([row_index[int(r)] for r in deg_int])

    # プール（plan §4.2 手順 4）
    sen_set = set(senmayo_names)
    pool_names, pool_freq = [], {}
    for name, cols in name_to_cols.items():
        if name in sen_set or name in blocked or name not in geneformer_symbols:
            continue
        f_ = cnt.get(int(cols[0]), 0) / float(n_int)
        if f_ <= 0:
            continue
        pool_names.append(name)
        pool_freq[name] = f_
    pool_names.sort()
    log(f"[{tag}] 無作為抽出のプール: {len(pool_names)} 個（A2/A1 の var にあり・辞書にあり・盲検と SenMayo を除く・区間で freq>0）")

    all_cols = sorted(set(gene122_cols) | {name_to_cols[g][0] for g in pool_names})
    D = detection_matrix(rows, read_rows, np.asarray(all_cols, dtype=np.int64))
    cpos = {int(c): i for i, c in enumerate(all_cols)}
    idx122 = [cpos[int(c)] for c in gene122_cols]

    def scores_of(cols):
        return D[:, np.asarray(cols)].sum(axis=1).astype(float) / float(len(senmayo_names))

    sc = scores_of(idx122)
    ref_sc, deg_sc = sc[ref_rows_i], sc[deg_rows_i]
    a5 = bool(np.all(np.isfinite(sc)))

    m = float(np.mean(ref_sc))
    sd = float(np.std(ref_sc, ddof=1))
    z_all = (sc - m) / sd if sd > 0 else np.full_like(sc, np.nan)

    d_matched = cohens_d_pooled(deg_sc, ref_sc)
    p_welch = welch_one_sided(deg_sc, ref_sc)
    p_student = float(stats.ttest_ind(deg_sc, ref_sc, equal_var=True, alternative="greater").pvalue) if min(n_ref, n_deg) > 1 else float("nan")
    try:
        p_mw = float(stats.mannwhitneyu(deg_sc, ref_sc, alternative="greater").pvalue)
    except Exception:
        p_mw = float("nan")
    log(f"[{tag}] S-1: d_matched={d_matched:.4f}、片側 p(Welch)={p_welch:.3e}（参照 {n_ref} / 比較対象 {n_deg}）")

    # S-3: そろえる前（区間を付けない全参照・全比較対象）
    have = set(int(r) for r in read_rows)
    ref_all = [r for r in np.where(is_ref)[0] if int(r) in have]
    deg_all = [r for r in np.where(is_deg)[0] if int(r) in have]
    d_unmatched = cohens_d_pooled(sc[[row_index[int(r)] for r in deg_all]],
                                  sc[[row_index[int(r)] for r in ref_all]]) if (len(ref_all) > 1 and len(deg_all) > 1) else float("nan")
    ratio = (d_matched / d_unmatched) if (d_unmatched and not np.isnan(d_unmatched) and d_unmatched != 0) else float("nan")
    log(f"[{tag}] S-3: d_unmatched={d_unmatched:.4f}（{len(ref_all)} 対 {len(deg_all)} 細胞）、比={ratio}")

    # S-2: 頻度をそろえた無作為 122 遺伝子集合
    n_sets = int(n_random_override) if n_random_override else PLAN["S2"]["n_random_sets"]
    if len(pool_names) >= len(senmayo_names):
        sets, info = draw_random_sets(pool_freq, senmayo_freq, PLAN["S2"]["random_key"], n_sets,
                                      len(senmayo_names), log)
        ds = np.empty(len(sets), dtype=float)
        for i, s in enumerate(sets):
            cols = [cpos[name_to_cols[g][0]] for g in s]
            v = scores_of(cols)
            ds[i] = cohens_d_pooled(v[deg_rows_i], v[ref_rows_i])
        pct = float(np.percentile(ds, PLAN["S2"]["percentile"]))
        emp_p = float((int(np.sum(ds >= d_matched)) + 1) / (len(ds) + 1))
        s2 = {"run": True, "n_sets": len(sets), "pool_size": len(pool_names),
              "percentile_99": pct, "senmayo_d": float(d_matched), "pass": bool(d_matched > pct),
              "empirical_one_sided_p": emp_p, "random_d_mean": float(np.mean(ds)),
              "random_d_sd": float(np.std(ds, ddof=1)), "random_d_min": float(np.min(ds)),
              "random_d_max": float(np.max(ds)), "info": info}
        log(f"[{tag}] S-2: SenMayo d={d_matched:.4f}、無作為 {len(sets)} 組の 99 パーセンタイル={pct:.4f}、"
            f"合格={s2['pass']}、経験的な片側割合={emp_p:.4f}")
    else:
        s2 = {"run": False, "reason": "プールが 122 個未満", "pool_size": len(pool_names)}
        log(f"[{tag}] S-2: 省略（プール {len(pool_names)} 個）")

    # LODO（提供者を 1 人ずつ除く。判定は変えない）
    lodo = {"unit": "donor", "scope": "interval_reference_and_degenerate", "changes_decision": False,
            "donors": [], "d_min": None, "d_median": None, "d_max": None,
            "n_pass_effect_and_p": 0, "n_reversed": 0}
    ds_lodo = []
    for dn in sorted(set(donor[interval_cells].tolist())):
        keep = donor[interval_cells] != dn
        rr = interval_cells[keep & is_ref[interval_cells]]
        dd = interval_cells[keep & is_deg[interval_cells]]
        if len(rr) < 2 or len(dd) < 2:
            lodo["donors"].append({"donor": str(dn), "n_reference": int(len(rr)), "n_degenerate": int(len(dd)),
                                   "d": None, "p": None, "note": "2 未満の群があり計算しない"})
            continue
        rsc = sc[[row_index[int(r)] for r in rr]]
        dsc = sc[[row_index[int(r)] for r in dd]]
        dl, pl = cohens_d_pooled(dsc, rsc), welch_one_sided(dsc, rsc)
        ds_lodo.append(dl)
        if dl >= PLAN["S1"]["d_threshold"] and pl < PLAN["S1"]["p_threshold"]:
            lodo["n_pass_effect_and_p"] += 1
        if dl < 0:
            lodo["n_reversed"] += 1
        lodo["donors"].append({"donor": str(dn), "n_reference": int(len(rr)), "n_degenerate": int(len(dd)),
                               "d": dl, "p": pl})
    if ds_lodo:
        lodo["d_min"], lodo["d_median"], lodo["d_max"] = (float(np.min(ds_lodo)), float(np.median(ds_lodo)), float(np.max(ds_lodo)))
    log(f"[{tag}] LODO: 提供者 {len(lodo['donors'])} 人。d の範囲 {lodo['d_min']}〜{lodo['d_max']}、"
        f"効果量と p を満たした残し方 {lodo['n_pass_effect_and_p']}、向きが逆 {lodo['n_reversed']}")

    def dist(x):
        return {"n": int(len(x)), "mean": float(np.mean(x)) if len(x) else None,
                "median": float(np.median(x)) if len(x) else None,
                "sd_unbiased": float(np.std(x, ddof=1)) if len(x) > 1 else None,
                "min": float(np.min(x)) if len(x) else None, "max": float(np.max(x)) if len(x) else None,
                "frac_at_0": float(np.mean(x == 0)) if len(x) else None,
                "frac_at_1": float(np.mean(x == 1)) if len(x) else None}

    desc = {
        "score_distribution": {"reference": dist(ref_sc), "degenerate": dist(deg_sc)},
        "nFeature_RNA": {"reference": dist(nfeat[ref_int]), "degenerate": dist(nfeat[deg_int]),
                         "degenerate_over_reference_median_ratio":
                             float(np.median(nfeat[deg_int]) / np.median(nfeat[ref_int])) if n_ref and n_deg else None},
        "donor_cell_counts": {"degenerate": {str(d): int((donor[deg_int] == d).sum()) for d in donors_deg},
                              "reference": {str(d): int((donor[ref_int] == d).sum()) for d in donors_ref}},
    }

    # 感度区間（記述のみ。plan §6.2）
    sensitivity = {}
    sm = (nfeat >= sens_lo) & (nfeat <= sens_hi)
    sr = [r for r in np.where(is_ref & sm)[0] if int(r) in have]
    sd_ = [r for r in np.where(is_deg & sm)[0] if int(r) in have]
    if len(sr) > 1 and len(sd_) > 1:
        sensitivity[sens_label] = {
            "interval": [sens_lo, sens_hi], "reference_cells": len(sr), "degenerate_cells": len(sd_),
            "d": cohens_d_pooled(sc[[row_index[int(r)] for r in sd_]], sc[[row_index[int(r)] for r in sr]]),
            "p_welch_one_sided": welch_one_sided(sc[[row_index[int(r)] for r in sd_]], sc[[row_index[int(r)] for r in sr]]),
            "note": "記述のみ。判定には使わない（plan §6.2）"}

    scores_rows = [[str(cell_id[r]), tag, "reference" if is_ref[r] else "degenerate", str(donor[r]),
                    int(nfeat[r]), float(sc[row_index[int(r)]]), float(z_all[row_index[int(r)]])]
                   for r in interval_cells]

    return {
        "tag": tag, "cells": cells_report, "a3": a3, "a4": a4,
        "standardization": {"method": "reference_only_within_interval", "m": m, "sd": sd, "n_reference": int(n_ref)},
        "S1": {"d_matched": d_matched, "p_welch_one_sided": p_welch,
               "mean_reference": float(np.mean(ref_sc)), "mean_degenerate": float(np.mean(deg_sc)),
               "sd_reference": float(np.std(ref_sc, ddof=1)), "sd_degenerate": float(np.std(deg_sc, ddof=1)),
               "p_student_one_sided": p_student, "p_mannwhitney_one_sided": p_mw,
               "direction": "degenerate_higher" if d_matched > 0 else ("degenerate_lower" if d_matched < 0 else "equal"),
               "pass_effect_and_p": bool(d_matched >= PLAN["S1"]["d_threshold"] and p_welch < PLAN["S1"]["p_threshold"])},
        "S2": s2,
        "S3": {"d_matched": d_matched, "d_unmatched": d_unmatched, "ratio": ratio,
               "all_reference_cells": int(is_ref.sum()), "all_degenerate_cells": int(is_deg.sum()),
               "threshold_ratio": PLAN["S3"]["ratio_threshold"]},
        "T6_too_few_cells": bool(t6),
        "LODO": lodo, "descriptive": desc, "sensitivity_interval": sensitivity,
        "scores_rows": scores_rows, "scores_finite": a5,
    }


# ---------------------------------------------------------------------------
def main() -> int:
    t_start = _dt.datetime.now()
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true", help="A2 の小さな試し実行（判定ではない）")
    ap.add_argument("--cells-per-state", type=int, default=20, help="試し実行で各状態から使う細胞数")
    ap.add_argument("--sets", type=int, default=None, help="試し実行での無作為集合の数（既定 200）")
    args = ap.parse_args()

    log = Log()
    out_dir = EXP / ("results-smoke" if args.smoke else "results")
    out_dir.mkdir(parents=True, exist_ok=True)
    log(f"開始（{'試し実行（判定ではない）' if args.smoke else '本実行'}）。出力先 {out_dir.relative_to(ROOT)}")
    log(f"plan.md sha256 {sha256(EXP/'plan.md')[:16]}… / config.yaml sha256 {sha256(EXP/'config.yaml')[:16]}…")

    cfg = yaml.safe_load((EXP / "config.yaml").read_text())
    errs = check_config(cfg)
    if errs:
        for e in errs:
            log(f"config 不一致: {e}")
        raise SystemExit("config.yaml の decision: 節が plan.md と違う。中止して報告する")
    log("config.yaml の decision: 節は plan.md §4 と一致した")

    blocked = set(cfg["blocked_genes"]["symbols"])
    geneformer_symbols = set(pickle.load((ROOT / "models/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl").open("rb")).keys())
    log(f"Geneformer の遺伝子記号の辞書: {len(geneformer_symbols)} 個")

    senmayo_names = build_122(cfg, log)

    # 循環性（結果を見る前に、記号一致だけで数える。plan §5.4）
    t12 = list(csv.DictReader((ROOT / cfg["genes"]["positive_control_file"]).open()))
    t12_all = set(r["gene_symbol"] for r in t12)
    t12_a2 = set(r["gene_symbol"] for r in t12 if r["in_A2_var"] == "True")
    sen_set = set(senmayo_names)
    cand = set(cfg["genes"]["candidate_genes"])
    circ = {
        "note": "結果を見る前に遺伝子記号の一致だけで数えた（発現値・得点は扱わない。plan §5.4）",
        "senmayo_122": len(sen_set), "candidate_8": len(cand),
        "positive_control_all": len(t12_all), "positive_control_A2": len(t12_a2), "blocked": len(blocked),
        "senmayo_122_vs_candidate": len(sen_set & cand),
        "senmayo_122_vs_positive_control_A2": len(sen_set & t12_a2),
        "candidate_vs_positive_control_A2": len(cand & t12_a2),
        "blocked_vs_senmayo_122": len(blocked & sen_set),
        "blocked_vs_positive_control_A2": len(blocked & t12_a2),
    }
    log(f"循環性（記号一致のみ）: SenMayo 122 ∩ 候補 8={circ['senmayo_122_vs_candidate']}、"
        f"∩ 陽性対照(A2)={circ['senmayo_122_vs_positive_control_A2']}、候補 ∩ 陽性対照={circ['candidate_vs_positive_control_A2']}、"
        f"盲検 ∩ SenMayo={circ['blocked_vs_senmayo_122']}、盲検 ∩ 陽性対照={circ['blocked_vs_positive_control_A2']}")

    # A2（主）
    a2 = analyze_dataset("A2", cfg["data"]["A2"], senmayo_names, blocked, geneformer_symbols, log,
                         smoke=args.smoke,
                         limit_per_state=(args.cells_per_state if args.smoke else None),
                         n_random_override=(args.sets or 200 if args.smoke else None))

    # S-4（A1。判定に使わない）
    a1, s4 = None, None
    if not args.smoke:
        log("--- S-4（A1。判定に使わない）---")
        a1 = analyze_dataset("A1", cfg["data"]["A1"], senmayo_names, blocked, geneformer_symbols, log)
        s4 = {"interval": list(_A1_INTERVAL), "random_key_descriptive": PLAN["S4"]["random_key"],
              "d": a1["S1"]["d_matched"], "p_one_sided": a1["S1"]["p_welch_one_sided"],
              "percentile_99": a1["S2"].get("percentile_99") if a1["S2"].get("run") else None,
              "pass_senmayo_gt_percentile": a1["S2"].get("pass") if a1["S2"].get("run") else None,
              "cells": a1["cells"],
              "LODO_summary": {"d_min": a1["LODO"]["d_min"], "d_max": a1["LODO"]["d_max"],
                               "n_pass_effect_and_p": a1["LODO"]["n_pass_effect_and_p"], "n_reversed": a1["LODO"]["n_reversed"]},
              "note": ("A1 では細胞老化の遺伝子の検出率が 1% 前後以下である（仮説 §3.2）。"
                       "S-4 が不合格でも段階 1 の判定には影響しない。判定には使わない。")}
        log(f"S-4: d={s4['d']:.4f}、片側 p={s4['p_one_sided']:.3e}、99 パーセンタイル={s4['percentile_99']}、"
            f"合格={s4['pass_senmayo_gt_percentile']}（判定には使わない）")
    else:
        log("試し実行: A1（S-4）は省略（plan §7: A2 だけを対象にした小さな試し実行）")

    # T-1〜T-6（plan §4.5）と評価順
    d_matched, p_w = a2["S1"]["d_matched"], a2["S1"]["p_welch_one_sided"]
    d_un, ratio = a2["S3"]["d_unmatched"], a2["S3"]["ratio"]
    s2_pass = bool(a2["S2"].get("pass")) if a2["S2"].get("run") else False
    T = {
        "T5_indicator_not_viable": bool(a2["a4"]["usable_in_var_and_dict"] < PLAN["T5_min_gene_count"]),
        "T6_too_few_cells": bool(a2["T6_too_few_cells"]),
        "T2_direction_reversed": bool((not np.isnan(d_matched)) and d_matched < 0),
        "T1_effect_below_0.5": bool((np.isnan(d_matched) or d_matched < PLAN["S1"]["d_threshold"])
                                    or (np.isnan(p_w) or p_w >= PLAN["S1"]["p_threshold"])),
        "T3_not_specific": bool(not s2_pass),
        "T4_depth_explains": bool((not np.isnan(d_un)) and d_un > 0 and (not np.isnan(ratio))
                                  and ratio < PLAN["S3"]["ratio_threshold"]),
    }
    triggered = "all_pass"
    for short in PLAN["evaluation_order"]:
        if T[T_KEY[short]]:
            triggered = f"T-{short[1]}"
            break
    log(f"評価順（T-5→T-6→T-2→T-1→T-3→T-4）で最初に当たった条件: {triggered}")

    dur = (_dt.datetime.now() - t_start).total_seconds()
    criteria = {
        "experiment_id": cfg["experiment_id"], "hypothesis_id": cfg["hypothesis_id"], "stage": 1,
        "smoke": bool(args.smoke),
        "generated_at": _dt.datetime.now().astimezone().isoformat(),
        "duration_seconds": round(dur, 1),
        "note": ("plan.md §4 の判定規則をそのまま適用した値。支持／否定／判定不能の確定は PI（ゲート G3）。"
                 + ("【試し実行。この値は判定ではない】" if args.smoke else "")),
        "hashes": {"plan.md": sha256(EXP / "plan.md"), "config.yaml": sha256(EXP / "config.yaml"),
                   "script": sha256(Path(__file__).resolve())},
        "gene_list": {"n": len(senmayo_names), "source": cfg["genes"]["senmayo_file"],
                      "excluded_not_in_A2_count": len(cfg["genes"]["exclude_not_in_A2"]),
                      "excluded_blind_count": len(cfg["genes"]["exclude_blind_from_list"]),
                      "n_in_A1_var": None},
        "blocked_genes": {"count": len(blocked), "names_omitted": True,
                          "note": "盲検中のため名前を値とともに出力に書かない（plan §3.5・§7 A-1）"},
        "circularity": circ,
        "A2": {k: v for k, v in a2.items() if k != "scores_rows"},
        "S4_A1": s4,
        "A1_secondary_raw": ({k: v for k, v in a1.items() if k != "scores_rows"} if a1 else None),
        "T_conditions": T,
        "decision_proposal": {"evaluation_order": PLAN["evaluation_order"], "first_triggered": triggered,
                              "S1_pass": a2["S1"]["pass_effect_and_p"], "S2_pass": s2_pass,
                              "S3_not_T4": (not T["T4_depth_explains"]),
                              "all_of_S1_S2_S3_pass": bool(a2["S1"]["pass_effect_and_p"] and s2_pass
                                                           and not T["T4_depth_explains"]),
                              "confirmed_by": "PI（ゲート G3）"},
        "A_checks": {
            "A-2_state_labels": "合格（指定した参照・比較対象のラベルがデータに実在）",
            "A-3_counts_integer_like": a2["a3"]["counts_integer_like"],
            "A-4_genes": a2["a4"],
            "A-5_scores_finite": a2["scores_finite"],
            "A-6_outputs_written": True,
            "A-1_blocked_symbol_occurrences": None,
            "A-1_pass": None,
            "A-7_duration_seconds": round(dur, 1),
        },
    }

    # 得点一覧
    all_rows = list(a2["scores_rows"]) + (list(a1["scores_rows"]) if a1 else [])
    scores_csv = out_dir / "scores_012.csv"
    with scores_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cell_id", "dataset", "state", "donor", "nFeature_RNA", "score", "z"])
        for row in all_rows:
            w.writerow(row)
    log(f"得点一覧を書いた: {scores_csv.relative_to(ROOT)}（{len(all_rows)} 行）")

    # A-1: 盲検 5 遺伝子が、出力・ログのどこにも 0 件か（機械的に確認）
    criteria_txt = json.dumps(criteria, ensure_ascii=False, indent=2)
    csv_txt = scores_csv.read_text()
    log_txt = log.text()
    hits = 0
    for sym in sorted(blocked):
        for blob in (criteria_txt, csv_txt, log_txt):
            hits += len(re.findall(re.escape(sym), blob))
    criteria["A_checks"]["A-1_blocked_symbol_occurrences"] = int(hits)
    criteria["A_checks"]["A-1_pass"] = bool(hits == 0)
    if hits != 0:
        log(f"A-1 不合格: 盲検の遺伝子記号が出力に {hits} 件現れた。中止して報告する")
        raise SystemExit("A-1: 盲検の混入。中止する")
    log("A-1 合格: 盲検 5 遺伝子の記号は、得点一覧・出力・ログのどこにも現れない")

    criteria_json = out_dir / "criteria_012.json"
    criteria_json.write_text(json.dumps(criteria, ensure_ascii=False, indent=2))
    log(f"判定値の JSON を書いた: {criteria_json.relative_to(ROOT)}")

    (out_dir / "README.md").write_text(_readme_text(criteria, scores_csv.name, args.smoke))
    log(f"説明を書いた: {(out_dir / 'README.md').relative_to(ROOT)}")
    log(f"終了（終了コード 0）。所要 {dur:.1f} 秒")
    (out_dir / "run.log").write_text(log.text())
    return 0


def _readme_text(criteria: dict, scores_name: str, smoke: bool) -> str:
    a2 = criteria["A2"]
    s1, s2, s3 = a2["S1"], (a2["S2"] if a2["S2"].get("run") else {}), a2["S3"]
    L = a2["LODO"]
    lines = [
        f"# EXP-20261007-002 段階 1 — 出力（{'試し実行。判定ではない' if smoke else '本実行'}）",
        "",
        f"- 生成: {criteria['generated_at']}、所要 {criteria['duration_seconds']} 秒",
        f"- 仮説: {criteria['hypothesis_id']} 段階 1。**摂動なし・Geneformer 不使用**（GPU 0 時間）",
        "- **判定（支持／否定／判定不能）の確定は PI（ゲート G3）。ここは値と提案まで。**",
        "",
        "## 何を出したか",
        "",
        "| ファイル | 中身 |",
        "|---|---|",
        "| `criteria_012.json` | 判定値と、判定に使ったすべての数値（A2 の主判定、特異性、深さ依存、反証条件、"
        "提供者を 1 人ずつ除く確認、循環性、S-4） |",
        f"| `{scores_name}` | 区間の細胞ごとの SenMayo 得点（細胞 ID・データ・状態・提供者・`nFeature_RNA`・得点・標準得点） |",
        "| `README.md` | この説明 |",
        "| `run.log` | 実行の記録 |",
        "",
        "## 主な数値（A2。区間 1,000〜2,000）",
        "",
        f"- 細胞: 参照 {a2['cells']['reference_cells']} / 比較対象 {a2['cells']['degenerate_cells']}"
        f"（提供者 {a2['cells']['reference_donors']} / {a2['cells']['degenerate_donors']}）。"
        f"仮説 §2.2 の実測との一致: {a2['cells']['matches_expected']}",
        f"- 標準化（参照だけ。案 B）: m={a2['standardization']['m']:.6f}, sd={a2['standardization']['sd']:.6f}, "
        f"n={a2['standardization']['n_reference']}",
        f"- **主判定 S-1**: Cohen's d={s1['d_matched']:.4f}、片側 p（Welch）={s1['p_welch_one_sided']:.3e}、"
        f"向き={s1['direction']}、効果量と p の合格={s1['pass_effect_and_p']}",
        f"  - 参考: Student 片側 p={s1['p_student_one_sided']:.3e}、Mann-Whitney 片側 p={s1['p_mannwhitney_one_sided']:.3e}（判定は Welch）",
        f"- 特異性 S-2: SenMayo の d={s2.get('senmayo_d')}、無作為 {s2.get('n_sets')} 組の 99 パーセンタイル="
        f"{s2.get('percentile_99')}、合格={s2.get('pass')}、経験的な片側割合={s2.get('empirical_one_sided_p')}、"
        f"プール {s2.get('pool_size')} 個",
        f"- 深さ依存 S-3: そろえる前 d={s3['d_unmatched']:.4f}（{s3['all_reference_cells']} 対 {s3['all_degenerate_cells']} 細胞）、"
        f"そろえた後 d={s3['d_matched']:.4f}、比={s3['ratio']}",
        f"- 提供者を 1 人ずつ除く確認（判定は変えない）: d の {L['d_min']} 〜 {L['d_max']}、"
        f"効果量と p を満たした残し方 {L['n_pass_effect_and_p']}、向きが逆 {L['n_reversed']}",
        f"- 循環性（記号一致のみ。結果を見る前に数えた）: SenMayo 122 ∩ 候補 8 = {criteria['circularity']['senmayo_122_vs_candidate']}、"
        f"∩ 陽性対照(A2) = {criteria['circularity']['senmayo_122_vs_positive_control_A2']}",
        "",
        "## 反証条件（plan §4.5）",
        "",
        "| 条件 | 当てはまるか |",
        "|---|---|",
    ]
    for k, v in criteria["T_conditions"].items():
        lines.append(f"| {k} | {v} |")
    lines += [
        "",
        f"- 評価順（T-5→T-6→T-2→T-1→T-3→T-4）で最初に当たった条件: **{criteria['decision_proposal']['first_triggered']}**",
        f"- S-1・S-2・S-3 のすべて合格の提案: {criteria['decision_proposal']['all_of_S1_S2_S3_pass']}",
        "",
        "## 副次（A1。判定に使わない）",
        "",
    ]
    if criteria.get("S4_A1"):
        q = criteria["S4_A1"]
        lines += [f"- A1（区間 {q['interval']}）: d={q['d']:.4f}、片側 p={q['p_one_sided']:.3e}、"
                  f"99 パーセンタイル={q['percentile_99']}、合格={q['pass_senmayo_gt_percentile']}",
                  f"- **{q['note']}**"]
    else:
        lines.append("- 試し実行のため省略。")
    lines += [
        "",
        "## 注意",
        "",
        "- in silico の計算である。「示唆される」までで、実験的事実としては書かない。",
        "- 盲検中の 5 遺伝子の名前・値は、出力・ログのどこにも書かない（plan §3.5・§7 A-1）。"
        "循環性の重なりは個数だけを記録する。",
        "- 判定基準・区間・乱数の鍵は plan.md と config.yaml に固定した値のまま。実行時に変えていない。",
        "",
    ]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    sys.exit(main())

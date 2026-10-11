#!/usr/bin/env python3
"""EXP-20261011-001 — 事前確認（precheck）。摂動は 1 回も使わない。

やること（HYP-20261011-01 §8・§9.3 を結果を見る前に数えて固定する）:
1. トークン化済み入力を作業用に複製（A2・A1）し、sha256 を記録する（guardrail 2026-09-26）
2. 循環性: 判定に関わるすべての集合の重なりを総当たりで数える（A2 は確定値の再確認、A1 は本実験で初めて数える）
3. A1 で新しく摂動する退化の陽性対照 30 個が、Geneformer の辞書にあり、A1 の var にあり、
   検出率の中央値が 0.06 以上であるかを確かめる（§9.3。**1 つでも欠けたら個数を記録して PI に諮る**）
4. 盲検 5 個（TP53・CDKN1A・LTBR・SPP1・LCN2）がどの集合にも入っていないことを確かめる（値は数えない）

出力: experiments/EXP-20261011-001/work/precheck_014.json（事前登録のハッシュに含める）
"""
from __future__ import annotations

import hashlib
import json
import pickle
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
EXP = ROOT / "experiments/EXP-20261011-001"
A2_SRC = ROOT / "experiments/EXP-20261007-001"
A1_SRC = ROOT / "experiments/EXP-20261004-001"
EXP003 = ROOT / "experiments/EXP-20261007-003"

sys.path.insert(0, str(A2_SRC / "work"))
try:
    from geneformer import TOKEN_DICTIONARY_FILE
except Exception:  # 辞書は後で遅延インポートする
    TOKEN_DICTIONARY_FILE = None

BLOCKED = {
    "ENSG00000111321": "LTBR", "ENSG00000118785": "SPP1", "ENSG00000124762": "CDKN1A",
    "ENSG00000141510": "TP53", "ENSG00000148346": "LCN2",
}

# 退化の陽性対照 30 個（EXP-20261007-003 config.yaml の positive_control_degenerate_30 から写した。不変）
DEG30 = {
    "ENSG00000085662": "AKR1B1", "ENSG00000135931": "ARMC9", "ENSG00000130724": "CHMP2A",
    "ENSG00000103811": "CTSH", "ENSG00000205795": "CYS1", "ENSG00000169738": "DCXR",
    "ENSG00000125868": "DSTN", "ENSG00000125971": "DYNLRB1", "ENSG00000173812": "EIF1",
    "ENSG00000156976": "EIF4A2", "ENSG00000127951": "FGL2", "ENSG00000206503": "HLA-A",
    "ENSG00000004776": "HSPB6", "ENSG00000142798": "HSPG2", "ENSG00000115461": "IGFBP5",
    "ENSG00000130522": "JUND", "ENSG00000169223": "LMAN2", "ENSG00000140545": "MFGE8",
    "ENSG00000111341": "MGP", "ENSG00000133392": "MYH11", "ENSG00000151366": "NDUFC2",
    "ENSG00000110717": "NDUFS8", "ENSG00000138735": "PDE5A", "ENSG00000123349": "PFDN5",
    "ENSG00000100911": "PSME2", "ENSG00000063177": "RPL18", "ENSG00000124614": "RPS10",
    "ENSG00000163191": "S100A11", "ENSG00000109610": "SOD3", "ENSG00000006327": "TNFRSF12A",
}


def sha256f(p: Path) -> str | None:
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None


def clone_tokenized(name: str, src: Path) -> dict:
    dst = EXP / f"work/tokenized-{name}.dataset"
    dst.mkdir(parents=True, exist_ok=True)
    names = sorted([p.name for p in src.glob("data-*.arrow")]) + ["dataset_info.json", "state.json"]
    rec: dict = {"複製元": str(src.relative_to(ROOT)), "複製先": str(dst.relative_to(ROOT)), "ファイル": {}}
    ok = True
    for n in names:
        s, d = src / n, dst / n
        if not d.exists() or sha256f(d) != sha256f(s):
            r = subprocess.run(["cp", "-c", str(s), str(d)], capture_output=True, text=True)
            if r.returncode != 0:                     # APFS クローンが不可なら通常コピー
                subprocess.run(["cp", str(s), str(d)], check=True)
        a, b = sha256f(s), sha256f(d)
        rec["ファイル"][n] = {"sha256_元": a, "sha256_複製": b, "一致": a == b}
        ok &= (a == b)
    rec["pass"] = bool(ok)
    return rec


def set_map(g: dict, key: str) -> dict[str, str]:
    v = g[key]
    return dict(v) if isinstance(v, dict) else {r["ensembl"]: r["symbol"] for r in v}


def main() -> None:
    g11 = json.loads((A2_SRC / "work/genes_011.json").read_text())
    cand = dict(g11["候補"])
    pos71 = dict(g11["陽性対照71"])
    hk3 = dict(g11["恒常発現3"])
    neu7 = dict(g11["中立7"])
    ctrl500 = {r["ensembl"]: r["symbol"] for r in g11["対照500"]}
    # すべて Ensembl ID の集合にそろえる（記号の集合と混ぜない。2026-10-11 に一度間違えた）
    cand_ens = set(cand.values())
    pos71_ens = set(pos71.values())
    hk3_ens = set(hk3.values())
    neu7_ens = set(neu7.values())
    poolA2 = set(ctrl500) | hk3_ens | neu7_ens
    sym_of = {**{v: k for k, v in cand.items()}, **{v: k for k, v in pos71.items()},
              **{v: k for k, v in hk3.items()}, **{v: k for k, v in neu7.items()},
              **ctrl500, **DEG30}

    # A1 の対照プール 549 と測る遺伝子 626（EXP-20261004-001 config.yaml の controls から機械的に再構成）
    cfgA1 = yaml_safe_load = __import__("yaml").safe_load((A1_SRC / "config.yaml").read_text())
    ctl = cfgA1["controls"]
    rand538 = [x for x in ctl["negative_sampler"]["genes"]][: ctl["negative_random_n"]]
    neg_fixed = list(ctl["negative_fixed"])
    extra = [x for grp in (ctl.get("extra_fixed") or []) for x in grp["genes"]]
    poolA1 = set(rand538) | set(neg_fixed) | set(extra)
    targets8 = list(cfgA1["perturbation"]["genes_to_perturb"])
    pos71A1 = list(ctl["positive"])
    measured626 = set(targets8) | set(pos71A1) | poolA1

    def overlap(desc: str, a: set[str], b: set[str]) -> dict:
        inter = sorted(a & b)
        return {"組": desc, "件数": len(inter),
                "遺伝子": [sym_of.get(x, x) for x in inter]}

    overlaps = [
        overlap("退化30 ∩ 対照プールA2(510)", set(DEG30), poolA2),
        overlap("退化30 ∩ 活性化PC71", set(DEG30), pos71_ens),
        overlap("退化30 ∩ 候補8", set(DEG30), cand_ens),
        overlap("退化30 ∩ 恒常発現3", set(DEG30), hk3_ens),
        overlap("退化30 ∩ 中立7", set(DEG30), neu7_ens),
        overlap("候補8 ∩ 活性化PC71", cand_ens, pos71_ens),
        overlap("退化30 ∩ 対照プールA1(549)", set(DEG30), poolA1),
        overlap("退化30 ∩ A1測定626", set(DEG30), measured626),
        overlap("候補8 ∩ 対照プールA1(549)", cand_ens, poolA1),
        overlap("活性化PC71 ∩ 対照プールA1(549)", set(pos71A1), poolA1),
        overlap("退化30 ∩ 盲検5", set(DEG30), set(BLOCKED)),
        overlap("候補8 ∩ 盲検5", cand_ens, set(BLOCKED)),
        overlap("対照プールA2 ∩ 盲検5", poolA2, set(BLOCKED)),
        overlap("対照プールA1 ∩ 盲検5", poolA1, set(BLOCKED)),
        overlap("活性化PC71 ∩ 盲検5", pos71_ens, set(BLOCKED)),
    ]

    # 退化 30 の辞書・var・検出率（A1）
    tokd = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb")) if TOKEN_DICTIONARY_FILE else {}
    import anndata as ad
    import numpy as np
    A1 = ad.read_h5ad(str(ROOT / "experiments/EXP-20260929-001/work/input.h5ad"), backed="r")
    var_idx = set(A1.var.index.astype(str))
    mask = A1.obs["subclass_l2"].astype(str).isin(["FIB", "M-FIB", "dFIB"]).to_numpy()
    det_cols = {ens: int(np.where(A1.var.index.astype(str) == ens)[0][0]) for ens in DEG30 if ens in var_idx}
    deg_ok = []
    n_missing_dict = n_missing_var = n_below_detect = 0
    per_gene = {}
    for ens, sym in DEG30.items():
        in_dict = ens in tokd
        in_var = ens in var_idx
        det = 0.0
        if in_var:
            col = A1[:, det_cols[ens]].X
            vals = col.toarray().ravel() if hasattr(col, "toarray") else np.asarray(col).ravel()
            det = float((vals[mask] > 0).mean())
        if not in_dict:
            n_missing_dict += 1
        if not in_var:
            n_missing_var += 1
        if det < 0.06:
            n_below_detect += 1
        deg_ok.append((ens, sym, in_dict, in_var, det))
        per_gene[sym] = {"辞書": in_dict, "A1 var": in_var, "検出率(start+goal)": round(det, 4),
                         "0.06 未満": bool(det < 0.06)}
    dets = [d for _, _, _, _, d in deg_ok if d >= 0]
    det_med = float(np.median(dets)) if dets else None

    # A2 の var にもあるか（退化 30 は EXP-001 の保存済み 587 に含まれる前提の再確認）
    A2in = __import__("anndata").read_h5ad(str(ROOT / "outputs/design/2026-10-07-A2-feasibility/work/input.h5ad"), backed="r")
    a2var = set(A2in.var.index.astype(str))
    n_missing_a2 = sum(1 for ens in DEG30 if ens not in a2var)

    out = {
        "実験": "EXP-20261011-001", "仮説": "HYP-20261011-01",
        "注": "摂動は 1 回も使っていない。HYP-20261011-01 §8（循環性）と §9.3（A1 の退化の陽性対照）の結果を見る前の数え上げ。",
        "1. トークン化済み入力の複製": {
            "A2": clone_tokenized("a2", A2_SRC / "work/tokenized.dataset"),
            "A1": clone_tokenized("a1", A1_SRC / "work/tokenized.dataset"),
        },
        "2. 集合の大きさ": {
            "候補8": len(cand), "活性化PC71": len(pos71), "退化PC30": len(DEG30),
            "恒常発現3": len(hk3), "中立7": len(neu7), "対照500": len(ctrl500),
            "対照プールA2(510)": len(poolA2), "対照プールA1(549)": len(poolA1),
            "A1 測定626": len(measured626),
        },
        "3. 重なり（総当たり）": overlaps,
        "3b. 固定値との一致": {
            "退化30 ∩ 対照プールA2 = 23": any(o["組"] == "退化30 ∩ 対照プールA2(510)" and o["件数"] == 23 for o in overlaps),
            "退化30 ∩ 活性化PC71 = 7": any(o["組"] == "退化30 ∩ 活性化PC71" and o["件数"] == 7 for o in overlaps),
            "退化30 ∩ 候補8 = 0": any(o["組"] == "退化30 ∩ 候補8" and o["件数"] == 0 for o in overlaps),
            "候補8 ∩ 活性化PC71 = 2": any(o["組"] == "候補8 ∩ 活性化PC71" and o["件数"] == 2 for o in overlaps),
            "盲検との重なり = 0（全組）": all(o["件数"] == 0 for o in overlaps if "盲検" in o["組"]),
        },
        "3c. 仮説 §3.3 の訂正（結果を見る前に発見した事実）": {
            "仮説 §3.3 の記述": "A1 の保存済み 626 遺伝子には、退化した状態の陽性対照 30 個が 1 個も入っていない",
            "実測": "10 個が既に A1 の 626 に含まれる（活性化 PC71 と重なる 7 個＋対照プール A1 と重なる 3 個）",
            "内訳": {
                "活性化 PC71 と重なる（新規摂動でも A1 保存済みでも測れる）": ["MGP", "IGFBP5", "DSTN", "ARMC9", "PDE5A", "MFGE8", "HSPG2"],
                "対照プール A1(549) と重なる": ["CTSH", "PFDN5", "CYS1"],
            },
            "設計への反映": "新規摂動はタスク指示どおり 30 個すべてを回す。既存の 10 個は（同じ機体・同じ bf16 のはずなので）"
                          "機械的な一致性の確認に使い、P-1 (ii) の退化の陽性対照 30 個の値は新規摂動の結果を主とする。",
        },
        "4. 退化30 の A1 での利用可能性（§9.3）": {
            "辞書に無い": n_missing_dict, "A1 var に無い": n_missing_var,
            "A2 var に無い（再確認）": n_missing_a2,
            "検出率の中央値": round(det_med, 4) if det_med is not None else None,
            "§9.3 の条件": "30 個すべてが辞書・A1 var にあり、検出率の中央値が 0.06 以上",
            "検出率 0.06 未満の個数（記述。止める条件ではない）": n_below_detect,
            "検出率 0.06 未満の遺伝子（記述）": [sym for _, sym, _, _, d in deg_ok if d < 0.06],
            "遺伝子ごと": per_gene,
            "注": "1 つでも辞書・var から欠けたら個数を記録して PI に諮る（§9.3）。"
                  "過剰発現型は未発現の遺伝子でも先頭に挿入できる（Geneformer の overexpress_tokens）ため、"
                  "検出率が低い遺伝子も摂動の対象にはできる。検出率は A1 の出発点（FIB+M-FIB）+ 退化（dFIB）で測った。",
        },
        "5. 判定": None,
    }
    cond_ok = (n_missing_dict == 0 and n_missing_var == 0 and n_missing_a2 == 0
               and (det_med is not None and det_med >= 0.06)
               and all(o["件数"] == 0 for o in overlaps if "盲検" in o["組"]))
    out["5. 判定"] = ("§9.3 の条件は満たす（辞書・var の欠け 0 個、検出率の中央値 "
                      + (f"{det_med:.4f}" if det_med is not None else "nan")
                      + " ≥ 0.06）。着手前チェックと G2 の判定は別途行う。"
                      if cond_ok else "**§9.3 の条件を満たさない。個数を記録して PI に諮る。**")
    (EXP / "work").mkdir(parents=True, exist_ok=True)
    (EXP / "work/precheck_014.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
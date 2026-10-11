#!/usr/bin/env python3
"""EXP-20260926-004 (a): VCAM1/HAVCR1 と同程度に aPT_CKD_TLShigh で上昇する「非ドライバー」比較遺伝子（DE マッチ）の選定。

**発現データのみ**（ISP / perturbation の結果は一切読まない）。入力は negctrl_candidates.py が
EXP-20260926-001 の input.h5ad（生カウント）から作った全遺伝子の統計表 genome_stats.csv（読み取りのみ）:
  - log2fc_start_vs_goal: donor pseudo-bulk（≥ 20 細胞の donor, log2(CPM+1) の donor 平均の差）start vs goal
  - welch_q_goal: donor 水準の Welch t の BH q
  - tok_cells_start: start（TLS-high aPT 7,533 細胞）で Geneformer 入力（4,094 トークン）に入る細胞数
  - relrank_start: start の入力内の相対順位（0 = 先頭）の平均
  - donors_ge25_start: 25 細胞以上で検出される start donor 数
規則（plan.md §A.2 に事前登録）:
  E1 Geneformer 辞書に収載、E2 log2FC ≥ 2 かつ q < 0.05（上昇）、E3 tok_cells_start ≥ 440（HAVCR1 659 の 1/1.5）
  かつ donors_ge25_start ≥ 8、E4 自動除外（EXCLUDE: negctrl の既知リスト + 本実験の追加カテゴリ）、E5 文献レビューで
  除外（MANUAL_REVIEW, 理由と不確実性を記録）
  距離 d_r(g) = max(|Δlog2FC|/1.0, |ln(cells_g/cells_r)|/ln 1.5, |Δrelrank|/0.10)（r = VCAM1, HAVCR1）
  各 r について d の小さい順に、E4/E5 で除外されない遺伝子を 5 つ（重複なし）→ 計 10 遺伝子
出力: experiments/EXP-20260926-004/de_matched/{candidates.csv, selected.csv, selection.json}
usage: .venv/bin/python experiments/EXP-20260926-004/analysis/de_matched_select.py
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
GS = ROOT / "outputs/literature/2026-09-26-negative-controls/genome_stats.csv"
OUT = ROOT / "experiments/EXP-20260926-004/de_matched"
N_PER_REF = 5
REFS = ["VCAM1", "HAVCR1"]
TOL = {"log2fc": 1.0, "lncells": np.log(1.5), "relrank": 0.10}
E2_MIN_LOG2FC, E2_MAX_Q = 2.0, 0.05
E3_MIN_CELLS, E3_MIN_DONORS = 440, 8
BLIND = {"TP53", "CDKN1A", "LTBR", "SPP1", "LCN2"}

# negctrl_candidates.py の既知除外リスト（傷害・aPT マーカー、TLS/免疫、解離ストレス、性染色体、細胞周期、ストレス応答）
_spec = importlib.util.spec_from_file_location("negctrl", ROOT / "scripts/analysis/negctrl_candidates.py")
_neg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_neg)  # main() は呼ばない（読み込むだけ）

# 本実験で追加する除外カテゴリ（候補を見る前に決めたカテゴリ。個々の遺伝子は文献知識から）
EXTRA_EXCLUDE = {
    # 傷害・aPT/failed-repair・老化の既知マーカー／メディエーター
    "TNFRSF12A", "CD44", "TIMP1", "TIMP2", "IGFBP7", "GDF15", "CDKN2B", "CDKN1B", "IL32", "TACSTD2", "KRT20",
    "CDH6", "PLA2R1", "GLB1", "LMNB1", "NNMT", "SOX4", "PAX2", "PAX8", "HNF1B", "NT5E",
    # 補体・IFN・炎症
    "CFB", "CFH", "SERPING1", "C1R", "C1S", "C4A", "C4B", "C7", "CD74", "CIITA", "B2M", "TNFSF10", "TNFAIP2",
    "TNFAIP3", "ICAM1", "STAT1", "IRF1", "IRF7",
    # 修復シグナル経路（HGF/MET・EGFR・Wnt・Notch・Hedgehog・TGF-β・YAP）
    "MET", "HGF", "EGFR", "AREG", "HBEGF", "EREG", "HES1", "HEY1", "HEYL", "SHH", "IHH", "GLI1", "GLI2",
    "YAP1", "WWTR1", "CTGF", "CCN2", "CYR61", "CCN1", "AXIN2", "LGR5", "ITGAV", "ITGB6",
    # ECM・線維化
    "FN1", "TNC", "POSTN", "SPARC", "VCAN", "THBS1", "LUM", "DCN", "ELN", "FBN1",
    # 酸化ストレス応答
    "PON2", "NQO1", "GCLM", "GCLC", "TXNRD1", "NFE2L2", "NFE2L3", "SRXN1", "GPX2", "AKR1C1", "AKR1B10",
    # 代謝ストレス・低酸素
    "VEGFA", "HIF1A", "EPAS1", "CA9", "NDRG1", "BNIP3", "PDK1",
}
EXTRA_PREFIX = ("WNT", "FZD", "NOTCH", "JAG", "DLL", "TGFB", "SMAD", "COL", "MMP", "LOX", "IFI", "IFIT", "ISG",
                "OAS", "GBP", "MX1", "MX2", "CXCL", "CCL", "IL1", "IL6", "IL1", "TNF", "SERPINE", "KRT")

# 文献レビュー（距離順に上から。E4 を通った遺伝子のうち選定までに現れたものは全て記録する）
# decision: include / exclude。uncertainty: 低 / 中 / 高（「PT 傷害・修復・老化・TLS の駆動因子の証拠がない」判断の確からしさ）
MANUAL_REVIEW = {
    # VCAM1 に近い帯
    "PPFIBP1": ("include", "liprin-β1（接着斑・細胞遊走の足場）。PT 傷害・修復・老化の駆動因子とする報告は把握していない", "中"),
    "AKAP12": ("include", "gravin（PKA/PKC の足場）。傷害 PT で上昇するとの記載はありうるが、駆動因子としての因果証拠は把握していない", "高"),
    "SELENOM": ("include", "セレノプロテイン M（ER の酸化還元）。腎傷害での役割の報告は把握していない", "中"),
    "DMKN": ("include", "dermokine（表皮の分泌タンパク質）。腎での機能報告は把握していない", "中"),
    "MACC1": ("exclude", "MET の転写調節因子。HGF/MET は PT 再生の修復経路（経路単位の保守的な除外）", "高"),
    "SCIN": ("include", "scinderin（アクチン切断）。腎傷害での役割の報告は把握していない", "中"),
    # HAVCR1 に近い帯
    "ANKRD2": ("include", "筋のストレス応答性 ankyrin repeat タンパク質。腎 PT での機能報告は把握していない（検出はアンビエント・ダブレットの可能性も）", "中"),
    "C1QTNF3": ("exclude", "CTRP3（アディポカイン）。糖尿病性腎症・腎線維化モデルで保護的との報告がある（記憶に基づく, 未検証）", "高"),
    "SCARA3": ("exclude", "別名 CSR1（cellular stress response 1）。酸化ストレスで誘導 → ストレス応答カテゴリとして除外", "中"),
    "WNT2B": ("exclude", "Wnt リガンド（Wnt/β-catenin は傷害・修復・線維化の主要経路）。EXTRA_PREFIX でも除外", "低"),
    "ARRDC4": ("include", "α-arrestin（TXNIP の類縁, 糖代謝応答）。腎傷害の駆動因子とする報告は把握していない", "中"),
    "NFE2L3": ("exclude", "NRF3（CNC-bZIP, NRF2 の類縁）。酸化ストレス応答カテゴリとして保守的に除外。EXTRA_EXCLUDE でも除外", "高"),
    "IDUA": ("include", "α-L-iduronidase（リソソーム酵素）。腎傷害の駆動因子とする報告は把握していない", "低"),
    "SLITRK4": ("include", "神経のシナプス接着分子。腎での機能報告は把握していない", "低"),
    "ADAMTS16": ("include", "血圧関連座位（ラットで腎に発現）。PT 傷害・修復の駆動因子とする報告は把握していない", "中"),
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(GS)
    ref = d.set_index("symbol").loc[REFS]
    excl_auto = d.symbol.isin(_neg.EXCLUDE_SYMBOLS | EXTRA_EXCLUDE | BLIND | set(REFS)) | \
        d.symbol.str.startswith(_neg.EXCLUDE_PREFIX + EXTRA_PREFIX) | d.known_excluded.astype(bool)
    d["E1"] = d.in_gf_vocab.astype(bool)
    d["E2"] = (d.log2fc_start_vs_goal >= E2_MIN_LOG2FC) & (d.welch_q_goal < E2_MAX_Q)
    d["E3"] = (d.tok_cells_start >= E3_MIN_CELLS) & (d.donors_ge25_start >= E3_MIN_DONORS)
    d["E4_auto_excluded"] = excl_auto
    pool = d[d.E1 & d.E2 & d.E3].copy()
    rows, chosen = [], []
    for r in REFS:
        fc, tc, rr = ref.loc[r, ["log2fc_start_vs_goal", "tok_cells_start", "relrank_start"]]
        x = pool.copy()
        x["d_fc"] = (x.log2fc_start_vs_goal - fc).abs() / TOL["log2fc"]
        x["d_cells"] = np.log(x.tok_cells_start / tc).abs() / TOL["lncells"]
        x["d_rank"] = (x.relrank_start - rr).abs() / TOL["relrank"]
        x["dist"] = x[["d_fc", "d_cells", "d_rank"]].max(1)
        x = x.sort_values("dist")
        n = 0
        for _, g in x.iterrows():
            if n >= N_PER_REF:
                break
            s = g.symbol
            if s in chosen:
                continue
            if g.E4_auto_excluded:
                dec, why, unc = "exclude_auto", "既知リスト / 接頭辞 / 陽性・盲検対象（E4）", "—"
            else:
                assert s in MANUAL_REVIEW, f"{s} は文献レビュー未記入（選定前に記入すること）"
                dec, why, unc = MANUAL_REVIEW[s]
                if dec == "include":
                    chosen.append(s); n += 1
            rows.append({"ref": r, "rank_by_dist": len([q for q in rows if q["ref"] == r]) + 1, "symbol": s,
                         "ensembl_id": g.ensembl_id, "decision": dec, "reason": why, "uncertainty": unc,
                         **{k: g[k] for k in ["dist", "d_fc", "d_cells", "d_rank", "log2fc_start_vs_goal", "welch_q_goal",
                                              "tok_cells_start", "det_start", "det_goal", "relrank_start",
                                              "relrank_goal", "donors_ge25_start", "donor_cv_start",
                                              "log2fc_start_vs_alt", "cp10k_start", "cp10k_goal"]}})
        assert n == N_PER_REF, f"{r}: {n} 遺伝子しか選べない"
    cand = pd.DataFrame(rows)
    cand.to_csv(OUT / "candidates.csv", index=False, float_format="%.4g")
    sel = cand[cand.decision == "include"].copy()
    refrows = d.set_index("symbol").loc[REFS].reset_index()
    sel.to_csv(OUT / "selected.csv", index=False, float_format="%.4g")
    summ = {"rule": {"E2_min_log2fc": E2_MIN_LOG2FC, "E2_max_q": E2_MAX_Q, "E3_min_cells": E3_MIN_CELLS,
                     "E3_min_donors_ge25": E3_MIN_DONORS, "tolerances": {k: float(v) for k, v in TOL.items()},
                     "n_per_ref": N_PER_REF},
            "n_pool_E1_E3": int(len(pool)), "n_pool_after_auto_exclusion": int((~pool.E4_auto_excluded).sum()),
            "selected": sel[["ref", "symbol", "ensembl_id", "dist"]].to_dict("records"),
            "refs": refrows[["symbol", "ensembl_id", "log2fc_start_vs_goal", "tok_cells_start", "det_start",
                             "relrank_start", "donors_ge25_start"]].to_dict("records"),
            "source": str(GS.relative_to(ROOT)), "uses_isp_results": False}
    (OUT / "selection.json").write_text(json.dumps(summ, ensure_ascii=False, indent=1, default=float) + "\n")
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(cand[["ref", "rank_by_dist", "symbol", "decision", "uncertainty", "dist", "log2fc_start_vs_goal",
                    "tok_cells_start", "det_start", "det_goal", "relrank_start", "donors_ge25_start"]].round(3).to_string(index=False))
        print(refrows[["symbol", "log2fc_start_vs_goal", "tok_cells_start", "det_start", "relrank_start"]].round(3))
    print(json.dumps(summ["selected"], ensure_ascii=False, default=float))


if __name__ == "__main__":
    main()

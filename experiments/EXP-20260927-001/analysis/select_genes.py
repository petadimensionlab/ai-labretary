"""EXP-20260927-001 事前登録の一部: (A) 発現をそろえた比較遺伝子群、(B) null 遺伝子プールの選定と、
(C) 主検定の実行可能性の確認。ISP の結果は一切読まない（発現統計とトークン化済み入力だけ）。

入力（読み取りのみ）
  - outputs/literature/2026-09-26-negative-controls/genome_stats.csv（EXP-004 と同じ。donor pseudo-bulk の log2FC 等）
  - work/pos_by_gene.parquet（precheck_pos.py。start 細胞での元の位置 pos の中央値・長さの中央値・提供者数）
(A) 比較遺伝子群（陽性対照 1 つにつき 5 遺伝子, plan.md §5.2）
  候補: Geneformer 辞書に収載、start で 300 細胞以上・10 提供者以上で入力に入る
  距離 d_c = max(|Δlog2FC|/0.5, |ln(細胞数比)|/ln1.5, |Δ相対順位|/0.10, |Δlog2(1+pos中央値)|/0.25)
  除外 E4（自動）: EXP-004 の既知リスト + 本実験の追加（PT の同一性・脂肪酸酸化・ミトコンドリア生合成の調節因子、
  NRF2 経路）、盲検、陽性対照（削除型・過剰発現型）、固定遺伝子。E5（文献知識）: MANUAL_REVIEW に理由と不確実性
  陽性対照ごとに d_c の小さい順に E4・E5 を通る 5 遺伝子（重複なし, 陽性対照の並び順に確定）
(B) null プール（固定 18 遺伝子それぞれに 3 遺伝子 = 54）
  候補: (A) と同じ収載・頻度条件、E4 で除外されない、|log2FC| < 1.0、比較遺伝子と重ならない
  距離 d_n = max(|ln(細胞数比)|/ln1.5, |Δlog2(1+pos中央値)|/0.25, |Δlog2(長さ中央値)|/0.25)（発現頻度・元の順位・入力の長さ）
  固定 18 遺伝子の並び順に、d_n の小さい順で 3 遺伝子ずつ（重複なし）
(C) seed 11/12/13 で run_isp.py と同じ規則で細胞を選び、(i) 全遺伝子×seed の細胞数、(ii) 主検定で使える null 遺伝子数
   （対象の細胞の 80% 以上に d ≤ 0.5 のマッチがある null 遺伝子）、(iii) EXP-005 の選抜細胞との重なり
出力: work/comparators.csv, work/comparator_candidates.csv, work/null_pool.json, work/precheck_coverage.json
実行: .venv/bin/python experiments/EXP-20260927-001/analysis/select_genes.py [--list]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts/geneformer"))
import run_isp as ri  # noqa: E402

E = ROOT / "experiments/EXP-20260927-001"
CFG = yaml.safe_load(open(E / "config.yaml"))
TOK = ROOT / CFG["data"]["tokenized_dataset"]
GS = ROOT / "outputs/literature/2026-09-26-negative-controls/genome_stats.csv"
N_COMP, N_NULL = 5, 3
TOL_C = {"log2fc": 0.5, "lncells": np.log(1.5), "relrank": 0.10, "log2pos": 0.25}
TOL_N = {"lncells": np.log(1.5), "log2pos": 0.25, "log2len": 0.25}
NULL_MAX_ABS_LOG2FC = 1.0
MIN_CELLS, MIN_DONORS = 300, 10
BLIND = {"TP53", "CDKN1A", "LTBR", "SPP1", "LCN2"}
POS_OE = ["HNF4A", "PPARGC1A", "PPARA", "ESRRA", "CPT1A", "NFE2L2"]
POS_DEL = ["HIF1A", "NOTCH1", "EGFR", "YAP1", "KEAP1"]
FIXED = POS_OE + ["ACTB", "GAPDH", "PPIA", "CERS4", "TATDN1", "IFT25", "VAPA", "PHPT1", "ANAPC16", "EEF1D",
                  "VCAM1", "HAVCR1"]

# EXP-004 の除外（既知リスト + 追加カテゴリ）をそのまま読む（EXP-004 のコードは変更しない）
_spec = importlib.util.spec_from_file_location("dm", ROOT / "experiments/EXP-20260926-004/analysis/de_matched_select.py")
_dm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_dm)
# 本実験の追加カテゴリ（候補を見る前に決めた。過剰発現型の陽性対照と同じ経路の調節因子 = 駆動因子になりうるもの）
EXTRA_EXCLUDE = {
    # PT の同一性を決める転写因子・核内受容体
    "HNF1A", "HNF1B", "HNF4G", "ONECUT1", "FOXA1", "FOXA2", "FOXA3", "GATA3", "TFCP2L1", "EMX2", "LHX1", "OSR1",
    "SALL1", "SIX2", "WT1", "PAX2", "PAX8", "KLF15", "KLF6", "KLF4", "KLF5", "TFEB",
    # 脂肪酸酸化・ミトコンドリア生合成・エネルギー代謝の調節因子
    "PPARGC1B", "TFAM", "NRF1", "GABPA", "SIRT1", "SIRT3", "PRKAA1", "PRKAA2", "MLYCD", "ACACA", "ACACB",
    "SLC25A20", "SLC27A2", "FABP1", "CD36", "LPL", "PDK4", "UCP2", "MFN2", "OPA1", "PINK1", "PRKN", "SREBF1",
    "SREBF2", "MLXIPL", "INSR", "FOXO1", "FOXO3",
    # NRF2 経路
    "KEAP1", "MAFF", "MAFG", "MAFK", "BACH1", "HMOX1", "SQSTM1", "G6PD", "TXN", "SLC7A11", "ME1",
}
EXTRA_PREFIX = ("NR1", "NR2", "NR4", "NR5", "NR0", "RXR", "RAR", "PPAR", "ESRR", "HNF", "ACAD", "HADH", "CPT",
                "ACOX", "ACSL", "NDUF", "SDH", "UQCR", "COX", "ATP5", "MRPL", "MRPS", "RPL", "RPS", "EIF",
                "SLC25A", "SIRT", "FOXO", "GST", "GPX", "PRDX")

# E5 文献知識による確認（距離順に現れたものを記録する。ISP 結果は見ていない）
# decision: include / exclude。uncertainty: 低 / 中 / 高（「PT の傷害・修復・老化の駆動因子の証拠がない」判断の確からしさ）
MANUAL_REVIEW: dict[str, tuple[str, str, str]] = {
    # HNF4A の帯
    "UGT3A1": ("include", "UDP 糖転移酵素 3A1（PT に発現する代謝酵素）。傷害・修復の駆動因子とする報告は把握していない", "中"),
    "METRN": ("include", "meteorin（神経栄養性の分泌タンパク質）。腎保護の報告があるのは別遺伝子の METRNL で、METRN の報告は把握していない", "高"),
    "CISD3": ("include", "ミトコンドリアの鉄硫黄タンパク質（構成要素で調節因子ではない）。腎傷害の報告は把握していない", "中"),
    "SNX5": ("include", "sorting nexin 5（エンドソーム輸送。PT のドーパミン受容体の再利用の報告はあるが傷害の駆動因子ではない）", "中"),
    "DDX10": ("include", "RNA ヘリカーゼ（リボソーム生合成）。腎での機能報告は把握していない", "低"),
    # PPARGC1A の帯
    "AKAP11": ("include", "A キナーゼ足場タンパク質 11。腎傷害の報告は把握していない", "中"),
    "SBNO1": ("exclude", "strawberry notch 相同遺伝子 1（Notch 経路の調節因子）。Notch は PT 傷害・線維化の経路 → 経路単位で保守的に除外", "中"),
    "TERF1": ("exclude", "テロメア結合タンパク質（shelterin）。テロメアは細胞老化の駆動経路 → 本研究の主題（老化）に関わるため除外", "中"),
    "ZNF24": ("include", "ジンクフィンガー転写因子 24。腎での機能報告は把握していない", "低"),
    "ARFGAP3": ("include", "ARF GTPase 活性化タンパク質 3（ゴルジ輸送）。腎傷害の報告は把握していない", "低"),
    "TMUB1": ("include", "膜貫通ユビキチン様タンパク質 1。腎での機能報告は把握していない", "低"),
    "HOXA9": ("include", "HOX 転写因子。腎の発生で主要なのは HOX10/11 群で、HOXA9 の PT 傷害への関与は把握していない", "中"),
    # PPARA の帯
    "FKBP4": ("include", "FKBP52（ステロイド受容体のシャペロン補因子）。PT 傷害の駆動因子とする報告は把握していない", "中"),
    "KLF13": ("exclude", "Krüppel 様転写因子。KLF4/5/6/15 に腎傷害の調節因子の報告が多く、ファミリー単位で保守的に除外", "中"),
    "WBP2": ("exclude", "YAP/TAZ の転写共役因子 → YAP 経路（EXP-004 の除外カテゴリ）", "中"),
    "MGLL": ("exclude", "モノアシルグリセロールリパーゼ。阻害で腎虚血再灌流傷害が軽減するとの報告がある（記憶に基づく, 未検証）", "高"),
    "MSH3": ("include", "DNA ミスマッチ修復。腎傷害の駆動因子とする報告は把握していない", "低"),
    "PDXDC1": ("include", "ピリドキサール依存性脱炭酸酵素ドメイン 1。腎での機能報告は把握していない", "低"),
    "PHGDH": ("exclude", "セリン合成の律速酵素。腎臓病でセリン代謝が保護的に働くとの報告がある（記憶に基づく, 未検証）", "高"),
    "HIPK3": ("exclude", "HIPK ファミリー（HIPK2 は腎線維化の駆動因子）→ ファミリー単位で保守的に除外", "中"),
    "ABHD17A": ("include", "脱パルミトイル化酵素。腎での機能報告は把握していない", "低"),
    "VPS51": ("include", "GARP/EARP 複合体（エンドソーム輸送）。腎傷害の報告は把握していない", "低"),
    # ESRRA の帯
    "PGPEP1": ("include", "ピログルタミルペプチダーゼ。腎での機能報告は把握していない", "低"),
    "DHTKD1": ("include", "ミトコンドリアの 2-オキソアジピン酸脱水素酵素（代謝酵素で調節因子ではない）", "中"),
    "ABHD10": ("include", "脱パルミトイル化・グルクロニド加水分解酵素。腎傷害の報告は把握していない", "低"),
    "PNPLA4": ("include", "パタチン様リパーゼ 4。腎傷害の報告は把握していない", "中"),
    "HDAC6": ("exclude", "HDAC6 の阻害が腎傷害・線維化を軽減する報告が複数ある", "低"),
    "AFMID": ("exclude", "キヌレニン経路（NAD のデノボ合成）の酵素。NAD 合成は AKI の保護経路として報告がある → 経路単位で除外", "中"),
    "PALM": ("include", "paralemmin（膜の形態）。腎での機能報告は把握していない", "低"),
    # CPT1A の帯
    "SH3RF1": ("include", "POSH（JNK の足場 E3 リガーゼ）。PT 傷害の駆動因子とする報告は把握していない", "中"),
    "KRCC1": ("include", "機能未詳のタンパク質。腎での報告は把握していない", "低"),
    "CHCHD7": ("include", "ミトコンドリア膜間腔タンパク質（構成要素）。腎傷害の報告は把握していない", "低"),
    "SMNDC1": ("include", "スプライシング因子 SPF30。腎での報告は把握していない", "低"),
    "TRAF4": ("exclude", "TNF 受容体結合因子 4（TGF-β・炎症シグナルの調節）→ 炎症・TGF-β カテゴリとして除外", "中"),
    "ETF1": ("include", "翻訳終結因子 eRF1。腎傷害の駆動因子とする報告は把握していない", "低"),
    # NFE2L2 の帯
    "SRI": ("include", "sorcin（Ca 結合タンパク質）。腎傷害の報告は把握していない", "低"),
    "PDIA6": ("exclude", "ER のジスルフィド異性化酵素（IRE1 を介した小胞体ストレス応答の調節）→ ストレス応答カテゴリとして除外", "中"),
    "ANXA11": ("include", "アネキシン A11。腎傷害の駆動因子とする報告は把握していない", "低"),
    "PHF3": ("include", "PHD フィンガー 3（転写伸長）。腎での報告は把握していない", "低"),
    "PRP4K": ("include", "スプライシングのキナーゼ PRPF4B。腎での報告は把握していない", "低"),
    "SELENOF": ("include", "セレノプロテイン F（ER の糖タンパク質の品質管理）。腎傷害の報告は把握していない（EXP-004 の SELENOM と同じ扱い）", "中"),
}


# null 側の追加除外（E5 で除外した遺伝子に加え、null 候補に現れた明らかな傷害メディエーター等。ISP 結果は見ていない）
NULL_EXCLUDE = {
    "HMGB1": "傷害関連分子パターン（DAMP）。腎傷害・炎症のメディエーター",
    "MGST1": "ミクロソームのグルタチオン S 転移酵素（NRF2 標的・酸化ストレス応答）",
    "MAPKAPK2": "MK2（p38 経路の炎症・ストレス応答キナーゼ）",
}


def load_table() -> pd.DataFrame:
    d = pd.read_csv(GS)
    p = pd.read_parquet(E / "work/pos_by_gene.parquet")
    d = d.merge(p[["ensembl", "n_start", "n_donor", "pos_med", "len_med", "relpos_med"]],
                left_on="ensembl_id", right_on="ensembl", how="inner")
    neg = _dm._neg
    auto = d.symbol.isin(neg.EXCLUDE_SYMBOLS | _dm.EXTRA_EXCLUDE | EXTRA_EXCLUDE | BLIND | set(POS_OE) | set(POS_DEL)
                         | set(FIXED)) | d.symbol.str.startswith(neg.EXCLUDE_PREFIX + _dm.EXTRA_PREFIX + EXTRA_PREFIX) \
        | d.known_excluded.astype(bool)
    d["E4_auto_excluded"] = auto
    d["ok"] = d.in_gf_vocab.astype(bool) & (d.n_start >= MIN_CELLS) & (d.n_donor >= MIN_DONORS)
    d["l2pos"] = np.log2(1 + d.pos_med)
    d["l2len"] = np.log2(d.len_med)
    return d


def comparators(d: pd.DataFrame, listing: bool):
    ref = d.set_index("symbol")
    pool = d[d.ok].copy()
    rows, chosen = [], []
    for r in POS_OE:
        x = pool.copy()
        x["d_fc"] = (x.log2fc_start_vs_goal - ref.loc[r, "log2fc_start_vs_goal"]).abs() / TOL_C["log2fc"]
        x["d_cells"] = np.log(x.n_start / ref.loc[r, "n_start"]).abs() / TOL_C["lncells"]
        x["d_rank"] = (x.relrank_start - ref.loc[r, "relrank_start"]).abs() / TOL_C["relrank"]
        x["d_pos"] = (x.l2pos - ref.loc[r, "l2pos"]).abs() / TOL_C["log2pos"]
        x["dist"] = x[["d_fc", "d_cells", "d_rank", "d_pos"]].max(1)
        x = x.sort_values("dist")
        n, k = 0, 0
        for _, g in x.iterrows():
            if n >= N_COMP or (listing and k >= 25):
                break
            s = g.symbol
            if s in chosen:
                continue
            k += 1
            if g.E4_auto_excluded:
                dec, why, unc = "exclude_auto", "既知リスト / 接頭辞 / 陽性・盲検・固定（E4）", "—"
            elif s in MANUAL_REVIEW:
                dec, why, unc = MANUAL_REVIEW[s]
            elif listing:
                dec, why, unc = "UNREVIEWED", "", ""
            else:
                raise SystemExit(f"{s}（{r} の候補）は文献確認が未記入")
            if dec == "include":
                chosen.append(s); n += 1
            rows.append({"ref": r, "symbol": s, "ensembl_id": g.ensembl_id, "decision": dec, "reason": why,
                         "uncertainty": unc, **{c: g[c] for c in ["dist", "d_fc", "d_cells", "d_rank", "d_pos",
                                                                  "log2fc_start_vs_goal", "welch_q_goal", "n_start",
                                                                  "det_start", "det_goal", "relrank_start", "pos_med",
                                                                  "len_med", "n_donor"]}})
        if not listing:
            assert n == N_COMP, f"{r}: {n} 遺伝子しか選べない"
    return pd.DataFrame(rows), chosen


def null_pool(d: pd.DataFrame, comp: list[str]) -> list[str]:
    ref = d.set_index("symbol")
    pool = d[d.ok & ~d.E4_auto_excluded & (d.log2fc_start_vs_goal.abs() < NULL_MAX_ABS_LOG2FC)
             & ~d.symbol.isin(comp) & ~d.symbol.isin(set(NULL_EXCLUDE))
             & ~d.symbol.isin({k for k, v in MANUAL_REVIEW.items() if v[0] == "exclude"})].copy()
    chosen = []
    for r in FIXED:
        x = pool[~pool.symbol.isin(chosen)].copy()
        x["dist"] = pd.concat([np.log(x.n_start / ref.loc[r, "n_start"]).abs() / TOL_N["lncells"],
                               (x.l2pos - ref.loc[r, "l2pos"]).abs() / TOL_N["log2pos"],
                               (x.l2len - ref.loc[r, "l2len"]).abs() / TOL_N["log2len"]], axis=1).max(1)
        take = x.sort_values(["dist", "symbol"]).head(N_NULL)
        chosen += take.symbol.tolist()
    return chosen


def coverage(d: pd.DataFrame, fixed: list[str], comp: list[str], nulls: list[str]) -> dict:
    from datasets import load_from_disk
    sym2id = dict(zip(d.symbol, d.ensembl_id))
    m, rows, donors = ri.start_matrix(TOK, CFG)
    gdict = ri.gene_token_dict()
    sel = CFG["cell_selection"]
    ids_all = load_from_disk(str(TOK)).select(rows.tolist())["input_ids"]

    def pick(sym, seed, nmax, per_d):
        idx = m[:, gdict[sym2id[sym]]].nonzero()[0].tolist()
        return ri.balanced_subset(idx, donors, nmax, per_d, random.Random(f"{seed}-{sym2id[sym]}"))

    def prof(sym, keep):
        tid = gdict[sym2id[sym]]
        return (np.log2(1 + np.array([list(ids_all[i]).index(tid) for i in keep])),
                np.log2(np.array([len(ids_all[i]) for i in keep])))

    seeds = CFG["robustness"]["seeds"]
    tol = CFG["analysis"]["pos_matched_permutation"]["within_gene_tolerance_log2"]
    frac = CFG["analysis"]["pos_matched_permutation"]["min_cover_frac"]
    genes = fixed + comp + nulls
    ncell, P = {}, {}
    for s in seeds:
        for g in genes:
            k = pick(g, s, sel["max_cells_per_gene"], sel["max_cells_per_donor"])
            ncell[f"{g}@{s}"] = len(k)
            P[(s, g)] = prof(g, k)
    elig = {}
    for s in seeds:
        for g in genes:
            lp, ll = P[(s, g)]
            n = 0
            for h in nulls:
                if h == g:
                    continue
                hp, hl = P[(s, h)]
                dd = np.maximum(np.abs(lp[:, None] - hp[None, :]), np.abs(ll[:, None] - hl[None, :]))
                n += (dd <= tol).any(axis=1).mean() >= frac
            elig[f"{g}@{s}"] = int(n)
    overlap = {}
    for g in fixed + comp:
        old = set().union(*[set(pick(g, s, 200, 25)) for s in (3, 4, 5)])
        overlap[g] = [round(len(set(pick(g, s, sel["max_cells_per_gene"], sel["max_cells_per_donor"])) & old)
                            / sel["max_cells_per_gene"], 3) for s in seeds]
    grp = {"fixed": fixed, "comparator": comp, "null": nulls}
    return {"seeds": seeds, "tolerance_log2": tol, "min_cover_frac": frac,
            "cells_min": int(min(ncell.values())),
            "cells_below_min": {k: v for k, v in ncell.items() if v < sel["min_cells_per_gene"]},
            "eligible_null_min_by_group": {k: int(min(elig[f"{g}@{s}"] for g in v for s in seeds)) for k, v in grp.items()},
            "eligible_null_median_by_group": {k: float(np.median([elig[f"{g}@{s}"] for g in v for s in seeds]))
                                              for k, v in grp.items()},
            "eligible_null_fixed": {g: [elig[f"{g}@{s}"] for s in seeds] for g in fixed},
            "eligible_null_below_20": {k: v for k, v in elig.items() if v < 20},
            "overlap_with_EXP005_cells_frac": overlap, "uses_isp_results": False}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true", help="文献確認のため候補を距離順に表示するだけ")
    a = ap.parse_args()
    d = load_table()
    cand, comp = comparators(d, a.list)
    with pd.option_context("display.width", 250, "display.max_rows", 500):
        print(cand[["ref", "symbol", "decision", "dist", "log2fc_start_vs_goal", "n_start", "relrank_start",
                    "pos_med", "n_donor"]].round(3).to_string(index=False))
    if a.list:
        return
    nulls = null_pool(d, comp)
    assert len(set(nulls)) == len(FIXED) * N_NULL and not set(nulls) & (set(comp) | set(FIXED) | BLIND)
    cand.to_csv(E / "work/comparator_candidates.csv", index=False, float_format="%.4g")
    cand[cand.decision == "include"].to_csv(E / "work/comparators.csv", index=False, float_format="%.4g")
    t = d.set_index("symbol")
    (E / "work/null_pool.json").write_text(json.dumps({
        "rule": "固定 18 遺伝子ごとに d_n（発現頻度・元の位置・入力の長さ）の近い 3 遺伝子, |log2FC|<1, E4・E5 除外＋NULL_EXCLUDE",
        "null_exclude": NULL_EXCLUDE,
        "tolerances": {k: float(v) for k, v in TOL_N.items()}, "n": len(nulls),
        "genes": [{"symbol": s, "ensembl": t.loc[s, "ensembl_id"], "matched_to": FIXED[i // N_NULL],
                   "n_start": int(t.loc[s, "n_start"]), "pos_med": float(t.loc[s, "pos_med"]),
                   "len_med": float(t.loc[s, "len_med"]), "log2fc": round(float(t.loc[s, "log2fc_start_vs_goal"]), 3)}
                  for i, s in enumerate(nulls)], "uses_isp_results": False}, ensure_ascii=False, indent=1) + "\n")
    cov = coverage(d, FIXED, comp, nulls)
    (E / "work/precheck_coverage.json").write_text(json.dumps(cov, ensure_ascii=False, indent=1) + "\n")
    print("comparators", comp)
    print("nulls", nulls)
    print({k: cov[k] for k in ["cells_min", "cells_below_min", "eligible_null_min_by_group",
                               "eligible_null_median_by_group", "eligible_null_below_20"]})
    print("fixed", cov["eligible_null_fixed"])


if __name__ == "__main__":
    main()

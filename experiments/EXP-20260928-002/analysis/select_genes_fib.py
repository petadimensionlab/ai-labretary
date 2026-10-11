"""EXP-20260927-003 事前登録の一部（摂動の実行前に確定する）。

ISP（摂動）の結果は一切読まない。使うのは work/gene_stats_fib.csv（トークン化済み入力と生カウントから
作った遺伝子統計表）と、トークン化済み入力の長さだけである。

(A) 発現をそろえた比較遺伝子群（主対象 1 つにつき 5 遺伝子。plan.md §5.2）
(B) null 遺伝子プール（固定遺伝子 1 つにつき 3 遺伝子。plan.md §5.1）
(C) 深さ（トークン長）をそろえた状態重心用の細胞（seed ごと。plan.md §4.3）
(D) 実行可能性の確認（遺伝子ごとの細胞数、主検定で「使える null 遺伝子」の数）

出力: work/{comparators.csv, comparator_candidates.csv, null_pool.json,
       state_rows-seed{S}.json, depth_match-seed{S}.json, precheck_coverage.json}
実行: .venv/bin/python experiments/EXP-20260927-003/analysis/select_genes_fib.py [--list]
"""
from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260927-003"
TOK = E / "work/tokenized.dataset"
STATS = E / "work/gene_stats_fib.csv"

SEEDS = [21, 22, 23]                 # 既出（0〜5, 11〜13）と重ならない
STATE_MAX, STATE_PER_DONOR = 800, 60   # EXP-20260927-002 と同じ
DEPTH_TOL = 0.10                     # |Δlog(トークン長)| ≤ 0.10 の最近傍 1:1（EXP-20260927-002 と同じ）
MAX_CELLS_PER_GENE, MAX_CELLS_PER_DONOR, MIN_CELLS_PER_GENE = 100, 15, 80

N_COMP, N_NULL = 5, 3
TOL_C = {"log2fc": 0.5, "lncells": np.log(1.5), "relpos": 0.10, "log2pos": 0.25}
TOL_N = {"lncells": np.log(1.5), "log2pos": 0.25, "log2len": 0.25}
NULL_MAX_ABS_LOG2FC = 1.0
MIN_CELLS, MIN_DONORS = 300, 10            # 高頻度側の候補条件
LOW_CELLS = (30, 150)                      # 低頻度側（CDKN1A 用）の候補条件
LOW_MIN_DONORS = 8

# 主検定の設定（config.yaml の analysis.pos_matched_permutation と同じ値。plan.md §6.1）
PERM_TOL, PERM_COVER, PERM_MIN_NULL = 0.5, 0.8, 20

# ---- 固定遺伝子（役割ごと。摂動対象は CCL19/CCL21/CXCL13 を含まない = 循環性の回避, plan.md §4.1）
TARGET_MAIN = {"LTBR": "ENSG00000111321", "RELB": "ENSG00000104856", "NFKB2": "ENSG00000077150"}
TARGET_REF = {"CDKN1A": "ENSG00000124762"}          # 参考（細胞数不足。主判定に数えない）
HK = {"ACTB": "ENSG00000075624", "GAPDH": "ENSG00000111640", "PPIA": "ENSG00000196262"}
PANEL = {"CERS4": "ENSG00000090661", "TATDN1": "ENSG00000147687", "IFT25": "ENSG00000081870",
         "VAPA": "ENSG00000101558", "PHPT1": "ENSG00000054148", "ANAPC16": "ENSG00000166295",
         "EEF1D": "ENSG00000104529"}
FIBRO_DIR = {"ACTA2": "ENSG00000107796", "TAGLN": "ENSG00000149591"}   # 線維化方向の対照（参考）
FIXED_HIGH = list(TARGET_MAIN) + list(HK) + list(PANEL) + list(FIBRO_DIR)   # 15 遺伝子（null の参照）

# ---- 自動除外（候補を見る前に決めた。本系（腎線維芽細胞・三次リンパ組織・細胞老化）の駆動因子になりうるもの）
EXCLUDE_SYMBOLS = {
    # 非古典的・古典的 NF-κB 経路とその上流受容体
    "NFKB1", "NFKB2", "RELA", "RELB", "REL", "CHUK", "IKBKB", "IKBKG", "MAP3K14", "NFKBIA", "NFKBIB",
    "NFKBIE", "NFKBID", "NFKBIZ", "TRAF1", "TRAF2", "TRAF3", "TRAF5", "TRAF6", "BIRC2", "BIRC3",
    "TNFRSF11A", "TNFRSF13B", "TNFRSF13C", "TNFRSF17", "CD40", "CD40LG", "LTBR", "LTA", "LTB",
    "TNFSF14", "TNF", "TNFRSF1A", "TNFRSF1B", "TNFAIP3", "TNIP1", "CYLD", "OTULIN", "RIPK1", "RIPK2",
    # 三次リンパ組織・オーガナイザーのケモカイン／接着分子
    "CCL19", "CCL21", "CXCL13", "CCL2", "CCL5", "CXCL9", "CXCL10", "CXCL11", "CXCL12", "CXCR4",
    "CXCR5", "CCR7", "IL7", "IL7R", "TNFSF13B", "PDPN", "NGFR", "VCAM1", "ICAM1", "MADCAM1",
    "PECAM1", "PROX1", "LYVE1", "PTGDS", "CR2", "FDCSP", "CLU",
    # 細胞老化・DNA 損傷応答（本研究の主題）
    "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B", "TP53", "MDM2", "RB1", "E2F1", "GLB1", "LMNB1",
    "SERPINE1", "IGFBP3", "IGFBP5", "IGFBP7", "ATM", "ATR", "CHEK1", "CHEK2", "TERT", "TERF1",
    "TERF2", "TP53BP1", "TRP53", "SIRT1", "SIRT6", "FOXO3", "HMGB1", "GDF15",
    # 線維化・筋線維芽細胞化
    "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3", "SMAD2", "SMAD3", "SMAD4", "SMAD7",
    "ACTA2", "TAGLN", "CTGF", "CCN2", "POSTN", "FN1", "YAP1", "WWTR1", "PDGFRA", "PDGFRB", "THY1",
    "LOX", "LOXL1", "LOXL2", "TIMP1", "MMP2", "MMP9", "SPP1", "LCN2",
    # インターフェロン・自然免疫シグナル
    "STAT1", "STAT2", "STAT3", "STAT6", "IRF1", "IRF3", "IRF4", "IRF7", "IFNG", "JAK1", "JAK2",
    "MYD88", "TICAM1", "STING1", "TMEM173", "CGAS", "MB21D1", "NLRP3", "CASP1", "PYCARD",
    # 低酸素・Notch・Wnt・Hippo（腎傷害で繰り返し出る経路）
    "HIF1A", "EPAS1", "VHL", "NOTCH1", "NOTCH2", "NOTCH3", "JAG1", "HES1", "HEY1", "RBPJ",
    "CTNNB1", "WNT4", "WNT5A", "LEF1", "TCF7", "TCF7L2", "AXIN2",
}
EXCLUDE_PREFIX = ("NFKB", "TRAF", "CCL", "CXCL", "CCR", "CXCR", "TNFSF", "TNFRSF", "IL1", "IL2", "IL6",
                  "IL10", "IL11", "IL17", "IL33", "COL", "SMAD", "IFI", "ISG", "OAS", "MX1", "HLA-",
                  "IGH", "IGK", "IGL", "MMP", "TIMP", "CDKN", "WNT", "FZD", "NOTCH", "JAK", "STAT",
                  "SERPIN", "LOXL", "TLR", "MT-", "RPL", "RPS", "MRPL", "MRPS")
BLIND = {"TP53", "CDKN1A", "LTBR", "SPP1", "LCN2"}   # 盲検（解除の指示があるまで摂動・解析しない）
LABEL_GENES = {"CCL19", "CCL21", "CXCL13"}           # goal ラベル定義遺伝子（循環性のため対象から外す）

# ---- 文献知識による確認（E5）。距離順に現れた候補すべてに判断・理由・不確実性を記録する。
# decision: include / exclude。uncertainty: 低 / 中 / 高（「本系の駆動因子である報告を把握していない」判断の確からしさ）
MANUAL_REVIEW: dict[str, tuple[str, str, str]] = {
    # --- LTBR の帯
    "TBRG1": ("exclude", "NIAM（ARF・MDM2 と結合する核小体タンパク質）。p53 経路＝細胞老化の主題に関わるため保守的に除外", "中"),
    "HOXA6": ("include", "HOX 転写因子。腎の発生で主要なのは HOX9〜11 群で、線維芽細胞の三次リンパ組織化への関与は把握していない", "中"),
    "TRIB2": ("exclude", "tribbles 偽キナーゼ。TRIB1/2/3 は NF-κB・MAPK シグナルの調節因子の報告がある → ファミリー単位で除外", "中"),
    "MTIF3": ("include", "ミトコンドリア翻訳開始因子 3（翻訳装置の構成因子）。腎・免疫での駆動因子の報告は把握していない", "低"),
    "ATP6V1B2": ("include", "液胞型 ATPase の V1 B2 サブユニット（構成因子）。三次リンパ組織・老化の駆動因子の報告は把握していない", "低"),
    "OAT": ("include", "オルニチンアミノトランスフェラーゼ（代謝酵素で調節因子ではない）", "低"),
    "SPG21": ("include", "maspardin（エンドソーム輸送）。腎線維芽細胞での駆動因子の報告は把握していない", "低"),
    # --- RELB の帯
    "AGTR1": ("exclude", "アンジオテンシン II 1 型受容体。腎の炎症・線維化の主要な駆動因子", "低"),
    "CD226": ("exclude", "DNAM-1（リンパ球の共刺激・接着受容体）。リンパ組織形成に関わりうるため保守的に除外", "中"),
    "SLC26A7": ("include", "陰イオン輸送体（腎遠位尿細管）。線維芽細胞の状態変化の駆動因子の報告は把握していない", "低"),
    "TPRG1": ("exclude", "tumor protein p63 regulated 1。p63 は p53 ファミリー → 老化・p53 経路として保守的に除外", "中"),
    "SHC3": ("include", "Shc アダプター 3（受容体型キナーゼのアダプター、主に神経系）。腎線維芽細胞での報告は把握していない", "中"),
    "MEGF9": ("include", "multiple EGF-like domains 9（膜タンパク質）。駆動因子の報告は把握していない", "低"),
    "OSBPL3": ("include", "オキシステロール結合タンパク質様 3（脂質輸送）。駆動因子の報告は把握していない", "低"),
    "VPS35": ("include", "レトロマー複合体の構成因子（小胞輸送）。駆動因子の報告は把握していない", "低"),
    # --- NFKB2 の帯
    "ECH1": ("include", "エノイル CoA ヒドラターゼ 1（脂肪酸酸化の酵素で調節因子ではない）", "低"),
    "ARL2": ("include", "ADP リボシル化因子様 2（チューブリン折りたたみ・ミトコンドリア）。駆動因子の報告は把握していない", "低"),
    "TGM2": ("exclude", "トランスグルタミナーゼ 2。細胞外基質の架橋と TGF-β 活性化を介して腎線維化を駆動する報告がある", "低"),
    "CCT8": ("include", "シャペロニン CCT の構成サブユニット。駆動因子の報告は把握していない", "低"),
    "KIAA0040": ("include", "機能未詳のタンパク質。腎・免疫での報告は把握していない", "低"),
    "LONRF2": ("include", "LON ペプチダーゼ N 末端・RING フィンガー 2。駆動因子の報告は把握していない", "低"),
    "DEPP1": ("exclude", "FOXO3 の標的でオートファジー・酸化ストレスを調節する → ストレス応答・老化に隣接するため保守的に除外", "中"),
    "FRZB": ("exclude", "分泌型 frizzled 関連タンパク質（Wnt 拮抗因子）→ Wnt 経路として除外", "中"),
    # --- CDKN1A の帯（参考対象。低頻度のため候補は免疫系・神経系の低発現遺伝子が中心になる）
    "PIK3R5": ("exclude", "PI3K 制御サブユニット p101。PI3K-AKT は増殖・老化の経路 → 経路単位で保守的に除外", "中"),
    "GAD1": ("include", "グルタミン酸脱炭酸酵素 1（神経系の代謝酵素）。腎線維芽細胞での報告は把握していない", "低"),
    "ASTN1": ("include", "アストロタクチン 1（神経の接着分子）。腎での報告は把握していない", "低"),
    "CYTIP": ("include", "cytohesin 結合タンパク質（免疫細胞）。線維芽細胞では周辺 RNA 由来の可能性が高く、駆動因子の報告は把握していない", "中"),
    "S100A1": ("include", "S100 カルシウム結合タンパク質 A1（心筋・神経）。炎症性 DAMP の S100A8/A9 とは別分子", "中"),
    "GPR34": ("include", "プリン受容体（主にミクログリア）。腎線維芽細胞での報告は把握していない", "中"),
}


# ---- null プールだけに効く追加除外（最初の抽出で現れた候補を、ISP 結果を見ずに文献知識で点検した結果）
NULL_EXCLUDE: dict[str, str] = {
    "RBCK1": "HOIL-1（LUBAC の構成因子）。直線状ユビキチン鎖を介した NF-κB 活性化の中心因子",
    "CKLF": "chemokine-like factor（ケモカイン様因子）。走化性因子としてリンパ球集積に関わりうる",
    "CCND1": "サイクリン D1。細胞周期停止＝細胞老化の主題に直結する",
    "KLF6": "KLF ファミリー転写因子。腎傷害・線維化の調節因子の報告が多い（EXP-20260927-001 と同じ扱い）",
    "B2M": "β2 ミクログロブリン。MHC クラス I の構成因子であり、老化促進因子としての報告もある",
    "SULF2": "ヘパラン硫酸 6-O 脱硫酸化酵素。Wnt・FGF など複数経路のシグナル強度を変える",
    "TNC": "テネイシン C。線維化・傷害関連の細胞外基質で DAMP としても働く",
    "LUM": "ルミカン。線維化に関わる細胞外基質プロテオグリカン",
    "ZFP36L2": "TTP ファミリーの RNA 結合タンパク質。炎症性転写産物の安定性を介して NF-κB 応答を調節する",
    "ELMO1": "Rac の GEF 複合体。糖尿病性腎症の感受性遺伝子で線維化促進の報告がある",
    "RASSF4": "RASSF ファミリー（RASSF1A は老化・腫瘍抑制）→ ファミリー単位で保守的に除外",
    "MTUS1": "ATIP（アンジオテンシン II 2 型受容体結合タンパク質）→ アンジオテンシン経路",
    "NPM1": "ヌクレオフォスミン。ARF・p53 と結合し老化経路に関わる",
    "FAM162A": "HIF-1 標的（低酸素誘導アポトーシス因子）→ 低酸素経路",
    "LBH": "Wnt 標的。関節リウマチの線維芽細胞様滑膜細胞の活性化因子としての報告がある（線維芽細胞の活性化に隣接）",
    "MSRA": "メチオニンスルホキシド還元酵素 A（酸化ストレス応答）",
    "PLAC8": "placenta specific 8。腎細胞癌・尿細管再生での報告がある",
}
NULL_EXCLUDE_PREFIX = ("PRDX", "KLF", "RASSF")


# ------------------------------------------------------------------ 読み込み
def load_stats() -> pd.DataFrame:
    d = pd.read_csv(STATS)
    d = d[d.symbol.astype(str) != ""].copy()
    d["l2pos"] = np.log2(1 + d.pos_med)
    d["l2len"] = np.log2(d.len_med)
    fixed_syms = set(FIXED_HIGH) | set(TARGET_REF)
    d["auto_excluded"] = (d.symbol.isin(EXCLUDE_SYMBOLS | BLIND | LABEL_GENES | fixed_syms)
                          | d.symbol.str.startswith(EXCLUDE_PREFIX))
    d["ok_high"] = (d.n_start >= MIN_CELLS) & (d.n_donor_start >= MIN_DONORS)
    d["ok_low"] = (d.n_start.between(*LOW_CELLS)) & (d.n_donor_start >= LOW_MIN_DONORS)
    return d


# ------------------------------------------------------------------ (A) 比較遺伝子群
def comparators(d: pd.DataFrame, listing: bool):
    ref = d.set_index("symbol")
    rows, chosen = [], []
    for r, pool_col in [(s, "ok_high") for s in TARGET_MAIN] + [(s, "ok_low") for s in TARGET_REF]:
        x = d[d[pool_col]].copy()
        x["d_fc"] = (x.log2fc_start_vs_goal - ref.loc[r, "log2fc_start_vs_goal"]).abs() / TOL_C["log2fc"]
        x["d_cells"] = np.log(x.n_start / ref.loc[r, "n_start"]).abs() / TOL_C["lncells"]
        x["d_relpos"] = (x.relpos_med - ref.loc[r, "relpos_med"]).abs() / TOL_C["relpos"]
        x["d_pos"] = (x.l2pos - ref.loc[r, "l2pos"]).abs() / TOL_C["log2pos"]
        x["dist"] = x[["d_fc", "d_cells", "d_relpos", "d_pos"]].max(1)
        x = x.sort_values(["dist", "symbol"])
        n, k = 0, 0
        for _, g in x.iterrows():
            if n >= N_COMP or (listing and k >= 30):
                break
            s = g.symbol
            if s in chosen:
                continue
            k += 1
            if g.auto_excluded:
                dec, why, unc = "exclude_auto", "自動除外（既知リスト / 接頭辞 / 盲検 / ラベル定義 / 固定遺伝子）", "—"
            elif s in MANUAL_REVIEW:
                dec, why, unc = MANUAL_REVIEW[s]
            elif listing:
                dec, why, unc = "UNREVIEWED", "", ""
            else:
                raise SystemExit(f"{s}（{r} の候補）は文献確認が未記入")
            if dec == "include":
                chosen.append(s)
                n += 1
            rows.append({"ref": r, "symbol": s, "ensembl": g.ensembl, "decision": dec, "reason": why,
                         "uncertainty": unc,
                         **{c: g[c] for c in ["dist", "d_fc", "d_cells", "d_relpos", "d_pos",
                                              "log2fc_start_vs_goal", "n_start", "n_donor_start",
                                              "relpos_med", "pos_med", "len_med", "frac_start", "frac_goal"]}})
        if not listing and n != N_COMP:
            raise SystemExit(f"{r}: {n} 遺伝子しか選べない")
    return pd.DataFrame(rows), chosen


# ------------------------------------------------------------------ (B) null プール
def null_pool(d: pd.DataFrame, comp: list[str]) -> list[tuple[str, str]]:
    ref = d.set_index("symbol")
    base = d[~d.auto_excluded & (d.log2fc_start_vs_goal.abs() < NULL_MAX_ABS_LOG2FC)
             & ~d.symbol.isin(comp)
             & ~d.symbol.isin(set(NULL_EXCLUDE))
             & ~d.symbol.str.startswith(NULL_EXCLUDE_PREFIX)
             & ~d.symbol.isin({k for k, v in MANUAL_REVIEW.items() if v[0] == "exclude"})]
    chosen: list[tuple[str, str]] = []
    for r, pool_col in [(s, "ok_high") for s in FIXED_HIGH] + [(s, "ok_low") for s in TARGET_REF]:
        taken = {s for s, _ in chosen}
        x = base[base[pool_col] & ~base.symbol.isin(taken)].copy()
        x["dist"] = pd.concat([
            np.log(x.n_start / ref.loc[r, "n_start"]).abs() / TOL_N["lncells"],
            (x.l2pos - ref.loc[r, "l2pos"]).abs() / TOL_N["log2pos"],
            (x.l2len - ref.loc[r, "l2len"]).abs() / TOL_N["log2len"]], axis=1).max(1)
        for s in x.sort_values(["dist", "symbol"]).head(N_NULL).symbol:
            chosen.append((s, r))
    return chosen


# ------------------------------------------------------------------ (C) 深さをそろえた状態重心の細胞
def depth_match(log_len: np.ndarray, pos: np.ndarray, neg: np.ndarray, tol: float) -> np.ndarray:
    """log(トークン長)の最近傍 1:1 マッチ（置換なし）。EXP-20260927-002 の precheck_ci.depth_match と同じ。"""
    pi, ni = np.where(pos)[0], np.where(neg)[0]
    ni = ni[np.argsort(log_len[ni], kind="stable")]
    nv = log_len[ni]
    used = np.zeros(len(ni), bool)
    keep: list[int] = []
    for i in pi:
        v = log_len[i]
        cand = np.where((~used) & (np.abs(nv - v) <= tol))[0]
        if len(cand):
            j = cand[np.argmin(np.abs(nv[cand] - v))]
            used[j] = True
            keep += [int(i), int(ni[j])]
    m = np.zeros(len(log_len), bool)
    m[keep] = True
    return m


def balanced(idx: list[int], donors: list[str], n_max: int, per_donor: int, rng: random.Random) -> list[int]:
    """run_isp.balanced_subset と同じ規則（donor あたり per_donor 以下・合計 n_max 以下）。"""
    order = list(idx)
    rng.shuffle(order)
    used: dict[str, int] = {}
    keep = []
    for i in order:
        dn = donors[i]
        if per_donor is not None and used.get(dn, 0) >= per_donor:
            continue
        used[dn] = used.get(dn, 0) + 1
        keep.append(i)
        if n_max is not None and len(keep) >= n_max:
            break
    return sorted(keep)


def state_rows(ds, seeds: list[int]) -> dict:
    """seed ごとに、状態重心に使う dataset の行番号を決める。
    start / goal / alt をそれぞれ 800 細胞・donor 60 で抽出したうえで、
    goal（と alt）は start と log(トークン長)の 1:1 マッチが付いた細胞だけを残す。"""
    key = "fib_tls_state"
    states = np.array([str(x) for x in ds[key]])
    donors = [str(x) for x in ds["donor_id"]]
    length = np.array(ds["length"], dtype=float)
    log_len = np.log(length)
    out = {}
    for s in seeds:
        rng = random.Random(s)
        cap = {}
        for st in ["FIB_neg_TLShigh", "FIB_org", "FIB_neg_TLSlow"]:
            cap[st] = balanced(np.flatnonzero(states == st).tolist(), donors,
                               STATE_MAX, STATE_PER_DONOR, rng)
        rows: dict[str, list[int]] = {}
        info: dict[str, dict] = {}
        for st in ["FIB_org", "FIB_neg_TLSlow"]:
            m0 = np.zeros(len(states), bool)
            m0[cap[st]] = True
            m1 = np.zeros(len(states), bool)
            m1[cap["FIB_neg_TLShigh"]] = True
            m = depth_match(log_len, m0, m1, DEPTH_TOL)
            rows[st] = sorted(np.flatnonzero(m & m0).tolist())
            rows.setdefault("_start_" + st, sorted(np.flatnonzero(m & m1).tolist()))
            info[st] = {"n_matched": len(rows[st]),
                        "median_length": {st: float(np.median(length[rows[st]])),
                                          "FIB_neg_TLShigh": float(np.median(length[rows["_start_" + st]]))},
                        "n_before_match": len(cap[st])}
        # start 側の重心は goal とのマッチで使った start 細胞（深さをそろえた対）を使う
        rows["FIB_neg_TLShigh"] = rows.pop("_start_FIB_org")
        rows.pop("_start_FIB_neg_TLSlow")
        info["FIB_neg_TLShigh"] = {"n_matched": len(rows["FIB_neg_TLShigh"]),
                                   "median_length": float(np.median(length[rows["FIB_neg_TLShigh"]])),
                                   "n_before_match": len(cap["FIB_neg_TLShigh"])}
        out[s] = {"rows": rows, "info": info,
                  "unmatched_goal_rows": sorted(cap["FIB_org"])}
    return out


# ------------------------------------------------------------------ (D) 実行可能性
def coverage(ds, d: pd.DataFrame, genes: dict[str, str]) -> dict:
    import scipy.sparse as sp
    key = "fib_tls_state"
    states = np.array([str(x) for x in ds[key]])
    donors_all = [str(x) for x in ds["donor_id"]]
    ids = list(ds["input_ids"])
    start_rows = np.flatnonzero(states == "FIB_neg_TLShigh")
    sel = [ids[i] for i in start_rows]
    donors = [donors_all[i] for i in start_rows]
    lens = np.array([len(x) for x in sel], float)
    tok = dict(zip(d.symbol, d.token))
    indptr = np.cumsum([0] + [len(x) for x in sel])
    col = np.concatenate([np.asarray(x) for x in sel])
    ntok = int(col.max()) + 1
    m = sp.csr_matrix((np.ones(len(col), bool), col, indptr), shape=(len(sel), ntok)).tocsc()

    def pick(sym: str, seed: int) -> list[int]:
        t = int(tok[sym])
        idx = m[:, t].nonzero()[0].tolist() if t < ntok else []
        return balanced(idx, donors, MAX_CELLS_PER_GENE, MAX_CELLS_PER_DONOR,
                        random.Random(f"{seed}-{genes[sym]}"))

    def prof(sym: str, keep: list[int]):
        t = int(tok[sym])
        return (np.log2(1 + np.array([list(sel[i]).index(t) for i in keep], float)),
                np.log2(lens[keep]))

    P, ncell, ndon = {}, {}, {}
    for s in SEEDS:
        for g in genes:
            k = pick(g, s)
            ncell[f"{g}@{s}"] = len(k)
            ndon[f"{g}@{s}"] = len({donors[i] for i in k})
            P[(s, g)] = prof(g, k)
    nulls = [g for g, r in genes.items() if r == "negative"] if isinstance(next(iter(genes.values())), str) else []
    return {"ncell": ncell, "ndon": ndon, "_P": P}


def eligible(P: dict, genes: list[str], nulls: list[str]) -> dict:
    out = {}
    for s in SEEDS:
        for g in genes:
            lp, ll = P[(s, g)]
            if len(lp) == 0:
                out[f"{g}@{s}"] = 0
                continue
            n = 0
            for h in nulls:
                if h == g:
                    continue
                hp, hl = P[(s, h)]
                if len(hp) == 0:
                    continue
                dd = np.maximum(np.abs(lp[:, None] - hp[None, :]), np.abs(ll[:, None] - hl[None, :]))
                n += (dd <= PERM_TOL).any(axis=1).mean() >= PERM_COVER
            out[f"{g}@{s}"] = int(n)
    return out


# ------------------------------------------------------------------ main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true", help="文献確認のため候補を距離順に表示するだけ")
    a = ap.parse_args()
    d = load_stats()
    cand, comp = comparators(d, a.list)
    with pd.option_context("display.width", 260, "display.max_rows", 400):
        print(cand[["ref", "symbol", "decision", "dist", "log2fc_start_vs_goal", "n_start", "n_donor_start",
                    "relpos_med", "pos_med"]].round(3).to_string(index=False))
    if a.list:
        return

    nulls = null_pool(d, comp)
    null_syms = [s for s, _ in nulls]
    assert len(set(null_syms)) == len(null_syms)
    assert not (set(null_syms) & (set(comp) | set(FIXED_HIGH) | set(TARGET_REF) | BLIND | LABEL_GENES))
    cand.to_csv(E / "work/comparator_candidates.csv", index=False, float_format="%.5g")
    cand[cand.decision == "include"].to_csv(E / "work/comparators.csv", index=False, float_format="%.5g")

    t = d.set_index("symbol")
    (E / "work/null_pool.json").write_text(json.dumps({
        "rule": "固定遺伝子ごとに d_n（発現頻度・元の位置・入力の長さ）の近い 3 遺伝子。|log2FC| < 1、自動除外・文献確認で除外したものを除く",
        "null_exclude": NULL_EXCLUDE, "null_exclude_prefix": list(NULL_EXCLUDE_PREFIX),
        "tolerances": {k: float(v) for k, v in TOL_N.items()},
        "n": len(nulls),
        "genes": [{"symbol": s, "ensembl": t.loc[s, "ensembl"], "matched_to": r,
                   "n_start": int(t.loc[s, "n_start"]), "n_donor_start": int(t.loc[s, "n_donor_start"]),
                   "pos_med": float(t.loc[s, "pos_med"]), "len_med": float(t.loc[s, "len_med"]),
                   "log2fc": round(float(t.loc[s, "log2fc_start_vs_goal"]), 3)} for s, r in nulls],
        "uses_isp_results": False}, ensure_ascii=False, indent=1) + "\n")

    from datasets import load_from_disk
    ds = load_from_disk(str(TOK))

    sr = state_rows(ds, SEEDS)
    for s in SEEDS:
        (E / f"work/state_rows-seed{s}.json").write_text(json.dumps(
            {"seed": s, "depth_match": {"feature": "length", "log_tolerance": DEPTH_TOL},
             "cap": {"max_cells_per_state": STATE_MAX, "max_cells_per_donor": STATE_PER_DONOR},
             "rows": sr[s]["rows"], "unmatched_goal_rows": sr[s]["unmatched_goal_rows"],
             "info": sr[s]["info"], "uses_isp_results": False}, ensure_ascii=False) + "\n")
    (E / "work/depth_match_summary.json").write_text(json.dumps(
        {str(s): sr[s]["info"] for s in SEEDS}, ensure_ascii=False, indent=1) + "\n")

    role = {}
    for g in TARGET_MAIN:
        role[g] = "target"
    for g in TARGET_REF:
        role[g] = "target_ref"
    for g in HK:
        role[g] = "housekeeping"
    for g in PANEL:
        role[g] = "negctrl_panel"
    for g in FIBRO_DIR:
        role[g] = "fibrosis_direction"
    for g in comp:
        role[g] = "comparator"
    for g in null_syms:
        role[g] = "negative"
    sym2ens = dict(zip(d.symbol, d.ensembl))
    genes = {g: sym2ens[g] for g in role}

    cov = coverage(ds, d, genes)
    P = cov.pop("_P")
    elig = eligible(P, list(role), null_syms)
    grp: dict[str, list[str]] = {}
    for g, r in role.items():
        grp.setdefault(r, []).append(g)
    out = {
        "seeds": SEEDS, "tolerance_log2": PERM_TOL, "min_cover_frac": PERM_COVER,
        "n_genes": len(role), "roles": {k: len(v) for k, v in grp.items()},
        "cells_by_gene_seed": cov["ncell"], "donors_by_gene_seed": cov["ndon"],
        "cells_min": int(min(cov["ncell"].values())),
        "cells_below_min": {k: v for k, v in cov["ncell"].items() if v < MIN_CELLS_PER_GENE},
        "eligible_null": elig,
        "eligible_null_min_by_role": {k: int(min(elig[f"{g}@{s}"] for g in v for s in SEEDS))
                                      for k, v in grp.items()},
        "eligible_null_median_by_role": {k: float(np.median([elig[f"{g}@{s}"] for g in v for s in SEEDS]))
                                         for k, v in grp.items()},
        "eligible_null_below_min": {k: v for k, v in elig.items() if v < PERM_MIN_NULL},
        "total_tokens": None, "uses_isp_results": False,
    }
    # 計算量見積もり用の総トークン数
    states = np.array([str(x) for x in ds["fib_tls_state"]])
    ids = list(ds["input_ids"])
    start_rows = np.flatnonzero(states == "FIB_neg_TLShigh")
    lens = np.array([len(ids[i]) for i in start_rows], float)
    tokd = dict(zip(d.symbol, d.token))
    import scipy.sparse as sp
    sel = [ids[i] for i in start_rows]
    indptr = np.cumsum([0] + [len(x) for x in sel])
    col = np.concatenate([np.asarray(x) for x in sel])
    m = sp.csr_matrix((np.ones(len(col), bool), col, indptr),
                      shape=(len(sel), int(col.max()) + 1)).tocsc()
    donors_all = [str(x) for x in ds["donor_id"]]
    donors = [donors_all[i] for i in start_rows]
    tot = 0
    for s in SEEDS:
        for g in role:
            idx = m[:, int(tokd[g])].nonzero()[0].tolist()
            k = balanced(idx, donors, MAX_CELLS_PER_GENE, MAX_CELLS_PER_DONOR, random.Random(f"{s}-{genes[g]}"))
            tot += float(lens[k].sum())
    out["total_tokens"] = tot
    (E / "work/precheck_coverage.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    print("comparators:", comp)
    print("nulls:", null_syms)
    print({k: out[k] for k in ["n_genes", "roles", "cells_min", "cells_below_min",
                               "eligible_null_min_by_role", "eligible_null_median_by_role",
                               "eligible_null_below_min", "total_tokens"]})
    print("depth match:", json.dumps({str(s): sr[s]["info"] for s in SEEDS}, ensure_ascii=False))


if __name__ == "__main__":
    main()

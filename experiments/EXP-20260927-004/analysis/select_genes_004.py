"""EXP-20260927-004 事前登録の一部（摂動の実行前に確定する）。

摂動（ISP）の結果は一切読まない。使うのは
  - work/gene_stats_fib.csv（EXP-20260927-003 が生カウントとトークン化済み入力から作った遺伝子統計表）
  - トークン化済み入力（遺伝子の位置・入力の長さ・提供者・状態ラベル）
  - 文献（outputs/literature/TLS-organizer-drivers-2026-09-27.md）
だけである。

(A) 目標の軸の陽性対照候補 10 遺伝子（文献側 5 ＋ データ側 5）を固定する
(B) 新しい陽性対照候補に、移動距離（元の位置・入力の長さ）と細胞数をそろえた対照遺伝子を 3 つずつ選ぶ
(C) 深さをそろえた状態重心に使う細胞（seed 31 / 32 / 33）
(D) 実行可能性（遺伝子ごとの細胞数・提供者数、主検定で「使える対照遺伝子」の数）

出力: work/{positive_controls_004.csv, null_pool_004.json, state_rows-seed{S}.json,
        depth_match_summary_004.json, precheck_coverage_004.json, genes_004.json}
実行: .venv/bin/python experiments/EXP-20260927-004/analysis/select_genes_004.py
"""
from __future__ import annotations

import importlib.util
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260927-004"
SRC3 = ROOT / "experiments/EXP-20260927-003"
STATS = SRC3 / "work/gene_stats_fib.csv"          # 読み取りのみ（EXP-003 の事前登録済みファイル）
TOK = E / "work/tokenized.dataset"                # EXP-003 のデータ部分を cp -c で複製したもの

SEEDS = [31, 32, 33]                              # 既出（0〜5, 11〜13, 21〜23）と重ならない
MAX_CELLS_PER_GENE, MAX_CELLS_PER_DONOR, MIN_CELLS_PER_GENE = 200, 20, 120
N_NULL_PER_PC = 3
# 主検定で「使える対照遺伝子」の下限。下回る遺伝子×seed は判定不能（NaN）になるので、
# 計画の段階で全遺伝子がこれを満たすまで対照遺伝子を補充する（2026-09-27 追加。結果は一切見ていない）
MIN_ELIG_NULLS = 20
N_AUG_SEARCH = 400                                # 補充の探索幅（不足している遺伝子ごとに、距離の近い順に何個まで調べるか）
N_AUG_ROUNDS = 8
TOL_N = {"lncells": np.log(1.5), "log2pos": 0.25, "log2len": 0.25}
NULL_MAX_ABS_LOG2FC = 1.0
MIN_CELLS, MIN_DONORS = 250, 12                   # 対照遺伝子の候補条件（1 遺伝子 200 細胞を確保できる水準）
PERM_TOL, PERM_COVER = 0.5, 0.8

# ---- EXP-20260927-003 の事前登録済みコードから、深さマッチと抽出規則をそのまま使う（実装のずれを避ける）
_spec = importlib.util.spec_from_file_location("sg3", SRC3 / "analysis/select_genes_fib.py")
sg3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sg3)
assert sg3.DEPTH_TOL == 0.10 and sg3.STATE_MAX == 800 and sg3.STATE_PER_DONOR == 60

# ---- この実験で摂動しない遺伝子（読み込み時に assert する）
BLOCKED = {"LTBR": "主対象（HYP-20260925-02）。本実験では回さない",
           "RELB": "主対象。本実験では回さない",
           "NFKB2": "主対象。本実験では回さない",
           "CDKN1A": "参考対象・盲検。本実験では回さない",
           "TP53": "盲検（未解除）", "SPP1": "盲検（未解除）", "LCN2": "盲検（未解除）",
           "CCL19": "目標ラベルを定義する遺伝子（循環性）",
           "CCL21": "目標ラベルを定義する遺伝子（循環性）",
           "CXCL13": "目標ラベルを定義する遺伝子（循環性）"}

# ---- (A) 目標の軸の陽性対照候補（結果を見る前に固定）
# tier=lit : 文献で三次リンパ組織オーガナイザー／線維芽細胞細網細胞の形質として報告のある分子（固定リスト）
# tier=data: データだけから機械的に選ぶマーカー操作型。選定規則は下の DATA_TIER_RULE のとおり。
#            文献に「線維芽細胞で単一遺伝子を過剰発現してオーガナイザー様へ転換させた」報告が無い
#            （outputs/literature/TLS-organizer-drivers-2026-09-27.md）ため、代替として置く。
DATA_TIER_RULE = (
    "log2 発現比 ≤ −1.0（目標の状態に偏る）・出発点 250 細胞以上・提供者 12 人以上を満たす遺伝子から、"
    "(i) 目標の状態での検出率が最大のもの 1 つ（＝目標の細胞の過半に載る代表的マーカー）と、"
    "(ii) 残りを log2 発現比の大きい順に並べ、実行可能性の条件を満たすものを上から 4 つ、を採る。"
    "実行可能性の条件 = 抽出規則で 120 細胞以上取れ、かつ対照遺伝子の補充後に主検定で使える対照遺伝子が 20 以上。"
)
DATA_TIER_FC_MAX, DATA_TIER_N_FC, DATA_TIER_N_FRAC = -1.0, 4, 1
POSITIVE_CANDIDATES_LIT: list[tuple[str, str, str, str]] = [
    ("VCAM1", "lit", "抗 LTβR アゴニストで前駆脂肪細胞がオーガナイザー様になるとき誘導される接着分子の 1 つ（PMID 22940098）",
     "非古典的 NF-κB の下流でもあるため、本仮説の軸と部分的に重なる"),
    ("ICAM1", "lit", "同上（PMID 22940098）。線維芽細胞細網細胞でも高い（PMID 18941209）",
     "同上（NF-κB 標的）"),
    ("CXCL12", "lit", "線維芽細胞細網細胞・リンパ組織間質の代表的ケモカイン（PMID 22807928, 33981032）",
     "目標ラベルを定義する 3 遺伝子とは別のケモカインだが、ケモカイン一般の共通性は残る"),
    ("IL7", "lit", "線維芽細胞細網細胞の必須因子。LTβR 欠損で産生が失われる（PMID 23623380）",
     "LTβR 下流であり本仮説の軸と部分的に重なる"),
    ("TNFSF13B", "lit", "三次リンパ組織の間質が出す B 細胞生存因子（BAFF）",
     "非古典的 NF-κB の標的でもある"),
]
# data 層の選定で落ちた遺伝子は work/positive_controls_004.csv の rejected 欄に残す。
# 参考（2026-09-27 時点の実測、摂動の結果は見ていない）:
#   CA8   log2 発現比 −1.97 で偏りは最大だが、抽出規則（1 遺伝子 200 細胞・1 提供者 20 細胞・
#         出発点 326 細胞 / 16 提供者）で 111 細胞しか取れず、下限 120 を満たさない。
#   HPSE2 log2 発現比 −1.35 で順位が 10 位のため、上位 5 つが埋まった時点で選ばれない。
#         仮に選んだとしても、細胞は取れる（193 細胞）が、元の位置の広がりが大きすぎて
#         （log2(1+順位) が 5〜95 パーセンタイルで 4.1〜11.4）、80% を許容 0.5 で覆える対照遺伝子が
#         補充しても 11〜15 個しか見つからず、主検定の下限 20 に届かない＝構造的に判定不能になる。

# ---- 既存の役割（EXP-20260927-003 と同じ遺伝子。比較可能性のため変えない）
HK = {"ACTB": "ENSG00000075624", "GAPDH": "ENSG00000111640", "PPIA": "ENSG00000196262"}
PANEL = {"CERS4": "ENSG00000090661", "TATDN1": "ENSG00000147687", "IFT25": "ENSG00000081870",
         "VAPA": "ENSG00000101558", "PHPT1": "ENSG00000054148", "ANAPC16": "ENSG00000166295",
         "EEF1D": "ENSG00000104529"}
POS_MYOF = {"ACTA2": "ENSG00000107796", "TAGLN": "ENSG00000149591"}   # 筋線維芽細胞の軸の陽性対照（参照用）


def load_stats() -> pd.DataFrame:
    d = pd.read_csv(STATS)
    d = d[d.symbol.astype(str) != ""].copy()
    d["l2pos"] = np.log2(1 + d.pos_med)
    d["l2len"] = np.log2(d.len_med)
    return d


def existing_nulls(d: pd.DataFrame) -> list[dict]:
    """EXP-20260927-003 の対照遺伝子 48 のうち、出発点で 250 細胞以上のものを引き継ぐ。
    低頻度の 3 遺伝子（CDKN1A 用）は、本実験に CDKN1A がないので外す（細胞数をそろえるため）。"""
    pool = json.loads((SRC3 / "work/null_pool.json").read_text())["genes"]
    keep = [g for g in pool if g["n_start"] >= MIN_CELLS and g["n_donor_start"] >= MIN_DONORS]
    dropped = [g["symbol"] for g in pool if g not in keep]
    print(f"既存の対照遺伝子 {len(pool)} → 引き継ぎ {len(keep)}（除外 {dropped}）")
    return keep


def new_nulls(d: pd.DataFrame, pcs: list[str], taken: set[str]) -> list[tuple[str, str]]:
    """陽性対照候補ごとに、細胞数・元の位置・入力の長さの近い対照遺伝子を 3 つ選ぶ。
    除外規則は EXP-20260927-003 と同じものをそのまま使う（sg3 の定数）。"""
    ref = d.set_index("symbol")
    auto = (d.symbol.isin(sg3.EXCLUDE_SYMBOLS | sg3.BLIND | sg3.LABEL_GENES | set(BLOCKED)
                          | set(HK) | set(PANEL) | set(POS_MYOF) | set(pcs))
            | d.symbol.str.startswith(sg3.EXCLUDE_PREFIX)
            | d.symbol.isin(set(sg3.NULL_EXCLUDE))
            | d.symbol.str.startswith(sg3.NULL_EXCLUDE_PREFIX)
            | d.symbol.isin({k for k, v in sg3.MANUAL_REVIEW.items() if v[0] == "exclude"}))
    base = d[~auto & (d.log2fc_start_vs_goal.abs() < NULL_MAX_ABS_LOG2FC)
             & (d.n_start >= MIN_CELLS) & (d.n_donor_start >= MIN_DONORS)]
    chosen: list[tuple[str, str]] = []
    for r in pcs:
        got = taken | {s for s, _ in chosen}
        x = base[~base.symbol.isin(got)].copy()
        x["dist"] = pd.concat([
            np.log(x.n_start / ref.loc[r, "n_start"]).abs() / TOL_N["lncells"],
            (x.l2pos - ref.loc[r, "l2pos"]).abs() / TOL_N["log2pos"],
            (x.l2len - ref.loc[r, "l2len"]).abs() / TOL_N["log2len"]], axis=1).max(1)
        sel = x.sort_values(["dist", "symbol"]).head(N_NULL_PER_PC)
        if len(sel) < N_NULL_PER_PC:
            raise SystemExit(f"{r}: 対照遺伝子が {len(sel)} 個しか選べない")
        for s in sel.symbol:
            chosen.append((s, r))
    return chosen


class Profiles:
    """出発点（FIB_neg_TLShigh）の細胞から、遺伝子 × seed ごとに抽出される細胞と
    （log2(1+元の位置), log2(入力の長さ)）を返す。抽出規則は EXP-20260927-003 の balanced() をそのまま使う。"""

    def __init__(self, ds, d: pd.DataFrame) -> None:
        import scipy.sparse as sp
        states = np.array([str(x) for x in ds["fib_tls_state"]])
        donors_all = [str(x) for x in ds["donor_id"]]
        ids = list(ds["input_ids"])
        start_rows = np.flatnonzero(states == "FIB_neg_TLShigh")
        self.sel = [list(ids[i]) for i in start_rows]
        self.donors = [donors_all[i] for i in start_rows]
        self.lens = np.array([len(x) for x in self.sel], float)
        self.tok = dict(zip(d.symbol, d.token))
        self.eid = dict(zip(d.symbol, d.ensembl))
        indptr = np.cumsum([0] + [len(x) for x in self.sel])
        col = np.concatenate([np.asarray(x) for x in self.sel])
        self.ntok = int(col.max()) + 1
        self.m = sp.csr_matrix((np.ones(len(col), bool), col, indptr),
                               shape=(len(self.sel), self.ntok)).tocsc()
        self._cells: dict[tuple[str, int], tuple] = {}
        self._pair: dict[tuple[str, str, int], bool] = {}

    def get(self, sym: str, seed: int) -> tuple[list[int], np.ndarray, np.ndarray]:
        key = (sym, seed)
        if key not in self._cells:
            t = int(self.tok[sym])
            idx = self.m[:, t].nonzero()[0].tolist() if t < self.ntok else []
            keep = sg3.balanced(idx, self.donors, MAX_CELLS_PER_GENE, MAX_CELLS_PER_DONOR,
                                random.Random(f"{seed}-{self.eid[sym]}"))
            lp = np.log2(1 + np.array([self.sel[i].index(t) for i in keep], float))
            ll = np.log2(self.lens[keep]) if len(keep) else np.zeros(0)
            self._cells[key] = (keep, lp, ll)
        return self._cells[key]

    def covers(self, gene: str, null: str, seed: int) -> bool:
        """null の細胞が gene の細胞の 80% 以上を、log2 距離 0.5 以内で覆えるか（主検定の使用条件）。"""
        key = (gene, null, seed)
        if key not in self._pair:
            _, lp, ll = self.get(gene, seed)
            _, hp, hl = self.get(null, seed)
            if len(lp) == 0 or len(hp) == 0:
                self._pair[key] = False
            else:
                dd = np.maximum(np.abs(lp[:, None] - hp[None, :]), np.abs(ll[:, None] - hl[None, :]))
                self._pair[key] = bool((dd <= PERM_TOL).any(axis=1).mean() >= PERM_COVER)
        return self._pair[key]

    def n_eligible(self, gene: str, nulls: list[str], seed: int) -> int:
        return sum(self.covers(gene, h, seed) for h in nulls if h != gene)


def augment_nulls(prof: Profiles, d: pd.DataFrame, all_genes: list[str],
                  nulls: list[str], taken: set[str], base: pd.DataFrame) -> list[tuple[str, str]]:
    """主検定で使える対照遺伝子が MIN_ELIG_NULLS を下回る遺伝子に、
    細胞数・元の位置・入力の長さの近い対照遺伝子を足す（結果は一切使わない規則）。"""
    ref = d.set_index("symbol")
    added: list[tuple[str, str]] = []
    for rnd in range(N_AUG_ROUNDS):
        short = {g: min(prof.n_eligible(g, nulls, s) for s in SEEDS) for g in all_genes}
        short = {g: v for g, v in short.items() if v < MIN_ELIG_NULLS}
        if not short:
            print(f"補充 {rnd} 巡目: 不足なし")
            return added
        print(f"補充 {rnd} 巡目: 不足 {len(short)} 遺伝子 " +
              ", ".join(f"{g}({v})" for g, v in sorted(short.items(), key=lambda x: x[1])))
        got = taken | {s for s, _ in added}
        progressed = False
        for g in sorted(short, key=lambda x: short[x]):
            x = base[~base.symbol.isin(got)].copy()
            if x.empty:
                break
            x["dist"] = pd.concat([
                np.log(x.n_start / ref.loc[g, "n_start"]).abs() / TOL_N["lncells"],
                (x.l2pos - ref.loc[g, "l2pos"]).abs() / TOL_N["log2pos"],
                (x.l2len - ref.loc[g, "l2len"]).abs() / TOL_N["log2len"]], axis=1).max(1)
            need = MIN_ELIG_NULLS - min(prof.n_eligible(g, nulls + [s for s, _ in added], s) for s in SEEDS)
            for sym in x.sort_values(["dist", "symbol"]).symbol.head(N_AUG_SEARCH):
                if need <= 0:
                    break
                if all(prof.covers(g, sym, s) for s in SEEDS) and len(prof.get(sym, SEEDS[0])[0]) >= MIN_CELLS_PER_GENE:
                    added.append((sym, g))
                    got.add(sym)
                    need -= 1
                    progressed = True
        nulls = nulls + [s for s, _ in added]
        if not progressed:
            print("補充できる対照遺伝子が尽きた")
            return added
    return added


def try_feasible(prof: Profiles, d: pd.DataFrame, sym: str, nulls: list[str],
                 base: pd.DataFrame, taken: set[str]) -> tuple[bool, list[str], str]:
    """陽性対照候補 1 つの実行可能性を調べる。
    返り値: (可否, その遺伝子のために足す必要のある対照遺伝子, 不可の理由)。摂動の結果は一切使わない。"""
    ns = [len(prof.get(sym, s)[0]) for s in SEEDS]
    if min(ns) < MIN_CELLS_PER_GENE:
        return False, [], f"抽出できる細胞が {min(ns)} で下限 {MIN_CELLS_PER_GENE} 未満"
    ref = d.set_index("symbol")
    pool, add = list(nulls), []
    got = set(taken) | set(nulls) | {sym}
    x = base[~base.symbol.isin(got)].copy()
    x["dist"] = pd.concat([
        np.log(x.n_start / ref.loc[sym, "n_start"]).abs() / TOL_N["lncells"],
        (x.l2pos - ref.loc[sym, "l2pos"]).abs() / TOL_N["log2pos"],
        (x.l2len - ref.loc[sym, "l2len"]).abs() / TOL_N["log2len"]], axis=1).max(axis=1)
    for cand in x.sort_values(["dist", "symbol"]).symbol.head(N_AUG_SEARCH):
        if min(prof.n_eligible(sym, pool, s) for s in SEEDS) >= MIN_ELIG_NULLS:
            break
        if all(prof.covers(sym, cand, s) for s in SEEDS) and \
                min(len(prof.get(cand, s)[0]) for s in SEEDS) >= MIN_CELLS_PER_GENE:
            pool.append(cand); add.append(cand)
    worst = min(prof.n_eligible(sym, pool, s) for s in SEEDS)
    if worst < MIN_ELIG_NULLS:
        return False, [], f"補充しても使える対照遺伝子が {worst} で下限 {MIN_ELIG_NULLS} 未満（構造的に判定不能）"
    return True, add, ""


def pick_data_tier(prof: Profiles, d: pd.DataFrame, lit: list[str], nulls: list[str],
                   base: pd.DataFrame, taken: set[str]) -> tuple[list[tuple[str, str]], list[str], list[dict]]:
    """DATA_TIER_RULE に従って、データだけからマーカー操作型の陽性対照を選ぶ。"""
    cand = d[(d.log2fc_start_vs_goal <= DATA_TIER_FC_MAX) & (d.n_start >= MIN_CELLS)
             & (d.n_donor_start >= MIN_DONORS)
             & ~d.symbol.isin(sg3.EXCLUDE_SYMBOLS | sg3.BLIND | sg3.LABEL_GENES | set(BLOCKED)
                              | set(HK) | set(PANEL) | set(POS_MYOF) | set(lit))
             & ~d.symbol.str.startswith(sg3.EXCLUDE_PREFIX)
             & ~d.symbol.isin({k for k, v in sg3.MANUAL_REVIEW.items() if v[0] == "exclude"})].copy()
    print(f"data 層の候補 {len(cand)}（log2 発現比 ≤ {DATA_TIER_FC_MAX}、出発点 {MIN_CELLS} 細胞以上、"
          f"提供者 {MIN_DONORS} 人以上）")
    chosen: list[tuple[str, str]] = []
    extra: list[str] = []
    rejected: list[dict] = []
    pool = list(nulls)
    took = set(taken)

    def consider(sym: str, why: str) -> bool:
        ok, add, reason = try_feasible(prof, d, sym, pool, base, took)
        if not ok:
            rejected.append({"symbol": sym, "reason": reason,
                             "log2fc_start_vs_goal": round(float(d.set_index("symbol").loc[sym,
                                                          "log2fc_start_vs_goal"]), 3)})
            print(f"  ✗ {sym}: {reason}")
            return False
        chosen.append((sym, why))
        pool.extend(add); extra.extend(add); took.add(sym); took.update(add)
        print(f"  ✓ {sym}（{why}）対照遺伝子の補充 {len(add)}")
        return True

    # (i) 目標の状態での検出率が最大のもの
    for sym in cand.sort_values("frac_goal", ascending=False).symbol:
        if consider(sym, f"目標の状態での検出率が最大（{float(cand.set_index('symbol').loc[sym, 'frac_goal']):.3f}）"):
            break
    # (ii) log2 発現比の大きい順
    for sym in cand.sort_values("log2fc_start_vs_goal").symbol:
        if len(chosen) >= DATA_TIER_N_FC + DATA_TIER_N_FRAC:
            break
        if sym in {s for s, _ in chosen}:
            continue
        consider(sym, f"目標の状態への偏りが大きい順（log2 発現比 "
                      f"{float(cand.set_index('symbol').loc[sym, 'log2fc_start_vs_goal']):.2f}）")
    return chosen, extra, rejected


def coverage(prof: Profiles, genes: dict[str, str], nulls: list[str]) -> dict:
    """遺伝子ごとに実際に抽出できる細胞数・提供者数と、主検定で使える対照遺伝子の数を数える。"""
    ncell, ndon, elig, ntoken = {}, {}, {}, 0.0
    for s in SEEDS:
        for g in genes:
            keep, _, _ = prof.get(g, s)
            ncell[f"{g}@{s}"] = len(keep)
            ndon[f"{g}@{s}"] = len({prof.donors[i] for i in keep})
            ntoken += float(prof.lens[keep].sum()) if len(keep) else 0.0
            elig[f"{g}@{s}"] = prof.n_eligible(g, nulls, s)
    return {"n_cells": ncell, "n_donors": ndon, "n_eligible_null_genes": elig,
            "total_cell_gene": int(sum(ncell.values())), "total_tokens": float(ntoken),
            "max_cells_per_gene": MAX_CELLS_PER_GENE, "max_cells_per_donor": MAX_CELLS_PER_DONOR,
            "min_cells_per_gene": MIN_CELLS_PER_GENE, "min_eligible_null_genes": MIN_ELIG_NULLS}


def main() -> None:
    (E / "work").mkdir(parents=True, exist_ok=True)
    d = load_stats()
    t = d.set_index("symbol")
    lit = [s for s, *_ in POSITIVE_CANDIDATES_LIT]
    assert not (set(lit) & set(BLOCKED)), "回さない遺伝子が陽性対照候補に入っている"

    def meta(sym: str, matched_to: str, source: str) -> dict:
        return {"symbol": sym, "ensembl": t.loc[sym, "ensembl"], "matched_to": matched_to, "source": source,
                "n_start": int(t.loc[sym, "n_start"]), "n_donor_start": int(t.loc[sym, "n_donor_start"]),
                "pos_med": float(t.loc[sym, "pos_med"]), "len_med": float(t.loc[sym, "len_med"]),
                "log2fc": round(float(t.loc[sym, "log2fc_start_vs_goal"]), 3)}

    # (B) 対照遺伝子（EXP-20260927-003 からの引き継ぎ）
    old = existing_nulls(d)
    null_meta: list[dict] = [
        {"symbol": g["symbol"], "ensembl": g["ensembl"], "matched_to": g["matched_to"],
         "source": "EXP-20260927-003", "n_start": g["n_start"], "n_donor_start": g["n_donor_start"],
         "pos_med": g["pos_med"], "len_med": g["len_med"], "log2fc": g["log2fc"]} for g in old]

    # (C) 深さをそろえた状態重心の細胞
    from datasets import load_from_disk
    ds = load_from_disk(str(TOK))
    sr = sg3.state_rows(ds, SEEDS)
    summary = {}
    for s in SEEDS:
        (E / f"work/state_rows-seed{s}.json").write_text(
            json.dumps(sr[s]["rows"], ensure_ascii=False) + "\n")
        summary[s] = sr[s]["info"]
        (E / f"work/unmatched_goal_rows-seed{s}.json").write_text(
            json.dumps(sr[s]["unmatched_goal_rows"]) + "\n")
    (E / "work/depth_match_summary_004.json").write_text(
        json.dumps({"depth_tolerance_log": sg3.DEPTH_TOL, "seeds": SEEDS, "info": summary},
                   ensure_ascii=False, indent=1, default=float) + "\n")
    print("深さマッチ:", {s: summary[s]["FIB_org"]["n_matched"] for s in SEEDS})

    # (A) 陽性対照候補: 文献層は固定、data 層は DATA_TIER_RULE で機械的に選ぶ
    prof = Profiles(ds, d)
    aug_base = d[~(d.symbol.isin(sg3.EXCLUDE_SYMBOLS | sg3.BLIND | sg3.LABEL_GENES | set(BLOCKED))
                   | d.symbol.str.startswith(sg3.EXCLUDE_PREFIX)
                   | d.symbol.isin(set(sg3.NULL_EXCLUDE))
                   | d.symbol.str.startswith(sg3.NULL_EXCLUDE_PREFIX)
                   | d.symbol.isin({k for k, v in sg3.MANUAL_REVIEW.items() if v[0] == "exclude"}))
                & (d.log2fc_start_vs_goal.abs() < NULL_MAX_ABS_LOG2FC)
                & (d.n_start >= MIN_CELLS) & (d.n_donor_start >= MIN_DONORS)]
    fixed_roles = set(HK) | set(PANEL) | set(POS_MYOF) | set(lit)
    nulls0 = [g["symbol"] for g in null_meta]
    dsel, dextra, drej = pick_data_tier(prof, d, lit, nulls0, aug_base,
                                        set(nulls0) | fixed_roles)
    for s in dextra:
        null_meta.append(meta(s, "data 層の選定", "screen"))
    pcs = lit + [s for s, _ in dsel]
    assert not (set(pcs) & set(BLOCKED)), "回さない遺伝子が陽性対照候補に入っている"
    tiers = {s: ("lit", why, lim) for s, _, why, lim in POSITIVE_CANDIDATES_LIT}
    tiers.update({s: ("data", why, "三次リンパ組織との機能的関係は知られていない。純粋なマーカー操作")
                  for s, why in dsel})
    rows = []
    for s in pcs:
        r = t.loc[s]
        tier, why, lim = tiers[s]
        rows.append({"symbol": s, "ensembl": r.ensembl, "tier": tier, "reason": why, "limitation": lim,
                     "n_start": int(r.n_start), "n_donor_start": int(r.n_donor_start),
                     "n_goal": int(r.n_goal), "frac_goal": round(float(r.frac_goal), 4),
                     "log2fc_start_vs_goal": round(float(r.log2fc_start_vs_goal), 3),
                     "pos_med": float(r.pos_med), "len_med": float(r.len_med),
                     "relpos_med": round(float(r.relpos_med), 3)})
    pd.DataFrame(rows).to_csv(E / "work/positive_controls_004.csv", index=False)
    (E / "work/positive_controls_rejected_004.json").write_text(
        json.dumps({"rule": DATA_TIER_RULE, "rejected": drej, "uses_isp_results": False},
                   ensure_ascii=False, indent=1) + "\n")
    print(pd.DataFrame(rows)[["symbol", "tier", "n_start", "n_donor_start",
                              "frac_goal", "log2fc_start_vs_goal", "pos_med", "len_med"]].to_string(index=False))

    # (B2) 陽性対照候補に細胞数・元の位置・入力の長さの近い対照遺伝子を 3 つずつ足す
    nulls1 = [g["symbol"] for g in null_meta]
    new = new_nulls(d, pcs, set(nulls1) | fixed_roles | set(pcs))
    for s, r in new:
        null_meta.append(meta(s, r, "new"))
    print(f"対照遺伝子 {len(null_meta)}（引き継ぎ {len(old)} ＋ 選定時の補充 {len(dextra)} ＋ 新規 {len(new)}）")

    # (B3) 使える対照遺伝子が 20 を下回る「判定対象の遺伝子」に対照遺伝子を補充する。
    # 対照遺伝子自身が 20 を下回る場合は、その遺伝子×seed を較正の欠測として扱う
    # （guardrails 2026-09-27「欠測を不合格と数えない」）。判定対象ではないので補充の対象にしない。
    judged = list(pcs) + list(POS_MYOF) + list(HK) + list(PANEL)
    roles_fixed = judged
    nulls = [g["symbol"] for g in null_meta]
    extra = augment_nulls(prof, d, judged, nulls, set(nulls) | set(roles_fixed), aug_base)
    for s, r in extra:
        null_meta.append(meta(s, r, "augment"))
    syms = [g["symbol"] for g in null_meta]
    assert len(set(syms)) == len(syms), "対照遺伝子が重複している"
    assert not (set(syms) & (set(pcs) | set(HK) | set(PANEL) | set(POS_MYOF) | set(BLOCKED)))
    (E / "work/null_pool_004.json").write_text(json.dumps(
        {"rule": "EXP-20260927-003 の対照遺伝子のうち出発点 250 細胞以上のものを引き継ぎ、"
                 "陽性対照候補それぞれに細胞数・元の位置・入力の長さの近い 3 遺伝子を足し、"
                 "そのうえで主検定で使える対照遺伝子が 20 未満になる遺伝子には、"
                 "同じ距離の近い順に対照遺伝子を補充する（摂動の結果は一切使わない）",
         "tolerances": {k: float(v) for k, v in TOL_N.items()},
         "min_eligible_null_genes": MIN_ELIG_NULLS,
         "n": len(null_meta), "n_from_exp003": len(old), "n_screen": len(dextra),
         "n_new": len(new), "n_augment": len(extra),
         "genes": null_meta, "uses_isp_results": False}, ensure_ascii=False, indent=1) + "\n")
    print(f"対照遺伝子 {len(null_meta)}（引き継ぎ {len(old)} ＋ 選定時 {len(dextra)} ＋ "
          f"新規 {len(new)} ＋ 補充 {len(extra)}）")

    # 全遺伝子（役割つき）
    genes: dict[str, str] = {}
    role: dict[str, str] = {}
    for s in pcs:
        genes[s] = t.loc[s, "ensembl"]; role[s] = "positive_goal"
    for s, eid in POS_MYOF.items():
        genes[s] = eid; role[s] = "positive_myof"
    for s, eid in HK.items():
        genes[s] = eid; role[s] = "housekeeping"
    for s, eid in PANEL.items():
        genes[s] = eid; role[s] = "negctrl_panel"
    for g in null_meta:
        genes[g["symbol"]] = g["ensembl"]; role[g["symbol"]] = "negative"
    assert not (set(genes) & set(BLOCKED)), "回さない遺伝子が遺伝子リストに入っている"
    print(f"摂動する遺伝子 {len(genes)}")

    # (D) 実行可能性
    cov = coverage(prof, genes, [g["symbol"] for g in null_meta])
    (E / "work/precheck_coverage_004.json").write_text(
        json.dumps(cov, ensure_ascii=False, indent=1, default=float) + "\n")
    nc = pd.Series(cov["n_cells"])
    ne = pd.Series(cov["n_eligible_null_genes"])
    print(f"細胞数 中央値 {nc.median():.0f}、最小 {nc.min():.0f}（下限 {MIN_CELLS_PER_GENE}）")
    print(f"使える対照遺伝子 最小 {ne.min()}（下限 {MIN_ELIG_NULLS}）、中央値 {ne.median():.0f}")
    print(f"延べ {cov['total_cell_gene']} 細胞×遺伝子、総トークン {cov['total_tokens']:.3e}")
    judged_keys = [f"{g}@{s}" for g in judged for s in SEEDS]
    ne_j = ne[judged_keys]
    null_keys = [f"{g}@{s}" for g in [x["symbol"] for x in null_meta] for s in SEEDS]
    n_valid_null = int((ne[null_keys] >= MIN_ELIG_NULLS).sum())
    print(f"判定対象 {len(judged)} 遺伝子: 使える対照遺伝子 最小 {ne_j.min()}、中央値 {ne_j.median():.0f}")
    print(f"較正に使える対照遺伝子×seed: {n_valid_null} / {len(null_keys)}（下限 60）")
    cov["n_judged_genes"] = len(judged)
    cov["min_eligible_among_judged"] = int(ne_j.min())
    cov["n_valid_null_gene_seed"] = n_valid_null
    (E / "work/precheck_coverage_004.json").write_text(
        json.dumps(cov, ensure_ascii=False, indent=1, default=float) + "\n")
    assert nc.min() >= MIN_CELLS_PER_GENE, f"細胞数が下限を下回る遺伝子がある: {nc[nc < MIN_CELLS_PER_GENE].to_dict()}"
    assert ne_j.min() >= MIN_ELIG_NULLS, \
        f"判定対象で使える対照遺伝子が下限を下回る: {ne_j[ne_j < MIN_ELIG_NULLS].to_dict()}"
    assert n_valid_null >= 60, f"較正に使える対照遺伝子×seed が {n_valid_null} で下限 60 未満"

    (E / "work/genes_004.json").write_text(json.dumps(
        {"seeds": SEEDS, "blocked": BLOCKED,
         "genes": [{"symbol": s, "ensembl": e, "role": role[s]} for s, e in genes.items()],
         "counts": {r: sum(1 for v in role.values() if v == r) for r in sorted(set(role.values()))}},
        ensure_ascii=False, indent=1) + "\n")
    print(json.dumps({r: sum(1 for v in role.values() if v == r) for r in sorted(set(role.values()))},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()

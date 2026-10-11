"""EXP-20260928-002 事前登録の一部（摂動の実行前に確定する）。

摂動（ISP）の結果は一切読まない。使うのはトークン化済み入力
（遺伝子の検出率・元の位置・入力の長さ・提供者・状態ラベル・細胞種ラベル）と、
analysis/pathway_candidates_002.tsv（文献・データベース由来の候補。本データと無関係）だけである。

段取り（HYP-20260928-01 §3.1 の S1〜S6 をそのまま実装する）:
 (A) 遺伝子の統計表を作る（EXP-20260928-001 と同じ規則・同じデータなので gene_stats_001.csv と一致するはず。検算する）
 (B) 経路セット K=120 を決める（S4 除外 → S5 実行可能性 → S6 切り出し）
 (C) 対照遺伝子の上位集合 980 と、入れ子の部分集合 600（案 B 相当）・520（案 C 相当）を決める
 (D) 副次的な部分セット（非マーカー／目標マーカー条件を除いた版／一般機構／経路特異的）を決める
 (E) 深さをそろえた状態重心の細胞（seed 51 / 52 / 53）
 (F) 実行可能性の記録（細胞数・提供者数・使える対照遺伝子の数）
 (G) config に貼る遺伝子リストと、段階実行の遺伝子リスト

出力: work/{gene_stats_002.csv, pathway_set_002.json, null_pool_002.json, subsets_002.json,
        genes_002.json, state_rows-seed{S}.json, depth_match_summary_002.json,
        precheck_coverage_002.json, config_gene_lists_002.yaml, stage_genes_002.json}
実行: .venv/bin/python experiments/EXP-20260928-002/analysis/select_genes_002.py
"""
from __future__ import annotations

import importlib.util
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260928-002"
SRC3 = ROOT / "experiments/EXP-20260927-003"
SRC4 = ROOT / "experiments/EXP-20260927-004"
SRC1 = ROOT / "experiments/EXP-20260928-001"
TOK = E / "work/tokenized.dataset"        # EXP-20260928-001 のデータ部分を cp -c で複製（sha256 一致を確認済み）
CAND = Path(__file__).resolve().parent / "pathway_candidates_002.tsv"

SEEDS = [51, 52, 53]                      # 既出（0〜5, 11〜13, 21〜23, 31〜33, 41〜43）と重ならない
SUBSET_SEED = 20260929                    # 対照の上位集合の「近い順に取った後の残り」を並べ替える seed（固定）

K_PATH = 120                              # 経路セットの大きさ（PI 指示 2026-09-28）
NULL_SUPERSET, NULL_SUB_B, NULL_SUB_C = 980, 600, 520   # 対照の上位集合と入れ子の部分集合

MAX_CELLS_PER_GENE, MAX_CELLS_PER_DONOR, MIN_CELLS_PER_GENE = 200, 20, 120
MIN_CELLS, MIN_DONORS = 250, 12           # S5(2) 出発点での細胞数・提供者数
MIN_ELIG_NULLS = 20                       # S5(4) 主検定で使える対照遺伝子の下限
PERM_TOL, PERM_COVER = 0.5, 0.8
TOL_N = {"lncells": np.log(1.5), "log2pos": 0.25, "log2len": 0.25}
N_NEAREST_PER_JUDGED = 2                  # 判定遺伝子ごとに最初に確保する近傍の対照遺伝子の数
MAX_SWAP_ROUNDS = 12                      # S5(4) で経路遺伝子を差し替える上限

NULL_ABS_BIAS_MAX = 0.5                   # 対照遺伝子の条件 |bias| < 0.5
MARKER_BIAS_MIN, MARKER_MARGIN_MIN = 0.8, 0.5   # 「目標マーカーの条件」（副次的な部分セットの定義に使う）

MYOF_SUBCLASS = ["MYOF", "cycMYOF"]
STATE_KEY, START, GOAL, ALT = "fib_tls_state", "FIB_neg_TLShigh", "FIB_org", "FIB_neg_TLSlow"

_spec = importlib.util.spec_from_file_location("sg3", Path(__file__).resolve().parent / "select_genes_fib.py")
sg3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sg3)
assert sg3.DEPTH_TOL == 0.10 and sg3.STATE_MAX == 800 and sg3.STATE_PER_DONOR == 60

# ---- S4 除外（HYP-20260928-01 §2.2。結果を見る前に固定）
BLIND = {"TP53", "CDKN1A", "LTBR", "SPP1", "LCN2"}         # 盲検（解除の指示があるまで摂動も解析もしない）
LABEL_GENES = {"CCL19", "CCL21", "CXCL13"}                 # 目標ラベルを定義する遺伝子（循環）
ORIG_TARGETS = {"LTBR", "RELB", "NFKB2"}                   # 元の主対象。セットに入れない
BLOCKED = BLIND | LABEL_GENES | ORIG_TARGETS

HK = {"ACTB": "ENSG00000075624", "GAPDH": "ENSG00000111640", "PPIA": "ENSG00000196262"}
PANEL = {"CERS4": "ENSG00000090661", "TATDN1": "ENSG00000147687", "IFT25": "ENSG00000081870",
         "VAPA": "ENSG00000101558", "PHPT1": "ENSG00000054148", "ANAPC16": "ENSG00000166295",
         "EEF1D": "ENSG00000104529"}
REF_MYOF = {"ACTA2": "ENSG00000107796", "TAGLN": "ENSG00000149591"}

PATHWAY_RULE = (
    "経路セットの選定に、本データの目標ラベル（CCL19/CCL21/CXCL13 の検出）を一切使わない。"
    "S1 経路の構成因子（GO:0038061 と LTbetaR->NIK->IKKalpha->p100/p52->RelB）、"
    "S2 経路の転写標的（一次文献 PMID 12387745 / 22940098 / 23623380 / 28580957）、"
    "S3 キュレーション済みセット（Reactome TNFR2 non-canonical NF-kB / GOBP / NF-kB 標的カタログ）"
    "の順に並べ、S4 で（ラベル定義 3・EXP-20260928-001 の目標マーカー 24・筋線維芽細胞マーカー 24・"
    "盲検 5・元の主対象 3・EXP-20260927-004 で摂動済み 106）を除き、"
    "S5 で実行可能性（辞書にある／出発点 250 細胞・12 提供者以上／抽出で 120 細胞以上／使える対照 20 以上）を当て、"
    "S6 で tier の順・同順位内は遺伝子記号の辞書順に並べて先頭 120 個を取る。無作為抽出も結果参照も使わない。")

NULL_RULE = (
    "対照遺伝子は |bias_goal| < 0.5 かつ |bias_myof| < 0.5、出発点 250 細胞・12 提供者以上、"
    "抽出で 120 細胞以上取れるもの。経路セット・アッセイ対照とは 1 遺伝子も重複させない。"
    "並び順は (i) 判定遺伝子 156 個の順に、元の位置・入力の長さ・細胞数が近いものを 2 個ずつ確保し、"
    "(ii) 残りを seed 20260929 の固定の無作為並べ替えで並べる。"
    f"先頭 {NULL_SUB_C} 個を案 C 相当、先頭 {NULL_SUB_B} 個を案 B 相当、"
    f"全 {NULL_SUPERSET} 個を主判定（案 A）とする。部分集合は入れ子である。")


# ------------------------------------------------------------------ (A) 遺伝子統計
def base_states(ds) -> np.ndarray:
    """状態ラベルを「付け替える前の形」に戻して返す（MYOF_ 接頭辞を外す）。"""
    st = np.array([str(x) for x in ds[STATE_KEY]])
    return np.array([x[len("MYOF_"):] if x.startswith("MYOF_") else x for x in st])


def gene_stats(ds) -> tuple[pd.DataFrame, dict]:
    import scipy.sparse as sp
    st = base_states(ds)
    sc = np.array([str(x) for x in ds["subclass_l2"]])
    ids = [list(x) for x in ds["input_ids"]]
    indptr = np.cumsum([0] + [len(x) for x in ids])
    col = np.concatenate([np.asarray(x) for x in ids])
    ntok = int(col.max()) + 1
    M = sp.csr_matrix((np.ones(len(col), bool), col, indptr), shape=(len(ids), ntok)).tocsc()
    myo = np.isin(sc, MYOF_SUBCLASS)
    sets = {"start": (st == START) & ~myo, "goal": (st == GOAL) & ~myo, "alt": (st == ALT) & ~myo,
            "myof": myo & np.isin(st, [START, ALT])}
    donors = np.array([str(x) for x in ds["donor_id"]])
    d = pd.DataFrame({"token": np.arange(ntok)})
    for k, v in sets.items():
        d[f"det_{k}"] = np.asarray(M[v].sum(0)).ravel() / v.sum()
        d[f"n_{k}"] = np.asarray(M[v].sum(0)).ravel()
    eps = 1.0 / int(sets["goal"].sum())
    d["bias_goal"] = np.log2((d.det_goal + eps) / (d.det_start + eps))
    d["bias_myof"] = np.log2((d.det_myof + eps) / (d.det_start + eps))
    Ms = M[sets["start"]].tocsc()
    codes = pd.Categorical(donors[sets["start"]]).codes
    nd = np.zeros(ntok, int)
    for t in range(ntok):
        r = Ms[:, t].nonzero()[0]
        nd[t] = len(np.unique(codes[r])) if len(r) else 0
    d["n_donor_start"] = nd
    ref = pd.read_csv(SRC3 / "work/gene_stats_fib.csv")[["token", "symbol", "ensembl"]]
    d = d.merge(ref, on="token", how="inner")
    d = d[d.symbol.astype(str) != ""].reset_index(drop=True)
    info = {k: {"n_cells": int(v.sum()), "n_donors": int(len(set(donors[v])))} for k, v in sets.items()}
    info["eps"] = eps
    return d, info


# ------------------------------------------------------------------ (E) 抽出と被覆
class Profiles:
    """出発点（MYOF を除いた FIB_neg_TLShigh）から、遺伝子 x seed ごとに抽出される細胞と
    (log2(1+元の位置), log2(入力の長さ)) を返す。抽出規則は run_isp.balanced_subset と同じ。"""

    def __init__(self, ds, d: pd.DataFrame) -> None:
        import scipy.sparse as sp
        st = base_states(ds)
        sc = np.array([str(x) for x in ds["subclass_l2"]])
        donors_all = [str(x) for x in ds["donor_id"]]
        ids = [list(x) for x in ds["input_ids"]]
        rows = np.flatnonzero((st == START) & ~np.isin(sc, MYOF_SUBCLASS))
        self.rows = rows
        self.sel = [ids[i] for i in rows]
        self.pos = [{t: j for j, t in enumerate(x)} for x in self.sel]
        self.donors = [donors_all[i] for i in rows]
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

    def get(self, sym: str, seed: int):
        key = (sym, seed)
        if key not in self._cells:
            t = int(self.tok[sym])
            idx = self.m[:, t].nonzero()[0].tolist() if t < self.ntok else []
            keep = sg3.balanced(idx, self.donors, MAX_CELLS_PER_GENE, MAX_CELLS_PER_DONOR,
                                random.Random(f"{seed}-{self.eid[sym]}"))
            lp = np.log2(1 + np.array([self.pos[i][t] for i in keep], float))
            ll = np.log2(self.lens[keep]) if len(keep) else np.zeros(0)
            self._cells[key] = (keep, lp, ll)
        return self._cells[key]

    def n_donors(self, sym: str, seed: int) -> int:
        keep, _, _ = self.get(sym, seed)
        return len({self.donors[i] for i in keep})

    def covers(self, gene: str, null: str, seed: int) -> bool:
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


# ------------------------------------------------------------------ main
def main() -> None:
    from datasets import load_from_disk
    (E / "work").mkdir(parents=True, exist_ok=True)
    ds = load_from_disk(str(TOK))
    d, info = gene_stats(ds)
    d.to_csv(E / "work/gene_stats_002.csv", index=False)
    print("状態ごとの細胞数（MYOF を除く）:", json.dumps(info, ensure_ascii=False))
    # 同じデータ・同じ規則なので EXP-20260928-001 の表と一致するはず（検算）
    ref1 = pd.read_csv(SRC1 / "work/gene_stats_001.csv")
    cmp_cols = ["token", "n_start", "n_goal", "n_alt", "n_myof", "bias_goal", "bias_myof", "n_donor_start"]
    a = ref1.sort_values("token").reset_index(drop=True)
    b = d.sort_values("token").reset_index(drop=True)
    assert len(a) == len(b) and (a.ensembl.to_numpy() == b.ensembl.to_numpy()).all(), \
        "gene_stats の遺伝子が EXP-20260928-001 と一致しない"
    dmax = float(np.abs(a[cmp_cols].to_numpy() - b[cmp_cols].to_numpy()).max())
    assert dmax < 1e-9, f"gene_stats の数値が EXP-20260928-001 と一致しない（最大差 {dmax}）"
    print(f"gene_stats 検算: EXP-20260928-001 と {len(d)} 行・{len(cmp_cols)} 列が一致（最大差 {dmax:.2e}。"
          "差は CSV への書き出しでの丸めによるもの）")

    ref = d.set_index("symbol")
    have = set(d.symbol)

    # ---------------------------------------------------------- S4 除外集合
    g1 = json.loads((SRC1 / "work/genes_001.json").read_text())["genes"]
    markers_goal_001 = {g["symbol"] for g in g1 if g["role"] == "marker_goal"}
    markers_myof_001 = {g["symbol"] for g in g1 if g["role"] == "marker_myof"}
    measured_004 = {g["symbol"] for g in json.loads((SRC4 / "work/genes_004.json").read_text())["genes"]}
    exclude_s4 = (LABEL_GENES | markers_goal_001 | markers_myof_001 | BLIND | ORIG_TARGETS
                  | measured_004 | set(HK) | set(PANEL) | set(REF_MYOF))
    print(f"\nS4 除外: ラベル定義 {len(LABEL_GENES)} / 目標マーカー {len(markers_goal_001)} / "
          f"筋線維芽細胞マーカー {len(markers_myof_001)} / 盲検 {len(BLIND)} / 元の主対象 {len(ORIG_TARGETS)} / "
          f"EXP-20260927-004 で摂動済み {len(measured_004)} → 合計 {len(exclude_s4)}（重複を除く）")

    # ---------------------------------------------------------- S1〜S3 候補の読み込み
    cand = pd.read_csv(CAND, sep="\t", comment="#")
    assert list(cand.columns) == ["symbol", "tier", "group", "note"]
    assert cand.symbol.is_unique, "候補に重複がある"
    print(f"候補（S1〜S3, 重複除去後）{len(cand)}")

    trace = {"n_candidates": int(len(cand))}
    cand["in_dict"] = cand.symbol.isin(have)
    cand["s4_excluded"] = cand.symbol.isin(exclude_s4)
    trace["s4_excluded"] = sorted(cand.symbol[cand.s4_excluded].tolist())
    trace["not_in_data"] = sorted(cand.symbol[~cand.in_dict].tolist())
    c = cand[cand.in_dict & ~cand.s4_excluded].copy()
    c = c.join(ref[["n_start", "n_donor_start", "bias_goal", "bias_myof", "ensembl"]], on="symbol")
    print(f"S4 除外 {int(cand.s4_excluded.sum())} / 本データに無い {int((~cand.in_dict).sum())} "
          f"→ 残り {len(c)}")

    # ---------------------------------------------------------- S5 (2) 細胞数・提供者数
    ok2 = (c.n_start >= MIN_CELLS) & (c.n_donor_start >= MIN_DONORS)
    trace["s5_reject_cells_donors"] = sorted(
        f"{r.symbol}({int(r.n_start)}細胞/{int(r.n_donor_start)}提供者)" for _, r in c[~ok2].iterrows())
    c = c[ok2].copy()
    print(f"S5(2) 出発点 {MIN_CELLS} 細胞・{MIN_DONORS} 提供者以上 → 残り {len(c)}")

    prof = Profiles(ds, d)

    # ---------------------------------------------------------- S5 (3) 抽出できる細胞数
    keep3 = []
    rej3 = []
    for sym in c.symbol:
        ns = [len(prof.get(sym, s)[0]) for s in SEEDS]
        (keep3 if min(ns) >= MIN_CELLS_PER_GENE else rej3).append((sym, min(ns)))
    trace["s5_reject_selected_cells"] = sorted(f"{s}({n})" for s, n in rej3)
    c = c[c.symbol.isin({s for s, _ in keep3})].copy()
    print(f"S5(3) 抽出で {MIN_CELLS_PER_GENE} 細胞以上（3 seed すべて） → 残り {len(c)}")

    # ---------------------------------------------------------- S6 並べ替え
    c["tier_rank"] = c.tier.map({"S1": 0, "S2": 1, "S3": 2})
    c = c.sort_values(["tier_rank", "symbol"]).reset_index(drop=True)
    assert len(c) >= K_PATH, f"S5 通過が {len(c)} 個で K={K_PATH} に届かない"
    print(f"S6 並べ替え後の通過候補 {len(c)}（S1 {int((c.tier=='S1').sum())} / "
          f"S2 {int((c.tier=='S2').sum())} / S3 {int((c.tier=='S3').sum())}）")

    # ---------------------------------------------------------- 対照遺伝子の母集団
    d_named = d[d.symbol.notna() & (d.symbol.astype(str) != "nan")]
    null_base = d_named[(d_named.bias_goal.abs() < NULL_ABS_BIAS_MAX) & (d_named.bias_myof.abs() < NULL_ABS_BIAS_MAX)
                  & (d_named.n_start >= MIN_CELLS) & (d_named.n_donor_start >= MIN_DONORS)
                  & ~d_named.symbol.isin(exclude_s4 | set(cand.symbol))
                  & ~d_named.symbol.isin(set(sg3.NULL_EXCLUDE))
                  & ~d_named.symbol.str.startswith(sg3.NULL_EXCLUDE_PREFIX)
                  & ~d_named.symbol.str.startswith(sg3.EXCLUDE_PREFIX)
                  & ~d_named.symbol.isin(sg3.EXCLUDE_SYMBOLS)].copy()
    print(f"\n対照遺伝子の母集団（条件を満たす遺伝子）{len(null_base)}")

    def profile_med(sym: str) -> tuple[float, float, int]:
        keep, lp, ll = prof.get(sym, SEEDS[0])
        return (float(np.median(lp)) if len(lp) else np.nan,
                float(np.median(ll)) if len(ll) else np.nan, len(keep))

    bm = {s: profile_med(s) for s in null_base.symbol}
    null_base["l2pos"] = [bm[s][0] for s in null_base.symbol]
    null_base["l2len"] = [bm[s][1] for s in null_base.symbol]
    null_base["n_sel"] = [bm[s][2] for s in null_base.symbol]
    null_base = null_base[(null_base.n_sel >= MIN_CELLS_PER_GENE) & null_base.l2pos.notna()].reset_index(drop=True)
    print(f"うち抽出で {MIN_CELLS_PER_GENE} 細胞以上取れるもの {len(null_base)}")
    assert len(null_base) >= NULL_SUPERSET, f"対照遺伝子が {len(null_base)} 個しかない"

    assay = list(HK) + list(PANEL) + list(REF_MYOF)   # アッセイ対照 36（うち筋線維芽細胞マーカー 24 は下で足す）
    assay = list(markers_myof_001) + assay
    assert len(assay) == 36, len(assay)

    # ---------------------------------------------------------- 対照の上位集合と入れ子の部分集合
    def build_nulls(judged: list[str]) -> list[str]:
        meds = {s: profile_med(s) for s in judged}
        got: set[str] = set()
        ordered: list[str] = []
        for g in judged:                      # (i) 判定遺伝子ごとに近い順に 2 個ずつ
            x = null_base[~null_base.symbol.isin(got)].copy()
            x["dist"] = pd.concat([
                np.log(x.n_sel / max(meds[g][2], 1)).abs() / TOL_N["lncells"],
                (x.l2pos - meds[g][0]).abs() / TOL_N["log2pos"],
                (x.l2len - meds[g][1]).abs() / TOL_N["log2len"]], axis=1).max(axis=1)
            for sym in x.sort_values(["dist", "symbol"]).symbol.head(N_NEAREST_PER_JUDGED):
                ordered.append(sym)
                got.add(sym)
        rest = sorted(set(null_base.symbol) - got)          # (ii) 残りを固定の無作為並べ替え
        random.Random(SUBSET_SEED).shuffle(rest)
        ordered += rest
        return ordered[:NULL_SUPERSET]

    # ---------------------------------------------------------- S5(4) 使える対照遺伝子 >= 20
    path_syms = list(c.symbol[:K_PATH])
    next_i = K_PATH
    swaps: list[dict] = []
    for rnd in range(MAX_SWAP_ROUNDS):
        judged = path_syms + assay
        nulls = build_nulls(judged)
        sub_c, sub_b = nulls[:NULL_SUB_C], nulls[:NULL_SUB_B]
        short = {}
        for g in judged:
            v = min(prof.n_eligible(g, sub_c, s) for s in SEEDS)   # 最も厳しい版（520）で判定する
            if v < MIN_ELIG_NULLS:
                short[g] = v
        print(f"S5(4) {rnd} 巡目: 使える対照が {MIN_ELIG_NULLS} 未満の判定遺伝子 {len(short)} 個"
              + (" " + ", ".join(f"{k}({v})" for k, v in sorted(short.items(), key=lambda x: x[1])[:8]) if short else ""))
        bad_path = [g for g in short if g in path_syms]
        if not bad_path:
            if short:
                raise SystemExit(f"アッセイ対照で使える対照が不足: {short}")
            break
        for g in bad_path:
            if next_i >= len(c):
                raise SystemExit("差し替える候補が尽きた")
            newg = c.symbol.iloc[next_i]
            swaps.append({"removed": g, "n_eligible": short[g], "added": newg, "round": rnd})
            path_syms[path_syms.index(g)] = newg
            next_i += 1
    else:
        raise SystemExit("S5(4) が規定回数で収束しなかった")

    path_syms = sorted(path_syms, key=lambda s: (int(c[c.symbol == s].tier_rank.iloc[0]), s))
    judged = path_syms + assay
    nulls = build_nulls(judged)
    sub = {"A_980": nulls[:NULL_SUPERSET], "B_600": nulls[:NULL_SUB_B], "C_520": nulls[:NULL_SUB_C]}
    assert set(sub["C_520"]) < set(sub["B_600"]) < set(sub["A_980"]), "対照の部分集合が入れ子になっていない"
    assert not (set(nulls) & set(judged)), "対照と判定遺伝子が重複している"
    print(f"\n対照の上位集合 {len(sub['A_980'])}（部分集合 {len(sub['B_600'])} / {len(sub['C_520'])}、入れ子を確認）")

    elig = {g: {s: prof.n_eligible(g, sub["C_520"], s) for s in SEEDS} for g in judged}
    worst = min(min(v.values()) for v in elig.values())
    print(f"使える対照遺伝子の最小（520 版）{worst}（下限 {MIN_ELIG_NULLS}）")

    # ---------------------------------------------------------- (D) 副次的な部分セット
    p = c.set_index("symbol").loc[path_syms]
    nonmarker = sorted(p.index[(p.bias_goal.abs() < 0.5) & (p.bias_myof.abs() < 0.5)])
    goalmarker = sorted(p.index[(p.bias_goal >= MARKER_BIAS_MIN)
                                & ((p.bias_goal - p.bias_myof) >= MARKER_MARGIN_MIN)])
    subsets = {
        "pathway_all": path_syms,
        "nonmarker": nonmarker,                                     # HYP §4.2 の副次的な予測
        "drop_goalmarker": [s for s in path_syms if s not in set(goalmarker)],
        "general_mechanism": [s for s in path_syms if p.loc[s, "group"] == "一般機構"],
        "pathway_specific": [s for s in path_syms if p.loc[s, "group"] == "経路特異的"],
        "tier_S1": [s for s in path_syms if p.loc[s, "tier"] == "S1"],
        "tier_S2": [s for s in path_syms if p.loc[s, "tier"] == "S2"],
        "tier_S3": [s for s in path_syms if p.loc[s, "tier"] == "S3"],
    }
    print("\n経路セットの部分セット: " + ", ".join(f"{k}={len(v)}" for k, v in subsets.items()))
    print(f"  目標マーカーの条件を満たす経路遺伝子 {len(goalmarker)}: {goalmarker}")

    # ---------------------------------------------------------- (E) 状態重心の細胞
    st = base_states(ds)
    sc = np.array([str(x) for x in ds["subclass_l2"]])
    donors = [str(x) for x in ds["donor_id"]]
    length = np.array(ds["length"], dtype=float)
    log_len = np.log(length)
    myo = np.isin(sc, MYOF_SUBCLASS)
    summary = {}
    for s in SEEDS:
        rng2 = random.Random(s)
        cap = {}
        for stt in (START, GOAL, ALT):
            idx = np.flatnonzero((st == stt) & ~myo).tolist()
            cap[stt] = sg3.balanced(idx, donors, sg3.STATE_MAX, sg3.STATE_PER_DONOR, rng2)
        rows: dict[str, list[int]] = {}
        info_s: dict[str, dict] = {}
        for stt in (GOAL, ALT):
            m0 = np.zeros(len(st), bool); m0[cap[stt]] = True
            m1 = np.zeros(len(st), bool); m1[cap[START]] = True
            m = sg3.depth_match(log_len, m0, m1, sg3.DEPTH_TOL)
            rows[stt] = sorted(np.flatnonzero(m & m0).tolist())
            rows["_start_" + stt] = sorted(np.flatnonzero(m & m1).tolist())
            info_s[stt] = {"n_matched": len(rows[stt]),
                           "median_length": {stt: float(np.median(length[rows[stt]])),
                                             START: float(np.median(length[rows["_start_" + stt]]))}}
        rows[START] = rows.pop("_start_" + GOAL)
        rows.pop("_start_" + ALT)
        info_s[START] = {"n_matched": len(rows[START]),
                         "median_length": float(np.median(length[rows[START]]))}
        # 検算: 行番号がその状態を指し、除いたはずの細胞種が 0 件であること（guardrails 2026-09-28）
        for stt, rr in rows.items():
            assert (st[rr] == stt).all(), f"seed {s}: {stt} の行番号が状態と一致しない"
            assert int(myo[rr].sum()) == 0, f"seed {s}: {stt} に筋線維芽細胞が混ざっている"
        (E / f"work/state_rows-seed{s}.json").write_text(
            json.dumps({"rows": rows}, ensure_ascii=False) + "\n")
        summary[s] = info_s
    (E / "work/depth_match_summary_002.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=1) + "\n")
    print("\n深さマッチ後の細胞数:", {s: {k: v["n_matched"] for k, v in summary[s].items()} for s in SEEDS})

    # ---------------------------------------------------------- (G) 出力
    ens = dict(zip(d.symbol, d.ensembl))
    genes = []
    for sym in path_syms:
        genes.append({"symbol": sym, "ensembl": ens[sym], "role": "pathway",
                      "tier": p.loc[sym, "tier"], "group": p.loc[sym, "group"],
                      "note": p.loc[sym, "note"],
                      "bias_goal": round(float(ref.loc[sym, "bias_goal"]), 3),
                      "bias_myof": round(float(ref.loc[sym, "bias_myof"]), 3),
                      "n_start": int(ref.loc[sym, "n_start"]),
                      "n_donor_start": int(ref.loc[sym, "n_donor_start"])})
    for sym in sorted(markers_myof_001):
        genes.append({"symbol": sym, "ensembl": ens[sym], "role": "marker_myof",
                      "bias_goal": round(float(ref.loc[sym, "bias_goal"]), 3),
                      "bias_myof": round(float(ref.loc[sym, "bias_myof"]), 3)})
    for sym, e in HK.items():
        genes.append({"symbol": sym, "ensembl": e, "role": "housekeeping"})
    for sym, e in PANEL.items():
        genes.append({"symbol": sym, "ensembl": e, "role": "negctrl_panel"})
    for sym, e in REF_MYOF.items():
        genes.append({"symbol": sym, "ensembl": e, "role": "reference_myof"})
    for rank, sym in enumerate(sub["A_980"]):
        genes.append({"symbol": sym, "ensembl": ens[sym], "role": "negative", "null_rank": rank,
                      "bias_goal": round(float(ref.loc[sym, "bias_goal"]), 3),
                      "bias_myof": round(float(ref.loc[sym, "bias_myof"]), 3)})
    assert len({g["symbol"] for g in genes}) == len(genes), "遺伝子が重複している"
    assert not ({g["symbol"] for g in genes} & BLOCKED), "回さない遺伝子が混入している"
    assert len(genes) == K_PATH + 36 + NULL_SUPERSET, len(genes)

    (E / "work/genes_002.json").write_text(json.dumps(
        {"rule_pathway": PATHWAY_RULE, "rule_null": NULL_RULE, "n": len(genes),
         "seeds": SEEDS, "genes": genes}, ensure_ascii=False, indent=1) + "\n")
    (E / "work/pathway_set_002.json").write_text(json.dumps(
        {"rule": PATHWAY_RULE, "K": K_PATH, "swaps": swaps, "trace": trace,
         "genes": [g for g in genes if g["role"] == "pathway"]}, ensure_ascii=False, indent=1) + "\n")
    (E / "work/null_pool_002.json").write_text(json.dumps(
        {"rule": NULL_RULE, "subset_seed": SUBSET_SEED,
         "n": {k: len(v) for k, v in sub.items()},
         "versions": {k: [{"symbol": s, "ensembl": ens[s]} for s in v] for k, v in sub.items()}},
        ensure_ascii=False, indent=1) + "\n")
    (E / "work/subsets_002.json").write_text(json.dumps(
        {"note": "経路セットの部分セット（結果を見る前に固定）。主判定は pathway_all。",
         "goal_marker_condition_genes": goalmarker,
         "subsets": {k: [{"symbol": s, "ensembl": ens[s]} for s in v] for k, v in subsets.items()}},
        ensure_ascii=False, indent=1) + "\n")

    # 段階実行の遺伝子リスト（PI 指示: 途中で止まっても副次的な版が評価できること）
    stage = {
        "stage1": {"note": f"経路 {K_PATH} + アッセイ対照 36 + 対照 {NULL_SUB_C} → 案 C 相当（520 版）が評価できる",
                   "negatives": sub["C_520"]},
        "stage2": {"note": f"対照を {NULL_SUB_B} まで足す → 案 B 相当（600 版）が評価できる",
                   "negatives": sub["B_600"]},
        "stage3": {"note": f"対照を {NULL_SUPERSET} まで足す → 主判定（案 A・980 版）",
                   "negatives": sub["A_980"]},
    }
    for k, v in stage.items():
        v["n_negatives"] = len(v["negatives"])
        v["n_genes_cumulative"] = K_PATH + 36 + len(v["negatives"])
        v["cell_gene_cumulative"] = v["n_genes_cumulative"] * MAX_CELLS_PER_GENE * len(SEEDS)
    (E / "work/stage_genes_002.json").write_text(json.dumps(stage, ensure_ascii=False, indent=1) + "\n")

    (E / "work/precheck_coverage_002.json").write_text(json.dumps(
        {"n_cells": {f"{g}@{s}": len(prof.get(g, s)[0]) for g in judged for s in SEEDS},
         "n_donors": {f"{g}@{s}": prof.n_donors(g, s) for g in judged for s in SEEDS},
         "n_eligible_null_genes_520": {f"{g}@{s}": elig[g][s] for g in judged for s in SEEDS},
         "min_eligible_among_judged_520": worst,
         "n_genes_total": len(genes), "n_nulls": len(sub["A_980"]),
         "state_cells": info,
         "total_cell_gene": len(genes) * MAX_CELLS_PER_GENE * len(SEEDS)},
        ensure_ascii=False, indent=1) + "\n")

    lines = ["# 自動生成: analysis/select_genes_002.py。config.yaml へそのまま貼る"]
    for role, key in (("pathway", "genes_to_perturb"), ("marker_myof", "controls.positive"),
                      ("negative", "controls.negative_sampler.genes")):
        lines.append(f"# --- {key} ({role})")
        for g in genes:
            if g["role"] == role:
                lines.append(f"    - {g['ensembl']}  # {g['symbol']}")
    (E / "work/config_gene_lists_002.yaml").write_text("\n".join(lines) + "\n")

    print(f"\n合計 {len(genes)} 遺伝子（経路 {K_PATH} / 筋線維芽細胞マーカー 24 / 恒常発現 3 / "
          f"中立パネル 7 / 参照 2 / 対照 {len(sub['A_980'])}）")
    print(f"延べ 細胞 x 遺伝子 = {len(genes) * MAX_CELLS_PER_GENE * len(SEEDS):,}")
    print("段階: " + ", ".join(f"{k} 累計 {v['n_genes_cumulative']} 遺伝子 / "
                               f"{v['cell_gene_cumulative']:,} 細胞x遺伝子" for k, v in stage.items()))


if __name__ == "__main__":
    main()

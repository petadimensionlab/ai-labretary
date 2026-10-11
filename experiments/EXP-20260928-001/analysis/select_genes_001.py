"""EXP-20260928-001 事前登録の一部（摂動の実行前に確定する）。

摂動（ISP）の結果は一切読まない。使うのはトークン化済み入力
（遺伝子の検出率・元の位置・入力の長さ・提供者・状態ラベル・細胞種ラベル）だけである。

(A) 遺伝子の統計表を作る（検出率で状態への偏りを測る。過剰発現型が動かすのは「入力に載っているか」なので、
    発現量ではなく検出率で偏りを定義する）
(B) 目標マーカーセット（K=24）と 筋線維芽細胞マーカーセット（K=24）を固定する
(C) 対照遺伝子（null）を選び、主検定で使える対照遺伝子が全判定遺伝子で 20 以上になるまで補充する
(D) 深さをそろえた状態重心の細胞（seed 41 / 42 / 43）。**MYOF・cycMYOF を 3 状態から除く**
    （EXP-20260927-004 validation.md §2.7-4: 出発点と筋線維芽細胞の重心が細胞を 9.8〜11.2% 共有していた）
(E) 実行可能性（遺伝子ごとの細胞数・提供者数、使える対照遺伝子の数）

出力: work/{gene_stats_001.csv, marker_sets_001.csv, null_pool_001.json, genes_001.json,
        state_rows-seed{S}.json, depth_match_summary_001.json, precheck_coverage_001.json,
        config_gene_lists_001.yaml}
実行: .venv/bin/python experiments/EXP-20260928-001/analysis/select_genes_001.py
"""
from __future__ import annotations

import importlib.util
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260928-001"
SRC3 = ROOT / "experiments/EXP-20260927-003"
SRC4 = ROOT / "experiments/EXP-20260927-004"
TOK = E / "work/tokenized.dataset"        # EXP-20260927-004 のデータ部分を cp -c で複製（cache は削除済み）

SEEDS = [41, 42, 43]                      # 既出（0〜5, 11〜13, 21〜23, 31〜33）と重ならない
SEL_SEED = 20260928                       # マーカーセットの無作為半分を選ぶときの seed（固定）

MAX_CELLS_PER_GENE, MAX_CELLS_PER_DONOR, MIN_CELLS_PER_GENE = 200, 20, 120
MIN_CELLS, MIN_DONORS = 250, 12
MIN_ELIG_NULLS = 20
PERM_TOL, PERM_COVER = 0.5, 0.8
TOL_N = {"lncells": np.log(1.5), "log2pos": 0.25, "log2len": 0.25}
N_AUG_SEARCH, N_AUG_ROUNDS = 400, 8
MAX_NULLS = 110

# ---- マーカーセットの選定規則（結果を見る前に固定）
MARKER_RULE = (
    "検出率の偏り bias_X = log2((det_X + eps)/(det_start + eps)) を使う。"
    "eps = 1 / 目標状態の細胞数。状態 X のマーカーの条件は "
    "(1) bias_X >= 0.8、(2) もう一方のマーカー状態 Y に対して bias_X - bias_Y >= 0.5、"
    "(3) 出発点で 250 細胞以上・提供者 12 人以上、(4) 抽出規則で 120 細胞以上取れる、"
    "(5) 対照遺伝子の補充後に主検定で使える対照遺伝子が 20 以上。"
    "条件を満たす遺伝子を bias_X の大きい順に並べ、上位 12 個（偏りの大きい側）と、"
    "残りから seed 20260928 で無作為に選んだ 12 個を合わせて 24 個とする。"
    "上位半分と無作為半分を分けるのは、結論が『偏りが極端な遺伝子だけの性質』でないことを"
    "同じ実験の中で確かめるためである（両半分の結果は記述値として別々に報告する）。")
MARKER_BIAS_MIN, MARKER_MARGIN_MIN, K_SET, K_TOP = 0.8, 0.5, 24, 12
NULL_ABS_BIAS_MAX = 0.5

MYOF_SUBCLASS = ["MYOF", "cycMYOF"]
STATE_KEY, START, GOAL, ALT = "fib_tls_state", "FIB_neg_TLShigh", "FIB_org", "FIB_neg_TLSlow"

_spec = importlib.util.spec_from_file_location("sg3", SRC3 / "analysis/select_genes_fib.py")
sg3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sg3)
assert sg3.DEPTH_TOL == 0.10 and sg3.STATE_MAX == 800 and sg3.STATE_PER_DONOR == 60

BLOCKED = {"LTBR", "RELB", "NFKB2",          # 主対象（HYP-20260925-02）。本実験では回さない
           "CDKN1A",                          # 参考対象・盲検
           "TP53", "SPP1", "LCN2",            # 盲検（未解除）
           "CCL19", "CCL21", "CXCL13"}        # 目標ラベルを定義する遺伝子（循環性）

HK = {"ACTB": "ENSG00000075624", "GAPDH": "ENSG00000111640", "PPIA": "ENSG00000196262"}
PANEL = {"CERS4": "ENSG00000090661", "TATDN1": "ENSG00000147687", "IFT25": "ENSG00000081870",
         "VAPA": "ENSG00000101558", "PHPT1": "ENSG00000054148", "ANAPC16": "ENSG00000166295",
         "EEF1D": "ENSG00000104529"}
REF_MYOF = {"ACTA2": "ENSG00000107796", "TAGLN": "ENSG00000149591"}   # EXP-004 と比べるための参照（既測定）


# ------------------------------------------------------------------ (A) 遺伝子統計
def base_states(ds) -> np.ndarray:
    """状態ラベルを「付け替える前の形」に戻して返す。
    work/tokenized.dataset は MYOF / cycMYOF の細胞の状態を `MYOF_<元の状態>` に付け替えてあるが
    （plan.md §3.1）、本スクリプトは細胞種の絞り込みを自前で行うので、どちらの形でも同じ結果になるように
    接頭辞を外して扱う。**付け替えの前に走らせても後に走らせても出力が変わらないことを保証するための処理。**"""
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
    # 出発点での提供者数
    Ms = M[sets["start"]].tocsc()
    ds_don = donors[sets["start"]]
    codes = pd.Categorical(ds_don).codes
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
    """出発点（MYOF を除いた FIB_neg_TLShigh）から、遺伝子 × seed ごとに抽出される細胞と
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
            lp = np.log2(1 + np.array([self.sel[i].index(t) for i in keep], float))
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


# ------------------------------------------------------------------ (B) マーカーセット
def excluded_mask(d: pd.DataFrame, extra: set[str]) -> pd.Series:
    return (d.symbol.isin(sg3.EXCLUDE_SYMBOLS | sg3.BLIND | sg3.LABEL_GENES | BLOCKED
                          | set(HK) | set(PANEL) | set(REF_MYOF) | extra)
            | d.symbol.str.startswith(sg3.EXCLUDE_PREFIX)
            | d.symbol.isin({k for k, v in sg3.MANUAL_REVIEW.items() if v[0] == "exclude"}))


def pick_markers(d: pd.DataFrame, axis: str, other: str, measured: set[str], taken: set[str]) -> pd.DataFrame:
    b, ob = f"bias_{axis}", f"bias_{other}"
    cand = d[(d[b] >= MARKER_BIAS_MIN) & ((d[b] - d[ob]) >= MARKER_MARGIN_MIN)
             & (d.n_start >= MIN_CELLS) & (d.n_donor_start >= MIN_DONORS)
             & ~excluded_mask(d, measured | taken)].copy()
    return cand.sort_values([b, "symbol"], ascending=[False, True]).reset_index(drop=True)


def main() -> None:
    from datasets import load_from_disk
    (E / "work").mkdir(parents=True, exist_ok=True)
    ds = load_from_disk(str(TOK))
    d, info = gene_stats(ds)
    d.to_csv(E / "work/gene_stats_001.csv", index=False)
    print("状態ごとの細胞数（MYOF を除く）:", json.dumps(info, ensure_ascii=False))

    measured = {g["symbol"] for g in json.loads((SRC4 / "work/genes_004.json").read_text())["genes"]}
    print(f"EXP-20260927-004 で既に摂動した遺伝子 {len(measured)} 個をマーカーセットから除く")

    prof = Profiles(ds, d)
    ref = d.set_index("symbol")

    # --- 対照遺伝子の母集団（どちらの状態にも偏らない遺伝子）
    null_base = d[(d.bias_goal.abs() < NULL_ABS_BIAS_MAX) & (d.bias_myof.abs() < NULL_ABS_BIAS_MAX)
                  & (d.n_start >= MIN_CELLS) & (d.n_donor_start >= MIN_DONORS)
                  & ~excluded_mask(d, set())
                  & ~d.symbol.isin(set(sg3.NULL_EXCLUDE))
                  & ~d.symbol.str.startswith(sg3.NULL_EXCLUDE_PREFIX)].copy()
    null_base["l2pos"] = np.log2(1 + null_base.det_start * 0)  # 下で seed ごとに実測を使うので placeholder
    print(f"対照遺伝子の母集団 {len(null_base)}")

    # 元の位置・入力の長さは seed 41 の実測の中央値を使う（抽出規則と同じ細胞から）
    def profile_med(sym: str) -> tuple[float, float, int]:
        keep, lp, ll = prof.get(sym, SEEDS[0])
        return (float(np.median(lp)) if len(lp) else np.nan,
                float(np.median(ll)) if len(ll) else np.nan, len(keep))

    # --- マーカーセットを選ぶ（実行可能性を確かめながら）
    rng = random.Random(SEL_SEED)
    chosen: dict[str, list[tuple[str, str]]] = {}
    rejected: list[dict] = []
    taken: set[str] = set()
    for axis, other in (("goal", "myof"), ("myof", "goal")):
        cand = pick_markers(d, axis, other, measured, taken)
        print(f"\n{axis} マーカー候補 {len(cand)}（規則: bias >= {MARKER_BIAS_MIN}, "
              f"差 >= {MARKER_MARGIN_MIN}, 出発点 {MIN_CELLS} 細胞・{MIN_DONORS} 提供者以上）")
        ok: list[tuple[str, str]] = []
        for half, pool in (("top", list(cand.symbol[:60])), ("random", None)):
            if half == "random":
                rest = [s for s in cand.symbol if s not in {x for x, _ in ok}]
                rng.shuffle(rest)
                pool = rest
            need = K_TOP
            for sym in pool:
                if need == 0:
                    break
                if sym in {x for x, _ in ok} or sym in taken:
                    continue
                ns = [len(prof.get(sym, s)[0]) for s in SEEDS]
                if min(ns) < MIN_CELLS_PER_GENE:
                    rejected.append({"symbol": sym, "axis": axis,
                                     "reason": f"抽出できる細胞が {min(ns)} で下限 {MIN_CELLS_PER_GENE} 未満"})
                    continue
                ok.append((sym, half))
                taken.add(sym)
                need -= 1
            if need:
                raise SystemExit(f"{axis} の {half} 半分が {K_TOP - need} 個しか選べない")
        chosen[axis] = ok
        print(f"  {axis}: {[s for s, _ in ok]}")

    judged = [s for v in chosen.values() for s, _ in v]

    # --- 対照遺伝子（null）: 判定遺伝子ごとに 元の位置・入力の長さ・細胞数の近いものから選ぶ
    meds = {s: profile_med(s) for s in list(judged) + list(HK) + list(PANEL) + list(REF_MYOF)}
    base = null_base[~null_base.symbol.isin(taken)].copy()
    bmeds = {}
    for s in base.symbol:
        bmeds[s] = profile_med(s)
    base["l2pos"] = [bmeds[s][0] for s in base.symbol]
    base["l2len"] = [bmeds[s][1] for s in base.symbol]
    base["n_sel"] = [bmeds[s][2] for s in base.symbol]
    base = base[(base.n_sel >= MIN_CELLS_PER_GENE) & base.l2pos.notna()].reset_index(drop=True)
    print(f"\n対照遺伝子の母集団（抽出可能なもの）{len(base)}")

    nulls: list[str] = []
    got = set(taken)
    for g in judged:                      # 判定遺伝子ごとに近い順に 2 個ずつ確保
        x = base[~base.symbol.isin(got)].copy()
        x["dist"] = pd.concat([
            np.log(x.n_sel / max(meds[g][2], 1)).abs() / TOL_N["lncells"],
            (x.l2pos - meds[g][0]).abs() / TOL_N["log2pos"],
            (x.l2len - meds[g][1]).abs() / TOL_N["log2len"]], axis=1).max(1)
        for sym in x.sort_values(["dist", "symbol"]).symbol.head(2):
            nulls.append(sym)
            got.add(sym)

    # 補充: 主検定で使える対照遺伝子が 20 未満の判定遺伝子がなくなるまで
    all_judged = judged + list(HK) + list(PANEL) + list(REF_MYOF)
    for rnd in range(N_AUG_ROUNDS):
        short = {g: min(prof.n_eligible(g, nulls, s) for s in SEEDS) for g in all_judged}
        short = {g: v for g, v in short.items() if v < MIN_ELIG_NULLS}
        if not short or len(nulls) >= MAX_NULLS:
            print(f"補充 {rnd} 巡目: 不足 {len(short)} 件、対照遺伝子 {len(nulls)} 個で終了")
            break
        print(f"補充 {rnd} 巡目: 不足 {len(short)} 遺伝子 " +
              ", ".join(f"{g}({v})" for g, v in sorted(short.items(), key=lambda x: x[1])[:8]))
        progressed = False
        for g in sorted(short, key=lambda x: short[x]):
            if len(nulls) >= MAX_NULLS:
                break
            x = base[~base.symbol.isin(got)].copy()
            if x.empty:
                break
            x["dist"] = pd.concat([
                np.log(x.n_sel / max(meds[g][2], 1)).abs() / TOL_N["lncells"],
                (x.l2pos - meds[g][0]).abs() / TOL_N["log2pos"],
                (x.l2len - meds[g][1]).abs() / TOL_N["log2len"]], axis=1).max(1)
            need = MIN_ELIG_NULLS - min(prof.n_eligible(g, nulls, s) for s in SEEDS)
            for sym in x.sort_values(["dist", "symbol"]).symbol.head(N_AUG_SEARCH):
                if need <= 0 or len(nulls) >= MAX_NULLS:
                    break
                if all(prof.covers(g, sym, s) for s in SEEDS):
                    nulls.append(sym)
                    got.add(sym)
                    need -= 1
                    progressed = True
        if not progressed:
            print("補充できる対照遺伝子が尽きた")
            break

    elig = {g: {s: prof.n_eligible(g, nulls, s) for s in SEEDS} for g in all_judged}
    worst = min(min(v.values()) for v in elig.values())
    print(f"\n対照遺伝子 {len(nulls)} 個。使える対照遺伝子の最小 {worst}（下限 {MIN_ELIG_NULLS}）")

    # --- (D) 状態重心の細胞（MYOF を除く）
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
        (E / f"work/state_rows-seed{s}.json").write_text(
            json.dumps({"rows": rows}, ensure_ascii=False) + "\n")   # トップレベルを rows で包む（guardrails 2026-09-27）
        summary[s] = info_s
    (E / "work/depth_match_summary_001.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=1) + "\n")
    print("\n深さマッチ後の細胞数:", {s: {k: v["n_matched"] for k, v in summary[s].items()} for s in SEEDS})

    # --- 出力
    ens = dict(zip(d.symbol, d.ensembl))
    genes = []
    for axis, lst in chosen.items():
        for sym, half in lst:
            genes.append({"symbol": sym, "ensembl": ens[sym], "role": f"marker_{axis}", "half": half,
                          "bias_goal": round(float(ref.loc[sym, "bias_goal"]), 3),
                          "bias_myof": round(float(ref.loc[sym, "bias_myof"]), 3),
                          "n_start": int(ref.loc[sym, "n_start"]),
                          "n_donor_start": int(ref.loc[sym, "n_donor_start"])})
    for sym, e in HK.items():
        genes.append({"symbol": sym, "ensembl": e, "role": "housekeeping"})
    for sym, e in PANEL.items():
        genes.append({"symbol": sym, "ensembl": e, "role": "negctrl_panel"})
    for sym, e in REF_MYOF.items():
        genes.append({"symbol": sym, "ensembl": e, "role": "reference_myof"})
    for sym in nulls:
        genes.append({"symbol": sym, "ensembl": ens[sym], "role": "negative",
                      "bias_goal": round(float(ref.loc[sym, "bias_goal"]), 3),
                      "bias_myof": round(float(ref.loc[sym, "bias_myof"]), 3)})
    assert len({g["symbol"] for g in genes}) == len(genes), "遺伝子が重複している"
    assert not ({g["symbol"] for g in genes} & BLOCKED), "回さない遺伝子が混入している"
    assert not ({g["symbol"] for g in genes} & sg3.BLIND), "盲検の遺伝子が混入している"
    (E / "work/genes_001.json").write_text(json.dumps(
        {"rule": MARKER_RULE, "n": len(genes), "genes": genes}, ensure_ascii=False, indent=1) + "\n")
    pd.DataFrame([g for g in genes if g["role"].startswith("marker")]).to_csv(
        E / "work/marker_sets_001.csv", index=False)
    (E / "work/null_pool_001.json").write_text(json.dumps(
        {"n": len(nulls), "genes": [{"symbol": s, "ensembl": ens[s]} for s in nulls]},
        ensure_ascii=False, indent=1) + "\n")
    (E / "work/precheck_coverage_001.json").write_text(json.dumps(
        {"n_cells": {f"{g}@{s}": len(prof.get(g, s)[0]) for g in all_judged + nulls for s in SEEDS},
         "n_donors": {f"{g}@{s}": prof.n_donors(g, s) for g in all_judged + nulls for s in SEEDS},
         "n_eligible_null_genes": {f"{g}@{s}": elig[g][s] for g in all_judged for s in SEEDS},
         "n_genes_total": len(genes), "n_nulls": len(nulls),
         "min_eligible_among_judged": worst,
         "state_cells": info,
         "total_cell_gene": len(genes) * MAX_CELLS_PER_GENE * len(SEEDS)},
        ensure_ascii=False, indent=1) + "\n")

    # config.yaml に貼る形（実行スクリプトが読める書き方。guardrails 2026-09-27）
    lines = ["# 自動生成: analysis/select_genes_001.py。config.yaml へそのまま貼る"]
    for role, key in (("marker_goal", "genes_to_perturb"), ("marker_myof", "controls.positive"),
                      ("negative", "controls.negative_sampler.genes")):
        lines.append(f"# --- {key} ({role})")
        for g in genes:
            if g["role"] == role:
                lines.append(f"    - {g['ensembl']}  # {g['symbol']}")
    (E / "work/config_gene_lists_001.yaml").write_text("\n".join(lines) + "\n")
    print(f"\n合計 {len(genes)} 遺伝子（目標マーカー 24 / 筋線維芽細胞マーカー 24 / "
          f"恒常発現 3 / 中立パネル 7 / 参照 2 / 対照 {len(nulls)}）")
    print(f"延べ 細胞 × 遺伝子 = {len(genes) * MAX_CELLS_PER_GENE * len(SEEDS):,}")


if __name__ == "__main__":
    main()

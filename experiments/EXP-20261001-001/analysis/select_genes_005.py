#!/usr/bin/env python3
"""EXP-20261001-001 — 摂動する遺伝子と対照遺伝子群を決める（**結果を見る前に固定する**）。

使うのは (i) 仮説・選定手順・文献であらかじめ決まっている遺伝子の一覧と、(ii) 出発点・目標の
細胞だけから計算できる量（トークン化済み入力での検出細胞数・提供者数・入力中の位置・入力長・
提供者ごとの検出率）だけである。**摂動の結果は 1 つも使わない**。

## このスクリプトが EXP-20260930-001 の select_genes_004.py から変える 3 点

1. **介入セット 119 個を加える。** 出どころは選定手順
   `outputs/design/2026-10-01-intervention-gene-set/procedure.md`（**数え上げを 1 行も走らせる前に
   固定**。2026-10-01T19:35:00+09:00）の出力 `counts.json`。
   本スクリプトは counts.json を読むだけでなく、**元の経路セット 120 個
   （`experiments/EXP-20260928-002/work/pathway_set_002.json`）から除外の手順を自分でやり直して
   一致を確かめる**（1 個でも違えば停止する）。
2. **対照遺伝子群 500 個を引き直す。** 介入セット 119 個のうち **NFKBIA・PDE4B・REL・TNIP1 の
   4 個が EXP-20260930-001 の対照 500 個と重なる**。選定手順 §5.1 で**結果を見る前に決めたとおり、
   介入セットは削らず対照の側を引き直す**（介入は文献で決まり、対照は規則で決まるため）。
   引き直しの規則は EXP-20260930-001 と**同一**（|bias_goal| < 0.5、出発点で n_cap20 ≥ 120、
   判定対象と 1 個も重ねない、判定遺伝子ごとに近傍を 2 個ずつ確保した後に固定 seed で並べ替え）。
   除外集合に**介入セット 119 個と経路の候補表 406 個すべて**を加える。
3. **中立の対照の採用規則 H-1〜H-3 を数値で当てる**（plan.md §6.1.2）。
   恒常発現 4 個（ACTB・GAPDH・B2M・PPIA）と中立対照パネル 7 個の**11 個すべてに機械的に当て**、
   通ったものだけを中立性の判定 P-5 に使う。落ちたものは**記述にとどめる**。
   **結果を見てから外すのは事後の選び直しなので、閾値はここで先に固定する。**

出力:
  work/gene_stats_005.csv          全トークンの統計（検出率・偏り・位置・長さ・細胞数・提供者変動係数）
  work/genes_005.json              役割ごとの遺伝子（摂動する全遺伝子 717 個）
  work/null_pool_005.json          対照遺伝子 500（先頭 300 が入れ子の感度分析）
  work/subsets_005.json            判定に使うセット（結果を見る前に固定）
  work/neutral_rule_005.json       **H-1〜H-3 の当てはめ結果**（通った遺伝子・落ちた遺伝子）
  work/overlap_check_005.json      **重なり 0 件の確認**（総当たり）
  work/precheck_coverage_005.json  遺伝子ごとの抽出細胞数・提供者数（3 seed）
"""
from __future__ import annotations

import json
import random
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20261001-001"
P30 = ROOT / "experiments/EXP-20260930-001"
P29 = ROOT / "experiments/EXP-20260929-001"
P28 = ROOT / "experiments/EXP-20260928-002"
DESIGN = ROOT / "outputs/design/2026-10-01-intervention-gene-set"
LIT = ROOT / "outputs/literature/2026-09-30-independent-goal-state/markers"
sys.path.insert(0, str(ROOT / "scripts/geneformer"))
import run_isp as ri  # noqa: E402

SEEDS = [81, 82, 83]
MIN_CELLS, MAX_CELLS, MAX_PER_DONOR = 120, 200, 20
GOAL, START = "aStr", "FIB_MFIB"

# 対照遺伝子の条件（EXP-20260930-001 と**同一**。選定手順 §5.1）
NULL_ABS_BIAS_MAX = 0.5
NULL_MIN_DETECT, NULL_MIN_DONORS = 250, 12
N_NULL_PRIMARY, N_NULL_TOTAL = 300, 500
N_RESERVE_PER_JUDGED = 2
NULL_ORDER_SEED = 20261001        # **並べ替えの seed。結果を見る前に固定**（引き直しなので更新する）
TOL_N = {"lncells": float(np.log(1.5)), "log2pos": 0.25, "log2len": 0.25}

# ---------------------------------------------------------------------------
# 中立の対照の採用規則 H-1〜H-3（plan.md §6.1.2。**結果を見る前に数値で固定する**）
#
# 根拠（なぜこの数値か）。閾値は「どの遺伝子が通るか」ではなく、
# **『どの細胞でも一定』という前提が、この単一核データで検査できる水準か**から決めた。
#
#  H-1 検出率 ≥ 0.20（出発点・目標の両方）
#      中立の対照の前提は「細胞の種類によらず入力に載っている」こと。
#      5 細胞に 1 つも載らない遺伝子では、その前提自体がこのデータで確かめられない。
#      提供者あたり 20 細胞という抽出の設計上、0.20 は 1 提供者あたり期待 4 細胞に当たり、
#      提供者ごとの検出率をばらつきの推定ごと扱える下限でもある。
#  H-2 |log2(目標の検出率 / 出発点の検出率)| ≤ 0.25
#      対照遺伝子群に入れる条件は |bias_goal| < 0.5 である（EXP-20260930-001 以来）。
#      **中立であることを示すために使う遺伝子**は、対照の母集団に入れてよい水準の
#      **2 倍きびしく**なければ役に立たない。log2 で 0.25 は検出率で約 19% の差。
#  H-3 提供者間の検出率の変動係数 ≤ 0.50（出発点・目標の悪いほう）
#      変動係数 0.5 は「提供者ごとの検出率が平均の ±50% ほどばらつく」水準。
#      これを超えると『どの細胞でも一定』は提供者の段階ですでに成り立っておらず、
#      提供者上限 20 細胞の抽出で提供者が交絡として効いてくる。
#
# **この 3 つは、現行の恒常発現 4 個にも中立対照パネル 7 個にも、新しい候補にも、
# まったく同じように当てる。** 落ちた遺伝子は記述にとどめ、判定には使わない。
H1_DETECT_MIN = 0.20
H2_ABS_LOG2_BIAS_MAX = 0.25
H3_DONOR_CV_MAX = 0.50
H3_MIN_CELLS_PER_DONOR = 20       # 提供者ごとの検出率を数えるのに必要な、その状態での細胞数
NEUTRAL_MIN_ADOPTED = 3           # 採用された中立対照がこれを下回ったら着手しない
#   3 は過去実験の中立性セット（ACTB・GAPDH・B2M の 3 個）と同じ大きさである。
#   **初版ではここに 5 と書いたが、根拠のない数字だったので過去実験と同じ 3 に直した。**
#   直したのは H-1〜H-3 ではなく「何個あれば検査できるか」の下限だけで、
#   **採用される遺伝子の顔ぶれを決める規則は 1 文字も変えていない**（plan.md §11 の変更ログ）。

# **H-3 の僅差に対する感度分析**（結果を見る前に固定）。
#   H-3 を 0.50 にすると B2M（0.5002）と VAPA（0.5015）が **0.0002 / 0.0015 の差で落ちる**。
#   閾値を動かせばこの 2 個は通る。**規則は動かさない**が、動かした場合の版を
#   入れ子の感度分析として先に登録しておく（後から選び直す余地をなくすため）。
H3_DONOR_CV_MAX_SENSITIVITY = 0.55


def read_marker_csv(p: Path) -> list[str]:
    return [x.strip() for x in p.read_text().splitlines()[1:] if x.strip()]


def donor_cv(Mstate, donor_codes, min_cells: int) -> np.ndarray:
    """状態 S の細胞（行）× トークン（列）の行列から、提供者ごとの検出率の変動係数を返す。

    提供者は、その状態で `min_cells` 細胞以上を持つものだけを数える。
    """
    n_don = int(donor_codes.max()) + 1
    cnt = np.bincount(donor_codes, minlength=n_don).astype(float)
    use = cnt >= min_cells
    # 提供者 × トークン の検出細胞数
    Dm = sp.csr_matrix((np.ones(len(donor_codes), float),
                        (donor_codes, np.arange(len(donor_codes)))),
                       shape=(n_don, len(donor_codes)))
    num = np.asarray((Dm @ Mstate).todense())           # 提供者ごとの検出細胞数
    rate = num[use] / cnt[use][:, None]
    mean = rate.mean(axis=0)
    sd = rate.std(axis=0, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        cv = np.where(mean > 0, sd / mean, np.nan)
    return cv, int(use.sum())


def main() -> None:
    from datasets import load_from_disk

    G = json.loads((E / "analysis/gene_sets_005.json").read_text())
    G30 = json.loads((P30 / "analysis/gene_sets_004.json").read_text())
    for k in G30:
        if k != "note" and G[k] != G30[k]:
            sys.exit(f"gene_sets_005.json の {k} が EXP-20260930-001 と違う（変えてはいけない）")
    sym2ens = json.loads((E / "work/symbol_to_ensembl_005.json").read_text())
    tokd = ri.gene_token_dict()
    ens2sym: dict[str, str] = {}
    for s, e in sym2ens.items():
        ens2sym.setdefault(e, s)

    # ---------------------------------------------------------------- 文献との照合（陽性対照）
    lake = json.loads((LIT / "lake2023_fibroblast_state_markers.json").read_text())
    mini42 = read_marker_csv(LIT / "adaptive_stromal_state_aFIB_plus_MYOF_minimal.csv")
    full76 = read_marker_csv(LIT / "adaptive_stromal_state_aFIB_plus_MYOF_full_conserved.csv")
    if not set(mini42) <= set(full76):
        sys.exit("最小 42 が保守 76 の部分集合になっていない")

    def usable(g: str) -> bool:
        return g in sym2ens and sym2ens[g] in tokd

    pc39 = [g for g in mini42 if usable(g)]
    pc71 = [g for g in full76 if usable(g)]
    if sorted(pc39) != sorted(G["positive_control_39"]) or len(pc71) != 71:
        sys.exit("陽性対照が EXP-20260930-001 と一致しない（変えてはいけない）")
    pc_extra32 = [g for g in pc71 if g not in set(pc39)]

    # ---------------------------------------------------------------- 入力から数える
    ds = load_from_disk(str(E / "work/tokenized.dataset"))
    st = np.array([str(x) for x in ds["fib_anno_state"]])
    donors = np.array([str(x) for x in ds["donor_id"]])
    ids = [np.asarray(x) for x in ds["input_ids"]]
    lens = np.array([len(x) for x in ids], float)

    indptr = np.cumsum([0] + [len(x) for x in ids])
    col = np.concatenate(ids)
    ntok = max(int(col.max()) + 1, max(tokd.values()) + 1)
    posdata = np.concatenate([np.arange(1, len(x) + 1) for x in ids]).astype(np.int32)
    Mb = sp.csr_matrix((np.ones(len(col), bool), col, indptr), shape=(len(ids), ntok)).tocsc()
    Mp = sp.csr_matrix((posdata, col, indptr), shape=(len(ids), ntok)).tocsc()

    sets = {"start": st == START, "goal": st == GOAL, "dfib": st == "dFIB", "dmfib": st == "dMFIB"}
    d = pd.DataFrame({"token": np.arange(ntok)})
    for k, v in sets.items():
        n = np.asarray(Mb[v].sum(0)).ravel()
        d[f"n_{k}"] = n
        d[f"det_{k}"] = n / v.sum()
    eps = 1.0 / int(sets["goal"].sum())
    for k in ("goal", "dfib", "dmfib"):
        d[f"bias_{k}"] = np.log2((d[f"det_{k}"] + eps) / (d.det_start + eps))

    sidx = np.flatnonzero(sets["start"])
    Ms = Mb[sidx].tocsc()
    Mps = Mp[sidx].tocsc()
    dcodes = pd.Categorical(donors[sidx]).codes
    ndon = np.zeros(ntok, int)
    ncap = np.zeros(ntok, int)
    lp_med = np.full(ntok, np.nan)
    ll_med = np.full(ntok, np.nan)
    slen = lens[sidx]
    for t in range(ntok):
        r = Ms.indices[Ms.indptr[t]:Ms.indptr[t + 1]]
        if len(r) == 0:
            continue
        cnt = np.bincount(dcodes[r])
        ndon[t] = int((cnt > 0).sum())
        ncap[t] = int(min(MAX_CELLS, np.minimum(cnt, MAX_PER_DONOR).sum()))
        p = Mps.data[Mps.indptr[t]:Mps.indptr[t + 1]].astype(float)
        lp_med[t] = float(np.median(np.log2(p)))
        ll_med[t] = float(np.median(np.log2(slen[r])))
    d["n_donor_start"], d["n_cap20"] = ndon, ncap
    d["lp_med"], d["ll_med"] = lp_med, ll_med

    # **H-3 用: 提供者ごとの検出率の変動係数**（出発点・目標のそれぞれ）
    cv_start, n_don_start_used = donor_cv(
        Mb[sidx].tocsr().astype(float), pd.Categorical(donors[sidx]).codes,
        H3_MIN_CELLS_PER_DONOR)
    gidx = np.flatnonzero(sets["goal"])
    cv_goal, n_don_goal_used = donor_cv(
        Mb[gidx].tocsr().astype(float), pd.Categorical(donors[gidx]).codes,
        H3_MIN_CELLS_PER_DONOR)
    d["donor_cv_start"], d["donor_cv_goal"] = cv_start, cv_goal

    tok2ens = {v: k for k, v in tokd.items()}
    d["ensembl"] = [tok2ens.get(int(t), "") for t in d.token]
    d["symbol"] = [ens2sym.get(e, "") for e in d.ensembl]
    d = d[d.ensembl != ""].reset_index(drop=True)
    d.to_csv(E / "work/gene_stats_005.csv", index=False)

    # EXP-20260930-001 と同じ入力・同じ集計なので、共通の列は 1 つも違わないはず
    prev = pd.read_csv(P30 / "work/gene_stats_004.csv")
    same = ["token", "n_start", "det_start", "n_goal", "det_goal", "bias_goal",
            "n_donor_start", "n_cap20", "lp_med", "ll_med"]
    a = d.set_index("token")[same[1:]].sort_index()
    b = prev.set_index("token")[same[1:]].sort_index()
    if not a.reindex(b.index).fillna(-1).round(6).equals(b.fillna(-1).round(6)):
        sys.exit("gene_stats が EXP-20260930-001 と一致しない（入力が変わっていないはずなのに）")

    ref = d.set_index("symbol")
    ref = ref[~ref.index.duplicated()]

    def cap20(g: str) -> int:
        return int(ref.loc[g].n_cap20) if g in ref.index else -1

    # ---------------------------------------------------------------- 介入セット（選定手順の出力）
    counts = json.loads((DESIGN / "counts.json").read_text())
    p1 = counts["候補"]["P1_経路セット120"]
    INT = [x["symbol"] for x in p1["採用した遺伝子"]]
    tier = {x["symbol"]: x["tier"] for x in p1["採用した遺伝子"]}
    group = {x["symbol"]: x["group"] for x in p1["採用した遺伝子"]}
    if len(INT) != 119 or p1["K"] != 119:
        sys.exit(f"介入セットが 119 個でない（{len(INT)}）")

    # **選定手順 §5 の除外を自分でやり直して、counts.json と一致することを確かめる**
    pathway120 = [x["symbol"] for x in json.loads((P28 / "work/pathway_set_002.json").read_text())["genes"]]
    if len(pathway120) != 120:
        sys.exit("経路セットが 120 個でない")
    step, remain = [], list(pathway120)
    for label, drop in (
            ("1. 盲検中の 3 個", set(G["blinded_still_closed_3"])),
            ("2. 陽性対照 71 と重なる", set(pc71)),
            ("3. 目標の定義と重なる", set()),            # 目標は著者注釈。遺伝子で定義していない
            ("4. 恒常発現・中立対照パネルと重なる",
             set(G["housekeeping_3"]) | set(G["housekeeping_extra_descriptive"])
             | set(G["negctrl_panel_7"])),
            ("5. gene_stats に無い（var / トークン辞書に無い）", None),
            ("6. 出発点で上限 20 後 120 細胞未満", None)):
        if drop is None:
            lost = ([g for g in remain if g not in ref.index] if "gene_stats" in label
                    else [g for g in remain if cap20(g) < MIN_CELLS])
        else:
            lost = [g for g in remain if g in drop]
        remain = [g for g in remain if g not in set(lost)]
        step.append({"手順": label, "落ちた数": len(lost), "残り": len(remain), "落ちた遺伝子": sorted(lost)})
    if remain != INT:
        sys.exit(f"除外の手順をやり直した結果が counts.json と一致しない: "
                 f"{sorted(set(remain) ^ set(INT))}")
    recomputed_steps = step

    INT_S1 = [g for g in INT if tier[g] == "S1"]
    INT_S2 = [g for g in INT if tier[g] == "S2"]
    INT_S3 = [g for g in INT if tier[g] == "S3"]
    INT_SPEC = [g for g in INT if group[g] == "経路特異的"]
    T1, T2 = G["perturb_T1_11"], G["perturb_T2_7"]
    INT_T1T2 = INT + [g for g in T1 + T2 if g not in set(INT)]
    sens = counts["入れ子の感度分析"]
    for name, got, want in (("S1", INT_S1, sens["感度a_S1経路の構成因子"]["K"]),
                            ("S2", INT_S2, sens["感度b_S2経路の転写標的"]["K"]),
                            ("S3", INT_S3, sens["感度c_S3キュレーション"]["K"]),
                            ("経路特異的", INT_SPEC, sens["感度d_経路特異的のみ"]["K"]),
                            ("∪T1T2", INT_T1T2, sens["感度e_主判定とT1T2の和"]["K"])):
        if len(got) != want:
            sys.exit(f"感度分析 {name} の個数が counts.json と違う: {len(got)} ≠ {want}")

    # ---------------------------------------------------------------- 摂動できるかの確認
    for name in ("perturb_T1_11", "perturb_T2_7", "housekeeping_3", "negctrl_panel_7"):
        bad = [g for g in G[name] if cap20(g) < MIN_CELLS]
        if bad:
            sys.exit(f"{name}: 下限 {MIN_CELLS} 細胞に届かない遺伝子がある {bad}")
    for g in pc71 + INT:
        if cap20(g) < MIN_CELLS:
            sys.exit(f"{g} が下限 {MIN_CELLS} 細胞に届かない")

    # ---------------------------------------------------------------- 中立の対照の採用規則 H-1〜H-3
    neutral_cand = (G["housekeeping_3"] + G["housekeeping_extra_descriptive"] + G["negctrl_panel_7"])
    cat = ({g: "恒常発現" for g in G["housekeeping_3"] + G["housekeeping_extra_descriptive"]}
           | {g: "中立対照パネル" for g in G["negctrl_panel_7"]})
    neutral: dict[str, dict] = {}
    for g in neutral_cand:
        r = ref.loc[g]
        h1 = bool(min(r.det_start, r.det_goal) >= H1_DETECT_MIN)
        h2 = bool(abs(r.bias_goal) <= H2_ABS_LOG2_BIAS_MAX)
        cv = float(np.nanmax([r.donor_cv_start, r.donor_cv_goal]))
        h3 = bool(cv <= H3_DONOR_CV_MAX)
        neutral[g] = {"区分": cat[g], "det_start": round(float(r.det_start), 4),
                      "det_goal": round(float(r.det_goal), 4),
                      "bias_goal": round(float(r.bias_goal), 4),
                      "donor_cv_start": round(float(r.donor_cv_start), 4),
                      "donor_cv_goal": round(float(r.donor_cv_goal), 4),
                      "donor_cv_max": round(cv, 4),
                      "H1_検出率": h1, "H2_偏り": h2, "H3_提供者変動係数": h3,
                      "採用": bool(h1 and h2 and h3)}
    adopted = [g for g in neutral_cand if neutral[g]["採用"]]
    rejected = [g for g in neutral_cand if not neutral[g]["採用"]]
    # H-3 を 0.55 に緩めた版（感度分析。**主判定には使わない**）
    adopted_cv055 = [g for g in neutral_cand
                     if neutral[g]["H1_検出率"] and neutral[g]["H2_偏り"]
                     and neutral[g]["donor_cv_max"] <= H3_DONOR_CV_MAX_SENSITIVITY]
    (E / "work/neutral_rule_005.json").write_text(json.dumps(
        {"規則": {"H-1 出発点・目標の両方での検出率": f"≥ {H1_DETECT_MIN}",
                  "H-2 |log2(目標の検出率/出発点の検出率)|": f"≤ {H2_ABS_LOG2_BIAS_MAX}",
                  "H-3 提供者間の検出率の変動係数（悪いほう）": f"≤ {H3_DONOR_CV_MAX}",
                  "提供者を数える下限": f"その状態で {H3_MIN_CELLS_PER_DONOR} 細胞以上",
                  "注": "**結果を見る前に固定した。** 恒常発現にも中立対照パネルにも新しい候補にも"
                        "同じように当てる。落ちた遺伝子は記述にとどめ、判定には使わない"},
         "数えた提供者数": {"出発点": n_don_start_used, "目標": n_don_goal_used},
         "当てはめ": neutral,
         "採用": adopted, "落ちた": rejected,
         "採用の内訳": {"恒常発現": [g for g in adopted if cat[g] == "恒常発現"],
                        "中立対照パネル": [g for g in adopted if cat[g] == "中立対照パネル"]},
         "H-3 を 0.55 にした版（感度分析。主判定には使わない）": {
             "閾値": H3_DONOR_CV_MAX_SENSITIVITY, "採用": adopted_cv055,
             "加わる遺伝子": [g for g in adopted_cv055 if g not in set(adopted)],
             "注": "H-3 0.50 では B2M が 0.5002、VAPA が 0.5015 で**ごくわずかに**落ちる。"
                   "**閾値は動かさない**が、動かした場合を結果を見る前に登録しておく"},
         "P-5 に使うセット": "neutral_adopted",
         "必要な最小個数": NEUTRAL_MIN_ADOPTED},
        indent=1, ensure_ascii=False) + "\n")
    if len(adopted) < NEUTRAL_MIN_ADOPTED:
        sys.exit(f"採用された中立対照が {len(adopted)} 個しかない（必要 {NEUTRAL_MIN_ADOPTED}）。"
                 "設計を見直すこと")

    # ---------------------------------------------------------------- 除外集合
    prev_groups = json.loads((P29 / "analysis/excluded_symbols_003.json").read_text())["groups"]
    prev_states = json.loads((P29 / "work/input_states.json").read_text())
    lake_all = sorted({g for v in lake["states"].values() for g in v["genes"]})
    pathway_all = [x["symbol"] for x in p1["採用した遺伝子"]] + \
                  [x["symbol"] for x in counts["候補"]["P2_候補表の全体"]["採用した遺伝子"]]
    judged = (INT + pc71 + T1 + T2 + G["housekeeping_3"]
              + G["housekeeping_extra_descriptive"] + G["negctrl_panel_7"])
    judged = list(dict.fromkeys(judged))         # CCL2・TIMP1 は介入と T1/T2 の両方にある

    excluded = (set(judged) | set(lake_all) | set(mini42) | set(full76)
                | set(pathway_all)                       # **経路の候補表すべて**（介入の近縁を外す）
                | set(G["perturb_T1_candidates_16"]) | set(G["perturb_T2_candidates_12"])
                | set(G["blinded_5"]) | set(G["reference_dFIB_markers_6_descriptive"])
                | set(prev_groups["secretory25"]) | set(prev_groups["org_markers7"])
                | set(prev_groups["myof_markers24"]) | set(prev_groups["housekeeping3"])
                | set(prev_groups["negctrl_panel7"]) | set(prev_groups["reference2"])
                | set(prev_states["ecm_candidates_order"]) | set(prev_states["ecm_program"]))

    # ---------------------------------------------------------------- 対照遺伝子を引き直す
    cand = d[(~d.symbol.isin(excluded)) & (d.symbol != "")
             & (d.n_start >= NULL_MIN_DETECT) & (d.n_donor_start >= NULL_MIN_DONORS)
             & (d.n_cap20 >= MIN_CELLS)
             & (d.bias_goal.abs() < NULL_ABS_BIAS_MAX)].copy()
    cand = cand.set_index("symbol")
    cand = cand[~cand.index.duplicated()]

    order: list[str] = []
    taken: set[str] = set()
    for s in judged:                       # 判定する遺伝子に性質の近い対照を 2 個ずつ確保する
        r = ref.loc[s]
        av = cand[~cand.index.isin(taken)]
        dist = pd.DataFrame({
            "a": (np.log(av.n_cap20) - np.log(r.n_cap20)).abs() / TOL_N["lncells"],
            "b": (av.lp_med - r.lp_med).abs() / TOL_N["log2pos"],
            "c": (av.ll_med - r.ll_med).abs() / TOL_N["log2len"]}).max(axis=1)
        for s2 in dist.nsmallest(N_RESERVE_PER_JUDGED).index:
            order.append(str(s2)); taken.add(str(s2))
    rest = [s for s in cand.index if s not in taken]
    random.Random(NULL_ORDER_SEED).shuffle(rest)
    order += rest
    if len(order) < N_NULL_TOTAL:
        sys.exit(f"対照遺伝子が足りない（{len(order)} < {N_NULL_TOTAL}）")
    null500 = order[:N_NULL_TOTAL]
    null300 = null500[:N_NULL_PRIMARY]

    # ---------------------------------------------------------------- 重なりの総当たり確認
    named = {"【介入】主判定 119": INT,
             "（内訳）S1 経路の構成因子 61": INT_S1,
             "（内訳）S2 経路の転写標的 13": INT_S2,
             "（内訳）S3 キュレーション 45": INT_S3,
             "（内訳）経路特異的 114": INT_SPEC,
             "介入 ∪ T1 ∪ T2 135": INT_T1T2,
             "陽性対照 71": pc71, "（内訳）最小 42 版の 39": pc39,
             "（内訳）保守版で加わる 32": pc_extra32,
             "前回の摂動対象 T1 11": T1, "前回の摂動対象 T2 7": T2,
             "恒常発現 3": G["housekeeping_3"],
             "恒常発現の追加 1": G["housekeeping_extra_descriptive"],
             "中立対照パネル 7": G["negctrl_panel_7"],
             "採用した中立対照": adopted,
             "採用した中立対照（H-3 0.55 版）": adopted_cv055,
             "対照遺伝子 500": null500, "対照遺伝子 300": null300}
    # **入れ子の族**（重なって当然なので総当たりの検査から外し、代わりに包含関係を検算する）。
    # 介入セットの内訳（S1/S2/S3 は層、経路特異的は群）は、同じ 119 個を 2 通りに切った
    # ものなので互いに重なる。重なりの有無ではなく「すべて INT の部分集合か」で確かめる。
    INT_FAMILY = ["【介入】主判定 119", "（内訳）S1 経路の構成因子 61",
                  "（内訳）S2 経路の転写標的 13", "（内訳）S3 キュレーション 45",
                  "（内訳）経路特異的 114", "介入 ∪ T1 ∪ T2 135"]
    PC_FAMILY = ["陽性対照 71", "（内訳）最小 42 版の 39", "（内訳）保守版で加わる 32"]
    NEU_FAMILY = ["恒常発現 3", "恒常発現の追加 1", "中立対照パネル 7",
                  "採用した中立対照", "採用した中立対照（H-3 0.55 版）"]
    NULL_FAMILY = ["対照遺伝子 500", "対照遺伝子 300"]
    FAMILIES = [INT_FAMILY, PC_FAMILY, NEU_FAMILY, NULL_FAMILY]
    nested_checks = {
        "S1 ∪ S2 ∪ S3 = 介入 119": sorted(set(INT_S1) | set(INT_S2) | set(INT_S3)) == sorted(INT),
        "S1・S2・S3 は互いに素": len(INT_S1) + len(INT_S2) + len(INT_S3) == len(INT),
        "経路特異的 ⊂ 介入": set(INT_SPEC) <= set(INT),
        "介入 ⊂ 介入∪T1∪T2": set(INT) <= set(INT_T1T2),
        "T1 ⊂ 介入∪T1∪T2": set(T1) <= set(INT_T1T2),
        "T2 ⊂ 介入∪T1∪T2": set(T2) <= set(INT_T1T2),
        "39 ⊂ 71": set(pc39) <= set(pc71),
        "65 ⊂ 71": set(subs_pc65 := [g for g in pc71
                                     if g not in set(G["positive_control71_opposite_6"])]) <= set(pc71),
        "採用した中立対照 ⊂ 候補 11": set(adopted) <= set(neutral_cand),
        "採用（H-3 0.55）⊃ 採用": set(adopted) <= set(adopted_cv055),
        "対照 300 ⊂ 対照 500": set(null300) <= set(null500)}
    if not all(nested_checks.values()):
        sys.exit(f"入れ子の包含関係が成り立たない: "
                 f"{[k for k, v in nested_checks.items() if not v]}")

    # **介入と T1/T2 の重なり（CCL2・TIMP1）は循環ではない**（どちらも過去の介入対象であって、
    # 目標の定義にも対照にも使われていない）。重なりを許すが、件数は必ず記録して報告する。
    ALLOWED = {("【介入】主判定 119", "前回の摂動対象 T1 11"): {"CCL2"},
               ("【介入】主判定 119", "前回の摂動対象 T2 7"): {"TIMP1"},
               ("（内訳）S2 経路の転写標的 13", "前回の摂動対象 T1 11"): {"CCL2"},
               ("（内訳）S2 経路の転写標的 13", "前回の摂動対象 T2 7"): {"TIMP1"},
               ("（内訳）経路特異的 114", "前回の摂動対象 T1 11"): {"CCL2"},
               ("（内訳）経路特異的 114", "前回の摂動対象 T2 7"): {"TIMP1"}}
    # 組の順序に依らないように frozenset で引く（並べ方を間違えて検査が抜けないため）
    NESTED_FS = {frozenset((a, b)) for fam in FAMILIES for a, b in combinations(fam, 2)}
    # 感度 (e) は定義上 T1・T2 を丸ごと含む（包含関係は nested_checks で検算済み）
    NESTED_FS |= {frozenset(("介入 ∪ T1 ∪ T2 135", "前回の摂動対象 T1 11")),
                  frozenset(("介入 ∪ T1 ∪ T2 135", "前回の摂動対象 T2 7"))}
    ALLOWED_FS = {frozenset(k): v for k, v in ALLOWED.items()}
    pairs, violations = {}, {}
    for a, b in combinations(named, 2):
        if frozenset((a, b)) in NESTED_FS:
            continue
        ov = sorted(set(named[a]) & set(named[b]))
        pairs[f"{a} × {b}"] = ov
        if set(ov) - ALLOWED_FS.get(frozenset((a, b)), set()):
            violations[f"{a} × {b}"] = ov
    blinded_hit = sorted(set(G["blinded_still_closed_3"]) & {g for v in named.values() for g in v})
    prev_null = {x["symbol"] for x in
                 json.loads((P30 / "work/null_pool_004.json").read_text())["versions"]["B_500"]}
    overlap = {"各セットの大きさ": {k: len(v) for k, v in named.items()},
               "入れ子の族（総当たりから外し、包含関係で検算した）":
                   {"介入": INT_FAMILY, "陽性対照": PC_FAMILY, "中立対照": NEU_FAMILY,
                    "対照遺伝子": NULL_FAMILY},
               "包含関係の検算": nested_checks,
               "総当たりの重なり": pairs,
               "重なりの合計": int(sum(len(v) for v in pairs.values())),
               "違反（許していない重なり）": violations,
               "許している重なり": {f"{a} × {b}": sorted(v) for (a, b), v in ALLOWED.items()},
               "許す理由": "CCL2・TIMP1 は介入セットと前回の摂動対象 T1・T2 の両方にあるが、"
                           "どちらも過去の介入対象であって目標の定義にも対照にも使われていない。"
                           "循環ではない（HYP-20261001-01 §6.3）。感度 (e) では重複して数えない",
               "盲検中の 3 遺伝子の混入": blinded_hit,
               "対照 500 × 介入 119 の重なり": sorted(set(null500) & set(INT)),
               "引き直す前（EXP-20260930-001 の対照 500）× 介入 119 の重なり":
                   sorted(prev_null & set(INT)),
               "引き直した対照 500 と EXP-20260930-001 の対照 500 の共通": len(prev_null & set(null500))}
    (E / "work/overlap_check_005.json").write_text(
        json.dumps(overlap, indent=1, ensure_ascii=False) + "\n")
    if violations or blinded_hit:
        sys.exit(f"許していない重なりがある: {violations} / 盲検混入: {blinded_hit}")

    # ---------------------------------------------------------------- セットと役割
    subsets = {
        "INT": INT,                                       # **主判定 P-1**
        "INT_S1": INT_S1, "INT_S2": INT_S2, "INT_S3": INT_S3,       # 入れ子の感度分析
        "INT_pathway_specific": INT_SPEC, "INT_plus_T1T2": INT_T1T2,
        "positive_control71": pc71,                       # **P-4 健全性チェック**
        "positive_control65": [g for g in pc71
                               if g not in set(G["positive_control71_opposite_6"])],
        "positive_control39": pc39,
        "positive_control70": [g for g in pc71 if g != "LTBP1"],
        "positive_control38": [g for g in pc39 if g != "LTBP1"],
        "positive_control35": [g for g in pc39
                               if g not in set(G["positive_control_opposite_4"])],
        "neutral_adopted": adopted,                       # **P-5 中立性（規則 H-1〜H-3 を通ったもの）**
        "neutral_adopted_cv055": adopted_cv055,           # 入れ子の感度分析（H-3 を 0.55 に）
        "housekeeping": G["housekeeping_3"],              # 記述のみ
        "negctrl_panel": G["negctrl_panel_7"],            # 記述のみ
        "T1": T1, "T2": T2,                               # 記述のみ（前回の摂動対象）
    }

    roles: dict[str, str] = {}
    for g in null500:
        roles[g] = "negative"
    for g in G["negctrl_panel_7"]:
        roles[g] = "negctrl_panel"
    for g in G["housekeeping_extra_descriptive"]:
        roles[g] = "housekeeping_extra"
    for g in G["housekeeping_3"]:
        roles[g] = "housekeeping"
    for g in pc71:
        roles[g] = "positive_control"
    for g in T2:
        roles[g] = "target_T2"
    for g in T1:
        roles[g] = "target_T1"
    for g in INT:
        roles[g] = "intervention"          # **介入が最優先**（CCL2・TIMP1 はここに入る）
    allg = judged + null500
    if len(allg) != len(set(allg)):
        sys.exit("摂動する遺伝子が重複している")
    if set(allg) & set(G["blinded_still_closed_3"]):
        sys.exit("盲検中の遺伝子が入っている")

    def rec(s: str) -> dict:
        r = ref.loc[s]
        o = {"symbol": s, "ensembl": sym2ens[s], "role": roles[s],
             "bias_goal": round(float(r.bias_goal), 3),
             "bias_dfib": round(float(r.bias_dfib), 3),
             "n_start": int(r.n_start), "n_goal": int(r.n_goal),
             "det_start": round(float(r.det_start), 4), "det_goal": round(float(r.det_goal), 4),
             "n_donor_start": int(r.n_donor_start), "n_cap20": int(r.n_cap20),
             "lp_med": round(float(r.lp_med), 4), "ll_med": round(float(r.ll_med), 4)}
        if s in tier:
            o |= {"tier": tier[s], "group": group[s]}
        return o

    (E / "work/genes_005.json").write_text(json.dumps(
        {"genes": [rec(s) for s in allg], "n_total": len(allg),
         "by_role": {r: sorted(g for g in allg if roles[g] == r) for r in sorted(set(roles.values()))},
         "介入セットの除外の手順（やり直して counts.json と一致を確認）": recomputed_steps},
        indent=1, ensure_ascii=False) + "\n")
    (E / "work/null_pool_005.json").write_text(json.dumps(
        {"order_seed": NULL_ORDER_SEED, "n_candidates": int(len(cand)),
         "reserved_per_judged": N_RESERVE_PER_JUDGED, "n_reserved": len(taken),
         "n_excluded_symbols": len(excluded),
         "conditions": {"abs_bias_goal_max": NULL_ABS_BIAS_MAX, "min_detect_start": NULL_MIN_DETECT,
                        "min_donors_start": NULL_MIN_DONORS, "min_cap20": MIN_CELLS},
         "versions": {"A_300": [rec(s) for s in null300], "B_500": [rec(s) for s in null500]},
         "nested": null500[:N_NULL_PRIMARY] == null300}, indent=1, ensure_ascii=False) + "\n")
    (E / "work/subsets_005.json").write_text(json.dumps(
        {"subsets": {k: [rec(s) for s in v] for k, v in subsets.items()},
         "sizes": {k: len(v) for k, v in subsets.items()}}, indent=1, ensure_ascii=False) + "\n")

    # ---------------------------------------------------------------- seed ごとの抽出細胞数
    cov = {}
    for s in allg:
        e = sym2ens[s]
        t = tokd[e]
        idx = Ms.indices[Ms.indptr[t]:Ms.indptr[t + 1]].tolist()
        per = {}
        for sd in SEEDS:
            keep = ri.balanced_subset(idx, list(donors[sidx]), MAX_CELLS, MAX_PER_DONOR,
                                      random.Random(f"{sd}-{e}"))
            per[str(sd)] = {"n_cells": len(keep), "n_donors": len({donors[sidx][i] for i in keep})}
        cov[s] = {"ensembl": e, "role": roles[s], "per_seed": per,
                  "min_cells": min(v["n_cells"] for v in per.values())}
    (E / "work/precheck_coverage_005.json").write_text(
        json.dumps(cov, indent=1, ensure_ascii=False) + "\n")
    bad = {s: v["min_cells"] for s, v in cov.items() if v["min_cells"] < MIN_CELLS}

    print(f"対照候補 {len(cand)} 個 → 対照 {N_NULL_TOTAL}（先頭 {N_NULL_PRIMARY} が入れ子）")
    print("除外した記号:", len(excluded))
    print("摂動する遺伝子:", json.dumps(
        {r: len([g for g in allg if roles[g] == r]) for r in sorted(set(roles.values()))},
        ensure_ascii=False), "= 合計", len(allg))
    print("セットの大きさ:", json.dumps({k: len(v) for k, v in subsets.items()}, ensure_ascii=False))
    print("\n中立の対照の採用規則 H-1〜H-3:")
    for g in neutral_cand:
        v = neutral[g]
        print(f"  {g:10s} {v['区分']:8s} det {v['det_start']:.3f}/{v['det_goal']:.3f} "
              f"bias {v['bias_goal']:+.3f} CV {v['donor_cv_max']:.3f} "
              f"H1={v['H1_検出率']} H2={v['H2_偏り']} H3={v['H3_提供者変動係数']} "
              f"→ {'採用' if v['採用'] else '**落選（記述のみ）**'}")
    print("採用:", adopted, f"（{len(adopted)} 個）／落選:", rejected)
    print("H-3 を 0.55 にした版（感度分析）:", adopted_cv055,
          "← 加わるのは", [g for g in adopted_cv055 if g not in set(adopted)])
    print("\n重なりの合計:", overlap["重なりの合計"], "（うち許している重なり",
          sum(len(v) for v in ALLOWED.values()), "件）／違反:", violations or "なし",
          "／盲検混入:", blinded_hit or "なし")
    print("対照 500 × 介入 119 の重なり:", overlap["対照 500 × 介入 119 の重なり"] or "0 件")
    print("引き直す前の対照 500 × 介入 119:",
          overlap["引き直す前（EXP-20260930-001 の対照 500）× 介入 119 の重なり"])
    print("引き直した対照 500 と前回の共通:",
          overlap["引き直した対照 500 と EXP-20260930-001 の対照 500 の共通"], "個")
    print("120 細胞を下回る遺伝子:", bad if bad else "なし")
    if bad:
        sys.exit("抽出細胞数が下限を下回る遺伝子がある。設計を見直すこと")


if __name__ == "__main__":
    main()

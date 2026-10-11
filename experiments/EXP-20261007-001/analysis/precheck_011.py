#!/usr/bin/env python3
"""EXP-20261007-001 の設計に必要な量を、**摂動を 1 回も回さずに**作る・数える。

作るもの（すべて work/ 配下。data/raw には書かない）:
  - work/tokenized.dataset            A2 の FIB / aStr / dFIB を raw.X からトークン化したもの
  - work/genes_011.json               候補 8・陽性対照 71・恒常発現 3・中立 7・対照 500（入れ子 300）
  - work/donor_groups_011.json        出発点の提供者 50 人を 3 群に分けた結果
  - work/state_rows-seed{81,82,83}.json  状態の重心に使う行（3 群で同一の内容）
  - work/precheck_011.json            検査 C-1〜C-9 の結果

**乱数の鍵はすべてこのファイルに書いてある。結果を見る前に固定する。**

走らせ方:
    .venv/bin/python experiments/EXP-20261007-001/analysis/precheck_011.py
"""
from __future__ import annotations

import hashlib
import json
import pickle
import random
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261007-001"
W = E / "work"
sys.path.insert(0, str(ROOT / "scripts/geneformer"))
import run_isp as ri  # noqa: E402

# ---------------------------------------------------------------- 事前に固定する値
BLINDED = {"TP53", "CDKN1A", "LTBR", "SPP1", "LCN2"}           # 盲検。数えない・使わない
CANDIDATES = ["MADCAM1", "CXCL13", "CCL21", "CA12", "ZFHX3", "SLC4A9", "GSN", "ELN"]
HOUSEKEEPING = ["ACTB", "GAPDH", "B2M"]
NEUTRAL = ["CERS4", "TATDN1", "IFT25", "VAPA", "PHPT1", "ANAPC16", "EEF1D"]
PC76_CSV = ROOT / ("outputs/literature/2026-09-30-independent-goal-state/markers/"
                   "adaptive_stromal_state_aFIB_plus_MYOF_full_conserved.csv")
N_CTRL_EXPR, N_CTRL_ZERO = 400, 100          # 対照 500 = 発現あり 400 + 出発点で検出 0 の 100
N_NESTED_EXPR, N_NESTED_ZERO = 240, 60       # 入れ子 300（同じ比率 4:1）
EXPR_MIN_FRAC = 0.01                         # 「発現あり」= 出発点 252 細胞の 1% 以上（3 細胞以上）で入力に載る
KEY_CTRL = 20261007                          # 対照を引く乱数の鍵
KEY_GROUP = 20261008                         # 提供者の群分けの同順位を崩す乱数の鍵（2026-10-08 改訂。疾患で層別）
STRATA_ORDER = ["acute kidney injury", "normal", "chronic kidney disease"]  # 提供者数の少ない層から
KEY_STATE = "20261007-state"                 # 重心の細胞を引く乱数の鍵
CAP = 20                                     # 提供者あたりの上限（摂動・重心とも）
SEEDS = [81, 82, 83]                         # 3 群の名前（run_isp.py の seed として使う）
STATES = {"start": "FIB", "goal": "aStr", "alt": ["dFIB"]}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def tokenize() -> Path:
    cfg = {"data": {"input_h5ad": str(W / "input.h5ad"), "tokenized_dataset": str(W / "tokenized.dataset"),
                    "custom_attr_name_dict": {"isp_state": "isp_state", "donor_id": "donor_id",
                                              "author_cell_type": "author_cell_type",
                                              "disease": "disease", "n_counts": "n_counts"}},
           "perturbation": {"nproc": 1}, "model": {"version": "V2"}}
    ri.single_process_datasets(cfg)
    return ri.tokenize(cfg, W)


def main() -> None:
    import anndata as ad
    from datasets import load_from_disk

    out: dict = {"実験": "EXP-20261007-001", "注": "摂動は 1 回も回していない"}
    tok = tokenize()
    ds = load_from_disk(str(tok))
    st = np.array(ds["isp_state"])
    don = np.array([str(x) for x in ds["donor_id"]])
    dis = np.array([str(x) for x in ds["disease"]])
    ln = np.array(ds["length"])
    ids = [list(x) for x in ds["input_ids"]]
    out["C-0 トークン化"] = {"細胞数": int(len(ds)), "状態ごと": dict(Counter(st.tolist())),
                            "入力の長さの中央（状態ごと）": {s: float(np.median(ln[st == s])) for s in np.unique(st)}}

    # ------------------------------------------------------------ 辞書と記号の対応
    from geneformer import TOKEN_DICTIONARY_FILE
    tokd = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
    A = ad.read_h5ad(str(W / "input.h5ad"), backed="r")
    sym = A.var["gene_name"].astype(str).to_numpy()
    ens = A.var_names.astype(str).to_numpy()
    s2e: dict[str, str] = {}
    for s_, e_ in zip(sym, ens):
        s2e.setdefault(s_, e_)
    e2s = {e_: s_ for s_, e_ in zip(sym, ens)}
    inv = {v: k for k, v in tokd.items()}

    def to_ens(symbols: list[str]) -> dict[str, str | None]:
        return {s_: (s2e.get(s_) if s2e.get(s_) in tokd else None) for s_ in symbols}

    blinded_ens = {s2e[s_] for s_ in BLINDED if s_ in s2e}
    out["C-1 盲検の Ensembl ID（除外にだけ使う）"] = sorted(blinded_ens)

    pc76 = [l.strip() for l in PC76_CSV.read_text().splitlines()[1:] if l.strip()]
    pc_var = [g for g in pc76 if g in s2e]
    pc71 = [g for g in pc_var if s2e[g] in tokd]
    cand = to_ens(CANDIDATES)
    hk = to_ens(HOUSEKEEPING)
    neu = to_ens(NEUTRAL)
    miss = [k for d in (cand, hk, neu) for k, v in d.items() if v is None]
    if miss:
        raise SystemExit(f"辞書または A2 の var に無い: {miss}")
    if len(pc71) != 71:
        raise SystemExit(f"陽性対照が 71 個でない: {len(pc71)}")
    out["C-1 辞書"] = {"陽性対照 76 → var": len(pc_var), "→ 辞書": len(pc71),
                       "var に無い": sorted(set(pc76) - set(pc_var)),
                       "辞書に無い": sorted(set(pc_var) - set(pc71)),
                       "候補・恒常発現・中立": "すべて辞書にある"}

    # ------------------------------------------------------------ 出発点での検出（入力に載る細胞の割合）
    fib_rows = np.where(st == STATES["start"])[0]
    cnt = Counter()
    for i in fib_rows:
        cnt.update(set(int(t) for t in ids[i]))
    lineage_cnt = Counter()
    for i in np.where(st != STATES["start"])[0]:
        lineage_cnt.update(set(int(t) for t in ids[i]))
    n_fib = len(fib_rows)

    # ------------------------------------------------------------ 対照 500（入れ子 300）
    used = ({s2e[g] for g in pc71} | set(cand.values()) | set(hk.values()) | set(neu.values()) | blinded_ens)
    special = {"<pad>", "<mask>", "<cls>", "<eos>"}
    universe = sorted(e for e in tokd if e not in special and e in e2s and e not in used)
    expr = [e for e in universe if cnt.get(tokd[e], 0) >= int(np.ceil(EXPR_MIN_FRAC * n_fib))]
    zero = [e for e in universe if cnt.get(tokd[e], 0) == 0 and lineage_cnt.get(tokd[e], 0) >= 1]
    rng = random.Random(KEY_CTRL)
    pick_e = rng.sample(expr, N_CTRL_EXPR)
    pick_z = rng.sample(zero, N_CTRL_ZERO)
    controls = pick_e + pick_z
    nested = pick_e[:N_NESTED_EXPR] + pick_z[:N_NESTED_ZERO]
    out["C-2 対照の引き方"] = {
        "母集団": "辞書にあり A2 の var にあり、候補・陽性対照・恒常発現・中立・盲検のどれでもない遺伝子",
        "母集団の数": len(universe),
        "発現あり（出発点 252 細胞の 1% 以上で入力に載る）の数": len(expr),
        "出発点で 0 だが目標か参照で 1 細胞以上に載るの数": len(zero),
        "引いた数": {"発現あり": N_CTRL_EXPR, "出発点で 0": N_CTRL_ZERO},
        "入れ子 300": {"発現あり": N_NESTED_EXPR, "出発点で 0": N_NESTED_ZERO, "作り方": "引いた順の先頭"},
        "乱数の鍵": KEY_CTRL,
        "出発点で 0 の対照の数（候補 MADCAM1・CXCL13 と同じ条件）": N_CTRL_ZERO}

    # ------------------------------------------------------------ 重なりの総当たり（盲検を含む）
    sets = {"候補": set(cand.values()), "陽性対照71": {s2e[g] for g in pc71}, "恒常発現3": set(hk.values()),
            "中立7": set(neu.values()), "対照500": set(controls), "盲検5": blinded_ens}
    ov = {}
    names = list(sets)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            x = sorted(e2s.get(e, e) for e in sets[a] & sets[b])
            ov[f"{a} × {b}"] = x
    allowed = {"候補 × 陽性対照71": ["ELN", "GSN"]}
    bad = {k: v for k, v in ov.items() if v and allowed.get(k) != v}
    out["C-3 重なりの総当たり"] = {"結果": ov, "許される重なり": allowed, "それ以外の重なり": bad,
                                   "合格": not bad}
    if bad:
        raise SystemExit(f"想定外の重なり: {bad}")
    measured = sorted(set().union(*[v for k, v in sets.items() if k != "盲検5"]))
    if set(measured) & blinded_ens:
        raise SystemExit("盲検の遺伝子が測る遺伝子に入っている")
    out["C-3 測る遺伝子の数"] = len(measured)

    # ------------------------------------------------------------ 提供者の 3 群
    fib_donors = don[fib_rows]
    capped = Counter({d: min(c, CAP) for d, c in Counter(fib_donors.tolist()).items()})
    dis_of = {d: Counter(dis[fib_rows][fib_donors == d].tolist()).most_common(1)[0][0] for d in capped}
    rg = random.Random(KEY_GROUP)
    groups = {s: [] for s in SEEDS}
    tot = {s: 0 for s in SEEDS}
    strata = sorted(set(dis_of.values()))
    if sorted(STRATA_ORDER) != strata:
        raise SystemExit(f"疾患の区分が想定と違う: {strata}")
    for stratum in STRATA_ORDER:                      # 提供者数の少ない層から
        order = sorted(d for d in capped if dis_of[d] == stratum)
        rg.shuffle(order)                             # 同順位を崩すための並べ替え
        order.sort(key=lambda d: -capped[d])          # 安定ソート: 細胞数の多い提供者から
        n_in = {s: 0 for s in SEEDS}
        for d in order:   # (1) この層の提供者が最も少ない群 (2) 細胞数の合計が最も少ない群 (3) 番号の小さい群
            s = min(SEEDS, key=lambda x: (n_in[x], tot[x], x))
            groups[s].append(d)
            n_in[s] += 1
            tot[s] += capped[d]
    gsum = {}
    for s in SEEDS:
        gsum[str(s)] = {"提供者数": len(groups[s]), "細胞数（上限 20 の後）": tot[s],
                        "提供者あたり最大": max(capped[d] for d in groups[s]),
                        "疾患の内訳（提供者数）": dict(Counter(dis_of[d] for d in groups[s]))}
    allg = [d for s in SEEDS for d in groups[s]]
    assert len(allg) == len(set(allg)) == len(capped), "群が重なっているか欠けている"
    (W / "donor_groups_011.json").write_text(json.dumps(
        {"規則": "疾患で層別（急性腎障害 → 正常 → 慢性腎臓病の順に処理）。層の中では上限 20 を当てた細胞数の多い順に、"
                 "(1) その層の提供者が最も少ない群 (2) 細胞数の合計が最も少ない群 (3) 番号の小さい群 へ入れる。"
                 "同じ細胞数の提供者の順は random.Random(20261008).shuffle（層の順に 1 本の乱数を使う）で決める",
         "groups": {str(s): sorted(groups[s]) for s in SEEDS}, "要約": gsum},
        ensure_ascii=False, indent=1) + "\n")

    # ------------------------------------------------------------ 各群で実際に摂動される細胞（run_isp.py と同じ手順を再現）
    cfg_sel = {"cell_selection": {"mode": "common", "donor_key": "donor_id",
                                  "donor_groups": {str(s): groups[s] for s in SEEDS}},
               "states": {"start_state": STATES["start"]}, "data": {"state_key": "isp_state"}}
    m, rows, donors = ri.start_matrix(tok, cfg_sel)
    pools, c4 = {}, {}
    for s in SEEDS:
        gset = set(groups[s])
        idx = [i for i in range(m.shape[0]) if donors[i] in gset]
        keep = ri.balanced_subset(idx, donors, 200, CAP, random.Random(f"{s}-common"))
        pools[s] = [int(rows[i]) for i in keep]
        L = ln[pools[s]]
        c4[str(s)] = {"細胞数": len(keep), "提供者数": len({donors[i] for i in keep}),
                      "入力の長さの中央": float(np.median(L)),
                      "入力の長さの第1〜第3四分位": [float(np.percentile(L, 25)), float(np.percentile(L, 75))]}
    inter = [len(set(pools[a]) & set(pools[b])) for a in SEEDS for b in SEEDS if a < b]
    out["C-4 3 群の細胞"] = {"群ごと": c4, "群どうしで重なる細胞": inter, "合格": all(x == 0 for x in inter),
                             "出発点全体の入力の長さの中央": float(np.median(ln[fib_rows]))}
    out["C-4 群の要約"] = gsum

    # ------------------------------------------------------------ 候補が各群で元から入力に載っている細胞の数
    c5 = {}
    for g, e in cand.items():
        t = tokd[e]
        c5[g] = {str(s): int(sum(t in set(ids[i]) for i in pools[s])) for s in SEEDS}
    out["C-5 候補が元から入力に載っている細胞（群ごと）"] = c5
    # 対照のうち、各群で元から載っている細胞が 0 の遺伝子の数（低い側の較正）
    zero_in_all = sum(all(tokd[e] not in set().union(*[set(ids[i]) for i in pools[s]]) for s in SEEDS)
                      for e in controls)
    out["C-5 対照のうち 3 群すべてで元から載っている細胞が 0 の遺伝子"] = int(zero_in_all)

    # ------------------------------------------------------------ 状態の重心に使う細胞（3 群で同一）
    rs = random.Random(KEY_STATE)
    srows = {}
    for name in [STATES["start"], STATES["goal"]] + STATES["alt"]:
        idx = [int(i) for i in np.where(st == name)[0]]
        srows[name] = ri.balanced_subset(idx, don.tolist(), None, CAP, rs)
    spec = {"rows": srows,
            "注": "提供者あたり上限 20。random.Random('20261007-state') で引いた。3 群で同一の内容（軸を群によらず同じにするため）"}
    for s in SEEDS:
        (W / f"state_rows-seed{s}.json").write_text(json.dumps(spec, ensure_ascii=False) + "\n")
    out["C-6 状態の重心の細胞"] = {k: {"細胞数": len(v), "入力の長さの中央": float(np.median(ln[v])),
                                       "提供者数": len(set(don[v].tolist()))} for k, v in srows.items()}
    out["C-6 state_rows の sha256"] = {str(s): sha(W / f"state_rows-seed{s}.json") for s in SEEDS}

    # ------------------------------------------------------------ 中立対照パネルの A2 での差（出発点 対 目標）
    goal_rows = np.where(st == STATES["goal"])[0]
    gcnt = Counter()
    for i in goal_rows:
        gcnt.update(set(int(t) for t in ids[i]))
    c7 = {}
    for g, e in {**hk, **neu}.items():
        a_, b_ = cnt.get(tokd[e], 0) / n_fib, gcnt.get(tokd[e], 0) / len(goal_rows)
        c7[g] = {"出発点の検出率": round(a_, 4), "目標の検出率": round(b_, 4),
                 "比（目標÷出発点）": round(b_ / a_, 3) if a_ else None}
    out["C-7 恒常発現・中立対照の検出率（出発点 対 目標）"] = c7

    genes = {"候補": cand, "陽性対照71": {g: s2e[g] for g in pc71}, "恒常発現3": hk, "中立7": neu,
             "対照500": [{"ensembl": e, "symbol": e2s[e], "層": "発現あり" if e in set(pick_e) else "出発点で0",
                         "出発点で入力に載る細胞": int(cnt.get(tokd[e], 0))} for e in controls],
             "対照300（入れ子）": nested, "測る遺伝子": measured,
             "盲検（除外にだけ使う）": sorted(blinded_ens)}
    (W / "genes_011.json").write_text(json.dumps(genes, ensure_ascii=False, indent=1) + "\n")
    (W / "precheck_011.json").write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps(out, ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    main()

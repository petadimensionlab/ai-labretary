#!/usr/bin/env python3
"""EXP-20261004-001 着手前チェック — **摂動を 1 回も回さずに**計算できるものだけを出す。

この実験の主目的は「遺伝子ごとに摂動した細胞の顔ぶれが違う」という弱点を直すことなので、
着手前チェックの中心は**共通の細胞集合が本当に作れるか・本当に共通か**である。

出力: experiments/EXP-20261004-001/work/precheck_009.json

| # | 見るもの | 不合格の条件 |
|---|---|---|
| C-1 | 測る遺伝子が Geneformer の辞書にあるか | 1 つでも無い |
| C-2 | 共通の細胞集合が 3 seed とも取れるか（提供者上限 20・合計 200） | 200 未満の seed がある |
| C-3 | 共通の細胞集合が**遺伝子によらず同一**か（設計上の保証の検算） | 1 つでも違う |
| C-4 | 盲検中の遺伝子・判定に使うセットどうしの重なり | 盲検が混ざる |
| C-5 | 共通細胞の中で各遺伝子が「すでに入力に載っている」割合 | （記述。較正の範囲を見る） |
| C-6 | 低発現側の対照が候補の範囲を覆うか | 覆えない |
| C-7 | 深さ（入力の長さ）が状態の重心とそろっているか | （記述） |

使い方: `.venv/bin/python experiments/EXP-20261004-001/analysis/precheck_009.py`
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261004-001"
sys.path.insert(0, str(ROOT / "scripts/geneformer"))

SEEDS = [81, 82, 83]
N_CELLS = 200
PER_DONOR = 20
# 盲検中（config.blocked_genes と同じ。LTBR・CDKN1A は 2026-09-27 に解除済みだが、
# この実験では使わないので対象外として扱う）
BLOCKED = {"ENSG00000141510": "TP53", "ENSG00000118785": "SPP1", "ENSG00000148346": "LCN2"}

# 上位候補（探索 EXP-20261001-001 + EXP-20261002-002 の 1,223 遺伝子の順位表から、
# 3 seed すべてで z > 2.73 を満たしたもののうち、PI が指定した 8 個）。
# **Ensembl ID は手で書かない。** 既存の遺伝子表から引いて、引けなければ止める
# （guardrails 2026-10-03: 記録に書く値には出どころを添える）。
CANDIDATE_SYMBOLS = ["MADCAM1", "CA12", "CXCL13", "CCL21", "ZFHX3", "SLC4A9", "GSN", "ELN"]


def load_gene_tables():
    g = json.loads((ROOT / "experiments/EXP-20261001-001/work/genes_005.json").read_text())["genes"]
    x = json.loads((ROOT / "experiments/EXP-20261002-002/work/extra_null_007.json").read_text())["追加する対照"]
    subs = json.loads((ROOT / "experiments/EXP-20261001-001/work/subsets_005.json").read_text())["subsets"]
    npool = json.loads((ROOT / "experiments/EXP-20261001-001/work/null_pool_005.json").read_text())
    return g, x, subs, npool


def main() -> None:
    import run_isp
    cfg = yaml.safe_load((E / "work/config-precheck.yaml").read_text())
    tok = ROOT / cfg["data"]["tokenized_dataset"]
    out: dict = {
        "実験": "EXP-20261004-001",
        "注": "摂動は 1 回も回していない。トークン化済み入力と遺伝子リストだけから計算した。",
        "トークン化済み入力": str(tok.relative_to(ROOT)),
    }

    g, x, subs, npool = load_gene_tables()
    allg = g + x
    ens = {r["symbol"]: r["ensembl"] for r in allg}
    sym = {r["ensembl"]: r["symbol"] for r in allg}
    n_start_tbl = {r["ensembl"]: r["n_start"] for r in allg}

    pos71 = [d["ensembl"] for d in subs["positive_control71"]]
    hk = [d["ensembl"] for d in subs["housekeeping"]]           # ACTB / GAPDH / B2M
    neutral = [d["ensembl"] for d in subs["neutral_adopted"]]   # ACTB / CERS4 / ANAPC16 / EEF1D
    hk_all = list(dict.fromkeys(hk + neutral))
    b500 = [ens[d["symbol"]] for d in npool["versions"]["B_500"]]
    miss_sym = [s for s in CANDIDATE_SYMBOLS if s not in ens]
    if miss_sym:
        raise SystemExit(f"候補の Ensembl ID を既存の遺伝子表から引けません: {miss_sym}")
    CANDIDATES = {s: ens[s] for s in CANDIDATE_SYMBOLS}
    cand = list(CANDIDATES.values())

    # ---- 対照プール: B_500 から候補を除く（候補が自分の属する対照プールに入っていてはいけない）
    null_main = [t for t in b500 if t not in set(cand)]
    out["対照プール"] = {
        "土台": "EXP-20261001-001 の null_pool_005.json / versions.B_500（入れ子: A_300 ⊂ B_500）",
        "B_500 の個数": len(b500),
        "候補として抜いたもの": [sym[t] for t in b500 if t in set(cand)],
        "主対照の個数": len(null_main),
    }

    # ---- start_matrix（摂動なし）
    m, rows, donors = run_isp.start_matrix(tok, cfg)
    tokdict = run_isp.gene_token_dict()
    out["出発点"] = {"状態": cfg["states"]["start_state"], "細胞数": int(m.shape[0]),
                     "提供者数": len(set(donors))}

    # ---- C-1 辞書
    def tid(e):
        t = tokdict.get(e)
        return int(t) if t is not None and t < m.shape[1] else None
    measured = list(dict.fromkeys(cand + pos71 + hk_all + null_main))
    missing = [sym.get(e, e) for e in measured if tid(e) is None]
    out["C-1 辞書"] = {"測る遺伝子": len(measured), "辞書に無いもの": missing, "pass": not missing}

    # ---- C-2 / C-3 共通の細胞集合
    sel = cfg["cell_selection"]
    assert sel["mode"] == "common" and sel["max_cells_per_gene"] == N_CELLS and sel["max_cells_per_donor"] == PER_DONOR
    lengths_all = np.diff(m.indptr) if m.format == "csr" else None
    common: dict[int, list[int]] = {}
    c2 = []
    for s in SEEDS:
        rng = random.Random(f"{s}-common")
        keep = run_isp.balanced_subset(list(range(m.shape[0])), donors, N_CELLS, PER_DONOR, rng)
        common[s] = [int(rows[i]) for i in keep]
        ds = [donors[i] for i in keep]
        c2.append({"seed": s, "細胞数": len(keep), "提供者数": len(set(ds)),
                   "提供者あたり最大": max(ds.count(d) for d in set(ds))})
    out["C-2 共通の細胞集合"] = {"seed ごと": c2,
                                 "pass": all(r["細胞数"] == N_CELLS and r["提供者あたり最大"] <= PER_DONOR for r in c2)}

    # C-3: 遺伝子を変えても同じ集合になることの検算（common モードは鍵に遺伝子を含めないので設計上同一）
    probe = [cand[0], pos71[0], null_main[0], hk_all[0]]
    same = True
    for s in SEEDS:
        for e in probe:
            rng = random.Random(f"{s}-common")   # select_cells と同じ鍵
            k = [int(rows[i]) for i in run_isp.balanced_subset(list(range(m.shape[0])), donors, N_CELLS, PER_DONOR, rng)]
            same &= (k == common[s])
    out["C-3 遺伝子によらず同一か"] = {
        "確かめ方": "select_cells と同じ手順を 4 つの役割の遺伝子で再現し、行番号の並びが完全一致するかを見る",
        "pass": bool(same),
        "注": "common モードは乱数の鍵を f'{seed}-common' とし遺伝子を含めないので設計上同一。"
              "本番では全遺伝子の cells.csv の cell_uid 列を突き合わせて機械的に確かめる（A-2）",
    }

    # ---- C-4 重なり
    sets = {"候補8": cand, "陽性対照71": pos71, "恒常発現・中立": hk_all, "主対照": null_main}
    ov = {}
    names = list(sets)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            inter = sorted(set(sets[a]) & set(sets[b]))
            if inter:
                ov[f"{a} ∩ {b}"] = [sym.get(t, t) for t in inter]
    blocked_hit = [BLOCKED[t] for t in measured if t in BLOCKED]
    out["C-4 重なり"] = {"セットどうしの重なり": ov, "盲検の混入": blocked_hit,
                         "pass": not blocked_hit}

    # ---- C-5 共通細胞の中での「すでに載っている」割合
    start_index = {int(r): i for i, r in enumerate(rows)}
    common_start_idx = {s: [start_index[r] for r in common[s]] for s in SEEDS}
    def frac(e):
        t = tid(e)
        if t is None:
            return None
        col = m[:, t].nonzero()[0]
        colset = set(int(v) for v in col)
        return {str(s): sum(1 for i in common_start_idx[s] if i in colset) for s in SEEDS}, int(len(col))

    rowsC5 = []
    for name, e in CANDIDATES.items():
        f, n = frac(e)
        rowsC5.append({"symbol": name, "ensembl": e, "出発点で入力に載る細胞": n,
                       "出発点に占める割合": round(n / m.shape[0], 4),
                       "共通200細胞のうち載っている数(seed81/82/83)": [f[str(s)] for s in SEEDS]})
    out["C-5 共通細胞の中で既に載っている割合（候補）"] = rowsC5

    ns_null = np.array([n_start_tbl[t] for t in null_main if t in n_start_tbl])
    ns_cand = np.array([n_start_tbl.get(t, np.nan) for t in cand], dtype=float)
    out["C-5 較正の範囲"] = {
        "主対照の出発点発現細胞数": {"最小": int(ns_null.min()), "1 パーセンタイル": float(np.percentile(ns_null, 1)),
                                     "中央": float(np.median(ns_null)), "最大": int(ns_null.max())},
        "候補の出発点発現細胞数": {sym.get(t, t): (None if np.isnan(v) else int(v))
                                   for t, v in zip(cand, ns_cand)},
        "主対照より少ない候補": [sym.get(t, t) for t, v in zip(cand, ns_cand) if v < ns_null.min()],
    }

    # ---- C-6 低発現の追加対照（候補の範囲を覆う）
    det = {k: int(np.asarray(m.sum(axis=0)).ravel()[v])
           for k, v in tokdict.items() if str(k).startswith("ENS") and v < m.shape[1]}
    judged = set(measured) | set(BLOCKED) | {d["ensembl"] for v in subs.values() for d in v} | \
             {r["ensembl"] for r in allg}
    lo, hi = int(cfg["controls"]["low_expression_controls"]["n_start_min"]), \
             int(cfg["controls"]["low_expression_controls"]["n_start_max"])
    n_extra = int(cfg["controls"]["low_expression_controls"]["n"])
    pool = sorted([e for e, c in det.items() if lo <= c <= hi and e not in judged])
    rng = random.Random(cfg["controls"]["low_expression_controls"]["seed"])
    extra = sorted(rng.sample(pool, n_extra)) if len(pool) >= n_extra else []
    out["C-6 低発現の追加対照"] = {
        "条件": f"出発点で入力に載る細胞が {lo}〜{hi}、既存のどの判定セットにも入っていない",
        "候補になりうる遺伝子": len(pool), "採る数": n_extra,
        "引き方": f"random.Random({cfg['controls']['low_expression_controls']['seed']}).sample",
        "採ったもの": extra,
        "発現細胞数": sorted(det[e] for e in extra)[:10] + ["..."] if extra else [],
        "pass": len(extra) == n_extra,
    }

    # ---- C-7 深さ
    L = np.diff(m.tocsr().indptr)
    out["C-7 深さ"] = {
        "出発点全体の入力の長さ 中央": float(np.median(L)),
        "共通細胞の入力の長さ 中央": {str(s): float(np.median(L[common_start_idx[s]])) for s in SEEDS},
        "注": "状態の重心は EXP-20261001-001 と同じ state_rows-seed*.json（深さをそろえた行）を使う",
    }

    # ---- 測る遺伝子の最終一覧
    final = list(dict.fromkeys(cand + pos71 + hk_all + null_main + extra))
    out["測る遺伝子の最終一覧"] = {
        "合計": len(final),
        "内訳": {"候補": len(cand), "陽性対照": len(pos71), "恒常発現・中立": len(hk_all),
                 "主対照(B_500 − 候補)": len(null_main), "低発現の追加対照": len(extra),
                 "重複を除いた合計": len(final)},
    }
    (E / "work").mkdir(parents=True, exist_ok=True)
    (E / "work/genes_009.json").write_text(json.dumps(
        {"候補": {k: v for k, v in CANDIDATES.items()},
         "陽性対照71": pos71, "恒常発現・中立": hk_all,
         "主対照": null_main, "低発現の追加対照": extra,
         "測る遺伝子（この順に回す）": final}, ensure_ascii=False, indent=1) + "\n")

    fails = [k for k, v in out.items() if isinstance(v, dict) and v.get("pass") is False]
    out["総合"] = {"不合格の項目": fails, "pass": not fails}
    (E / "work/precheck_009.json").write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps(out, ensure_ascii=False, indent=1, default=float))
    sys.exit(0 if not fails else 2)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""介入する遺伝子セットの数え上げ（procedure.md 2026-10-01T19:35 固定版をそのまま実装）。

**摂動は 1 回も回さない。** 読むのは遺伝子リストと、既に作ってある遺伝子ごとの集計だけ。
すべて読み取りのみ。書き込むのは --out の JSON だけ。

出力: outputs/design/2026-10-01-intervention-gene-set/counts.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E30 = ROOT / "experiments/EXP-20260930-001"
E28 = ROOT / "experiments/EXP-20260928-002"
OUT = ROOT / "outputs/design/2026-10-01-intervention-gene-set"

MIN_CAP20 = 120           # procedure §2 / F-1
K_MIN = 30                # procedure §4 条件 K
FLOOR_MAX = 1.29          # procedure §6 合格の条件
BLIND = {"TP53", "SPP1", "LCN2"}   # procedure §4 F-3（LTBR・CDKN1A は解除済み）

# ---- 検出下限（power_004.py と同一の規則・同一の対照分布・同一の τ）
N_DRAW, PCT, TAU, POWER_TARGET, SEED = 200_000, 0.99, 2.051, 0.80, 20260930


def znull() -> np.ndarray:
    g = pd.read_csv(ROOT / "experiments/EXP-20260929-001/results/eval/gene_seed_003.csv")
    g = g[(g.col == "m_goal") & (g.role == "negative")]
    seeds = sorted(g.seed.unique())
    p = g.pivot_table(index="gene", columns="seed", values="D_raw")
    z = np.column_stack([(p[s] - p[s].mean()) / p[s].std(ddof=1) for s in seeds])
    return z[~np.isnan(z).any(axis=1)]


class Power:
    def __init__(self) -> None:
        self.z = znull()
        self.rng = np.random.default_rng(SEED)

    def floor(self, k: int, abs_min: float) -> dict:
        z, rng = self.z, self.rng
        N = z.shape[0]
        thr = np.percentile(z[rng.integers(0, N, size=(N_DRAW, k))].mean(axis=1), PCT * 100, axis=0)
        thr = np.atleast_1d(thr)
        if thr.size == 1:                     # per-seed にする
            thr = np.array([np.percentile(z[rng.integers(0, N, size=(N_DRAW, k)), s].mean(axis=1),
                                          PCT * 100) for s in range(z.shape[1])])

        def power(delta: float, n: int = 40_000) -> float:
            idx = rng.integers(0, N, size=(n, k))
            eff = rng.normal(delta, TAU, size=(n, k))
            obs = z[idx].mean(axis=1) + eff.mean(axis=1)[:, None]
            return float(((obs > thr[None, :]).all(axis=1) & (obs >= abs_min).all(axis=1)).mean())

        lo, hi = 0.0, 6.0
        for _ in range(28):
            mid = (lo + hi) / 2
            if power(mid) >= POWER_TARGET:
                hi = mid
            else:
                lo = mid
        return {"K": k, "set_abs_min": abs_min,
                "percentile_99": [round(float(x), 4) for x in thr],
                "detection_floor": round((lo + hi) / 2, 3)}


def abs_min_for(k: int) -> float:
    return round(4.5 / np.sqrt(k), 2)


def main() -> None:
    # ---- 既存のセット（EXP-20260930-001。読み取りのみ）
    gj = json.loads((E30 / "work/genes_004.json").read_text())
    by_role: dict[str, set] = {}
    sym_of = {}
    for g in gj["genes"]:
        by_role.setdefault(g["role"], set()).add(g["symbol"])
        sym_of[g["ensembl"]] = g["symbol"]
    npool = json.loads((E30 / "work/null_pool_004.json").read_text())["versions"]
    null500 = {x["symbol"] for x in npool["B_500"]}
    null300 = {x["symbol"] for x in npool["A_300"]}

    existing = {
        "陽性対照 71": by_role["positive_control"],
        "目標の定義（著者注釈なので遺伝子なし）": set(),
        "対照遺伝子 500": null500,
        "対照遺伝子 300（入れ子）": null300,
        "恒常発現 3": by_role["housekeeping"],
        "恒常発現の追加 1": by_role["housekeeping_extra"],
        "中立対照パネル 7": by_role["negctrl_panel"],
        "前回の摂動対象 T1 11": by_role["target_T1"],
        "前回の摂動対象 T2 7": by_role["target_T2"],
    }

    # ---- 新しい出発点での遺伝子ごとの集計（読み取りのみ）
    st = pd.read_csv(E30 / "work/gene_stats_004.csv")
    st = st[st.symbol.notna()]
    stat = {r.symbol: r for r in st.itertuples()}

    # ---- 候補 (P1) 経路セット 120（EXP-20260928-002 で固定）
    p1 = json.loads((E28 / "work/pathway_set_002.json").read_text())["genes"]
    p1_rows = [{"symbol": g["symbol"], "tier": g["tier"], "group": g["group"], "note": g["note"]}
               for g in p1]

    # ---- 候補 (P2) 候補表の全体
    cand = pd.read_csv(E28 / "analysis/pathway_candidates_002.tsv", sep="\t", comment="#")
    p2_rows = [{"symbol": r.symbol, "tier": r.tier, "group": r.group, "note": r.note}
               for r in cand.itertuples()]

    def screen(rows: list[dict], label: str) -> dict:
        seen, uniq = set(), []
        for r in rows:
            if r["symbol"] not in seen:
                seen.add(r["symbol"])
                uniq.append(r)
        steps, kept = [], uniq
        drops: dict[str, list[str]] = {}

        def drop(name: str, pred) -> None:
            nonlocal kept
            out = [r for r in kept if pred(r["symbol"])]
            kept = [r for r in kept if not pred(r["symbol"])]
            drops[name] = sorted(x["symbol"] for x in out)
            steps.append({"手順": name, "落ちた数": len(out), "残り": len(kept)})

        drop("1. 盲検中の 3 個", lambda s: s in BLIND)
        drop("2. 陽性対照 71 と重なる", lambda s: s in existing["陽性対照 71"])
        drop("3. 目標の定義と重なる", lambda s: s in existing["目標の定義（著者注釈なので遺伝子なし）"])
        drop("4. 恒常発現・中立対照パネルと重なる",
             lambda s: s in (existing["恒常発現 3"] | existing["恒常発現の追加 1"]
                             | existing["中立対照パネル 7"]))
        drop("5. gene_stats に無い（var / トークン辞書に無い）", lambda s: s not in stat)
        drop("6. 出発点で上限 20 後 120 細胞未満",
             lambda s: int(stat[s].n_cap20) < MIN_CAP20)

        det = [{"symbol": r["symbol"],
                "n_start": int(stat[r["symbol"]].n_start),
                "n_donor_start": int(stat[r["symbol"]].n_donor_start),
                "n_cap20": int(stat[r["symbol"]].n_cap20),
                "det_start": round(float(stat[r["symbol"]].det_start), 4),
                "det_goal": round(float(stat[r["symbol"]].det_goal), 4),
                "tier": r["tier"], "group": r["group"]} for r in kept]
        return {"候補": label, "重複を除いた候補数": len(uniq), "除外の手順": steps,
                "落ちた遺伝子": drops, "K": len(kept), "採用した遺伝子": det}

    res: dict = {"tool": "outputs/design/2026-10-01-intervention-gene-set/count_sets.py",
                 "note": "摂動は 1 回も回していない。遺伝子リストと事前集計だけから計算した。",
                 "条件": {"K_min": K_MIN, "n_cap20_min": MIN_CAP20,
                          "検出下限の上限": FLOOR_MAX, "盲検": sorted(BLIND)},
                 "候補": {}}
    res["候補"]["P1_経路セット120"] = screen(p1_rows, "P1 経路セット 120")
    res["候補"]["P2_候補表の全体"] = screen(p2_rows, "P2 候補表の全体")

    # ---- 主判定の決定（procedure §3 の順番）
    pw = Power()
    chosen = None
    for key in ["P1_経路セット120", "P2_候補表の全体"]:
        k = res["候補"][key]["K"]
        am = abs_min_for(k) if k >= 2 else None
        f = pw.floor(k, am) if k >= 2 else None
        res["候補"][key]["set_abs_min"] = am
        res["候補"][key]["検出下限"] = f["detection_floor"] if f else None
        res["候補"][key]["対照の99パーセンタイル"] = f["percentile_99"] if f else None
        res["候補"][key]["条件K"] = "合格" if k >= K_MIN else "不合格"
        res["候補"][key]["条件_検出下限"] = ("合格" if f and f["detection_floor"] <= FLOOR_MAX
                                             else "不合格")
        if chosen is None and k >= K_MIN and f and f["detection_floor"] <= FLOOR_MAX:
            chosen = key
    res["主判定に採った候補"] = chosen

    if chosen:
        main_syms = [g["symbol"] for g in res["候補"][chosen]["採用した遺伝子"]]
        rows = res["候補"][chosen]["採用した遺伝子"]
        # ---- 入れ子の感度分析（procedure §5.2）
        subs = {
            "感度a_S1経路の構成因子": [r["symbol"] for r in rows if r["tier"] == "S1"],
            "感度b_S2経路の転写標的": [r["symbol"] for r in rows if r["tier"] == "S2"],
            "感度c_S3キュレーション": [r["symbol"] for r in rows if r["tier"] == "S3"],
            "感度d_経路特異的のみ": [r["symbol"] for r in rows if r["group"] == "経路特異的"],
            "感度e_主判定とT1T2の和": sorted(set(main_syms) | existing["前回の摂動対象 T1 11"]
                                             | existing["前回の摂動対象 T2 7"]),
        }
        res["入れ子の感度分析"] = {}
        for name, syms in subs.items():
            k = len(syms)
            am = abs_min_for(k) if k >= 2 else None
            f = pw.floor(k, am) if k >= 2 else None
            res["入れ子の感度分析"][name] = {
                "K": k, "set_abs_min": am,
                "検出下限": f["detection_floor"] if f else None,
                "検出下限が 1.29 以下か": bool(f and f["detection_floor"] <= FLOOR_MAX),
                "遺伝子": syms}

        # ---- F-2 総当たりの重なり
        sets = dict(existing)
        sets["【介入】主判定"] = set(main_syms)
        names = list(sets)
        pair = {}
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                ov = sorted(sets[a] & sets[b])
                pair[f"{a} × {b}"] = ov
        res["F2_総当たりの重なり"] = {
            "各セットの大きさ": {k: len(v) for k, v in sets.items()},
            "重なり": pair,
            "介入セットが関わる重なりの合計":
                sum(len(v) for k, v in pair.items() if "【介入】" in k),
            "全体の重なりの合計": sum(len(v) for v in pair.values()),
            "盲検中の 3 個の混入": sorted(set(main_syms) & BLIND),
        }
        # 検出率
        d_start = [r["det_start"] for r in rows]
        d_goal = [r["det_goal"] for r in rows]
        cap = [r["n_cap20"] for r in rows]
        res["主判定セットの実測"] = {
            "K": len(rows),
            "出発点での検出率": {"中央値": round(float(np.median(d_start)), 4),
                                 "最小": min(d_start), "最大": max(d_start)},
            "目標での検出率": {"中央値": round(float(np.median(d_goal)), 4),
                               "最小": min(d_goal), "最大": max(d_goal)},
            "上限20後の細胞数": {"中央値": int(np.median(cap)), "最小": min(cap), "最大": max(cap)},
            "tier の内訳": {t: sum(1 for r in rows if r["tier"] == t) for t in ["S1", "S2", "S3"]},
            "群の内訳": {g: sum(1 for r in rows if r["group"] == g)
                         for g in sorted({r["group"] for r in rows})},
        }

    (OUT / "counts.json").write_text(json.dumps(res, ensure_ascii=False, indent=2))
    print(json.dumps({k: v for k, v in res.items() if k != "候補"}, ensure_ascii=False)[:2000])
    for k, v in res["候補"].items():
        print(k, "K=", v["K"], "floor=", v.get("検出下限"), "abs_min=", v.get("set_abs_min"))
        for s in v["除外の手順"]:
            print("   ", s)


if __name__ == "__main__":
    main()

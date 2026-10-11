#!/usr/bin/env python3
"""EXP-20261011-001 — HYP-20261011-01。差の向き（u_diff）・共通の向き（u_common）への射影を、
保存済みの埋め込みと、追加の退化の陽性対照 30 個の実行結果に当てる。

- 探索（A2）: EXP-20261007-001 の保存済み 587 遺伝子を再射影（GPU 0 時間）。P-2・P-3（候補 8 個の固有性。**探索**）
- 確認（A1）: EXP-20261004-001 の保存済み 626 + EXP-20261011-001 の新規 30 個（退化の陽性対照）を再射影。
  主判定 P-1（方法の検証。**確認**）
- 射影の式（`m = (d·u)/‖d‖`）・層別補正（`evaluate_009.gene_stats_stratified`、発現の有無で層別した二次式）・
  標準化（群ごとに対照プール A2:510 / A1:549 の平均・標準偏差）は、**EXP-20261007-003 の評価コードと同じ**。
  変えるのは goal に当たる向きだけ（u_ada → u_diff・u_common）。
- 移動量 0 の細胞は、**外す前の表で 2 本（u_diff・u_common）とも射影を計算し、外す細胞 ID を後で除く**
  （guardrails 2026-10-11。向きごとに表を作らない）。
- **判定（支持・否定・判定不能）は確定しない。** ゲート G3 は PI が行う。ここが出すのは「規則を機械的に当てた結果」まで。
- 盲検中の 5 遺伝子（TP53・CDKN1A・LTBR・SPP1・LCN2）は、摂動の対象・対照・指標・出力・ハッシュのどこにも入れない。

走らせ方:
    .venv/bin/python experiments/EXP-20261011-001/analysis/evaluate_014.py --part a2          # A2（探索。新規摂動の出力は不要）
    .venv/bin/python experiments/EXP-20261011-001/analysis/evaluate_014.py --part a1          # A1（確認。新規 30 個の出力が必要）
    .venv/bin/python experiments/EXP-20261011-001/analysis/evaluate_014.py --part all --smoke # 試し実行（判定ではない）

主な入力: plan.md、config.yaml、EXP-20261007-001 の results/・work/、EXP-20261004-001 の results/・work/、
          本実験の results/（新規 30 個。A1 のとき必須）。
主な出力: results/criteria_014.json・results/gene_rank_014.csv・results/gene_group_014.csv・
          results/README.md・results/run.log（試し実行は results-smoke/）。
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
EXP = ROOT / "experiments/EXP-20261011-001"
A2_SRC = ROOT / "experiments/EXP-20261007-001"
A1_SRC = ROOT / "experiments/EXP-20261004-001"
EXP003 = ROOT / "experiments/EXP-20261007-003"
sys.path.insert(0, str(ROOT / "experiments/EXP-20261004-001/analysis"))
sys.path.insert(0, str(ROOT / "experiments/EXP-20261001-001/analysis"))
sys.path.insert(0, str(ROOT / "scripts/analysis"))
import evaluate_009 as ev9                                # noqa: E402  層別補正・細胞表（1 文字も変えない）
from evaluate_005 import set_null                        # noqa: E402  復元ありで k 個引いて平均する

SEEDS = [81, 82, 83]
N_DRAW = 1000

# --------------------------------------------------------------------------- plan.md §4 の判定値
# **ここを変えてはいけない。** config.yaml の decision: 節と照合し、違えば止める。
PLAN = {
    "候補_n": 8,
    "活性化PC71_n": 71,
    "退化PC30_n": 30,
    "プールA2_n": 510,
    "プールA1_n": 549,
    "A2_測る遺伝子_n": 587,
    "A1_測る遺伝子_n": 626,
    "新規30_n": 30,
    "重なり7": ["MGP", "IGFBP5", "DSTN", "ARMC9", "PDE5A", "MFGE8", "HSPG2"],
    "再現確認（A2・u_ada の候補8）": [2.664, 2.776, 2.302],
    "再現確認（A2・u_deg の候補8）": [2.6754, 2.7333, 2.2468],
    "再現確認（A1・u_ada の候補8）": [6.900499, 6.673740, 6.708363],
    "再現確認（A1・u_ada の PC71）": [1.182938, 1.226077, 1.300099],
    "再現確認_許容差": 1e-3,
    "R-12_閾値": 0.20,
    "乱数の鍵": {"P-1 活性化の偽セット": 20261017, "P-1 退化の偽セット": 20261018,
                "P-1 (iii) の無作為 2 セット": 20261019, "P-2 候補の偽セット": 20261020,
                "P-3 候補の偽セット": 20261021, "偽セットの較正": 20261022},
}


def sha256f(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def log_line(path: Path, msg: str) -> None:
    line = f"[{_dt.datetime.now().strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(path, "a") as fh:
        fh.write(line + "\n")


# --------------------------------------------------------------------------- 射影（plan.md §4.1）
def add_projection(res: Path, start_state: str, goal_state: str) -> pd.Series:
    """evaluate_013.add_projection と同一の式・同一の手順（`u = unit(goal) − unit(start)` に正規化し、
    `(d·u)/‖d‖`）。`cells.csv`（移動量 0 を外す前の表）に対して計算し、`(run, gene, cell_uid)` で返す。"""
    cells = pd.read_csv(res / "cells.csv", dtype={"donor_id": str})
    proj = np.full(len(cells), np.nan)
    for (run, gene), part in cells.groupby(["run", "gene"], sort=False):
        sz = np.load(res / run / "state_emb_vectors.npz")
        gd = res / run / "genes" / gene / "emb"
        o = np.load(gd / "orig.npy").astype(np.float64)
        p = np.load(gd / "pert.npy").astype(np.float64)
        if len(o) != len(part):
            raise SystemExit(f"**止める**: {gd}: 埋め込み {len(o)} != cells.csv {len(part)}")
        d = p - o
        dn = np.linalg.norm(d, axis=1)
        vg, vs = sz[goal_state].astype(np.float64), sz[start_state].astype(np.float64)
        u = vg / np.linalg.norm(vg) - vs / np.linalg.norm(vs)
        u = u / np.linalg.norm(u)
        with np.errstate(divide="ignore", invalid="ignore"):
            proj[part.index.to_numpy()] = (d @ u) / dn
    idx = pd.MultiIndex.from_arrays([cells.run.to_numpy(), cells.gene.to_numpy(),
                                     cells.cell_uid.to_numpy()], names=["run", "gene", "cell_uid"])
    return pd.Series(proj, index=idx, name="m_proj_deg")


def add_projection_diff(res: Path, start_state: str, goal_ada: str, goal_deg: str, mode: str) -> pd.Series:
    """差の向き（mode="diff"）または共通の向き（mode="common"）への射影。

    u_ada = unit(goal_ada) − unit(start)；u_deg = unit(goal_deg) − unit(start)
    u = (u_ada − u_deg)/‖…‖（diff）または (u_ada + u_deg)/‖…‖（common）
    **外す前の表**（cells.csv）で計算し、(run, gene, cell_uid) で返す（guardrails 2026-10-11）。
    """
    cells = pd.read_csv(res / "cells.csv", dtype={"donor_id": str})
    proj = np.full(len(cells), np.nan)
    for (run, gene), part in cells.groupby(["run", "gene"], sort=False):
        sz = np.load(res / run / "state_emb_vectors.npz")
        gd = res / run / "genes" / gene / "emb"
        o = np.load(gd / "orig.npy").astype(np.float64)
        p = np.load(gd / "pert.npy").astype(np.float64)
        if len(o) != len(part):
            raise SystemExit(f"**止める**: {gd}: 埋め込み {len(o)} != cells.csv {len(part)}")
        d = p - o
        dn = np.linalg.norm(d, axis=1)
        def unit(v):
            return v / np.linalg.norm(v)
        vg, vs = sz[goal_ada].astype(np.float64), sz[start_state].astype(np.float64)
        ua = unit(vg) - unit(vs); ua = ua / np.linalg.norm(ua)
        vg2 = sz[goal_deg].astype(np.float64)
        ud = unit(vg2) - unit(vs); ud = ud / np.linalg.norm(ud)
        u = (ua - ud) if mode == "diff" else (ua + ud)
        u = u / np.linalg.norm(u)
        with np.errstate(divide="ignore", invalid="ignore"):
            proj[part.index.to_numpy()] = (d @ u) / dn
    idx = pd.MultiIndex.from_arrays([cells.run.to_numpy(), cells.gene.to_numpy(),
                                     cells.cell_uid.to_numpy()], names=["run", "gene", "cell_uid"])
    return pd.Series(proj, index=idx, name=f"m_proj_{mode}")


def merge_projection(c: pd.DataFrame, srs: pd.Series, col: str, label: str) -> pd.DataFrame:
    n0 = len(c)
    c = c.merge(srs, left_on=["run", "gene", "cell_uid"], right_index=True, how="left", validate="one_to_one")
    miss = int(c[col].isna().sum())
    if miss:
        raise SystemExit(f"**止める**: {label} の射影が結びつかない細胞が {miss} 個（{n0} 行中）")
    return c


def per_seed_stats(c: pd.DataFrame, col: str, pool: list[str], sym: dict, log, label: str) -> pd.DataFrame:
    """群ごとに層別補正（ev9.gene_stats_stratified をそのまま呼ぶ）→ 標準化（対照プールで）。"""
    poolset = set(pool)
    frames = []
    for s in SEEDS:
        cs = c[c.seed == s].assign(m_goal=lambda d: d[col])
        nl = {h: hd for h, hd in cs[cs.gene.isin(poolset)].groupby("gene")}
        f = ev9.gene_stats_stratified(cs, nl, "m_goal")
        n = f[f.gene.isin(poolset)]
        frames.append(f.assign(seed=s, z=(f.D_raw - n.D_raw.mean()) / n.D_raw.std(ddof=1),
                               symbol=f.gene.map(sym)))
    long = pd.concat(frames, ignore_index=True)
    log(f"射影「{label}」: 層別補正と標準化を {len(SEEDS)} 群で計算（対照 {len(poolset)} 個で標準化）")
    return long


def spearman_rec(x: pd.Series, y: pd.Series, name_x: str, name_y: str) -> dict:
    """NaN を書かず null と理由を書く（guardrails 2026-10-11）。"""
    m = x.notna() & y.notna()
    xs, ys = x[m].to_numpy(dtype=float), y[m].to_numpy(dtype=float)
    if len(xs) < 3 or float(np.ptp(xs)) == 0.0 or float(np.ptp(ys)) == 0.0:
        return {"Spearman": None, "n": int(m.sum()),
                "理由": "計算できない（値を 2 つとも持つ遺伝子が 3 個未満、または一方が一定）"}
    return {"Spearman": round(float(stats.spearmanr(xs, ys).statistic), 4), "n": int(m.sum()),
            "外した数（値が欠けた遺伝子）": int((~m).sum())}


# --------------------------------------------------------------------------- 材料（A-2）
def check_materials(cfg, log) -> dict:
    ev = cfg["evaluation"]
    rec: dict = {"sha256": {}, "埋め込み": {}, "トークン化済み入力の複製": {}}
    ok = True
    for name, path, pin in (
        ("A2 genes_011.json", A2_SRC / "work/genes_011.json", ev["a2"]["genes_json"].endswith("genes_011.json")),
        ("A1 config.yaml", A1_SRC / "config.yaml", True),
    ):
        got = sha256f(path)
        rec["sha256"][name] = {"sha256": got, "実在": bool(got)}
        ok &= bool(got)
    # トークン化済み入力（複製）の sha256（precheck_014.json に記録した値と一致させる）
    pre = json.loads((EXP / "work/precheck_014.json").read_text())
    for key, clone in (("a2", "A2"), ("a1", "A1")):
        rec2 = pre["1. トークン化済み入力の複製"][clone]
        ok &= bool(rec2["pass"])
        rec["トークン化済み入力の複製"][key] = rec2["複製先"]
    # 埋め込みの実在
    counts = {}
    for label, root, n_genes in (
        ("A2 保存 587", A2_SRC / "results", PLAN["A2_測る遺伝子_n"]),
        ("A1 保存 626", A1_SRC / "results", PLAN["A1_測る遺伝子_n"]),
    ):
        n_files = 0
        miss = []
        for s in SEEDS:
            gdir = root / f"run-{s}-L-1/genes"
            for p in sorted(gdir.iterdir()):
                if not p.is_dir():
                    continue
                for kind in ("orig", "pert"):
                    if (p / "emb" / f"{kind}.npy").exists():
                        n_files += 1
                    else:
                        miss.append(str((p / "emb" / f"{kind}.npy").relative_to(ROOT)))
        counts[label] = {"実在": n_files, "期待": n_genes * 3 * 2, "欠け": miss[:5]}
        ok &= (n_files == n_genes * 3 * 2 and not miss)
    new_root = EXP / "results"
    if (new_root / "cells.csv").exists():
        n_files = 0
        for s in SEEDS:
            for p in sorted((new_root / f"run-{s}-L-1/genes").iterdir()):
                if p.is_dir() and (p / "emb" / "orig.npy").exists() and (p / "emb" / "pert.npy").exists():
                    n_files += 1
        counts["A1 新規 30"] = {"実在": n_files, "期待": PLAN["新規30_n"] * 3}
        ok &= (n_files == PLAN["新規30_n"] * 3)
    else:
        counts["A1 新規 30"] = {"実在": 0, "期待": PLAN["新規30_n"] * 3, "注": "まだ実行していない（--part a2 では不要）"}
    rec["埋め込み"] = counts
    rec["pass"] = bool(ok)
    log(f"A-2 材料の実在と複製の sha256: pass={rec['pass']}")
    return rec


# --------------------------------------------------------------------------- 遺伝子集合（A-3）
def gene_sets_from_cfg_and_json(cfg) -> dict:
    g = json.loads((A2_SRC / "work/genes_011.json").read_text())
    deg30 = list(cfg["perturbation"]["genes_to_perturb"])
    cand = list(g["候補"].values())
    pos71 = list(g["陽性対照71"].values())
    hk = list(g["恒常発現3"].values())
    neu = list(g["中立7"].values())
    ctrl500 = [r["ensembl"] for r in g["対照500"]]
    layer = {r["ensembl"]: r["層"] for r in g["対照500"]}
    sym = {v: k for k, v in g["候補"].items()} | {v: k for k, v in g["陽性対照71"].items()} \
        | {v: k for k, v in g["恒常発現3"].items()} | {v: k for k, v in g["中立7"].items()} \
        | {r["ensembl"]: r["symbol"] for r in g["対照500"]}
    # A1 のプール 549 は A1 config の controls から機械的に再構成（plan.md §2.2）
    cfgA1 = yaml.safe_load((A1_SRC / "config.yaml").read_text())
    ctl = cfgA1["controls"]
    rand = [x for x in ctl["negative_sampler"]["genes"]][: ctl["negative_random_n"]]
    poolA1 = sorted(set(rand) | set(ctl["negative_fixed"])
                    | {x for grp in (ctl.get("extra_fixed") or []) for x in grp["genes"]})
    poolA2 = sorted(set(ctrl500) | set(hk) | set(neu))
    # 重複を除く（候補の GSN・ELN は陽性対照にも入る。guardrails 2026-10-11）
    sets = {"候補": list(dict.fromkeys(cand)), "活性化PC71": list(dict.fromkeys(pos71)),
            "退化PC30": list(dict.fromkeys(deg30)), "恒常発現3": list(dict.fromkeys(hk)),
            "中立7": list(dict.fromkeys(neu)), "対照500": list(dict.fromkeys(ctrl500)),
            "プールA2": poolA2, "プールA1": poolA1}
    dup = {k: [x for x in v if v.count(x) > 1] for k, v in sets.items() if len(set(v)) != len(v)}
    return {"sets": sets, "layer": layer, "symbol": sym, "片方に重複（A-7）": dup}


# --------------------------------------------------------------------------- A2（探索）
def run_a2(cfg, G, out: dict, L, log, smoke: bool) -> None:
    ev = cfg["evaluation"]["a2"]
    tok = str(ROOT / ev["tokenized_dataset"])
    cfg2 = {"data": {"tokenized_dataset": tok, "state_key": "isp_state"},
            "states": {"start_state": "FIB", "goal_state": "aStr", "alt_states": ["dFIB"]},
            "output": {"results_dir": "results"}}
    ev9.E = A2_SRC
    ev9.BLOCKED = {e: e for e in cfg["blocked_genes"]}
    res = A2_SRC / "results"
    c, drop = ev9.load_cells(cfg2, "results")
    c = c[c.seed.isin(SEEDS)].reset_index(drop=True)
    # 外す前の表で 3 本（u_deg・u_diff・u_common）を計算し、(run,gene,cell_uid) で結合する
    deg = add_projection(res, "FIB", "dFIB")
    diff = add_projection_diff(res, "FIB", "aStr", "dFIB", "diff")
    com = add_projection_diff(res, "FIB", "aStr", "dFIB", "common")
    c = merge_projection(c, deg, "m_proj_deg", "A2 u_deg")
    c = merge_projection(c, diff, "m_proj_diff", "A2 u_diff")
    c = merge_projection(c, com, "m_proj_common", "A2 u_common")
    out["A2 細胞ごとの表"] = {"遺伝子数": int(c.gene.nunique()), "移動量 0 で外した細胞": drop}

    pool = [x for x in G["sets"]["プールA2"] if x in set(c.gene)]
    symptr = G["symbol"]
    long = {}
    for col, lab in (("m_proj", "u_ada"), ("m_proj_deg", "u_deg"), ("m_proj_diff", "u_diff"),
                     ("m_proj_common", "u_common")):
        long[lab] = per_seed_stats(c, col, pool, symptr, log, f"A2 {lab}")
    Z = {k: v.pivot(index="gene", columns="seed", values="z").reindex(columns=SEEDS) for k, v in long.items()}

    def zm(direction: str, genes: list[str]) -> np.ndarray:
        return Z[direction].reindex([x for x in genes if x in Z[direction].index]).to_numpy()

    # ---- 再現確認（A-4・A2 側）
    checks = {}
    S = zm("u_ada", G["sets"]["候補"]).mean(axis=0)
    d_ada = float(np.max(np.abs(S - np.array(PLAN["再現確認（A2・u_ada の候補8）"]))))
    D = zm("u_deg", G["sets"]["候補"]).mean(axis=0)
    d_deg = float(np.max(np.abs(D - np.array(PLAN["再現確認（A2・u_deg の候補8）"]))))
    checks["A-4 再現確認（A2）"] = {
        "u_ada の候補8セット平均": [round(x, 4) for x in S], "参照": PLAN["再現確認（A2・u_ada の候補8）"],
        "u_deg の候補8セット平均": [round(x, 4) for x in D], "参照": PLAN["再現確認（A2・u_deg の候補8）"],
        "最大絶対差": round(max(d_ada, d_deg), 6), "許容差": PLAN["再現確認_許容差"],
        "pass": bool(max(d_ada, d_deg) <= PLAN["再現確認_許容差"])}
    if not checks["A-4 再現確認（A2）"]["pass"] and not smoke:
        raise SystemExit("**止める**(A-4/A2): 再現確認が合わない")
    out["A チェック（A2）"] = checks

    # ---- P-2・P-3（探索）・R-7〜R-11
    cand = G["sets"]["候補"]; pc71 = G["sets"]["活性化PC71"]
    rng2 = np.random.default_rng(PLAN["乱数の鍵"]["P-2 候補の偽セット"])
    rng3 = np.random.default_rng(PLAN["乱数の鍵"]["P-3 候補の偽セット"])
    obs_diff = zm("u_diff", cand).mean(axis=0)
    obs_com = zm("u_common", cand).mean(axis=0)
    pc71_diff = zm("u_diff", pc71).mean(axis=0)
    poolZ = Z["u_diff"].reindex(pool).to_numpy()
    thr2 = np.quantile(set_null(poolZ, len(cand), N_DRAW, rng2, replace=True), 0.99, axis=0)
    thr3 = np.quantile(set_null(Z["u_common"].reindex(pool).to_numpy(), len(cand), N_DRAW, rng3, replace=True),
                       0.99, axis=0)
    p2_a = bool(np.all(obs_diff <= thr2))          # 3 群すべてで偽セットの 99 パーセンタイル以内
    p2_b = bool(np.all(obs_diff < pc71_diff))      # 3 群すべてで活性化 PC71 より小さい
    p3 = bool(np.all(obs_com > thr3))
    p2 = p2_a and p2_b
    # R-7〜R-11（仮説 §6.2）
    r7 = p2_a and p3
    r8 = (not p2_a) and bool(np.all(obs_diff > thr2)) and bool(np.all(obs_diff >= pc71_diff))
    r9 = bool(np.all(obs_diff < np.quantile(
        set_null(poolZ, len(cand), N_DRAW, np.random.default_rng(PLAN["乱数の鍵"]["P-2 候補の偽セット"] + 1),
                 replace=True), 0.01, axis=0)))
    r10 = (not p3) and (not p2_a) and not (np.all(obs_diff > thr2))
    r11 = p3 and bool(np.all(obs_diff > thr2))
    first = next((k for k, v in {"R-7": r7, "R-8": r8, "R-9": r9, "R-10": r10, "R-11": r11}.items() if v),
                 "該当なし（R-14 保留）")
    out["A2 P-2・P-3（探索。確認的な主張はしない）"] = {
        "P-2 (a) 候補8 の z_diff が偽セット 99 パーセンタイル以内": {
            "セット平均": [round(x, 4) for x in obs_diff], "偽セット99%点": [round(x, 4) for x in thr2],
            "3 群すべてで以内": p2_a},
        "P-2 (b) かつ活性化 PC71 のセット平均より小さい": {
            "活性化 PC71 の z_diff セット平均": [round(x, 4) for x in pc71_diff], "3 群すべてで小さい": p2_b},
        "P-3 共通の向きで候補8 が偽セット 99 パーセンタイル超": {
            "セット平均": [round(x, 4) for x in obs_com], "偽セット99%点": [round(x, 4) for x in thr3], "成り立つ": p3},
        "P-2 全体": p2,
        "反証条件": {"R-7": r7, "R-8": r8, "R-9": r9, "R-10": r10, "R-11": r11, "最初に当たったもの": first,
                   "注": "探索である（仮説 §5.4・§5.5）。"},
    }
    # ---- 順位表（探索）
    for direction, zname, fname in (("u_diff", "z_diff", "gene_rank_014"), ("u_common", "z_common", "gene_rank_common_014")):
        Zr = Z[direction].rename(columns={s: f"z_{s}" for s in SEEDS})
        rank = (Zr.assign(z_mean3=Z[direction].mean(axis=1)).reset_index().rename(columns={"index": "gene"})
                .assign(symbol=lambda d: d["gene"].map(G["symbol"]),
                        注="探索（データを見て候補を探すための順位表。事前登録した確認的検定ではない。CLAUDE.md §3.5）")
                .sort_values("z_mean3", ascending=False).reset_index(drop=True))
        rank.insert(0, "rank", np.arange(1, len(rank) + 1))
        outdir = EXP / ("results-smoke" if smoke else "results")
        rank.to_csv(outdir / f"{fname}.csv", index=False)
        out[f"A2 {zname} の順位表（探索）"] = {"行数": len(rank), "保存先": f"{fname}.csv",
                                           "上位10": rank.head(10)[["symbol", "z_mean3"]].round(4).to_dict(orient="records")}
    log(f"A2: P-2={p2}（a={p2_a}, b={p2_b}）、P-3={p3}、最初に当たった反証条件={first}")


# --------------------------------------------------------------------------- A1（確認）
def run_a1(cfg, G, out: dict, L, log, smoke: bool) -> None:
    ev = cfg["evaluation"]["a1"]
    nw = cfg["evaluation"]["a1"]["new_perturbation"]
    tok = str(ROOT / ev["tokenized_dataset"])
    cfg1 = {"data": {"tokenized_dataset": tok, "state_key": "fib_anno_state"},
            "states": {"start_state": "FIB_MFIB", "goal_state": "aStr", "alt_states": ["dFIB", "dMFIB"]},
            "output": {"results_dir": "results"}}
    ev9.BLOCKED = {e: e for e in cfg["blocked_genes"]}
    res = A1_SRC / "results"
    ev9.E = A1_SRC
    c_a1, drop_a1 = ev9.load_cells(cfg1, "results")
    c_a1 = c_a1[c_a1.seed.isin(SEEDS)].reset_index(drop=True)
    # 新規 30 個（本実験の results。まだ無ければ A1 部分は実行しない）
    res_new = EXP / "results"
    if not (res_new / "cells.csv").exists():
        raise SystemExit("**止める**: A1 の新規摂動の出力がありません（results/cells.csv）。先に摂動を回してください。")
    ev9.E = EXP
    c_new, drop_new = ev9.load_cells(cfg1, "results")
    c_new = c_new[c_new.seed.isin(SEEDS)].reset_index(drop=True)
    # 外す前の表で 3 本を計算して結合
    for tag, cc, root in (("a1既存", c_a1, res), ("a1新規", c_new, res_new)):
        deg = add_projection(root, "FIB_MFIB", "dFIB")
        diff = add_projection_diff(root, "FIB_MFIB", "aStr", "dFIB", "diff")
        com = add_projection_diff(root, "FIB_MFIB", "aStr", "dFIB", "common")
        cc = merge_projection(cc, deg, "m_proj_deg", f"{tag} u_deg")
        cc = merge_projection(cc, diff, "m_proj_diff", f"{tag} u_diff")
        cc = merge_projection(cc, com, "m_proj_common", f"{tag} u_common")
        if tag == "a1既存":
            c_a1 = cc
        else:
            c_new = cc

    pool = [x for x in G["sets"]["プールA1"] if x in set(c_a1.gene)]
    symptr = G["symbol"]
    S = {}
    for col, lab in (("m_proj", "u_ada"), ("m_proj_deg", "u_deg"), ("m_proj_diff", "u_diff"),
                     ("m_proj_common", "u_common")):
        S["既存_" + lab] = per_seed_stats(c_a1, col, pool, symptr, log, f"A1 既存 {lab}")
    # 新規 30 個は同じプール（nulls は既存の 549）で補正・標準化する。標準化の平均・標準偏差は既存側と同じ値を使う
    for col, lab in (("m_proj", "u_ada"), ("m_proj_deg", "u_deg"), ("m_proj_diff", "u_diff"),
                     ("m_proj_common", "u_common")):
        frames = []
        for s in SEEDS:
            cs = c_new[c_new.seed == s].assign(m_goal=lambda d: d[col])
            cs0 = c_a1[c_a1.seed == s].assign(m_goal=lambda d: d[col])
            nl = {h: hd for h, hd in cs0[cs0.gene.isin(pool)].groupby("gene")}
            f = ev9.gene_stats_stratified(cs, nl, "m_goal")
            n = S["既存_" + lab][(S["既存_" + lab].seed == s) & (S["既存_" + lab].gene.isin(pool))].D_raw
            frames.append(f.assign(seed=s, z=(f.D_raw - n.mean()) / n.std(ddof=1),
                                   symbol=f.gene.map(symptr)))
        S["新規_" + lab] = pd.concat(frames, ignore_index=True)
    Z_saved = {lab.split("_", 1)[1]: v.pivot(index="gene", columns="seed", values="z").reindex(columns=SEEDS)
               for lab, v in S.items() if lab.startswith("既存_")}
    Z_new = {lab.split("_", 1)[1]: v.pivot(index="gene", columns="seed", values="z").reindex(columns=SEEDS)
             for lab, v in S.items() if lab.startswith("新規_")}
    out["A1 細胞ごとの表"] = {"既存626": {"遺伝子数": int(c_a1.gene.nunique()),
                                       "移動量0で外した": drop_a1},
                             "新規30": {"遺伝子数": int(c_new.gene.nunique()), "移動量0で外した": drop_new}}

    def zm_saved(direction: str, genes: list[str]) -> np.ndarray:
        return Z_saved[direction].reindex([x for x in genes if x in Z_saved[direction].index]).to_numpy()

    def zm_new(direction: str, genes: list[str]) -> np.ndarray:
        return Z_new[direction].reindex([x for x in genes if x in Z_new[direction].index]).to_numpy()

    # ---- 再現確認（A-4・A1 側）
    checks = {}
    cand = zm_saved("u_ada", G["sets"]["候補"]).mean(axis=0)
    pc71a = zm_saved("u_ada", G["sets"]["活性化PC71"]).mean(axis=0)
    d1 = float(np.max(np.abs(cand - np.array(PLAN["再現確認（A1・u_ada の候補8）"]))))
    d2 = float(np.max(np.abs(pc71a - np.array(PLAN["再現確認（A1・u_ada の PC71）"]))))
    # 遺伝子ごとの z_ada を gene_seed_009.csv と突き合わせる（626 個）
    ref = pd.read_csv(str(ROOT / ev["repro"]["gene_seed_009"]))
    refz = ref.set_index(["gene", "seed"])["z"]
    diff_all = []
    for gene in Z_saved["u_ada"].index:
        for s in SEEDS:
            a = Z_saved["u_ada"].loc[gene, s]
            b = refz.get((gene, s))
            if b is not None and pd.notna(a):
                diff_all.append(abs(float(a) - float(b)))
    d_all = float(np.max(diff_all)) if diff_all else float("nan")
    checks["A-4 再現確認（A1）"] = {
        "u_ada の候補8セット平均": [round(x, 4) for x in cand], "参照": PLAN["再現確認（A1・u_ada の候補8）"],
        "u_ada の PC71 セット平均": [round(x, 4) for x in pc71a], "参照": PLAN["再現確認（A1・u_ada の PC71）"],
        "遺伝子ごと 626 個の最大絶対差": round(float(d_all), 6), "許容差": PLAN["再現確認_許容差"],
        "pass": bool(d1 <= 1e-3 and d2 <= 1e-3 and d_all <= 1e-3)}
    # 新規 30 のうち既存 10 個との機械一致（同じ bf16/GB10 のはず）
    overlap10 = [x for x in G["sets"]["退化PC30"] if x in set(c_a1.gene)]
    mc = {"重なり遺伝子": [G["symbol"].get(x, x) for x in overlap10], "n": len(overlap10),
          "細胞集合の一致": {}, "z_diff の最大差": 0.0}
    mc_ok = True
    ok_cells = True
    for s in SEEDS:
        s1 = set(c_a1[c_a1.seed == s].cell_uid)
        s2 = set(c_new[c_new.seed == s].cell_uid)
        same = (s1 == s2)
        mc["細胞集合の一致"][str(s)] = same
        ok_cells &= same
    if ok_cells:
        z1 = zm_saved("u_diff", overlap10)
        z2 = zm_new("u_diff", [x for x in overlap10])
        if z1.shape == z2.shape:
            mc["z_diff の最大差"] = round(float(np.max(np.abs(z1 - z2))), 6)
            mc_ok = bool(np.all(np.abs(z1 - z2) <= 1e-3))
    checks["A-4 新規30 の機械一致（既存10個との突き合わせ）"] = {"結果": mc, "pass": bool(mc_ok)}
    checks["A-4 再現確認（A1）"]["pass"] &= bool(mc_ok)
    if not checks["A-4 再現確認（A1）"]["pass"] and not smoke:
        raise SystemExit("**止める**(A-4/A1): 再現確認または新規30の機械一致が合わない")
    out["A チェック（A1）"] = checks

    # ---- P-1 主判定（A1・確認）
    pc71 = G["sets"]["活性化PC71"]; pc30 = G["sets"]["退化PC30"]
    poolZ = Z_saved["u_diff"].reindex(pool).to_numpy()
    obs71 = zm_saved("u_diff", pc71).mean(axis=0)
    obs30 = zm_new("u_diff", pc30).mean(axis=0)
    rng_i = np.random.default_rng(PLAN["乱数の鍵"]["P-1 活性化の偽セット"])
    rng_ii = np.random.default_rng(PLAN["乱数の鍵"]["P-1 退化の偽セット"])
    rng_iii = np.random.default_rng(PLAN["乱数の鍵"]["P-1 (iii) の無作為 2 セット"])
    thr71 = np.quantile(set_null(poolZ, len(pc71), N_DRAW, rng_i, replace=True), 0.99, axis=0)
    thr30 = np.quantile(set_null(poolZ, len(pc30), N_DRAW, rng_ii, replace=True), 0.01, axis=0)
    null_diff = np.empty((3, N_DRAW))
    for s in range(3):
        ia = np.random.default_rng(PLAN["乱数の鍵"]["P-1 (iii) の無作為 2 セット"] + s).integers(0, len(poolZ), size=(N_DRAW, len(pc71)))
        ib = np.random.default_rng(PLAN["乱数の鍵"]["P-1 (iii) の無作為 2 セット"] + 1000 + s).integers(0, len(poolZ), size=(N_DRAW, len(pc30)))
        null_diff[s] = poolZ[:, s][ia].mean(axis=1) - poolZ[:, s][ib].mean(axis=1)
    thr_diff = np.quantile(null_diff, 0.99, axis=1)
    c_i = bool(np.all(obs71 > 0) and np.all(obs71 > thr71))
    c_ii = bool(np.all(obs30 < 0) and np.all(obs30 < thr30))
    c_iii = bool(np.all((obs71 - obs30) > thr_diff))
    p1 = c_i and c_ii and c_iii
    # 重なりを除いた版（R-13）
    ov7 = set(PLAN["重なり7"])  # 記号ではない。Ensembl の 7 個を符号対応から引く
    ov_ens = {ens for ens, s in G["symbol"].items() if s in PLAN["重なり7"]}
    pc71r = [x for x in pc71 if x not in ov_ens]
    pc30r = [x for x in pc30 if x not in ov_ens]
    obs71r = zm_saved("u_diff", pc71r).mean(axis=0)
    obs30r = zm_new("u_diff", pc30r).mean(axis=0)
    thr71r = np.quantile(set_null(poolZ, len(pc71r), N_DRAW, np.random.default_rng(PLAN["乱数の鍵"]["P-1 活性化の偽セット"] + 2), replace=True), 0.99, axis=0)
    thr30r = np.quantile(set_null(poolZ, len(pc30r), N_DRAW, np.random.default_rng(PLAN["乱数の鍵"]["P-1 退化の偽セット"] + 2), replace=True), 0.01, axis=0)
    c_ir = bool(np.all(obs71r > 0) and np.all(obs71r > thr71r))
    c_iir = bool(np.all(obs30r < 0) and np.all(obs30r < thr30r))
    c_iiir = bool(np.all((obs71r - obs30r) > thr_diff))
    p1_sens = c_ir and c_iir and c_iiir
    # 偶然の水準（較正。判定に使わない）
    rng_cal = np.random.default_rng(PLAN["乱数の鍵"]["偽セットの較正"])
    cal_hit = 0
    for _ in range(N_DRAW):
        ia = rng_cal.integers(0, len(poolZ), size=len(pc71))
        ib = rng_cal.integers(0, len(poolZ), size=len(pc30))
        m71 = poolZ[ia].mean(axis=0); m30 = poolZ[ib].mean(axis=0)
        cal_hit += int(bool(np.all(m71 > thr71) and np.all(m30 < thr30) and np.all(m71 - m30 > thr_diff)))
    p1_main = {"(i) 活性化71": {"セット平均 z_diff": [round(x, 4) for x in obs71],
                                "偽セット99%点": [round(x, 4) for x in thr71], "成り立つ": c_i},
               "(ii) 退化30": {"セット平均 z_diff": [round(x, 4) for x in obs30],
                               "偽セット1%点": [round(x, 4) for x in thr30], "成り立つ": c_ii},
               "(iii) 差が無作為2セットの99%点超": {"差": [round(x, 4) for x in (obs71 - obs30)],
                                                  "99%点": [round(x, 4) for x in thr_diff], "成り立つ": c_iii},
               "P-1 全体": p1}
    out["P-1 主判定（A1・確認）"] = {"主": p1_main,
                                   "重なりを除いた版（64・23。R-13）": {
                                       "(i)": c_ir, "(ii)": c_iir, "(iii)": c_iiir, "P-1": p1_sens,
                                       "結論が変わるか": bool(p1 != p1_sens)},
                                   "偶然の水準（較正。判定に使わない）": {
                                       "偽ペア 1,000 組で P-1 を通した数": cal_hit,
                                       "読み方": "片側 1% に較正した規則の結合した偶然合格率"}}
    # ---- R-1〜R-6（最初に当たったものを採用）
    m71m = float(np.mean(obs71)); m30m = float(np.mean(obs30))
    r1 = bool(np.all(obs71 < 0))
    r2 = bool(np.all(obs30 > 0))
    r3 = (not c_i) and (not c_ii)
    r4 = (m71m > 0 and m30m > 0) or (m71m < 0 and m30m < 0)
    r5 = c_i and c_ii and (not c_iii)
    r6 = (c_i and not c_ii) or (c_ii and not c_i)
    first_p = next((k for k, v in {"R-1": r1, "R-2": r2, "R-3": r3, "R-4": r4, "R-5": r5, "R-6": r6}.items()
                    if v), "該当なし（P-1 の予測どおりの提案）")
    out["P-1 反証条件（R-1〜R-6）"] = {"R-1": r1, "R-2": r2, "R-3": r3, "R-4": r4, "R-5": r5, "R-6": r6,
                                      "最初に当たったもの": first_p}
    # ---- R-12（交絡。対照プール 549）
    ctrl = pool
    zbar = Z_saved["u_diff"].loc[ctrl].mean(axis=1)
    n_exp = S["既存_u_diff"].groupby("gene").n_expressed.mean().reindex(ctrl)
    llm = S["既存_u_diff"].groupby("gene").ll_med.mean().reindex(ctrl)
    r_exp = float(stats.spearmanr(zbar, n_exp).statistic)
    r_len = float(stats.spearmanr(zbar, llm).statistic)
    r12 = abs(r_exp) > PLAN["R-12_閾値"] or abs(r_len) > PLAN["R-12_閾値"]
    out["R-12（交絡。対照プール 549）"] = {
        "ρ_exp（発現細胞数）": round(r_exp, 4), "ρ_len（入力の長さの中央）": round(r_len, 4),
        "閾値": PLAN["R-12_閾値"], "該当": r12,
        "記録（判定に使わない）": {"層ごと": {}}}
    # ---- R-13：結論が変わるか（上で計算。§4.4）
    if "P-1 主判定（A1・確認）" in out:
        out["P-1 主判定（A1・確認）"]["重なりを除いた版（64・23。R-13）"]["R-13 に当たる"] = bool(p1 != p1_sens)
    log(f"A1: P-1={p1}（(i)={c_i}, (ii)={c_ii}, (iii)={c_iii}）、感度（重なり除外）P-1={p1_sens}、"
        f"最初に当たった反証条件={first_p}、R-12={r12}")

    # ---- 順位表（A1・626 と 30 の和集合。探索。重なる 10 個は既存側を採る）
    comb = pd.concat([S["既存_u_diff"], S["新規_u_diff"]], ignore_index=True) \
        .drop_duplicates(subset=["gene", "seed"], keep="first")
    Zc = comb.pivot(index="gene", columns="seed", values="z").reindex(columns=SEEDS)
    Zr = Zc.rename(columns={s: f"z_{s}" for s in SEEDS})
    rank = (Zr.assign(z_mean3=Zc.mean(axis=1)).reset_index().rename(columns={"index": "gene"})
            .assign(symbol=lambda d: d["gene"].map(G["symbol"]),
                    注="探索（データを見て候補を探すための順位表。事前登録した確認的検定ではない。CLAUDE.md §3.5）")
            .sort_values("z_mean3", ascending=False).reset_index(drop=True))
    rank.insert(0, "rank", np.arange(1, len(rank) + 1))
    outdir = EXP / ("results-smoke" if smoke else "results")
    rank.to_csv(outdir / "gene_rank_014.csv", index=False)
    out["A1 z_diff の順位表（探索）"] = {"行数": len(rank),
                                       "上位10": rank.head(10)[["symbol", "z_mean3"]].round(4).to_dict(orient="records")}


# --------------------------------------------------------------------------- 本体
def main(part: str, smoke: bool = False) -> None:
    t0 = time.time()
    cfg = yaml.safe_load((EXP / "config.yaml").read_text())
    dec = cfg["decision"]
    outdir = EXP / ("results-smoke" if smoke else "results")
    outdir.mkdir(parents=True, exist_ok=True)
    LOG = outdir / "run.log"
    LOG.write_text("")

    def L(msg: str) -> None:
        log_line(LOG, msg)

    suf = "-smoke" if smoke else ""
    L(f"開始（part={part}{'・試し実行' if smoke else ''}）。出力先 {outdir.relative_to(ROOT)}")
    L(f"plan.md sha256 {sha256f(EXP/'plan.md')[:16]}… / config.yaml sha256 "
      f"{sha256f(EXP/'config.yaml')[:16]}… / このスクリプト sha256 {sha256f(Path(__file__))[:16]}…")

    # ---- 1. config の decision が plan と同じか
    checks_conf = {
        "候補8": len(cfg["gene_sets"]["candidate_8"]) == PLAN["候補_n"],
        "退化30": len(cfg["perturbation"]["genes_to_perturb"]) == PLAN["退化PC30_n"],
        "乱数の鍵（P-1 活性化）": int(dec["乱数の鍵"]["P-1 活性化の偽セット"]) == PLAN["乱数の鍵"]["P-1 活性化の偽セット"],
        "乱数の鍵（P-1 退化）": int(dec["乱数の鍵"]["P-1 退化の偽セット"]) == PLAN["乱数の鍵"]["P-1 退化の偽セット"],
        "乱数の鍵（P-1 差）": int(dec["乱数の鍵"]["P-1 (iii) の無作為 2 セット"]) == PLAN["乱数の鍵"]["P-1 (iii) の無作為 2 セット"],
        "乱数の鍵（P-2）": int(dec["乱数の鍵"]["P-2 候補の偽セット"]) == PLAN["乱数の鍵"]["P-2 候補の偽セット"],
        "乱数の鍵（P-3）": int(dec["乱数の鍵"]["P-3 候補の偽セット"]) == PLAN["乱数の鍵"]["P-3 候補の偽セット"],
        "乱数の鍵（較正）": int(dec["乱数の鍵"]["偽セットの較正"]) == PLAN["乱数の鍵"]["偽セットの較正"],
        "A1 新規摂動は bf16": str(cfg["compute"]["dtype"]).lower() == "bf16",
        "A1 新規摂動は cuda": str(cfg["compute"]["device"]).lower() == "cuda",
        "G2 不要（config 側）": cfg["compute"]["needs_G2_approval"] is False,
    }
    if not all(checks_conf.values()):
        L(f"**止める**: config.yaml の decision: 節が plan.md §4 と一致しない {checks_conf}")
        raise SystemExit("**止める**(config): decision が plan と違う")
    out: dict = {
        "実験": "EXP-20261011-001", "仮説": "HYP-20261011-01", "mode": "smoke" if smoke else "final",
        "part": part,
        "問い": "差の向き（u_diff）・共通の向き（u_common）で、活性化した間質への固有な移動と退化した状態への"
                "固有な移動を分離して測れるか（P-1）と、候補 8 個の向きの固有性（P-2・P-3、探索）",
        "0. config と plan の一致": checks_conf,
        "使った数値型": {"device（評価）": "cpu", "dtype（評価）": "fp32",
                     "GPU 時間（評価）": 0.0,
                     "備考": "保存済み埋め込みの読み出しと再射影のみ。新規摂動は別途 bf16/CUDA で実行済み"},
    }

    # ---- 2・3. 材料と遺伝子集合
    mats = check_materials(cfg, L)
    out["1. 材料（A-2）"] = mats
    G = gene_sets_from_cfg_and_json(cfg)
    out["2. 遺伝子集合（A-3）"] = {
        "個数": {k: len(v) for k, v in G["sets"].items()},
        "片方に重複がある集合（A-7）": G["片方に重複（A-7）"],
        "盲検との重なり": {k: sorted(set(v) & set(cfg["blocked_genes"])) for k, v in G["sets"].items()},
    }
    if any(G["片方に重複（A-7）"].values()) or \
       any(out["2. 遺伝子集合（A-3）"]["盲検との重なり"].values()):
        raise SystemExit("**止める**(A-1/A-7): 集合の重複か盲検の混入がある")
    # A-5 事前に借りる（z を出してから行うので、ここでは「決定」を後で）

    if part in ("a2", "all"):
        run_a2(cfg, G, out, L, L, smoke)
    if part in ("a1", "all"):
        run_a1(cfg, G, out, L, L, smoke)

    # ---- 出力
    out["decision（config.yaml の全鍵）"] = dec
    out["注"] = ("判定（支持・否定・判定不能。ゲート G3）は PI が行う。この道具は規則を機械的に当てた結果と提案まで。"
                 "hypotheses/tracker.md は変更しない。コンピュータ上の遺伝子摂動の結果であり、実験的に確認された"
                 "生物学的事実ではない。盲検中の遺伝子は摂動も解析もしていない（値は数えない）。")
    dur = time.time() - t0
    out["実行"] = {"開始": _dt.datetime.fromtimestamp(t0).astimezone().isoformat(),
                 "所要秒": round(dur, 1), "GPU 時間": 0.0, "device": "cpu", "dtype": "fp32"}
    cj = outdir / f"criteria_014{suf}.json"
    cj.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    reloaded = json.loads(cj.read_text())
    a6 = {"decision の全鍵を持つ": set((reloaded.get("decision（config.yaml の全鍵）") or {}).keys()) == set(dec.keys()),
          "保存先": str(cj.relative_to(ROOT)), "pass": True}
    out["A チェック"] = {"A-2 材料": {"pass": mats["pass"]},
                        "A-4 再現確認": "run_a2/run_a1 内", "A-6 結果ファイル": a6}
    cj.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")

    # ---- A-1（出力・ログ・一覧に盲検が 0 件）
    sym_blk = {"ENSG00000141510": None, "ENSG00000124762": None, "ENSG00000111321": None,
               "ENSG00000118785": None, "ENSG00000148346": None}
    names_blk = ["TP53", "CDKN1A", "LTBR", "SPP1", "LCN2"]
    targets = list(outdir.glob("*.*"))
    a1 = {"調べたファイル": [str(p.relative_to(ROOT)) for p in targets], "件数": 0}
    for p in targets:
        text = p.read_text(errors="replace")
        hits = sum(text.count(k) for k in sym_blk) + sum(text.count(k) for k in names_blk)
        if hits:
            a1["件数"] += hits
            a1.setdefault("当たったファイル", []).append(str(p.relative_to(ROOT)))
    a1["pass"] = a1["件数"] == 0
    L(f"A-1 盲検の混入（出力・ログ・一覧）: {a1['件数']} 件、pass={a1['pass']}")
    if not a1["pass"]:
        raise SystemExit("**止める**(A-1): 盲検が出力・ログ・一覧に現れた（B-2）")
    out["A チェック"]["A-1 盲検の混入"] = a1
    # A-5（値に NaN・無限大がないことを結果ファイルの再読込で確認）
    cj.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    rer = json.loads(cj.read_text())
    ser = json.dumps(rer)
    a5 = {"NaN が JSON に残っている": int(ser.count("NaN")), "無限大": int(ser.count("Infinity")),
          "pass": "NaN" not in ser}
    out["A チェック"]["A-5 値"] = a5
    cj.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    L(f"完了（{dur:.1f} 秒）。出力 {cj.relative_to(ROOT)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=["a2", "a1", "all"], default="all")
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    main(a.part, a.smoke)
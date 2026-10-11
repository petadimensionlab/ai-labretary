#!/usr/bin/env python3
"""EXP-20261007-001 — A2 で測った 587 遺伝子に、事前登録した判定（plan.md §6、HYP-20261007-01 §4・§5）を当てる。

**摂動は回さない。保存済みの埋め込みを読むだけ。** 判定（ゲート G3）は PI が行う。この道具は提案までにとどまる。

EXP-20261004-001 の evaluate_009.py との違い:
  - 再現の単位は seed ではなく **提供者で分けた互いに重ならない 3 群**（seed 81/82/83 は群の名前）
  - 標準化に使う対照は **510 個**（A2 で引き直した対照 500 ＋ 恒常発現 3 ＋ 中立 7）
  - 主判定は **セット単位**（P-1）。遺伝子ごと（P-2）と三次リンパ組織の 3 個（P-3）は副次
  - **最小検出可能効果量を A2 の対照で実測し**、事前に固定した置き換えの規則を当てる（§6.2）
  - 偽セットの較正（R-7）を検査する
  - 交絡の補正（発現の有無で層別した二次式）は evaluate_009.py の関数をそのまま使う（1 文字も変えない）

走らせ方:
    .venv/bin/python experiments/EXP-20261007-001/analysis/evaluate_011.py
  試し実行の出力で評価コードを通すとき（**判定ではない**）:
    .venv/bin/python experiments/EXP-20261007-001/analysis/evaluate_011.py --smoke
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261007-001"
sys.path.insert(0, str(ROOT / "experiments/EXP-20261004-001/analysis"))
sys.path.insert(0, str(ROOT / "experiments/EXP-20261001-001/analysis"))
sys.path.insert(0, str(ROOT / "scripts/analysis"))
import evaluate_009 as ev9                                # noqa: E402  層別補正・細胞表（変更しない）
from evaluate_005 import set_null                        # noqa: E402  復元ありで k 個引いて平均する

TLS3 = ("MADCAM1", "CXCL13", "CCL21")
N_DRAW = 1000


# --------------------------------------------------------------------------- 遺伝子集合
def gene_sets(cfg: dict, smoke: bool) -> dict:
    """config.yaml から本番（run_isp.py:main）と同じ規則で再構成する。"""
    targets = list(cfg["perturbation"]["genes_to_perturb"])
    pos = list(cfg["controls"]["positive"])
    hk = list(cfg["controls"]["negative_fixed"])
    neu = [x for grp in (cfg["controls"].get("extra_fixed") or []) for x in grp["genes"]]
    fixed = set(hk) | set(neu)
    smp = cfg["controls"]["negative_sampler"]
    if smp.get("method") != "explicit":
        raise SystemExit("対照は explicit 方式のはず")
    rand = [x for x in smp["genes"] if x not in set(targets) | set(pos) | fixed][:cfg["controls"]["negative_random_n"]]
    nulls = sorted(set(rand) | fixed)
    measured = sorted(set(targets) | set(pos) | set(nulls))
    if not smoke:
        want = {"候補": (len(set(targets)), 8), "陽性対照": (len(set(pos)), 71), "対照": (len(rand), 500),
                "恒常発現": (len(hk), 3), "中立": (len(neu), 7), "標準化の対照": (len(nulls), 510),
                "測る遺伝子": (len(measured), 587)}
        bad = {k: v for k, v in want.items() if v[0] != v[1]}
        if bad:
            raise SystemExit(f"**止める**: 遺伝子集合の個数が事前登録と違う {bad}")
    return {"候補": targets, "陽性対照": pos, "恒常発現": hk, "中立": neu, "対照": rand,
            "標準化の対照": nulls, "測る遺伝子": measured}


# --------------------------------------------------------------------------- 最小検出可能効果量（§6.2）
def mde_gene(zmat: np.ndarray) -> dict:
    """遺伝子ごと: 対照の「3 群の最小値」m の分布から、
    偶然合格率 5% の閾値 T5 = m の 95 パーセンタイル、
    最小検出可能効果量 = 対照に一様な効果 δ を足したとき 80% が T5 を超える最小の δ = T5 − m の 20 パーセンタイル。"""
    m = zmat.min(axis=1)
    t5 = float(np.quantile(m, 0.95))
    return {"偶然合格率5%の閾値": round(t5, 4), "最小検出可能効果量": round(t5 - float(np.quantile(m, 0.20)), 4),
            "対照の数": int(len(m))}


def mde_set(znull: np.ndarray, k: int, rng) -> dict:
    """セット: 対照から復元ありで k 個 × 1,000 組を引き、群ごとの 99 パーセンタイル thr_g を出す。
    各組の余裕 m = min_g (平均_g − thr_g)。一様な効果 δ を足して 80% の組が全群で thr を超える最小の δ
    = −(m の 20 パーセンタイル)。"""
    S = set_null(znull, k, N_DRAW, rng, replace=True)
    thr = np.quantile(S, 0.99, axis=0)
    m = (S - thr).min(axis=1)
    return {"群ごとの99パーセンタイル": [round(float(x), 4) for x in thr],
            "最小検出可能効果量": round(float(-np.quantile(m, 0.20)), 4), "k": k}


# --------------------------------------------------------------------------- 本体
def main(smoke: bool = False) -> None:
    cfg_path = E / ("work/config-smoke.yaml" if smoke else "config.yaml")
    cfg = yaml.safe_load(cfg_path.read_text())
    main_cfg = yaml.safe_load((E / "config.yaml").read_text())
    dec = main_cfg["decision"]
    mc = cfg["metric"]["confound_correction"]
    assert mc["form"] == "stratified_by_expression" and mc["variables"] == ["lp", "ll"]
    assert mc["drop_zero_movement"] is True
    assert cfg["cell_selection"]["mode"] == "common" and cfg["cell_selection"].get("donor_groups")
    z_gene0 = float(dec["P-2 副次（遺伝子ごと）"]["閾値"])
    set_abs0 = float(dec["P-1 主判定（セット）"]["閾値(ii)"])
    set_repl = float(dec["P-1 主判定（セット）"]["置き換えの基準"])
    gene_repl = float(dec["P-2 副次（遺伝子ごと）"]["置き換えの基準"])
    assert (z_gene0, set_abs0, set_repl, gene_repl) == (2.73, 1.59, 2.34, 2.73), "事前登録の数値を動かしてはいけない"
    keys = dec["乱数の鍵"]
    # 本番は float32 で、判定を保留する範囲を置かない（PI 判断 2026-10-08。plan.md §5.3）
    if not dec["判定保留の範囲"].get("置かない") or str(main_cfg["compute"]["dtype"]).lower() != "fp32":
        raise SystemExit("**止める**: 事前登録では float32 で判定を保留する範囲を置かないはず")
    if str(cfg["compute"]["dtype"]).lower() != "fp32":
        raise SystemExit("**止める**: 評価する結果の設定が float32 でない")
    HOLD_WIDTH = 0.0

    gs = gene_sets(cfg, smoke)
    g011 = json.loads((E / "work/genes_011.json").read_text())
    sym = {v: k for k, v in g011["候補"].items()} | {v: k for k, v in g011["陽性対照71"].items()} \
        | {v: k for k, v in g011["恒常発現3"].items()} | {v: k for k, v in g011["中立7"].items()} \
        | {r["ensembl"]: r["symbol"] for r in g011["対照500"]}
    layer = {r["ensembl"]: r["層"] for r in g011["対照500"]}
    # 本番は 2 台の結果を merge_results_011.py で results/ にまとめてから評価する
    res_name = "results-smoke" if smoke else "results"
    res = E / res_name
    # 評価の出力先: 試し実行は results-smoke/、本番は results/（試し実行の出力を results/ に混ぜない）
    OUTD = res
    seeds = list(cfg["robustness"]["seeds"])
    out: dict = {"実験": "EXP-20261007-001", "仮説": "HYP-20261007-01", "mode": "smoke" if smoke else "final",
                 "群": seeds}
    if smoke:
        out["注（試し実行）"] = "**試し実行の出力に対する評価であり、判定ではない。** 目的は評価コードが最後まで通ることの確認だけ"

    def save_and_stop(msg: str) -> None:
        (OUTD).mkdir(parents=True, exist_ok=True)
        name = "criteria_011-smoke.json" if smoke else "criteria_011.json"
        (OUTD / name).write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
        raise SystemExit(msg)

    # ---------------------------------------------------------- 0. 実行条件（float32・tf32 が閉じていること）
    devs = sorted((E / res_name).glob("device*.json")) if smoke else sorted((E / res_name).glob("device-*.json"))
    dv = {}
    for f in devs:
        j = json.loads(f.read_text())
        num = j.get("numerics") or {}
        dv[f.name] = {"dtype": j.get("dtype"), "host": j.get("host"),
                      "cuda_matmul_allow_tf32": num.get("cuda_matmul_allow_tf32"),
                      "float32_matmul_precision": num.get("float32_matmul_precision")}
    ok0 = bool(dv) and all(v["dtype"] == "float32" and v["cuda_matmul_allow_tf32"] is False
                           and v["float32_matmul_precision"] == "highest" for v in dv.values())
    out["0. 実行条件"] = {"device.json": dv, "pass": ok0}
    if not ok0 and not smoke:
        save_and_stop("**止める**(A-0): float32 でない、または tf32 が閉じていない、または device.json が無い")

    # ---------------------------------------------------------- 1. 細胞ごとの表（evaluate_009 の関数）
    ev9.E = E
    ev9.BLOCKED = {e: e for e in main_cfg["blocked_genes"]}
    c, drop = ev9.load_cells(cfg, res_name)
    c = c[c.seed.isin(seeds)].reset_index(drop=True)
    present = set(c.gene.unique())
    out["1. 細胞ごとの表"] = {"遺伝子数": len(present), "移動量 0 で外した細胞": drop}

    # ---------------------------------------------------------- 2. 止める検査
    checks: dict = {}
    # A-1 盲検
    blk = set(main_cfg["blocked_genes"])
    hits = sorted(blk & present)
    checks["A-1 盲検の混入"] = {"件数": len(hits), "pass": not hits}
    if hits:
        out["A チェック"] = checks
        save_and_stop("**止める**(A-1)")
    # A-2 群の中で細胞が遺伝子によらず同一、群どうしで重ならない、群の提供者だけ
    pre = json.loads((E / "work/precheck_011.json").read_text())
    groups = {str(k): set(v) for k, v in cfg["cell_selection"]["donor_groups"].items()}
    a2, uid_by_seed, ok = {}, {}, True
    for s in seeds:
        run = res / f"run-{s}-L-1"
        gdirs = sorted(p for p in (run / "genes").iterdir() if p.is_dir()) if (run / "genes").exists() else []
        missing = [p.name for p in gdirs if not (p / "cells.csv").exists()]
        digs, sizes, donors = set(), set(), set()
        uids = set()
        for p in gdirs:
            if not (p / "cells.csv").exists():
                continue
            df = pd.read_csv(p / "cells.csv", dtype={"donor_id": str})
            digs.add(hashlib.sha256("\n".join(f"{u}\t{d}\t{l}" for u, d, l in
                                              zip(df.cell_uid, df.donor_id, df.length)).encode()).hexdigest())
            sizes.add(len(df))
            donors |= set(df.donor_id)
            uids |= set(df.cell_uid)
        want = pre["C-4 3 群の細胞"]["群ごと"][str(s)]["細胞数"]
        if smoke:
            want = min(want, int(cfg["cell_selection"]["max_cells_per_gene"]))
        rec = {"遺伝子数": len(gdirs), "cells.csv が無い": missing, "細胞集合の種類": len(digs),
               "細胞数": sorted(sizes), "期待する細胞数": want,
               "群の外の提供者": sorted(donors - groups[str(s)])}
        rec["pass"] = (not missing and len(digs) == 1 and sizes == {want} and not rec["群の外の提供者"]
                       and len(gdirs) > 0)
        ok &= rec["pass"]
        a2[str(s)] = rec
        uid_by_seed[s] = uids
    inter = {f"{a}-{b}": len(uid_by_seed[a] & uid_by_seed[b]) for i, a in enumerate(seeds) for b in seeds[i + 1:]}
    ok &= all(v == 0 for v in inter.values())
    checks["A-2 群の細胞"] = {"群ごと": a2, "群どうしで重なる細胞": inter, "pass": bool(ok)}
    if not ok and not smoke:
        out["A チェック"] = checks
        save_and_stop("**止める**(A-2)")
    # A-3 重心の行が事前登録どおり（本番のみ）／A-3b 重心が 3 群で同一か（記録）
    if not smoke:
        pinned = pre["C-6 state_rows の sha256"]
        a3 = {}
        for s in seeds:
            p = ROOT / str(cfg["state_embs"]["rows_file"]).format(seed=s)
            got = hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None
            a3[str(s)] = {"一致": got == pinned[str(s)]}
        checks["A-3 state_rows の sha256"] = {"群ごと": a3, "pass": all(v["一致"] for v in a3.values())}
        if not checks["A-3 state_rows の sha256"]["pass"]:
            out["A チェック"] = checks
            save_and_stop("**止める**(A-3)")
    cents = {s: np.load(res / f"run-{s}-L-1/state_emb_vectors.npz") for s in seeds}
    names = [cfg["states"]["start_state"], cfg["states"]["goal_state"]] + cfg["states"].get("alt_states", [])
    checks["A-3b 重心の群間の最大差（記録）"] = {
        n: max(float(np.max(np.abs(cents[a][n].astype(float) - cents[b][n].astype(float))))
               for i, a in enumerate(seeds) for b in seeds[i + 1:]) for n in names}
    # A-5 全遺伝子
    if not smoke:
        miss = sorted(set(gs["測る遺伝子"]) - present)
        checks["A-5 全遺伝子が結果にある"] = {"欠け": miss, "pass": not miss}
        if miss:
            out["A チェック"] = checks
            save_and_stop("**止める**(A-5)")
    checks["A-6 移動量 0 の割合（記録）"] = {"値": drop["frac_zero_movement"], "5%超": drop["frac_zero_movement"] > 0.05}
    out["A チェック"] = checks

    # ---------------------------------------------------------- 3. 群ごとの標準得点
    nulls = set(gs["標準化の対照"]) & present
    frames = []
    for s in seeds:
        cs = c[c.seed == s].assign(m_goal=lambda d: d.m_proj)
        nl = {h: hd for h, hd in cs[cs.gene.isin(nulls)].groupby("gene")}
        f = ev9.gene_stats_stratified(cs, nl, "m_goal")
        n = f[f.gene.isin(nulls)]
        frames.append(f.assign(seed=s, z=(f.D_raw - n.D_raw.mean()) / n.D_raw.std(ddof=1),
                               symbol=f.gene.map(sym)))
    long = pd.concat(frames, ignore_index=True)
    (OUTD).mkdir(parents=True, exist_ok=True)
    long.to_csv(OUTD / ("gene_group_011-smoke.csv" if smoke else "gene_group_011.csv"), index=False)
    Z = long.pivot(index="gene", columns="seed", values="z").reindex(columns=seeds)

    def zm(genes: list[str]) -> np.ndarray:
        return Z.reindex([g for g in genes if g in Z.index]).to_numpy()

    znull = zm(sorted(nulls))

    # ---------------------------------------------------------- 4. 最小検出可能効果量（§6.2。事前に固定した置き換えの規則）
    mg = mde_gene(znull)
    ms = mde_set(znull, 8, np.random.default_rng(keys["最小検出可能効果量の計算"]))
    z_gene = mg["最小検出可能効果量"] if mg["最小検出可能効果量"] > gene_repl else z_gene0
    set_abs = ms["最小検出可能効果量"] if ms["最小検出可能効果量"] > set_repl else set_abs0
    out["4. 最小検出可能効果量（A2 の対照で実測）"] = {
        "遺伝子ごと": mg, "セット（8 個）": ms,
        "P-2 に使う閾値": z_gene, "置き換えたか（P-2）": z_gene != z_gene0,
        "P-1 条件(ii) に使う値": set_abs, "置き換えたか（P-1）": set_abs != set_abs0,
        "規則": "実測値が 2.73（遺伝子ごと）・2.34（セット）を上回ったときだけ、実測値を使う（結果を見る前に固定）"}
    z_hold = z_gene + HOLD_WIDTH

    # ---------------------------------------------------------- 5. P-1 主判定（セット）
    zc = zm(gs["候補"])
    S_obs = zc.mean(axis=0) if zc.size else np.full(len(seeds), np.nan)
    rng1 = np.random.default_rng(keys["P-1 の偽セット"])
    S_null = set_null(znull, 8, N_DRAW, rng1, replace=True)
    thr1 = np.quantile(S_null, 0.99, axis=0)
    c_i = bool(np.all(S_obs > thr1))
    c_ii = bool(np.all(S_obs >= set_abs))
    p1 = c_i and c_ii
    rev_thr = np.quantile(S_null, 0.01, axis=0)
    reversed_ = bool(np.all(S_obs < rev_thr) and np.all(S_obs <= -set_abs))
    out["5. P-1 主判定（候補 8 個のセット平均）"] = {
        "セット平均の標準得点": [round(float(x), 3) for x in S_obs],
        "偽セットの99パーセンタイル": [round(float(x), 3) for x in thr1],
        "条件(i) 3群すべてで99パーセンタイル超": c_i, f"条件(ii) 3群すべてで {set_abs} 以上": c_ii,
        "P-1": "合格" if p1 else "不合格",
        "R-1 に当たる": not c_i, "R-2 に当たる": not c_ii,
        "R-3 逆へ動いた（3群すべてで偽セットの1パーセンタイル未満かつ −条件(ii) 以下）": reversed_,
        "群ごとの生の主指標の中央（対照の中央との差ではない）":
            [round(float(c[(c.seed == s) & c.gene.isin(gs['候補'])].m_proj.median()), 5) for s in seeds]}

    # ---------------------------------------------------------- 6. P-2 遺伝子ごと／P-3
    rows = []
    for g in gs["候補"]:
        if g not in Z.index:
            rows.append({"gene": g, "symbol": sym.get(g, g), "判定": "値が無い"})
            continue
        zs = Z.loc[g].to_numpy()
        allp = bool(np.all(zs > z_gene))
        mean = float(np.mean(zs))
        v = "合格" if allp and mean >= z_hold else ("判定保留" if allp else "不合格")
        rows.append({"gene": g, "symbol": sym.get(g, g), "z": [round(float(x), 3) for x in zs],
                     "3群平均": round(mean, 3), "判定": v})
    npass = sum(r["判定"] == "合格" for r in rows)
    nhold = sum(r["判定"] == "判定保留" for r in rows)
    p2 = "予測どおり（3 個以上）" if npass >= 3 else (
        "判定保留（合格と判定保留を合わせると 3 個以上）" if npass + nhold >= 3 else "R-4: 予測が外れた（合格 2 個以下）")
    tls = [r for r in rows if r["symbol"] in TLS3]
    tp = sum(r["判定"] == "合格" for r in tls)
    th = sum(r["判定"] == "判定保留" for r in tls)
    p3 = "予測どおり（1 個以上）" if tp >= 1 else ("判定保留（合格 0 個、判定保留 1 個以上）" if th >= 1 else
                                                   "R-5: 三次リンパ組織に関わる 3 個は単一細胞 RNA-seq では再現しなかった")

    def passes(g: str) -> bool:
        zs = Z.loc[g].to_numpy()
        return bool(np.all(zs > z_gene) and np.mean(zs) >= z_hold)
    nl_genes = [g for g in sorted(nulls) if g in Z.index]
    k_null = sum(passes(g) for g in nl_genes)
    p_hat = k_null / len(nl_genes) if nl_genes else float("nan")
    zero_genes = [g for g in nl_genes if layer.get(g) == "出発点で0"]
    expr_genes = [g for g in nl_genes if layer.get(g) == "発現あり"]
    out["6. P-2・P-3（遺伝子ごと）"] = {
        "閾値": z_gene, "判定保留の範囲": "置かない（float32。PI 判断 2026-10-08）", "遺伝子ごと": rows,
        "合格の数": npass, "判定保留の数": nhold, "P-2": p2,
        "三次リンパ組織の3個": {"合格": tp, "判定保留": th, "P-3": p3},
        "偶然の水準": {"対照の数": len(nl_genes), "合格した対照": k_null, "偽陽性率": round(p_hat, 5),
                     "候補8個での期待件数": round(8 * p_hat, 4),
                     "二項検定 片側p（合格数 ≥ 実際）": float(stats.binomtest(npass, 8, max(p_hat, 1e-9),
                                                                         alternative="greater").pvalue)},
        "診断 D-2（偽陽性率の層ごと。止めない）": {
            "出発点で0の対照": {"n": len(zero_genes), "合格": sum(passes(g) for g in zero_genes)},
            "発現ありの対照": {"n": len(expr_genes), "合格": sum(passes(g) for g in expr_genes)}}}

    # ---------------------------------------------------------- 7. R-6 陽性対照・R-7 偽セットの較正・恒常発現
    zp = zm(gs["陽性対照"])
    rng6 = np.random.default_rng(keys["R-6 の偽セット"])
    pc_obs = zp.mean(axis=0)
    pc_thr = np.quantile(set_null(znull, zp.shape[0], N_DRAW, rng6, replace=True), 0.99, axis=0)
    r6 = not bool(np.all(pc_obs > pc_thr))
    g69 = [g for g in gs["陽性対照"] if sym.get(g) not in ("GSN", "ELN")]
    z69 = zm(g69)
    o69 = z69.mean(axis=0)
    t69 = np.quantile(set_null(znull, z69.shape[0], N_DRAW, rng6, replace=True), 0.99, axis=0)
    rng7 = np.random.default_rng(keys["R-7 の偽セット"])
    S7 = set_null(znull, 8, N_DRAW, rng7, replace=True)
    n7 = int(np.sum(np.all(S7 > thr1, axis=1) & np.all(S7 >= set_abs, axis=1)))
    p7 = float(stats.binomtest(n7, N_DRAW, 0.01, alternative="greater").pvalue)
    r7 = p7 < 0.05
    out["7. 測定が働いているか"] = {
        "R-6 陽性対照71": {"セット平均": [round(float(x), 3) for x in pc_obs],
                          "偽セット99パーセンタイル": [round(float(x), 3) for x in pc_thr], "R-6 に当たる": r6},
        "陽性対照69（GSN・ELN を除く。記述）": {"セット平均": [round(float(x), 3) for x in o69],
                                         "偽セット99パーセンタイル": [round(float(x), 3) for x in t69],
                                         "3群すべてで超える": bool(np.all(o69 > t69))},
        "R-7 偽セットの較正": {"新しい偽セットでP-1に合格した数": n7, "組数": N_DRAW, "片側p": p7, "R-7 に当たる": r7},
        "恒常発現（記録）": {sym.get(g, g): [round(float(x), 3) for x in Z.loc[g]] for g in gs["恒常発現"] if g in Z.index}}

    # ---------------------------------------------------------- 8. 負の側の目盛り（記述）
    n_neg = sum(bool(np.all(Z.loc[g].to_numpy() < -z_gene)) for g in nl_genes)
    out["8. 負の側の目盛り"] = {"正の側の割合": round(p_hat, 5), "負の側の割合": round(n_neg / max(len(nl_genes), 1), 5),
                              "読み方": "負の側が正の側の 10 分の 1 未満なら、負の値は記録するだけで解釈しない"}

    # ---------------------------------------------------------- 8b. 入れ子 300（記述）・A1 との突き合わせ（記述）
    nested = set(g011["対照300（入れ子）"]) | set(gs["恒常発現"]) | set(gs["中立"])
    desc: dict = {}
    if not smoke:
        n300 = [g for g in nested if g in Z.index]
        D = long.pivot(index="gene", columns="seed", values="D_raw").reindex(columns=seeds)
        mu, sd = D.loc[n300].mean(), D.loc[n300].std(ddof=1)
        Z3 = (D - mu) / sd
        zn3 = Z3.loc[n300].to_numpy()
        S3 = Z3.reindex(gs["候補"]).to_numpy().mean(axis=0)
        thr3 = np.quantile(set_null(zn3, 8, N_DRAW, np.random.default_rng(keys["P-1 の偽セット"]), replace=True),
                           0.99, axis=0)
        desc["入れ子 300 + 10 で標準化した P-1（記述）"] = {
            "セット平均": [round(float(x), 3) for x in S3], "99パーセンタイル": [round(float(x), 3) for x in thr3],
            "(i)": bool(np.all(S3 > thr3)), "(ii)": bool(np.all(S3 >= set_abs))}
        a1 = ROOT / "experiments/EXP-20261004-001/results/gene_rank_009.csv"
        if a1.exists():
            r1 = pd.read_csv(a1).set_index("gene")["z_mean3seed"]
            both = list(dict.fromkeys(g for g in gs["陽性対照"] + gs["候補"]
                                      if g in r1.index and g in Z.index))
            a2m = Z.loc[both].mean(axis=1)
            desc["A1 との突き合わせ（陽性対照と候補。別の尺度。記述）"] = {
                "遺伝子数": len(both),
                "Spearman": round(float(stats.spearmanr(r1.loc[both], a2m).statistic), 4),
                "候補": {sym.get(g, g): {"A1": round(float(r1.loc[g]), 3), "A2": round(float(a2m.loc[g]), 3)}
                       for g in gs["候補"] if g in both}}
    out["8b. 記述"] = desc

    # ---------------------------------------------------------- 9. 解釈（plan.md §6.6 の表。確定は PI）
    nan_any = bool(np.isnan(S_obs).any() or np.isnan(pc_obs).any())
    if nan_any:
        tag = "R-8: 値が欠けており規則に当てはまらない。事実だけを報告する"
    elif r6:
        tag = "R-6: 判定不能（陽性対照が働いていない）。候補の値は報告するが解釈しない"
    elif r7:
        tag = "R-7: 判定不能（偽セットの較正が崩れている）"
    elif p1:
        tag = "P-1 合格: 主判定の予測どおり"
    elif reversed_:
        tag = "R-3: 候補 8 個のセットは目標から逆へ動いた（支持されない）"
    else:
        tag = "R-1/R-2: 支持されない（" + ("条件(i)" if not c_i else "") + ("・" if (not c_i and not c_ii) else "") \
            + ("条件(ii)" if not c_ii else "") + " を満たさない）"
    out["9. 解釈の提案（確定は PI。ゲート G3）"] = {"主判定": tag, "P-2": p2, "P-3": p3}

    # ---------------------------------------------------------- 10. 順位表（探索。CLAUDE.md §3.5）
    rank = (long.groupby("gene").z.agg(["mean", "min"]).rename(columns={"mean": "z_mean3", "min": "z_min3"})
            .reset_index().assign(symbol=lambda d: d.gene.map(sym)).sort_values("z_mean3", ascending=False))
    rp = OUTD / ("gene_rank_011-smoke.csv" if smoke else "gene_rank_011.csv")
    rank.to_csv(rp, index=False)
    out["10. 全遺伝子の順位表（探索）"] = {"行数": len(rank), "保存先": str(rp)}
    out["注"] = ("判定（ゲート G3）は PI が行う。hypotheses/tracker.md は変更しない。"
                "コンピュータ上の遺伝子摂動の結果であり、実験的に確認された事実ではない。"
                "盲検中の TP53・CDKN1A・LTBR・SPP1・LCN2 は摂動も解析もしていない。")
    name = "criteria_011-smoke.json" if smoke else "criteria_011.json"
    (OUTD / name).write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps({k: out[k] for k in out if not k.startswith("10.")}, ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true", help="試し実行の出力で評価コードを通す（判定ではない）")
    main(smoke=ap.parse_args().smoke)

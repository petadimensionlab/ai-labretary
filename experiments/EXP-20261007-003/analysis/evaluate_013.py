#!/usr/bin/env python3
"""EXP-20261007-003 — HYP-20261007-02 段階 2。保存済み埋め込みを、2 本の方向へ射影し直す（評価のみ）。

**摂動は回さない。Geneformer も使わない。GPU 0 時間。** EXP-20261007-001 が保存した
587 遺伝子の摂動前後埋め込み・`cells.csv`・状態の重心（FIB・aStr・dFIB）を読み、
**射影する方向だけを 2 本**（活性化した間質＝aStr / 退化した状態＝dFIB）にして計算する。

- 射影の式・交絡の補正（発現の有無で層別した二次式 `lp`・`ll`）・標準化（対照 510 個）・
  移動量 0 の除外・偽セット `set_null`（復元あり）は、**EXP-20261007-001 の評価コードと同一**
  （`evaluate_009.load_cells`・`evaluate_009.gene_stats_stratified`・`evaluate_005.set_null` を
  そのまま呼ぶ。**関数の中身は 1 文字も変えていない**）。変えるのは goal に当たる重心だけ。
- 本ファイルは plan.md §5 の判定規則をそのまま実装し、`config.yaml` の `decision:` 節が
  plan.md と同じ値かを検査する（違えば止める）。
- **判定（支持／否定／判定不能）は確定しない。** ゲート G3 は PI が行う。ここが出すのは
  「規則を機械的に当てた結果」まで。
- 盲検中の 5 遺伝子（TP53・CDKN1A・LTBR・SPP1・LCN2）は、摂動の対象・対照・指標・出力の
  どこにも入れない（plan.md §6.1。§7 A-1 で 0 件を確認する）。

走らせ方:
    .venv/bin/python experiments/EXP-20261007-003/analysis/evaluate_013.py            # 本実行
    .venv/bin/python experiments/EXP-20261007-003/analysis/evaluate_013.py --smoke    # 試し実行（判定ではない）

試し実行は、保存済み埋め込みの一部（候補 4・退化の陽性対照 4・対照 30・恒常発現 3・中立 7）だけを
symlink で見せる作業ディレクトリを作り、評価コードを最後まで通す（guardrail 2026-10-01）。
**対照を 30 個入れるのは、層別補正（自分自身を抜く回帰）が成り立つため**（guardrail 2026-10-01）。
値は判定ではない。出力は `results-smoke/` に分ける。

主な入力: plan.md、config.yaml、EXP-20261007-001 の `results/`（読み取りのみ）と `work/`
          （`genes_011.json`・`precheck_011.json`・`state_rows-seed*.json`・`tokenized.dataset`）。
主な出力: `results/criteria_013.json`・`results/gene_rank_013.csv`・`results/gene_group_013.csv`・
          `results/README.md`・`results/run.log`（試し実行は `results-smoke/`）。
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
EXP = ROOT / "experiments/EXP-20261007-003"
SRC = ROOT / "experiments/EXP-20261007-001"
sys.path.insert(0, str(ROOT / "experiments/EXP-20261004-001/analysis"))
sys.path.insert(0, str(ROOT / "experiments/EXP-20261001-001/analysis"))
sys.path.insert(0, str(ROOT / "scripts/analysis"))
import evaluate_009 as ev9                                # noqa: E402  層別補正・細胞表（1 文字も変えない）
from evaluate_005 import set_null                        # noqa: E402  復元ありで k 個引いて平均する

N_DRAW = 1000
SEEDS = [81, 82, 83]

# --------------------------------------------------------------------------- plan.md §5 の判定値
# **ここを変えてはいけない。** config.yaml の decision: 節と照合し、違えば止める。
PLAN = {
    "閾値_候補": 1.59,                       # 候補 8 個（K=8、4.5/√8）
    "閾値_遺伝子ごと": 2.73,                  # 記録のみ（段階 2 は遺伝子ごとの合否を判定に使わない）
    "置き換え": "なし",                        # 段階 2 は保存済み埋め込みの再射影（plan.md §5.1）
    "候補_n": 8,
    "陽性対照30_n": 30,
    "対照500_n": 500,
    "プール_n": 510,                          # 対照 500 + 恒常発現 3 + 中立 7
    "測る遺伝子_n": 587,
    "再現確認_参照": [2.664, 2.776, 2.302],     # criteria_011.json の候補 8 セット平均
    "再現確認_許容差": 1e-3,
    "G-c_閾値": 0.20,
    "乱数の鍵": {"候補の偽セット": 20261014, "陽性対照の偽セット": 20261015, "偽セットの較正": 20261016},
    "群": SEEDS,
    "state_rows_sha256": {str(s): "3432812c65583a9902204915af32bc474c75b699e49251a62e3e3c5a79ff4de8"
                          for s in SEEDS},
    # 再現確認に使う参照値（criteria_011.json「5. P-1 主判定」の 99 パーセンタイル。記録のみ）
    "再現確認_参照99": [0.855, 0.884, 0.894],
}


# --------------------------------------------------------------------------- 小道具
def sha256f(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def log_line(path: Path, msg: str) -> None:
    line = f"[{_dt.datetime.now().strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(path, "a") as fh:
        fh.write(line + "\n")


# --------------------------------------------------------------------------- 材料（plan.md §7 手順 2・3、A-2）
def clone_tokenized(log) -> dict:
    """トークン化済み入力を作業用に複製する（guardrail 2026-09-26）。

    HF datasets は `select`/`map` のキャッシュをデータセットのディレクトリ自体に書くため、
    **EXP-001 の `work/tokenized.dataset` を直接指さない。** データ部分だけを APFS クローンで
    `EXP-003/work/tokenized.dataset` に複製し、sha256 の一致を記録する。
    """
    src = SRC / "work/tokenized.dataset"
    dst = EXP / "work/tokenized.dataset"
    dst.mkdir(parents=True, exist_ok=True)
    names = sorted([p.name for p in src.glob("data-*.arrow")]) + ["dataset_info.json", "state.json"]
    rec: dict = {"複製元": str(src.relative_to(ROOT)), "複製先": str(dst.relative_to(ROOT)), "ファイル": {}}
    ok = True
    for n in names:
        s, d = src / n, dst / n
        if not d.exists() or sha256f(d) != sha256f(s):
            r = subprocess.run(["cp", "-c", str(s), str(d)], capture_output=True, text=True)
            if r.returncode != 0:                     # APFS クローンが不可なら通常コピー
                subprocess.run(["cp", str(s), str(d)], check=True)
        a, b = sha256f(s), sha256f(d)
        rec["ファイル"][n] = {"sha256_元": a, "sha256_複製": b, "一致": a == b}
        ok &= (a == b)
    rec["pass"] = bool(ok)
    log(f"トークン化済み入力を複製し sha256 を照合した（{len(names)} 件、一致={ok}）")
    return rec


def check_materials(cfg, pre, log) -> dict:
    """A-2: 材料の sha256 が固定値と一致し、埋め込み 587 × 3 が実在すること。"""
    se = cfg["source_experiment"]
    rec: dict = {"sha256": {}, "埋め込み": {}}
    ok = True
    for name, path, pin in (
        ("genes_011.json", SRC / "work/genes_011.json", se["genes_json_sha256"]),
        ("precheck_011.json", SRC / "work/precheck_011.json", se["precheck_json_sha256"]),
        ("EXP-001/config.yaml", SRC / "config.yaml", se["config_sha256"]),
    ):
        got = sha256f(path)
        m = bool(got) and got.startswith(str(pin))
        rec["sha256"][name] = {"sha256": got, "固定値（先頭16桁）": pin, "一致": m}
        ok &= m
    rows = {}
    for s in SEEDS:
        p = SRC / f"work/state_rows-seed{s}.json"
        got = sha256f(p)
        pin = pre["C-6 state_rows の sha256"][str(s)]
        m = got == pin
        rows[str(s)] = {"sha256": got, "一致": m}
        ok &= m
    rec["state_rows"] = rows
    # 記録のみ（判定には使わない）: 使ったコードの版
    rec["使ったコード（記録）"] = {str(p.relative_to(ROOT)): sha256f(p) for p in (
        ROOT / "experiments/EXP-20261004-001/analysis/evaluate_009.py",
        ROOT / "experiments/EXP-20261001-001/analysis/evaluate_005.py",
        SRC / "analysis/evaluate_011.py",
        ROOT / "scripts/analysis/movement_correction.py",
    )}
    # 埋め込みの実在（587 遺伝子 × 3 群 × orig/pert）
    want = set(json.loads((SRC / "work/genes_011.json").read_text())["測る遺伝子"])
    n_files = 0
    for s in SEEDS:
        gdir = SRC / f"results/run-{s}-L-1/genes"
        genes, missing = [], []
        for p in sorted(gdir.iterdir()):
            if not p.is_dir():
                continue
            genes.append(p.name)
            for kind in ("orig", "pert"):
                f = p / "emb" / f"{kind}.npy"
                if f.exists():
                    n_files += 1
                else:
                    missing.append(str(f.relative_to(ROOT)))
        m = (len(genes) == PLAN["測る遺伝子_n"] and set(genes) == want and not missing)
        rec["埋め込み"][str(s)] = {"遺伝子ディレクトリ": len(genes), "orig/pert が無い": missing,
                                   "遺伝子集合が一致": set(genes) == want, "pass": m}
        ok &= m
    rec["埋め込みファイル数"] = {"実在": n_files, "期待": PLAN["測る遺伝子_n"] * 3 * 2}
    ok &= (n_files == PLAN["測る遺伝子_n"] * 3 * 2)
    rec["pass"] = bool(ok)
    log(f"A-2 材料の sha256 と埋め込み（{n_files} 個）の実在: pass={rec['pass']}")
    return rec


# --------------------------------------------------------------------------- 遺伝子集合（A-3）
def load_gene_sets(cfg, log, smoke: bool) -> dict:
    g = json.loads((SRC / "work/genes_011.json").read_text())
    cand = list(g["候補"].values())                                     # 8
    pos_ada = list(g["陽性対照71"].values())                            # 71（記述にだけ使う）
    pos_deg = list(cfg["gene_sets"]["positive_control_degenerate_30"]["genes"])   # 30
    hk = list(g["恒常発現3"].values())                                  # 3
    neu = list(g["中立7"].values())                                    # 7
    ctrl500 = [r["ensembl"] for r in g["対照500"]]                      # 500
    layer = {r["ensembl"]: r["層"] for r in g["対照500"]}
    sym = {v: k for k, v in g["候補"].items()} | {v: k for k, v in g["陽性対照71"].items()} \
        | {v: k for k, v in g["恒常発現3"].items()} | {v: k for k, v in g["中立7"].items()} \
        | {r["ensembl"]: r["symbol"] for r in g["対照500"]}
    pool = sorted(set(ctrl500) | set(hk) | set(neu))
    out = {"候補": cand, "陽性対照_活性化71": pos_ada, "陽性対照_退化30": pos_deg,
           "恒常発現": hk, "中立": neu, "対照500": ctrl500, "標準化の対照": pool, "層": layer, "symbol": sym}
    want = {"候補": (len(cand), PLAN["候補_n"]), "陽性対照_退化30": (len(pos_deg), PLAN["陽性対照30_n"]),
            "対照500": (len(ctrl500), PLAN["対照500_n"]), "標準化の対照": (len(pool), PLAN["プール_n"])}
    bad = {k: v for k, v in want.items() if v[0] != v[1]}
    # config の候補 8 と genes_011.json の候補 8 が一致すること
    if sorted(cfg["gene_sets"]["candidate_8"]) != sorted(cand):
        bad["config 候補8 != genes_011 候補"] = (sorted(cfg["gene_sets"]["candidate_8"]), cand)
    # 重なりの総当たり（plan.md §6.2 の固定値の確認。判定に使わない）
    blk = set(cfg["blocked_genes"])
    ov = {
        "陽性対照30 ∩ 対照500": sorted(set(pos_deg) & set(ctrl500)),
        "陽性対照30 ∩ 陽性対照71": sorted(set(pos_deg) & set(pos_ada)),
        "陽性対照30 ∩ 候補8": sorted(set(pos_deg) & set(cand)),
        "陽性対照30 ∩ 恒常発現3": sorted(set(pos_deg) & set(hk)),
        "陽性対照30 ∩ 中立7": sorted(set(pos_deg) & set(neu)),
        "陽性対照30 ∩ 盲検5": sorted(set(pos_deg) & blk),
        "候補8 ∩ 盲検5": sorted(set(cand) & blk),
        "対照500 ∩ 盲検5": sorted(set(ctrl500) & blk),
        "陽性対照71 ∩ 盲検5": sorted(set(pos_ada) & blk),
    }
    out["重なり"] = {k: {"件数": len(v), "遺伝子": [sym.get(x, x) for x in v]} for k, v in ov.items()}
    out["重なりの固定値と一致"] = {
        "陽性対照30 ∩ 対照500 = 23": len(ov["陽性対照30 ∩ 対照500"]) == 23,
        "陽性対照30 ∩ 陽性対照71 = 7": len(ov["陽性対照30 ∩ 陽性対照71"]) == 7,
        "陽性対照30 ∩ 候補8 = 0": len(ov["陽性対照30 ∩ 候補8"]) == 0,
        "陽性対照30 ∩ 盲検5 = 0": len(ov["陽性対照30 ∩ 盲検5"]) == 0,
        "候補8 ∩ 盲検5 = 0": len(ov["候補8 ∩ 盲検5"]) == 0,
        "対照500 ∩ 盲検5 = 0": len(ov["対照500 ∩ 盲検5"]) == 0,
        "陽性対照71 ∩ 盲検5 = 0": len(ov["陽性対照71 ∩ 盲検5"]) == 0,
    }
    out["pass"] = (not bad) and all(out["重なりの固定値と一致"].values())
    if bad and not smoke:
        raise SystemExit(f"**止める**(A-3): 遺伝子集合の個数が事前登録と違う {bad}")
    log(f"A-3 遺伝子集合: 候補 {len(cand)}・退化の陽性対照 {len(pos_deg)}・対照 {len(ctrl500)}・"
        f"標準化プール {len(pool)}。重なりの固定値と一致={out['pass']}")
    return out


# --------------------------------------------------------------------------- 射影（plan.md §4.2）
def add_projection(res: Path, start_state: str, goal_state: str) -> pd.Series:
    """`evaluate_009.load_cells` の射影ブロック（`u = unit(goal) − unit(start)` に正規化し、
    `(d·u)/‖d‖`）と**同じ式・同じ手順**で、goal に当たる重心だけを差し替えて細胞ごとの主指標を計算する。

    `load_cells` と同じく **移動量 0 の細胞を外す前の表**（`results/cells.csv`）に対して、
    同じ `(run, gene)` のまとまりで計算する（埋め込みの行はそのまとまりの行と 1 対 1）。
    返すのは `(run, gene, cell_uid)` を鍵にした Series で、`load_cells` が返した表（外した後）へ
    結合すれば、外す細胞は 1 本目（活性化の方向）と同一になる。
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
        vg, vs = sz[goal_state].astype(np.float64), sz[start_state].astype(np.float64)
        u = vg / np.linalg.norm(vg) - vs / np.linalg.norm(vs)
        u = u / np.linalg.norm(u)
        with np.errstate(divide="ignore", invalid="ignore"):
            proj[part.index.to_numpy()] = (d @ u) / dn
    idx = pd.MultiIndex.from_arrays(
        [cells.run.to_numpy(), cells.gene.to_numpy(), cells.cell_uid.to_numpy()],
        names=["run", "gene", "cell_uid"])
    return pd.Series(proj, index=idx, name="m_proj_deg")


def per_seed_stats(c: pd.DataFrame, col: str, pool: list[str], sym: dict,
                   log, label: str) -> pd.DataFrame:
    """群ごとに層別補正（evaluate_009.gene_stats_stratified をそのまま呼ぶ）→ 標準化（対照 510 個）。"""
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
    """Spearman を記録用に出す。値が欠けた遺伝子は外し、外した数を記録する。

    発現細胞が 0 の遺伝子は `lp_med` が定義できず（`lp` は載っている細胞でしか定義できない）、
    層「出発点で0」の 100 個は発現細胞数がすべて 0 で一定になる。**計算できないものは
    NaN を書かずに null と理由を書く**（JSON を素直に読めるようにするため）。
    """
    m = x.notna() & y.notna()
    xs, ys = x[m].to_numpy(dtype=float), y[m].to_numpy(dtype=float)
    if len(xs) < 3 or float(np.ptp(xs)) == 0.0 or float(np.ptp(ys)) == 0.0:
        return {"Spearman": None, "n": int(m.sum()),
                "理由": "計算できない（値を 2 つとも持つ遺伝子が 3 個未満、または一方が一定）"}
    return {"Spearman": round(float(stats.spearmanr(xs, ys).statistic), 4), "n": int(m.sum()),
            "外した数（値が欠けた遺伝子）": int((~m).sum())}


def expressed_length_median(c: pd.DataFrame, genes: list[str], log) -> pd.DataFrame:
    """G-c の ρ_len に使う「その遺伝子が発現していた細胞の入力の長さの中央値」（群ごと）。

    発現 0 の遺伝子は、その群の全細胞の入力の長さの中央値で補完し、補完数を記録する
    （plan.md §5.4）。common モードなので群の全細胞は遺伝子によらず同一。
    """
    sub = c[c.gene.isin(set(genes))]
    rows = []
    for (s, gene), part in sub.groupby(["seed", "gene"], sort=True):
        e = part[part.expressed]
        rows.append({"seed": int(s), "gene": gene, "n_expressed": int(len(e)),
                     "ll_exp_med": float(e.ll.median()) if len(e) else np.nan})
    df = pd.DataFrame(rows)
    fallback = {s: float(c[c.seed == s].drop_duplicates("cell_uid").ll.median()) for s in SEEDS}
    df["補完したか"] = df.ll_exp_med.isna()
    for s in SEEDS:
        m = df.seed == s
        df.loc[m, "ll_exp_med"] = df.loc[m, "ll_exp_med"].fillna(fallback[s])
    log(f"G-c の ρ_len: 発現していた細胞の入力の長さの中央値を計算（補完 {int(df['補完したか'].sum())} 件 / "
        f"{len(df)} 遺伝子×群）")
    return df


# --------------------------------------------------------------------------- 本体
def main(smoke: bool = False) -> None:
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
    L(f"開始（{'試し実行' if smoke else '本実行'}）。出力先 {outdir.relative_to(ROOT)}")
    L(f"plan.md sha256 {sha256f(EXP/'plan.md')[:16]}… / config.yaml sha256 "
      f"{sha256f(EXP/'config.yaml')[:16]}… / このスクリプト sha256 {sha256f(Path(__file__))[:16]}…")

    # ---------------------------------------------------------------- 1. config の decision: 節が plan と同じか
    gs_cfg = cfg["gene_sets"]
    checks_conf = {
        "config の候補8": len(gs_cfg["candidate_8"]) == PLAN["候補_n"],
        "config の陽性対照30": len(gs_cfg["positive_control_degenerate_30"]["genes"]) == PLAN["陽性対照30_n"],
        "標準化プール = 510": gs_cfg["standardize_pool_n"] == PLAN["プール_n"],
        "U-1 の閾値 = 1.59": float(dec["U-1 主判定"]["閾値"]) == PLAN["閾値_候補"],
        "U-1 の合格 = 99 パーセンタイル超 かつ 絶対値 1.59 以上":
            ("99 パーセンタイル超" in dec["U-1 主判定"]["合格"]) and ("1.59 以上" in dec["U-1 主判定"]["合格"]),
        "U-1 の置き換え = なし": str(dec["U-1 主判定"]["置き換え"]).startswith(PLAN["置き換え"]),
        "退化の方向が主対象": "合格しないこと（主対象）" in str(dec["U-1 主判定"]["方向 2 (退化した状態)"]),
        "G-b の合格 = 99 パーセンタイル超のみ":
            ("99 パーセンタイル" in dec["G-b 陽性対照（退化の方向・30 個）"]["合格"])
            and ("絶対値の条件は課さない" in dec["G-b 陽性対照（退化の方向・30 個）"]["合格"]),
        "G-c の閾値 = 0.20":
            ("0.20" in str(dec["G-c 交絡"]["相関 1"])) and ("0.20" in str(dec["G-c 交絡"]["相関 2"])),
        "乱数の鍵":
            int(dec["乱数の鍵"]["U-1 候補の偽セット"]) == PLAN["乱数の鍵"]["候補の偽セット"]
            and int(dec["乱数の鍵"]["G-b 陽性対照の偽セット"]) == PLAN["乱数の鍵"]["陽性対照の偽セット"]
            and int(dec["乱数の鍵"]["偽セットの較正"]) == PLAN["乱数の鍵"]["偽セットの較正"],
        "数値型 = fp32": str(cfg["compute"]["dtype"]).lower() == "fp32",
        "実行先 = cpu": str(cfg["compute"]["device"]).lower() == "cpu",
        "GPU 見積もり = 0": float(cfg["compute"]["gpu_hours_estimate"]) == 0.0,
        "ゲート G2 不要": cfg["compute"]["needs_G2_approval"] is False,
    }
    if not all(checks_conf.values()):
        L(f"**止める**: config.yaml の decision: 節が plan.md §5 と一致しない {checks_conf}")
        raise SystemExit("**止める**(config): decision: 節が plan.md §5 と違う")
    L(f"config.yaml の decision: 節は plan.md §5 と一致した（{len(checks_conf)} 項目）")

    out: dict = {
        "実験": "EXP-20261007-003",
        "仮説": "HYP-20261007-02（段階 2）",
        "mode": "smoke" if smoke else "final",
        "問い": "保存済み埋め込みの再射影で、候補 8 個のセット平均は活性化した間質の方向では合格し、"
                "退化した状態の方向では合格しないか。退化の方向の陽性対照 30 個は合格するか",
        "0. config と plan の一致": checks_conf,
        "使った数値型": {"device": cfg["compute"]["device"], "dtype": cfg["compute"]["dtype"],
                    "GPU 時間": 0.0, "備考": "保存済み埋め込みの読み出しと再射影のみ。Geneformer を動かさない"},
    }
    if smoke:
        out["注（試し実行）"] = ("**試し実行の出力に対する評価であり、判定ではない。** "
                              "保存済み埋め込みの一部（遺伝子 48 個）だけを見せた作業ディレクトリで、"
                              "評価コードが最後まで通ることだけを確かめる（guardrail 2026-10-01）。")

    # ---------------------------------------------------------------- 2. 材料（A-2）と複製
    pre = json.loads((SRC / "work/precheck_011.json").read_text())
    clone = clone_tokenized(L)
    mats = check_materials(cfg, pre, L)
    out["1. 材料（A-2）"] = {"複製（トークン化済み入力）": clone, "sha256 と実在": mats}
    if not mats["pass"] and not smoke:
        out["A チェック"] = {"A-2 材料": {"pass": False}}
        (outdir / f"criteria_013{suf}.json").write_text(
            json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
        raise SystemExit("**止める**(A-2): 材料の sha256 が合わない、または埋め込みが欠けている")

    # ---------------------------------------------------------------- 3. 遺伝子集合（A-3）
    G = load_gene_sets(cfg, L, smoke)
    out["2. 遺伝子集合（A-3）"] = {
        "候補8": [G["symbol"].get(x, x) for x in G["候補"]],
        "陽性対照_退化30": [G["symbol"].get(x, x) for x in G["陽性対照_退化30"]],
        "個数": {"候補": len(G["候補"]), "陽性対照_退化": len(G["陽性対照_退化30"]),
               "対照500": len(G["対照500"]), "標準化プール": len(G["標準化の対照"])},
        "対照500 の層": {k: sum(1 for x in G["対照500"] if G["層"][x] == k)
                     for k in sorted({G["層"][x] for x in G["対照500"]})},
        "重なり（総当たり・plan.md §6.2 の固定値の確認）": G["重なり"],
        "重なりの固定値と一致": G["重なりの固定値と一致"],
        "pass": G["pass"],
    }
    # 盲検（A-1）: 遺伝子集合に入っていないこと（値は数えない）
    blk = set(cfg["blocked_genes"])
    sets_hit = sorted((blk & set(G["候補"])) | (blk & set(G["陽性対照_退化30"]))
                      | (blk & set(G["対照500"])) | (blk & set(G["恒常発現"])) | (blk & set(G["中立"])))
    if sets_hit:
        L("**止める**(A-1): 盲検の遺伝子が遺伝子集合に混入している")
        raise SystemExit("**止める**(A-1): 盲検の遺伝子が遺伝子集合に混入している")
    if not G["pass"] and not smoke:
        raise SystemExit("**止める**(A-3): 遺伝子集合の個数または重なりが事前登録と違う")

    # ---------------------------------------------------------------- 4. 細胞ごとの表と 2 本の射影
    # 試し実行は、一部の遺伝子だけを symlink で見せる作業ディレクトリを作る
    if smoke:
        IN = EXP / "work/smoke-input"
        if IN.exists():
            shutil.rmtree(IN)
        smoke_genes = (G["候補"][:4] + G["陽性対照_退化30"][:4]
                       + [x for x in G["対照500"] if G["層"][x] == "発現あり"][:25]
                       + [x for x in G["対照500"] if G["層"][x] == "出発点で0"][:5]
                       + G["恒常発現"] + G["中立"])
        smoke_genes = list(dict.fromkeys(smoke_genes))
        src_cells = pd.read_csv(SRC / "results/cells.csv", dtype={"donor_id": str})
        (IN / "run-81-L-1/genes").mkdir(parents=True, exist_ok=True)
        for s in SEEDS:
            gdir = IN / f"run-{s}-L-1/genes"
            gdir.mkdir(parents=True, exist_ok=True)
            (IN / f"run-{s}-L-1/state_emb_vectors.npz").symlink_to(
                SRC / f"results/run-{s}-L-1/state_emb_vectors.npz")
            for g in smoke_genes:
                (gdir / g).symlink_to(SRC / f"results/run-{s}-L-1/genes" / g)
        src_cells[src_cells.gene.isin(smoke_genes)].to_csv(IN / "cells.csv", index=False)
        ev9.E = EXP
        res_rel = "work/smoke-input"
        res = EXP / res_rel
        L(f"試し実行の入力を作った（{len(smoke_genes)} 遺伝子 × 3 群、{IN.relative_to(ROOT)}）")
    else:
        ev9.E = SRC
        res_rel = "results"
        res = SRC / res_rel

    ev9.BLOCKED = {e: e for e in cfg["blocked_genes"]}
    cfg001 = yaml.safe_load((SRC / "config.yaml").read_text())
    cfg001["data"]["tokenized_dataset"] = str((EXP / "work/tokenized.dataset").relative_to(ROOT))
    c, drop = ev9.load_cells(cfg001, res_rel)               # ← EXP-001 の関数をそのまま使う（活性化の方向）
    c = c[c.seed.isin(SEEDS)].reset_index(drop=True)
    # 退化した状態の方向（移動量 0 を外す前の表で計算し、(run, gene, cell_uid) で結合する）
    deg = add_projection(res, "FIB", "dFIB")
    n0 = len(c)
    c = c.merge(deg, left_on=["run", "gene", "cell_uid"], right_index=True, how="left",
                validate="one_to_one")
    miss = int(c.m_proj_deg.isna().sum())
    if miss:
        raise SystemExit(f"**止める**: 退化の方向の射影が結びつかない細胞が {miss} 個（{n0} 行中）")
    present = set(c.gene.unique())
    out["3. 細胞ごとの表"] = {"遺伝子数": len(present), "移動量 0 で外した細胞": drop}
    L(f"細胞ごとの表: 遺伝子 {len(present)} 個、移動量 0 で外した細胞 {drop['n_zero_movement']} 個 "
      f"（{drop['frac_zero_movement']:.6f}）")

    # ---------------------------------------------------------------- 5. 重心の幾何（記録。plan.md §2.2 の確認）
    z = np.load(SRC / "results/run-81-L-1/state_emb_vectors.npz")
    u_ada = z["aStr"].astype(float) / np.linalg.norm(z["aStr"].astype(float)) \
        - z["FIB"].astype(float) / np.linalg.norm(z["FIB"].astype(float))
    u_deg = z["dFIB"].astype(float) / np.linalg.norm(z["dFIB"].astype(float)) \
        - z["FIB"].astype(float) / np.linalg.norm(z["FIB"].astype(float))
    geo = {
        "（記録。判定に使わない。plan.md §2.2 の値の確認）": True,
        "コサイン（FIB と aStr）": float(np.dot(z["FIB"], z["aStr"]) / (np.linalg.norm(z["FIB"]) * np.linalg.norm(z["aStr"]))),
        "コサイン（FIB と dFIB）": float(np.dot(z["FIB"], z["dFIB"]) / (np.linalg.norm(z["FIB"]) * np.linalg.norm(z["dFIB"]))),
        "コサイン（活性化の方向ベクトル と 退化の方向ベクトル）": float(
            np.dot(u_ada, u_deg) / (np.linalg.norm(u_ada) * np.linalg.norm(u_deg))),
        "単位重心間の距離 FIB→aStr": float(np.linalg.norm(
            z["aStr"].astype(float) / np.linalg.norm(z["aStr"]) - z["FIB"].astype(float) / np.linalg.norm(z["FIB"]))),
        "単位重心間の距離 FIB→dFIB": float(np.linalg.norm(
            z["dFIB"].astype(float) / np.linalg.norm(z["dFIB"]) - z["FIB"].astype(float) / np.linalg.norm(z["FIB"]))),
    }
    out["3b. 重心の幾何"] = geo
    L(f"重心の幾何（記録）: cos(FIB,aStr)={geo['コサイン（FIB と aStr）']:.6f}、"
      f"cos(FIB,dFIB)={geo['コサイン（FIB と dFIB）']:.6f}、2 方向のコサイン="
      f"{geo['コサイン（活性化の方向ベクトル と 退化の方向ベクトル）']:.6f}")

    # ---------------------------------------------------------------- 6. 層別補正 → 標準化 → z（2 方向）
    pool = [g for g in G["標準化の対照"] if g in present]
    long_ada = per_seed_stats(c, "m_proj", pool, G["symbol"], L, "活性化した間質（aStr）")
    long_deg = per_seed_stats(c, "m_proj_deg", pool, G["symbol"], L, "退化した状態（dFIB）")
    long_ada["方向"] = "adaptive"
    long_deg["方向"] = "degenerate"
    long = pd.concat([long_ada, long_deg], ignore_index=True)
    long.to_csv(outdir / f"gene_group_013{suf}.csv", index=False)
    longs = {"adaptive": long_ada, "degenerate": long_deg}
    Z = {d: longs[d].pivot(index="gene", columns="seed", values="z").reindex(columns=SEEDS)
         for d in longs}

    def zm(direction: str, genes: list[str]) -> np.ndarray:
        return Z[direction].reindex([g for g in genes if g in Z[direction].index]).to_numpy()

    # ---------------------------------------------------------------- 7. 止める検査（A-4・A-5）
    checks: dict = {}
    cand_ada = zm("adaptive", G["候補"])
    S_ada = cand_ada.mean(axis=0)
    ref = np.array(PLAN["再現確認_参照"])
    a4_diff = float(np.max(np.abs(S_ada - ref)))
    checks["A-4 再現確認（活性化の方向）"] = {
        "再計算した候補8のセット平均": [round(float(x), 4) for x in S_ada],
        "criteria_011.json の参照値": PLAN["再現確認_参照"],
        "最大絶対差": a4_diff, "許容差": PLAN["再現確認_許容差"],
        "pass": bool(a4_diff <= PLAN["再現確認_許容差"]),
    }
    if smoke:
        checks["A-4 再現確認（活性化の方向）"]["注（試し実行）"] = (
            "試し実行は標準化の対照が 510 個ではなく一部（43 個）なので、参照値とは一致しない。"
            "A-4 の合否は本実行でのみ判定する。")
    L(f"A-4 再現確認（活性化の方向の候補 8 セット平均）: {[round(float(x), 4) for x in S_ada]}、"
      f"参照値との最大絶対差 {a4_diff:.2e}（許容 1e-3）、pass={checks['A-4 再現確認（活性化の方向）']['pass']}"
      + ("（試し実行。対照プールが違うので参照値とは一致しない。判定は本実行で行う）" if smoke else ""))
    if not checks["A-4 再現確認（活性化の方向）"]["pass"] and not smoke:
        out["A チェック"] = checks
        (outdir / f"criteria_013{suf}.json").write_text(
            json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
        raise SystemExit("**止める**(A-4): 活性化の方向の再現確認が合わない（B-1）")

    # A-5 値（NaN・無限大なし）
    finite = {}
    for d in ("adaptive", "degenerate"):
        arr = Z[d].to_numpy()
        finite[d] = {"NaN": int(np.isnan(arr).sum()), "無限大": int(np.isinf(arr).sum()),
                     "行数": int(arr.shape[0]), "列数": int(arr.shape[1])}
    checks["A-5 値（標準得点）"] = {"方向ごと": finite,
                                 "pass": all(v["NaN"] == 0 and v["無限大"] == 0 for v in finite.values())}
    L(f"A-5 標準得点に NaN・無限大なし: {[f'{d}={finite[d]}' for d in finite]}")
    if not checks["A-5 値（標準得点）"]["pass"]:
        out["A チェック"] = checks
        (outdir / f"criteria_013{suf}.json").write_text(
            json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
        raise SystemExit("**止める**(A-5): 標準得点に NaN・無限大がある")
    out["A チェック"] = checks

    # ---------------------------------------------------------------- 8. U-1 主判定（plan.md §5.3）
    def judge(direction: str, genes: list[str], key: int, abs_thr: float | None) -> dict:
        obs = zm(direction, genes).mean(axis=0)
        thr = np.quantile(set_null(zm(direction, [g for g in pool]), len(genes), N_DRAW,
                                   np.random.default_rng(key), replace=True), 0.99, axis=0)
        ci = bool(np.all(obs > thr))
        cii = bool(np.all(np.abs(obs) >= abs_thr)) if abs_thr is not None else None
        return {"セットの大きさ": len(genes), "セット平均の標準得点": [round(float(x), 4) for x in obs],
                "偽セットの99パーセンタイル": [round(float(x), 4) for x in thr],
                "条件(i) 3 群すべてで 99 パーセンタイル超": ci,
                "条件(ii) 3 群すべてで 絶対値 >= 閾値": cii,
                "閾値（絶対値）": abs_thr, "乱数の鍵": key,
                "合格": bool(ci and (cii is not False)) if abs_thr is not None else bool(ci)}

    u1a = judge("adaptive", G["候補"], PLAN["乱数の鍵"]["候補の偽セット"], PLAN["閾値_候補"])
    u1b = judge("degenerate", G["候補"], PLAN["乱数の鍵"]["候補の偽セット"], PLAN["閾値_候補"])
    u1c = judge("degenerate", G["陽性対照_退化30"], PLAN["乱数の鍵"]["陽性対照の偽セット"], None)
    out["4. U-1 主判定（config の decision をそのまま適用）"] = {
        "合格の定義（候補8個）": "偽セット 1,000 組（対照 510 個から復元ありで 8 個）の 3 群それぞれの "
                              "99 パーセンタイル超 かつ 3 群すべてで絶対値 1.59 以上",
        "合格の定義（陽性対照30個）": "偽セット 1,000 組（30 個）の 3 群それぞれの 99 パーセンタイル超のみ"
                                 "（絶対値の条件は課さない。PI 判断 2026-10-11）",
        "U-1a 候補8（活性化した間質の方向。再現確認）": u1a,
        "U-1b 候補8（退化した状態の方向。主対象）": u1b,
        "U-1c 陽性対照30（退化した状態の方向）": u1c,
    }
    if smoke:
        out["4. U-1 主判定（config の decision をそのまま適用）"]["注（試し実行）"] = (
            "**試し実行の値であり判定ではない。** 標準化の対照が 510 個ではなく一部（43 個）で、"
            "遺伝子も 48 個だけなので、値は本実行と比較できない。")
    L(f"U-1a 候補8・活性化の方向: セット平均 {u1a['セット平均の標準得点']}、99%点 "
      f"{u1a['偽セットの99パーセンタイル']}、合格={u1a['合格']}")
    L(f"U-1b 候補8・退化の方向: セット平均 {u1b['セット平均の標準得点']}、99%点 "
      f"{u1b['偽セットの99パーセンタイル']}、合格={u1b['合格']}")
    L(f"U-1c 陽性対照30・退化の方向: セット平均 {u1c['セット平均の標準得点']}、99%点 "
      f"{u1c['偽セットの99パーセンタイル']}、合格={u1c['合格']}")

    # 陽性対照 71（活性化の方向。EXP-001 の R-6 の再現。記述。判定に使わない）
    pc_ada = zm("adaptive", G["陽性対照_活性化71"]).mean(axis=0)
    out["5. 陽性対照71（活性化の方向。記述。判定に使わない）"] = {
        "セット平均の標準得点": [round(float(x), 4) for x in pc_ada],
        "criteria_011.json の参照値（R-6）": [1.449, 1.398, 1.262],
        "読み方": "EXP-001 の R-6 と同じ量。U-1 の判定には使わない（plan.md §5.3）",
    }
    L(f"陽性対照71・活性化の方向（記述）: {[round(float(x), 4) for x in pc_ada]}（EXP-001 の参照 1.449/1.398/1.262）")

    # ---------------------------------------------------------------- 9. G-c（plan.md §5.4）
    ctrl = [g for g in G["対照500"] if g in Z["degenerate"].index]
    zbar = Z["degenerate"].loc[ctrl].mean(axis=1)
    n_exp = long_deg.pivot(index="gene", columns="seed", values="n_expressed").reindex(
        columns=SEEDS).loc[ctrl].mean(axis=1)
    ll_med_by_gene = expressed_length_median(c, ctrl, L)
    llm = ll_med_by_gene.groupby("gene").ll_exp_med.mean().reindex(ctrl)
    lp_med = long_deg.pivot(index="gene", columns="seed", values="lp_med").reindex(
        columns=SEEDS).loc[ctrl].mean(axis=1)
    r_exp = float(stats.spearmanr(zbar, n_exp).statistic)
    r_len = float(stats.spearmanr(zbar, llm).statistic)
    gc_pass = abs(r_exp) <= PLAN["G-c_閾値"] and abs(r_len) <= PLAN["G-c_閾値"]
    layer_corr = {}
    for lname in sorted({G["層"][g] for g in ctrl}):
        gg = [g for g in ctrl if G["層"][g] == lname]
        layer_corr[lname] = {
            "n": len(gg),
            "z と発現細胞数": spearman_rec(zbar.loc[gg], n_exp.loc[gg], "z", "発現細胞数"),
            "z と入力の長さの中央": spearman_rec(zbar.loc[gg], llm.loc[gg], "z", "入力の長さ"),
        }
    out["6. G-c（対照500 の退化の方向の標準得点と交絡）"] = {
        "対象": {"遺伝子数": len(ctrl), "層": {k: sum(1 for g in ctrl if G["層"][g] == k)
                                          for k in sorted({G["層"][g] for g in ctrl})}},
        "ρ_exp（発現細胞数との Spearman）": round(r_exp, 4),
        "ρ_len（入力の長さの中央値との Spearman）": round(r_len, 4),
        "閾値": PLAN["G-c_閾値"],
        "合格条件": "|ρ_exp| <= 0.20 かつ |ρ_len| <= 0.20（3 群平均の値で判定）",
        "合格": bool(gc_pass),
        "不合格のとき": "U-R4（判定不能）",
        "記録（判定に使わない）": {
            "ρ_lp（元の順位との Spearman）": spearman_rec(zbar, lp_med, "z", "元の順位の中央"),
            "層ごと": layer_corr,
            "発現 0 で補完した 遺伝子×群 の数": int(ll_med_by_gene["補完したか"].sum()),
            "入力の長さの中央値の記述": {
                "遺伝子ごとの 3 群平均の中央": round(float(llm.median()), 2),
                "遺伝子ごとの 3 群平均の最小": round(float(llm.min()), 2),
                "遺伝子ごとの 3 群平均の最大": round(float(llm.max()), 2)},
        },
    }
    L(f"G-c: |ρ_exp|={abs(r_exp):.4f}、|ρ_len|={abs(r_len):.4f}（閾値 0.20）、合格={gc_pass}")

    # ---------------------------------------------------------------- 10. 反証条件 U-R1〜U-R5（plan.md §5.5）
    u_r1 = bool(u1b["合格"])                       # 退化の方向でも候補 8 が合格する
    u_r2 = not bool(u1a["合格"])                   # 活性化の方向で候補 8 が合格しない
    u_r3 = not bool(u1c["合格"])                   # 退化の方向の陽性対照が合格しない
    u_r4 = not bool(gc_pass)                       # 対照で |ρ| > 0.20
    if u_r4:
        first = "U-R4"
    elif u_r3:
        first = "U-R3"
    elif u_r1:
        first = "U-R1"
    elif u_r2:
        first = "U-R2"
    else:
        first = "該当なし（U-1 の予測どおり。plan §5.5 評価順 6）"
    out["7. 反証条件 U-R1〜U-R5（config をそのまま適用）"] = {
        "評価順（固定）": "1 再現確認 A-4 → 2 G-c（U-R4）→ 3 陽性対照（U-R3）→ 4 候補8の退化（U-R1）→ 5 U-R2",
        "U-R1 候補8が退化した状態の方向でも合格する（2 つの方向が分離しない）": u_r1,
        "U-R2 候補8が活性化した間質の方向でも合格しない": u_r2,
        "U-R3 退化した状態の方向の陽性対照が合格しない（判定不能）": u_r3,
        "U-R4 対照遺伝子で |ρ| > 0.20（判定不能）": u_r4,
        "U-R5 上のどれにも当たらない（結果が事前登録の型に当てはまらない）": None,
        "評価順で最初に当たったもの": first,
        "注": "plan.md §5.5 の表の U-R5 は「保留」と書き、同 §5.5 の評価順 6 は「どれにも当たらない → "
              "U-1 の予測どおり（PI に提案）」と書いている。**この 2 つの文言は食い違う。** ここでは"
              "規則を後から作らず、上に書いたとおり『評価順で最初に当たったもの』と各条件の真偽だけを出す。"
              "U-R1〜U-R4 がすべて偽なら『U-1 の予測どおり』であり、その場合は U-R5 を true とも false とも"
              "確定しない（判定は PI＝ゲート G3）。",
    }
    out["7b. U-1 の予測との一致（提案。確定は PI）"] = {
        "予測": "U-1a 合格 かつ U-1b 不合格 かつ U-1c 合格",
        "U-1a 合格": bool(u1a["合格"]), "U-1b 不合格": not bool(u1b["合格"]), "U-1c 合格": bool(u1c["合格"]),
        "予測どおり": bool(u1a["合格"] and not u1b["合格"] and u1c["合格"]),
        "G-c 合格（U-R4 に当たらない）": bool(gc_pass),
    }
    if smoke:
        out["7. 反証条件 U-R1〜U-R5（config をそのまま適用）"]["注（試し実行）"] = (
            "**試し実行の値であり判定ではない。** 対照プールと遺伝子数が本実行と違う。")
    L(f"反証条件: U-R1={u_r1} U-R2={u_r2} U-R3={u_r3} U-R4={u_r4} / 最初に当たったもの: {first}")

    # ---------------------------------------------------------------- 11. 退化の方向の順位表（探索。CLAUDE.md §3.5）
    Zr = Z["degenerate"].rename(columns={s: f"z_{s}" for s in SEEDS})
    rank = (Zr.assign(z_mean3=Z["degenerate"].mean(axis=1), z_min3=Z["degenerate"].min(axis=1))
            .reset_index().rename(columns={"index": "gene"})
            .assign(symbol=lambda d: d["gene"].map(G["symbol"]),
                    注="探索（データを見て候補を探すための順位表。事前登録した確認的検定ではない。"
                       "確認は別データか別 seed で行う。CLAUDE.md §3.5）")
            .sort_values("z_mean3", ascending=False).reset_index(drop=True))
    rank.insert(0, "rank", np.arange(1, len(rank) + 1))
    rp = outdir / f"gene_rank_013{suf}.csv"
    rank.to_csv(rp, index=False)
    out["8. 退化した状態の方向の順位表（探索）"] = {
        "行数": len(rank), "保存先": str(rp.relative_to(ROOT)),
        "注": "**探索である。** 事前登録した確認的検定ではない。ファイル内の `注` 列にも明記した",
        "列": list(rank.columns),
        "上位10": rank.head(10)[["symbol", "z_mean3", "z_min3"]].round(4).to_dict(orient="records"),
        "下位10": rank.tail(10)[["symbol", "z_mean3", "z_min3"]].round(4).to_dict(orient="records"),
    }
    L(f"順位表（探索）を書いた: {rp.relative_to(ROOT)}（{len(rank)} 行）")

    # ---------------------------------------------------------------- 12. 出力（A-6）・README・所要時間
    out["decision（config.yaml の全鍵）"] = dec
    out["注"] = ("判定（支持・否定・判定不能。ゲート G3）は PI が行う。この道具は規則を機械的に当てた"
                "結果と提案まで。hypotheses/tracker.md は変更しない。コンピュータ上の遺伝子摂動の"
                "結果であり、実験的に確認された生物学的事実ではない。盲検中の遺伝子は摂動も解析も"
                "していない（値は数えない）。")
    dur = time.time() - t0
    out["実行"] = {"開始": _dt.datetime.fromtimestamp(t0).astimezone().isoformat(),
                 "所要秒": round(dur, 1), "GPU 時間": 0.0, "device": cfg["compute"]["device"],
                 "dtype": cfg["compute"]["dtype"]}
    cj = outdir / f"criteria_013{suf}.json"
    cj.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    reloaded = json.loads(cj.read_text())
    nested = reloaded.get("decision（config.yaml の全鍵）") or {}
    a6 = {"decision の全鍵を持つ": set(nested.keys()) == set(dec.keys()),
          "decision の鍵": sorted(dec.keys()), "書き込めた": True, "保存先": str(cj.relative_to(ROOT)),
          "置き場所": "「decision（config.yaml の全鍵）」の節に、config.yaml の decision: の全鍵をそのまま入れた"}
    a6["pass"] = bool(a6["decision の全鍵を持つ"] and a6["書き込めた"])
    L(f"A-6 結果ファイル: {cj.relative_to(ROOT)} を書いた。decision: の全鍵を持つ={a6['pass']}")
    write_readme(outdir, suf, out)
    L(f"説明を書いた: {(outdir / 'README.md').relative_to(ROOT)}")

    # ---------------------------------------------------------------- 13. A-1（出力・ログ・一覧に盲検が 0 件。値は数えない）
    sym_blk = {"ENSG00000141510": None, "ENSG00000124762": None, "ENSG00000111321": None,
               "ENSG00000118785": None, "ENSG00000148346": None}
    names_blk = ["TP53", "CDKN1A", "LTBR", "SPP1", "LCN2"]
    targets = [cj, rp, outdir / f"gene_group_013{suf}.csv", outdir / "README.md", LOG]
    a1 = {"調べたファイル": [str(p.relative_to(ROOT)) for p in targets if p.exists()], "件数": 0}
    for p in targets:
        if not p.exists():
            continue
        text = p.read_text(errors="replace")
        hits = sum(text.count(k) for k in sym_blk) + sum(text.count(k) for k in names_blk)
        if hits:
            a1["件数"] += hits
            a1.setdefault("当たったファイル", []).append(str(p.relative_to(ROOT)))
    a1["pass"] = a1["件数"] == 0
    L(f"A-1 盲検の混入（出力・ログ・一覧）: {a1['件数']} 件、pass={a1['pass']}")
    if not a1["pass"]:
        raise SystemExit("**止める**(A-1): 盲検が出力・ログ・一覧に現れた（B-2）")

    # ---------------------------------------------------------------- 14. 検査の結果を書き戻す（A-6 の再確認）
    out["A チェック"]["A-1 盲検の混入"] = a1
    out["A チェック"]["A-6 結果ファイル"] = a6
    out["A チェック"]["A-2 材料"] = {"pass": mats["pass"]}
    out["A チェック"]["A-3 遺伝子集合"] = {"pass": G["pass"]}
    cj.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    re2 = json.loads(cj.read_text())
    a6b = set((re2.get("decision（config.yaml の全鍵）") or {}).keys()) == set(dec.keys())
    L(f"A-6 追記後も criteria_013{suf}.json は decision: の全鍵を持つ={a6b}。所要 {dur:.1f} 秒")
    L("終了（終了コード 0）")


def write_readme(outdir: Path, suf: str, out: dict) -> None:
    u = out["4. U-1 主判定（config の decision をそのまま適用）"]
    gc = out["6. G-c（対照500 の退化の方向の標準得点と交絡）"]
    a4 = out["A チェック"]["A-4 再現確認（活性化の方向）"]
    txt = f"""# EXP-20261007-003 段階 2 — 出力（{'試し実行' if suf else '本実行'}）

- 仮説: HYP-20261007-02 段階 2。**保存済み埋め込みの再射影のみ。摂動なし・Geneformer 不使用**（GPU 0 時間）
- 所要: {out['実行']['所要秒']} 秒。数値型 fp32、device cpu
- **判定（支持／否定／判定不能）の確定は PI（ゲート G3）。ここは値と、規則を機械的に当てた結果まで。**
- ゲート G2 不要（PI 判断 2026-10-11）。`Approved-by: PI 2026-10-11`
{''
 + '''- **これは試し実行である。値は判定ではない。** 保存済み埋め込みの一部（遺伝子 48 個）だけを見せ、
  標準化の対照も 510 個ではなく一部（43 個）で回した。**下の U-1a/U-1b/U-1c と反証条件は本実行と
  比較できない。** A-4（再現確認）の合否は本実行でのみ判定する。目的は評価コードが最後まで通ることの
  確認だけ（guardrail 2026-10-01）。
''' if suf else ''}

## 何を出したか

| ファイル | 中身 |
|---|---|
| `criteria_013{suf}.json` | 判定値と、判定に使ったすべての数値（config の decision: の全鍵を含む） |
| `gene_rank_013{suf}.csv` | 587 遺伝子の**退化した状態の方向**の標準得点の順位表（**探索**であることを `注` 列に明記） |
| `gene_group_013{suf}.csv` | 遺伝子 × 群 × 方向の `D_raw`・`z`・`n_expressed`・`lp_med`・`ll_med` など（補助） |
| `README.md` | この説明 |
| `run.log` | 実行の記録 |

## 主な数値

- 再現確認 A-4（活性化した間質の方向の候補 8 セット平均）: {a4['再計算した候補8のセット平均']}、
  参照（criteria_011.json）{a4['criteria_011.json の参照値']}、最大絶対差 {a4['最大絶対差']:.2e}、合格={a4['pass']}
- U-1a 候補 8・活性化の方向: セット平均 {u['U-1a 候補8（活性化した間質の方向。再現確認）']['セット平均の標準得点']}、
  99 パーセンタイル {u['U-1a 候補8（活性化した間質の方向。再現確認）']['偽セットの99パーセンタイル']}、
  合格={u['U-1a 候補8（活性化した間質の方向。再現確認）']['合格']}
- U-1b 候補 8・**退化した状態の方向**（主対象）: セット平均 {u['U-1b 候補8（退化した状態の方向。主対象）']['セット平均の標準得点']}、
  99 パーセンタイル {u['U-1b 候補8（退化した状態の方向。主対象）']['偽セットの99パーセンタイル']}、
  合格={u['U-1b 候補8（退化した状態の方向。主対象）']['合格']}
- U-1c 陽性対照 30・退化の方向: セット平均 {u['U-1c 陽性対照30（退化した状態の方向）']['セット平均の標準得点']}、
  99 パーセンタイル {u['U-1c 陽性対照30（退化した状態の方向）']['偽セットの99パーセンタイル']}、
  合格={u['U-1c 陽性対照30（退化した状態の方向）']['合格']}
- G-c: |ρ_exp|={abs(gc['ρ_exp（発現細胞数との Spearman）']):.4f}、|ρ_len|={abs(gc['ρ_len（入力の長さの中央値との Spearman）']):.4f}（閾値 0.20）、合格={gc['合格']}
- 反証条件で最初に当たったもの: **{out['7. 反証条件 U-R1〜U-R5（config をそのまま適用）']['評価順で最初に当たったもの']}**

## 注意

- in silico の計算である。「示唆される」までで、実験的事実としては書かない。
- **`gene_rank_013{suf}.csv` は探索である。** 事前登録した確認的検定ではない（CLAUDE.md §3.5）。
- 盲検中の 5 遺伝子の名前・値は、出力・ログのどこにも書かない（plan.md §6.1・§7 A-1）。
- 判定基準・乱数の鍵・遺伝子集合は plan.md と config.yaml に固定した値のまま。実行時に変えていない。
- 限界は plan.md §11（とくに §11-6 陽性対照 30 個のうち 23 個が標準化プールと重なる、§11-7 2 つの方向は
  埋め込み上ほぼ同じ向きにある、§11-9 評価できるのは保存済み 587 遺伝子だけ）。
"""
    (outdir / "README.md").write_text(txt)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true", help="保存済み埋め込みの一部で評価コードを通す（判定ではない）")
    main(smoke=ap.parse_args().smoke)

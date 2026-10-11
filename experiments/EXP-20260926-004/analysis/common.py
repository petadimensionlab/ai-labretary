"""EXP-20260926-004 共通関数（CPU のみ）。EXP-002/003 の a2 と判定（C1〜C3, C3-panel）を再実装し、
細胞単位・donor 単位の両方に使う。a2 の定義は experiments/EXP-20260926-002/analysis/reanalysis_cosine.py の
corrected() の m_a2 と同一（長さ 5 分位ビン［null 細胞で境界］の null 平均を細胞ごとに引く。null 遺伝子自身は
leave-one-gene-out）。一致は donor_level.py の先頭で EXP-003 の gene_seed_metrics.csv と照合する。
盲検: 主要対象・参照（BLIND）は読み込み直後に除外する。
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BLIND = {"ENSG00000141510", "ENSG00000124762", "ENSG00000111321", "ENSG00000118785", "ENSG00000148346"}
POS = {"ENSG00000162692": "VCAM1", "ENSG00000113249": "HAVCR1"}
HK = {"ENSG00000075624": "ACTB", "ENSG00000111640": "GAPDH", "ENSG00000196262": "PPIA"}
PANEL = {"ENSG00000090661": "CERS4", "ENSG00000147687": "TATDN1", "ENSG00000081870": "IFT25",
         "ENSG00000101558": "VAPA", "ENSG00000054148": "PHPT1", "ENSG00000166295": "ANAPC16",
         "ENSG00000104529": "EEF1D"}
RP = {"ENSG00000110700": "RPS13", "ENSG00000144713": "RPL32"}
DE = {"ENSG00000110841": "PPFIBP1", "ENSG00000131016": "AKAP12", "ENSG00000198832": "SELENOM",
      "ENSG00000161249": "DMKN", "ENSG00000006747": "SCIN", "ENSG00000165887": "ANKRD2",
      "ENSG00000140450": "ARRDC4", "ENSG00000127415": "IDUA", "ENSG00000179542": "SLITRK4",
      "ENSG00000145536": "ADAMTS16"}
DE_BAND = {**{g: "VCAM1" for g in list(DE)[:5]}, **{g: "HAVCR1" for g in list(DE)[5:]}}
SYM = {**POS, **HK, **PANEL, **RP, **DE}
CRIT = {"C1_q": 0.90, "min_seeds": 2, "C2_frac": 15 / 18, "C3_range": (0.10, 0.90), "panel_min_genes": 5}


def load_cells(res: Path, logits: bool = False) -> pd.DataFrame:
    """results ディレクトリの cells.csv（全 run）。role は negatives.json で補う（seed を別プロセスで回した場合の対策）。
    logits=True なら genes/<g>/emb/logits_{orig,pert}.npy から dP（PT_healthy の確率変化）を付ける。"""
    c = pd.read_csv(res / "cells.csv", dtype={"donor_id": str})
    c = c[~c.gene.isin(BLIND)].reset_index(drop=True)
    assert not set(c.gene) & BLIND
    for run in c.run.unique():
        nj = res / run / "negatives.json"
        if nj.exists():
            negs = set(json.loads(nj.read_text()))
            c.loc[(c.run == run) & c.gene.isin(negs) & (c.role == "screen"), "role"] = "negative"
    if logits:
        parts = []
        for (run, gene), _ in c.groupby(["run", "gene"], sort=False):
            gd = res / run / "genes" / gene
            per = pd.read_csv(gd / "cells.csv", dtype={"donor_id": str})
            lo = np.load(gd / "emb" / "logits_orig.npy").astype(np.float64)
            lp = np.load(gd / "emb" / "logits_pert.npy").astype(np.float64)
            assert len(lo) == len(lp) == len(per)
            p = lambda z: np.exp(z[:, 1] - np.logaddexp(z[:, 0], z[:, 1]))  # softmax(z)[1]
            parts.append(pd.DataFrame({"run": run, "gene": gene, "cell_uid": per.cell_uid.to_numpy(),
                                       "dP": p(lp) - p(lo)}))
        d = pd.concat(parts, ignore_index=True)
        n0 = len(c)
        c = c.merge(d, on=["run", "gene", "cell_uid"], how="left", validate="one_to_one")
        assert len(c) == n0 and c.dP.notna().all()
    return c


def bin_correct(run_df: pd.DataFrame, base: str, nbins: int = 5) -> np.ndarray:
    """1 run の細胞表で base 列を長さビンごとの null 平均で補正（a2 と同じ手順）。"""
    df = run_df
    out = np.full(len(df), np.nan)
    null_genes = sorted(df.loc[df.role == "negative", "gene"].unique())
    groups = [(df.gene == g, (df.role == "negative") & (df.gene != g)) for g in null_genes]
    groups.append((df.role != "negative", df.role == "negative"))
    qs = np.linspace(0, 1, nbins + 1)[1:-1]
    for tgt, fit in groups:
        F, T = df[fit], df[tgt]
        if len(T) == 0:
            continue
        edges = np.quantile(F.length, qs)
        bf, bt = np.digitize(F.length, edges), np.digitize(T.length, edges)
        mu = F.groupby(bf)[base].mean()
        out[tgt.to_numpy()] = T[base].to_numpy() - mu.reindex(bt).to_numpy()
    return out


def add_metrics(c: pd.DataFrame) -> pd.DataFrame:
    parts = []
    for _, s in c.groupby("run", sort=True):
        s = s.copy()
        s["m_raw"] = s.shift_goal
        s["m_a2"] = bin_correct(s, "shift_goal")
        if "dP" in s:
            s["dP_goal_a2"] = bin_correct(s, "dP")
        parts.append(s)
    return pd.concat(parts, ignore_index=True)


def null_pct(x: float, null) -> float:
    return float((np.asarray(null) < x).mean())


def criteria(g: pd.DataFrame, col: str = "v") -> dict:
    """g: run, gene, role, <col>（遺伝子×seed のスコア）。EXP-002/003 と同じ C1〜C3・C3-panel。"""
    runs = sorted(g.run.unique())
    need = CRIT["min_seeds"] if len(runs) >= 3 else len(runs)
    pct = {}
    nfn = []
    for run, s in g.groupby("run"):
        null = s[s.role == "negative"][col].to_numpy()
        nfn.append(float((null < 0).mean()))
        v = s.set_index("gene")[col]
        for gene in v.index:
            if gene in SYM:
                pct[(run, gene)] = null_pct(v[gene], null)
    res = {"n_seeds": len(runs), "null_frac_negative": round(float(np.mean(nfn)), 3)}
    res["pct"] = {SYM[k[1]] + "@" + k[0]: round(v, 3) for k, v in pct.items()}
    c1 = {POS[p]: sum(pct.get((r, p), 0) >= CRIT["C1_q"] for r in runs) for p in POS}
    lo, hi = CRIT["C3_range"]
    c3 = {HK[h]: sum(lo <= pct.get((r, h), -1) <= hi for r in runs) for h in HK}
    res["C1_pct"] = {POS[p]: [round(pct.get((r, p), np.nan), 3) for r in runs] for p in POS}
    res["C1"] = all(v >= need for v in c1.values())
    res["C3_pct"] = {HK[h]: [round(pct.get((r, h), np.nan), 3) for r in runs] for h in HK}
    res["C3"] = all(v >= need for v in c3.values())
    wins, tot = 0, 0
    for run, s in g.groupby("run"):
        v = s.set_index("gene")[col]
        for p in POS:
            for h in HK:
                if p in v and h in v:
                    wins += v[p] > v[h]; tot += 1
    res["C2_wins"] = f"{wins}/{tot}"
    res["C2"] = tot > 0 and wins >= CRIT["C2_frac"] * tot
    per = {PANEL[p]: sum(lo <= pct.get((r, p), -1) <= hi for r in runs) for p in PANEL}
    ok = sum(v >= need for v in per.values())
    res["C3_panel_n_ok"] = f"{ok}/{len(PANEL)}"
    res["C3_panel"] = ok >= CRIT["panel_min_genes"]
    res["all_C1_C3_panel"] = res["C1"] and res["C2"] and res["C3"] and res["C3_panel"]
    return res

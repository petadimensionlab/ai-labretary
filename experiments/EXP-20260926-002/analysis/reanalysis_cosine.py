"""EXP-20260926-002 step 1: EXP-20260926-001 の cells.csv（コサインのみ）で削除ドリフト補正指標を比較する。
盲検: 主要対象 TP53/CDKN1A/LTBR と参照 SPP1/LCN2 は読み込み直後に除外する（plan.md §6.0）。
実行: .venv/bin/python experiments/EXP-20260926-002/analysis/reanalysis_cosine.py
出力: experiments/EXP-20260926-002/results/reanalysis/{criteria.csv,gene_seed_metrics.csv,criteria.json}
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "experiments/EXP-20260926-001"
OUT = ROOT / "experiments/EXP-20260926-002/results/reanalysis"
BLIND = {"ENSG00000141510", "ENSG00000124762", "ENSG00000111321", "ENSG00000118785", "ENSG00000148346"}
POS = {"ENSG00000162692": "VCAM1", "ENSG00000113249": "HAVCR1"}
HK = {"ENSG00000075624": "ACTB", "ENSG00000111640": "GAPDH", "ENSG00000196262": "PPIA"}
START = {"tls_high": "aPT_CKD_TLShigh", "tls_low": "aPT_CKD_TLSlow"}
N_SPLIT = 20


def token_positions(cells: pd.DataFrame) -> pd.Series:
    """削除した遺伝子の入力内の位置（0 = <cls> の次ではなく input_ids の index）。tokenized dataset は読み取りのみ。"""
    import pickle
    from datasets import load_from_disk
    from geneformer import TOKEN_DICTIONARY_FILE
    tokd = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
    ds = load_from_disk(str(ROOT / "experiments/EXP-20260926-002/work/tokenized.dataset"))  # EXP-001 の複製（EXP-001 に cache を書かない）
    uids = sorted(cells.cell_uid.unique())
    ids = dict(zip(uids, ds.select(uids)["input_ids"]))
    pos = [ids[u].index(tokd[g]) for u, g in zip(cells.cell_uid, cells.gene)]
    return pd.Series(pos, index=cells.index)


def ols_fit_predict(Xfit, yfit, Xall):
    Xf = np.column_stack([np.ones(len(Xfit)), Xfit])
    Xa = np.column_stack([np.ones(len(Xall)), Xall])
    beta, *_ = np.linalg.lstsq(Xf, yfit, rcond=None)
    return Xa @ beta


def design(df: pd.DataFrame, kind: str, donors: list[str]) -> np.ndarray:
    ll = np.log(df.length.to_numpy(float))
    if kind == "len_donor":
        D = np.column_stack([(df.donor_id == d).to_numpy(float) for d in donors[1:]])
        return np.column_stack([ll, D])
    if kind == "len_pos":
        rp = df.pos.to_numpy(float) / df.length.to_numpy(float)
        return np.column_stack([ll, rp, np.log1p(df.pos.to_numpy(float))])
    if kind == "start_alt":
        return df[["shift_start", "alt"]].to_numpy(float)
    raise ValueError(kind)


def corrected(run_df: pd.DataFrame) -> pd.DataFrame:
    """1 run（=1 seed）の細胞表に候補指標を付ける。null から推定する補正は、null 遺伝子自身には
    leave-one-gene-out（その遺伝子を除いた null で推定）を使う。"""
    df = run_df.copy()
    df["m_raw"] = df.shift_goal
    df["m_d"] = df.shift_goal - df.shift_start
    donors = sorted(df.donor_id.unique())
    null_genes = sorted(df.loc[df.role == "negative", "gene"].unique())
    groups = [(g, df.gene == g, (df.role == "negative") & (df.gene != g)) for g in null_genes]
    groups.append((None, df.role != "negative", df.role == "negative"))
    cols = {k: np.full(len(df), np.nan) for k in ["m_a2", "m_c", "m_a3", "m_a4", "m_b", "m_d2", "m_d4"]}
    for _, tgt, fit in groups:
        F, T = df[fit], df[tgt]
        # 長さ 5 分位ビン（fit 側の null 細胞で境界を決める）
        edges = np.quantile(F.length, [0.2, 0.4, 0.6, 0.8])
        bf, bt = np.digitize(F.length, edges), np.digitize(T.length, edges)
        for col, base in [("m_a2", "shift_goal"), ("m_d2", "m_d")]:
            mu = F.groupby(bf)[base].mean()
            cols[col][tgt.to_numpy()] = T[base].to_numpy() - mu.reindex(bt).to_numpy()
        mu, sd = F.groupby(bf).shift_goal.mean(), F.groupby(bf).shift_goal.std()
        cols["m_c"][tgt.to_numpy()] = (T.shift_goal.to_numpy() - mu.reindex(bt).to_numpy()) / sd.reindex(bt).to_numpy()
        for col, kind, base in [("m_a3", "len_donor", "shift_goal"), ("m_a4", "len_pos", "shift_goal"),
                                ("m_b", "start_alt", "shift_goal"), ("m_d4", "len_pos", "m_d")]:
            pred = ols_fit_predict(design(F, kind, donors), F[base].to_numpy(float), design(T, kind, donors))
            cols[col][tgt.to_numpy()] = T[base].to_numpy() - pred
    for k, v in cols.items():
        df[k] = v
    return df


METRICS = {"m_raw": "raw shift_goal", "m_d": "d: goal−start", "m_a2": "a2: −null(長さビン)",
           "m_c": "c: z(長さビン)", "m_a3": "a3: −null(log長+donor)", "m_a4": "a4: −null(log長+削除位置)",
           "m_b": "b: 残差(start,alt)", "m_d2": "d2: d−null(長さビン)", "m_d4": "d4: d−null(log長+削除位置)"}


def null_pct(x, null):
    return float((np.asarray(null) < x).mean())


def evaluate(cells: pd.DataFrame, det: dict, metric: str) -> dict:
    g = cells.groupby(["run", "gene", "role"]).agg(v=(metric, "mean"), L=("length", "mean")).reset_index()
    g["det"] = g.gene.map(det)
    res = {"metric": metric}
    # C1 / C3 / 付記
    pct = {}
    for run, s in g.groupby("run"):
        null = s[s.role == "negative"].set_index("gene").v
        for gene in list(POS) + list(HK):
            if gene in set(s.gene):
                pct[(run, gene)] = null_pct(s.set_index("gene").v[gene], null)
        res.setdefault("null_frac_negative", []).append(float((null < 0).mean()))
    runs = sorted(g.run.unique())
    need = 2 if len(runs) >= 3 else len(runs)  # 3 seed 中 2（phase 0 の 1 seed では 1/1）
    res["n_seeds"] = len(runs)
    c1 = {POS[p]: sum(pct.get((r, p), 0) >= 0.90 for r in runs) for p in POS}
    c3 = {HK[h]: sum(0.10 <= pct.get((r, h), -1) <= 0.90 for r in runs) for h in HK}
    res["C1_seeds"] = c1
    res["C1_pct_mean"] = {POS[p]: round(np.mean([pct[(r, p)] for r in runs]), 3) for p in POS}
    res["C1"] = all(v >= need for v in c1.values())
    res["C3_seeds"] = c3
    res["C3_pct_mean"] = {HK[h]: round(np.mean([pct[(r, h)] for r in runs]), 3) for h in HK}
    res["C3"] = all(v >= need for v in c3.values())
    # C2
    wins, tot, margin = 0, 0, []
    for run, s in g.groupby("run"):
        v = s.set_index("gene").v
        nsd = s[s.role == "negative"].v.std()
        for p in POS:
            for h in HK:
                wins += v[p] > v[h]; tot += 1; margin.append((v[p] - v[h]) / nsd)
    res["C2_wins"] = f"{wins}/{tot}"
    res["C2_margin_nullSD"] = round(float(np.mean(margin)), 2)
    res["C2"] = wins >= 15 / 18 * tot
    res["null_frac_negative"] = round(float(np.mean(res["null_frac_negative"])), 2)
    # C4 (i) split-half（donor を 2 群に分け、25 遺伝子の順位一致）
    rng = np.random.default_rng(0)
    rhos = []
    for run, s in cells.groupby("run"):
        donors = np.array(sorted(s.donor_id.unique()))
        for _ in range(N_SPLIT):
            half = set(rng.permutation(donors)[: len(donors) // 2])
            h = s.donor_id.isin(half)
            a = s[h].groupby("gene")[metric].mean(); b = s[~h].groupby("gene")[metric].mean()
            idx = a.index.intersection(b.index)
            rhos.append(stats.spearmanr(a[idx], b[idx]).statistic)
    res["C4i_splithalf_rho"] = round(float(np.nanmean(rhos)), 3)
    # C4 (ii) seed 間, 固定対照 5 遺伝子
    fixed = list(POS) + list(HK)
    piv = g[g.gene.isin(fixed)].pivot(index="gene", columns="run", values="v")
    rr = [stats.spearmanr(piv[a], piv[b]).statistic for i, a in enumerate(runs) for b in runs[i + 1:]]
    res["C4ii_seed_rho_fixed5"] = round(float(np.mean(rr)), 3) if rr else None
    res["C4"] = res["C4i_splithalf_rho"] >= 0.7 and (not rr or res["C4ii_seed_rho_fixed5"] >= 0.7)
    # C5
    r_len = stats.spearmanr(g.L, g.v).statistic
    r_det = stats.spearmanr(g.det, g.v).statistic
    rc = [stats.spearmanr(s.length, s[metric]).statistic for _, s in cells.groupby(["run", "gene"]) if s.length.nunique() > 5]
    res["C5_rho_len_gene"] = round(float(r_len), 3)
    res["C5_rho_det_gene"] = round(float(r_det), 3)
    res["C5_rho_len_cell_median"] = round(float(np.nanmedian(rc)), 3)
    res["C5"] = abs(r_len) <= 0.3 and abs(r_det) <= 0.3 and abs(np.nanmedian(rc)) <= 0.2
    res["n_pass_C1_C3"] = int(res["C1"]) + int(res["C2"]) + int(res["C3"])
    res["n_pass_all"] = res["n_pass_C1_C3"] + int(res["C4"]) + int(res["C5"])
    return res, g


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    allres, allg = [], []
    for r in ["tls_high", "tls_low"]:
        c = pd.read_csv(SRC / f"results/{r}/cells.csv", dtype={"donor_id": str})
        c = c[~c.gene.isin(BLIND)].reset_index(drop=True)  # 盲検: 主要対象・参照を除外
        assert not set(c.gene) & BLIND
        alt = [x for x in c.columns if x.startswith("shift_alt")][0]
        c["alt"] = c[alt]
        c["pos"] = token_positions(c)
        det = json.loads((SRC / f"work/detection-{START[r]}.json").read_text())
        c = pd.concat([corrected(s) for _, s in c.groupby("run")], ignore_index=True)
        for m in METRICS:
            res, g = evaluate(c, det, m)
            res = {"run": r, "label": METRICS[m], **res}
            allres.append(res)
            g["run_name"] = r; g["metric"] = m
            allg.append(g)
        # 参考: 削除位置と shift の関係（null 細胞）
        n = c[c.role == "negative"]
        print(f"[{r}] null cells: spearman(shift_goal, pos)={stats.spearmanr(n.pos, n.shift_goal).statistic:.3f}, "
              f"(shift_goal, pos/len)={stats.spearmanr(n.pos / n.length, n.shift_goal).statistic:.3f}, "
              f"(shift_goal, length)={stats.spearmanr(n.length, n.shift_goal).statistic:.3f}; "
              f"control mean pos: " + ", ".join(f"{(POS | HK)[k]} {c[c.gene == k].pos.median():.0f}" for k in (POS | HK))
              + f"; null median pos {n.pos.median():.0f}")
    df = pd.DataFrame(allres)
    pd.concat(allg).to_csv(OUT / "gene_seed_metrics.csv", index=False)
    flat = df.copy()
    for col in ["C1_seeds", "C1_pct_mean", "C3_seeds", "C3_pct_mean"]:
        flat[col] = flat[col].map(lambda d: " ".join(f"{k}:{v}" for k, v in d.items()))
    flat.to_csv(OUT / "criteria.csv", index=False)
    (OUT / "criteria.json").write_text(json.dumps(allres, ensure_ascii=False, indent=1, default=str) + "\n")
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print(flat[["run", "label", "C1_pct_mean", "C2_wins", "C2_margin_nullSD", "C3_pct_mean", "C4i_splithalf_rho",
                "C4ii_seed_rho_fixed5", "C5_rho_len_gene", "C5_rho_det_gene", "C5_rho_len_cell_median",
                "null_frac_negative", "C1", "C2", "C3", "C4", "C5", "n_pass_all"]].to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())

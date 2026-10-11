"""EXP-20260927-003 事前登録の一部（摂動の実行前）: 線維芽細胞の遺伝子統計表を作る。

ISP（摂動）の結果は一切読まない。使うのは
  - トークン化済み入力（EXP-20260927-002 の tokenized.dataset を APFS clone したもの。guardrails 2026-09-26）
  - 生カウントの h5ad（EXP-20260927-002 の work/input.h5ad。読み取りのみ）
だけである。

出力 work/gene_stats_fib.csv の列
  ensembl, symbol（分かるもののみ）, token
  n_start        : start 状態（FIB_neg_TLShigh）の細胞のうち、この遺伝子が Geneformer 入力（rank 符号化後）に載る細胞数
  frac_start     : n_start / start 細胞数
  n_donor_start  : 同じく提供者数
  pos_med        : 入力内の位置（index_in_input_ids, <cls> を含む）の中央値（載っている細胞のみ）
  len_med        : その細胞の入力の長さの中央値
  relpos_med     : pos/length の中央値
  n_goal, frac_goal : goal 状態（FIB_org）での同じ検出細胞数・割合
  log2fc_start_vs_goal : 提供者ごとの疑似バルク CPM（生カウント）から log2((start+1)/(goal+1))
  n_alt, frac_alt   : alt 状態（FIB_neg_TLSlow）での検出細胞数・割合

実行: .venv/bin/python experiments/EXP-20260927-003/analysis/gene_stats_fib.py
"""
from __future__ import annotations

import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260927-003"
TOK = E / "work/tokenized.dataset"
H5 = ROOT / "experiments/EXP-20260927-002/work/input.h5ad"   # 読み取りのみ
START, GOAL, ALT = "FIB_neg_TLShigh", "FIB_org", "FIB_neg_TLSlow"
KEY = "fib_tls_state"


def token_dict() -> dict[str, int]:
    from geneformer import TOKEN_DICTIONARY_FILE
    with open(TOKEN_DICTIONARY_FILE, "rb") as f:
        return pickle.load(f)


def main() -> None:
    from datasets import load_from_disk

    ds = load_from_disk(str(TOK))
    states = np.array(ds[KEY])
    donors = np.array([str(x) for x in ds["donor_id"]])
    ids = list(ds["input_ids"])
    tokd = token_dict()
    ntok = max(max(tokd.values()) + 1, max(max(x) for x in ids) + 1)

    def mat(mask: np.ndarray):
        rows = np.flatnonzero(mask)
        sel = [ids[i] for i in rows]
        indptr = np.cumsum([0] + [len(x) for x in sel])
        col = np.concatenate([np.asarray(x) for x in sel])
        m = sp.csr_matrix((np.ones(len(col), bool), col, indptr), shape=(len(rows), ntok)).tocsc()
        return m, rows, sel

    m_s, rows_s, sel_s = mat(states == START)
    m_g, _, _ = mat(states == GOAL)
    m_a, _, _ = mat(states == ALT)
    n_s, n_g, n_a = m_s.shape[0], m_g.shape[0], m_a.shape[0]
    print(f"cells start={n_s} goal={n_g} alt={n_a}")

    # start 細胞の位置・長さ（token -> list）
    lens = np.array([len(x) for x in sel_s], dtype=float)
    pos_sum: dict[int, list[int]] = {}
    len_of: dict[int, list[float]] = {}
    for r, x in enumerate(sel_s):
        L = lens[r]
        for p, t in enumerate(x):
            pos_sum.setdefault(t, []).append(p)
            len_of.setdefault(t, []).append(L)

    dons_s = donors[rows_s]
    cnt_s = np.asarray(m_s.sum(0)).ravel()
    cnt_g = np.asarray(m_g.sum(0)).ravel()
    cnt_a = np.asarray(m_a.sum(0)).ravel()

    # 提供者数（start）
    don_codes, don_idx = np.unique(dons_s, return_inverse=True)
    D = sp.csr_matrix((np.ones(len(don_idx), bool), (don_idx, np.arange(len(don_idx)))),
                      shape=(len(don_codes), len(don_idx)))
    ndon = np.asarray(((D @ m_s) > 0).sum(0)).ravel()

    # ---- 生カウントから提供者ごとの疑似バルク CPM（start / goal）
    import anndata as ad
    A = ad.read_h5ad(H5)
    st = A.obs[KEY].astype(str).to_numpy()
    dn = A.obs["donor_id"].astype(str).to_numpy()
    var_ids = A.var_names.to_numpy().astype(str)
    X = A.X.tocsc() if sp.issparse(A.X) else sp.csc_matrix(A.X)

    def pseudobulk(mask):
        out = []
        for d in np.unique(dn[mask]):
            sub = X[np.flatnonzero(mask & (dn == d))]
            if sub.shape[0] < 5:
                continue
            v = np.asarray(sub.sum(0)).ravel()
            tot = v.sum()
            if tot <= 0:
                continue
            out.append(v / tot * 1e6)
        return np.mean(np.vstack(out), axis=0) if out else np.zeros(X.shape[1])

    cpm_s = pseudobulk(st == START)
    cpm_g = pseudobulk(st == GOAL)
    l2fc = np.log2((cpm_s + 1.0) / (cpm_g + 1.0))
    fc_by_id = dict(zip(var_ids, l2fc))
    cpm_s_by_id = dict(zip(var_ids, cpm_s))
    sym_by_id = dict(zip(var_ids, A.var["gene_name"].astype(str).to_numpy()))

    rows = []
    for ens, t in tokd.items():
        if not str(ens).startswith("ENS") or t >= len(cnt_s) or cnt_s[t] == 0:
            continue
        p = np.asarray(pos_sum.get(t, []), float)
        L = np.asarray(len_of.get(t, []), float)
        rows.append({
            "ensembl": ens, "symbol": sym_by_id.get(ens, ""), "token": int(t),
            "n_start": int(cnt_s[t]), "frac_start": float(cnt_s[t] / n_s), "n_donor_start": int(ndon[t]),
            "pos_med": float(np.median(p)), "len_med": float(np.median(L)),
            "relpos_med": float(np.median(p / L)),
            "n_goal": int(cnt_g[t]) if t < len(cnt_g) else 0,
            "frac_goal": float(cnt_g[t] / n_g) if t < len(cnt_g) else 0.0,
            "n_alt": int(cnt_a[t]) if t < len(cnt_a) else 0,
            "frac_alt": float(cnt_a[t] / n_a) if t < len(cnt_a) else 0.0,
            "log2fc_start_vs_goal": float(fc_by_id.get(ens, np.nan)),
            "cpm_start": float(cpm_s_by_id.get(ens, np.nan)),
        })
    d = pd.DataFrame(rows).sort_values("ensembl")
    d.to_csv(E / "work/gene_stats_fib.csv", index=False, float_format="%.6g")
    (E / "work/gene_stats_fib.meta.json").write_text(json.dumps({
        "tokenized": str(TOK.relative_to(ROOT)), "h5ad": str(H5.relative_to(ROOT)),
        "n_cells": {"start": int(n_s), "goal": int(n_g), "alt": int(n_a)},
        "n_genes_in_start_input": int(len(d)),
        "median_length": {"start": float(np.median(lens))},
        "uses_isp_results": False,
    }, ensure_ascii=False, indent=1) + "\n")
    print(d.shape, "→ work/gene_stats_fib.csv")


if __name__ == "__main__":
    main()

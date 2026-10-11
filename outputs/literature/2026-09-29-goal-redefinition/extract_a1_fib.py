#!/usr/bin/env python3
"""A1（単一核）の線維芽細胞について、候補遺伝子パネルの生カウントを抽出する（読み取り専用）。

出力: fib_A1_panel.csv.gz（細胞 × 遺伝子の生カウント + メタ情報）
      panel_meta.json（遺伝子の Ensembl ID、辞書の有無、欠けた遺伝子）
"""
import json, sys, pickle, warnings
from pathlib import Path
import h5py, numpy as np, pandas as pd, scipy.sparse as sp

warnings.filterwarnings("ignore")
ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
RAW = ROOT / "data/raw/kidney_tls/A1_KPMP_v1.5_snRNA_Lake2023.h5ad"
OUT = ROOT / "outputs/literature/2026-09-29-goal-redefinition"
TOKDICT = Path("/Users/petadimensionlab/workspace/research/Geneformer/geneformer_hf/"
               "geneformer/token_dictionary_gc104M.pkl")

FIB_LABELS = ["FIB", "M-FIB", "aFIB", "dFIB", "dM-FIB"]

# --- 目標の定義に使いうるマーカー候補（ケモカイン 3 遺伝子を含まない） ---
GOAL_MARKERS = [
    # リンパ組織様間質・オーガナイザー系
    "PDPN", "NGFR", "ICAM1", "VCAM1", "TNFSF13B", "TNFSF13", "CR2", "MADCAM1",
    "ENPP2", "IL7", "LTA", "TNFSF11", "GLYCAM1", "BST1", "ALDH1A2", "WNT5A",
    # 抗原提示
    "CD74", "HLA-DRA", "HLA-DRB1", "CIITA", "HLA-DPA1",
    # 補体・分泌
    "C3", "C7", "CFD", "SERPING1", "C1S", "C1R", "PTX3", "CFH",
    # 線維芽細胞のサブタイプ／汎マーカー
    "PI16", "CD34", "CLU", "GREM1", "SFRP2", "SFRP4", "FAP", "THY1", "TNC",
    "POSTN", "COL1A1", "COL1A2", "COL3A1", "LUM", "DCN", "PDGFRA", "PDGFRB",
    "FN1", "ITGA8", "MEG3", "NR4A1", "PTGDS", "APOE", "S100A4", "ELN", "COL15A1",
    # 筋線維芽細胞・血管周囲（対照の軸）
    "ACTA2", "TAGLN", "MYH11", "RGS5", "NOTCH3", "DES",
]

# --- 摂動の対象にしうる遺伝子（分泌プログラム。盲検 5 遺伝子は含めない） ---
PERT_CANDIDATES = [
    "CXCL13", "CCL19", "CCL21", "LTB", "TNFSF14", "CXCL8", "IL6", "CCL2",
    "CXCL9", "CXCL10", "CXCL11", "CCL5", "CCL20", "IL33", "IL1B", "CSF1",
    "IL15", "IL34", "CXCL14", "CXCL1", "CXCL2", "CXCL3", "CXCL16", "CXCL12",
    "CCL8", "CCL13", "CCL7", "CCL11", "CCL28", "LIF", "IL11", "IL32",
    "TNFSF10", "TNFSF15", "IL1A", "IL18", "CSF2", "CSF3", "CXCL5", "CXCL6",
    "XCL1", "IL16", "IL34", "VEGFA", "PGF", "ANGPT1",
]

# --- 細胞老化に関係する遺伝子（盲検の TP53/CDKN1A/LTBR/SPP1/LCN2 は除く） ---
SENESCENCE = [
    "CDKN2A", "CDKN2B", "SERPINE1", "IGFBP3", "IGFBP7", "GDF15", "B2M",
    "MMP1", "MMP2", "MMP3", "TIMP1", "TIMP2", "CCND1", "LMNB1", "E2F1",
    "RB1", "GLB1", "TNFRSF10C", "NFKB1", "NFKB2", "RELB", "TRAF1", "BIRC3",
]

# 盲検中（摂動も解析もしない）。抽出もしない。
BLINDED = ["TP53", "CDKN1A", "LTBR", "SPP1", "LCN2"]

PANEL = []
for g in GOAL_MARKERS + PERT_CANDIDATES + SENESCENCE:
    if g not in PANEL and g not in BLINDED:
        PANEL.append(g)


def rd(h, c):
    g = h["obs"][c]
    if isinstance(g, h5py.Group):
        cats = np.array(g["categories"].asstr()[...])
        codes = g["codes"][...]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "nan")
    v = g[...]
    return np.array(g.asstr()[...]) if v.dtype.kind in "OS" else v


def main():
    tok = pickle.loads(TOKDICT.read_bytes())
    with h5py.File(RAW, "r") as h:
        ct = rd(h, "subclass.l2")
        don = rd(h, "donor_id")
        dis = rd(h, "disease")
        n = len(ct)
        fib = np.isin(ct, FIB_LABELS)
        v = h["raw/var"]
        ens = np.array(v["_index"].asstr()[...])
        fn = v["feature_name"]
        sym = (np.array(fn["categories"].asstr()[...])[fn["codes"][...]]
               if isinstance(fn, h5py.Group) else np.array(fn.asstr()[...]))
        s2i = {}
        for i, s in enumerate(sym):
            s2i.setdefault(s, i)
        present = [g for g in PANEL if g in s2i]
        missing = [g for g in PANEL if g not in s2i]
        cols = np.array([s2i[g] for g in present])
        meta = {
            "n_cells_total": int(n),
            "n_fibroblast": int(fib.sum()),
            "fib_labels": FIB_LABELS,
            "panel_size": len(PANEL),
            "genes_missing_in_var": missing,
            "blinded_excluded": BLINDED,
            "gene_ensembl": {g: str(ens[s2i[g]]) for g in present},
            "gene_in_geneformer_dict": {g: bool(str(ens[s2i[g]]) in tok) for g in present},
            "gene_token": {g: int(tok[str(ens[s2i[g]])]) for g in present
                           if str(ens[s2i[g]]) in tok},
            "token_dict": str(TOKDICT),
            "token_dict_size": len(tok),
        }
        Xg = h["raw/X"]
        indptr = Xg["indptr"][...]
        a, b = indptr[0], indptr[min(200, n)]
        chk = Xg["data"][a:b][:100000]
        meta["counts_integer_like"] = bool(np.allclose(chk, np.round(chk)))

        blocks, tot, ng, rows = [], [], [], []
        step = 20000
        for s in range(0, n, step):
            e = min(n, s + step)
            if not fib[s:e].any():
                continue
            a, b = indptr[s], indptr[e]
            m = sp.csr_matrix((Xg["data"][a:b], Xg["indices"][a:b],
                               indptr[s:e + 1] - a), shape=(e - s, len(ens)))
            loc = np.where(fib[s:e])[0]
            m = m[loc]
            tot.append(np.asarray(m.sum(1)).ravel())
            ng.append(np.diff(m.indptr))
            blocks.append(m[:, cols].toarray().astype(np.int32))
            rows.append(loc + s)
            print(f"  {e}/{n}", file=sys.stderr, end="\r")
        rows = np.concatenate(rows)
        D = pd.DataFrame(np.vstack(blocks), columns=present)
        D.insert(0, "row", rows)
        D.insert(1, "donor", don[rows])
        D.insert(2, "disease", dis[rows])
        D.insert(3, "subclass", ct[rows])
        D["n_counts"] = np.concatenate(tot)
        D["n_genes"] = np.concatenate(ng)

    D.to_csv(OUT / "fib_A1_panel.csv.gz", index=False, compression="gzip")
    (OUT / "panel_meta.json").write_text(
        json.dumps(meta, indent=1, ensure_ascii=False))
    print(f"\n線維芽細胞 {len(D)} 細胞 / 遺伝子 {len(present)}（欠け {len(missing)}）", file=sys.stderr)


if __name__ == "__main__":
    main()

"""EXP-20260926-005: 陽性対照候補の (i) 実行可能性（選抜できる細胞数）と (ii) マーカー度（start vs goal の発現差）。

事前登録の一部。ISP の結果は一切使わない（トークン化済み入力の発現情報だけ）。
  det_start : start（aPT_CKD_TLShigh, 7,533 細胞）で Geneformer 入力に含まれる細胞の割合
  det_goal  : goal（PT_healthy のうち held-out 9 donor）で同じ割合
  ddet      : det_start − det_goal（＝ aPT 側マーカーなら正, 健常側マーカーなら負）
  ddet_pct  : |ddet| の背景分布（どちらかの状態で 300 細胞以上検出される遺伝子）での分位
  rank_pct  : 発現細胞での rank 位置の中央値（0 = 最上位）。overexpress の効きしろの目安
実行: .venv/bin/python experiments/EXP-20260926-005/analysis/gene_annotation.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts/geneformer"))
import run_isp as ri  # noqa: E402

E = ROOT / "experiments/EXP-20260926-005"
CFG = yaml.safe_load(open(E / "config-delete.yaml"))
TOK = ROOT / CFG["data"]["tokenized_dataset"]


def main() -> None:
    import pickle
    from datasets import load_from_disk
    from geneformer import TOKEN_DICTIONARY_FILE
    gdict = ri.gene_token_dict()
    with open(Path(TOKEN_DICTIONARY_FILE).with_name("gene_name_id_dict_gc104M.pkl"), "rb") as f:
        name2id = pickle.load(f)
    ds = load_from_disk(str(TOK))
    state = ds[CFG["data"]["state_key"]]
    donor = [str(x) for x in ds["donor_id"]]
    goal_donors = set(CFG["state_embs"]["include_donors"]["PT_healthy"])
    st, go = CFG["states"]["start_state"], CFG["states"]["goal_state"]
    i_start = [i for i, s in enumerate(state) if s == st]
    i_goal = [i for i, s in enumerate(state) if s == go and donor[i] in goal_donors]
    ids = ds["input_ids"]

    ntok = max(gdict.values()) + 1
    cnt_s, cnt_g = np.zeros(ntok, int), np.zeros(ntok, int)
    rank_sum, rank_n = np.zeros(ntok), np.zeros(ntok, int)
    for i in i_start:
        x = np.asarray(ids[i])
        cnt_s[x] += 1
        rank_sum[x] += np.arange(len(x)) / len(x)
        rank_n[x] += 1
    for i in i_goal:
        cnt_g[np.asarray(ids[i])] += 1

    fs, fg = cnt_s / len(i_start), cnt_g / len(i_goal)
    ddet = fs - fg
    bg = (cnt_s >= 300) | (cnt_g >= 300)
    bg_abs = np.sort(np.abs(ddet[bg]))

    # 実行可能性: 遺伝子ごとに seed 3/4/5 で実際に選抜できる細胞数
    m, rows, donors = ri.start_matrix(TOK, CFG)
    sel = CFG["cell_selection"]
    import random as _r

    genes: dict[str, str] = {}
    for grp in ["positive_delete", "positive_overexpress", "prev_positive", "housekeeping", "negctrl_panel"]:
        for n in GROUPS[grp]:
            genes[n] = grp

    out = []
    for name, grp in genes.items():
        gid = name2id.get(name)
        tid = gdict.get(gid)
        idx = m[:, tid].nonzero()[0].tolist() if tid is not None else []
        ns = []
        for seed in CFG["robustness"]["seeds"]:
            keep = ri.balanced_subset(idx, donors, sel["max_cells_per_gene"], sel["max_cells_per_donor"],
                                      _r.Random(f"{seed}-{gid}"))
            ns.append(len(keep))
        nd = len({donors[i] for i in idx})
        out.append({
            "gene": name, "ensembl": gid, "group": grp,
            "det_start": int(cnt_s[tid]), "det_goal": int(cnt_g[tid]),
            "f_start": round(float(fs[tid]), 4), "f_goal": round(float(fg[tid]), 4),
            "ddet": round(float(ddet[tid]), 4),
            "ddet_pct": round(float(np.searchsorted(bg_abs, abs(ddet[tid])) / len(bg_abs)), 3),
            "rank_pct_start": round(float(rank_sum[tid] / max(rank_n[tid], 1)), 3),
            "n_donor_start": nd, "n_cells_seeds": ns,
        })
    out.sort(key=lambda r: (r["group"], -abs(r["ddet"])))
    (E / "work" / "gene_annotation.json").write_text(json.dumps(
        {"n_start": len(i_start), "n_goal": len(i_goal), "n_bg_genes": int(bg.sum()), "genes": out},
        ensure_ascii=False, indent=1) + "\n")
    w = f"{'gene':10s} {'group':21s} {'detS':>6s} {'detG':>6s} {'fS':>6s} {'fG':>6s} {'ddet':>7s} {'pct':>5s} {'rank':>5s} {'nd':>3s} cells"
    print(f"start {len(i_start)} / goal {len(i_goal)} 細胞、背景 {int(bg.sum())} 遺伝子")
    print(w)
    for r in out:
        print(f"{r['gene']:10s} {r['group']:21s} {r['det_start']:6d} {r['det_goal']:6d} {r['f_start']:6.3f} "
              f"{r['f_goal']:6.3f} {r['ddet']:7.3f} {r['ddet_pct']:5.2f} {r['rank_pct_start']:5.2f} "
              f"{r['n_donor_start']:3d} {r['n_cells_seeds']}")


GROUPS = {
    "positive_delete": ["HIF1A", "NOTCH1", "EGFR", "YAP1", "KEAP1"],
    "positive_overexpress": ["HNF4A", "PPARGC1A", "PPARA", "ESRRA", "CPT1A", "NFE2L2"],
    "prev_positive": ["VCAM1", "HAVCR1"],
    "housekeeping": ["ACTB", "GAPDH", "PPIA"],
    "negctrl_panel": ["CERS4", "TATDN1", "IFT25", "VAPA", "PHPT1", "ANAPC16", "EEF1D"],
}

if __name__ == "__main__":
    main()

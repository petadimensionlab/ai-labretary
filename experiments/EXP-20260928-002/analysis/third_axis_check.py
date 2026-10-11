"""EXP-20260928-002 — 3 本目の軸（対照の状態 FIB_neg_TLSlow）でのセット検定。**事後の検証用**。

作成: Labretary（validator） / 2026-09-29。`validation.md` §1.3 / §3.1 / §3.4 の数値の出所。

**このファイルは事前登録（work/SHA256SUMS-preregistered.txt の 38 ファイル）に含まれない。**
事前登録は結果を見る前に固定したものであり、本ファイルは結果を見た後に書いた検証用のコードである。
**判定・閾値・事前登録ファイルは一切変更しない。** 本ファイルの出力は判定に使わない。

なぜ必要か
----------
plan.md §6.1 は「軸を入れ替えて 3 本つくる（目標 / 筋線維芽細胞 / 置き換え）」と書いているが、
evaluate_002.py のセット検定は目標の軸と筋線維芽細胞の軸の 2 本にしか回っていない。
3 本目（対照の状態＝三次リンパ組織が少ない提供者の腎線維芽細胞）を入れると、
経路セットは目標の状態より対照の状態のほうへわずかに大きく動いていることが分かる。
これは交差（G-2）の合格の読み方を変える所見なので、検算できる形で残す。

**二次回帰の補正は省略できない。** 補正を省いた素の平均では、
どの重心へも一様に近づく成分（元の位置・入力の長さで決まる）が残り、
3 本の軸の大小関係が逆転する。本ファイルはその比較（`補正なしとの比較`）も出力する。

検算のしかた
------------
本ファイルは公式の評価コード analysis/evaluate_002.py の
metric_columns() / design() / gene_stats() と**同じ式**を使っている。
同一性は出力の `m_goal_の再現` で確かめられる（公式の
results/eval/gene_seed_002.csv の D_raw との差を 3 軸すべてについて出す）。

実行: .venv/bin/python experiments/EXP-20260928-002/analysis/third_axis_check.py
出力: experiments/EXP-20260928-002/results/eval/third_axis_002.json
所要: 約 1 分（埋め込み 1,136 遺伝子 × 3 seed を読み直す）
"""
from __future__ import annotations

import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20260928-002"
OUT = E / "results/eval/third_axis_002.json"

# 摂動していない（してはいけない）遺伝子。1 行でも残っていれば停止する。evaluate_002.py と同一
BLOCKED = {"ENSG00000111321", "ENSG00000104856", "ENSG00000077150",
           "ENSG00000124762", "ENSG00000141510", "ENSG00000118785", "ENSG00000148346",
           "ENSG00000172724", "ENSG00000137077", "ENSG00000156234"}

AXES = {"目標の軸": "m_goal", "筋線維芽細胞の軸": "m_myofaxis", "対照の状態の軸": "m_placebo"}


def cos(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """evaluate_002.py と同一。"""
    return (a @ b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b) + 1e-12)


def design(lp, ll) -> np.ndarray:
    """元の位置 log2(1+index) と入力の長さ log2(length) の二次の計画行列。evaluate_002.py と同一。"""
    lp, ll = np.asarray(lp, float), np.asarray(ll, float)
    return np.column_stack([np.ones_like(lp), lp, ll, lp ** 2, ll ** 2, lp * ll])


def load(cfg: dict) -> pd.DataFrame:
    """細胞ごとの 4 重心へのコサイン差と移動量を、埋め込みから作り直す。evaluate_002.py の load() と同一。"""
    res = E / cfg["output"]["results_dir"]
    c = pd.read_csv(res / "cells.csv", dtype={"donor_id": str})
    n0 = len(c)
    c = c[~c.gene.isin(BLOCKED)].reset_index(drop=True)
    if n0 != len(c):
        raise SystemExit(f"回さないはずの遺伝子が {n0 - len(c)} 行あった")
    c["seed"] = c.run.astype(str).str.extract(r"run-(\d+)", expand=False).astype(int)

    from datasets import load_from_disk
    from geneformer import TOKEN_DICTIONARY_FILE
    tokd = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
    ds = load_from_disk(str(E / "work/tokenized.dataset"))
    uids = sorted(c.cell_uid.unique())
    ids = {u: np.asarray(x) for u, x in zip(uids, ds.select(uids)["input_ids"])}
    pos = np.empty(len(c), dtype=np.int64)
    for gene, part in c.groupby("gene", sort=False):
        t = tokd[gene]
        pos[part.index.to_numpy()] = [int(np.flatnonzero(ids[u] == t)[0]) for u in part.cell_uid]
    c["pos"] = pos
    c["lp"] = np.log2(1 + c.pos.astype(float))
    c["ll"] = np.log2(c.length.astype(float))

    names = ["goal", "start", "alt", "myof"]
    vals = {n: np.full(len(c), np.nan) for n in names}
    dnorm = np.full(len(c), np.nan)
    s_ = cfg["states"]
    max_cells_csv_diff = 0.0
    for (run, gene), part in c.groupby(["run", "gene"], sort=False):
        seed = int(str(run).split("-")[1])
        sz = np.load(res / run / "state_emb_vectors.npz")
        xz = np.load(E / f"work/extra_centroids-seed{seed}.npz")
        cen = {"goal": sz[s_["goal_state"]], "start": sz[s_["start_state"]],
               "alt": sz[s_["alt_states"][0]], "myof": xz["MYOF_neg"]}
        gd = res / run / "genes" / gene / "emb"
        o = np.load(gd / "orig.npy").astype(np.float64)
        p = np.load(gd / "pert.npy").astype(np.float64)
        assert len(o) == len(part), f"{gd}: 埋め込み {len(o)} ≠ cells.csv {len(part)}"
        idx = part.index.to_numpy()
        for n in names:
            vals[n][idx] = cos(p, cen[n]) - cos(o, cen[n])
        dnorm[idx] = np.linalg.norm(p - o, axis=1)
        # cells.csv に記録済みの shift_* と、埋め込みからの再計算が一致することの確認
        max_cells_csv_diff = max(
            max_cells_csv_diff,
            float(np.abs(vals["goal"][idx] - part.shift_goal.to_numpy()).max()),
            float(np.abs(vals["start"][idx] - part.shift_start.to_numpy()).max()),
            float(np.abs(vals["alt"][idx] - part["shift_alt:FIB_neg_TLSlow"].to_numpy()).max()))
    for n in names:
        c[f"s_{n}"] = vals[n]
    c["dnorm"] = dnorm
    c.attrs["max_cells_csv_diff"] = max_cells_csv_diff

    # 軸 3 本。evaluate_002.py の metric_columns() と同一の式（同じ関数形・同じ尺度）
    g, st, al, my, dn = c.s_goal, c.s_start, c.s_alt, c.s_myof, c.dnorm + 1e-30
    form = lambda t, o: (t - sum(o) / len(o)) / dn
    c["m_goal"] = form(g, [st, al, my])
    c["m_myofaxis"] = form(my, [st, al, g])
    c["m_placebo"] = form(al, [st, g, my])          # ← 3 本目。対照の状態に特異的な成分
    # 補正を省いた素の平均（比較用。判定には使わない）
    c["raw_goal"] = g
    c["raw_myofaxis"] = my
    c["raw_placebo"] = al
    return c


def gene_D_raw(cs: pd.DataFrame, null_genes: set[str], col: str) -> pd.Series:
    """遺伝子ごとに、対照遺伝子（自分自身は除く）で当てた二次回帰からの残差の平均 D_raw。

    evaluate_002.py の gene_stats() と**同じ式**（正規方程式版）。
    自分が対照遺伝子のときは、自分の寄与 X'X・X'y を全体から引いてから解く。"""
    nulls = {h: hd for h, hd in cs[cs.gene.isin(null_genes)].groupby("gene")}
    G = {h: (X := design(hd.lp, hd.ll)).T @ X for h, hd in nulls.items()}
    B = {h: design(hd.lp, hd.ll).T @ hd[col].to_numpy() for h, hd in nulls.items()}
    G_all, B_all = sum(G.values()), sum(B.values())
    out = {}
    for gene, tg in cs.groupby("gene"):
        Gm, Bm = (G_all - G[gene], B_all - B[gene]) if gene in G else (G_all, B_all)
        beta = np.linalg.solve(Gm, Bm)
        out[gene] = float(tg[col].to_numpy().mean() - design(tg.lp, tg.ll).mean(axis=0) @ beta)
    return pd.Series(out, name=col)


def main() -> None:
    cfg = yaml.safe_load(open(E / "config.yaml"))
    cr = cfg["analysis"]["criteria"]
    seeds = sorted(cfg["robustness"]["seeds"])

    gj = json.loads((E / "work/genes_002.json").read_text())["genes"]
    ens = {x["symbol"]: x["ensembl"] for x in gj}
    sym = {v: k for k, v in ens.items()}
    roles = {x["symbol"]: x["role"] for x in gj}
    npool = json.loads((E / "work/null_pool_002.json").read_text())["versions"]
    NULL = set(ens[x["symbol"]] for x in npool["A_980"])          # 主判定の対照 980 個
    subs = json.loads((E / "work/subsets_002.json").read_text())["subsets"]
    sets_def = {k: [ens[x["symbol"]] for x in v] for k, v in subs.items()}
    sets_def |= {"筋線維芽細胞マーカー": [ens[s] for s, r in roles.items() if r == "marker_myof"],
                 "恒常発現遺伝子": [ens[s] for s, r in roles.items() if r == "housekeeping"],
                 "中立対照パネル": [ens[s] for s, r in roles.items() if r == "negctrl_panel"]}
    JP = {"pathway_all": "経路セット", "nonmarker": "非マーカー部分セット",
          "drop_goalmarker": "目標マーカー 8 個を抜いた版", "pathway_specific": "経路特異的な版",
          "general_mechanism": "一般機構", "tier_S1": "層: 経路の構成因子",
          "tier_S2": "層: 一次文献の転写標的", "tier_S3": "層: キュレーション済みセット"}
    order = ["pathway_all", "nonmarker", "drop_goalmarker", "pathway_specific", "general_mechanism",
             "tier_S1", "tier_S2", "tier_S3", "筋線維芽細胞マーカー", "恒常発現遺伝子", "中立対照パネル"]

    c = load(cfg)
    assert sorted(c.seed.unique()) == seeds, "seed が config と違う"

    # --- 遺伝子ごとの値と、対照 980 個での標準化（z）。3 軸すべてについて
    Zt: dict[tuple[str, int], pd.Series] = {}
    Rt: dict[tuple[str, int], pd.Series] = {}      # 補正なし（比較用）
    Dt: dict[tuple[str, int], pd.Series] = {}      # 補正後の生の値（公式との照合用）
    for s in seeds:
        cs = c[c.seed == s]
        for jp, col in AXES.items():
            D = gene_D_raw(cs, NULL, col)
            Dt[(col, s)] = D
            nu = D[D.index.isin(NULL)]
            Zt[(jp, s)] = (D - nu.mean()) / nu.std(ddof=1)
            R = cs.groupby("gene")[col.replace("m_", "raw_")].mean()
            nr = R[R.index.isin(NULL)]
            Rt[(jp, s)] = (R - nr.mean()) / nr.std(ddof=1)

    # --- (3) 公式の指標と同じ補正を使っていることの確認
    off = pd.read_csv(E / "results/eval/gene_seed_002.csv")
    repro = {}
    for col in AXES.values():
        diffs = []
        for s in seeds:
            o = off[(off.col == col) & (off.seed == s)].set_index("gene").D_raw
            j = Dt[(col, s)].to_frame("mine").join(o)
            assert len(j) == 1136 and j.D_raw.notna().all(), f"{col} seed {s}: 遺伝子数が合わない"
            diffs.append(float(np.abs(j.mine - j.D_raw).max()))
        repro[col] = {"公式 gene_seed_002.csv の D_raw との最大絶対差（seed 別）": diffs,
                      "最大": max(diffs), "n_遺伝子": 1136}
    repro["cells.csv の shift_* と埋め込みからの再計算の最大絶対差"] = c.attrs["max_cells_csv_diff"]
    repro["注"] = ("公式と同じ正規方程式・同じ自己除外を使っているので、差は浮動小数点の丸めの水準に収まる。"
                   "validation.md §2.1 が報告した 2.1e-6 は、自己除外を省いた簡易版での照合値であり、"
                   "本ファイルはそれより厳密な照合になっている。")

    # --- 無作為な対照セットの分布（参考値。判定には使わない）
    rng = np.random.default_rng(cr["null_draw_seed"])
    znull = {jp: np.column_stack([Zt[(jp, s)][sorted(NULL)].to_numpy() for s in seeds]) for jp in AXES}
    cache: dict[tuple[str, int], np.ndarray] = {}

    def null_dist(jp: str, k: int) -> np.ndarray:
        if (jp, k) not in cache:
            zn = znull[jp]
            cache[(jp, k)] = zn[rng.integers(0, zn.shape[0], size=(cr["n_permutation"], k))].mean(axis=1)
        return cache[(jp, k)]

    # --- (2) セットごと × 3 seed × 3 軸
    sets_out = {}
    for name in order:
        genes = sets_def[name]
        row = {"K": len(genes), "seeds": seeds}
        for jp in AXES:
            z = np.column_stack([Zt[(jp, s)][genes].to_numpy() for s in seeds])
            obs = z.mean(axis=0)
            nd = null_dist(jp, len(genes))
            row[jp] = {
                "セット平均_seed別": [round(float(x), 4) for x in obs],
                "セット平均_3seed平均": round(float(obs.mean()), 4),
                "無作為セットの99パーセンタイル_参考": [round(float(np.quantile(nd[:, i], 0.99)), 4)
                                                       for i in range(len(seeds))],
                "無作為セットの中での分位_参考": [round(float((nd[:, i] < obs[i]).mean()), 4)
                                                  for i in range(len(seeds))],
                "遺伝子間SD_3seed平均": round(float(z.mean(axis=1).std(ddof=1)), 4),
                "95%信頼区間_遺伝子サンプリング_3seed平均": [
                    round(float(z.mean() - 1.96 * z.mean(axis=1).std(ddof=1) / np.sqrt(len(genes))), 4),
                    round(float(z.mean() + 1.96 * z.mean(axis=1).std(ddof=1) / np.sqrt(len(genes))), 4)],
            }
        g_ = np.column_stack([Zt[("目標の軸", s)][genes].to_numpy() for s in seeds]).mean(axis=0)
        m_ = np.column_stack([Zt[("筋線維芽細胞の軸", s)][genes].to_numpy() for s in seeds]).mean(axis=0)
        p_ = np.column_stack([Zt[("対照の状態の軸", s)][genes].to_numpy() for s in seeds]).mean(axis=0)
        row["目標の軸 − 対照の状態の軸"] = {
            "seed別": [round(float(x), 4) for x in g_ - p_],
            "3seed平均": round(float((g_ - p_).mean()), 4),
            "3seedとも正か": bool(all(g_ - p_ > 0))}
        row["目標の軸 − 筋線維芽細胞の軸（交差 G-2 と同じ量）"] = {
            "seed別": [round(float(x), 4) for x in g_ - m_],
            "3seed平均": round(float((g_ - m_).mean()), 4),
            "3seedとも正か": bool(all(g_ - m_ > 0))}
        # 補正を省いた場合（比較用）
        row["補正なしとの比較"] = {
            jp: [round(float(np.column_stack([Rt[(jp, s)][genes].to_numpy()
                                              for s in seeds]).mean(axis=0)[i]), 4)
                 for i in range(len(seeds))] for jp in AXES}
        sets_out[JP.get(name, name)] = row

    # --- (4) 偽経路セットの中心が 0 からずれていること
    import ast
    sn = pd.read_csv(E / "results/eval/set_null_002.csv")
    sn["zs"] = sn.z.apply(lambda x: np.array(ast.literal_eval(x)))
    sham_out = {}
    for label, mask, k in [("経路セット 120 に合わせた偽セット", ~sn.set.str.contains("nonmarker"), 120),
                           ("非マーカー 76 に合わせた偽セット", sn.set.str.contains("nonmarker"), 76)]:
        a = np.vstack(sn[mask].zs.to_numpy())
        nd = null_dist("目標の軸", k)
        sham_out[label] = {
            "組数": int(mask.sum()),
            "偽セットの平均_seed別": [round(float(x), 4) for x in a.mean(axis=0)],
            "偽セットの平均_全体": round(float(a.mean()), 4),
            "無制限の無作為セットの平均": round(float(nd.mean()), 4),
            "中心のずれ（偽セット − 無作為）": round(float(a.mean() - nd.mean()), 4),
            "偽セットの散らばり（組間SD）": round(float(a.mean(axis=1).std(ddof=1)), 4),
            "無作為セットの散らばり": round(float(nd.mean(axis=1).std()), 4),
            "散らばりの比": round(float(a.mean(axis=1).std(ddof=1) / nd.mean(axis=1).std()), 4)}
    sham_out["読み方"] = ("元の位置・入力の長さ・細胞数・目標側への偏りをそろえる手続きそのものが、"
                          "偽セットの分布の中心と散らばりを動かしている。"
                          "中心が負にずれる非マーカー版では『観測値が 50 組中 49 組より上』という順位が過大になる。"
                          "本実験では主判定が不合格なので偽経路セットの検査（G-6）は結論に効かないが、"
                          "合格したときの裏づけに使う設計だったので、次の実験では中心と散らばりを検査すること。")

    M = {
        "実験": "EXP-20260928-002",
        "作成": "Labretary（validator） / 2026-09-29",
        "位置づけ": ("**事後の検証**。事前登録（work/SHA256SUMS-preregistered.txt の 38 ファイル）に含まれない。"
                     "判定・閾値・事前登録ファイルは変更していない。本ファイルの数値は判定に使わない。"
                     "結論の確定は PI（結論確定の関門 G3）が行う。"),
        "なぜ必要か": ("plan.md §6.1 は軸を 3 本つくると書いているが、evaluate_002.py のセット検定は "
                       "目標の軸と筋線維芽細胞の軸の 2 本にしか回っていない。3 本目（対照の状態）を入れると、"
                       "交差（G-2）の合格の読み方が変わる。"),
        "軸の定義": {
            "目標の軸": "(s_goal − (s_start + s_alt + s_myof)/3) / ‖Δ‖",
            "筋線維芽細胞の軸": "(s_myof − (s_start + s_alt + s_goal)/3) / ‖Δ‖",
            "対照の状態の軸": "(s_alt − (s_start + s_goal + s_myof)/3) / ‖Δ‖  ← 3 本目",
            "注": ("s_X は摂動後と摂動前の『重心 X とのコサイン』の差、‖Δ‖ は埋め込みの移動量。"
                   "3 本は同じ関数形・同じ尺度なので互いに比べてよい"
                   "（evaluate_002.py の compare_same_scale が許す組み合わせ）。"),
            "重心": {"目標": "FIB_org（深さをそろえた）", "出発点": "FIB_neg_TLShigh",
                     "対照の状態": "FIB_neg_TLSlow", "筋線維芽細胞": "MYOF_neg（手元の Mac で作成）"}},
        "補正": ("遺伝子ごとの値 D_raw は、元の位置 log2(1+index) と入力の長さ log2(length) の二次回帰で "
                 "対照遺伝子 980 個から予測した値からの残差の平均。z はその値を同じ seed の対照 980 個で標準化したもの。"
                 "**補正を省くと、どの重心へも一様に近づく成分が残って 3 軸の大小関係が逆転する**"
                 "（各セットの『補正なしとの比較』を参照）。"),
        "m_goal_の再現": repro,
        "セットごとの結果": sets_out,
        "偽経路セットの中心": sham_out,
        "主な所見": [
            "経路セット 120 は、目標の軸（3 seed 平均 %.3f）より対照の状態の軸（%.3f）のほうへ大きく動いている。"
            % (sets_out["経路セット"]["目標の軸"]["セット平均_3seed平均"],
               sets_out["経路セット"]["対照の状態の軸"]["セット平均_3seed平均"]),
            "非マーカー部分セット 76 は 3 seed とも「目標の軸 − 対照の状態の軸」が負である。",
            "したがって交差（G-2）が合格したのは「筋線維芽細胞から離れた」ためで、"
            "「目標の状態へ寄った」ためではない。",
            "層のうち一次文献の転写標的 13 個だけが「目標へ寄り、対照の状態から離れる」形を示すが、"
            "K が小さすぎて検定にならない（信頼区間は 0 を含む）。記述値にとどめる。"],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(M, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"書き出し: {OUT}")
    print(f"m_goal の再現（公式との最大差）: {repro['m_goal']['最大']:.3e}")
    for n in ("経路セット", "非マーカー部分セット"):
        r = sets_out[n]
        print(f"{n}（K={r['K']}）: 目標 {r['目標の軸']['セット平均_seed別']} / "
              f"筋線維芽 {r['筋線維芽細胞の軸']['セット平均_seed別']} / "
              f"対照の状態 {r['対照の状態の軸']['セット平均_seed別']} / "
              f"目標−対照の状態 {r['目標の軸 − 対照の状態の軸']['seed別']}")


if __name__ == "__main__":
    main()

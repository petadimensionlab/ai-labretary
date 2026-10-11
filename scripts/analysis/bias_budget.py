#!/usr/bin/env python3
"""交絡による歪みの直接見積もり（歪みの予算 / bias budget）。

**何をする道具か**

セット（遺伝子の集まり）の平均 z を判定に使うとき、そのセットが交絡の候補
（移動量・元の順位・入力の長さ・相対位置・一様成分）について対照から偏っていると、
**中身がまったく中立でもセット平均 z がずれる**。この道具は、そのずれの大きさを
**標準偏差の単位で直接見積もる**。

    交絡の候補 v ごとに
        b_v = （対照だけで測った「z と v の対照内分位」の傾き） × （セットの v の平均分位 − 0.5）
    歪みの合計 B = Σ_v b_v（同時回帰の偏回帰係数を使う）

合格条件は「**守りたい判断**（セットの絶対値条件 0.41 標準偏差）を歪みが覆さないこと」
から決める。既定では **|B| の 95% 信頼区間の上側 ≤ 0.10 標準偏差**（判断の線の 24.4%）。

**この道具が見ないもの**

- セットの主指標の値（目標の向きへの射影）を**一度も使わない**。使うのは
  (1) 対照遺伝子だけで測った傾きと、(2) セットの遺伝子の交絡の候補の値、の 2 つだけ。
- したがって**主判定の値を見る前に合否が決まる**。
- ただし交絡の候補のうち「移動量」と「一様成分」は摂動の出力から作られる量である
  （入力だけからは作れない）。**主指標の向きには依存しない**が、「入力側の性質だけ」
  ではない。この区別は報告で省略しない。

**周辺回帰と同時回帰**

交絡の候補どうしは強く重なっている（対照 500 個で移動量と元の順位の順位相関 0.878、
移動量と相対位置 0.787）。候補を 1 つずつ当てた**周辺回帰**の b_v を単純に足すと、
**同じ 1 つの歪みを 5 回数える**ことになり、歪みが無い指標でも合計が大きく出る。
既定の判定は**同時回帰**（5 つを一度に入れた重回帰）の偏回帰係数を使う。
周辺回帰版の Σ|b_v| も必ず併記する（重複して数えた上限として読む）。

使い方（呼び出し側）:

    from bias_budget import bias_budget, bias_budget_seeds
    r = bias_budget(gene_table, set_genes, control_genes)

`gene_table` は**遺伝子ごと** 1 行の表で、少なくとも次の列を持つこと:
    gene, z, dnorm_med, lp_med, ll_med, relpos_med, uniform
（`z` は対照で標準化した残差。`uniform` は一様成分の遺伝子ごとの残差。）

自己検査:  python3 scripts/analysis/bias_budget.py --self-test
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- 既定値
#: 呼び出し側が用意しなければならない列（2026-10-01 に列の写し漏れで評価が止まった教訓）
REQUIRED_COLUMNS: tuple[str, ...] = ("gene", "z")

#: 交絡の候補（既定）。表の列名で指定する
CANONICAL_COVARIATES: tuple[str, ...] = (
    "dnorm_med",    # 移動量 ‖Δ‖ の遺伝子ごとの中央値
    "lp_med",       # 元の順位 log2(1+pos)
    "ll_med",       # 入力の長さ log2(length)
    "relpos_med",   # 相対位置 pos / length
    "uniform",      # 一様成分（遺伝子ごとの残差）
)

#: 守りたい判断（セットの絶対値条件、標準偏差の単位）
DECISION_LINE = 0.41

#: 歪みの予算（標準偏差の単位）。DECISION_LINE の 24.4%
BUDGET = 0.10

#: 信頼区間の片側係数（正規近似の 95%）
CI_Z = 1.959964


# --------------------------------------------------------------------------- 分位
def control_quantile(values: np.ndarray, control_values: np.ndarray) -> np.ndarray:
    """各遺伝子の値が、対照遺伝子の分布の中で何分位に当たるかを返す（0〜1）。

    対照自身についても「自分より小さい対照の割合」で測る（自分を除かない。
    500 個あれば 1/500 の違いしか生まないため）。
    """
    v = np.asarray(values, dtype=float)
    c = np.sort(np.asarray(control_values, dtype=float))
    return np.searchsorted(c, v, side="left") / len(c)


def _quantile_frame(gene_table: pd.DataFrame, covariates, control_genes) -> pd.DataFrame:
    ctrl = gene_table[gene_table.gene.isin(set(control_genes))]
    if len(ctrl) < 30:
        raise ValueError(f"対照遺伝子が少なすぎる: {len(ctrl)} 個")
    out = pd.DataFrame({"gene": gene_table.gene.to_numpy()})
    for v in covariates:
        if v not in gene_table.columns:
            raise KeyError(f"交絡の候補の列がない: {v}")
        out[v] = control_quantile(gene_table[v].to_numpy(), ctrl[v].to_numpy())
    return out


# --------------------------------------------------------------------------- 回帰
def _ols(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """最小二乗。係数と係数の共分散行列（定数項込み）を返す。"""
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = XtX_inv @ (X.T @ y)
    resid = y - X @ beta
    dof = max(len(y) - X.shape[1], 1)
    s2 = float(resid @ resid) / dof
    return beta, s2 * XtX_inv


def _fit(q_ctrl: pd.DataFrame, z_ctrl: np.ndarray, covariates) -> dict:
    """対照だけで、同時回帰と周辺回帰の傾きを当てる。"""
    n = len(z_ctrl)
    Xj = np.column_stack([np.ones(n)] + [q_ctrl[v].to_numpy() for v in covariates])
    bj, Vj = _ols(Xj, z_ctrl)
    marg = {}
    for v in covariates:
        Xm = np.column_stack([np.ones(n), q_ctrl[v].to_numpy()])
        bm, Vm = _ols(Xm, z_ctrl)
        marg[v] = (float(bm[1]), float(np.sqrt(Vm[1, 1])))
    return {"beta_joint": bj[1:], "cov_joint": Vj[1:, 1:], "marginal": marg, "n": n}


# --------------------------------------------------------------------------- 本体
def bias_budget(gene_table: pd.DataFrame, set_genes, control_genes,
                covariates=CANONICAL_COVARIATES, budget: float = BUDGET,
                decision_line: float = DECISION_LINE, ci_z: float = CI_Z) -> dict:
    """1 つのセットについて、交絡による歪みを見積もる。

    返り値の主な項目:
      B_joint       同時回帰による歪みの点推定（標準偏差の単位）
      upper_joint   |B_joint| の信頼区間の上側。**これが判定に使う値**
      pass          upper_joint <= budget
      sum_abs_marginal / upper_marginal  周辺回帰版（重複して数えた上限）
      covariates    候補ごとの内訳
    """
    missing = [c for c in REQUIRED_COLUMNS if c not in gene_table.columns]
    if missing:
        raise KeyError(f"呼び出し側が必要とする列が足りない: {missing}")
    covariates = list(covariates)
    sg, cg = set(set_genes), set(control_genes)
    q = _quantile_frame(gene_table, covariates, cg)
    in_ctrl = q.gene.isin(cg).to_numpy()
    in_set = q.gene.isin(sg).to_numpy()
    if in_set.sum() == 0:
        raise ValueError("セットの遺伝子が表に 1 つもない")
    z_ctrl = gene_table.z.to_numpy()[in_ctrl]
    fit = _fit(q[in_ctrl].reset_index(drop=True), z_ctrl, covariates)

    dq = np.array([float(q[v].to_numpy()[in_set].mean()) - 0.5 for v in covariates])
    bj = fit["beta_joint"] * dq
    B = float(bj.sum())
    se_B = float(np.sqrt(max(dq @ fit["cov_joint"] @ dq, 0.0)))
    upper = abs(B) + ci_z * se_B

    det = {}
    sum_abs_m, upper_m = 0.0, 0.0
    for i, v in enumerate(covariates):
        sm, se_m = fit["marginal"][v]
        bm = sm * dq[i]
        sum_abs_m += abs(bm)
        upper_m += abs(bm) + ci_z * se_m * abs(dq[i])
        det[v] = {
            "セットの平均分位": round(float(dq[i] + 0.5), 4),
            "分位のずれ": round(float(dq[i]), 4),
            "同時回帰の傾き": round(float(fit["beta_joint"][i]), 4),
            "歪み b_v（同時）": round(float(bj[i]), 4),
            "周辺回帰の傾き": round(sm, 4),
            "傾きの標準誤差（周辺）": round(se_m, 4),
            "歪み b_v（周辺）": round(float(bm), 4),
        }
    return {
        "n_set": int(in_set.sum()), "n_controls": fit["n"],
        "covariates": det,
        "B_joint": round(B, 4), "se_B": round(se_B, 4), "upper_joint": round(upper, 4),
        "sum_abs_marginal": round(sum_abs_m, 4), "upper_marginal": round(upper_m, 4),
        "budget": budget, "decision_line": decision_line,
        "budget_frac_of_decision": round(budget / decision_line, 3),
        "pass": bool(upper <= budget),
        "pass_marginal_strict": bool(upper_m <= budget),
    }


def bias_budget_seeds(tables: dict, set_genes, control_genes, min_seeds=None, **kw) -> dict:
    """seed ごとの表に当てて、全 seed（既定）で合格かを返す。"""
    per = {int(s): bias_budget(t, set_genes, control_genes, **kw) for s, t in tables.items()}
    need = len(per) if min_seeds is None else min_seeds
    return {"per_seed": per,
            "n_pass": int(sum(r["pass"] for r in per.values())),
            "n_seeds": len(per), "min_seeds": need,
            "upper_joint": [per[s]["upper_joint"] for s in sorted(per)],
            "B_joint": [per[s]["B_joint"] for s in sorted(per)],
            "pass": bool(sum(r["pass"] for r in per.values()) >= need)}


# --------------------------------------------------------------------------- 自己検査
def _synth(rng, n_ctrl=500, n_extra=220, rho=0.85, effect=0.0, noise=1.0,
           driver="dnorm_med"):
    """交絡の候補が強く重なっている状況の合成データを作る。

    候補 5 つはすべて 1 つの潜在変数 u から作る（互いの順位相関はおよそ 0.85 で、
    実測の「移動量と元の順位 0.878 / 移動量と相対位置 0.787」に近い）。

    `effect` は「交絡 1 単位（分位で 1）あたり z が何標準偏差動くか」。
    `driver` が候補の列名なら**見えている交絡**、"_u" なら**見えていない潜在の交絡**
    （候補はその雑音まじりの代理にすぎない）。
    """
    n = n_ctrl + n_extra
    u = rng.normal(size=n)
    cols = {v: rho * u + np.sqrt(1 - rho ** 2) * rng.normal(size=n)
            for v in CANONICAL_COVARIATES}
    truth = np.zeros(n)
    if effect:
        d = u if driver == "_u" else cols[driver]
        qd = (np.argsort(np.argsort(d)) + 0.5) / n
        truth = effect * (qd - 0.5)
    z = truth + noise * rng.normal(size=n)
    t = pd.DataFrame({"gene": [f"G{i:04d}" for i in range(n)], "z": z,
                      "_truth": truth, **cols})
    t["_u"] = u
    return t, list(t.gene[:n_ctrl])


def _standardize(t, ctrl):
    """対照で標準化する。`_truth`（雑音を抜いた本当の期待値）も同じ変換をかける。
    対照だけで作った塊では、この `_truth` の平均が**真の歪み**である。"""
    c = t[t.gene.isin(set(ctrl))]
    m, sd = c.z.mean(), c.z.std(ddof=1)
    t = t.copy()
    t["z"] = (t.z - m) / sd
    t["_truth"] = (t._truth - c._truth.mean()) / sd
    return t


def _skewed_set(t, ctrl, k, p, by="_u"):
    """対照だけから、性質 `by` の分位が平均 p になる k 個の塊を取る。"""
    c = t[t.gene.isin(set(ctrl))].sort_values(by).reset_index(drop=True)
    st = int(round(p * len(c) - k / 2))
    st = max(0, min(len(c) - k, st))
    return list(c.gene[st:st + k])


def self_test(verbose=True) -> dict:
    rng = np.random.default_rng(20261002)
    out = {}
    K = 119

    # ---------------------------------------------------------------- (0) 列
    t, ctrl = _synth(rng)
    t = _standardize(t, ctrl)
    r = bias_budget(t, _skewed_set(t, ctrl, K, 0.5), ctrl)
    need = ("B_joint", "se_B", "upper_joint", "pass", "sum_abs_marginal",
            "upper_marginal", "covariates", "budget")
    miss = [k for k in need if k not in r]
    assert not miss, f"返り値に必要な項目が足りない: {miss}"
    for v in CANONICAL_COVARIATES:
        assert v in r["covariates"], f"内訳に {v} がない"
    try:
        bias_budget(t.drop(columns=["z"]), ctrl[:K], ctrl)
        raise AssertionError("列が足りないのに通ってしまった")
    except KeyError:
        pass
    out["(0) 列の検査"] = "OK（必要な項目がそろい、列が欠けると KeyError で止まる）"
    if verbose:
        print(f"(0) 列と返り値  OK  ({len(need)} 項目 / 候補 {len(CANONICAL_COVARIATES)} 個)")

    # ------------------------------------- (1) きれいなデータ（交絡の無い指標 60 本）の合格率
    #   指標が交絡の候補とまったく関係しない（effect=0）とき、歪みは真に 0。
    #   セットは「偏っていない無作為な 119 個」と「分位 0.75 に偏った 119 個」の 2 通り。
    clean = {"無作為な 119 個": [], "分位 0.75 に偏った 119 個": []}
    clean_m = {k: [] for k in clean}
    for _ in range(60):
        t, ctrl = _synth(rng, effect=0.0)
        t = _standardize(t, ctrl)
        for label, p in (("無作為な 119 個", None), ("分位 0.75 に偏った 119 個", 0.75)):
            gs = (list(rng.choice(ctrl, K, replace=False)) if p is None
                  else _skewed_set(t, ctrl, K, p))
            r = bias_budget(t, gs, ctrl)
            clean[label].append(r["pass"])
            clean_m[label].append(r["pass_marginal_strict"])
    out["(1) きれいなデータ（60 本）の合格率"] = {
        k: {"同時回帰（判定に使う版）": round(float(np.mean(v)), 3),
            "周辺回帰 Σ|b_v|（参考）": round(float(np.mean(clean_m[k])), 3)}
        for k, v in clean.items()}
    if verbose:
        print("(1) きれいなデータ 60 本の合格率 "
              + json.dumps(out["(1) きれいなデータ（60 本）の合格率"], ensure_ascii=False))

    # --------------------------------- (2) 既知の歪みを混ぜたとき、見積もりが実際と一致するか
    #   対照だけで作った塊は**中身が中立**なので、観測されるセット平均 z がそのまま
    #   「実際の歪み」である。見積もり B と比べる。
    def recover(driver, reps=40):
        rows = []
        for eff in (0.0, 0.4, 0.8, 1.6):
            for p in (0.1, 0.25, 0.5, 0.75, 0.9):
                est, act, est_m = [], [], []
                for _ in range(reps):
                    t, ctrl = _synth(rng, effect=eff, driver=driver)
                    t = _standardize(t, ctrl)
                    gs = _skewed_set(t, ctrl, K, p)
                    r = bias_budget(t, gs, ctrl)
                    est.append(r["B_joint"])
                    est_m.append(r["sum_abs_marginal"])
                    act.append(float(t.set_index("gene")._truth.reindex(gs).mean()))
                a, e = np.mean(act), np.mean(est)
                rows.append({"混ぜた歪みの強さ": eff, "セットの露出 p": p,
                             "実際の歪み（平均）": round(float(a), 4),
                             "見積もり B（平均）": round(float(e), 4),
                             "食い違い（平均）": round(float(e - a), 4),
                             "食い違いの標準偏差":
                                 round(float(np.std(np.array(est) - np.array(act))), 4),
                             "周辺回帰 Σ|b_v| ÷ 実際の歪み":
                                 round(float(np.mean(est_m) / abs(a)), 2)
                                 if abs(a) > 0.02 else None})
        return rows

    rows = recover("dnorm_med")
    out["(2a) 見えている交絡での再現"] = rows
    big = [r for r in rows if abs(r["実際の歪み（平均）"]) > 0.02]
    worst = max(abs(r["食い違い（平均）"]) for r in big)
    assert worst <= 0.02, f"見えている交絡で見積もりが合わない（最大 {worst:.4f}）"
    over = [r["周辺回帰 Σ|b_v| ÷ 実際の歪み"] for r in big
            if r["周辺回帰 Σ|b_v| ÷ 実際の歪み"]]
    out["(2a) 判定"] = (f"OK（歪みが 0.02 を超える {len(big)} 条件で、食い違いの最大 {worst:.4f}。"
                        f"周辺回帰 Σ|b_v| は実際の歪みの {min(over):.1f}〜{max(over):.1f} 倍）")
    if verbose:
        print(f"(2a) 見えている交絡  OK  食い違いの最大 {worst:.4f}（{len(big)} 条件）"
              f" / 周辺回帰 Σ|b_v| は実際の {min(over):.1f}〜{max(over):.1f} 倍")

    rows2 = recover("_u")
    out["(2b) 見えていない潜在の交絡での再現"] = rows2
    big2 = [r for r in rows2 if abs(r["実際の歪み（平均）"]) > 0.02]
    worst2 = max(abs(r["食い違い（平均）"]) for r in big2)
    ratio = [abs(r["見積もり B（平均）"]) / abs(r["実際の歪み（平均）"]) for r in big2]
    out["(2b) 回収率"] = {"最小": round(min(ratio), 3), "最大": round(max(ratio), 3),
                          "食い違いの絶対値の最大": round(float(worst2), 4)}
    out["(2b) 判定"] = (
        f"候補が雑音まじりの代理（潜在の交絡との順位相関 0.85）にすぎない場合、"
        f"見積もりは実際の歪みの {min(ratio):.0%}〜{max(ratio):.0%} しか回収しない"
        f"（取りこぼしの最大 {worst2:.4f} 標準偏差）。"
        f"**この道具は「列挙した候補で説明できる歪み」しか見積もらない**という限界。")
    assert min(ratio) >= 0.55, f"潜在の交絡の回収率が低すぎる（{min(ratio):.3f}）"
    if verbose:
        print(f"(2b) 見えていない潜在の交絡  回収率 {min(ratio):.0%}〜{max(ratio):.0%}、"
              f"取りこぼしの最大 {worst2:.4f}（限界として記録）")

    # ------------------------------------------------- (3) 信頼区間が本当に 95% を覆うか
    cov = []
    for _ in range(300):
        t, ctrl = _synth(rng, effect=0.8)
        t = _standardize(t, ctrl)
        gs = _skewed_set(t, ctrl, K, 0.8)
        r = bias_budget(t, gs, ctrl)
        act = float(t.set_index("gene")._truth.reindex(gs).mean())
        cov.append(abs(r["B_joint"] - act) <= CI_Z * r["se_B"])
    out["(3) 信頼区間の被覆率（目標 0.95）"] = round(float(np.mean(cov)), 3)
    if verbose:
        print(f"(3) 信頼区間の被覆率 {np.mean(cov):.3f}（目標 0.95）")

    # ------------------------------------------------- (4) 感度: 本物の歪みで落ちるか
    sens = []
    for eff in (0.0, 0.2, 0.4, 0.8, 1.6):
        for p in (0.6, 0.75, 0.9):
            ps, bs = [], []
            for _ in range(40):
                t, ctrl = _synth(rng, effect=eff)
                t = _standardize(t, ctrl)
                gs = _skewed_set(t, ctrl, K, p)
                r = bias_budget(t, gs, ctrl)
                ps.append(r["pass"])
                bs.append(abs(float(t.set_index("gene")._truth.reindex(gs).mean())))
            sens.append({"混ぜた歪みの強さ": eff, "露出 p": p,
                         "実際の歪みの大きさ（平均）": round(float(np.mean(bs)), 4),
                         "合格率": round(float(np.mean(ps)), 3)})
    out["(4) 感度"] = sens
    # 実際の歪みが予算 0.10 の 2 倍を超える条件では、ほぼ必ず落ちること
    bad = [r for r in sens if r["実際の歪みの大きさ（平均）"] > 2 * BUDGET and r["合格率"] > 0.1]
    assert not bad, f"大きな歪みを見逃している: {bad}"
    out["(4) 判定"] = "OK（実際の歪みが予算の 2 倍を超える条件の合格率はすべて 0.1 以下）"
    if verbose:
        print("(4) 感度  OK  " + json.dumps(sens, ensure_ascii=False))

    # ------------------- (5) 「どこまで偏ったセットなら、きれいな指標で合格と言えるか」
    #   交絡が無い指標でも、セットが大きく偏っていると**歪みが小さいと証明できない**
    #   （対照 500 個から測った傾きの不確かさが予算を食い尽くす）。その限界を測る。
    lim = []
    for p in (0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9):
        ps, dqs = [], []
        for _ in range(60):
            t, ctrl = _synth(rng, effect=0.0)
            t = _standardize(t, ctrl)
            gs = _skewed_set(t, ctrl, K, p)
            r = bias_budget(t, gs, ctrl)
            ps.append(r["pass"])
            dqs.append(np.mean([abs(d["分位のずれ"]) for d in r["covariates"].values()]))
        lim.append({"露出 p": p, "分位のずれの平均絶対値": round(float(np.mean(dqs)), 3),
                    "合格率": round(float(np.mean(ps)), 3)})
    ok = [r for r in lim if r["合格率"] >= 0.80]
    out["(5) 証明できる偏りの限界（きれいな指標）"] = lim
    out["(5) 判定"] = (
        f"交絡が無い指標でも、分位のずれの平均絶対値が約 "
        f"{max(r['分位のずれの平均絶対値'] for r in ok):.2f} を超えると合格率が 80% を割る。"
        f"無作為な 119 個のセットのずれは標準偏差 0.0265 なので十分内側だが、"
        f"**偏ったセットは「歪みが小さい」と証明できない**。")
    if verbose:
        print("(5) 証明できる偏りの限界  " + json.dumps(lim, ensure_ascii=False))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="交絡による歪みの直接見積もり")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--out", default=None, help="自己検査の結果を書き出す JSON")
    a = ap.parse_args()
    if not a.self_test:
        ap.error("この道具は呼び出して使う。単体では --self-test のみ")
    r = self_test()
    if a.out:
        with open(a.out, "w") as f:
            json.dump(r, f, ensure_ascii=False, indent=1, default=float)
        print(f"書き出し: {a.out}")
    print("自己検査 すべて合格")


if __name__ == "__main__":
    sys.exit(main())

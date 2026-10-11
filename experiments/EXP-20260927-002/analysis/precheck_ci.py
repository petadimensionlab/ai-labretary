#!/usr/bin/env python3
"""EXP-20260927-002: 事前チェックの AUC・信頼区間・交絡チェック（plan.md §4, §5）。

run_isp.py --precheck が出した CLS 埋め込み（results-precheck/seed{S}/precheck/cell_embs.csv）を
再利用する。**追加の forward は行わない。**

  --seed S   : seed ごとの計算 → results-precheck/seed{S}/precheck_ci.json
  --summary  : 3 seed をまとめ、config.yaml の criteria を機械的に当てはめる
               → results-precheck/summary.json

判定は事前登録（plan.md §5）のとおり。仮説の Status は変えない（G1 は PI）。
"""
from __future__ import annotations
import argparse, json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[3]


def rel(p: str) -> Path:
    p = Path(p)
    return p if p.is_absolute() else ROOT / p


def res_dir(cfg: dict, seed: int) -> Path:
    return rel(cfg["output"]["experiment_dir"]) / cfg["output"]["results_dir"] / f"seed{seed}"


def emb_dataset_path(cfg: dict, seed: int) -> Path:
    se = cfg["state_embs"]
    base = f"state-embs-seed{seed}-n{se['max_cells_per_state']}-d{se['max_cells_per_donor']}"
    # 結合の不具合修正（2026-09-27, plan.md 末尾「逸脱の記録」参照）:
    # EmbExtractor.extract_embs は内部で downsample_and_sort（length 降順・安定ソート）により行を
    # 並べ替えるため、(state, donor) だけの位置合わせは同一 donor 内の複数細胞を取り違える恐れがある
    # （実際、length・donor が同じでも n_counts/subclass_l2 が異なる細胞が 3〜4% 存在することを確認済み）。
    # run_isp.py 側で row_id 列を付けた dataset（*-rowid.dataset）を使い、cell_embs.csv 側の row_id と
    # ID 結合する。
    rowid_path = rel(cfg["output"]["experiment_dir"]) / "work" / f"{base}-rowid.dataset"
    if rowid_path.exists():
        return rowid_path
    raise SystemExit(f"{rowid_path} がありません。run_isp.py --precheck を row_id 付きで再実行してください。")


# ------------------------------------------------------------------ 読み込み
def load_cells(cfg: dict, seed: int) -> tuple[np.ndarray, pd.DataFrame]:
    """埋め込み（cell_embs.csv）と、同じ細胞のメタ情報（*-rowid.dataset）を row_id で結合して返す。
    (state, donor) だけの位置合わせは同一 donor 内の複数細胞を一意に対応付けられないため使わない
    （2026-09-27 の結合の不具合修正。plan.md 末尾「逸脱の記録」参照）。"""
    key, dkey = cfg["data"]["state_key"], cfg["cell_selection"]["donor_key"]
    emb = pd.read_csv(res_dir(cfg, seed) / "precheck" / "cell_embs.csv", index_col=0)
    if "row_id" not in emb.columns:
        raise SystemExit("cell_embs.csv に row_id 列がありません。run_isp.py --precheck を再実行してください。")
    X = emb.drop(columns=[key, dkey, "row_id"]).to_numpy(dtype=np.float64)

    from datasets import load_from_disk
    ds = load_from_disk(str(emb_dataset_path(cfg, seed)))   # 読むだけ（map/select/filter はしない）
    want = [key, dkey, "length", "n_counts", "subclass_l2", "region", "disease", "row_id"]
    full = pd.DataFrame({c: list(ds[c]) for c in want if c in ds.column_names})
    if full["row_id"].duplicated().any():
        raise SystemExit("dataset 側の row_id が重複しています")
    if emb["row_id"].duplicated().any():
        raise SystemExit("cell_embs.csv 側の row_id が重複しています")
    missing = set(emb["row_id"]) - set(full["row_id"])
    if missing:
        raise SystemExit(f"cell_embs.csv の row_id が dataset にありません: {sorted(missing)[:10]}…")

    # row_id で ID 結合（位置ではなく ID で突き合わせる）
    meta = full.set_index("row_id").loc[emb["row_id"].to_numpy()].reset_index()
    if len(meta) != len(emb):
        raise SystemExit(f"細胞数が一致しません: cell_embs {len(emb)} / dataset {len(meta)}")
    same = (meta[key].astype(str).to_numpy() == emb[key].astype(str).to_numpy()).all() and \
           (meta[dkey].astype(str).to_numpy() == emb[dkey].astype(str).to_numpy()).all()
    if not same:
        raise SystemExit("row_id で結合した後も cell_embs.csv と dataset の (状態, donor) が一致しません。"
                         "結合ロジックを見直してください")
    # 件数の一致を明示的に assert（コーディネーター指示。donor_id の dtype 差異（str/int）を吸収するため
    # 比較前に文字列化する）
    assert len(meta) == len(emb) == X.shape[0], "行数が一致しません"
    meta_counts = meta.assign(**{dkey: meta[dkey].astype(str)}).groupby([key, dkey]).size().sort_index()
    emb_counts = emb.assign(**{dkey: emb[dkey].astype(str)}).groupby([key, dkey]).size().sort_index()
    assert meta_counts.index.equals(emb_counts.index) and (meta_counts.to_numpy() == emb_counts.to_numpy()).all(), \
        "(状態, donor) ごとの件数が一致しません"
    meta = meta.reset_index(drop=True)
    meta["log_len"] = np.log(meta["length"].to_numpy(dtype=float))
    meta["log_counts"] = np.log1p(meta["n_counts"].to_numpy(dtype=float))
    return X, meta


# ------------------------------------------------------------------ AUC
def auc_cv(X: np.ndarray, y: np.ndarray, g: np.ndarray, folds: int, C: float) -> dict:
    """donor 単位 GroupKFold の fold 平均 AUC と out-of-fold 予測。"""
    aucs, oof = [], np.full(len(y), np.nan)
    used = 0
    for tr, te in GroupKFold(n_splits=folds).split(X, y, groups=g):
        if len(set(y[te])) < 2 or len(set(y[tr])) < 2:
            continue
        clf = make_pipeline(StandardScaler(), LogisticRegression(C=C, max_iter=5000))
        clf.fit(X[tr], y[tr])
        p = clf.predict_proba(X[te])[:, 1]
        oof[te] = p
        aucs.append(roc_auc_score(y[te], p))
        used += 1
    return {"auc_mean": round(float(np.mean(aucs)), 4) if aucs else None,
            "auc_sd": round(float(np.std(aucs, ddof=1)), 4) if len(aucs) > 1 else None,
            "fold_aucs": [round(float(a), 4) for a in aucs], "n_folds_used": used,
            "n_cells": int(len(y)), "n_pos": int(y.sum()), "n_donors": int(len(set(g))),
            "_oof": oof}


def boot_ci(oof: np.ndarray, y: np.ndarray, g: np.ndarray, n: int, seed: int) -> dict:
    """out-of-fold 予測に対する donor 単位ブートストラップ 95% CI。"""
    ok = ~np.isnan(oof)
    oof, y, g = oof[ok], y[ok], g[ok]
    donors = np.array(sorted(set(g)))
    idx = {d: np.where(g == d)[0] for d in donors}
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n):
        pick = rng.choice(donors, size=len(donors), replace=True)
        sel = np.concatenate([idx[d] for d in pick])
        if len(set(y[sel])) < 2:
            continue
        vals.append(roc_auc_score(y[sel], oof[sel]))
    if not vals:
        return {"pooled_auc": None, "ci95": None, "n_boot_used": 0}
    return {"pooled_auc": round(float(roc_auc_score(y, oof)), 4),
            "ci95": [round(float(np.percentile(vals, 2.5)), 4), round(float(np.percentile(vals, 97.5)), 4)],
            "n_boot_used": len(vals)}


def perm_test(X, y, g, folds, C, n_perm: int, seed: int, obs_auc: float) -> dict:
    """donor の中でラベルを入れ替える並べ替え検定（donor 構成の効果は残す = 保守的）。"""
    rng = np.random.default_rng(seed)
    donors = sorted(set(g))
    idx = {d: np.where(g == d)[0] for d in donors}
    null = []
    for _ in range(n_perm):
        yp = y.copy()
        for d in donors:
            i = idx[d]
            yp[i] = rng.permutation(y[i])
        r = auc_cv(X, yp, g, folds, C)
        if r["auc_mean"] is not None:
            null.append(r["auc_mean"])
    if not null:
        return {"p": None, "n_perm_used": 0}
    p = (int(np.sum(np.array(null) >= obs_auc)) + 1) / (len(null) + 1)
    return {"p": round(float(p), 4), "n_perm_used": len(null),
            "null_mean": round(float(np.mean(null)), 4), "null_p95": round(float(np.percentile(null, 95)), 4)}


# ------------------------------------------------------------------ 深さマッチ
def depth_match(meta: pd.DataFrame, pos: np.ndarray, neg: np.ndarray, tol: float) -> np.ndarray:
    """log(トークン長)の最近傍 1:1 マッチ（置換なし）。埋め込みは使わない。"""
    pi = np.where(pos)[0]
    ni = np.where(neg)[0]
    order = np.argsort(meta["log_len"].to_numpy()[ni], kind="stable")
    ni = ni[order]
    nv = meta["log_len"].to_numpy()[ni]
    used = np.zeros(len(ni), dtype=bool)
    keep = []
    for i in pi:
        v = meta["log_len"].iat[i]
        cand = np.where((~used) & (np.abs(nv - v) <= tol))[0]
        if len(cand):
            j = cand[np.argmin(np.abs(nv[cand] - v))]
            used[j] = True
            keep += [i, int(ni[j])]
    m = np.zeros(len(meta), dtype=bool)
    m[keep] = True
    return m


# ------------------------------------------------------------------ 1 seed
def run_seed(cfg: dict, seed: int) -> dict:
    key, dkey = cfg["data"]["state_key"], cfg["cell_selection"]["donor_key"]
    cr, pc = cfg["criteria"], cfg["precheck"]
    folds, C = pc["folds"], pc["C"]
    tol = cr["depth_match"]["log_tolerance"]
    X, meta = load_cells(cfg, seed)
    st, gl, alt = cfg["states"]["start_state"], cfg["states"]["goal_state"], cfg["states"]["alt_states"][0]
    state = meta[key].astype(str).to_numpy()
    donor = meta[dkey].astype(str).to_numpy()
    out = {"seed": seed, "n_cells_total": int(len(meta)),
           "n_by_state": pd.Series(state).value_counts().to_dict(),
           "median_length_by_state": meta.groupby(key)["length"].median().to_dict()}

    def pair(a: str, b: str, sub: np.ndarray | None = None, matched: bool = True, label: str = "") -> dict:
        m = np.isin(state, [a, b])
        if sub is not None:
            m &= sub
        if matched:
            m = m & depth_match(meta, m & (state == b), m & (state == a), tol)
        y = (state[m] == b).astype(int)          # goal = 1
        if len(set(y)) < 2:
            return {"label": label, "error": "片方のクラスしかありません"}
        r = auc_cv(X[m], y, donor[m], folds, C)
        oof = r.pop("_oof")
        r |= boot_ci(oof, y, donor[m], cr["bootstrap_n"], cr["bootstrap_seed"])
        d = np.c_[meta.loc[m, "log_counts"].to_numpy(), meta.loc[m, "log_len"].to_numpy()]
        rd = auc_cv(d, y, donor[m], folds, C)
        rd.pop("_oof")
        r["depth_only_auc"] = rd["auc_mean"]
        r["depth_only_folds"] = rd["fold_aucs"]
        r["median_length"] = {a: float(meta.loc[m & (state == a), "length"].median()),
                              b: float(meta.loc[m & (state == b), "length"].median())}
        r["label"] = label
        r["_mask"] = m
        r["_y"] = y
        return r

    # C1: 主判定（深さマッチ後の start vs goal）
    prim = pair(st, gl, matched=True, label="goal_vs_start_depth_matched")
    m, y = prim.pop("_mask"), prim.pop("_y")
    prim["permutation"] = perm_test(X[m], y, donor[m], folds, C,
                                    cr["permutations"], cr["permutation_seed"], prim["auc_mean"])
    out["primary"] = prim

    # 参考: マッチ前 / start vs alt
    for label, kw in [("goal_vs_start_unmatched", dict(matched=False)),
                      ("start_vs_alt_depth_matched", dict(matched=True))]:
        a, b = (st, gl) if "goal" in label else (alt, st)
        r = pair(a, b, label=label, **kw)
        r.pop("_mask", None)
        r.pop("_y", None)
        out.setdefault("reference", {})[label] = r

    # C3: 陽性対照（FIB vs M-FIB、全状態をまとめて）
    p = cr["C3_positive_control"]
    sc = meta["subclass_l2"].astype(str).to_numpy()
    m3 = np.isin(sc, [p["a"], p["b"]])
    r3 = auc_cv(X[m3], (sc[m3] == p["b"]).astype(int), donor[m3], folds, C)
    r3.pop("_oof")
    out["positive_control"] = {"label": f"{p['a']}_vs_{p['b']}", "min_auc": p["min_auc"], **r3}

    # 感度分析（参考値。判定には使わない）
    sens = {}
    s = cr.get("sensitivity") or {}
    for name, rule in s.items():
        if name == "within_start_donors":
            if not rule:
                continue
            sub = np.isin(donor, list(set(donor[state == st])))   # start と同じ donor の細胞だけ
        else:
            col, vals = next(iter(rule.items()))                  # 例: {region: [C]}, {disease: [AKI, CKD]}
            sub = np.isin(meta[col].astype(str).to_numpy(), [str(v) for v in vals])
        r = pair(st, gl, sub=sub, matched=True, label=name)
        r.pop("_mask", None)
        r.pop("_y", None)
        sens[name] = r
    out["sensitivity"] = sens

    p = res_dir(cfg, seed) / "precheck_ci.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=2, ensure_ascii=False, default=float) + "\n")
    print(f"→ {p}")
    return out


# ------------------------------------------------------------------ まとめ
def summarize(cfg: dict) -> dict:
    cr = cfg["criteria"]
    seeds = cfg["robustness"]["seeds"]
    per = {}
    for s in seeds:
        f = res_dir(cfg, s) / "precheck_ci.json"
        if f.exists():
            per[s] = json.loads(f.read_text())
    if not per:
        raise SystemExit("precheck_ci.json がありません（--seed を先に実行してください）")
    a = [v["primary"]["auc_mean"] for v in per.values()]
    lo = [v["primary"]["ci95"][0] for v in per.values() if v["primary"].get("ci95")]
    pv = [v["primary"]["permutation"]["p"] for v in per.values()]
    dep = [v["primary"]["depth_only_auc"] for v in per.values()]
    pos = [v["positive_control"]["auc_mean"] for v in per.values()]
    c1 = bool(np.mean(a) >= cr["C1_min_auc"] and sum(x >= cr["C1_min_auc"] for x in a) >= cr["C1_min_seeds_passing"])
    c1b = bool(lo and np.median(lo) >= cr["C1b_min_ci_lower"])
    c1c = bool(sum((x is not None and x < cr["C1c_max_p"]) for x in pv) >= cr["C1_min_seeds_passing"])
    c2 = bool(np.mean(dep) <= cr["C2_max_depth_auc"] and (np.mean(a) - np.mean(dep)) >= cr["C2_min_margin"])
    c3 = bool(np.mean(pos) >= cr["C3_positive_control"]["min_auc"])
    if c1 and c1b and c1c and c2 and c3:
        verdict = "合格（ISP の前提を満たす。G1・G2 の判断は PI）"
    elif not c1:
        verdict = "不合格（ISP は行わない。goal の定義を見直す案を PI に提出）"
    else:
        verdict = "判定不能（C2 または C3 を満たさない）"
    out = {"seeds": list(per), "primary_auc_per_seed": a, "primary_auc_mean": round(float(np.mean(a)), 4),
           "ci95_lower_per_seed": lo, "permutation_p_per_seed": pv, "depth_only_auc_per_seed": dep,
           "positive_control_auc_per_seed": pos,
           "criteria": {"C1": c1, "C1b": c1b, "C1c": c1c, "C2": c2, "C3": c3},
           "verdict": verdict,
           "note": "仮説の Status は変えない。G1（採択）と G3（結論）の判断は PI。"}
    f = rel(cfg["output"]["experiment_dir"]) / cfg["output"]["results_dir"] / "summary.json"
    f.write_text(json.dumps(out, indent=2, ensure_ascii=False, default=float) + "\n")
    print(json.dumps(out, indent=2, ensure_ascii=False))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--seed", type=int)
    ap.add_argument("--summary", action="store_true")
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    if args.summary:
        summarize(cfg)
    elif args.seed is not None:
        run_seed(cfg, args.seed)
    else:
        raise SystemExit("--seed S か --summary のどちらかを指定してください")


if __name__ == "__main__":
    main()

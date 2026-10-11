#!/usr/bin/env python3
"""着手前チェック — 摂動を 1 回も回さずに「その実験の目標の軸が測れる設計か」を判定する。

## なぜあるか

EXP-20260927-004 / EXP-20260928-001 / EXP-20260928-002 / EXP-20260929-001 の 4 実験が、
**目標の軸の目盛りが立っていない**（目標の定義に使った遺伝子を動かしても目標の軸で動かない）
という同じ理由で止まった。EXP-20260929-001 ではこれが判明するまでに機体 11.3 時間を使ったが、
**判定に必要な量はどれも摂動を 1 回も使わずに計算できた**（validation.md §2.3.2 / §2.4）。

本スクリプトは、その量を計算して **合格 / 判定保留 / 不合格** を返す。

- **不合格**: この設計では目標の軸が測れない。**摂動の本番に着手しない。**
  計画を作り直して再度このチェックを通す。
- **判定保留**: 測れる可能性はあるが、弱点が特定されている。
  **弱点と、それが結論に与える影響を `plan.md` に明記した上でのみ着手してよい**（PI の判断を仰ぐ）。
- **合格**: 弱点は検出されなかった。

**このチェックは事前登録の一部である。** 結果を見る前に実行し、出力 JSON を事前登録のハッシュに含める。
判定が出たあとで閾値を変えてはならない。

## 測るもの（摂動を使わない）

| # | 量 | 何を捕まえるか |
|---|---|---|
| 1 | **陽性対照の遺伝子数 `K_positive_control`（最も強い条件）** | 目標の軸の健全性チェックに検出力があるか |
| 2 | マーカーの検出率 `marker_detection_median` | そもそも測定法（単一核 RNA-seq 等）で拾えているか |
| 3 | マーカー群のまとまり `marker_corr_median` | 目標を定義するマーカーが 1 つの細胞状態を指しているか |
| 4 | 個々のマーカーの偏り `skew_*` | セットの値が 1 個の遺伝子に依存していないか |
| 5 | 軸の長さの比 `axis_length_ratio` | 目標が、うまく測れている参照状態に比べてどれだけ近いか |

1 は遺伝子リストだけ、2・3・4 は発現データだけ、5 は状態の重心だけから出る。

**順番には意味がある。** 過去 4 実験を分離したのは 1（K）であり、5（軸の長さの比）は
うまくいった実験と失敗した実験を分けなかった（0.442〜0.637 で重なる）。
当初 Labretary は 5 を第 1 の指標に挙げたが、実測では止める力が無かった（2026-09-30 訂正）。

## 使い方

    .venv/bin/python scripts/analysis/axis_feasibility.py --config <設定 yaml> [--out <json>]

設定ファイルの雛形は `templates/experiment/config.yaml` の `axis_feasibility:` 節。
実験の `config.yaml` をそのまま渡してもよい（`axis_feasibility:` 節だけを読む）。

終了コード: 0 = 合格 / 1 = 判定保留 / 2 = 不合格 / 3 = 実行エラー。
`--no-exit-code` を付けると判定によらず 0 で終わる。

`data/raw/` を含め、入力ファイルは**読み取りのみ**。書き込むのは `--out` の 1 ファイルだけ。
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]

# ---------------------------------------------------------------------------
# 既定の合否の目安
#
# **この値は EXP-20260927-004 / EXP-20260928-001 / EXP-20260928-002 / EXP-20260929-001 の
# 実測値に合わせて決めた（2026-09-30）。** 当てはめた実測値は次のとおり。
#
#   | 軸 / セット | 軸の長さの比 | 相関の中央値 | 検出率の中央値 | K | 実際の結果 |
#   |---|---|---|---|---|---|
#   | 筋線維芽細胞の軸（毎回働いた） | 1.00（基準） | 0.068 / 0.016 | 0.112 / 0.095 | 24 | 合格 |
#   | EXP-20260928-001 目標の軸 | 0.452〜0.498 | 0.019 | 0.088 | 24 | 限定つき合格 |
#   | EXP-20260927-004 目標の軸 | 0.608〜0.637 | 0.004 | 0.039 | 1 | 偶然と区別できず |
#   | EXP-20260928-002 目標の軸 | 0.448〜0.503 | 0.004 | 0.039 | 0 | 陽性対照を置けず |
#   | EXP-20260929-001 目標の軸 | 0.442〜0.450 | 0.009 | 0.042 | 7 | 不合格 |
#
# 要点（詳しい理由は judge() の説明を読むこと）:
#   - **止める条件は「実際に検定に使う陽性対照のセット」に当てる。**
#     目標を定義したマーカーの値は 4 実験を分離しなかった（28-001 と 28-002 は同じ定義で結果が違う）。
#   - 目標の軸で陽性対照が働いたのはこのプロジェクトで 1 回だけ（EXP-20260928-001、限定つき）。
#     その実験の実測値を「判定保留で通る」水準に置き、失敗した 3 実験を「不合格」に落とす。
#   - **軸の長さの比は 4 実験を分離しない**（0.442〜0.637）。うまくいった EXP-20260928-001 も 0.45〜0.50 で、
#     当初案の 0.70 はこのプロジェクトで一度も到達していない。**単独で止める条件には使わない。**
#   - **検出率 0.15 も到達例がない**（常に働いている筋線維芽細胞の軸でさえ 0.095〜0.112）。
#   - **相関の中央値は 0.009（不合格）対 0.019（限定つき合格）で差が 2 倍しかない。**
#     単独では止める根拠にならないので、検出率と 2 つそろったときだけ止める条件にする。
#   - **最も強い条件は陽性対照の遺伝子数 K である。** これは検出力の議論から出ており、
#     閾値の位置に恣意性がない（K = 7 → 検出下限 1.71 標準偏差。本プロジェクトで到達例 0）。
# ---------------------------------------------------------------------------
DEFAULT_THRESHOLDS = {
    # --- 不合格（着手しない）。1 つでも当たれば不合格 ---
    "block_K_min": 20,                  # 陽性対照の遺伝子数がこれ未満
    "block_detection_median_min": 0.06,  # マーカー検出率の中央値がこれ未満
    "block_axis_ratio_min": 0.35,        # 軸の長さの比がこれ未満
    # 相関と検出率が 2 つそろって低いとき（相関は単独では止めない）
    "block_pair_corr_max": 0.010,
    "block_pair_detection_max": 0.10,
    # --- 判定保留（弱点を明記した上でのみ着手） ---
    "hold_corr_median_min": 0.030,
    "hold_detection_median_min": 0.12,
    "hold_axis_ratio_min": 0.60,
    "hold_score_share_max": 0.25,        # 1 個のマーカーが得点の何割を占めるか
    "hold_detection_spread_max": 5.0,    # 検出率の 最大 / 最小
    "hold_circular_overlap_max": 0.50,   # 陽性対照と目標の定義マーカーの重なり
}


# ---------------------------------------------------------------------------
# 入出力
# ---------------------------------------------------------------------------
def rel(p) -> Path:
    q = Path(str(p))
    return q if q.is_absolute() else ROOT / q


def load_cfg(path: Path) -> dict:
    raw = yaml.safe_load(path.read_text())
    if not isinstance(raw, dict):
        raise SystemExit(f"設定を読めません: {path}")
    cfg = raw.get("axis_feasibility", raw)
    if "axis_feasibility" not in raw and "expression" not in raw:
        raise SystemExit(f"{path} に axis_feasibility 節がありません")
    return cfg


# ---------------------------------------------------------------------------
# 1. 発現データの指標（摂動を使わない）
# ---------------------------------------------------------------------------
def expression_metrics(cfg: dict) -> dict:
    """マーカー群のまとまり・検出率・偏りを、発現データだけから計算する。"""
    import anndata as ad
    import scipy.sparse as sp

    ex = cfg["expression"]
    h5 = rel(ex["h5ad"])
    if not h5.exists():
        raise SystemExit(f"発現データがありません: {h5}")
    A = ad.read_h5ad(str(h5))

    state_col = ex["state_col"]
    if state_col not in A.obs.columns:
        raise SystemExit(f"obs に列 {state_col} がありません（あるのは {list(A.obs.columns)}）")
    states = cfg["states"]
    pop_keys = cfg.get("population", ["start", "goal"])
    labels: list[str] = []
    for k in pop_keys:
        labels += [str(x) for x in states[k]["labels"]]
    mask = A.obs[state_col].astype(str).isin(labels).to_numpy()
    if mask.sum() == 0:
        raise SystemExit(f"細胞が 0 個です（{state_col} in {labels}）")

    X = A[mask].X
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(np.asarray(X))
    name_col = ex.get("gene_name_col", "gene_name")
    gene_names = (A.var[name_col].astype(str).to_numpy() if name_col in A.var.columns
                  else A.var.index.astype(str).to_numpy())

    # 目標の状態にあたる細胞（得点の偏りを測る対象）
    goal_labels = [str(x) for x in states["goal"]["labels"]]
    goal_mask = A.obs[state_col].astype(str).to_numpy()[mask]
    goal_mask = np.isin(goal_mask, goal_labels)

    total = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
    total[total == 0] = 1.0

    def one_set(symbols: list[str], label: str) -> dict:
        idx, found, missing = [], [], []
        for s in symbols:
            w = np.where(gene_names == s)[0]
            if len(w):
                idx.append(int(w[0]))
                found.append(s)
            else:
                missing.append(s)
        if len(idx) < 2:
            return {"set": label, "n_genes_requested": len(symbols), "n_genes_found": len(idx),
                    "missing": missing, "error": "遺伝子が 2 個未満で相関を計算できない"}
        C = np.asarray(X[:, idx].todense(), dtype=np.float64)
        det = (C > 0).mean(axis=0)
        n_cells_per_gene = (C > 0).sum(axis=0).astype(int)
        L = np.log1p(C / total[:, None] * 1e4)          # 1 万カウントあたり + log1p
        R = np.corrcoef(L.T)
        iu = np.triu_indices(len(idx), 1)
        off = R[iu]
        # 共検出のリフト（まばらなデータでの相関を補う指標）
        lifts = []
        D = C > 0
        for i, j in zip(*iu):
            exp = det[i] * det[j]
            lifts.append(float((D[:, i] & D[:, j]).mean() / exp) if exp > 0 else np.nan)
        # 目標の細胞での得点寄与（z 得点の和に対する 1 個ごとの取り分）
        sd = L.std(axis=0)
        Z = (L - L.mean(axis=0)) / np.where(sd == 0, 1.0, sd)
        contrib = np.clip(Z[goal_mask].sum(axis=0), 0, None) if goal_mask.any() else np.zeros(len(idx))
        share = contrib / contrib.sum() if contrib.sum() > 0 else contrib
        opposite = [found[i] for i in range(len(idx))
                    if goal_mask.any() and Z[goal_mask, i].mean() < 0]
        order = np.argsort(-share)
        return {
            "set": label,
            "n_cells": int(X.shape[0]),
            "n_genes_requested": len(symbols), "n_genes_found": len(idx), "missing": missing,
            "detection_median": float(np.median(det)),
            "detection_min": float(det.min()), "detection_max": float(det.max()),
            "detection_spread": float(det.max() / det.min()) if det.min() > 0 else float("inf"),
            "corr_median": float(np.median(off)),
            "corr_min": float(off.min()), "corr_max": float(off.max()),
            "codetection_lift_median": float(np.nanmedian(lifts)),
            "score_share_max": float(share.max()) if len(share) else float("nan"),
            "score_share_top_gene": found[int(order[0])] if len(share) else None,
            "n_opposite_direction": len(opposite), "opposite_direction": opposite,
            "per_gene": [{"symbol": found[i], "detection": float(det[i]),
                          "n_cells_detected": int(n_cells_per_gene[i]),
                          "score_share": float(share[i])} for i in np.argsort(-det)],
        }

    mk = cfg["markers"]
    gd_genes = [str(s) for s in (mk.get("goal_defining") or [])]
    # 目標を**著者の注釈など、遺伝子リストによらない形**で定義した場合は goal_defining が空になる。
    # これは正しい状態である（循環が構造的に起きない）。代わりの遺伝子リストを詰めてはいけない
    # ——詰めると circular_overlap = 0 が「重なりが無い」ではなく「別のリストを入れた」の意味になり、
    # 循環していないことの根拠にならなくなる（2026-09-30、HYP-20260930-01 の設計で実際に起きた）。
    out = {"population_labels": labels, "n_cells": int(X.shape[0])}
    out["goal_defining"] = (one_set(gd_genes, "目標を定義するマーカー") if len(gd_genes) >= 2 else
                            {"set": "目標を定義するマーカー", "n_genes_requested": len(gd_genes),
                             "genes": gd_genes,
                             "note": "目標は遺伝子リストで定義されていない（著者の注釈など）。"
                                     "定義マーカーの指標は計算しない"})
    # 陽性対照は「書いていない」と「空（置かなかった）」を区別する
    pc = [str(s) for s in mk["goal_positive_control"]] if "goal_positive_control" in mk \
        else [str(s) for s in mk["goal_defining"]]
    out["goal_positive_control"] = (one_set(pc, "目標の軸の陽性対照") if len(pc) >= 2 else
                                    {"set": "目標の軸の陽性対照", "n_genes_requested": len(pc),
                                     "n_genes_found": len(pc), "genes": pc,
                                     "error": "遺伝子が 2 個未満。まとまり・偏りを計算できない"})
    if mk.get("reference_positive_control"):
        out["reference_positive_control"] = one_set(
            [str(s) for s in mk["reference_positive_control"]], "参照の軸の陽性対照（対比）")
    if mk.get("contrast_set"):
        out["contrast_set"] = one_set([str(s) for s in mk["contrast_set"]], "対比のためのマーカー群")

    # 循環（陽性対照が目標の定義マーカーそのものか）
    # 目標を遺伝子リストで定義していないときは「重なりが 0」ではなく「そもそも問えない」。
    # 0.0 を入れると循環していない根拠のように読めるので、明示的に区別する。
    if len(gd_genes) == 0:
        out["circular_overlap"] = None
        out["circular_overlap_note"] = (
            "目標は遺伝子リストで定義されていないため、定義マーカーとの重なりは定義されない。"
            "**ただし注釈は発現のクラスタリングから作られているので、循環が完全に無いわけではない。**")
    else:
        a, b = set(gd_genes), set(pc)
        out["circular_overlap"] = len(a & b) / len(b) if b else 0.0

    # 摂動できる細胞数の下限を満たさない遺伝子
    floor = int(cfg.get("perturbation", {}).get("min_cells_per_gene", 0))
    if floor:
        out["below_cell_floor"] = [g["symbol"] for g in out["goal_positive_control"].get("per_gene", [])
                                   if g["n_cells_detected"] < floor]
        out["min_cells_per_gene"] = floor
    return out


# ---------------------------------------------------------------------------
# 2. 重心の指標（摂動を使わない）
# ---------------------------------------------------------------------------
def load_centroids(cfg: dict) -> dict[int, dict[str, np.ndarray]]:
    """seed ごとの状態重心を読む。

    mode: precomputed — 既に作ってある .npz を読む（`{seed}` を展開）
    mode: extract     — Geneformer の EmbExtractor で作る（摂動は行わない）
    """
    cc = cfg["centroids"]
    mode = cc.get("mode", "precomputed")
    seeds = [int(s) for s in cc.get("seeds", [0])]
    if mode == "extract":
        return extract_centroids(cfg, seeds)
    if mode != "precomputed":
        raise SystemExit(f"centroids.mode は precomputed か extract です（{mode}）")
    out: dict[int, dict[str, np.ndarray]] = {}
    for s in seeds:
        vec: dict[str, np.ndarray] = {}
        for f in cc["files"]:
            p = rel(str(f).format(seed=s))
            if not p.exists():
                raise SystemExit(f"重心ファイルがありません: {p}")
            z = np.load(str(p))
            vec.update({k: z[k].astype(np.float64) for k in z.files})
        out[s] = vec
    return out


def extract_centroids(cfg: dict, seeds: list[int]) -> dict[int, dict[str, np.ndarray]]:
    """トークン化済みデータから状態重心を作る（摂動は行わない）。

    `experiments/EXP-20260929-001/analysis/extra_centroids_003.py` と同じ手順。
    重心の細胞は `rows_file`（seed ごとの行番号 JSON）で指定する。
    **モデルの前向き計算だけを行い、摂動は 1 回も行わない。**
    """
    import os
    cc = cfg["centroids"]
    sys.path.insert(0, str(ROOT / "scripts/geneformer"))
    import run_isp as ri                                   # noqa: E402
    from datasets import load_from_disk                    # noqa: E402
    from geneformer import EmbExtractor                    # noqa: E402

    base = ri.load_config(rel(cc["run_config"]))
    os.environ["GENEFORMER_DEVICE"] = ri.resolve_device(base)
    ri.apply_dtype(base)
    ri.single_process_datasets(base)
    ds = load_from_disk(str(rel(cc.get("tokenized_dataset") or base["data"]["tokenized_dataset"])))
    m, p = base["model"], base.get("perturbation", {})
    out: dict[int, dict[str, np.ndarray]] = {}
    for s in seeds:
        spec = json.loads(rel(str(cc["rows_file"]).format(seed=s)).read_text())
        rows_by_state = {k: [int(i) for i in v] for k, v in spec["rows"].items()}
        for k, v in spec.get("extra_rows", {}).items():
            rows_by_state[k] = [int(i) for i in v]
        rows = sorted({i for v in rows_by_state.values() for i in v})
        work = rel(cc["work_dir"])
        work.mkdir(parents=True, exist_ok=True)
        dpath = work / f"axisfeas-cells-seed{s}.dataset"
        if not dpath.exists():
            ds.select(rows).map(lambda ex, i: {"row_id": rows[i]}, with_indices=True).save_to_disk(str(dpath))
        ee = ri.call_compat(EmbExtractor, {
            "model_type": m["type"], "num_classes": m.get("num_classes", 0), "filter_data": None,
            "max_ncells": None, "emb_mode": p.get("emb_mode", "cls").replace("_and_gene", ""),
            "emb_layer": m.get("emb_layer", -1), "emb_label": ["row_id"], "summary_stat": None,
            "forward_batch_size": p.get("forward_batch_size", 8), "nproc": p.get("nproc", 1),
            "model_version": m.get("version")}, "EmbExtractor")
        odir = work / f"axisfeas-embs-seed{s}"
        odir.mkdir(parents=True, exist_ok=True)
        embs = ee.extract_embs(str(rel(m["directory"])), str(dpath), str(odir), "cell_embs")
        X = embs.drop(columns=["row_id"]).to_numpy(dtype=np.float64)
        by_row = dict(zip(embs["row_id"].astype(int).to_numpy(), range(len(embs))))
        out[s] = {name: X[[by_row[i] for i in keep]].mean(axis=0)
                  for name, keep in rows_by_state.items()}
    return out


def centroid_metrics(cfg: dict) -> dict:
    """軸の長さを 2 通りで測る。

    - `distance_ratio`: 出発点から目標までの距離 ÷ 出発点から参照状態までの距離
      （**主の指標**。PI が定義した量）
    - `direction_ratio`: 軸の向きベクトル（自分の重心 − 他の重心の平均）の長さの比
      （摂動の指標が実際に使う向き。参考値）

    重心はすべて長さ 1 に直してから測る（指標がコサインの差なので長さは効かない）。
    """
    cents = load_centroids(cfg)
    nm = cfg["centroids"]["names"]
    start, goal, ref = nm["start"], nm["goal"], nm["reference"]
    others = [str(x) for x in nm.get("others", [])]
    per_seed = []
    for s in sorted(cents):
        v = cents[s]
        need = [start, goal, ref] + others
        miss = [n for n in need if n not in v]
        if miss:
            raise SystemExit(f"seed {s}: 重心に {miss} がありません（あるのは {sorted(v)}）")
        U = {n: v[n] / np.linalg.norm(v[n]) for n in need}
        d_goal = float(np.linalg.norm(U[start] - U[goal]))
        d_ref = float(np.linalg.norm(U[start] - U[ref]))

        def direction(n: str) -> np.ndarray:
            rest = [U[m] for m in need if m != n]
            return U[n] - np.mean(rest, axis=0)

        lg = float(np.linalg.norm(direction(goal)))
        lr = float(np.linalg.norm(direction(ref)))
        per_seed.append({
            "seed": s,
            "dist_start_to_goal": d_goal, "dist_start_to_reference": d_ref,
            "distance_ratio": d_goal / d_ref if d_ref else float("nan"),
            "dirlen_goal": lg, "dirlen_reference": lr,
            "direction_ratio": lg / lr if lr else float("nan"),
        })
    return {
        "states": {"start": start, "goal": goal, "reference": ref, "others": others},
        "per_seed": per_seed,
        "distance_ratio_min": min(x["distance_ratio"] for x in per_seed),
        "distance_ratio_median": float(np.median([x["distance_ratio"] for x in per_seed])),
        "direction_ratio_min": min(x["direction_ratio"] for x in per_seed),
        "direction_ratio_median": float(np.median([x["direction_ratio"] for x in per_seed])),
    }


# ---------------------------------------------------------------------------
# 3. 判定
# ---------------------------------------------------------------------------
def judge(expr: dict, cent: dict, K: int, th: dict) -> dict:
    """止める条件は「実際に検定に使う陽性対照のセット」に当てる。

    **なぜ定義マーカーではなく陽性対照なのか（過去 4 実験の実測から）。**
    EXP-20260928-001 は目標の状態を 3 つのケモカイン（検出率の中央値 0.039、相関 0.004）で定義したが、
    目標の軸の陽性対照には**別に選んだ 24 遺伝子**（検出率 0.088、相関 0.019）を使い、
    その健全性チェックに（限定つきで）合格した。一方 EXP-20260929-001 は定義マーカー 7 個を
    そのまま陽性対照にし（検出率 0.042、相関 0.009、K = 7）、不合格になった。
    **つまり定義マーカーの値は 4 実験を分離せず、陽性対照のセットの値が分離する。**
    定義マーカーの値は記録し、判定保留の材料にするが、単独では止めない。
    """
    pc = expr["goal_positive_control"]
    gd = expr["goal_defining"]
    corr = pc.get("corr_median", float("nan"))      # 陽性対照のまとまり（止める条件に使う）
    det = pc.get("detection_median", float("nan"))  # 陽性対照の検出率（止める条件に使う）
    dcorr = gd.get("corr_median", float("nan"))     # 定義マーカー（記録・判定保留のみ）
    ddet = gd.get("detection_median", float("nan"))
    ratio = cent["distance_ratio_min"]
    computable = "error" not in pc                  # K < 2 のときは相関・検出率を計算できない

    blocks, holds = [], []

    def add(lst, cond, name, msg):
        if cond:
            lst.append({"id": name, "reason": msg})

    # --- 不合格（1 つでも当たれば着手しない） ---
    add(blocks, K < th["block_K_min"], "K",
        f"陽性対照にできる遺伝子が {K} 個（下限 {th['block_K_min']} 個）。"
        f"K がこれより小さいと、目標の軸の健全性チェックの検出下限が、"
        f"本プロジェクトで観測されたどの陽性対照の効果よりも大きくなる"
        f"（K = 7 で 1.71 標準偏差。到達例 0 回）")
    add(blocks, computable and det < th["block_detection_median_min"], "detection",
        f"陽性対照の検出率の中央値が {det:.4f}"
        f"（下限 {th['block_detection_median_min']}）。測定法で拾えていない")
    add(blocks, ratio < th["block_axis_ratio_min"], "axis_length",
        f"軸の長さの比が {ratio:.4f}（下限 {th['block_axis_ratio_min']}）。"
        f"参照状態に比べて目標が近すぎる")
    add(blocks, computable and (corr < th["block_pair_corr_max"])
        and (det < th["block_pair_detection_max"]), "corr_and_detection",
        f"陽性対照のまとまり {corr:.4f}（< {th['block_pair_corr_max']}）と"
        f"検出率 {det:.4f}（< {th['block_pair_detection_max']}）が 2 つそろって低い。"
        f"目標は一貫した細胞状態ではなく「どれか 1 個がたまたま検出された細胞の寄せ集め」に近い")

    # --- 判定保留（弱点を明記した上でのみ着手） ---
    add(holds, computable and corr < th["hold_corr_median_min"], "corr",
        f"陽性対照のまとまりが {corr:.4f}（目安 {th['hold_corr_median_min']} 以上）")
    add(holds, computable and det < th["hold_detection_median_min"], "detection",
        f"陽性対照の検出率の中央値が {det:.4f}（目安 {th['hold_detection_median_min']} 以上）")
    add(holds, dcorr < th["hold_corr_median_min"] or ddet < th["hold_detection_median_min"],
        "goal_definition",
        f"目標を定義するマーカーのまとまりが {dcorr:.4f}・検出率の中央値が {ddet:.4f}"
        f"（目安 {th['hold_corr_median_min']} / {th['hold_detection_median_min']} 以上）。"
        f"**目標の状態そのものの根拠が弱い。** 陽性対照が働いても、"
        f"それは「目標のラベルが正しい」ことの証拠にはならない")
    add(holds, ratio < th["hold_axis_ratio_min"], "axis_length",
        f"軸の長さの比が {ratio:.4f}（目安 {th['hold_axis_ratio_min']} 以上）。"
        f"同じ大きさの摂動でも読み取れる信号がこの割合まで小さくなる")
    add(holds, pc.get("score_share_max", 0) >= th["hold_score_share_max"], "score_share",
        f"1 個のマーカー（{pc.get('score_share_top_gene')}）が得点の "
        f"{pc.get('score_share_max', float('nan')):.1%} を占める"
        f"（目安 {th['hold_score_share_max']:.0%} 未満）")
    add(holds, pc.get("detection_spread", 0) >= th["hold_detection_spread_max"], "detection_spread",
        f"陽性対照の検出率が {pc.get('detection_min', float('nan')):.4f}〜"
        f"{pc.get('detection_max', float('nan')):.4f} と "
        f"{pc.get('detection_spread', float('nan')):.1f} 倍ばらついている"
        f"（目安 {th['hold_detection_spread_max']:.0f} 倍未満）。"
        f"最も検出率の低い 1 個にセットの値が引きずられやすい")
    # 目標を遺伝子リストで定義していないときは None（問えない）。0 と同じに扱わない
    _ov = expr.get("circular_overlap")
    add(holds, _ov is not None and _ov > th["hold_circular_overlap_max"], "circular",
        f"陽性対照の {(_ov if _ov is not None else float('nan')):.0%} が目標の定義マーカーそのもの。"
        f"健全性チェックとしては必要だが（guardrails 2026-09-30）、"
        f"独立な陽性対照を別に置かないと循環になる")
    add(holds, bool(expr.get("below_cell_floor")), "cell_floor",
        f"摂動できる細胞数が下限 {expr.get('min_cells_per_gene')} を下回る遺伝子がある: "
        f"{expr.get('below_cell_floor')}")
    add(holds, pc.get("n_opposite_direction", 0) > 0, "opposite",
        f"目標の細胞で平均が負のマーカーが {pc.get('n_opposite_direction')} 個ある: "
        f"{pc.get('opposite_direction')}")

    verdict = "不合格" if blocks else ("判定保留" if holds else "合格")
    action = {
        "不合格": "**摂動の本番に着手しない。** 計画を作り直し、このチェックを通してから G2 に諮る。",
        "判定保留": "着手してよいが、**上の弱点と、それが結論に与える影響を plan.md に明記し、PI の判断を仰ぐ。**",
        "合格": "着手してよい。出力 JSON を事前登録のハッシュに含めること。",
    }[verdict]
    return {"verdict": verdict, "blocks": blocks, "holds": holds, "action": action,
            "K_positive_control": K}


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description="着手前チェック（摂動なしで軸が立つかを判定する）")
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", default=None, help="出力 JSON（既定は設定の output）")
    ap.add_argument("--no-exit-code", action="store_true", help="判定によらず 0 で終わる")
    a = ap.parse_args()

    cfg = load_cfg(rel(a.config))
    th = dict(DEFAULT_THRESHOLDS)
    th.update(cfg.get("thresholds", {}) or {})

    expr = expression_metrics(cfg)
    cent = centroid_metrics(cfg)
    mk = cfg["markers"]
    K = len(mk["goal_positive_control"]) if "goal_positive_control" in mk else len(mk["goal_defining"])
    res = judge(expr, cent, K, th)

    out = {
        "tool": "scripts/analysis/axis_feasibility.py",
        "generated_at": _dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "experiment_id": cfg.get("experiment_id"),
        "note": "摂動を 1 回も行っていない。状態の重心と発現データだけから計算した。",
        "thresholds": th,
        "centroid_metrics": cent,
        "expression_metrics": expr,
        "judgement": res,
    }
    dst = rel(a.out or cfg.get("output") or "outputs/axis_feasibility.json")
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")

    gd = expr["goal_defining"]
    pc = expr["goal_positive_control"]
    print(f"# 着手前チェック {cfg.get('experiment_id') or ''}")
    print(f"| 指標 | 値 | 目安 |")
    print(f"|---|---|---|")
    print(f"| 軸の長さの比（最小 / 3 回） | {cent['distance_ratio_min']:.4f} | "
          f"不合格 < {th['block_axis_ratio_min']} / 判定保留 < {th['hold_axis_ratio_min']} |")
    print(f"| 陽性対照の遺伝子数 K | {K} | 不合格 < {th['block_K_min']} |")
    print(f"| 陽性対照のまとまり（相関の中央値） | {pc.get('corr_median', float('nan')):+.4f} | "
          f"判定保留 < {th['hold_corr_median_min']} |")
    print(f"| 陽性対照の検出率の中央値 | {pc.get('detection_median', float('nan')):.4f} | "
          f"不合格 < {th['block_detection_median_min']} / 判定保留 < {th['hold_detection_median_min']} |")
    print(f"| （記録）目標を定義するマーカーのまとまり | {gd.get('corr_median', float('nan')):+.4f} | "
          f"判定保留 < {th['hold_corr_median_min']} |")
    print(f"| （記録）同じく検出率の中央値 | {gd.get('detection_median', float('nan')):.4f} | "
          f"判定保留 < {th['hold_detection_median_min']} |")
    print(f"| 1 個の得点寄与の最大 | {pc.get('score_share_max', float('nan')):.1%}"
          f"（{pc.get('score_share_top_gene')}） | 判定保留 ≥ {th['hold_score_share_max']:.0%} |")
    print(f"| 検出率のばらつき（最大/最小） | {pc.get('detection_spread', float('nan')):.1f} 倍 | "
          f"判定保留 ≥ {th['hold_detection_spread_max']:.0f} 倍 |")
    print()
    print(f"**判定: {res['verdict']}**")
    for b in res["blocks"]:
        print(f"- 不合格の理由（{b['id']}）: {b['reason']}")
    for h in res["holds"]:
        print(f"- 弱点（{h['id']}）: {h['reason']}")
    print(res["action"])
    print(f"\n出力: {dst}")

    if a.no_exit_code:
        return 0
    return {"合格": 0, "判定保留": 1, "不合格": 2}[res["verdict"]]


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SystemExit as e:
        if isinstance(e.code, int) or e.code is None:
            raise
        print(f"設定・入力のエラー: {e.code}", file=sys.stderr)   # 判定以外の停止は 3
        sys.exit(3)
    except Exception as e:  # noqa: BLE001
        print(f"実行エラー: {type(e).__name__}: {e}", file=sys.stderr)
        sys.exit(3)

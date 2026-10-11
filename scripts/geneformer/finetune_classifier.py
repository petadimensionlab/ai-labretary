#!/usr/bin/env python3
"""
Labretary: Geneformer CellClassifier の fine-tuning（config 駆動・デバイス非依存・再開可能）。EXP-20260926-003 で追加。

stages（--stage）:
  prepare  : [ローカル, CPU] config の donor 分割を検証し、学習・検証 donor の細胞だけを bundle ディレクトリに書く
             （train.dataset / eval.dataset / id_class.json / split_manifest.json / token dictionary / config / 本スクリプト /
             リモート用スクリプト / SHA256SUMS）。held-out（test）donor の細胞はローカル（bundle.local_test_dataset）にだけ書く
  train    : [リモート GPU または Mac] bundle から学習。CUDA / MPS / CPU（compute.device または GENEFORMER_DEVICE）。
             bf16 は対応を確認してから使う（非対応なら compute.bf16_unsupported に従う。既定は停止）。
             step checkpoint から自動再開。<out>/env.json, run_record.json, train_log.jsonl, training_summary.json,
             final/（最良モデル + SHA256SUMS）を書く。--probe N で N optimizer step だけ回し s/step とメモリを測る
  evaluate : [ローカル] held-out donor で FT 分類器と「凍結 Pretrained CLS 埋め込み + ロジスティック回帰」を比較し、
             事前登録したゲート（config evaluate.gate）で go/no-go を判定
  smoke    : [CPU] 小さなランダム BERT で prepare 済み bundle の一部を使い、train → 再開 → evaluate を通す（GPU 不要）

依存: torch, transformers, datasets, accelerate, scikit-learn, pyyaml と geneformer パッケージのソース
（collator_for_classification.py と token dictionary を直接読む。geneformer/__init__.py の重い依存は import しない）。
パスはプロジェクトルート相対。bundle 内のコピーから実行する場合（リモート）は --bundle / --model-dir / --out を渡す。
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.metadata as md
import importlib.util
import json
import math
import os
import pickle
import platform
import random
import shutil
import subprocess
import sys
import time
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
IN_BUNDLE = (HERE / "split_manifest.json").exists()
ROOT = HERE if IN_BUNDLE else HERE.parents[1]
BUNDLE_DATA = ("train.dataset", "eval.dataset")


# ----------------------------------------------------------------- utilities
def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def die(msg: str) -> None:
    raise SystemExit(f"ERROR: {msg}")


def rel(p) -> Path | None:
    if p is None:
        return None
    p = Path(p)
    return p if p.is_absolute() else ROOT / p


def sha256_file(p: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def write_sha256sums(d: Path, name: str = "SHA256SUMS") -> Path:
    """d 以下の全ファイル（SHA256SUMS 自身を除く）の sha256 を 'hash  relpath' 形式で書く。"""
    lines = [f"{sha256_file(p)}  {p.relative_to(d).as_posix()}"
             for p in sorted(d.rglob("*")) if p.is_file() and p.name != name]
    out = d / name
    out.write_text("\n".join(lines) + "\n")
    return out


def verify_sha256sums(d: Path, name: str = "SHA256SUMS") -> int:
    f = d / name
    if not f.exists():
        die(f"{f} がありません")
    bad, n = [], 0
    for line in f.read_text().splitlines():
        if not line.strip():
            continue
        h, rp = line.split("  ", 1)
        p = d / rp
        n += 1
        if not p.exists() or sha256_file(p) != h:
            bad.append(rp)
    if bad:
        die(f"{d}: sha256 不一致 / 欠落 {len(bad)} 件: {bad[:5]}")
    return n


def load_config(path: Path) -> dict:
    cfg = yaml.safe_load(path.read_text())
    cfg["_path"] = str(path)
    cfg["_sha256"] = sha256_file(path)
    return cfg


def apply_profile(cfg: dict, name: str | None) -> dict:
    if not name:
        return cfg
    prof = (cfg.get("profiles") or {}).get(name)
    if prof is None:
        die(f"profiles に {name} がありません（{list(cfg.get('profiles') or {})}）")
    cfg = copy.deepcopy(cfg)
    cfg["train"].update(prof)
    cfg["_profile"] = name
    return cfg


def geneformer_pkg_dir() -> Path:
    """geneformer パッケージのディレクトリ（__init__ を実行せずに探す）。GENEFORMER_PKG_DIR で上書き可。"""
    env = os.environ.get("GENEFORMER_PKG_DIR")
    if env:
        return Path(env)
    spec = importlib.util.find_spec("geneformer")
    if spec is None or not spec.submodule_search_locations:
        die("geneformer パッケージが見つかりません（../Geneformer/geneformer_hf を pip install -e するか GENEFORMER_PKG_DIR を指定）")
    return Path(list(spec.submodule_search_locations)[0])


def load_collator_cls():
    """geneformer/collator_for_classification.py を単体で読み込む（__init__ の ray / peft / scanpy などを避ける）。"""
    f = geneformer_pkg_dir() / "collator_for_classification.py"
    spec = importlib.util.spec_from_file_location("gf_collator_for_classification", f)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.DataCollatorForCellClassification


def load_token_dict(cfg: dict, bundle: Path | None = None) -> tuple[dict, Path]:
    name = cfg["model"].get("token_dictionary", "token_dictionary_gc104M.pkl")
    p = (bundle / name) if bundle and (bundle / name).exists() else geneformer_pkg_dir() / name
    with open(p, "rb") as f:
        return pickle.load(f), p


def pkg_versions() -> dict:
    out = {}
    for p in ["torch", "transformers", "datasets", "accelerate", "geneformer", "numpy", "scikit-learn",
              "safetensors", "pyyaml", "pyarrow"]:
        try:
            out[p] = md.version(p)
        except md.PackageNotFoundError:
            out[p] = None
    return out


def geneformer_digest() -> dict:
    """パッチ済み geneformer の同一性確認用: パッケージ内 .py の sha256 をまとめた digest と主要ファイルの sha256。"""
    d = geneformer_pkg_dir()
    h = hashlib.sha256()
    files = sorted(d.glob("*.py"))
    for f in files:
        h.update(f.name.encode())
        h.update(sha256_file(f).encode())
    key = {f.name: sha256_file(f) for f in files if f.name in ("collator_for_classification.py", "device.py",
                                                              "perturber_utils.py", "classifier.py")}
    return {"path": str(d), "py_digest": h.hexdigest(), "files": key}


def env_record(cfg: dict, device: str, extra: dict | None = None) -> dict:
    import torch
    rec = {
        "time": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "hostname": platform.node(), "platform": platform.platform(), "machine": platform.machine(),
        "python": sys.version.split()[0], "packages": pkg_versions(),
        "torch_cuda": torch.version.cuda, "cudnn": torch.backends.cudnn.version() if torch.backends.cudnn.is_available() else None,
        "mps_available": bool(getattr(torch.backends, "mps", None) and torch.backends.mps.is_available()),
        "device": device, "config": cfg["_path"], "config_sha256": cfg["_sha256"], "seed": cfg.get("seed"),
        "profile": cfg.get("_profile"), "geneformer": geneformer_digest(),
        "env_GENEFORMER_DEVICE": os.environ.get("GENEFORMER_DEVICE"),
        "env_CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES"),
    }
    if device == "cuda" and torch.cuda.is_available():
        props = torch.cuda.get_device_properties(0)
        rec["gpu"] = {"name": props.name, "total_memory_gb": round(props.total_memory / 2**30, 2),
                      "capability": f"{props.major}.{props.minor}", "count_visible": torch.cuda.device_count(),
                      "bf16_supported": torch.cuda.is_bf16_supported()}
        try:
            q = subprocess.run(["nvidia-smi", "--query-gpu=name,driver_version,memory.total", "--format=csv,noheader"],
                               capture_output=True, text=True, timeout=20)
            rec["nvidia_smi"] = q.stdout.strip().splitlines()
        except Exception as e:  # noqa: BLE001
            rec["nvidia_smi"] = f"unavailable: {e}"
    if extra:
        rec |= extra
    return rec


def check_versions(cfg: dict, strict: bool) -> list[str]:
    msgs = []
    have = pkg_versions()
    for p, want in (cfg.get("pinned_versions") or {}).items():
        got = have.get(p)
        if got is None or str(got).split("+")[0] != str(want):
            msgs.append(f"{p}: pinned {want} / installed {got}")
    for m in msgs:
        log(("ERROR " if strict else "WARN ") + "バージョン不一致 " + m)
    if msgs and strict:
        die("pinned_versions と一致しません（--strict-versions）")
    return msgs


def resolve_device(cfg: dict, section: str = "compute") -> str:
    import torch
    want = (os.environ.get("GENEFORMER_DEVICE") or str(cfg.get(section, {}).get("device", "auto"))).lower()
    avail = (["cuda"] if torch.cuda.is_available() else []) + \
            (["mps"] if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available() else []) + ["cpu"]
    dev = avail[0] if want == "auto" else want
    if dev not in avail:
        die(f"device={want} はこの環境で使えません（利用可能: {avail}）。黙って切り替えません")
    return dev


def resolve_dtype(cfg: dict, device: str, section: str = "compute") -> tuple[str, str | None]:
    """返り値: (実際に使う dtype 名, 切り替えた場合の理由)。bf16 非対応は compute.bf16_unsupported に従う。"""
    import torch
    want = str(cfg.get(section, {}).get("dtype", "bf16")).lower()
    if want != "bf16":
        return want, None
    if device == "cpu":
        return "fp32", "CPU のため fp32（bf16 は CUDA / MPS のみで使う）"
    if device == "cuda" and not torch.cuda.is_bf16_supported():
        fb = str(cfg.get("compute", {}).get("bf16_unsupported", "error")).lower()
        if fb == "error":
            die("この GPU は bf16 非対応です（compute.bf16_unsupported を fp16 か fp32 にするか、PI に確認）")
        return fb, f"bf16 非対応の GPU のため {fb} に変更（compute.bf16_unsupported）"
    return "bf16", None


def allow_numpy_in_checkpoints() -> None:
    """torch>=2.6 の torch.load(weights_only=True) は、transformers 4.46 が checkpoint に保存する rng_state.pth
    （numpy の乱数状態）を読めず、再開が UnpicklingError で止まる（EXP-20260926-003 の CPU smoke で確認）。
    自分で作った checkpoint に限って numpy の型だけを許可する（任意コード実行を許す weights_only=False にはしない）。"""
    import numpy as np
    import torch
    try:
        from numpy._core.multiarray import _reconstruct
    except ImportError:  # numpy<2
        from numpy.core.multiarray import _reconstruct
    safe = [_reconstruct, np.ndarray, np.dtype, type(np.dtype(np.uint32)), type(np.dtype(np.int64)),
            type(np.dtype(np.float64)), type(np.dtype(np.float32))]
    torch.serialization.add_safe_globals(safe)


def set_seed(seed: int) -> None:
    import numpy as np
    import torch
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


# ------------------------------------------------------------------- prepare
def expand_split(cfg: dict, donors_by_state: dict[str, set]) -> dict[str, dict[str, list[str]]]:
    """split セクションを {split: {state: [donor...]}} に展開（'all' はその状態の全 donor）。"""
    out = {}
    for sp in ("train", "eval", "test"):
        out[sp] = {}
        for state, ds in (cfg["split"].get(sp) or {}).items():
            out[sp][state] = sorted(donors_by_state.get(state, set())) if ds == "all" else [str(x) for x in ds]
    return out


def check_split(cfg: dict, split: dict, donors_by_state: dict[str, set]) -> list[str]:
    errs = []
    seen: dict[str, str] = {}
    for sp, st in split.items():
        for state, ds in st.items():
            if state not in cfg["data"]["label_map"]:
                errs.append(f"split.{sp}: 状態 {state} が data.label_map にありません")
            for d in ds:
                if d not in donors_by_state.get(state, set()):
                    errs.append(f"split.{sp}.{state}: donor {d} はこの状態の細胞を持ちません")
                if d in seen and seen[d] != sp:
                    errs.append(f"donor {d} が {seen[d]} と {sp} の両方にあります（漏洩）")
                seen[d] = sp
    # 対象状態の donor はすべて、どれか 1 つの分割に入っていること
    for state in cfg["data"]["label_map"]:
        miss = sorted(donors_by_state.get(state, set()) - {d for st in split.values() for d in st.get(state, [])})
        if miss:
            errs.append(f"状態 {state} の donor が分割に割り当てられていません: {miss}")
    start = cfg["data"]["isp_start_state"]
    if any(start in split[sp] for sp in ("train", "eval")):
        errs.append(f"ISP の start 状態 {start} が train/eval に入っています")
    if set(split["test"].get(start, [])) != donors_by_state.get(start, set()):
        errs.append(f"ISP の start 状態 {start} の全 donor を test（held-out）にしてください")
    for grp in cfg["split"].get("batch_groups") or []:
        sps = {seen.get(str(d)) for d in grp}
        if len(sps) > 1:
            errs.append(f"処理バッチ {grp} が複数の分割にまたがっています: {sps}")
    classes = cfg["data"]["classes"]
    for sp in ("train", "eval", "test"):
        labs = {cfg["data"]["label_map"][s] for s in split[sp]}
        if set(classes) - labs:
            errs.append(f"split.{sp} に含まれないクラス: {sorted(set(classes) - labs)}")
    return errs


def subsample_rows(rows: list[int], donors: list[str], cap: int | None, seed_key: str) -> list[int]:
    """donor ごとに cap 以下へ抽出（donor 名と seed から決まる乱数。順序に依存しない）。"""
    by: dict[str, list[int]] = {}
    for r in rows:
        by.setdefault(donors[r], []).append(r)
    keep = []
    for d in sorted(by):
        rr = sorted(by[d])
        if cap and len(rr) > cap:
            rr = sorted(random.Random(f"{seed_key}-{d}").sample(rr, cap))
        keep += rr
    return sorted(keep)


def stage_prepare(cfg: dict, args) -> None:
    import numpy as np
    from datasets import load_from_disk
    d, sp_cfg = cfg["data"], cfg["split"]
    tok = rel(d["tokenized_dataset"])
    ds = load_from_disk(str(tok), keep_in_memory=True)  # in-memory: 元の dataset ディレクトリに cache を書かない
    states = list(ds[d["state_key"]])
    donors = [str(x) for x in ds[d["donor_key"]]]
    lengths = np.asarray(list(ds["length"]))
    donors_by_state: dict[str, set] = {}
    for s, dn in zip(states, donors):
        donors_by_state.setdefault(s, set()).add(dn)
    split = expand_split(cfg, donors_by_state)
    errs = check_split(cfg, split, donors_by_state)
    if errs:
        for e in errs:
            log(f"ERROR {e}")
        die("split の検証に失敗しました")
    maxlen = cfg["train"].get("max_input_length", 4096)
    if lengths.max() > maxlen:
        die(f"入力長 {lengths.max()} > max_input_length {maxlen}")
    classes = d["classes"]
    cid = {c: i for i, c in enumerate(classes)}
    caps = sp_cfg.get("max_cells_per_donor") or {}

    rows_by_split, summary = {}, {}
    for sp, st in split.items():
        rows = []
        for state, dl in st.items():
            dset = set(dl)
            cand = [i for i, (s, dn) in enumerate(zip(states, donors)) if s == state and dn in dset]
            rows += subsample_rows(cand, donors, (caps.get(sp) or {}).get(state), f"{sp_cfg.get('seed', 0)}-{sp}-{state}")
        rows = sorted(rows)
        rows_by_split[sp] = rows
        per_donor: dict[str, dict] = {}
        for r in rows:
            k = per_donor.setdefault(donors[r], {"state": states[r], "n_cells": 0, "n_available": 0})
            k["n_cells"] += 1
        for dn, k in per_donor.items():
            k["n_available"] = int(sum(1 for s, x in zip(states, donors) if x == dn and s == k["state"]))
            k["median_length"] = float(np.median(lengths[[r for r in rows if donors[r] == dn]]))
        labs = [d["label_map"][states[r]] for r in rows]
        summary[sp] = {"n_cells": len(rows), "n_donors": len(per_donor),
                       "per_class": {c: labs.count(c) for c in classes},
                       "per_state": {s: sum(1 for r in rows if states[r] == s) for s in st},
                       "length_mean": round(float(lengths[rows].mean()), 1),
                       "length_median": float(np.median(lengths[rows])),
                       "tokens_total": int(lengths[rows].sum()), "per_donor": per_donor}

    test_donors = {dn for dl in split["test"].values() for dn in dl}
    bundle = rel(args.bundle or cfg["bundle"]["dir"])
    if bundle.exists():
        if not args.force:
            die(f"{bundle} は既にあります（作り直すなら --force）")
        shutil.rmtree(bundle)
    bundle.mkdir(parents=True)
    keep_cols = ["input_ids", "length", d["donor_key"], d["state_key"]]

    def build(rows):
        sub = ds.select(rows)
        sub = sub.remove_columns([c for c in sub.column_names if c not in keep_cols])
        sub = sub.add_column("cell_uid", [int(r) for r in rows])
        return sub.add_column("label", [cid[d["label_map"][states[r]]] for r in rows])

    for sp in ("train", "eval"):
        build(rows_by_split[sp]).save_to_disk(str(bundle / f"{sp}.dataset"))
    test_path = rel(cfg["bundle"]["local_test_dataset"])
    if test_path.exists():
        shutil.rmtree(test_path)
    test_path.parent.mkdir(parents=True, exist_ok=True)
    build(rows_by_split["test"]).save_to_disk(str(test_path))

    tdict, tpath = load_token_dict(cfg)
    shutil.copy2(tpath, bundle / tpath.name)
    (bundle / "id_class.json").write_text(json.dumps({str(i): c for i, c in enumerate(classes)}, indent=1) + "\n")
    for f in cfg["bundle"].get("copy_files") or []:
        src = rel(f)
        dst = bundle / ("remote" if "remote/" in f else "") / src.name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    shutil.copy2(cfg["_path"], bundle / "config-finetune.yaml")
    mw = rel(cfg["model"]["base_directory"]) / "model.safetensors"
    model_sha = sha256_file(mw) if mw.exists() else None  # リモートで同じ重みか照合する（run_remote.sh）
    if model_sha:
        (bundle / "remote").mkdir(exist_ok=True)
        (bundle / "remote" / "MODEL_SHA256").write_text(f"{model_sha}  Geneformer-V2-316M/model.safetensors\n")
    manifest = {
        "experiment_id": cfg["experiment_id"], "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "config_sha256": cfg["_sha256"], "source_tokenized": str(d["tokenized_dataset"]),
        "source_tokenized_data_sha256": sha256_file(next(tok.glob("data-*.arrow"))),
        "classes": {str(i): c for i, c in enumerate(classes)}, "label_map": d["label_map"],
        "splits": {sp: {k: v for k, v in summary[sp].items()} for sp in ("train", "eval")},
        "donors": {sp: split[sp] for sp in ("train", "eval")},
        "heldout": {"n_donors_excluded": len(test_donors), "n_cells_excluded_after_cap": summary["test"]["n_cells"],
                    "note": "held-out（test）donor の ID と細胞はローカル専用（bundle には含めない）"},
        "token_dictionary": {"file": tpath.name, "sha256": sha256_file(bundle / tpath.name)},
        "base_model": {"name": Path(cfg["model"]["base_directory"]).name, "model_safetensors_sha256": model_sha},
    }
    (bundle / "split_manifest.json").write_text(json.dumps(manifest, indent=1, ensure_ascii=False) + "\n")
    # 漏洩チェック（書いた bundle を読み直して確認）
    for sp in BUNDLE_DATA:
        b = load_from_disk(str(bundle / sp), keep_in_memory=True)
        bd = {str(x) for x in b[d["donor_key"]]}
        if bd & test_donors:
            die(f"bundle/{sp} に held-out donor が含まれています: {sorted(bd & test_donors)}")
        if d["isp_start_state"] in set(b[d["state_key"]]):
            die(f"bundle/{sp} に ISP の start 状態の細胞が含まれています")
        if set(b["cell_uid"]) & set(rows_by_split["test"]):
            die(f"bundle/{sp} と test で細胞が重複しています")
    n = write_sha256sums(bundle)
    full = {"split": split, "summary": summary, "bundle": str(bundle), "bundle_sha256sums_sha256": sha256_file(n),
            "test_dataset": str(test_path), "test_data_sha256": sha256_file(next(test_path.glob("data-*.arrow"))),
            "config_sha256": cfg["_sha256"]}
    lm = rel(cfg["bundle"]["local_manifest"])
    lm.write_text(json.dumps(full, indent=1, ensure_ascii=False) + "\n")
    for sp in ("train", "eval", "test"):
        s = summary[sp]
        log(f"{sp}: {s['n_cells']} 細胞 / {s['n_donors']} donor / クラス {s['per_class']} / 長さ中央値 {s['length_median']:.0f}"
            f" 平均 {s['length_mean']:.0f} / tokens {s['tokens_total']}")
    size = sum(p.stat().st_size for p in bundle.rglob("*") if p.is_file())
    log(f"bundle → {bundle}（{size / 2**20:.1f} MiB, SHA256SUMS {len(n.read_text().splitlines())} 件）")
    log(f"held-out test（ローカル専用）→ {test_path}、全分割の記録 → {lm}")


# --------------------------------------------------------------------- train
def resolve_bundle(cfg: dict, args) -> Path:
    if args.bundle:
        return Path(args.bundle).resolve()
    if IN_BUNDLE:
        return HERE
    return rel(cfg["bundle"]["dir"])


def load_bundle(bundle: Path, cfg: dict, verify: bool = True):
    from datasets import load_from_disk
    if verify:
        n = verify_sha256sums(bundle)
        log(f"bundle の sha256 を確認: {n} ファイル一致")
    man = json.loads((bundle / "split_manifest.json").read_text())
    tr = load_from_disk(str(bundle / "train.dataset"), keep_in_memory=True)
    ev = load_from_disk(str(bundle / "eval.dataset"), keep_in_memory=True)
    return tr, ev, man


def subsample_dataset(ds, per_donor: dict | None, cfg: dict, seed: int):
    if not per_donor:
        return ds
    dk, sk = cfg["data"]["donor_key"], cfg["data"]["state_key"]
    donors, states = [str(x) for x in ds[dk]], list(ds[sk])
    rows = []
    for state, cap in per_donor.items():
        cand = [i for i, s in enumerate(states) if s == state]
        rows += subsample_rows(cand, donors, cap, f"sub{seed}-{state}")
    return ds.select(sorted(rows))


def compute_metrics_fn(eval_pred):
    import numpy as np
    from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
    logits, labels = eval_pred
    logits = logits[0] if isinstance(logits, tuple) else logits
    p = np.exp(logits - logits.max(1, keepdims=True))
    p /= p.sum(1, keepdims=True)
    pred = p.argmax(1)
    out = {"accuracy": accuracy_score(labels, pred), "macro_f1": f1_score(labels, pred, average="macro")}
    if p.shape[1] == 2 and len(set(labels.tolist())) == 2:
        out["auroc"] = roc_auc_score(labels, p[:, 1])
    return out


def stage_train(cfg: dict, args) -> dict:
    import numpy as np
    import torch
    from transformers import (BertForSequenceClassification, EarlyStoppingCallback, Trainer, TrainerCallback,
                              TrainingArguments)
    from transformers.trainer_utils import get_last_checkpoint

    tcfg = cfg["train"]
    bundle = resolve_bundle(cfg, args)
    out = Path(args.out).resolve() if args.out else rel(tcfg["output_dir"])
    if args.probe:
        out = out.parent / (out.name + "-probe")
    model_dir = Path(args.model_dir or os.environ.get("GF_BASE_MODEL_DIR") or rel(cfg["model"]["base_directory"]))
    device = resolve_device(cfg)
    os.environ["GENEFORMER_DEVICE"] = device
    dtype, dtype_note = resolve_dtype(cfg, device)
    vmsgs = check_versions(cfg, args.strict_versions)
    if device == "cuda" and cfg["compute"].get("single_gpu", True) and torch.cuda.device_count() > 1:
        die(f"GPU が {torch.cuda.device_count()} 枚見えています。CUDA_VISIBLE_DEVICES で 1 枚に絞ってください")
    if device == "cuda" and cfg["compute"].get("allow_tf32", True):
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True
    if not (model_dir / "config.json").exists():
        die(f"事前学習モデルが見つかりません: {model_dir}（--model-dir または GF_BASE_MODEL_DIR）")

    tr, ev, man = load_bundle(bundle, cfg, verify=not args.skip_verify)
    if cfg["data"]["isp_start_state"] in set(tr[cfg["data"]["state_key"]]) | set(ev[cfg["data"]["state_key"]]):
        die("bundle に ISP の start 状態の細胞があります（held-out の漏洩）")
    seed = int(cfg.get("seed", 42))
    tr = subsample_dataset(tr, tcfg.get("subsample_train_per_donor"), cfg, seed)
    if args.limit_cells:  # smoke 用
        tr = tr.shuffle(seed=seed).select(range(min(args.limit_cells, len(tr))))
        ev = ev.shuffle(seed=seed).select(range(min(max(4, args.limit_cells // 2), len(ev))))
    maxlen = int(tcfg.get("max_input_length", 4096))
    if max(tr["length"]) > maxlen or max(ev["length"]) > maxlen:
        die(f"max_input_length {maxlen} を超える細胞があります")
    id_class = json.loads((bundle / "id_class.json").read_text())
    n_cls = len(id_class)
    tdict, tpath = load_token_dict(cfg, bundle)
    if sha256_file(tpath) != man["token_dictionary"]["sha256"]:
        die("token dictionary の sha256 が bundle の記録と一致しません")

    out.mkdir(parents=True, exist_ok=True)
    ckdir = out / "checkpoints"
    rec_path = out / "run_record.json"
    last = get_last_checkpoint(str(ckdir)) if ckdir.exists() else None
    if last and rec_path.exists():
        prev = json.loads(rec_path.read_text())
        if prev.get("config_sha256") != cfg["_sha256"] and not args.allow_config_change:
            die(f"再開しようとした checkpoint は別の config（sha {prev.get('config_sha256', '')[:12]}）で作られています"
                "（--allow-config-change で強制）")
        if prev.get("profile") != cfg.get("_profile"):
            log(f"WARN profile が変わりました: {prev.get('profile')} → {cfg.get('_profile')}（実効バッチは同じか確認）")

    eff = tcfg["per_device_train_batch_size"] * tcfg["gradient_accumulation_steps"]
    steps_per_epoch = math.ceil(len(tr) / tcfg["per_device_train_batch_size"] / tcfg["gradient_accumulation_steps"])
    max_steps = int(args.probe or args.max_steps or tcfg.get("max_steps", -1) or -1)
    env = env_record(cfg, device, {"dtype": dtype, "dtype_note": dtype_note, "version_mismatch": vmsgs,
                                   "bundle": str(bundle), "model_dir": str(model_dir), "out": str(out),
                                   "argv": sys.argv, "resume_from": last, "probe": bool(args.probe)})
    (out / ("env-probe.json" if args.probe else "env.json")).write_text(json.dumps(env, indent=1, ensure_ascii=False) + "\n")
    rec = {"config_sha256": cfg["_sha256"], "profile": cfg.get("_profile"), "bundle_manifest_config_sha256": man["config_sha256"],
           "n_train": len(tr), "n_eval": len(ev), "effective_batch": eff, "steps_per_epoch": steps_per_epoch,
           "train_args": tcfg, "device": device, "dtype": dtype, "seed": seed, "started": env["time"]}
    rec_path.write_text(json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    log(f"device {device} / dtype {dtype} / train {len(tr)} 細胞, eval {len(ev)} 細胞 / 実効バッチ {eff} / "
        f"{steps_per_epoch} step/epoch / resume: {last}")

    set_seed(seed)
    model = BertForSequenceClassification.from_pretrained(
        str(model_dir), num_labels=n_cls, id2label={int(k): v for k, v in id_class.items()},
        label2id={v: int(k) for k, v in id_class.items()}, output_hidden_states=False, output_attentions=False)
    nl = model.config.num_hidden_layers
    exp_nl = cfg["model"].get("expected_num_hidden_layers")
    if exp_nl and nl != exp_nl and not args.limit_cells:
        die(f"num_hidden_layers {nl} ≠ expected {exp_nl}")
    fl = int(tcfg.get("freeze_layers", 0))
    for layer in model.bert.encoder.layer[:fl]:
        for p in layer.parameters():
            p.requires_grad = False
    if tcfg.get("freeze_embeddings"):
        for p in model.bert.embeddings.parameters():
            p.requires_grad = False
    n_train_p = sum(p.numel() for p in model.parameters() if p.requires_grad)
    n_all_p = sum(p.numel() for p in model.parameters())
    log(f"凍結 {fl}/{nl} 層、学習パラメータ {n_train_p / 1e6:.1f}M / {n_all_p / 1e6:.1f}M")

    keep = ["input_ids", "length", "label"]
    tr_m = tr.remove_columns([c for c in tr.column_names if c not in keep])
    ev_m = ev.remove_columns([c for c in ev.column_names if c not in keep])
    Collator = load_collator_cls()
    collator = Collator(token_dictionary=tdict)

    eval_on = not args.probe
    targs = TrainingArguments(
        output_dir=str(ckdir), seed=seed, data_seed=seed,
        per_device_train_batch_size=tcfg["per_device_train_batch_size"],
        per_device_eval_batch_size=tcfg.get("per_device_eval_batch_size", 8),
        gradient_accumulation_steps=tcfg["gradient_accumulation_steps"],
        gradient_checkpointing=bool(tcfg.get("gradient_checkpointing", True)),
        gradient_checkpointing_kwargs={"use_reentrant": False} if tcfg.get("gradient_checkpointing", True) else None,
        learning_rate=float(tcfg["learning_rate"]), lr_scheduler_type=tcfg.get("lr_scheduler_type", "linear"),
        warmup_ratio=float(tcfg.get("warmup_ratio", 0.0)), weight_decay=float(tcfg.get("weight_decay", 0.0)),
        num_train_epochs=float(tcfg.get("num_train_epochs", 1)), max_steps=max_steps,
        eval_strategy="steps" if eval_on else "no", eval_steps=int(tcfg["eval_steps"]),
        save_strategy="steps" if eval_on else "no", save_steps=int(tcfg["eval_steps"]),
        save_total_limit=int(tcfg.get("save_total_limit", 3)),
        load_best_model_at_end=eval_on, metric_for_best_model=tcfg.get("metric_for_best_model", "eval_loss"),
        greater_is_better=bool(tcfg.get("greater_is_better", False)),
        logging_steps=int(tcfg.get("logging_steps", 5)), logging_first_step=True,
        group_by_length=bool(tcfg.get("group_by_length", True)), length_column_name="length",
        dataloader_num_workers=int(tcfg.get("dataloader_num_workers", 0)),
        bf16=(dtype == "bf16" and device != "cpu"), fp16=(dtype == "fp16"),
        use_cpu=(device == "cpu"), report_to="none", disable_tqdm=True,
    )

    class JsonlLogger(TrainerCallback):
        def __init__(self):
            self.t0 = time.time()
            self.t_last, self.s_last = time.time(), 0

        def on_log(self, a, state, control, logs=None, **kw):
            now = time.time()
            ds = state.global_step - self.s_last
            row = {"ts": time.strftime("%Y-%m-%dT%H:%M:%S"), "step": state.global_step, "epoch": state.epoch,
                   "elapsed_s": round(now - self.t0, 1), **(logs or {})}
            if ds > 0:
                row["s_per_opt_step"] = round((now - self.t_last) / ds, 3)
                self.t_last, self.s_last = now, state.global_step
            if device == "cuda":
                row["cuda_max_mem_gb"] = round(torch.cuda.max_memory_allocated() / 2**30, 2)
            elif device == "mps":
                row["mps_alloc_gb"] = round(torch.mps.driver_allocated_memory() / 2**30, 2)
            with open(out / ("train_log-probe.jsonl" if args.probe else "train_log.jsonl"), "a") as f:
                f.write(json.dumps(row, default=float) + "\n")

    callbacks = [JsonlLogger()]
    if eval_on and tcfg.get("early_stopping_patience"):
        callbacks.append(EarlyStoppingCallback(early_stopping_patience=int(tcfg["early_stopping_patience"])))
    trainer = Trainer(model=model, args=targs, data_collator=collator, train_dataset=tr_m,
                      eval_dataset=ev_m if eval_on else None, compute_metrics=compute_metrics_fn, callbacks=callbacks)
    if last:
        allow_numpy_in_checkpoints()
        log(f"checkpoint から再開: {last}")
    t0 = time.time()
    result = trainer.train(resume_from_checkpoint=last)
    wall = time.time() - t0
    st = trainer.state
    summ = {"finished": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "global_step": st.global_step, "epoch": st.epoch,
            "wall_seconds_this_session": round(wall, 1), "train_runtime": result.metrics.get("train_runtime"),
            "s_per_opt_step_session": round(wall / max(1, st.global_step - (int(Path(last).name.split("-")[-1]) if last else 0)), 3),
            "micro_batch": tcfg["per_device_train_batch_size"], "grad_accum": tcfg["gradient_accumulation_steps"],
            "effective_batch": eff, "steps_per_epoch": steps_per_epoch, "max_steps_cap": max_steps if max_steps > 0 else None,
            "device": device, "dtype": dtype, "dtype_note": dtype_note, "profile": cfg.get("_profile"),
            "n_train": len(tr), "n_eval": len(ev), "freeze_layers": fl, "trainable_params": n_train_p}
    if device == "cuda":
        summ["cuda_peak_mem_gb"] = round(torch.cuda.max_memory_allocated() / 2**30, 2)
        summ["gpu_total_mem_gb"] = env.get("gpu", {}).get("total_memory_gb")
    if args.probe:
        summ["probe"] = True
        (out / "probe_summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
        log(f"probe: {summ['s_per_opt_step_session']} s/optimizer step（実効バッチ {eff}）, peak "
            f"{summ.get('cuda_peak_mem_gb')} GB → 1 epoch ≈ {summ['s_per_opt_step_session'] * steps_per_epoch / 3600:.2f} h")
        return summ
    early = st.global_step < (max_steps if max_steps > 0 else steps_per_epoch * float(tcfg.get("num_train_epochs", 1))) - 1
    summ |= {"best_metric": st.best_metric, "best_model_checkpoint": st.best_model_checkpoint,
             "stopped_early_or_capped": bool(early),
             "truncated_by_max_steps": bool(max_steps > 0 and st.global_step >= max_steps)}
    final = out / "final"
    if final.exists():
        shutil.rmtree(final)
    trainer.save_model(str(final))  # load_best_model_at_end=True なので最良 checkpoint の重み
    shutil.copy2(bundle / "id_class.json", final / "id_class.json")
    shutil.copy2(bundle / "split_manifest.json", final / "split_manifest.json")
    ev_metrics = trainer.evaluate()
    summ["eval_metrics_best"] = ev_metrics
    (final / "training_summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    shutil.copy2(out / "env.json", final / "env.json")
    shutil.copy2(cfg["_path"], final / "config-finetune.yaml")
    write_sha256sums(final)
    (out / "training_summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    log(f"完了: best {summ['best_metric']} @ {summ['best_model_checkpoint']} / step {st.global_step} / "
        f"eval {json.dumps({k: round(v, 4) for k, v in ev_metrics.items() if isinstance(v, float)})}")
    log(f"final → {final}（SHA256SUMS 付き。ローカルへ戻して --stage evaluate）")
    return summ


# ------------------------------------------------------------------ evaluate
def forward_batches(model, ds, pad_id: int, device: str, dtype: str, bs: int, emb_layer: int | None):
    """返り値: (logits or None, cls_emb or None)。長さ順に並べて padding を減らし、元の順に戻す。"""
    import numpy as np
    import torch
    ids = list(ds["input_ids"])
    order = np.argsort([len(x) for x in ids])
    logits, embs = [None] * len(ids), [None] * len(ids)
    ac = torch.autocast(device_type=device, dtype=torch.bfloat16) if (dtype == "bf16" and device in ("cuda", "mps")) \
        else torch.autocast(device_type="cpu", enabled=False)
    model.eval()
    with torch.no_grad(), ac:
        for s in range(0, len(order), bs):
            idx = order[s:s + bs]
            L = max(len(ids[i]) for i in idx)
            x = torch.full((len(idx), L), pad_id, dtype=torch.long)
            m = torch.zeros((len(idx), L), dtype=torch.long)
            for j, i in enumerate(idx):
                x[j, :len(ids[i])] = torch.tensor(ids[i])
                m[j, :len(ids[i])] = 1
            o = model(input_ids=x.to(device), attention_mask=m.to(device), output_hidden_states=emb_layer is not None)
            if hasattr(o, "logits") and o.logits.shape[-1] < 100:
                lg = o.logits.float().cpu().numpy()
                for j, i in enumerate(idx):
                    logits[i] = lg[j]
            if emb_layer is not None:
                hs = o.hidden_states  # (embeddings, layer1, ..., layerN)。Geneformer: layer_to_quant = N + emb_layer
                n = len(hs) - 1
                e = hs[n + emb_layer][:, 0, :].float().cpu().numpy()
                for j, i in enumerate(idx):
                    embs[i] = e[j]
    L_ = np.vstack(logits) if logits[0] is not None else None
    E_ = np.vstack(embs) if embs[0] is not None else None
    return L_, E_


def donor_bootstrap(df, fn, n: int, seed: int = 0):
    import numpy as np
    rng = np.random.default_rng(seed)
    donors = df.donor.unique()
    vals = []
    for _ in range(n):
        pick = rng.choice(donors, len(donors), replace=True)
        sub = df.set_index("donor").loc[pick].reset_index()
        try:
            vals.append(fn(sub))
        except ValueError:
            continue
    return [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))] if vals else None


def stage_evaluate(cfg: dict, args) -> dict:
    import numpy as np
    import pandas as pd
    import torch
    from datasets import load_from_disk
    from scipy.stats import spearmanr
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import accuracy_score, f1_score, log_loss, roc_auc_score
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from transformers import BertForSequenceClassification, BertModel

    ecfg = cfg["evaluate"]
    ck = Path(args.checkpoint).resolve() if args.checkpoint else rel(ecfg["checkpoint"])
    base_dir = Path(args.baseline_dir).resolve() if args.baseline_dir else rel(ecfg["baseline"]["model_directory"])
    out = Path(args.eval_out).resolve() if args.eval_out else rel(ecfg["out_dir"])
    out.mkdir(parents=True, exist_ok=True)
    if not args.skip_verify:
        log(f"checkpoint の sha256 を確認: {verify_sha256sums(ck)} ファイル一致")
    device = resolve_device(cfg, "evaluate")
    dtype, _ = resolve_dtype(cfg, device, "evaluate")
    bundle = resolve_bundle(cfg, args)
    tr, _, man = load_bundle(bundle, cfg, verify=not args.skip_verify)
    test = load_from_disk(str(Path(args.test_dataset).resolve() if args.test_dataset
                              else rel(cfg["bundle"]["local_test_dataset"])), keep_in_memory=True)
    if args.limit_cells:
        tr = tr.shuffle(seed=0).select(range(min(args.limit_cells, len(tr))))
        test = test.shuffle(seed=0).select(range(min(args.limit_cells, len(test))))
    dk, sk = cfg["data"]["donor_key"], cfg["data"]["state_key"]
    train_donors = {str(x) for x in tr[dk]}
    if train_donors & {str(x) for x in test[dk]}:
        die("test と train で donor が重複しています")
    tdict, _ = load_token_dict(cfg, bundle)
    pad = tdict["<pad>"]
    id_class = json.loads((ck / "id_class.json").read_text())
    goal = cfg["data"]["classes"].index("PT_healthy") if "PT_healthy" in cfg["data"]["classes"] else 1
    bs = int(ecfg.get("batch_size", 8))
    env = env_record(cfg, device, {"stage": "evaluate", "checkpoint": str(ck), "dtype": dtype, "argv": sys.argv})
    (out / "env.json").write_text(json.dumps(env, indent=1, ensure_ascii=False) + "\n")
    tm = {}

    t = time.time()
    ft = BertForSequenceClassification.from_pretrained(str(ck)).to(device)
    lg, _ = forward_batches(ft, test, pad, device, dtype, bs, None)
    tm["ft_test_s"] = round(time.time() - t, 1)
    del ft
    p_ft = np.exp(lg - lg.max(1, keepdims=True))
    p_ft /= p_ft.sum(1, keepdims=True)

    t = time.time()
    base = BertModel.from_pretrained(str(base_dir), add_pooling_layer=False).to(device)
    _, e_tr = forward_batches(base, tr, pad, device, dtype, bs, int(ecfg["baseline"].get("emb_layer", -1)))
    _, e_te = forward_batches(base, test, pad, device, dtype, bs, int(ecfg["baseline"].get("emb_layer", -1)))
    tm["baseline_emb_s"] = round(time.time() - t, 1)
    del base
    if device == "mps":
        torch.mps.empty_cache()
    y_tr, y_te = np.asarray(tr["label"]), np.asarray(test["label"])
    clf = make_pipeline(StandardScaler(), LogisticRegression(C=float(ecfg["baseline"].get("C", 1.0)), max_iter=5000))
    clf.fit(e_tr.astype(np.float64), y_tr)
    p_base = clf.predict_proba(e_te.astype(np.float64))

    df = pd.DataFrame({"cell_uid": test["cell_uid"], "donor": [str(x) for x in test[dk]], "state": test[sk],
                       "length": test["length"], "label": y_te, "p_goal_ft": p_ft[:, goal], "p_goal_base": p_base[:, goal]})
    df.to_csv(out / "test_predictions.csv", index=False)

    def metrics(col: str, sub: pd.DataFrame) -> dict:
        p = sub[col].to_numpy()
        pred = (p >= 0.5).astype(int) if goal == 1 else (p < 0.5).astype(int)
        m = {"auroc": float(roc_auc_score(sub.label, p if goal == 1 else 1 - p)),
             "macro_f1": float(f1_score(sub.label, pred, average="macro")),
             "accuracy": float(accuracy_score(sub.label, pred)),
             "log_loss": float(log_loss(sub.label, np.column_stack([1 - p, p]) if goal == 1 else np.column_stack([p, 1 - p]),
                                        labels=[0, 1]))}
        sub = sub.assign(correct=(pred == sub.label.to_numpy()))
        dacc = sub.groupby("donor").agg(n=("correct", "size"), acc=("correct", "mean"), state=("state", "first"))
        big = dacc[dacc.n >= int(ecfg["gate"].get("min_donor_cells", 20))]
        m["worst_donor_accuracy"] = float(big.acc.min()) if len(big) else None
        m["worst_donor"] = str(big.acc.idxmin()) if len(big) else None
        m["per_donor_accuracy"] = {k: round(float(v), 4) for k, v in dacc.acc.items()}
        rho = []
        for lab in sorted(sub.label.unique()):
            s = sub[sub.label == lab]
            if len(s) > 10:
                rho.append(float(spearmanr(s[col], s.length).statistic))
        m["within_class_spearman_p_vs_length"] = rho
        return m

    res = {"n_test_cells": len(df), "n_test_donors": int(df.donor.nunique()),
           "n_per_state": df.state.value_counts().to_dict(), "timing_s": tm,
           "ft": metrics("p_goal_ft", df), "baseline": metrics("p_goal_base", df)}
    nb = int(ecfg.get("bootstrap_n", 1000)) if not args.limit_cells else 20
    for key, col in [("ft", "p_goal_ft"), ("baseline", "p_goal_base")]:
        res[key]["auroc_ci95_donor_bootstrap"] = donor_bootstrap(df, lambda s, c=col: roc_auc_score(s.label, s[c]), nb)
    res["diff_auroc_ci95_donor_bootstrap"] = donor_bootstrap(
        df, lambda s: roc_auc_score(s.label, s.p_goal_ft) - roc_auc_score(s.label, s.p_goal_base), nb)

    g = ecfg["gate"]
    checks = {}
    for mname in ("auroc", "macro_f1"):
        b, f = res["baseline"][mname], res["ft"][mname]
        margin = float(g.get("ceiling_margin", 0.01)) if b >= float(g.get("ceiling", 0.99)) else 0.0
        checks[f"{mname}_ft_ge_baseline"] = {"ft": f, "baseline": b, "margin": margin, "pass": bool(f >= b - margin)}
    w = res["ft"]["worst_donor_accuracy"]
    checks["worst_donor_accuracy"] = {"ft": w, "min": g.get("min_worst_donor_accuracy", 0.8),
                                      "pass": bool(w is not None and w >= float(g.get("min_worst_donor_accuracy", 0.8)))}
    res["gate"] = {"checks": checks, "go": all(c["pass"] for c in checks.values()),
                   "rule": "AUROC と macro-F1 がともにベースライン以上（ベースライン ≥ ceiling の指標は ceiling_margin の非劣性）"
                           "かつ held-out donor（≥ min_donor_cells 細胞）の最小正解率 ≥ min_worst_donor_accuracy"}
    res["smoke"] = bool(args.limit_cells)
    (out / "gate.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float) + "\n")
    log(f"FT: AUROC {res['ft']['auroc']:.4f} / macro-F1 {res['ft']['macro_f1']:.4f} / 最小 donor 正解率 {w}")
    log(f"baseline: AUROC {res['baseline']['auroc']:.4f} / macro-F1 {res['baseline']['macro_f1']:.4f}")
    log(f"ゲート: {'GO' if res['gate']['go'] else 'NO-GO'} → {out / 'gate.json'}")
    return res


# --------------------------------------------------------------------- smoke
def stage_smoke(cfg: dict, args) -> None:
    """CPU だけで train（途中停止 → 再開）と evaluate を通す。小さなランダム BERT を作業ディレクトリに作る。"""
    import torch
    from datasets import load_from_disk
    from transformers import BertConfig, BertForMaskedLM
    work = Path(args.smoke_dir or (ROOT / "experiments" / cfg["experiment_id"] / "work" / "smoke")).resolve()
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    src = resolve_bundle(cfg, args)
    verify_sha256sums(src)
    # 1) 小さな bundle（先頭 64 トークンに切り詰め）
    sb = work / "bundle"
    sb.mkdir()
    for name in BUNDLE_DATA:
        ds = load_from_disk(str(src / name), keep_in_memory=True).shuffle(seed=0).select(range(24))
        ds = ds.map(lambda e: {"input_ids": e["input_ids"][:64], "length": min(64, e["length"])}, keep_in_memory=True)
        ds.save_to_disk(str(sb / name))
    for f in ("id_class.json", "split_manifest.json", cfg["model"].get("token_dictionary", "token_dictionary_gc104M.pkl")):
        shutil.copy2(src / f, sb / f)
    write_sha256sums(sb)
    test_src = rel(cfg["bundle"]["local_test_dataset"])
    tds = load_from_disk(str(test_src), keep_in_memory=True).shuffle(seed=0).select(range(24))
    tds.map(lambda e: {"input_ids": e["input_ids"][:64], "length": min(64, e["length"])}, keep_in_memory=True) \
        .save_to_disk(str(work / "test.dataset"))
    # 2) 小さなランダム「事前学習」モデル
    tdict, _ = load_token_dict(cfg, src)
    bc = BertConfig(vocab_size=cfg["model"].get("expected_vocab_size", len(tdict)), hidden_size=32, num_hidden_layers=2,
                    num_attention_heads=2, intermediate_size=64, max_position_embeddings=4096, pad_token_id=tdict["<pad>"])
    torch.manual_seed(0)
    tiny = work / "tiny-base"
    BertForMaskedLM(bc).save_pretrained(str(tiny))
    write_sha256sums(tiny)
    c = copy.deepcopy(cfg)
    c["compute"] |= {"device": "cpu", "dtype": "fp32"}
    c["evaluate"] |= {"device": "cpu", "dtype": "fp32", "batch_size": 4}
    c["train"] |= {"freeze_layers": 1, "per_device_train_batch_size": 2, "gradient_accumulation_steps": 2,
                   "eval_steps": 2, "logging_steps": 1, "num_train_epochs": 1, "subsample_train_per_donor": None,
                   "early_stopping_patience": 50}
    os.environ["GENEFORMER_DEVICE"] = "cpu"
    ns = argparse.Namespace(bundle=str(sb), out=str(work / "ft"), model_dir=str(tiny), probe=0, max_steps=4,
                            strict_versions=False, skip_verify=False, allow_config_change=False, limit_cells=24)
    s1 = stage_train(c, ns)                      # 4 step で打ち切り（checkpoint-2, -4 を保存）
    ns.max_steps = 6
    s2 = stage_train(c, ns)                      # checkpoint-4 から再開して 6 step
    rr = json.loads((work / "ft" / "env.json").read_text())["resume_from"]
    assert rr and rr.endswith("checkpoint-4"), f"再開していません: {rr}"
    assert s2["global_step"] == 6, s2
    ns.probe = 3
    stage_train(c, ns)                           # probe の経路
    ne = argparse.Namespace(checkpoint=str(work / "ft" / "final"), baseline_dir=str(tiny), eval_out=str(work / "eval"),
                            test_dataset=str(work / "test.dataset"), bundle=str(sb), skip_verify=False, limit_cells=24)
    r = stage_evaluate(c, ne)
    assert "go" in r["gate"]
    log(f"smoke OK: train {s1['global_step']} → 再開 {s2['global_step']} step、evaluate（CPU・ランダム小モデル, 数値は無意味）"
        f" → {work}")


# ---------------------------------------------------------------------- main
def estimate(cfg: dict, n_train: int, n_eval: int, mean_len: float) -> dict:
    """MPS 実測（V2-316M, batch 2, 凍結 6 層, GC + bf16, 中央値 2,767 トークン ≈ 10 s/micro-step）からの目安。"""
    t = cfg["train"]
    micro = math.ceil(n_train / t["per_device_train_batch_size"]) * t["per_device_train_batch_size"] / 2  # batch-2 換算
    scale = mean_len / 2767
    ep = float(t.get("num_train_epochs", 1))
    return {"micro_steps_per_epoch_at_batch2": int(micro), "length_scale_vs_ref": round(scale, 2),
            "mps_hours_per_epoch_range": [round(micro * 10 * scale / 3600, 2), round(micro * 10 / 3600, 2)],
            "mps_hours_max_epochs_range": [round(micro * 10 * scale * ep / 3600, 2), round(micro * 10 * ep / 3600, 2)],
            "note": "下限 = 長さに比例と仮定、上限 = 10 s/step のまま。GPU は --probe で実測して補正"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--stage", required=True, choices=["prepare", "train", "evaluate", "smoke"])
    ap.add_argument("--profile", default=None, help="config の profiles から 1 つ（VRAM に合わせたバッチ構成）")
    ap.add_argument("--dry-run", action="store_true", help="検証と見積もりだけ（学習しない）")
    ap.add_argument("--bundle", default=None)
    ap.add_argument("--model-dir", default=None, help="事前学習モデル（train）")
    ap.add_argument("--out", default=None, help="学習の出力ディレクトリ（train）")
    ap.add_argument("--probe", type=int, default=0, help="N optimizer step だけ回して s/step とメモリを測る（train）")
    ap.add_argument("--max-steps", type=int, default=0)
    ap.add_argument("--strict-versions", action="store_true")
    ap.add_argument("--skip-verify", action="store_true", help="sha256 照合を省く（非推奨）")
    ap.add_argument("--allow-config-change", action="store_true")
    ap.add_argument("--force", action="store_true", help="prepare: 既存 bundle を作り直す")
    ap.add_argument("--checkpoint", default=None)
    ap.add_argument("--baseline-dir", default=None)
    ap.add_argument("--eval-out", default=None)
    ap.add_argument("--test-dataset", default=None)
    ap.add_argument("--smoke-dir", default=None)
    ap.add_argument("--limit-cells", type=int, default=0, help=argparse.SUPPRESS)
    args = ap.parse_args()
    cfg = apply_profile(load_config(Path(args.config).resolve()), args.profile)

    if args.stage == "prepare":
        if args.dry_run:
            from datasets import load_from_disk
            ds = load_from_disk(str(rel(cfg["data"]["tokenized_dataset"])), keep_in_memory=True)
            dbs: dict[str, set] = {}
            for s, dn in zip(ds[cfg["data"]["state_key"]], ds[cfg["data"]["donor_key"]]):
                dbs.setdefault(s, set()).add(str(dn))
            errs = check_split(cfg, expand_split(cfg, dbs), dbs)
            for e in errs:
                log(f"ERROR {e}")
            if errs:
                sys.exit(1)
            log("prepare dry-run OK（split の検証のみ）")
            return
        stage_prepare(cfg, args)
    elif args.stage == "train":
        if args.dry_run:
            bundle = resolve_bundle(cfg, args)
            tr, ev, man = load_bundle(bundle, cfg, verify=not args.skip_verify)
            tr = subsample_dataset(tr, cfg["train"].get("subsample_train_per_donor"), cfg, int(cfg.get("seed", 42)))
            import numpy as np
            est = estimate(cfg, len(tr), len(ev), float(np.mean(tr["length"])))
            dev = resolve_device(cfg)
            log(f"device {dev} / dtype {resolve_dtype(cfg, dev)[0]} / profile {cfg.get('_profile')} / train {len(tr)} / "
                f"eval {len(ev)} / 見積もり {json.dumps(est, ensure_ascii=False)}")
            check_versions(cfg, args.strict_versions)
            log("train dry-run OK")
            return
        stage_train(cfg, args)
    elif args.stage == "evaluate":
        stage_evaluate(cfg, args)
    elif args.stage == "smoke":
        stage_smoke(cfg, args)


if __name__ == "__main__":
    main()

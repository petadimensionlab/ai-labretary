#!/usr/bin/env python3
"""
Labretary: Geneformer in silico perturbation runner.

config.yaml を唯一のパラメータ源として、以下を実行する:
  1. tokenize   : h5ad -> Geneformer 用 .dataset（未作成の場合のみ）
  2. subsample  : seed ごとに細胞をサブサンプル
  3. state_embs : start/goal/alt 状態の埋め込み重心を抽出
  4. perturb    : InSilicoPerturber
  5. stats      : InSilicoPerturberStats
  6. summarize  : 全 run を results/summary.csv に統合（role = target/positive/negative）

モード:
  --dry-run  : config とファイルの存在・計算量を確認するだけ（GPU 不要）
  --smoke    : max_ncells=50, 遺伝子 3 個, min_cells_per_gene<=5 に縮小して全工程を通す
  --mock     : geneformer を使わず合成結果を作る（検証・図のパイプライン確認用）
  --approved : 夜間モードで大規模ジョブを許可する印（hook が参照）
  --precheck : 状態（start vs goal / alt）の分離 AUC だけを計算する（config の precheck セクション）
  --run NAME : config の runs から 1 つ選ぶ（states を上書き。結果は <results_dir>/<NAME>/）

任意の config セクション（EXP-20260926-001 で追加。なければ従来の挙動）:
  cell_selection : 遺伝子ごとに「その遺伝子が入力に含まれる start 細胞」から donor 上限つきで抽出して perturb
                   （Geneformer の max_ncells は発現細胞で絞った後に shuffle(seed=42) で切るだけで donor 上限がない）
  cell_selection.mode : "expressing"（既定・従来どおり）/ "common"
                   **common** は「その遺伝子が入力に載っている細胞」で絞らず、start 状態の全細胞から
                   seed ごとに 1 つだけ細胞集合を作り、**すべての遺伝子を同じ細胞で摂動する**。
                   過剰発現型（perturb_type: overexpress）でのみ意味がある。Geneformer 本体は未発現の
                   遺伝子の過剰発現を支援しており（perturber_utils.overexpress_tokens: 載っていなければ
                   rank-value encoding の先頭に挿入、in_silico_perturber.py:525 で overexpress は
                   細胞の絞り込みを行わない）、細胞を絞っていたのは run_isp.py 側の選択だった
                   （guardrails 2026-10-03、PI 指摘。EXP-20261004-001 で既定を変えずに選べるようにした）。
  state_embs     : 状態重心の細胞を状態ごと・donor ごとの上限つきで seed ごとに抽出
  controls.negative_sampler: {method: freq_matched, ...} 発現頻度マッチのランダム遺伝子
  controls.fixed_role      : negative_fixed の role 名（既定 negative）
  output.results_dir / output.experiment_dir, compute.cache_model, compute.sec_per_cell_gene
EXP-20260926-002 で追加（既定 off。なければ従来の挙動）:
  capture.embeddings: true → 遺伝子ごとに元/摂動後の CLS 埋め込み（fp32）を genes/<ENSG>/emb/{orig,pert}.npy、
                   状態重心を run-*/state_emb_vectors.npz に保存（cells.csv と同じ細胞順）。geneformer 本体は変更せず
                   pu.quant_cos_sims をラップして取り出す
  controls.negative_sampler: {method: explicit, genes: [...]} null 遺伝子を固定リストで与える（全 seed 共通）
  perturbation.controls_only: true → genes_to_perturb を空にして対照だけを回す（方法検証用）
  controls.negative_sampler.ref_genes: freq_matched の頻度の参照遺伝子（回さない・null からも除く。controls_only 用）
  環境変数 LABRETARY_DEVICE: compute.device を上書き（リモート GPU で同じ config を使う。device.json に記録）
  controls.extra_fixed: [{role: negctrl_panel, genes: [...]}, ...] negative_fixed とは別の role で回す固定遺伝子群
                   （陰性対照パネル・探索枠など。null からは除く）
EXP-20260926-003 で追加（既定 off。なければ従来の挙動）:
  state_embs.include_donors: {状態: [donor, ...]} → その状態の重心を指定 donor の細胞だけで作る（fine-tuning に使っていない
                   held-out donor で goal 重心を作るため）。書かない状態は従来どおり全 donor
  capture.classifier_logits: true（model.type: CellClassifier のみ）→ 分類ヘッドの logits を元/摂動後で取り出し
                   genes/<ENSG>/emb/logits_{orig,pert}.npy（fp32, cells.csv と同じ細胞順）、クラス名を run-*/classifier_labels.json に保存

Geneformer は版によって引数が異なるため、各クラスの signature を調べて
受け付けられる引数だけを渡す（渡せなかった引数は警告として表示）。
"""
from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import math
import os
import random
import sys
import time
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("pyyaml が必要です: pip install pyyaml")

ROOT = Path(__file__).resolve().parents[2]
VALID_PERTURB = {"delete", "overexpress", "inhibit", "activate"}
VALID_MODEL_TYPES = {"Pretrained", "CellClassifier", "GeneClassifier", "MTLCellClassifier"}
VALID_STATS = {"goal_state_shift", "vs_null", "mixture_model", "aggregate_data", "aggregate_gene_shifts"}
VALID_DEVICES = {"cuda", "mps", "cpu", "auto"}
VALID_DTYPES = {"bf16", "fp32"}


# ----------------------------------------------------------------- utilities
def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_config(path: Path) -> dict:
    with open(path) as f:
        cfg = yaml.safe_load(f)
    cfg["_path"] = str(path)
    cfg["_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return cfg


def rel(p: str | None) -> Path | None:
    if p is None:
        return None
    p = Path(p)
    return p if p.is_absolute() else ROOT / p


def exp_dir(cfg: dict) -> Path:
    """実験ディレクトリ（work/ と results*/ の親）。既定は config の置き場所。
    output.experiment_dir があればそれ（phase 0 用の派生 config を別フォルダに置く場合）。"""
    d = cfg.get("output", {}).get("experiment_dir")
    return rel(d) if d else Path(cfg["_path"]).resolve().parent


def call_compat(cls_or_fn, kwargs: dict, label: str):
    """signature にある引数だけ渡す。落とした引数は表示する。"""
    target = cls_or_fn.__init__ if inspect.isclass(cls_or_fn) else cls_or_fn
    params = inspect.signature(target).parameters
    accepts_var_kw = any(p.kind == p.VAR_KEYWORD for p in params.values())
    used = {k: v for k, v in kwargs.items() if accepts_var_kw or k in params}
    dropped = sorted(set(kwargs) - set(used))
    if dropped:
        log(f"WARN {label}: この版の geneformer は次の引数を受け付けません → {dropped}")
    return cls_or_fn(**used)


def set_seed(seed: int) -> None:
    random.seed(seed)
    try:
        import numpy as np
        np.random.seed(seed)
    except ImportError:
        pass
    try:
        import torch
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
        if torch.backends.mps.is_available():
            torch.mps.manual_seed(seed)
    except ImportError:
        pass


def available_devices() -> list[str]:
    try:
        import torch
    except ImportError:
        return ["cpu"]
    devs = ["cuda"] if torch.cuda.is_available() else []
    if torch.backends.mps.is_available():
        devs.append("mps")
    return devs + ["cpu"]


def resolve_device(cfg: dict) -> str:
    """compute.device を実デバイス名に解決する。auto は CUDA > MPS > CPU。
    指定デバイスが使えない場合は黙って切り替えず停止する（再現性のため）。"""
    # LABRETARY_DEVICE: config を変えずに実行先だけ変える（例: Mac 用 config をリモート CUDA で実行）。device.json に記録
    want = (os.environ.get("LABRETARY_DEVICE") or str(cfg.get("compute", {}).get("device", "auto"))).lower()
    avail = available_devices()
    dev = avail[0] if want == "auto" else want
    if dev not in avail:
        raise SystemExit(f"compute.device={want} はこの環境で使えません（利用可能: {avail}）。"
                         "config.yaml の compute.device を変更してください")
    return dev


def single_process_datasets(cfg: dict) -> None:
    """perturbation.nproc ≤ 1 のとき、datasets の map/filter に num_proc=None を渡して子プロセスを作らない。
    datasets 4.8 系は num_proc=1 でも fork したプール上で実行し、CUDA 初期化後の fork は
    "Cannot re-initialize CUDA in forked subprocess" で落ちる（EXP-20260926-003, spark-cf5b）。MPS でも同種の停止がある（guardrails）。"""
    if int(cfg.get("perturbation", {}).get("nproc", 4) or 1) > 1:
        return
    import datasets
    for name in ("map", "filter"):
        orig = getattr(datasets.Dataset, name)
        if getattr(orig, "_labretary_single", False):
            continue

        def wrapped(self, *a, _orig=orig, **kw):
            if kw.get("num_proc") is not None and kw["num_proc"] <= 1:
                kw["num_proc"] = None
            return _orig(self, *a, **kw)
        wrapped._labretary_single = True
        setattr(datasets.Dataset, name, wrapped)


def _numerics_record() -> dict:
    """結果の数値を左右する条件を、実行時の実測値として残す（2026-10-04）。

    **突き合わせずに「機体が違うから」と書かないための記録。** 2 つの結果を比べるときは、
    まずここを比べる。ここが違えば、埋め込みが違うのは当たり前であって発見ではない。
    """
    import torch
    import transformers
    rec = {
        "torch": torch.__version__,
        "transformers": transformers.__version__,
        "attn_implementation": sorted(set(_ATTN_IMPL)) or None,
        "cuda_matmul_allow_tf32": bool(torch.backends.cuda.matmul.allow_tf32),
        "cudnn_allow_tf32": bool(torch.backends.cudnn.allow_tf32),
        "float32_matmul_precision": torch.get_float32_matmul_precision(),
    }
    if torch.cuda.is_available():
        rec["cuda_device"] = torch.cuda.get_device_name(0)
        rec["cuda_version"] = torch.version.cuda
    return rec


def apply_dtype(cfg: dict) -> str:
    """compute.dtype=bf16 のとき、geneformer が読み込むモデルを bfloat16 にキャストする。
    ../Geneformer/analysis/_isp_common.py の GF_DTYPE と同じ方式（load_model をラップ）。
    bf16 は fp32 比 ~5x、埋め込み cosine 0.99995 以上、ISP 順位 Spearman 0.98
    （../Geneformer/docs/quantization/EXPERIMENT.md）。numpy 変換の修正（patches/bf16）が前提。"""
    name = str(cfg.get("compute", {}).get("dtype", "fp32")).lower()
    if name == "fp32":
        # **fp32 は「ものさしを作る」用途（CLAUDE.md §6.1）。機体をまたいで同じ値が出ることが条件。**
        # torch の既定に頼らず、精度を落としうる経路を明示的に閉じる（2026-10-04）。
        # これを閉じないと CUDA 側が tf32（仮数 10 bit）で積和し、Mac の真の fp32 と一致しない。
        import torch
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
        torch.set_float32_matmul_precision("highest")
        return "float32"
    import torch
    from geneformer import perturber_utils as pu
    original = pu.load_model

    def load_model_bf16(*args, **kwargs):
        return original(*args, **kwargs).to(torch.bfloat16)

    pu.load_model = load_model_bf16

    # shift は ≈1 のコサイン 2 つの差（~1e-4）で、bf16（有効桁 ~3 桁）では丸めに埋もれる。
    # forward は bf16 のまま、コサイン比較だけ fp32 に上げる（EXP-20260925-001 で確認）
    original_quant = pu.quant_cos_sims

    def quant_cos_sims_fp32(perturbation_emb, original_emb, cell_states_to_model, state_embs_dict, *a, **kw):
        if state_embs_dict is not None:
            state_embs_dict = {k: v.float() for k, v in state_embs_dict.items()}
        return original_quant(perturbation_emb.float(), original_emb.float(), cell_states_to_model,
                              state_embs_dict, *a, **kw)

    pu.quant_cos_sims = quant_cos_sims_fp32
    os.environ["GF_DTYPE"] = "bf16"
    return "bfloat16"


_CAPTURE: dict = {"on": False, "orig": [], "pert": [], "states": None,
                  "logits": False, "logits_on": False, "logits_orig": [], "logits_pert": [], "n_get": 0, "tag": "orig",
                  "id2label": None}


def apply_logit_capture(cfg: dict) -> bool:
    """capture.classifier_logits=true のとき、CellClassifier の分類ヘッド（model.classifier）の出力を取り出す。
    in_silico_perturber の CLS モードは 1 バッチごとに get_embs(元) → get_embs(摂動後) の順に呼ぶので、
    in_silico_perturber.get_embs を包んで呼び出しの偶奇で元/摂動後を区別する（geneformer 本体は変更しない）。
    状態重心（EmbExtractor）は emb_extractor.get_embs を直接呼び、logits_on=False の間なので記録されない。"""
    if not (cfg.get("capture") or {}).get("classifier_logits", False):
        return False
    from geneformer import in_silico_perturber as ispm
    from geneformer import perturber_utils as pu
    inner_get, inner_load = ispm.get_embs, pu.load_model

    def get_embs_tagged(*a, **kw):
        _CAPTURE["n_get"] += 1
        _CAPTURE["tag"] = "orig" if _CAPTURE["n_get"] % 2 == 1 else "pert"
        return inner_get(*a, **kw)

    def hook(_mod, _inp, out):
        if _CAPTURE["logits_on"]:
            _CAPTURE["logits_" + _CAPTURE["tag"]].append(out.detach().float().cpu().numpy())

    def load_model_hooked(*a, **kw):
        m = inner_load(*a, **kw)
        if hasattr(m, "classifier") and not getattr(m, "_labretary_logit_hook", False):
            m.classifier.register_forward_hook(hook)
            m._labretary_logit_hook = True
            _CAPTURE["id2label"] = {str(k): v for k, v in (m.config.id2label or {}).items()}
        return m

    ispm.get_embs = get_embs_tagged
    pu.load_model = load_model_hooked
    _CAPTURE["logits"] = True
    return True


def apply_capture(cfg: dict) -> bool:
    """capture.embeddings=true のとき、CLS モードの pu.quant_cos_sims 呼び出し（emb_mode="cell", 2 次元）から
    元・摂動後の埋め込みと状態重心を fp32 で取り出す。apply_dtype の後に呼ぶ。既定 off（従来と同じ）。
    in_silico_perturber.isp_perturb_set_special は pu.quant_cos_sims を属性参照で呼ぶので、ここでの差し替えが効く。"""
    if not (cfg.get("capture") or {}).get("embeddings", False):
        return False
    from geneformer import perturber_utils as pu
    inner = pu.quant_cos_sims

    def quant_cos_sims_capture(perturbation_emb, original_emb, cell_states_to_model, state_embs_dict, *a, **kw):
        if kw.get("emb_mode", a[0] if a else "gene") == "cell" and perturbation_emb.dim() == 2:
            _CAPTURE["pert"].append(perturbation_emb.detach().float().cpu().numpy())
            _CAPTURE["orig"].append(original_emb.detach().float().cpu().numpy())
            if state_embs_dict is not None and _CAPTURE["states"] is None:
                _CAPTURE["states"] = {str(k): v.detach().float().cpu().numpy() for k, v in state_embs_dict.items()}
        return inner(perturbation_emb, original_emb, cell_states_to_model, state_embs_dict, *a, **kw)

    pu.quant_cos_sims = quant_cos_sims_capture
    _CAPTURE["on"] = True
    return True


def capture_reset() -> None:
    _CAPTURE.update(orig=[], pert=[], states=None, logits_orig=[], logits_pert=[], n_get=0, tag="orig",
                    logits_on=_CAPTURE["logits"])


def capture_save(gd: Path, run_dir: Path, n_expected: int | None) -> None:
    """取り出した埋め込みを保存する。細胞数が cells.csv と合わなければ停止（対応付けを誤らない）。"""
    import numpy as np
    if _CAPTURE["logits"]:
        _CAPTURE["logits_on"] = False
        if not _CAPTURE["logits_pert"]:
            raise SystemExit(f"{gd}: capture.classifier_logits=true だが logits が取り出せませんでした（CellClassifier か確認）")
        lo, lp = np.concatenate(_CAPTURE["logits_orig"]), np.concatenate(_CAPTURE["logits_pert"])
        if len(lo) != len(lp) or (n_expected is not None and len(lo) != n_expected):
            raise SystemExit(f"{gd}: logits の細胞数が合いません（元 {len(lo)} / 摂動後 {len(lp)} / cells.csv {n_expected}）")
        (gd / "emb").mkdir(parents=True, exist_ok=True)
        np.save(gd / "emb" / "logits_orig.npy", lo.astype(np.float32))
        np.save(gd / "emb" / "logits_pert.npy", lp.astype(np.float32))
        lab = run_dir / "classifier_labels.json"
        if not lab.exists():
            lab.write_text(json.dumps(_CAPTURE["id2label"], ensure_ascii=False) + "\n")
    if not _CAPTURE["on"]:
        return
    if not _CAPTURE["pert"]:
        raise SystemExit(f"{gd}: capture.embeddings=true だが埋め込みが取り出せませんでした（emb_mode が cls か確認）")
    orig, pert = np.concatenate(_CAPTURE["orig"]), np.concatenate(_CAPTURE["pert"])
    if n_expected is not None and len(orig) != n_expected:
        raise SystemExit(f"{gd}: 取り出した埋め込み {len(orig)} 細胞 ≠ cells.csv {n_expected} 細胞")
    (gd / "emb").mkdir(parents=True, exist_ok=True)
    np.save(gd / "emb" / "orig.npy", orig.astype(np.float32))
    np.save(gd / "emb" / "pert.npy", pert.astype(np.float32))
    sv = run_dir / "state_emb_vectors.npz"
    if _CAPTURE["states"] is not None and not sv.exists():
        np.savez(sv, **_CAPTURE["states"])


# ------------------------------------------------------------------ validate
def validate(cfg: dict) -> list[str]:
    errs, warns = [], []
    for sec in ("model", "data", "states", "perturbation", "stats", "robustness"):
        if sec not in cfg:
            errs.append(f"セクション '{sec}' がありません")
    if errs:
        return errs
    m, d, s, p = cfg["model"], cfg["data"], cfg["states"], cfg["perturbation"]
    if m.get("type") not in VALID_MODEL_TYPES:
        errs.append(f"model.type は {VALID_MODEL_TYPES} のいずれか")
    if m.get("type") == "CellClassifier" and not m.get("num_classes"):
        errs.append("CellClassifier には model.num_classes が必要")
    if p.get("perturb_type") not in VALID_PERTURB:
        errs.append(f"perturbation.perturb_type は {VALID_PERTURB}")
    if p.get("perturb_type") in {"inhibit", "activate"} and p.get("perturb_rank_shift") is None:
        errs.append("inhibit/activate には perturb_rank_shift が必要")
    if cfg["stats"].get("mode") not in VALID_STATS:
        errs.append(f"stats.mode は {VALID_STATS}")
    if s.get("start_state") == s.get("goal_state"):
        errs.append("start_state と goal_state が同じです")
    g = p.get("genes_to_perturb")
    if p.get("controls_only") and g == []:
        pass  # 方法検証用: 対照だけを回す
    elif g != "all" and not (isinstance(g, list) and g):
        errs.append("genes_to_perturb は 'all' か Ensembl ID のリスト（対照だけなら perturbation.controls_only: true と []）")
    if isinstance(g, list):
        bad = [x for x in g if not str(x).startswith("ENS")]
        if bad:
            errs.append(f"Ensembl ID でない遺伝子があります: {bad[:5]}")
    dev = str(cfg.get("compute", {}).get("device", "auto")).lower()
    if dev not in VALID_DEVICES:
        errs.append(f"compute.device は {VALID_DEVICES}")
    elif dev != "auto" and dev not in available_devices():
        warns.append(f"compute.device={dev} はこの環境で使えません（利用可能: {available_devices()}）")
    dt = str(cfg.get("compute", {}).get("dtype", "fp32")).lower()
    if dt not in VALID_DTYPES:
        errs.append(f"compute.dtype は {VALID_DTYPES}")
    if dt == "bf16" and str(cfg.get("compute", {}).get("device", "auto")).lower() == "cpu":
        warns.append("bf16 は CPU では遅くなることが多い")
    sel_mode = str((cfg.get("cell_selection") or {}).get("mode", "expressing")).lower()
    if sel_mode not in {"expressing", "common"}:
        errs.append("cell_selection.mode は 'expressing'（既定）か 'common'")
    elif sel_mode == "common" and p.get("perturb_type") != "overexpress":
        errs.append("cell_selection.mode: common は perturb_type: overexpress のときだけ使えます")
    dg = (cfg.get("cell_selection") or {}).get("donor_groups")
    if dg:
        if sel_mode != "common":
            errs.append("cell_selection.donor_groups は mode: common のときだけ使えます")
        seen: dict[str, str] = {}
        for k, v in dg.items():
            if not (isinstance(v, list) and v):
                errs.append(f"cell_selection.donor_groups.{k} は提供者のリスト")
                continue
            for d_ in v:
                if str(d_) in seen:
                    errs.append(f"cell_selection.donor_groups: 提供者 {d_} が群 {seen[str(d_)]} と {k} に重複")
                seen[str(d_)] = str(k)
        miss = [x for x in cfg["robustness"].get("seeds", []) if str(x) not in {str(k) for k in dg}]
        if miss:
            errs.append(f"cell_selection.donor_groups に seed {miss} の群がありません")
    emb_mode = p.get("emb_mode", "cls")
    if m.get("version") == "V2" and "cls" not in emb_mode:
        errs.append("V2 は入力先頭が <cls> のため perturbation.emb_mode は cls か cls_and_gene")
    if not cfg.get("controls", {}).get("positive"):
        warns.append("ポジティブコントロールが空です（validator で Fail になります）")
    smp = cfg.get("controls", {}).get("negative_sampler") or {}
    if p.get("controls_only") and g == [] and smp.get("method") == "freq_matched" and not smp.get("ref_genes"):
        errs.append("controls_only で freq_matched を使うときは negative_sampler.ref_genes が必要")
    if smp.get("method") == "explicit" and not smp.get("genes"):
        errs.append("negative_sampler.method=explicit には genes（Ensembl ID のリスト）が必要")
    if (cfg.get("capture") or {}).get("embeddings") and "cls" not in emb_mode:
        errs.append("capture.embeddings は emb_mode: cls のときだけ使えます")
    if (cfg.get("capture") or {}).get("classifier_logits"):
        if m.get("type") != "CellClassifier":
            errs.append("capture.classifier_logits は model.type: CellClassifier のときだけ使えます")
        if emb_mode != "cls":
            errs.append("capture.classifier_logits は emb_mode: cls のときだけ使えます")
    inc = (cfg.get("state_embs") or {}).get("include_donors") or {}
    known = {s.get("start_state"), s.get("goal_state"), *s.get("alt_states", [])}
    for st, dl in inc.items():
        if st not in known:
            errs.append(f"state_embs.include_donors: {st} は start/goal/alt のどれでもありません")
        if not (isinstance(dl, list) and dl):
            errs.append(f"state_embs.include_donors.{st} は donor のリスト")
    if len(cfg["robustness"].get("seeds", [])) < 3:
        warns.append("seeds が 3 未満です（CLAUDE.md §6）")
    md = rel(m.get("directory"))
    if not (md and md.exists()):
        warns.append(f"モデルディレクトリが見つかりません: {md}")
    elif m.get("type") == "CellClassifier" and (md / "config.json").exists():
        n_lab = len(json.loads((md / "config.json").read_text()).get("id2label") or {})
        if n_lab and n_lab != m.get("num_classes"):
            errs.append(f"model.num_classes={m.get('num_classes')} ≠ {md}/config.json の クラス数 {n_lab}")
    tok, h5 = rel(d.get("tokenized_dataset")), rel(d.get("input_h5ad"))
    if not ((tok and tok.exists()) or (h5 and h5.exists())):
        warns.append(f"入力データが見つかりません: {tok} / {h5}")
    for w in warns:
        log(f"WARN {w}")
    return errs


def estimate(cfg: dict) -> dict:
    """ざっくりした forward 回数と GPU 時間の見積もり（A100 相当、実測で要補正。MPS/CPU ではさらに遅い）。"""
    sel = cfg.get("cell_selection") or {}
    n_cells = sel.get("max_cells_per_gene") or cfg["data"].get("max_ncells") or 10000
    g = cfg["perturbation"]["genes_to_perturb"]
    n_seeds = len(cfg["robustness"].get("seeds", [0]))
    n_layers = len(cfg["robustness"].get("emb_layers", [-1]))
    n_runs = len(cfg.get("runs") or [None])
    input_len = 4096 if cfg["model"].get("version") == "V2" else 2048
    if g == "all":
        n_genes = None
        forwards = n_cells * min(input_len, 2000)  # 1 細胞あたり発現遺伝子数 ~ 2000 と仮定
    else:
        c = cfg.get("controls", {})
        n_genes = len(g) + len(c.get("positive", [])) + len(c.get("negative_fixed", [])) + c.get("negative_random_n", 0) \
            + sum(len(x.get("genes", [])) for x in c.get("extra_fixed", []) or [])
        forwards = n_cells * n_genes
    forwards *= n_seeds * n_layers
    per_sec = 400 if input_len == 2048 else 150
    hours = forwards / per_sec / 3600
    out = {"forward_passes": int(forwards), "gpu_hours_estimate": round(hours, 2)}
    # 実測値（compute.sec_per_cell_gene: 1 細胞 × 1 遺伝子の秒数, モデルロード除く）があればそれで見積もる
    spc = cfg.get("compute", {}).get("sec_per_cell_gene")
    if spc and n_genes:
        se = cfg.get("state_embs") or {}
        n_states = 2 + len(cfg["states"].get("alt_states", []))
        emb_cells = (se.get("max_cells_per_state") or cfg["data"].get("max_ncells") or 500) * n_states
        per_run = (forwards * spc + n_seeds * n_layers * emb_cells * spc / 2) / 3600  # 重心は 1 forward/細胞
        out |= {"measured_hours_per_run": round(per_run, 2), "n_runs": n_runs,
                "measured_hours_total": round(per_run * n_runs, 2)}
        hours = per_run * n_runs
    out["needs_G2_approval"] = hours > 24 or n_cells > 10000 or g == "all"  # GPU 24 h（PI 2026-09-27）
    return out


# ------------------------------------------------------------------ pipeline
def tokenize(cfg: dict, work: Path) -> Path:
    tok = rel(cfg["data"]["tokenized_dataset"])
    if tok.exists():
        log(f"tokenized dataset を再利用: {tok}")
        return tok
    from geneformer import TranscriptomeTokenizer
    h5 = rel(cfg["data"]["input_h5ad"])
    src = work / "h5ad_in"
    src.mkdir(parents=True, exist_ok=True)
    link = src / h5.name
    if not link.exists():
        link.symlink_to(h5)  # data/raw は触らない
    tk = call_compat(TranscriptomeTokenizer, {
        "custom_attr_name_dict": cfg["data"].get("custom_attr_name_dict"),
        "nproc": cfg["perturbation"].get("nproc", 4),
        "model_version": cfg["model"].get("version"),
    }, "TranscriptomeTokenizer")
    tk.tokenize_data(str(src), str(tok.parent), tok.stem, file_format="h5ad")
    log(f"tokenized → {tok}")
    return tok


def subsample(tok: Path, seed: int, frac: float, work: Path) -> Path:
    from datasets import load_from_disk
    out = work / f"subsample-seed{seed}.dataset"
    if out.exists():
        return out
    ds = load_from_disk(str(tok)).shuffle(seed=seed)
    ds = ds.select(range(int(len(ds) * frac)))
    ds.save_to_disk(str(out))
    return out


def sample_negative_controls(tok: Path, n: int, exclude: set[str], seed: int) -> list[str]:
    """発現頻度でビン分けし、ターゲットと同程度の頻度帯からランダム遺伝子を選ぶ（簡易版）。"""
    if n <= 0:
        return []
    from collections import Counter
    from datasets import load_from_disk
    import pickle
    try:
        from geneformer import TOKEN_DICTIONARY_FILE
    except ImportError:
        log("WARN TOKEN_DICTIONARY_FILE が見つからないためネガティブコントロールを自動選択できません")
        return []
    with open(TOKEN_DICTIONARY_FILE, "rb") as f:
        tok2id = pickle.load(f)
    id2ens = {v: k for k, v in tok2id.items() if str(k).startswith("ENS")}
    ds = load_from_disk(str(tok))
    cnt = Counter()
    for ids in ds.select(range(min(len(ds), 2000)))["input_ids"]:
        cnt.update(ids)
    genes = [id2ens[i] for i, _ in cnt.most_common() if i in id2ens and id2ens[i] not in exclude]
    mid = genes[len(genes) // 10: len(genes) // 2]  # 極端な高発現・低発現を避ける
    rng = random.Random(seed)
    return rng.sample(mid, min(n, len(mid)))


def gene_token_dict() -> dict[str, int]:
    import pickle
    from geneformer import TOKEN_DICTIONARY_FILE
    with open(TOKEN_DICTIONARY_FILE, "rb") as f:
        return pickle.load(f)


_START_CACHE: dict = {}


def start_matrix(tok: Path, cfg: dict):
    """start 状態の細胞 × トークンの有無（CSC 疎行列）。dataset 上の行番号と donor も返す（キャッシュ）。"""
    start, key = cfg["states"]["start_state"], cfg["data"]["state_key"]
    ck = (str(tok), start)
    if ck in _START_CACHE:
        return _START_CACHE[ck]
    import numpy as np
    import scipy.sparse as sp
    from datasets import load_from_disk
    dkey = (cfg.get("cell_selection") or {}).get("donor_key", "donor_id")
    ds = load_from_disk(str(tok))
    rows = np.array([i for i, s in enumerate(ds[key]) if s == start])
    if len(rows) == 0:
        raise SystemExit(f"start_state={start} の細胞が {tok} にありません")
    ids = ds.select(rows)["input_ids"]
    indptr = np.cumsum([0] + [len(x) for x in ids])
    col = np.concatenate([np.asarray(x) for x in ids])
    m = sp.csr_matrix((np.ones(len(col), dtype=bool), col, indptr),
                      shape=(len(rows), max(int(col.max()) + 1, len(gene_token_dict())))).tocsc()
    donors = [str(x) for x in ds.select(rows)[dkey]] if dkey in ds.column_names else ["NA"] * len(rows)
    _START_CACHE[ck] = (m, rows, donors)
    return _START_CACHE[ck]


def start_detection(tok: Path, cfg: dict, work: Path) -> dict[str, int]:
    """start 状態の細胞のうち、各遺伝子が Geneformer 入力（rank 符号化後の input_ids）に
    含まれる細胞数。delete が効くのはこの細胞だけなので、発現頻度の定義にもこれを使う。"""
    import numpy as np
    start = cfg["states"]["start_state"]
    cache = work / f"detection-{start}.json"
    if cache.exists():
        return json.loads(cache.read_text())
    m, _, _ = start_matrix(tok, cfg)
    cnt = np.asarray(m.sum(axis=0)).ravel()
    det = {k: int(cnt[v]) for k, v in gene_token_dict().items()
           if str(k).startswith("ENS") and v < len(cnt) and cnt[v] > 0}
    cache.write_text(json.dumps(det) + "\n")
    return det


def sample_negative_matched(det: dict[str, int], n: int, refs: list[str], exclude: set[str],
                            seed: int, opts: dict) -> list[str]:
    """発現頻度マッチのランダム遺伝子（HYP-20260925-01 の設計）。
    候補 = start 細胞のうち min_cells 以上で発現している遺伝子（exclude を除く）。
    参照遺伝子（既定は対象遺伝子）を順番に回し、発現細胞数の比が ratio_tolerance 以内の候補から
    重複なしで 1 つずつ抽出する（該当がなければ log 比が最も近い候補）。seed ごとに独立に抽出。"""
    if n <= 0:
        return []
    rng = random.Random(seed)
    tol = math.log(opts.get("ratio_tolerance", 1.5))
    pool = sorted(g for g, c in det.items() if c >= opts.get("min_cells", 300) and g not in exclude)
    refs = [r for r in refs if det.get(r, 0) > 0]
    if not refs:
        raise SystemExit("negative_sampler: 参照遺伝子が start 細胞で 1 つも検出されていません")
    picked: list[str] = []
    for i in range(n):
        lr = math.log(det[refs[i % len(refs)]])
        cand = [g for g in pool if g not in picked and abs(math.log(det[g]) - lr) <= tol]
        if not cand:
            rest = [g for g in pool if g not in picked]
            cand = [min(rest, key=lambda g: abs(math.log(det[g]) - lr))]
        picked.append(rng.choice(cand))
    return picked


def balanced_subset(idx: list[int], donors: list[str], n_max: int | None, per_donor: int | None,
                    rng: random.Random) -> list[int]:
    """idx をシャッフルし、donor あたり per_donor 以下・合計 n_max 以下で先頭から採用する。"""
    order = list(idx)
    rng.shuffle(order)
    used: dict[str, int] = {}
    keep = []
    for i in order:
        d = donors[i]
        if per_donor and used.get(d, 0) >= per_donor:
            continue
        used[d] = used.get(d, 0) + 1
        keep.append(i)
        if n_max and len(keep) >= n_max:
            break
    return sorted(keep)


def rows_file_path(cfg: dict, seed: int) -> Path | None:
    """state_embs.rows_file（{seed} を含むパス）を解決する。指定がなければ None。"""
    rf = (cfg.get("state_embs") or {}).get("rows_file")
    return rel(str(rf).format(seed=seed)) if rf else None


def emb_name(cfg: dict, seed: int) -> str:
    se = cfg["state_embs"]  # 抽出条件をファイル名に入れ、条件の違う dataset を誤って再利用しない
    rf = rows_file_path(cfg, seed)
    if rf:  # 行番号を明示指定した場合は、その内容の sha256 を名前に入れる
        h = hashlib.sha256(rf.read_bytes()).hexdigest()[:8]
        return f"state-embs-seed{seed}-rows{h}.dataset"
    inc = se.get("include_donors")
    tag = "-inc" + hashlib.sha256(json.dumps(inc, sort_keys=True).encode()).hexdigest()[:8] if inc else ""
    return f"state-embs-seed{seed}-n{se.get('max_cells_per_state')}-d{se.get('max_cells_per_donor')}{tag}.dataset"


def state_emb_dataset(tok: Path, cfg: dict, seed: int, out: Path) -> Path:
    """状態重心用に、状態ごと max_cells_per_state・donor ごと max_cells_per_donor で抽出した dataset。
    state_embs.rows_file（{"rows": {状態: [dataset の行番号]}}）がある場合は、その行だけを使う
    （深さをそろえた重心など、抽出規則を実験側の解析スクリプトで事前登録したいとき。EXP-20260927-003）。"""
    if out.exists():
        return out
    from datasets import load_from_disk
    se, key = cfg["state_embs"], cfg["data"]["state_key"]
    dkey = cfg.get("cell_selection", {}).get("donor_key", "donor_id")
    ds = load_from_disk(str(tok))
    rf = rows_file_path(cfg, seed)
    if rf:
        spec = json.loads(rf.read_text())["rows"]
        states_col = [str(x) for x in ds[key]]
        keep = []
        for s in [cfg["states"]["start_state"], cfg["states"]["goal_state"]] + cfg["states"].get("alt_states", []):
            rows = spec.get(s)
            if not rows:
                raise SystemExit(f"state_embs.rows_file: 状態 {s} の行がありません（{rf}）")
            bad = [int(i) for i in rows if states_col[int(i)] != s]
            if bad:
                raise SystemExit(f"state_embs.rows_file: 状態 {s} と一致しない行 {bad[:5]}…（{rf}）")
            keep += [int(i) for i in rows]
        if len(set(keep)) != len(keep):
            raise SystemExit(f"state_embs.rows_file: 行番号が重複しています（{rf}）")
        ds.select(sorted(keep)).save_to_disk(str(out))
        return out
    states, donors = ds[key], [str(x) for x in ds[dkey]]
    rng = random.Random(seed)
    keep = []
    inc = {k: {str(x) for x in v} for k, v in (se.get("include_donors") or {}).items()}
    for s in [cfg["states"]["start_state"], cfg["states"]["goal_state"]] + cfg["states"].get("alt_states", []):
        idx = [i for i, x in enumerate(states) if x == s and (s not in inc or donors[i] in inc[s])]
        if s in inc:
            miss = inc[s] - {donors[i] for i in idx}
            if miss or not idx:
                raise SystemExit(f"state_embs.include_donors.{s}: 細胞のない donor {sorted(miss)}")
        keep += balanced_subset(idx, donors, se.get("max_cells_per_state"), se.get("max_cells_per_donor"), rng)
    ds.select(sorted(keep)).save_to_disk(str(out))
    return out


def state_emb_dataset_with_rowid(data_file: Path) -> Path:
    """precheck 専用: 一意な row_id 列（並べ替え前の行番号）を付けた dataset を用意する。
    EmbExtractor.extract_embs は内部で downsample_and_sort（length 降順・安定ソート）により
    行を並べ替えて返すため、(state, donor) だけでは同一 donor 内の複数細胞を一意に対応付けられない
    （EXP-20260927-002 で発覚。plan.md 末尾の「逸脱の記録」参照）。row_id を emb_label に含めることで、
    cell_embs.csv 側にも row_id 列が残り、analysis 側で ID 結合できる。追加の forward は発生しない
    （このデータセット自体は変わらず列が増えるだけ）。"""
    out = data_file.parent / (data_file.name.removesuffix(".dataset") + "-rowid.dataset")
    if out.exists():
        return out
    from datasets import load_from_disk
    ds = load_from_disk(str(data_file))
    ds = ds.map(lambda ex, idx: {"row_id": idx}, with_indices=True)
    ds.save_to_disk(str(out))
    return out


def common_cell_mode(cfg: dict) -> bool:
    """cell_selection.mode == "common"（全遺伝子を同じ細胞集合で摂動する）か。既定は False（従来どおり）。"""
    return str((cfg.get("cell_selection") or {}).get("mode", "expressing")).lower() == "common"


def cells_dir_name(cfg: dict, seed: int) -> str:
    """抽出条件をディレクトリ名に入れ、条件の違う細胞集合を誤って再利用しない。"""
    sel = cfg["cell_selection"]
    suffix = "-common" if common_cell_mode(cfg) else ""
    grp = donor_group(cfg, seed)
    if grp is not None:  # 提供者の群で絞るときは、群の中身をディレクトリ名に入れる（EXP-20261007-001）
        suffix += "-g" + hashlib.sha256(json.dumps(sorted(grp)).encode()).hexdigest()[:8]
    return f"cells-seed{seed}-n{sel.get('max_cells_per_gene')}-d{sel.get('max_cells_per_donor')}{suffix}"


def donor_group(cfg: dict, seed: int) -> list[str] | None:
    """cell_selection.donor_groups（{seed: [donor, ...]}）があれば、その seed で摂動に使う提供者の一覧。
    無ければ None（従来どおり。出発点の全提供者から選ぶ）。

    EXP-20261007-001 で追加（2026-10-07）。出発点の提供者を互いに重ならない群に分け、
    seed ごとに別の群の細胞だけを摂動する。**状態の重心はこの設定の影響を受けない**
    （重心は state_embs の設定だけで決まる）。common モードでのみ使える（validate で検査）。"""
    groups = (cfg.get("cell_selection") or {}).get("donor_groups")
    if not groups:
        return None
    g = groups.get(str(seed), groups.get(seed))
    if g is None:
        raise SystemExit(f"cell_selection.donor_groups に seed {seed} の群がありません")
    return [str(x) for x in g]


def select_cells(tok: Path, cfg: dict, gene: str, seed: int, out: Path):
    """遺伝子ごとに、start 細胞のうちその遺伝子が入力に含まれる細胞だけを、donor あたり
    max_cells_per_donor 以下・合計 max_cells_per_gene 以下で抽出して保存する。

    cell_selection.mode == "common" のときは**遺伝子で絞らず**、start 状態の全細胞から
    同じ donor 上限・合計上限で 1 つだけ細胞集合を作る（乱数の鍵は f"{seed}-common" で
    遺伝子に依存しない）。呼び出し側は out を遺伝子によらない 1 つのパスにするので、
    **すべての遺伝子がまったく同じ細胞・同じ順序で摂動される**。
    InSilicoPerturber は start 状態 → 遺伝子を含む細胞 → max_ncells（shuffle seed=42 固定）の順に絞るため、
    donor 上限と seed ごとのサブサンプルはここで先に作り、perturber には max_ncells=None で渡す。
    返り値の表は perturber 内部と同じ pu.downsample_and_sort 後の細胞順（pickle のリスト順と一致）。"""
    import pandas as pd
    from datasets import load_from_disk
    from geneformer import perturber_utils as pu
    sel = cfg["cell_selection"]
    dkey = sel.get("donor_key", "donor_id")
    meta_path = out.with_suffix(".cells.csv")
    if out.exists() and meta_path.exists():
        meta = pd.read_csv(meta_path, dtype={"donor_id": str})
        return out, meta
    m, rows, donors = start_matrix(tok, cfg)
    if common_cell_mode(cfg):
        if cfg["perturbation"].get("perturb_type") != "overexpress":
            raise SystemExit("cell_selection.mode: common は perturb_type: overexpress でのみ使える"
                             "（削除型は発現していない細胞では定義できない）")
        idx = list(range(m.shape[0]))
        grp = donor_group(cfg, seed)
        if grp is not None:                     # 提供者の群で絞る（遺伝子には依存しない）
            gs = set(grp)
            idx = [i for i in idx if donors[i] in gs]
            miss = gs - {donors[i] for i in idx}
            if miss:
                raise SystemExit(f"cell_selection.donor_groups[{seed}]: 出発点に細胞の無い提供者 {sorted(miss)[:5]}")
        rng = random.Random(f"{seed}-common")   # 遺伝子に依存しない鍵
    else:
        tid = gene_token_dict().get(gene)
        idx = m[:, tid].nonzero()[0].tolist() if tid is not None and tid < m.shape[1] else []
        rng = random.Random(f"{seed}-{gene}")
    keep = balanced_subset(idx, donors, sel.get("max_cells_per_gene"), sel.get("max_cells_per_donor"), rng)
    keep = [int(rows[i]) for i in keep]  # start 内の番号 → dataset の行番号
    if len(keep) < sel.get("min_cells_per_gene", 1):
        log(f"WARN {gene}: 発現細胞 {len(idx)} から抽出 {len(keep)} < min_cells_per_gene → skip")
        return None, None
    out.parent.mkdir(parents=True, exist_ok=True)
    sub = load_from_disk(str(tok)).select(keep)
    sub.add_column("cell_uid", keep).save_to_disk(str(out))
    srt = pu.downsample_and_sort(load_from_disk(str(out)), None)  # perturber と同じ並べ替え
    meta = pd.DataFrame({"cell_uid": srt["cell_uid"], "donor_id": [str(x) for x in srt[dkey]],
                         "length": srt["length"]})
    meta.to_csv(meta_path, index=False)
    return out, meta


def per_cell_shifts(isp_dir: Path, meta, cfg: dict):
    """raw pickle（状態 → {(token, 'cell_emb'): [細胞ごとの cos 変化]}）を細胞表にする。"""
    import pickle
    f = sorted(isp_dir.glob("*_raw.pickle"))
    if not f:
        return None
    d = pickle.load(open(f[0], "rb"))
    s = cfg["states"]
    df = meta.copy()
    names = {s["goal_state"]: "shift_goal", s["start_state"]: "shift_start"} | \
            {a: f"shift_alt:{a}" for a in s.get("alt_states", [])}
    for state, col in names.items():
        vals = next(iter(d[state].values())) if d.get(state) else []
        if len(vals) != len(df):
            raise SystemExit(f"{isp_dir}: 細胞数が一致しません（pickle {len(vals)} / 抽出 {len(df)}）")
        df[col] = vals
    return df


_LOAD_TIMES: list[float] = []
_ATTN_IMPL: list[str] = []   # 実際に使われた注意機構の実装（2026-10-04）


def wrap_model_loading(cfg: dict) -> None:
    """モデルロード時間を記録し、compute.cache_model=true なら同一モデルを使い回す（遺伝子ごとの再ロードを省く）。
    apply_dtype の後に呼ぶ（bf16 化済みのモデルがキャッシュされる）。"""
    from geneformer import perturber_utils as pu
    original, cache = pu.load_model, {}
    use_cache = bool(cfg.get("compute", {}).get("cache_model", False))

    def load_model_timed(*args, **kwargs):
        k = repr((args, sorted(kwargs.items())))
        if use_cache and k in cache:
            return cache[k]
        t = time.time()
        m = original(*args, **kwargs)
        _LOAD_TIMES.append(time.time() - t)
        # **実際に使われた注意機構の実装を記録する**（2026-10-04）。transformers の版によって
        # 既定が eager / sdpa と変わり、**同じ装置・同じ重みでも bf16 で 0.3125 ずれる**。
        _ATTN_IMPL.append(getattr(getattr(m, "config", None), "_attn_implementation", "unknown"))
        if use_cache:
            cache[k] = m
        return m

    pu.load_model = load_model_timed


def run_one(cfg: dict, data_file: Path, genes, out_dir: Path, emb_layer: int, state_embs, max_ncells="cfg"):
    from geneformer import InSilicoPerturber
    m, p, s, d = cfg["model"], cfg["perturbation"], cfg["states"], cfg["data"]
    isp = call_compat(InSilicoPerturber, {
        "perturb_type": p["perturb_type"],
        "perturb_rank_shift": p.get("perturb_rank_shift"),
        "genes_to_perturb": genes,
        "combos": p.get("combos", 0),
        "anchor_gene": p.get("anchor_gene"),
        "model_type": m["type"],
        "num_classes": m.get("num_classes", 0),
        "emb_mode": p.get("emb_mode", "cls"),
        "cell_emb_style": "mean_pool",
        "filter_data": d.get("filter_data"),
        "cell_states_to_model": cell_states(cfg),
        "state_embs_dict": state_embs,
        "max_ncells": d.get("max_ncells") if max_ncells == "cfg" else max_ncells,
        "emb_layer": emb_layer,
        "forward_batch_size": p.get("forward_batch_size", 100),
        "nproc": p.get("nproc", 4),
        "model_version": m.get("version"),
    }, "InSilicoPerturber")
    out_dir.mkdir(parents=True, exist_ok=True)
    isp.perturb_data(str(rel(m["directory"])), str(data_file), str(out_dir), "isp")


def cell_states(cfg: dict) -> dict:
    s = cfg["states"]
    return {"state_key": cfg["data"]["state_key"], "start_state": s["start_state"],
            "goal_state": s["goal_state"], "alt_states": s.get("alt_states", [])}


def get_state_embs(cfg: dict, data_file: Path, out_dir: Path, emb_layer: int, max_ncells="cfg"):
    from geneformer import EmbExtractor
    m, d, p = cfg["model"], cfg["data"], cfg["perturbation"]
    ee = call_compat(EmbExtractor, {
        "model_type": m["type"], "num_classes": m.get("num_classes", 0),
        "filter_data": d.get("filter_data"),
        "max_ncells": d.get("max_ncells") if max_ncells == "cfg" else max_ncells,
        # 状態重心は perturbation 側と同じ埋め込み種別で作る（cls_and_gene → cls, cell_and_gene → cell）
        "emb_mode": p.get("emb_mode", "cls").replace("_and_gene", ""),
        "emb_layer": emb_layer, "summary_stat": "exact_mean",
        "forward_batch_size": p.get("forward_batch_size", 100), "nproc": p.get("nproc", 4),
        "model_version": m.get("version"),
    }, "EmbExtractor")
    out_dir.mkdir(parents=True, exist_ok=True)
    return ee.get_state_embs(cell_states(cfg), str(rel(m["directory"])), str(data_file), str(out_dir), "state_embs")


def precheck(cfg: dict, tok: Path, work: Path, results: Path) -> dict:
    """事前チェック: start/goal/alt の細胞ごとの埋め込みでロジスティック回帰の k-fold AUC を出す。
    donor 単位の分割（GroupKFold）を主、層化分割を参考値とする。閾値は config の precheck セクション。"""
    import numpy as np
    from geneformer import EmbExtractor
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import GroupKFold, StratifiedKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    pc, m, d, p = cfg["precheck"], cfg["model"], cfg["data"], cfg["perturbation"]
    key = d["state_key"]
    dkey = (cfg.get("cell_selection") or {}).get("donor_key", "donor_id")
    data_file = state_emb_dataset(tok, cfg, pc.get("seed", 0), work / emb_name(cfg, pc.get("seed", 0)))
    data_file = state_emb_dataset_with_rowid(data_file)   # 結合の不具合修正（emb_label に row_id を追加）
    out = results / "precheck"
    out.mkdir(parents=True, exist_ok=True)
    ee = call_compat(EmbExtractor, {
        "model_type": m["type"], "num_classes": m.get("num_classes", 0), "filter_data": None, "max_ncells": None,
        "emb_mode": p.get("emb_mode", "cls").replace("_and_gene", ""), "emb_layer": pc.get("emb_layer", -1),
        "emb_label": [key, dkey, "row_id"], "summary_stat": None,
        "forward_batch_size": p.get("forward_batch_size", 100), "nproc": p.get("nproc", 4),
        "model_version": m.get("version"),
    }, "EmbExtractor")
    t = time.time()
    embs = ee.extract_embs(str(rel(m["directory"])), str(data_file), str(out), "cell_embs")
    t_emb = time.time() - t
    X_all = embs.drop(columns=[key, dkey, "row_id"]).to_numpy(dtype=np.float64)
    s = cfg["states"]
    res = {"n_cells": embs[key].value_counts().to_dict(),
           "n_donors": embs.groupby(key)[dkey].nunique().to_dict(),
           "embedding_seconds": round(t_emb, 1), "folds": pc.get("folds", 5), "comparisons": {}}
    for other, thr in [(s["goal_state"], pc.get("min_auc_goal", 0.8))] + \
                      [(a, pc.get("min_auc_alt", 0.6)) for a in s.get("alt_states", [])]:
        mask = embs[key].isin([s["start_state"], other]).to_numpy()
        X, y = X_all[mask], (embs.loc[mask, key] == s["start_state"]).astype(int).to_numpy()
        g = embs.loc[mask, dkey].astype(str).to_numpy()
        r = {"threshold": thr}
        for name, cv in [("group_by_donor", GroupKFold(n_splits=pc.get("folds", 5))),
                         ("stratified", StratifiedKFold(n_splits=pc.get("folds", 5), shuffle=True,
                                                        random_state=pc.get("seed", 0)))]:
            aucs = []
            for tr, te in cv.split(X, y, groups=g):
                if len(set(y[te])) < 2:
                    continue  # donor 分割で片方のクラスしかない fold は AUC を定義できない
                clf = make_pipeline(StandardScaler(), LogisticRegression(C=pc.get("C", 1.0), max_iter=5000))
                clf.fit(X[tr], y[tr])
                aucs.append(roc_auc_score(y[te], clf.predict_proba(X[te])[:, 1]))
            r[name] = {"auc_mean": round(float(np.mean(aucs)), 4) if aucs else None,
                       "auc_sd": round(float(np.std(aucs, ddof=1)), 4) if len(aucs) > 1 else None,
                       "n_folds_used": len(aucs), "fold_aucs": [round(float(a), 4) for a in aucs]}
        main_auc = r["group_by_donor"]["auc_mean"]
        r["pass"] = bool(main_auc is not None and main_auc >= thr)
        res["comparisons"][f"{s['start_state']}_vs_{other}"] = r
    (out / "precheck.json").write_text(json.dumps(res, indent=2, ensure_ascii=False) + "\n")
    return res


def run_stats(cfg: dict, isp_dir: Path, out_dir: Path, genes):
    from geneformer import InSilicoPerturberStats
    st = call_compat(InSilicoPerturberStats, {
        "mode": cfg["stats"]["mode"], "genes_perturbed": genes,
        "combos": cfg["perturbation"].get("combos", 0), "anchor_gene": cfg["perturbation"].get("anchor_gene"),
        "cell_states_to_model": cell_states(cfg),
    }, "InSilicoPerturberStats")
    out_dir.mkdir(parents=True, exist_ok=True)
    st.get_stats(str(isp_dir), None, str(out_dir), "stats")


def summarize(results: Path, roles: dict[str, str]) -> Path:
    """stats CSV 群を標準化した summary.csv にまとめる。"""
    import pandas as pd
    rows = []
    for csv in results.glob("run-*/**/stats*.csv"):
        run = csv.relative_to(results).parts[0]
        df = pd.read_csv(csv)
        colmap = {}
        for c in df.columns:
            lc = c.lower()
            if lc in ("ensembl_id", "gene") and "gene" not in colmap.values():
                colmap[c] = "gene"
            elif lc == "gene_name":
                colmap[c] = "gene_name"
            elif "shift_to_goal" in lc:
                colmap[c] = "shift"
            elif "goal" in lc and "fdr" in lc:
                colmap[c] = "fdr"
            elif "goal" in lc and "pval" in lc:
                colmap[c] = "pval"
            elif lc == "n_detections":
                colmap[c] = "n_detections"
        df = df.rename(columns=colmap)
        parts = csv.relative_to(results).parts
        if "gene" not in df.columns and "genes" in parts:
            # 1 遺伝子ずつ perturb した場合、stats CSV に遺伝子列がないのでディレクトリ名から補う
            df["gene"] = parts[parts.index("genes") + 1]
        keep = [c for c in ("gene", "gene_name", "shift", "pval", "fdr", "n_detections") if c in df.columns]
        df = df[keep].copy()
        df["run"] = run
        rows.append(df)
    if not rows:
        raise SystemExit("stats CSV が見つかりません")
    out = pd.concat(rows, ignore_index=True)
    # seed を別プロセスで回した場合も、各 run の negatives.json から null の role を補う（EXP-20260926-002 案 B）
    for nj in results.glob("run-*/negatives.json"):
        roles = {g: "negative" for g in json.loads(nj.read_text())} | roles
    out["role"] = out["gene"].map(roles).fillna("screen")
    # cell_selection を使った場合: 細胞ごとの表（cells.csv）を結合し、細胞数・donor 数・alt 方向のシフトを付ける
    cells = []
    for csv in results.glob("run-*/genes/*/cells.csv"):
        c = pd.read_csv(csv, dtype={"donor_id": str})
        c.insert(0, "gene", csv.parent.name)
        c.insert(0, "run", csv.relative_to(results).parts[0])
        cells.append(c)
    if cells:
        cells = pd.concat(cells, ignore_index=True)
        cells["role"] = cells["gene"].map(roles).fillna("screen")
        cells.to_csv(results / "cells.csv", index=False)
        cols = [c for c in cells.columns if c.startswith("shift_")]
        agg = cells.groupby(["run", "gene"]).agg(n_cells=("cell_uid", "size"), n_donors=("donor_id", "nunique"),
                                                 **{f"{c}_mean": (c, "mean") for c in cols}).reset_index()
        out = out.merge(agg, on=["run", "gene"], how="left")
        out["seed"] = out["run"].str.extract(r"run-(\d+)-", expand=False).astype("Int64")
    path = results / "summary.csv"
    out.to_csv(path, index=False)
    return path


def mock(cfg: dict, results: Path) -> Path:
    """geneformer なしで合成データを作り、validate/plot を試せるようにする。"""
    import numpy as np
    import pandas as pd
    rng = np.random.default_rng(0)
    pos = cfg.get("controls", {}).get("positive") or [f"ENSGPOS{i:05d}" for i in range(5)]
    targets = cfg["perturbation"]["genes_to_perturb"]
    targets = targets if isinstance(targets, list) else [f"ENSGTGT{i:05d}" for i in range(200)]
    neg = [f"ENSGNEG{i:05d}" for i in range(cfg.get("controls", {}).get("negative_random_n", 50))]
    true = {g: rng.normal(0.02, 0.01) for g in pos}
    true |= {g: rng.normal(0, 0.006) for g in targets}
    true |= {g: rng.normal(0, 0.004) for g in neg}
    rows = []
    for seed in cfg["robustness"].get("seeds", [0, 1, 2]):
        for g, mu in true.items():
            shift = mu + rng.normal(0, 0.003)
            z = shift / 0.005
            p = float(math.erfc(abs(z) / math.sqrt(2)))  # 両側 p
            rows.append({"gene": g, "shift": shift, "pval": p, "n_detections": int(rng.integers(50, 2000)),
                         "run": f"run-{seed}"})
    df = pd.DataFrame(rows)
    # BH
    for run, sub in df.groupby("run"):
        pv = sub["pval"].to_numpy(); order = np.argsort(pv); n = len(pv)
        fdr = np.empty(n); fdr[order] = np.minimum.accumulate((pv[order] * n / np.arange(1, n + 1))[::-1])[::-1]
        df.loc[sub.index, "fdr"] = np.clip(fdr, 0, 1)
    roles = {g: "positive" for g in pos} | {g: "negative" for g in neg} | {g: "target" for g in targets}
    df["role"] = df["gene"].map(roles)
    results.mkdir(parents=True, exist_ok=True)
    path = results / "summary.csv"
    df.to_csv(path, index=False)
    (results / "MOCK_DATA_DO_NOT_REPORT").write_text("このディレクトリの結果は合成データです\n")
    return path


# ---------------------------------------------------------------------- main
def apply_run(cfg: dict, name: str | None) -> str | None:
    """runs: {名前: {states: {...}}} があれば、指定 run の states を上書きする（同じ設計を start だけ替えて回す）。"""
    runs = cfg.get("runs")
    if not runs:
        if name:
            raise SystemExit(f"--run {name} が指定されましたが config に runs がありません")
        return None
    if name is None:
        raise SystemExit(f"config に runs {list(runs)} があります。--run で 1 つ選んでください")
    if name not in runs:
        raise SystemExit(f"runs に {name} がありません（{list(runs)}）")
    cfg["states"] = {**cfg["states"], **runs[name].get("states", {})}
    return name


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--mock", action="store_true")
    ap.add_argument("--precheck", action="store_true", help="状態の分離 AUC だけを計算する（config の precheck）")
    ap.add_argument("--run", default=None, help="config の runs から 1 つ選ぶ")
    ap.add_argument("--approved", action="store_true")
    args = ap.parse_args()

    cfg = load_config(Path(args.config))
    errs = validate(cfg)
    est = estimate(cfg)
    log(f"見積もり: {json.dumps(est, ensure_ascii=False)}")
    if errs:
        for e in errs:
            log(f"ERROR {e}")
        sys.exit(1)
    if args.dry_run:
        for name, r in (cfg.get("runs") or {}).items():
            log(f"run {name}: states={json.dumps({**cfg['states'], **r.get('states', {})}, ensure_ascii=False)}")
        log("dry-run OK")
        return

    ed = exp_dir(cfg)
    base = cfg.get("output", {}).get("results_dir") or ("results-smoke" if args.smoke else "results")
    run_name = None if (args.mock or args.precheck) else apply_run(cfg, args.run)
    results = ed / base / run_name if run_name else ed / base
    work = ed / "work"
    work.mkdir(parents=True, exist_ok=True)

    if args.mock:
        log(f"mock summary → {mock(cfg, ed / 'results-mock')}")
        return

    if os.environ.get("LABRETARY_MODE") == "night" and est["needs_G2_approval"] and not args.approved:
        sys.exit("夜間モード: G2 承認が必要なジョブです（--approved が必要）")

    # geneformer.device は GENEFORMER_DEVICE 環境変数を参照する
    device = resolve_device(cfg)
    os.environ["GENEFORMER_DEVICE"] = device
    dtype = apply_dtype(cfg)
    single_process_datasets(cfg)
    if apply_capture(cfg):
        log("capture.embeddings: 元/摂動後の CLS 埋め込みを genes/<ENSG>/emb/ に保存します（fp32）")
    if apply_logit_capture(cfg):
        log("capture.classifier_logits: 分類ヘッドの logits（元/摂動後）を genes/<ENSG>/emb/ に保存します（fp32）")
    wrap_model_loading(cfg)
    results.mkdir(parents=True, exist_ok=True)
    (results / "device.json").write_text(json.dumps(
        {"requested": cfg.get("compute", {}).get("device", "auto"), "resolved": device,
         "override_env": os.environ.get("LABRETARY_DEVICE"), "host": os.uname().nodename,
         "available": available_devices(), "dtype": dtype,
         # **ここから下は「比べる前に突き合わせる条件」**（2026-10-04）。
         # Mac の機体照合が落ちたとき、原因を MPS と断定して止めたが、実際は bf16 で、
         # 注意機構の実装が違うだけで同じ装置内でも 0.3125 ずれていた。
         # 記録していなかったので、後から突き合わせることもできなかった。
         "numerics": _numerics_record()}, indent=2) + "\n")
    log(f"device: {device} / dtype: {dtype}")

    if args.smoke:
        cfg["data"]["max_ncells"] = 50
        cfg["robustness"]["seeds"] = cfg["robustness"]["seeds"][:1]
        cfg["robustness"]["emb_layers"] = cfg["robustness"].get("emb_layers", [-1])[:1]
        cfg["controls"]["negative_random_n"] = min(2, cfg.get("controls", {}).get("negative_random_n", 0))
        if cfg.get("cell_selection"):
            sel_s = cfg["cell_selection"]
            sel_s["max_cells_per_gene"] = min(50, sel_s.get("max_cells_per_gene") or 50)
            # smoke は「全工程が通るか」を見るためのもの。min_cells_per_gene を本番の値の
            # まま残すと、max を 50 に絞った結果すべての遺伝子が skip され、摂動が 1 回も
            # 動かないまま「成功」して見える（2026-09-28 修正）。
            min_before = sel_s.get("min_cells_per_gene", 1)
            sel_s["min_cells_per_gene"] = min(min_before, 5)
            if sel_s["min_cells_per_gene"] != min_before:
                log(f"smoke: min_cells_per_gene {min_before} -> {sel_s['min_cells_per_gene']}（本番の下限では全遺伝子が skip されるため）")

    tok = tokenize(cfg, work)
    if args.precheck:
        res = precheck(cfg, tok, work, results)
        log(f"precheck: {json.dumps(res['comparisons'], ensure_ascii=False)}")
        return

    g = cfg["perturbation"]["genes_to_perturb"]
    ctl = cfg.get("controls", {})
    pos = list(ctl.get("positive", []))
    neg = list(ctl.get("negative_fixed", []))
    fixed_role = ctl.get("fixed_role", "negative")
    extra = [(x, grp["role"]) for grp in ctl.get("extra_fixed", []) or [] for x in grp["genes"]]
    neg_all = set(neg) | {x for x, _ in extra}  # null の抽出から除く固定遺伝子
    sampler = ctl.get("negative_sampler") or {}
    sel = cfg.get("cell_selection")
    se = cfg.get("state_embs")
    roles: dict[str, str] = {}
    timing = results / "timing.jsonl"
    wd = work / run_name if run_name else work

    for seed in cfg["robustness"]["seeds"]:
        set_seed(seed)
        data_file = subsample(tok, seed, cfg["robustness"].get("subsample_frac", 1.0), work)
        for layer in cfg["robustness"].get("emb_layers", [-1]):
            run_dir = results / f"run-{seed}-L{layer}"
            log(f"=== {run_dir.name} ===")
            t0 = time.time()
            if se:  # donor 上限つきで抽出した状態重心用の細胞（seed ごと）
                emb_file = state_emb_dataset(tok, cfg, seed, wd / emb_name(cfg, seed))
                state_embs = get_state_embs(cfg, emb_file, run_dir / "state_embs", layer, max_ncells=None)
            else:
                state_embs = get_state_embs(cfg, data_file, run_dir / "state_embs", layer)
            with open(timing, "a") as f:
                f.write(json.dumps({"seed": seed, "step": "state_embs", "seconds": round(time.time() - t0, 2)}) + "\n")
            if g == "all":
                run_one(cfg, data_file, "all", run_dir / "isp", layer, state_embs)
                run_stats(cfg, run_dir / "isp", run_dir / "stats", "all")
                roles |= {x: "positive" for x in pos} | {x: fixed_role for x in neg} | dict(extra)
            else:
                targets = g[:3] if args.smoke else g
                n_rand = ctl.get("negative_random_n", 0)
                if sampler.get("method") == "explicit":
                    rand = [x for x in sampler["genes"] if x not in set(targets) | set(pos) | neg_all][:n_rand or None]
                elif sampler.get("method") == "freq_matched":
                    det = start_detection(tok, cfg, work)
                    refs = targets if sampler.get("match_to", "targets") == "targets" else targets + pos
                    ref_genes = list(sampler.get("ref_genes") or [])  # controls_only 用: 頻度の参照（回さず null からも除く）
                    refs = ref_genes or refs
                    rand = sample_negative_matched(det, n_rand, refs, set(targets) | set(pos) | neg_all | set(ref_genes),
                                                   seed + sampler.get("seed_offset", 0), sampler)
                else:
                    rand = sample_negative_controls(tok, n_rand, set(targets) | set(pos), seed)
                (run_dir).mkdir(parents=True, exist_ok=True)
                (run_dir / "negatives.json").write_text(json.dumps(rand) + "\n")
                plan = [(x, "target") for x in targets] + [(x, "positive") for x in pos] + \
                       [(x, fixed_role) for x in neg] + extra + [(x, "negative") for x in rand]
                for gene, role in plan:
                    roles[gene] = role
                    gd = run_dir / "genes" / gene
                    if (gd / "stats").exists() and (not sel or (gd / "cells.csv").exists()) and \
                            (not _CAPTURE["on"] or (gd / "emb" / "pert.npy").exists()) and \
                            (not _CAPTURE["logits"] or (gd / "emb" / "logits_pert.npy").exists()):
                        continue  # 再開時: 済んだ遺伝子は飛ばす
                    t0, nl = time.time(), len(_LOAD_TIMES)
                    capture_reset()
                    if sel:
                        cdir = cells_dir_name(cfg, seed)
                        # common モードでは遺伝子によらない 1 つのファイルを使い回す（＝顔ぶれが同一であることが設計で保証される）
                        stem = "common" if common_cell_mode(cfg) else gene
                        f_cells, meta = select_cells(tok, cfg, gene, seed, wd / cdir / f"{stem}.dataset")
                        if f_cells is None:
                            continue
                        run_one(cfg, f_cells, [gene], gd / "isp", layer, state_embs, max_ncells=None)
                    else:
                        run_one(cfg, data_file, [gene], gd / "isp", layer, state_embs)
                    t_isp = time.time() - t0
                    run_stats(cfg, gd / "isp", gd / "stats", [gene])
                    n_cells = None
                    if sel:
                        cells = per_cell_shifts(gd / "isp", meta, cfg)
                        cells.to_csv(gd / "cells.csv", index=False)
                        n_cells = len(cells)
                    if _CAPTURE["on"] or _CAPTURE["logits"]:
                        capture_save(gd, run_dir, n_cells)
                    with open(timing, "a") as f:
                        f.write(json.dumps({"seed": seed, "gene": gene, "role": role, "n_cells": n_cells,
                                            "isp_seconds": round(t_isp, 2),
                                            "model_load_seconds": round(sum(_LOAD_TIMES[nl:]), 2)}) + "\n")
    path = summarize(results, roles)
    (results / "config.sha256").write_text(cfg["_sha256"] + "\n")
    # **device.json を書き直す。** 最初に書いた時点ではモデルをまだ読んでおらず、
    # attn_implementation（結果を左右する条件）が空のままになるため（2026-10-04）。
    d = json.loads((results / "device.json").read_text())
    d["numerics"] = _numerics_record()
    (results / "device.json").write_text(json.dumps(d, indent=2) + "\n")
    log(f"完了: {path} / 数値条件: {json.dumps(d['numerics'], ensure_ascii=False)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env bash
# 実験の再現情報を <EXP>/results/env.json に保存する
set -euo pipefail
EXP="${1:?usage: capture_env.sh experiments/EXP-...}"
mkdir -p "$EXP/results"
# 実験を走らせる環境の版を記録するため、プロジェクトの .venv を優先
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${PY:-$([[ -x "$ROOT/.venv/bin/python" ]] && echo "$ROOT/.venv/bin/python" || echo python3)}"
"$PY" - "$EXP" <<'PY'
import json, os, platform, subprocess, sys, hashlib, datetime, pathlib
exp = pathlib.Path(sys.argv[1])
def sh(c):
    try: return subprocess.run(c, shell=True, capture_output=True, text=True, timeout=20).stdout.strip()
    except Exception as e: return f"error: {e}"
pkgs = {}
for p in ["geneformer", "torch", "transformers", "datasets", "anndata", "scanpy", "numpy", "pandas", "scipy"]:
    try:
        from importlib.metadata import version; pkgs[p] = version(p)
    except Exception: pkgs[p] = None
cfg = exp / "config.yaml"
dev_json = exp / "results" / "device.json"  # run_isp.py が書く実使用デバイス
def accel():
    try:
        import torch
        return {"cuda": torch.cuda.is_available(), "mps": torch.backends.mps.is_available()}
    except Exception as e:
        return {"error": str(e)}
def cfg_compute(key):
    try:
        import yaml
        return (yaml.safe_load(cfg.read_text()).get("compute") or {}).get(key)
    except Exception:
        return None
env = {
  "captured_at": datetime.datetime.now().astimezone().isoformat(),
  "host": platform.node(), "platform": platform.platform(), "python": platform.python_version(),
  "git_commit": sh("git rev-parse HEAD"), "git_dirty": bool(sh("git status --porcelain")),
  "packages": pkgs,
  "gpu": sh("nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader") or "none",
  "accelerators": accel(),
  "cpu_model": sh("sysctl -n machdep.cpu.brand_string") if platform.system() == "Darwin" else sh("grep -m1 'model name' /proc/cpuinfo | cut -d: -f2"),
  "config_device": cfg_compute("device"),
  "config_dtype": cfg_compute("dtype") or "fp32",
  "geneformer_device_env": os.environ.get("GENEFORMER_DEVICE"),
  "run_device": json.loads(dev_json.read_text()) if dev_json.exists() else None,
  "config_sha256": hashlib.sha256(cfg.read_bytes()).hexdigest() if cfg.exists() else None,
}
(exp / "results" / "env.json").write_text(json.dumps(env, indent=2, ensure_ascii=False))
print(f"saved {exp/'results'/'env.json'}")
PY

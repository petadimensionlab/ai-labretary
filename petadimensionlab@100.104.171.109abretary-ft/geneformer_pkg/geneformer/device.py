"""
Geneformer device detection.

Routes compute selection across the best available accelerator with an explicit
override via the ``GENEFORMER_DEVICE`` environment variable.

Supported device names:

- ``cuda`` : NVIDIA CUDA **or** AMD ROCm. ROCm exposes the CUDA API, so a ROCm
  build reports ``torch.cuda.is_available() == True`` and is driven through the
  same code path.
- ``mps``  : Apple Silicon Metal Performance Shader (Mac MPS).
- ``cpu``  : CPU fallback.

Auto-detection priority is CUDA/ROCm > MPS > CPU. Set ``GENEFORMER_DEVICE`` to
force a specific backend (e.g. ``GENEFORMER_DEVICE=mps`` on a Mac that also has
CUDA available via a remote/egpu setup).
"""
from __future__ import annotations

import os

import torch

_VALID = {"cuda", "mps", "cpu"}


def _detect() -> str:
    """Return the best available device name without consulting the override."""
    if torch.cuda.is_available():
        return "cuda"
    mps = getattr(torch.backends, "mps", None)
    if mps is not None and mps.is_available():
        return "mps"
    return "cpu"


def get_device() -> str:
    """Return the resolved device name (``cuda``, ``mps``, or ``cpu``).

    Honors the ``GENEFORMER_DEVICE`` environment variable override.
    """
    override = os.environ.get("GENEFORMER_DEVICE", "").strip().lower()
    if override:
        if override not in _VALID:
            raise ValueError(
                f"GENEFORMER_DEVICE={override!r} is invalid; "
                f"choose from {sorted(_VALID)}"
            )
        return override
    return _detect()


def get_device_obj() -> torch.device:
    """Return a ``torch.device`` for the resolved accelerator."""
    return torch.device(get_device())


def move_to_device(model: torch.nn.Module) -> torch.nn.Module:
    """Move ``model`` to the resolved accelerator if it is not already there."""
    device = get_device_obj()
    model_device = next(model.parameters()).device
    if model_device.type != device.type:
        model.to(device)
    return model


def empty_cache() -> None:
    """Best-effort accelerator cache clear. No-op on MPS/CPU."""
    if torch.cuda.is_available():
        torch.cuda.empty_cache()


def manual_seed_all(seed: int) -> None:
    """Seed all available generators (CUDA/ROCm aware)."""
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
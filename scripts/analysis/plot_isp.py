#!/usr/bin/env python3
"""summary.csv から標準図を作る: (1) 効果量ランク図 (2) run 間散布図 (3) コントロール分布"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

COL = {"positive": "#d62728", "negative": "#7f7f7f", "target": "#1f77b4", "screen": "#1f77b4"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", required=True)
    ap.add_argument("--results", default="results")
    a = ap.parse_args()
    exp = Path(a.exp); res = exp / a.results; fig_dir = exp / "figures"; fig_dir.mkdir(exist_ok=True)
    df = pd.read_csv(res / "summary.csv")
    agg = df.groupby(["gene", "role"], as_index=False)["shift"].agg(["mean", "std"]).reset_index()
    agg = agg.sort_values("mean", ascending=False).reset_index(drop=True)
    tag = " (MOCK)" if (res / "MOCK_DATA_DO_NOT_REPORT").exists() else ""

    fig, ax = plt.subplots(figsize=(7, 4))
    for role, sub in agg.groupby("role"):
        ax.errorbar(sub.index, sub["mean"], yerr=sub["std"], fmt="o", ms=3, lw=0.5, color=COL.get(role, "k"), label=role)
    ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("rank"); ax.set_ylabel("shift toward goal state (mean ± SD)")
    ax.set_title("Perturbation effect ranking" + tag); ax.legend(); fig.tight_layout()
    fig.savefig(fig_dir / "rank.png", dpi=200); plt.close(fig)

    wide = df.pivot_table(index="gene", columns="run", values="shift")
    if wide.shape[1] >= 2:
        r1, r2 = wide.columns[:2]
        fig, ax = plt.subplots(figsize=(4, 4))
        ax.scatter(wide[r1], wide[r2], s=5)
        ax.set_xlabel(r1); ax.set_ylabel(r2); ax.set_title("Run-to-run reproducibility" + tag)
        fig.tight_layout(); fig.savefig(fig_dir / "reproducibility.png", dpi=200); plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 3.5))
    for role, sub in agg.groupby("role"):
        ax.hist(sub["mean"], bins=30, alpha=0.5, color=COL.get(role, "k"), label=role, density=True)
    ax.set_xlabel("shift"); ax.set_title("Controls vs targets" + tag); ax.legend(); fig.tight_layout()
    fig.savefig(fig_dir / "controls.png", dpi=200); plt.close(fig)
    print("\n".join(str(p) for p in sorted(fig_dir.glob("*.png"))))


if __name__ == "__main__":
    main()

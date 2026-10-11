#!/usr/bin/env python3
"""SHA256SUMS の作成・照合（標準ライブラリのみ。リモートでもローカルでも同じ）。EXP-20260926-003。

  python verify_sha256.py check DIR        # DIR/SHA256SUMS と全ファイルを照合（不一致・欠落があれば exit 1）
  python verify_sha256.py make DIR         # DIR/SHA256SUMS を作り直す（学習結果など）
  python verify_sha256.py check DIR --no-extra   # SHA256SUMS に載っていないファイルがあっても失敗にする

形式は `sha256sum` と同じ（'<hash>  <相対パス>'）。`sha256sum -c SHA256SUMS`（Linux）/ `shasum -a 256 -c`（macOS）でも照合できる。
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

NAME = "SHA256SUMS"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def make(d: Path) -> int:
    files = [p for p in sorted(d.rglob("*")) if p.is_file() and p.name != NAME]
    (d / NAME).write_text("".join(f"{sha256(p)}  {p.relative_to(d).as_posix()}\n" for p in files))
    print(f"{d / NAME}: {len(files)} files")
    return 0


def check(d: Path, no_extra: bool) -> int:
    f = d / NAME
    if not f.exists():
        print(f"NG: {f} がありません")
        return 1
    listed, bad = set(), []
    for line in f.read_text().splitlines():
        if not line.strip():
            continue
        h, rp = line.split("  ", 1)
        listed.add(rp)
        p = d / rp
        if not p.exists():
            bad.append(f"missing {rp}")
        elif sha256(p) != h:
            bad.append(f"mismatch {rp}")
    extra = [p.relative_to(d).as_posix() for p in d.rglob("*") if p.is_file() and p.name != NAME
             and p.relative_to(d).as_posix() not in listed]
    if extra:
        print(f"{'NG' if no_extra else 'WARN'}: SHA256SUMS にないファイル {len(extra)} 件: {extra[:5]}")
        if no_extra:
            bad += [f"extra {x}" for x in extra]
    if bad:
        print(f"NG: {len(bad)} 件: {bad[:10]}")
        return 1
    print(f"OK: {len(listed)} files match ({d})")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1] not in ("check", "make"):
        sys.exit(__doc__)
    d = Path(sys.argv[2]).resolve()
    sys.exit(make(d) if sys.argv[1] == "make" else check(d, "--no-extra" in sys.argv))

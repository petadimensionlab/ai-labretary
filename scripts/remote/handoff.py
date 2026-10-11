#!/usr/bin/env python3
"""ある seed が**正常に完了したときだけ**、同じ機体で次の seed を始める。

## なぜあるか

1 台で複数の seed を順に回すとき、`isp_remote.sh` は次の seed を自動で始めない。
2026-09-30、実行担当の引き継ぎメモだけが頼りだったため、seed82 が 23:03 に終わったのに
seed83 の開始が 23:47 になり **43 分遅れた**。しかもそのとき動いていた完了の見張りは
「3 つの seed がすべて終わったら知らせる」作りで、**引き継ぎが死ぬと永久に鳴らない**。
気づけたのは PI が状況を尋ねたからで、見落としに上限が無かった。

さらにその引き継ぎを macOS の `/bin/bash`（3.2）で書いたため、
`$VAR` の直後に日本語が来ると変数名を取り違える不具合で落ちた。
**この Mac のシェルは zsh、プロジェクトの計算は Python。bash 3.2 を選ぶ理由は無かった。**

## 考え方

- **正常に完了したときだけ次を始める。** 遺伝子数が足りない・エラーがある・tmux が残っている
  のいずれかなら**始めずに理由を書いて終わる**（中断規則を飛び越えない）
- **判断はしない。** 異常なら人（Labretary → PI）に戻す
- 終了コード: 0 = 次を開始した / 2 = 開始しなかった（異常）/ 3 = 見張り自身の異常

## 使い方

    .venv/bin/python scripts/remote/handoff.py \
        --exp EXP-20261001-001 --host coi-next1@100.94.57.17 \
        --after 82 --start 83 --expect 717 \
        --remote-gf-home '$HOME/workspace/research/Geneformer' \
        --log logs/handoff-20261001.log

読み取りは ssh ごしの `ls` と `grep` だけ。開始は `isp_remote.sh run` を呼ぶ。
"""
from __future__ import annotations

import argparse
import datetime as dt
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def now() -> str:
    return dt.datetime.now().strftime("%F %T")


class Log:
    def __init__(self, path: Path):
        self.path = path
        path.parent.mkdir(parents=True, exist_ok=True)

    def __call__(self, msg: str) -> None:
        line = f"[{now()}] {msg}"
        print(line, flush=True)
        with self.path.open("a", encoding="utf-8") as f:
            f.write(line + "\n")


def ssh(host: str, cmd: str, timeout: int = 40) -> tuple[int, str]:
    try:
        r = subprocess.run(["ssh", "-o", "ConnectTimeout=15", "-o", "BatchMode=yes", host, cmd],
                           capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return 124, ""
    return r.returncode, r.stdout.strip()


def gene_count(host: str, exp: str, seed: str) -> int | None:
    rc, out = ssh(host, f"ls ~/labretary-remote/experiments/{exp}/results-seed{seed}"
                        f"/run-{seed}-L-1/genes 2>/dev/null | wc -l")
    if rc != 0:
        return None
    return int(out) if out.isdigit() else None


def session_alive(host: str, seed: str) -> bool | None:
    rc, out = ssh(host, "tmux ls 2>/dev/null")
    if rc not in (0, 1):
        return None
    return f"seed{seed}" in out


def finished_clean(host: str, exp: str, seed: str, expect: int, log: Log) -> bool:
    """遺伝子数・エラー件数・tmux の 3 つで、正常に終わったかを確かめる。"""
    n = gene_count(host, exp, seed)
    if n != expect:
        log(f"★ seed{seed} の遺伝子数が {n} / {expect}。正常完了ではない")
        return False
    rc, out = ssh(host, f"grep -ac 'ERROR\\|Traceback' "
                        f"~/labretary-remote/experiments/{exp}/logs/*seed{seed}* 2>/dev/null | "
                        f"awk -F: '{{s+=$NF}} END {{print s+0}}'")
    if rc == 0 and out.isdigit() and int(out) > 0:
        log(f"★ seed{seed} のログにエラーが {out} 件ある")
        return False
    if session_alive(host, seed):
        log(f"★ seed{seed} の tmux セッションがまだ残っている")
        return False
    log(f"seed{seed} は正常に完了（{n} / {expect}、エラー 0 件、tmux 終了）")
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", required=True)
    ap.add_argument("--host", required=True)
    ap.add_argument("--after", required=True, help="この seed の完了を待つ")
    ap.add_argument("--start", required=True, help="完了したら始める seed")
    ap.add_argument("--expect", type=int, required=True, help="1 seed あたりの遺伝子数")
    ap.add_argument("--remote-gf-home", default=None)
    ap.add_argument("--interval", type=int, default=300)
    ap.add_argument("--log", default="logs/handoff.log")
    a = ap.parse_args()

    log = Log(Path(a.log) if Path(a.log).is_absolute() else ROOT / a.log)
    log(f"監視開始: {a.host} の seed{a.after} が {a.expect} 遺伝子で正常完了したら seed{a.start} を始める")

    # --- 完了を待つ（遺伝子数が達する、または tmux が消える）
    while True:
        n = gene_count(a.host, a.exp, a.after)
        alive = session_alive(a.host, a.after)
        if n == a.expect:
            log(f"seed{a.after} が {a.expect} 遺伝子に到達")
            break
        if alive is False:
            log(f"seed{a.after} の tmux セッションが消えた（そのとき {n} / {a.expect}）")
            break
        log(f"待機中: seed{a.after} {n} / {a.expect}")
        time.sleep(a.interval)

    time.sleep(20)   # 書き込みの取りこぼしを避ける

    # --- 正常完了でなければ始めない
    if not finished_clean(a.host, a.exp, a.after, a.expect, log):
        log(f"★★ seed{a.start} を開始しない。中断規則に当たる可能性がある。")
        log("★★ Labretary が状況を確かめて PI に報告すること。")
        return 2

    # --- 次を始める
    env = dict(os.environ)
    if a.remote_gf_home:
        env["REMOTE_GF_HOME"] = a.remote_gf_home
    cfg = f"experiments/{a.exp}/config-seed{a.start}.yaml"
    log(f"seed{a.start} を開始する: {cfg}")
    r = subprocess.run(["bash", "scripts/remote/isp_remote.sh", "run", a.host, cfg],
                       cwd=str(ROOT), env=env, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        log(f"★★ 開始に失敗（終了コード {r.returncode}）: {r.stderr.strip()[:400]}")
        return 2
    log(f"開始した: {r.stdout.strip()[:200]}")

    # --- 本当に動き出したかを確かめる
    for _ in range(12):
        time.sleep(30)
        n = gene_count(a.host, a.exp, a.start)
        if n and n >= 1:
            log(f"seed{a.start} が動き出した（{n} 遺伝子）")
            return 0
    log(f"★★ seed{a.start} を開始したが、6 分たっても出力が出ない。Labretary が確かめること")
    return 2


if __name__ == "__main__":
    sys.exit(main())

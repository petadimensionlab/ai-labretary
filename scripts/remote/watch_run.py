#!/usr/bin/env python3
"""リモートで走っている摂動の見張り — **進んでいないことを検知する**。

## なぜあるか

2026-09-30、seed82 → seed83 の引き継ぎが「次を開始する行」で落ちた。
そのとき動いていた見張りは **「3 つの seed がすべて 600 になったら知らせる」** という作りで、
引き継ぎが死ねば条件が永久に満たされない。**異常を検知できない仕組みだった。**
気づけたのは PI が状況を尋ねたからで、放置すれば見落としに上限が無かった。

さらに、その見張りを **macOS の /bin/bash（3.2）** で書いたために、
連想配列が使えない・`$VAR` の直後に日本語が来ると変数名を取り違える、という
古い bash 固有の不具合で 2 度落とした。**この Mac のシェルは zsh であり、
そもそも bash 3.2 を選ぶ理由が無かった。** 道具の選び方の誤りである。

## 考え方

1. **「期待どおり進んだこと」ではなく「進んでいないこと」を見る。**
   合計が一定時間まったく増えなければ、原因が何であれ異常として鳴る
   （引き継ぎの死亡・計算の異常終了・通信断・機体の停止）。
2. **判断はしない。** 鳴ったら止まって知らせるだけ。再開するかどうかは Labretary が決め、
   中断規則に当たるなら PI に報告する。
3. **終了コードで区別する**: 0 = 正常に完了 / 9 = 止まっている / 3 = 見張り自身の異常。

## 使い方

    .venv/bin/python scripts/remote/watch_run.py --exp EXP-20260930-001 \
        --seed 81:petadimensionlab@100.104.171.109 \
        --seed 82:coi-next1@100.94.57.17 \
        --seed 83:coi-next1@100.94.57.17 \
        --expect 600 --stall-min 25 --log logs/watch-run.log

読み取りだけを行う（ssh で `ls | wc -l` を叩く）。リモートにも手元にも何も書かない。
"""
from __future__ import annotations

import argparse
import datetime as dt
import subprocess
import sys
import time
from pathlib import Path


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


RESULTS_TEMPLATE = "results-seed{seed}"   # --results-template で変えられる（既定は従来どおり）


def count(host: str, exp: str, seed: str, timeout: int = 30) -> int | None:
    """その seed の遺伝子ディレクトリ数を数える。数えられなければ None（通信断と区別する）。"""
    rdir = RESULTS_TEMPLATE.format(seed=seed, host=host)
    remote = (f"ls ~/labretary-remote/experiments/{exp}/{rdir}"
              f"/run-{seed}-L-1/genes 2>/dev/null | wc -l")
    try:
        r = subprocess.run(["ssh", "-o", "ConnectTimeout=15", "-o", "BatchMode=yes", host, remote],
                           capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None
    if r.returncode != 0:
        return None
    s = r.stdout.strip()
    return int(s) if s.isdigit() else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", required=True)
    ap.add_argument("--seed", action="append", required=True,
                    help="SEED:HOST の形。複数指定できる")
    ap.add_argument("--expect", type=int, default=600, help="1 seed あたりの遺伝子数")
    ap.add_argument("--stall-min", type=int, default=25,
                    help="合計がこれだけの分数まったく増えなければ異常とみなす")
    ap.add_argument("--interval", type=int, default=300)
    ap.add_argument("--log", default="logs/watch-run.log")
    ap.add_argument("--results-template", default="results-seed{seed}",
                    help="リモートの結果ディレクトリ名（{seed} を展開）。既定は従来どおり results-seed{seed}"
                         "（2026-10-08、EXP-20261007-001 で追加。機体ごとに results-<機体> へ書く実験のため）")
    a = ap.parse_args()
    global RESULTS_TEMPLATE
    RESULTS_TEMPLATE = a.results_template

    pairs = []
    for s in a.seed:
        if ":" not in s:
            print(f"--seed は SEED:HOST の形で指定する: {s}", file=sys.stderr)
            return 3
        seed, host = s.split(":", 1)
        pairs.append((seed, host))

    root = Path(__file__).resolve().parents[2]
    log = Log(root / a.log if not Path(a.log).is_absolute() else Path(a.log))
    total = a.expect * len(pairs)
    log(f"監視開始: {a.exp} / seed {[p[0] for p in pairs]} / 目標 {total} "
        f"（合計が {a.stall_min} 分増えなければ異常として鳴る）")

    prev, last_change = -1, time.time()
    unreachable_since: float | None = None

    while True:
        counts = {seed: count(host, a.exp, seed) for seed, host in pairs}
        missing = [s for s, v in counts.items() if v is None]
        got = {s: (v or 0) for s, v in counts.items()}
        ssum = sum(got.values())
        t = time.time()

        if missing:
            # 通信断そのものは異常だが、一時的なこともある。続いたときだけ鳴らす
            unreachable_since = unreachable_since or t
            idle = int((t - unreachable_since) / 60)
            log(f"数えられない seed がある: {missing}（{idle} 分）")
            if idle >= a.stall_min:
                log(f"★★ {idle} 分にわたって {missing} の状態を数えられない。通信か機体の異常。")
                log("★★ Labretary が確かめること。**自分で再開せず、中断規則に当たるかを先に判断する。**")
                return 9
        else:
            unreachable_since = None

        if ssum >= total:
            log(f"合計 {ssum} / {total} に到達。正常に完了したので監視を終える")
            log("内訳: " + " / ".join(f"seed{s} {v}" for s, v in got.items()))
            return 0

        if ssum != prev:
            prev, last_change = ssum, t
        else:
            idle = int((t - last_change) / 60)
            if idle >= a.stall_min:
                log(f"★★ 止まっている: 合計 {ssum} / {total} が {idle} 分まったく増えていない")
                log("★★ 内訳: " + " / ".join(f"seed{s} {v}" for s, v in got.items()))
                log("★★ 考えられる原因: 次の seed の開始が行われていない／計算が異常終了した／通信断")
                log("★★ Labretary が原因を確かめること。**自分で再開せず、中断規則に当たるかを先に判断する。**")
                return 9

        log(f"合計 {ssum} / {total}（" + " / ".join(f"seed{s} {v}" for s, v in got.items())
            + f"）最後に増えてから {int((t - last_change) / 60)} 分")
        time.sleep(a.interval)


if __name__ == "__main__":
    sys.exit(main())

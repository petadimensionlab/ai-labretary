#!/usr/bin/env bash
# EXP-20260930-001: 3 つの seed がすべて 600 遺伝子で終わるまで待ち、節目だけをログに書く。
#
# **判断はしない。** 異常を見つけても自分で直さず、ログに書いて終わる。
# 見込みを過ぎても終わらないときに気づけるよう、経過時間も一緒に記録する。
#
# 注意: macOS の /bin/bash は 3.2 で **連想配列（declare -A）が使えない**。
# 使うと起動直後に落ちる（2026-09-30 に実際に落とした）。ここでは使わない。
# また bash 3.2 は **$VAR の直後に日本語などの多バイト文字が来ると変数名を取り違える**。
# 必ず ${VAR} と書く（2026-09-30 に "N?: unbound variable" で落とした）。
set -u

CF5B="petadimensionlab@100.104.171.109"
C083="coi-next1@100.94.57.17"
EXP="EXP-20260930-001"
ROOT="/Users/petadimensionlab/workspace/research/ai-labretary"
LOG="$ROOT/logs/wait-exp004-2026-09-30.log"
N=600
START=$(date +%s)

say() { printf '[%s +%dm] %s\n' "$(date '+%F %T')" "$((($(date +%s)-START)/60))" "$*" >> "$LOG"; }

cnt() {  # cnt HOST SEED
  ssh -o ConnectTimeout=15 -o BatchMode=yes "$1" \
    "ls ~/labretary-remote/experiments/$EXP/results-seed$2/run-$2-L-1/genes 2>/dev/null | wc -l" 2>/dev/null \
    | tr -d ' \r'
}

d81=0; d82=0; d83=0
say "監視開始（3 seed が ${N} 遺伝子になるまで）"

while :; do
  a=$(cnt "$CF5B" 81); b=$(cnt "$C083" 82); c=$(cnt "$C083" 83)
  [ -z "$a" ] && a="?"; [ -z "$b" ] && b="?"; [ -z "$c" ] && c="?"

  if [ "$a" = "$N" ] && [ "$d81" = "0" ]; then d81=1; say "seed81 が完了（${N} / ${N}）"; fi
  if [ "$b" = "$N" ] && [ "$d82" = "0" ]; then d82=1; say "seed82 が完了（${N} / ${N}）"; fi
  if [ "$c" = "$N" ] && [ "$d83" = "0" ]; then d83=1; say "seed83 が完了（${N} / ${N}）"; fi

  if [ "$d81$d82$d83" = "111" ]; then
    say "3 seed すべて完了。回収・照合・評価へ進める"
    exit 0
  fi
  say "経過: seed81 $a / seed82 $b / seed83 $c （各 ${N}）"
  sleep 900
done

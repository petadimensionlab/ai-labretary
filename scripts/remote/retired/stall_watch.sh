#!/usr/bin/env bash
# 止まっていることを検知する見張り（EXP-20260930-001）。
#
# ## なぜあるか
#
# 2026-09-30、seed82 → seed83 の引き継ぎスクリプトが「次を開始する行」で落ちた。
# そのとき動いていた完了の見張りは **「3 seed がすべて 600 になったら知らせる」** という作りで、
# 引き継ぎが死ねば条件が永久に満たされない。**つまり異常を検知できない仕組みだった。**
# 気づけたのは PI が状況を尋ねたからで、放置すれば見落としに上限が無かった。
#
# ## 考え方
#
# **「期待どおり進んだこと」ではなく「進んでいないこと」を見る。**
# 3 つの seed の遺伝子数の合計が一定時間まったく増えなければ、
# 原因が何であれ（引き継ぎの死亡・計算の異常終了・通信断・機体の停止）異常として鳴る。
#
# **判断はしない。** 鳴ったら止まって知らせるだけ。再開するかどうかは Labretary が決め、
# 中断規則に当たるなら PI に報告する。
#
# 注意: macOS の bash 3.2 は連想配列が使えず、`$VAR` の直後に日本語が来ると
# 変数名を取り違える。**必ず ${VAR} と書く。**
set -u

CF5B="petadimensionlab@100.104.171.109"
C083="coi-next1@100.94.57.17"
EXP="EXP-20260930-001"
ROOT="/Users/petadimensionlab/workspace/research/ai-labretary"
LOG="${ROOT}/logs/stall-watch-2026-09-30.log"
TOTAL=1800                 # 600 遺伝子 × 3 seed
STALL_MIN=25               # これだけ増えなければ異常とみなす（1 遺伝子 ≈ 24 秒なので十分な余裕）
INTERVAL=300

say() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >> "${LOG}"; }

cnt() {  # cnt HOST SEED
  ssh -o ConnectTimeout=15 -o BatchMode=yes "$1" \
    "ls ~/labretary-remote/experiments/${EXP}/results-seed$2/run-$2-L-1/genes 2>/dev/null | wc -l" 2>/dev/null \
    | tr -d ' \r'
}

num() { case "$1" in ''|*[!0-9]*) echo 0 ;; *) echo "$1" ;; esac; }

prev=-1
last_change=$(date +%s)
say "監視開始: 合計が ${STALL_MIN} 分まったく増えなければ異常として鳴る（目標 ${TOTAL}）"

while :; do
  a=$(num "$(cnt "${CF5B}" 81)")
  b=$(num "$(cnt "${C083}" 82)")
  c=$(num "$(cnt "${C083}" 83)")
  sum=$((a + b + c))
  now=$(date +%s)

  if [ "${sum}" -ge "${TOTAL}" ]; then
    say "合計 ${sum} / ${TOTAL} に到達。正常に終わったので監視を終える"
    exit 0
  fi

  if [ "${sum}" -ne "${prev}" ]; then
    prev=${sum}; last_change=${now}
  else
    idle=$(( (now - last_change) / 60 ))
    if [ "${idle}" -ge "${STALL_MIN}" ]; then
      say "★★ 止まっている: 合計 ${sum} / ${TOTAL} が ${idle} 分まったく増えていない"
      say "★★ 内訳 seed81 ${a} / seed82 ${b} / seed83 ${c}"
      say "★★ Labretary が原因を確かめること。引き継ぎの死亡・異常終了・通信断のいずれか。"
      say "★★ **自分で再開せず、中断規則に当たるかを判断してから動く。**"
      exit 9
    fi
  fi
  say "合計 ${sum} / ${TOTAL}（seed81 ${a} / seed82 ${b} / seed83 ${c}）最後に増えてから $(( (now - last_change) / 60 )) 分"
  sleep "${INTERVAL}"
done

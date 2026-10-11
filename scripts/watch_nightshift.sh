#!/usr/bin/env bash
# 最新の nightshift ログを追いかける。夜間でなければキューを表示
cd "$(dirname "$0")/.." || exit 1
while true; do
  latest="$(ls -d logs/nightshift/*/ 2>/dev/null | sort | tail -1)"
  if [[ -n "$latest" && -f "$latest/shift.log" ]] && pgrep -f labretary-nightshift.sh >/dev/null; then
    clear; echo "== nightshift: $latest =="; tail -n 40 "$latest/shift.log"
  else
    clear; echo "== nightshift: 停止中 =="; echo; echo "[night-queue]"
    ls tasks/night-queue/*.md 2>/dev/null | grep -v README | sed 's/^/  /'
    [[ -n "$latest" ]] && { echo; echo "[前回: $latest]"; tail -n 8 "$latest/shift.log" 2>/dev/null; }
  fi
  sleep 15
done

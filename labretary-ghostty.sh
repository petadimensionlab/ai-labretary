#!/usr/bin/env bash
# Labretary を Ghostty で可視化する（tmux の代わり）
#
#   ./labretary-ghostty.sh              # ウィンドウ方式: 役割ごとに Ghostty ウィンドウを開く（macOS / Linux）
#   ./labretary-ghostty.sh --splits     # 分割方式: 1 ウィンドウを 4 分割（macOS のみ・アクセシビリティ許可が必要）
#   ./labretary-ghostty.sh --job EXP-…  # 実行中ジョブのログウィンドウを追加
#
# 画面構成
#   ┌─────────────────────┬─────────────────────┐
#   │ Labretary (claude)  │ status (10 秒更新)  │
#   ├─────────────────────┼─────────────────────┤
#   │ actions (行動ログ)  │ nightshift / job    │
#   └─────────────────────┴─────────────────────┘
set -euo pipefail
PROJECT="$(cd "$(dirname "$0")" && pwd)"
CONF="$PROJECT/ghostty/labretary.conf"
OS="$(uname)"

# 各ペインのコマンド（どちらの方式でも共通）
CMD_MAIN="claude"
CMD_STATUS="while true; do clear; bash scripts/lab_status.sh; sleep 10; done"
CMD_ACTIONS="bash scripts/watch_actions.sh"
CMD_NIGHT="bash scripts/watch_nightshift.sh"
job_cmd() { echo "tail -n 100 -F experiments/$1/logs/run.log"; }

open_window() { # $1 title, $2 command
  local title="$1" cmd="cd '$PROJECT' && $2"
  if [[ "$OS" == "Darwin" ]]; then
    # macOS: 新しい Ghostty インスタンスとして開く（-e は最後に置く）
    open -na Ghostty.app --args --config-file="$CONF" --title="$title" \
      --working-directory="$PROJECT" -e bash -lc "$cmd"
  else
    command -v ghostty >/dev/null || { echo "ghostty が PATH にありません"; exit 1; }
    nohup ghostty --config-file="$CONF" --title="$title" --working-directory="$PROJECT" \
      -e bash -lc "$cmd" >/dev/null 2>&1 &
  fi
  sleep 0.6
}

splits_macos() {
  # Ghostty の keybind（labretary.conf で定義）をキーストロークで送って分割・コマンド入力する
  open -a Ghostty.app --args --config-file="$CONF"
  sleep 1.5
  osascript <<OSA
on typeLine(t)
  tell application "System Events" to keystroke t
  tell application "System Events" to key code 36
end typeLine
tell application "Ghostty" to activate
delay 0.5
tell application "System Events"
  keystroke "n" using {command down}
  delay 0.8
end tell
typeLine("cd '$PROJECT' && $CMD_MAIN")
delay 0.5
tell application "System Events" to keystroke "d" using {command down}
delay 0.6
typeLine("cd '$PROJECT' && $CMD_STATUS")
tell application "System Events" to keystroke "d" using {command down, shift down}
delay 0.6
typeLine("cd '$PROJECT' && $CMD_NIGHT")
tell application "System Events" to keystroke "[" using {command down}
delay 0.3
tell application "System Events" to keystroke "[" using {command down}
delay 0.3
tell application "System Events" to keystroke "d" using {command down, shift down}
delay 0.6
typeLine("cd '$PROJECT' && $CMD_ACTIONS")
OSA
}

case "${1:-}" in
  --splits)
    [[ "$OS" == "Darwin" ]] || { echo "--splits は macOS のみ。Linux はウィンドウ方式を使ってください"; exit 1; }
    splits_macos ;;
  --job)
    open_window "Labretary · job ${2:?EXP-ID}" "$(job_cmd "$2")" ;;
  ""|--windows)
    open_window "Labretary"               "$CMD_MAIN"
    open_window "Labretary · status"      "$CMD_STATUS"
    open_window "Labretary · actions"     "$CMD_ACTIONS"
    open_window "Labretary · nightshift"  "$CMD_NIGHT" ;;
  *) sed -n '2,15p' "$0"; exit 1 ;;
esac

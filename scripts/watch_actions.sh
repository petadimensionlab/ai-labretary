#!/usr/bin/env bash
# 行動ログ（logs/actions-*.log, JSON Lines）を見やすく追いかける。
# このプロジェクトには行動ログを書く hook を設定していない（2026-09-27 確認）。hook を設定しない限りログは空のまま。
# 行動ログを書く hook（例: log-action.sh）が logs/actions-current.log を当日のファイルへの symlink に張り替えるので、
# tail -F（名前で追跡）で日付をまたいでも追従できる。
cd "$(dirname "$0")/.." || exit 1
mkdir -p logs; f="logs/actions-$(date +%F).log"; touch "$f"
[[ -L logs/actions-current.log ]] || ln -sf "$(basename "$f")" logs/actions-current.log
tail -n 30 -F logs/actions-current.log 2>/dev/null | while IFS= read -r line; do
  if command -v jq >/dev/null; then
    printf '%s\n' "$line" | jq -r '"\(.ts[11:19] // "")  \(.mode // "")  \((.tool // .event // "") | .[0:14])  \(.input // "")"' 2>/dev/null || echo "$line"
  else echo "$line"; fi
done

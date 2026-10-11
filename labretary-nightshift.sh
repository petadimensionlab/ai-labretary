#!/usr/bin/env bash
# Labretary 夜間自律シフト
#   tasks/night-queue/*.md（または GitHub Issues）を 1 件ずつ Claude Code のヘッドレス実行で消化する。
#   最後に学習セッションで learnings をメモリへ転記し、朝の要約を書く（状況の一覧は lab-reporter 担当に委託）。
#
#   使い方:  ./labretary-nightshift.sh            # ファイルキュー
#            SOURCE=github ./labretary-nightshift.sh   # label: labretary-night の Issue
#            DRY=1 ./labretary-nightshift.sh      # 何を実行するか表示だけ
set -uo pipefail
cd "$(dirname "$0")"
PROJECT="$(pwd)"

INTERVAL="${INTERVAL:-600}"      # タスク間の休憩（秒）= 10 分（レート制限保護）
MAX_RUNS="${MAX_RUNS:-20}"       # 1 晩の最大実行回数
END_AT="${END_AT:-05:45}"        # この時刻を過ぎたら新しいタスクを始めない（朝の作業枠を残す）
SOURCE="${SOURCE:-files}"        # files | github
MODEL="${MODEL:-}"               # 空ならデフォルトモデル
DRY="${DRY:-0}"

export LABRETARY_MODE=night
DATE="$(date +%F)"
LOGDIR="$PROJECT/logs/nightshift/$DATE"; mkdir -p "$LOGDIR"
MAIN_LOG="$LOGDIR/shift.log"
log() { echo "[$(date '+%T')] $*" | tee -a "$MAIN_LOG"; }

# --- スリープ防止（スクリプト終了で自動解除）
if [[ "$(uname)" == "Darwin" ]] && command -v caffeinate >/dev/null; then
  caffeinate -i -w $$ & log "caffeinate on (pid $!)"
elif command -v systemd-inhibit >/dev/null && systemd-inhibit --list >/dev/null 2>&1 && [[ -z "${LABRETARY_INHIBITED:-}" ]]; then
  export LABRETARY_INHIBITED=1
  exec systemd-inhibit --what=idle:sleep --why="Labretary nightshift" "$0" "$@"
fi

command -v claude >/dev/null || { log "claude CLI が見つかりません"; exit 1; }

past_end() { [[ "$(date +%H:%M)" > "$END_AT" && "$(date +%H:%M)" < "12:00" ]]; }

fm() { # frontmatter の値: fm <file> <key> <default>
  awk -v k="$2" 'NR==1&&$0=="---"{f=1;next} f&&$0=="---"{exit} f&&$1==k":"{sub(/^[^:]*:[ ]*/,"");print;exit}' "$1" | grep . || echo "$3"
}

next_task() {
  if [[ "$SOURCE" == "github" ]]; then
    gh issue list --label labretary-night --state open --limit 1 --json number,title,body \
      --jq '.[0] | select(.) | "\(.number)\t\(.title)"'
  else
    ls tasks/night-queue/*.md 2>/dev/null | grep -v README.md | sort | head -1
  fi
}

build_prompt() { # $1 = task file
  local f="$1" type approved body
  type="$(fm "$f" type custom)"; approved="$(fm "$f" approved false)"
  body="$(awk 'NR==1&&$0=="---"{f=1;next} f&&$0=="---"{f=0;next} !f' "$f")"
  cat <<P
あなたは Labretary です。現在は夜間自律シフト（LABRETARY_MODE=night）で、PI は不在です。
CLAUDE.md §8 のファイルを読んでから、次のタスクを 1 件だけ実行してください。

- タスク種別: $type
- 対象: $(fm "$f" target "-")
- 大規模ジョブの承認: $approved  （true の場合のみ run_isp.py に --approved を付けてよい。false なら dry-run / smoke まで）

## 依頼内容
$body

## 夜間ルール
- PI への確認が必要な判断は行わず、tasks/queue.md の「承認待ち」に理由付きで書いて終了する
- 発見・失敗は logs/learnings-$(date +%F).jsonl に必ず 1 行ずつ記録する
- 最後に、作成・更新したファイル一覧（ls で実在確認済み）と 3 行サマリーを出力する
P
}

run_claude() { # $1 prompt, $2 out json, $3 max_turns
  local args=(-p "$1" --output-format json --max-turns "$3" --permission-mode acceptEdits)
  [[ -n "$MODEL" ]] && args+=(--model "$MODEL")
  claude "${args[@]}" > "$2" 2>> "$MAIN_LOG"
}

log "=== nightshift start (source=$SOURCE, max=$MAX_RUNS, interval=${INTERVAL}s, end=$END_AT) ==="
runs=0
while (( runs < MAX_RUNS )); do
  past_end && { log "終了時刻 $END_AT を過ぎたので停止"; break; }
  task="$(next_task)"; [[ -z "$task" ]] && { log "キューが空です"; break; }
  runs=$((runs+1)); id="$(printf '%03d' $runs)"

  if [[ "$SOURCE" == "github" ]]; then
    num="${task%%$'\t'*}"; tf="$LOGDIR/issue-$num.md"
    gh issue view "$num" --json title,body --jq '"---\ntype: custom\ntarget: \(.title)\napproved: false\nmax_turns: 40\n---\n\(.body)"' > "$tf"
  else tf="$task"; fi

  log "[$id] start: $tf"
  if [[ "$DRY" == "1" ]]; then build_prompt "$tf"; break; fi

  out="$LOGDIR/$id-$(basename "$tf" .md).json"
  start=$(date +%s)
  run_claude "$(build_prompt "$tf")" "$out" "$(fm "$tf" max_turns 40)"; rc=$?
  log "[$id] done rc=$rc in $(( $(date +%s) - start ))s → $out"

  if grep -Eqi 'rate.?limit|usage limit|quota' "$out" 2>/dev/null; then
    log "レート制限を検知。タスクはキューに残して停止"; break
  fi
  if [[ "$SOURCE" == "github" ]]; then
    gh issue comment "$num" --body "Labretary nightshift: rc=$rc, log=$(basename "$out")" >/dev/null
    [[ $rc -eq 0 ]] && gh issue edit "$num" --remove-label labretary-night --add-label labretary-done >/dev/null
  else
    [[ $rc -eq 0 ]] && mv "$tf" "tasks/done/$(date +%F)-$(basename "$tf")" || mv "$tf" "${tf%.md}.failed.md"
  fi
  sleep "$INTERVAL"
done

# --- 学習セッション: learnings → memory、日報
if [[ "$DRY" != "1" ]] && ! past_end; then
  log "learning session"
  run_claude "夜間モードの学習セッションです（CLAUDE.md §7）。状況の一覧は lab-reporter 担当に委託して作らせてください。今夜の logs/nightshift/$DATE/ の結果と logs/learnings-$(date +%F).jsonl を読み、guardrails / labretary-level のドメインパターンに昇格すべきものがあれば更新し、朝 PI が最初に読む要約を outputs/reports/daily/$(date +%F)-night.md に書いてください。" \
    "$LOGDIR/zz-learning.json" 40
fi
log "=== nightshift end (runs=$runs) ==="

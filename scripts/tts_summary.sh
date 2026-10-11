#!/usr/bin/env bash
# Labretary: レポートの要約を音声にする。読み上げ原稿も同じ場所に保存する。
#
#   bash scripts/tts_summary.sh NAME TEXT_FILE [ja|en]
#   echo "本文" | bash scripts/tts_summary.sh NAME - ja
#
# 言語ごとの生成方法（PI 指示 2026-09-26）:
#   ja … OpenRouter の /audio/speech（scripts/tts_openrouter.py、既定 Fish Audio S2.1 Pro Free）。
#        鍵はプロジェクト直下の .env の OPENROUTER_API_KEY。**原稿本文を外部 API に送る**ので、
#        未公開の生データやヒト由来データの識別子を原稿に入れない。
#   en … macOS の say（Ava (Premium)）。ローカル生成のみ。既定。
#
# 原稿の書き方: 話し言葉で書く。記号・表は避け、用語は docs/用語集.md に合わせる。
#   英語例: VCAM1 → "V-CAM one"、0.45 → "zero point four five"
#   日本語例: VCAM1 → ブイカム・ワン、0.45 → れい点よんご
# 環境変数: TTS_VOICE / TTS_RATE（en のみ）、TTS_MODEL / TTS_VOICE_ID（ja のみ）
set -euo pipefail
NAME="${1:?usage: tts_summary.sh NAME TEXT_FILE|- [ja|en]}"
SRC="${2:?TEXT_FILE か - を指定}"
LANG_CODE="${3:-en}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/outputs/audio/$(date +%F)"
mkdir -p "$OUT"
if [[ "$SRC" == "-" ]]; then cat > "$OUT/$NAME.txt"; else cp "$SRC" "$OUT/$NAME.txt"; fi

case "$LANG_CODE" in
  ja)
    ARGS=(--file "$OUT/$NAME.txt" -o "$OUT/$NAME.mp3" --format mp3)
    [[ -n "${TTS_MODEL:-}" ]] && ARGS+=(--model "$TTS_MODEL")
    [[ -n "${TTS_VOICE_ID:-}" ]] && ARGS+=(--voice "$TTS_VOICE_ID")
    "$ROOT/.venv/bin/python" "$ROOT/scripts/tts_openrouter.py" "${ARGS[@]}"
    FILE="$OUT/$NAME.mp3"
    ;;
  en)
    TMP="$(mktemp -t labretary-tts).aiff"
    say -v "${TTS_VOICE:-Ava (Premium)}" -r "${TTS_RATE:-175}" -f "$OUT/$NAME.txt" -o "$TMP"
    afconvert -f m4af -d aac "$TMP" "$OUT/$NAME.m4a"
    rm -f "$TMP"
    FILE="$OUT/$NAME.m4a"
    ;;
  *) echo "言語は ja か en"; exit 1 ;;
esac

dur="$(afinfo "$FILE" | awk '/estimated duration/{printf "%.0f", $3}')"
echo "${FILE}（${dur} 秒）"

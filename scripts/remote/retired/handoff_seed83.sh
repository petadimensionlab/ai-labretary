#!/usr/bin/env bash
# EXP-20260930-001: c083 の seed82 が**正常に完了したときだけ** seed83 を開始する。
#
# なぜあるか: c083 は 2 つの seed を順に回すが、自動で連結されていない。
# 実行担当の引き継ぎメモに頼ると、過去に完了を 3 時間見落とした前例がある（guardrails 2026-09-28）。
#
# **中断規則を無視して先へ進めないこと。** seed82 が 600 遺伝子で終わっていなければ
# seed83 を開始せず、理由を書いて終わる。判断は人（Labretary → PI）に戻す。
#
# 注意: macOS の bash 3.2 は **$VAR の直後に日本語が来ると変数名を取り違える**。
# 必ず ${VAR} と書く（2026-09-30 にこの不具合で seed83 の開始が 43 分遅れた）。
set -u

HOST="coi-next1@100.94.57.17"
EXP="EXP-20260930-001"
ROOT="/Users/petadimensionlab/workspace/research/ai-labretary"
LOG="$ROOT/logs/handoff-seed83-2026-09-30.log"
N_EXPECT=600

say() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOG"; }

count_genes() {
  ssh -o ConnectTimeout=15 -o BatchMode=yes "$HOST" \
    "ls ~/labretary-remote/experiments/$EXP/results-seed82/run-82-L-1/genes 2>/dev/null | wc -l" 2>/dev/null \
    | tr -d ' \r'
}
session_alive() {
  ssh -o ConnectTimeout=15 -o BatchMode=yes "$HOST" "tmux ls 2>/dev/null" 2>/dev/null \
    | grep -q 'isp-config-seed82'
}

say "監視開始: $HOST の seed82（$N_EXPECT 遺伝子で完了）"

# --- seed82 の終了を待つ（tmux が消えるか、600 に達するまで）
while :; do
  n="$(count_genes)"
  [ -z "$n" ] && n="?"
  if [ "$n" = "$N_EXPECT" ]; then
    say "seed82 が $N_EXPECT 遺伝子に到達"
    break
  fi
  if ! session_alive; then
    say "seed82 の tmux セッションが消えた（そのとき $n / $N_EXPECT 遺伝子）"
    break
  fi
  sleep 120
done

# --- 正常完了でなければ開始しない
n="$(count_genes)"
if [ "$n" != "$N_EXPECT" ]; then
  say "★ seed82 が $n / $N_EXPECT で終わっている。**seed83 を開始しない。**"
  say "★ 中断規則に当たる可能性がある。Labretary が状況を確認して PI に報告すること。"
  exit 2
fi

say "seed82 は正常に完了（$n / ${N_EXPECT}）。seed83 を開始する"

cd "$ROOT" || { say "★ 作業ディレクトリへ移動できない"; exit 3; }
REMOTE_GF_HOME='$HOME/workspace/research/Geneformer' \
  bash scripts/remote/isp_remote.sh run "$HOST" "experiments/$EXP/config-seed83.yaml" >>"$LOG" 2>&1
rc=$?
if [ $rc -ne 0 ]; then
  say "★ seed83 の開始に失敗（終了コード ${rc}）。Labretary が確認すること"
  exit 4
fi

sleep 60
n83="$(ssh -o ConnectTimeout=15 -o BatchMode=yes "$HOST" \
  "tmux ls 2>/dev/null | grep -c seed83; ls ~/labretary-remote/experiments/$EXP/results-seed83/run-83-L-1/genes 2>/dev/null | wc -l" 2>/dev/null | tr '\n' ' ')"
say "seed83 を開始した（tmux 数 / 遺伝子数: ${n83}）"
exit 0

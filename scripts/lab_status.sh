#!/usr/bin/env bash
# 研究状況の一覧（lab-reporter 担当・make status・ghostty ダッシュボードで使用）
cd "$(dirname "$0")/.." || exit 1
printf '\033[1m== Labretary status  %s ==\033[0m\n' "$(date '+%F %T')"
printf '\n\033[1m[Experiments]\033[0m\n'
printf '%-20s %-8s %-8s %-10s %-10s\n' EXP plan results validation report
for d in experiments/EXP-*/; do
  [ -d "$d" ] || continue
  e=$(basename "$d")
  a=$(grep -q '^Approved-by: *[^ ]' "$d/plan.md" 2>/dev/null && echo ok || echo -)
  r=$([ -f "$d/results/summary.csv" ] && echo done || { [ -f "$d/logs/pid" ] && kill -0 "$(cat "$d/logs/pid")" 2>/dev/null && echo RUNNING || echo -; })
  v=$([ -f "$d/validation.md" ] && echo ok || echo -)
  s=$(grep -m1 -o 'Status: [a-z]*' "$d/report.md" 2>/dev/null | cut -d' ' -f2); s=${s:--}
  printf '%-20s %-8s %-8s %-10s %-10s\n' "$e" "$a" "$r" "$v" "$s"
done
printf '\n\033[1m[Hypotheses]\033[0m\n'
awk -F'|' 'NR>4 && NF>3 {gsub(/ /,"",$4); c[$4]++} END {for (k in c) printf "  %-14s %d\n", k, c[k]}' hypotheses/tracker.md
printf '\n\033[1m[Night queue]\033[0m\n'
ls tasks/night-queue/*.md 2>/dev/null | grep -v README | sed 's/^/  /' || echo "  (empty)"
printf '\n\033[1m[GPU]\033[0m\n'
if command -v nvidia-smi >/dev/null; then
  nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used,memory.total --format=csv,noheader | sed 's/^/  /'
elif [ "$(uname -s)" = Darwin ]; then
  PY=$([ -x .venv/bin/python ] && echo .venv/bin/python || echo python3)
  mps=$("$PY" -c 'import torch; print("available" if torch.backends.mps.is_available() else "unavailable")' 2>/dev/null || echo "unknown (torch なし)")
  echo "  $(sysctl -n machdep.cpu.brand_string)  MPS: $mps"
else echo "  nvidia-smi なし"; fi

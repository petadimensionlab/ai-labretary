#!/usr/bin/env bash
# 優先度 A データ（outputs/literature/2026-09-25-kidney-TLS-datasets.md）を取得し、サイズと sha256 を記録する
set -uo pipefail
D="$(cd "$(dirname "$0")/.." && pwd)/data/raw/kidney_tls"
while IFS=$'\t' read -r name url; do
  out="$D/$name.h5ad"
  expect=$(curl -sIL "$url" | awk 'tolower($1)=="content-length:"{v=$2} END{gsub(/\r/,"",v); print v}')
  if [ -f "$out" ] && [ "$(stat -f%z "$out")" = "$expect" ]; then echo "SKIP $name (完了済み)"; continue; fi
  echo "GET  $name ($expect B)"
  curl -sSL --retry 5 --retry-delay 10 -C - -o "$out.part" "$url" || { echo "FAIL $name"; continue; }
  got=$(stat -f%z "$out.part")
  if [ "$got" = "$expect" ]; then mv "$out.part" "$out"; echo "OK   $name $got B $(shasum -a 256 "$out" | cut -d' ' -f1)"
  else echo "SIZE_MISMATCH $name got=$got expect=$expect"; fi
done < "$D/urls.tsv"
echo ALL_DONE

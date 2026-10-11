# phase 0 の結果（EXP-20260926-001）

本番 ISP ではない。計画（plan.md §8）の根拠。

- `env.json` — 実行環境（`scripts/capture_env.sh`。同じ内容が `../results/env.json`）
- `precheck/precheck.json` — 状態の分離 AUC（start vs goal 1.00 / start vs alt 0.627、donor 分割 5-fold）、`precheck/cell_embs.csv` は CLS 埋め込み（1,405 細胞）
- `precision-bf16/tls_high/`, `precision-fp32/tls_high/` — dtype だけを変えた matched pair（陽性対照 2 + HK 3 + 発現頻度マッチのランダム 5、各 100 細胞、seed 0）。対象遺伝子（TP53/CDKN1A/LTBR/SPP1/LCN2）は本番前に結果を見ないため含めていない
- `precision_bf16_vs_fp32.json` — 比較（Spearman 0.879 < 0.98 → 本番は fp32）と 1 細胞×1 遺伝子の実測秒数
- ログ: `../../../logs/EXP-20260926-001-phase0-*.log`

> **注記（2026-09-26 Labretary）**: 上の「本番は fp32」は phase 0 時点の提案。PI の決定（「bf16 を標準に」）により **本番は bf16**（forward bf16 + コサイン fp32）で実行した（plan.md §10、report.md §2.2）。

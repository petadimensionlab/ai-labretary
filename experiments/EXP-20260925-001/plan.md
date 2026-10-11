# EXP-20260925-001 実験計画 — パイプライン技術検証（MPS smoke）

- 仮説: なし（技術検証）。生物学的な解釈・結論は出さない
- 作成: Labretary / 2026-09-25

## 1. 問い
ダウンロードした CELLxGENE データ → 前処理 → tokenize → ISP（`--smoke`）が MPS 上で最後まで通るか。1 forward あたりの実測時間を得て、A100 基準の見積もりを補正する。

## 2. データ
| 項目 | 値 |
|---|---|
| データセット | A5a ccRCC + A5b 正常腎（Zhang 2021 PNAS, GSE159115, CELLxGENE 版） |
| 細胞数（start / goal） | macrophage のみ、状態ごと最大 300（ccRCC 300 / normal ≤262） |
| 状態ラベル列 | disease（clear cell renal carcinoma → normal） |
| フィルタ | cell_type = macrophage |

## 3. モデル
Geneformer-V2-316M（Pretrained）、device: mps

## 4. Perturbation 設計
`--smoke`: 50 細胞、1 seed、delete × 3 遺伝子（CXCL13, CCL19, CDKN1A）＋ コントロール

## 5. コントロール
ポジ（仮置き）: CD163 / ネガ: ACTB, GAPDH, ランダム 2

## 6. 評価指標と判定基準
完走すること、`results-smoke/summary.csv` と `device.json` が生成されること、1 forward あたり時間の記録。**shift の値は解釈しない**（50 細胞・1 seed）

## 8. 計算量見積もり
smoke は 50 細胞 × 約 8 遺伝子。G2 対象外（< 1 万細胞、< 2 h）

## 10. 承認
Approved-by: （技術検証・G2 対象外。PI 指示「2 も並行して進めて」2026-09-25 に基づき実行）

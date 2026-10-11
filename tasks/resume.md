# tasks/resume.md — 復帰用チェックポイント

頻繁に止まるため、各ステップの完了ごとにここへ状態を書き、いつでもここから再開できるようにする。
**復帰時は先頭の「次の一手」だけを見て作業を再開する。**

## 現在のスレッド: A6（EXP-20261011-003）埋め込みのみの確認実験

更新日時: 2026-10-11

## 状況（1 行）

計画（EXP-20261011-003）は **PI 回答 2026-10-11 を反映して第 2 回固定**（候補 5・管状系陽性対照 23 個・G2 承認）。
**実行・検証・報告の修正・G3 確定まですべて完了**（2026-10-11）: 埋め込み本番（spark-cf5b・fp32・46,957 核 1 回）→ evaluate_a6.py 本番 → PI 指示で報告を学術用語に全面書き換え（plan/config/evaluate・エージェント定義・検証文書）→ SHA256SUMS 再固定（9 件・全 OK）。
**【G3 確定 2026-10-11】PI が「判定不能」で確定。** 記録: `experiments/EXP-20261011-003/report.md`（状態=確定）・tracker・HYP-20261011-01.md。

## 次の一手（順番どおり）

1. ✅ **実行（simulation-engineer）完了**: prepare_h5ad／smoke／評価コード固定 → 埋め込み本番（GB10/CUDA fp32・46,957 核・1 回）→ axis_feasibility 2 系統（記録のみ）→ evaluate_a6.py 本番（A-1〜A-7 合格）。
2. ✅ **validator 批判的レビュー完了**（validation.md）: 判定不能の提案。主要値を再計算し criteria と一致。P-2(a) の sin 誤り・P-3 の帰属誤りを指摘。
3. ✅ **PI 指示 2026-10-11 の反映（完了）**: 報告を全面的に書き直し・専門用語は学術用語で書く・造語を作らない・遵守事項は総称で。CLAUDE.md §9・用語集・エージェント定義 20 ファイルに反映。plan.md/config.yaml/evaluate_a6.py/embed_a6.py を学術用語に全面置換（`rg` 0 件・判定値不変）→ SHA256SUMS 再固定（変更ログ第 6・7 回）。
4. ✅ **PI へゲート G3 の報告 → 確定**: 「判定不能」で確定（PI 2026-10-11）。`report.md`（状態=確定）・tracker・仮説ファイルに記録。
5. **次回の一手**: 終了。次のテーマ（間質系の分布のそろえ方を総 UMI 数まで含めて揃えた再評価、または GSN・ZFHX3 の探索）は PI 判断で新規に開始。復帰時はここから。

## 完了したこと

1. **G3 確定の記録**: HYP-20261011-01 は判定不能（inconclusive）。EXP-20261011-001/002/003 の report.md（状態=確定）・tracker・仮説ファイルへ記録済み。EXP-003（A6 独立データ確認）も判定不能で確定（間質系の分離が総 UMI 数と検出遺伝子数で説明できる）。
2. **A6 データの読み取り調査**（46,957 核 × 29,702 遺伝子。Myofibroblast 585・Fibroblast 141・PT Injured 2,019・PT Inflammatory 613。候補の 8 遺伝子の検出率を実測）。
3. **PI 決定 2026-10-11 を EXP-20261011-001/plan.md §5.6 と EXP-20261011-003/plan.md §0.3 に記録**: (a) 間質系＋管状系の 2 系統 (b) CCL21 を入れて SLC4A9 は検出が少なく位置を推定できなかった (c) 埋め込みのみ（摂動なし・fp32・1 回）。
4. **EXP-20261011-003 の計画を experiment-planner が作成・事前登録**。SHA256 9 件で固定（`shasum -c` OK）。
5. **plan §12 の 6 項目を PI に諮り、すべて回答済み**（候補 5・管状系陽性対照 23・G2 承認）。plan §10.2 に `Approved-by: PI 2026-10-11` を記録。
6. **literature-scout 完了**: 管状系退化の陽性対照 40 個（材料 Reck2025 補足データ 5・Lake2023・Gisch2024）。

## 参照

- 実験: `experiments/EXP-20261011-003/`（plan.md・config.yaml・analysis/・results/・validation.md・report.md）
- 判定の確定: `experiments/EXP-20261011-003/report.md`（状態=確定、判定不能 G3 2026-10-11）
- ハッシュ固定: `analysis/SHA256SUMS-preregistered.txt`（**9 件**。`shasum -c` で全件 OK）
- 文献成果物: `outputs/literature/2026-10-11-tubular-degenerative-positive-control/`・`2026-10-11-single-cell-depth-terminology.md`
- A6 データ: `data/raw/kidney_tls/A6_UUO_Reck2025.h5ad`（読み取り専用、sha256 8d9cca13…）
- **G2 承認の記録**: plan.md §10.2 の `Approved-by: PI 2026-10-11（G2 保守側・埋め込み ≈1.13 h 承認）`
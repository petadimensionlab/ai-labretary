# EXP-20260926-001 実験計画 — CKD TLS-high 腎の aPT で TP53/CDKN1A/LTBR を delete すると PT_healthy へ近づくか（TLS-low と比べる）

- 仮説: [HYP-20260925-01](../../hypotheses/HYP-20260925-01.md)（G1 採択 2026-09-26, PI）
- 作成: Labretary（experiment-planner）/ 2026-09-26
- 状態: **計画のみ。本番 ISP は未実行。** phase 0（入力作成・dry-run・事前チェック・精度ペア・時間計測）まで実施済み（結果は `results-phase0/`）

## 1. 問い
1. ヒト CKD 腎の aPT のうち、推定 TLS-high の donor に属する細胞（start）で TP53 / CDKN1A / LTBR を delete すると、CLS 埋め込みは PT_healthy（goal）のセントロイドへ近づくか。近づく量は、発現頻度をマッチさせたランダム遺伝子の null より大きいか（run 1: `tls_high`）
2. 同じ摂動を TLS-low の donor の aPT（alt）を start にしたとき（run 2: `tls_low`）と比べ、シフトは 1.2 倍以上か（文脈依存性）

in silico の結果であり、支持されても「Geneformer の埋め込み空間でそう示唆される」以上のことは言わない。

## 2. データ
| 項目 | 値 |
|---|---|
| データセット | A2 KPMP v1.5 scRNA（Lake 2023, `data/raw/kidney_tls/A2_KPMP_v1.5_scRNA_Lake2023.h5ad`, 読み取り専用。raw.X の整数カウントを使用） |
| 前処理 | `scripts/geneformer/prepare_h5ad.py` → `work/input.h5ad`（17,684 細胞 × 30,635 遺伝子）→ tokenize（V2, 4,096 トークン）→ `work/tokenized.dataset` |
| 細胞数（start / goal / alt） | **7,533（14 donor）/ 6,149（22 donor）/ 4,002（8 donor）** — 仮説の値と一致（`work/input.summary.json`） |
| 状態ラベル列 | `isp_state`（**推定ラベル**。下表） |
| フィルタ | author_cell_type ∈ {aPT, PT-S1/S2, PT-S3}、disease ∈ {CKD, normal}。CKD の TLS-mid donor の aPT と CKD の PT-S1/S2・S3 は除外 |
| 分割（FT する場合） | FT なし（Pretrained） |

**ラベルの由来（どれが推定か）**
| ラベル | 定義 | 由来 |
|---|---|---|
| `author_cell_type`（aPT, PT-S1/S2, PT-S3） | 著者のアノテーション | Lake 2023（著者） |
| `disease`（CKD / normal） | CELLxGENE メタデータ | データ提供者 |
| `TLS_group`（high / mid / low） | donor 単位 pseudo-bulk の TLS 8 遺伝子スコアの三分位（基準 = 300 細胞以上の donor） | **推定**（`outputs/literature/2026-09-25-tls-feasibility/donor_A2.csv` の `TLS_group` 列をそのまま写した。再定義していない） |
| `isp_state` | aPT_CKD_TLShigh = aPT ∧ CKD ∧ TLS high / aPT_CKD_TLSlow = aPT ∧ CKD ∧ TLS low / PT_healthy = (PT-S1/S2 ∨ PT-S3) ∧ normal（TLS 群は問わない。hyp_state_counts.csv と同じ 6,149 細胞） | **推定**（上の組み合わせ。組織学的な TLS ではない） |

## 3. モデル
- Geneformer-V2-316M（Pretrained）、`emb_mode: cls`（V2 必須）、emb_layer −1
- device: mps（Apple M4 Pro 64 GB）、**dtype: bf16**（forward bf16 + コサイン fp32）。**PI 決定（2026-09-26「bf16 を標準に」）**。phase 0c では fp32 と Spearman 0.879（../Geneformer の基準 0.98 未達）、上位 5 一致 5/5・符号 10/10。順位の細部には bf16 由来の誤差（mean|Δ| = mean|shift| の 7%）が乗ることを限界としてレポートに明記する
- `compute.cache_model: true`（遺伝子ごとのモデル再ロードをやめた。phase 0 で使用）
- `perturbation.nproc: 1`（MPS 使用後の fork で datasets.map が落ちるため。guardrails 2026-09-26）

## 4. Perturbation 設計
- `perturb_type: delete`、1 遺伝子ずつ（`genes/<ENSG>/`）
- 対象 5: TP53, CDKN1A, LTBR（新規予測）, SPP1, LCN2（データ駆動の参照）
- **細胞の選び方（遺伝子ごと・seed ごと, `cell_selection`）**
  - Geneformer の InSilicoPerturber は「start 状態 → 対象遺伝子が入力に含まれる細胞 → `max_ncells`（shuffle seed=42 固定）」の順に絞る。つまり発現細胞だけが残り、`max_ncells` は発現細胞での絞り込みの**後**にかかるが、**donor 上限はなく、seed で変わらない**（`in_silico_perturber.py` `apply_additional_filters`, `perturber_utils.downsample_and_sort`）
  - そこで run_isp.py 側で、start 細胞のうち遺伝子が tokenized 入力に含まれる細胞から、**donor あたり ≤15、合計 ≤100** を seed ごとに抽出して小さな dataset を作り、perturber には `max_ncells=None` で渡す（`select_cells`）。30 細胞未満しか取れない遺伝子は実行せず記録
  - 細胞の順は perturber と同じ `pu.downsample_and_sort` で並べて `*.cells.csv` に記録し、raw pickle（細胞ごとのコサイン変化）と対応付けて `genes/<ENSG>/cells.csv`（donor_id, shift_goal, shift_start, shift_alt）を保存。→ donor ブートストラップが可能
  - 実測（TLS-high, seed 0）: 対象・対照 10 遺伝子すべてで 100 細胞 / 12〜14 donor を確保（入力に含まれる細胞数: TP53 1,058, CDKN1A 533, LTBR 1,594, SPP1 5,713, LCN2 449, VCAM1 4,674, HAVCR1 659）
- **状態セントロイド（`state_embs`）**: seed ごとに状態あたり ≤500 細胞・donor あたり ≤60 で抽出（seed 0: 500 / 500 / 405）
- run 1 `tls_high`: start = aPT_CKD_TLShigh, goal = PT_healthy, alt = [aPT_CKD_TLSlow]
- run 2 `tls_low`: start = aPT_CKD_TLSlow, goal = PT_healthy, alt = [aPT_CKD_TLShigh]（同じ遺伝子・同じ手順）
- 実行: `run_isp.py --config experiments/EXP-20260926-001/config.yaml --run tls_high`（次に `--run tls_low`）

## 5. コントロール
- ポジティブ: VCAM1, HAVCR1（delete で aPT → 健常側、既知の再確認）
- ネガティブ（null）: **発現頻度マッチのランダム 20 遺伝子を seed ごとに抽出**（`controls.negative_sampler: freq_matched`）。候補 = start 細胞のうち 300 細胞以上で入力に含まれる遺伝子（9,530 個）から対象・陽性・HK を除く。対象 5 遺伝子を順に参照し、発現細胞数の比が 1/1.5〜1.5 倍の候補から各 4 個（seed 0 の発現細胞数: 353〜6,478）
- ハウスキーピング: ACTB, GAPDH, PPIA（role = `housekeeping`。非特異的なシフトの参照で、null には入れない）
- 状態のコントロール: alt 方向・start 方向のシフトも細胞ごとに記録（「別の aPT 状態へ動いただけ」でないかの確認）

## 6. 評価指標と判定基準（事前登録。仮説からそのまま転記）
- 指標: goal_state_shift = PT_healthy セントロイドへのコサイン類似度の増加の細胞平均（`shift_goal_mean`, stats の `Shift_to_goal_end` と一致を確認済み）
- **主要**: TP53, CDKN1A, LTBR のそれぞれで、shift が同 seed のランダム 20 遺伝子の null の上位 5% に入り（経験的 p）、対象 5 遺伝子で BH-FDR < 0.05。3 seed のうち 2 以上で成立
- **文脈依存性（run 2）**: shift(TLS-high) / shift(TLS-low) ≥ 1.2 かつ donor 単位ブートストラップ 95% CI が 1 を含まない
- **参照**: SPP1 / LCN2 は VCAM1 / HAVCR1 と同程度（null の上位 10% 以内）
- **反証**: TP53, CDKN1A, LTBR が 3 seed すべてで null の 80 パーセンタイル以下、または向きが逆 → 主要予測は棄却。比の CI が 1 を含む → 文脈依存性は棄却
- **inconclusive**: 陽性対照（VCAM1, HAVCR1）が null の上位 10% に入らない → パイプライン・状態定義を見直す
- 計算: `validate_isp.py --exp experiments/EXP-20260926-001 --results results/tls_high --config experiments/EXP-20260926-001/config.yaml --alt-results results/tls_low`（seed ごとの経験的 p・null 内分位・BH、seed 一致数、陽性対照、比と donor ブートストラップ CI（2,000 回, 同一細胞は seed 間で平均してから donor を run ごとに復元抽出））

**⚠ G2 前に PI の判断が必要な点（設計上の問題。基準は変更していない）**
1. **経験的 p の下限**: null が 20 遺伝子だと p の最小値は 1/21 = 0.048。対象 5 遺伝子で BH をかけると、**5 遺伝子すべてが null の最大値を超えた場合にしか FDR < 0.05 にならない**（1 遺伝子だけなら q ≥ 0.24）。選択肢: (a) null を ≥100 遺伝子に増やす（run 1 で +8,000 細胞×遺伝子/seed... 計 +24,000、+4.1 h）、(b) null の正規近似 p を FDR に使う（validate_isp.py は `p_norm`/`fdr_norm` も出力済み）、(c) FDR をやめ「上位 5% ＋ 2/3 seed」だけにする。推奨は (b)＋上位 5% を併記（追加コストなし）
2. **非特異的な削除ドリフト**: phase 0c（陽性対照・HK・ランダム 10 遺伝子）では、**全遺伝子の delete が goal・start・alt の 3 セントロイドすべてへのコサインを同程度下げた**（shift_goal −1.0×10⁻⁴〜−4.6×10⁻⁴, 1 seed, n = 100 細胞/遺伝子）。VCAM1（−1.1×10⁻⁴）・HAVCR1（−1.4×10⁻⁴）は HK・ランダムと区別できなかった。このままでは陽性対照が上位 10% に入らず **inconclusive** になる可能性が高い。副次指標として「特異的シフト」= shift_goal − shift_start（または − shift_alt）を**本番前に**事前登録するか判断してほしい（cells.csv に必要な列はすべて保存される。主要基準は変えない）

## 7. 頑健性チェック
- 3 seed（seed ごとに遺伝子ごとの細胞抽出・ランダム 20 遺伝子・セントロイドの細胞が変わる）。run 間 Spearman（validate_isp.py 既存）
- 発現頻度と shift の相関（交絡チェック, 既存）
- donor 数: TLS-low は 8 donor のみ → run 2 は 1 遺伝子あたり最大 8×15 = 120 細胞、LCN2（TLS-low CKD aPT で 87 細胞）は 100 に届かない
- 精度: bf16 vs fp32 の matched pair（§8, phase 0c）

## 8. 計算量見積もり（phase 0 の実測, MPS）
**phase 0 の結果（`results-phase0/`）**
| 項目 | 結果 |
|---|---|
| a. 入力作成 + dry-run | 17,684 細胞（5 秒）、dry-run OK（`logs/EXP-20260926-001-phase0-{prepare,dryrun}.log`） |
| b. 事前チェック（CLS, LR, 5-fold） | start vs goal: donor 分割 **AUC 1.00 ± 0.00**（層化 1.00）→ 0.8 以上で ISP 可。start vs alt: donor 分割 **AUC 0.627 ± 0.206**（fold 0.34〜0.83, 層化 0.899 ± 0.025）→ 0.6 以上で run 2 可だが境界的（n = 500 / 500 / 405 細胞、14 / 22 / 8 donor。`results-phase0/precheck/precheck.json`） |
| c. bf16 vs fp32（10 遺伝子 × 100 細胞, 1 seed, セントロイド 200 細胞/状態） | **Spearman 0.879（基準 0.98 未達）**、上位 5 一致 5/5、符号一致 10/10（中央値超でも 5/5）、mean|Δ| 1.23×10⁻⁵ / mean|shift| 1.77×10⁻⁴ = 7.0%。細胞単位では Spearman 0.677（n = 1,000）。速度差は 1.33 倍 → PI 決定により **本番は bf16**（`results-phase0/precision_bf16_vs_fp32.json`） |
| d. 時間（モデルロード除く） | 1 細胞×1 遺伝子あたり **fp32 0.822 s / bf16 0.616 s**（モデル再ロードは cache_model で 0 s）。状態セントロイドは 0.244 s/細胞（fp32） |

**本番の見積もり**
| run | 遺伝子 | 1 遺伝子あたり細胞（実測の抽出数） | seed | 細胞×遺伝子 | ISP | セントロイド | 計 |
|---|---|---|---|---|---|---|---|
| 1 `tls_high` | 30（対象 5 + 陽性 2 + HK 3 + ランダム 20） | 100（全遺伝子で 100 / 12〜14 donor） | 3 | 9,000 | 2.05 h | 0.29 h | **2.34 h** |
| 2 `tls_low` | 30 | 平均 88（TP53 89, CDKN1A 69, LTBR 100, SPP1 100, LCN2 48, VCAM1 100, HAVCR1 62, HK 100, ランダム平均 92。donor は 6〜8） | 3 | 約 7,900 | 1.81 h | 0.29 h | **2.09 h** |
| 合計 | | | | 約 16,900 | | | **約 4.4 h（fp32 換算）→ bf16 で約 3.3 h（MPS）** |

- 摂動対象の細胞は延べ 16,900 細胞×遺伝子（ユニーク細胞は run 1 で最大 7,533、run 2 で 4,002 → **どちらも 1 万細胞未満**）
- bf16（PI 決定）で約 3.3 h（0.616/0.822 s で換算。dry-run の見積もり 5.1 h（fp32）との差は未確認）。いずれにせよ 2 h 超のため **G2 必須**
- 精度の補強策: 主要な結論（対象 5 遺伝子が null 上位 5% に入るか）について、境界付近（null 分位 90〜99%）の遺伝子のみ fp32 で再計算して判定が変わらないことを確認する（追加 ~0.3 h）
- セントロイドの抽出は run 1 と run 2 で同じ細胞（同一 seed）なので、実装上は再利用できれば −0.29 h

**G2 判定**: **必要**。基準のうち「GPU 2 時間超」に該当する（実測ベース 4.4 h、run 1 単独でも 2.34 h）。1 万細胞超・全遺伝子スクリーニングには該当しない。`run_isp.py --dry-run` の判定も `needs_G2_approval: true`

## 9. 想定される失敗と対策
| 失敗 | 対策 |
|---|---|
| 陽性対照が上位 10% に入らない（削除ドリフトが支配的, §6-2） | inconclusive として報告。事前登録した副次指標（採用された場合）を併記。状態定義・emb_layer の見直しは別実験 |
| TLS-high と TLS-low の分離が donor 依存（AUC 0.63、fold 0.34〜0.83） | run 2 の比は donor ブートストラップ CI を主に見る。CI が広ければ「判定不能」と書く（棄却とは区別） |
| 発現細胞が足りない遺伝子（run 2 の LCN2 は TLS-low で 87 細胞、`min_cells_per_gene: 30` 未満なら skip） | 実際の n を summary.csv / cells.csv に記録し、n の異なる遺伝子を同列に比較しない |
| null 20 遺伝子では FDR が原理的に出ない（§6-1） | G2 で方針を確定してから本番を回す |
| MPS で 2 ジョブ同時実行 | run 1 → run 2 を逐次実行する。夜間キューに入れる場合も直列 |
| 途中で落ちた場合 | 遺伝子ディレクトリ単位で再開できる（stats と cells.csv が揃った遺伝子は skip） |
| ambient RNA（IL7R/LTB/CCL19 の aPT 陽性） | 今回の対象遺伝子には含めていない。解釈時に言及 |

## 10. 承認
- G1（仮説採択）: 済（2026-09-26, PI）
- G2（計画承認）: 済（2026-09-26, PI）

### G2 で確定した事項（2026-09-26、本番データを見る前に登録）
1. **精度**: bf16（forward bf16 + コサイン fp32）。PI 決定。phase 0c で fp32 と Spearman 0.879（基準 0.98 未達）→ 限界として明記し、判定境界付近（null 分位 90〜99%）の遺伝子は fp32 で再計算して判定が変わらないか確認する
2. **FDR**: null 20 遺伝子の正規近似 p（`p_norm`）に BH（対象 5 遺伝子）をかける。経験的分位（null 上位 5%）も併記する（§6-1 の選択肢 b）
3. **指標**: `shift_goal` に加えて **`shift_goal_minus_start`（細胞ごとに goal 方向シフト − start 方向シフト）を主要指標として追加**（§6-2）。
   - 経緯: phase 0c で、全遺伝子の delete が goal/start/alt の 3 セントロイドへ同程度の負シフトを示し、陽性対照 VCAM1/HAVCR1 が HK・ランダムと区別できなかった（削除ドリフト）。差分はこの一様成分を差し引く。phase 0（本番とは別の抽出・1 seed・10 遺伝子）を見て追加したことをレポートに明記する
   - 判定規則は §6 と同じ（null 上位 5%・FDR < 0.05・3 seed 中 2 以上・陽性対照は上位 10%）を両指標に適用し、両指標で結論が分かれた場合はその旨を報告する
4. 実行: run 1（tls_high）→ run 2（tls_low）を逐次。見積もり約 3.3〜3.9 h（bf16, MPS）

Approved-by: PetadimensionLab (PI), 2026-09-26 — G2 承認（上記 1〜4 を含む）

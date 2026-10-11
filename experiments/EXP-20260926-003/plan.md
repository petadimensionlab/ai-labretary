# EXP-20260926-003 実験計画 — fine-tuned Geneformer CellClassifier による in silico perturbation（方法検証, 対照先行）

- 仮説: [HYP-20260926-02](../../hypotheses/HYP-20260926-02.md)（方法仮説, **Status: proposed。G1 未承認**）
- 起点: PI 指示（2026-09-26）「A の後に、そもそもファインチューニングをしていないと思うので、ファインチューニングをしたトレインドのモデルを使ったインシリコパーターベーションの実験を試みてください」。追加指示（同日）: **fine-tuning だけをリモートの NVIDIA GPU（Tailscale 経由）で行い、ISP は Mac（MPS）**
- 作成: Labretary（experiment-planner）/ 2026-09-26
- 状態: **計画のみ**。phase 0（CPU だけの準備・検証）まで実施。GPU/MPS は未使用（EXP-20260926-002 案 A が MPS で実行中のため）
- 事前登録の時点: 本計画の基準（§6）は EXP-002 案 A の結果を見る前（2026-09-26 15:50 頃, 案 A 実行中）に記入した

## 1. 問い
1. 腎 PT の障害軸で fine-tuning した CellClassifier の埋め込みでは、1 遺伝子 delete の「削除ドリフト」が Pretrained より小さくなるか
2. その埋め込みで、陽性対照（VCAM1, HAVCR1）が HK・陰性対照パネル・ランダム null から、EXP-002 と同じ事前登録基準（主指標 a2 の C1〜C3, C3-panel）で分離するか
3. 同じ細胞・同じ遺伝子・同じ goal セントロイド定義の Pretrained と比べて改善するか

方法の検証であり、生物学的な仮説の検証ではない。主要対象（TP53/CDKN1A/LTBR）と参照（SPP1/LCN2）は、対照が合格し PI が承認するまで回さない（盲検, `controls_only: true`）。

## 2. データ
| 項目 | 値 |
|---|---|
| データセット | A2 KPMP v1.5 scRNA（Lake 2023, PMID 37468583）。EXP-001/002 と同じ入力（`work/input.h5ad` → EXP-001 への symlink, 読み取りのみ） |
| tokenized | `work/tokenized.dataset` = EXP-001 の tokenized のデータ部分だけを APFS clone（sha256 `ba3a14c0…` 一致を確認）。他実験の dataset を直接指さない（guardrails） |
| 細胞数（start / goal / alt） | 7,533（14 donor）/ 6,149（22 donor）/ 4,002（8 donor）。3 状態は donor が重ならない（計 44 donor） |
| 状態ラベル列 | `isp_state`（**推定ラベル**。TLS 群は donor 単位 pseudo-bulk の TLS スコア三分位で、組織学的な TLS ではない） |
| フィルタ | EXP-001 と同じ |
| 分割（FT する場合） | **ドナー単位**（§2.2）。held-out donor の細胞はリモートへ送らない |

### 2.1 分類課題の選択: 2 クラス（aPT_CKD vs PT_healthy）
- **採用**: `aPT_CKD`（学習は TLS-low の CKD donor の aPT のみ）vs `PT_healthy` の 2 クラス。TLS 群は ISP の start/alt の定義にだけ使う
- 3 クラス（TLS-high / TLS-low / 健常）を採らない理由:
  1. TLS-high vs TLS-low は凍結埋め込みで分離が弱い（事前チェック AUC 0.627 ± 0.206, セントロイド間距離は start–goal の 1/12）
  2. TLS 群は donor 単位のラベル（14 + 8 donor）なので、donor 分割で学習すると donor の個性を学ぶだけになりやすく、held-out で評価できる donor が数人しか残らない
  3. ISP の読み出しは start→goal（aPT→健常）方向なので、その軸を学習させるのが目的に合う
- 広い PT 障害軸（他データの追加）は、同じ前処理・同じ donor 規約のデータがないため今回は採らない
- 循環性: goal 方向の移動を、start と goal を分けるよう学習したモデルで読むことになる（Theodoris 2023 の ISP も同じ構造）。対策は HYP-20260926-02 §「循環性」: (i) 埋め込みの移動（主）と分類確率の変化（副）を両方出す、(ii) 特異性は HK・パネル・null との比較で判定、(iii) 同じ細胞での Pretrained 比較、(iv) ISP の細胞は FT に一切使わない donor、(v) マーカー除去効果の探索枠（§4.5, PI 判断）

### 2.2 donor 分割（事前固定。config-finetune.yaml `split`, CPU で集計した実数）
処理バッチ（`Sample####-EO1/2/3` は性別・年齢が異なる別人だが同じバッチ）は同じ分割に置く。ISP の start 状態（TLS-high）は **14 donor 全員を held-out**。

| 分割 | 用途 | aPT | 健常 PT | 細胞（上限適用後） |
|---|---|---|---|---|
| train | 学習 | TLS-low 6 donor（31-10000, 31-10013, 29-10010, 28-11515, 29-10406, 28-10051）、donor ≤ 500 → **1,565** | 8 donor（164-10, 165-9, 164-3, Sample1157-EO1/2/3, PRE018-1, PRE062-1）、donor ≤ 250 → **1,571** | **3,136 / 14 donor**（長さ中央値 1,132, 平均 1,395 トークン, 計 4.37 M トークン） |
| eval | early stopping・モデル選択 | TLS-low 2 donor（27-10053, 29-10386）≤ 60 → 120 | 5 donor（PRE019, PRE055-1, Sample1158-EO1/2/3）≤ 60 → 254 | 374 / 7 donor |
| test（held-out） | ゲート評価・ISP | TLS-high **14 donor 全員**（ゲート用は ≤ 100 → 1,271。ISP は全 7,533 細胞から抽出） | 9 donor（PRE038, 164-6, PRE027, Sample1162-EO1/2, PRE98sc, Sample1153-EO1/2/3）。ゲート用 ≤ 100 → 811、利用可能 3,106 細胞 | ゲート用 2,082 / 23 donor |

- ISP の細胞数: 案 A と同じ規則（200 細胞/遺伝子, donor ≤ 25）で、start は全員 held-out なので**案 A と同一の細胞**になる（CPU で確認: VCAM1・HAVCR1・GAPDH・PHPT1 の seed 3 で cell_uid が 200/200 一致, 13〜14 donor）
- goal セントロイド: held-out の健常 9 donor から 500 細胞（donor ≤ 60。8 donor × 60 + 31 = 511 細胞が利用可能）
- alt（TLS-low）の 8 donor は全員 FT の train/eval に使われている → alt 方向のシフトは参考値（a2 は goal 方向だけを使う）
- bundle（リモートへ送る唯一のディレクトリ `transfer/bundle/`, 19.2 MiB）: train/eval の細胞だけ。prepare が書いた後に読み直して、held-out donor・start 状態・test と同じ細胞が 1 つもないことを assert。manifest には held-out の件数だけを書き ID は書かない
- 注意（交絡）: 健常 PT の細胞は短い（長さ中央値 ~1,000 トークン）、aPT は donor により 710〜2,540。分類器が細胞長を手がかりにする恐れ → ゲートで「クラス内の P と長さの Spearman」を報告（合否には使わない）、ISP では C5 で確認

## 3. モデル
- 事前学習: Geneformer-V2-316M（18 層, hidden 1152, vocab 20,275, max 4,096 トークン）
- fine-tuning: `BertForSequenceClassification`（2 クラス, 0 = aPT_CKD, 1 = PT_healthy）、下位 6 層凍結（埋め込みは学習, Geneformer Classifier と同じ）、実効バッチ 16、lr 5e-5 linear・warmup 10%、weight decay 0.001、最大 3 epoch、20 optimizer step ごとに eval（eval_loss）、early stopping patience 4、最良 checkpoint を採用、gradient checkpointing、bf16（PI 標準）、seed 42
- 実行場所: **リモート NVIDIA GPU**（remote/README.md）。VRAM に合わせて micro-batch × 勾配累積だけを変える profile（vram16: 2×8, vram24: 4×4, vram48: 8×2, vram80: 16×1）。bf16 非対応の GPU では停止（`compute.bf16_unsupported: error`。fp16/fp32 は PI 判断）
- ISP: FT 済みモデル（`models/ft-main/final`, `model.type: CellClassifier`, `num_classes: 2`）、`emb_mode: cls`、emb_layer −1、MPS・bf16 forward + fp32 コサイン、`cache_model: true`、`nproc: 1`（guardrails）
- 比較: Pretrained（EXP-002 案 A と同じ設定）

## 4. Perturbation 設計
### 4.1 phase 0（実施済み, CPU のみ）
- donor 分割の検証と bundle 作成（`finetune_classifier.py --stage prepare`）、train の dry-run（既定・vram24・mps_fallback）
- CPU smoke: 小さなランダム BERT で train → 途中停止 → checkpoint から再開 → probe → evaluate を通した。**torch 2.14 + transformers 4.46.3 では checkpoint からの再開が `torch.load(weights_only=True)` の UnpicklingError で止まる**ことを発見し修正（numpy の型だけを許可。guardrails に記録）
- run_isp.py の CellClassifier 経路を CPU smoke（小さな FT モデル, 3 遺伝子 × 6 細胞, forward_batch_size 4）: 取り出した logits が直接計算と一致（最大差 5.6×10⁻⁹、元と摂動後の差 ~10⁻⁶ より十分小さい）、goal セントロイドは held-out donor だけ
- `make demo` 通過、既存 config（EXP-20260925-001, EXP-20260926-001/002 と 002 の案 B・phase 0）の dry-run 通過

### 4.2 phase 1: fine-tuning（リモート GPU, G1・G2 承認後）
1. remote/README.md 手順 0〜4（環境・転送・sha256 照合・dry-run）
2. probe（20 optimizer step）で s/step とピークメモリを実測 → profile を確定
3. 本番（`run_remote.sh train <profile>`）。中断時は同じコマンドで再開
4. 結果を `models/ft-main/` へ戻し、`final/SHA256SUMS` を照合
- ISP（Mac）とは別の機械なので、G2 承認後なら EXP-002 案 A の実行中でも開始できる

### 4.3 phase 2: ゲート評価（Mac MPS, §6.1）
- `--stage evaluate`: held-out 2,082 細胞で FT の P(PT_healthy) と、凍結 Pretrained の CLS 埋め込み（emb_layer −1）＋標準化ロジスティック回帰（C = 1, train と同じ 3,136 細胞で学習）を比較
- NO-GO なら phase 3 以降は行わない

### 4.4 phase 3–4: 対照の ISP（ゲート GO 後）
- **phase 3（FT）**: `config.yaml`。案 A と同一の設計: delete、対照のみ 34 遺伝子（陽性 2 + HK 3 + パネル 7 + RP 2 + 発現頻度マッチ null 20）× 200 細胞（donor ≤ 25）× seed 3/4/5。null は案 A と同じ規則・同じ seed → 同じ 20 遺伝子。`capture.embeddings`（CLS）と新しい `capture.classifier_logits`（分類ヘッドの logits）を保存
- **phase 4（Pretrained の対応比較）**: 案 A は同じ細胞・同じ遺伝子で元/摂動後の CLS 埋め込みを保存している（capture on）→ goal セントロイドだけを held-out 健常 9 donor の細胞で作り直し（Pretrained の EmbExtractor, seed 3–5 × 500 細胞 ≈ 0.1 h MPS）、shift_goal を埋め込みから再計算（CPU。EXP-002 で埋め込みからの再計算と perturber の値の差は最大 3.1×10⁻⁷）。再計算スクリプトは phase 4 の前に作る
  - 代替: 案 A が未完了・capture 欠落なら `config-pretrained-matched.yaml`（同じ設計を Pretrained で再実行, 4.3 h, 別途 G2）
  - 感度分析: 案 A 本来の goal セントロイド（全 22 donor）での Pretrained の結果も並べる（セントロイドの定義の影響）

### 4.5 探索（任意・PI 判断。合否に使わない）
- マーカー除去効果の確認: start で陽性対照と同程度に上昇する（aPT vs 健常の発現差と発現細胞数でマッチした）非ドライバー 10 遺伝子を FT と Pretrained の両方で回す（各 10 × 200 × 3 seed ≈ 1.1 h, 計 2.3 h）。陽性対照がこれらを超えなければ、C1 の合格は「分類器がマーカーを使っている」ことの反映にすぎない
- emb_layer 0（最終層）での FT の対照（1 seed, 34 × 200 ≈ 1.3 h）。Geneformer は訓練目的に近い問いには最終層を推奨

### 4.6 主要対象（対照合格・PI 承認後）
- 5 遺伝子 × 200 細胞 × seed 3–5 を FT モデルで追加（≈ 0.6 h）。判定は EXP-002 §4.3 と同じ規則を a2 に適用

## 5. コントロール
- ポジティブ: VCAM1、HAVCR1（EXP-001/002 と同じ）
- ネガティブ: 発現頻度マッチのランダム null 20 遺伝子（seed ごと, 案 A と同一）
- ハウスキーピング: ACTB, GAPDH, PPIA（**C2・C3 に含める**, PI 方針）
- 陰性対照パネル: CERS4, TATDN1, IFT25, VAPA, PHPT1, ANAPC16, EEF1D（C3-panel。選定は発現データのみ, `outputs/literature/2026-09-26-negative-controls.md`）
- 探索枠: RPS13, RPL32（報告のみ）
- モデルの対照: 同じ細胞・同じ遺伝子・同じ goal セントロイドの Pretrained（phase 4）
- 分類器の対照: 凍結埋め込み＋ロジスティック回帰（ゲート）

## 6. 評価指標と判定基準（事前登録, 2026-09-26, 案 A の結果を見る前）
### 6.1 ゲート（phase 2）: FT モデルを ISP に使ってよいか
held-out 23 donor（TLS-high aPT 1,271 細胞 + 健常 811 細胞）で、閾値 0.5:
| # | 基準 | 合格ライン |
|---|---|---|
| G-a | AUROC | FT ≥ ベースライン（ベースライン ≥ 0.99 なら FT ≥ ベースライン − 0.01） |
| G-b | macro-F1 | 同上 |
| G-c | donor 別の最小正解率（≥ 20 細胞の donor） | FT ≥ 0.80 |
- G-a〜G-c をすべて満たせば GO。1 つでも不合格なら FT モデルでの ISP は行わない（HYP-20260926-02 は inconclusive）
- ベースラインは天井（事前チェックで start vs goal の donor 分割 AUC 1.00）に近いと予想されるため、「上回る」ではなく非劣性（margin 0.01）で判定する。PI 指示の「ベースラインに勝つ」を厳密に適用すると天井では判定不能になるため（§10 で PI に確認）
- 報告のみ: log-loss、donor ブートストラップ 95% CI（1,000 回）、クラス内の P と細胞長の Spearman、eval（TLS-low）での値
- 学習が max_steps で打ち切られた場合は `training_summary.json` の `truncated_by_max_steps` とともに報告（../Geneformer の 31% epoch 打ち切りの教訓）

### 6.2 ISP（phase 3–4）: EXP-002 §6.2 / §6.4 と同一
- **主指標 a2**（長さ 5 分位ビン［null 細胞で境界］ごとの null 平均 shift_goal を細胞ごとに引き、遺伝子で平均。null 遺伝子自身は leave-one-gene-out）
- C1: VCAM1・HAVCR1 がそれぞれ同 seed の null の上位 10%（分位 ≥ 0.90）に 3 seed 中 2 以上
- C2: 陽性 2 × HK 3 × 3 seed = 18 組のうち陽性 > HK が 15 組以上
- C3: ACTB・GAPDH・PPIA がすべて null 分位 0.10〜0.90 に 3 seed 中 2 以上
- C3-panel: パネル 7 遺伝子中 5 以上が 3 seed 中 2 以上で null 分位 0.10〜0.90
- C4（参考）: split-half（donor 2 群, 20 回）の遺伝子順位 Spearman ≥ 0.7、seed 間（固定対照）≥ 0.7
- C5（参考）: 遺伝子×seed 単位の指標と平均細胞長・発現細胞数の |ρ| ≤ 0.3、細胞単位の長さとの |ρ| 中央値 ≤ 0.2
- 判定: **C1・C2・C3・C3-panel をすべて満たせば HYP-20260926-02 の P1 を支持**。1 つでも不合格なら棄却（§HYP 反証条件）
- 副指標（比較のみ）: raw, d, d2, c, V2, V2L（EXP-002 と同じ）、**dP_goal_a2** = 細胞ごとの ΔP(PT_healthy) = softmax(logits_pert)[1] − softmax(logits_orig)[1] に a2 と同じ長さビン補正。C1〜C3・C3-panel を同様に計算して報告
- P2（ドリフト）: 補正前 shift_goal の null 遺伝子×seed の負の割合、null 平均 shift_goal / (1 − cos(c_start, c_goal)) を FT と Pretrained（phase 4, 同じ細胞・同じ goal セントロイド）で比較
- P3: C2 の勝ち数 FT ≥ Pretrained
- 解析コード: EXP-002 の `analysis/evaluate_controls.py`（`reanalysis_cosine.py` の a2 実装）を EXP-003 のパスに向けて使う。dP_goal_a2 の追加以外は変更しない（phase 3 の前にコードを固定し sha256 を記録）
- 盲検: 主要対象・参照は回さない（`controls_only: true`）。解析スクリプトでも読み込み直後に除外（assert）

## 7. 頑健性チェック
- seed 3–5（案 A と同じ細胞・同じ null。細胞は同じなので、FT と Pretrained の差は細胞の抽出ではなくモデルの差）
- split-half（donor 単位）と seed 間の順位一致（C4）
- goal セントロイドの定義（held-out 9 donor vs 全 22 donor）の感度（Pretrained 側, phase 4）
- 長さビン数 3 / 5 / 10 の感度
- 精度: bf16（PI 標準）。a2 の判定が境界（分位 0.85〜0.95）の対照は fp32 で再計算（最大 5 遺伝子 × 200 細胞 × 1 seed ≈ 0.25 h, MPS）
- 学習の再現性: 学習 seed・データ順序は固定だが GPU の非決定性は残る。env.json（GPU・driver・CUDA・各版・geneformer の .py digest）で記録

## 8. 計算量見積もり
**リモート GPU の時間も G2 の対象**（CLAUDE.md §3 の「GPU 2 時間超」）。

### 8.1 fine-tuning（phase 1）
学習 1 epoch = 3,136 細胞・4.37 M トークン。演算量の目安: 1 トークンあたり順伝播 ~0.69 GFLOP（非埋め込み 287 M パラメータ + 長さ ~1,400 の注意）、逆伝播（上位 12 層は重み勾配あり、埋め込みを学習するので下位 6 層も入力勾配が通る）~1.15、gradient checkpointing の再計算 ~0.69 → **~2.5 GFLOP/トークン、~11 PFLOP/epoch**。MPS の実測（../Geneformer: batch 2・2,767 トークンで ~10 s/step → 実効 ~1.4 TFLOPS）から逆算した MPS 1 epoch 2.2 h は dry-run の下限と一致。

| 実行場所 | 仮定（実効性能） | 1 epoch | 最大 3 epoch + eval | 備考 |
|---|---|---|---|---|
| RTX 4090（24 GB, vram24） | bf16 165 TFLOPS の 10〜40% | 3〜18 min | **0.2〜0.9 h** | |
| A100 / H100 級（vram48/80） | 312 / ~990 TFLOPS の 10〜40% | 1〜6 min | **0.1〜0.4 h** | |
| 未知の GPU（aarch64 統合メモリ機など） | 30〜100 TFLOPS 相当の 10〜40% | 5〜60 min | 0.3〜3 h | probe で実測して確定 |
| 環境構築・転送・probe・照合 | — | — | ~0.3 h（うち GPU ~0.1 h） | |
| **MPS フォールバック（本構成）** | 実測 10 s/step | 2.2〜4.4 h | 6.6〜13 h | 非推奨 |
| **MPS フォールバック（profile mps_fallback）** | 2,091 細胞、下位 12 層＋埋め込み凍結、1 epoch | — | **0.8〜1.6 h** | 逆伝播が上位 6 層だけになり ~45% 減。精度低下の恐れ |

- 不確実性: GPU が未確定（Tailscale 上の候補機は 2026-09-26 時点でオフライン）。実効性能の幅（10〜40%）が最大の不確実性。**probe（20 step, 数分）の実測で本番の時間を確定してから開始**し、見積もりの 2 倍を超えたら中断して報告
- 早期終了は eval_loss の改善が 80 step（~0.4 epoch）止まった時点。天井に近い課題なので 1 epoch 前後で止まる可能性が高い

### 8.2 Mac（MPS）の時間
| phase | 内容 | 計算 | 時間 |
|---|---|---|---|
| 2 | ゲート評価 | FT 予測 2,082 細胞 + Pretrained 埋め込み 5,218 細胞（~0.2 s/細胞） | ~0.5 h |
| 3 | FT の対照 ISP | 34 遺伝子 × 200 × 3 seed = 20,400 細胞×遺伝子 × 0.684 s（案 A 実測, capture 込み）+ セントロイド | **~4.3 h**（dry-run 4.30 h） |
| 4 | Pretrained の goal セントロイドの再抽出 + 再計算 | 1,500 細胞 + CPU | ~0.1 h |
| 4 代替 | Pretrained で再実行（案 A が使えない場合） | 20,400 | 4.3 h |
| 4.5 任意 | DE マッチ探索（FT + Pretrained） | 2 × 6,000 | 2.3 h |
| 4.5 任意 | emb_layer 0（1 seed） | 6,800 | 1.3 h |
| 4.6 | 主要対象（対照合格・PI 承認後） | 3,000 | 0.6 h |

### 8.3 合計（G2 の申請対象）
- **必須: リモート GPU 0.1〜0.9 h（未知の GPU なら最大 ~3 h）+ Mac MPS 約 4.9 h（ゲート 0.5 + ISP 4.3 + セントロイド 0.1）= 約 5.0〜5.8 h**（GPU によっては最大 ~8 h）→ **G2 必要**（ISP 単独で 2 h 超）
- リモートが使えず MPS フォールバックの場合: 0.8〜1.6 h + 4.9 h = 約 5.7〜6.5 h
- 1 万細胞超・全遺伝子スクリーニングには該当しない（ユニーク細胞は ISP 最大 7,533、学習 3,136）
- MPS は 1 ジョブで飽和するため Mac 側は逐次（EXP-002 案 A の完了後）。リモート FT は案 A と並行可

## 9. 想定される失敗と対策
| 失敗 | 対策 |
|---|---|
| ゲート不合格（FT がベースラインに劣る） | ISP は行わず inconclusive として報告。学習率・epoch の再検討は PI 判断（結果を見てからの調整は探索扱い） |
| 打ち切り学習でベースライン以下（../Geneformer の 31% epoch の例） | max_steps は既定で使わない。打ち切った場合は truncated と明記 |
| 分類器が細胞長を手がかりにする | ゲートでクラス内の P–長さ相関を報告、ISP で C5 と a2（長さビン補正）で確認 |
| FT で alt（TLS-low）が学習データ | alt 方向は参考値。a2 は goal のみ |
| 陽性対照の合格がマーカー除去効果 | §4.5 の DE マッチ探索（PI 判断）。結論の書き方を「分類器が依存する遺伝子」にとどめる |
| リモートの接続断・OOM | step checkpoint（eval ごと）から同じコマンドで再開（CPU smoke で確認済み）。OOM は 1 段小さい profile（実効バッチ同じ） |
| torch 2.14 の weights_only で再開不能 | 修正済み（numpy の型だけ許可） |
| bf16 非対応 GPU | 停止して PI 判断（fp16 / fp32） |
| torch 2.14.0 の wheel がリモートの arch/CUDA にない | 版を変えず報告し、逸脱として記録してから進める |
| 転送の破損・取り違え | bundle・重み・返送モデルの sha256 を両側で照合（verify_sha256.py, run_remote.sh） |
| held-out データの流出 | bundle は train/eval のみ（prepare で assert）。終了後にリモートの bundle と out を削除 |
| 案 A が未完了・埋め込み欠落 | phase 4 は `config-pretrained-matched.yaml` の再実行（4.3 h, 別途 G2） |
| MPS 使用後の fork 停止 | `nproc: 1`（guardrails） |

## 10. 承認
- G1（HYP-20260926-02 の採択）: **未**（PI）
- G2: **済（PI, 2026-09-26「結果が出たら、ファインチューニングしたモデルで perturbation を行う」）**。ゲート GO を条件に phase 3（4.3 h）・phase 4（0.1 h）を実行。申請内容: リモート GPU の fine-tuning（probe 後に本番, 見積もり 0.1〜0.9 h・未知の GPU なら最大 ~3 h）＋ Mac のゲート評価 0.5 h・FT の対照 ISP 4.3 h・Pretrained セントロイド 0.1 h（計 約 5.0〜5.8 h）。phase 3（ISP）はゲート GO が条件
- PI に判断してほしい点:
  1. HYP-20260926-02 の採択（G1）と 2 クラス課題（aPT_CKD vs PT_healthy, TLS 群は ISP のみ）でよいか
  2. ゲートの基準: ベースラインが天井（~1.0）の指標は「0.01 以内の非劣性」で GO としてよいか（厳密に「上回る」を求めると判定できない）
  3. リモート GPU の機種・VRAM・OS（未確定）。torch 2.14.0 の wheel がない場合の扱い、bf16 非対応だった場合に fp16 か fp32 か
  4. G2（上記の合計）。リモートが当面使えない場合に MPS フォールバック（縮小構成 0.8〜1.6 h）へ切り替えるか
  5. 任意の探索（§4.5: DE マッチ 2.3 h、emb_layer 0 1.3 h）を加えるか

### fine-tuning の実行記録（2026-09-26）
- PI の実行指示（「実行して下さい」）に基づき、fine-tuning（リモート GPU, 見積もり ≤ 0.8 h）を実行。ゲート評価と ISP（Mac, 約 4.9 h）は G1（HYP-20260926-02）・G2 をゲート通過後に改めて確認する
- 非劣性基準（0.01 以内）と任意の探索（なし）は Labretary 推奨どおり（PI 未回答のため既定値）
- リモート: spark-cf5b（NVIDIA GB10, aarch64, 128 GB 統合メモリ, driver 595.58.03）。当初の coi-next1@100.94.57.17 は TCP 22 に到達できず断念
- 環境の逸脱: リモートの既存 venv（~/workspace/Geneformer/.venv）を変更せずに使用 → torch 2.12.0+cu130 / transformers 4.57.6 / datasets 4.8.5 / accelerate 1.13.0（Mac: 2.14.0 / 4.46.3 / 4.0.0 / 1.15.0）。geneformer は Mac のパッチ済みパッケージを ~/labretary-ft/geneformer_pkg に置き PYTHONPATH で使用。モデル重みは Mac と sha256 一致
- probe: 4.648 s/optimizer step（実効バッチ 16, vram80 = 16×1）、ピーク 12.15 GB / 121.7 GB、196 step/epoch → 1 epoch ≈ 0.25 h。本番開始 17:00（tmux）

### 実行先の変更（2026-09-26, PI 方針「重い計算はリモートの DGX Spark 系」）
- ゲート: **GO**（held-out 23 donor, 2,082 細胞。FT AUROC 1.0000 / macro-F1 0.9980 / 最小 donor 正解率 0.99、baseline AUROC 1.0000 / macro-F1 0.9930）
- phase 3（FT の対照 ISP）は Mac ではなく spark-cf5b（CUDA, bf16）で実行。config は同一（LABRETARY_DEVICE=cuda で実行先だけ変更, device.json に記録）
- **Pretrained の対応比較も同じ spark-cf5b で `config-pretrained-matched.yaml` を実行**（PI 承認 2026-09-26）。デバイスの違い（CUDA vs MPS）を比較から除くため。phase 4（Mac で案 A を再計算）はデバイス間の参考値として残す
- 逸脱: リモートは datasets 4.8.5 のため、nproc ≤ 1 の map/filter を単一プロセスに強制（run_isp.py `single_process_datasets`）

Approved-by: PetadimensionLab (PI), 2026-09-26 — G1（HYP-20260926-02）・G2（fine-tuning 実施済み、ゲート GO を条件に phase 3–4）

## 11. コード変更
- **新規 `scripts/geneformer/finetune_classifier.py`**（config-finetune.yaml 駆動）: stage `prepare`（donor 分割の検証・bundle 作成・sha256・漏洩 assert）/ `train`（CUDA・MPS・CPU、bf16 対応確認、profile、gradient checkpointing・勾配累積、step checkpoint から自動再開、config sha256 の一致確認、env.json・run_record.json・train_log.jsonl・training_summary.json・final/SHA256SUMS、`--probe`）/ `evaluate`（held-out でのゲート判定）/ `smoke`（CPU・ランダム小モデルで全経路）。geneformer は collator と token dictionary のファイルだけを直接読む（__init__ の重い依存を避け、リモートでは `--no-deps` で入れられる）
- **`scripts/geneformer/run_isp.py`**（既定 off・後方互換）: `state_embs.include_donors`（セントロイドを指定 donor だけで作る。ファイル名に条件のハッシュ）、`capture.classifier_logits`（CellClassifier の分類ヘッドの logits を元/摂動後で保存。`in_silico_perturber.get_embs` の呼び出し順で区別）、CellClassifier の `num_classes` と config.json のクラス数の照合
- `remote/`: README.md（手順・VRAM 別の profile 表・rsync のひな形。ホスト名・資格情報はプレースホルダ）、run_remote.sh（照合 → probe/train → 返送前の照合）、verify_sha256.py（標準ライブラリのみ）
- 確認: prepare・train dry-run・CPU smoke（train/再開/probe/evaluate）・run_isp の CellClassifier CPU smoke（logits 一致）・`make demo`・既存 config の dry-run

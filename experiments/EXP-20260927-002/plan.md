# EXP-20260927-002 実験計画 — HYP-20260925-02 の事前チェック（状態の識別可能性）

- 仮説: [HYP-20260925-02](../../hypotheses/HYP-20260925-02.md)（Status: **proposed**。本 EXP は **G1 採択の前**に行う事前チェック）
- 作成: Labretary（experiment-planner） / 2026-09-27
- 種別: **事前チェックのみ**。摂動（in silico perturbation, ISP）は本 EXP では行わない
- 事前登録: 判定基準・細胞の選び方・解析手順は、埋め込みを取る前に本ファイルで固定する。結果を見てから変えない

---

## 1. 問い

HYP-20260925-02 の goal（TLS オーガナイザー様線維芽細胞, `FIB_org`）と start（TLS 高 donor の 3 遺伝子陰性線維芽細胞, `FIB_neg_TLShigh`）は、**Geneformer の細胞埋め込み上で区別できるか**。

区別できなければ、ISP の目的変数である「goal セントロイドへのシフト」は解釈できない（セントロイドが start と実質同じ場所にあるため）。したがって本チェックは HYP-20260925-02 の ISP を実行してよいかの前提条件であり、仮説自体の検証ではない。

補助的に、**その識別が「細胞の検出された遺伝子数と総 UMI 数」だけで説明されないか**を同時に確かめる。理由は §2.3。

## 2. データ

### 2.1 データセットと細胞

| 項目 | 値 |
|---|---|
| データセット | A1 `data/raw/kidney_tls/A1_KPMP_v1.5_snRNA_Lake2023.h5ad`（KPMP v1.5 snRNA, Lake 2023）**読み取り専用** |
| 細胞種の絞り込み | `subclass.l2` ∈ {FIB, aFIB, dFIB, M-FIB, dM-FIB, MYOF, cycMYOF}（著者ラベル）= 21,769 細胞 |
| 状態ラベル列 | `fib_tls_state`（**推定ラベル**。細胞種は著者の `subclass.l2`、TLS 状態は本解析の遺伝子検出による推定） |
| donor 列 | `donor_id`（80 donor）。TLS 群は `outputs/literature/2026-09-25-tls-feasibility/donor_A1.csv` の `TLS_group`（donor 単位の推定。ここで再定義しない） |
| 作業コピー | `work/input_fib.h5ad` →（ラベル付与）→ `work/input.h5ad`。元データは変更しない |

### 2.2 goal ラベルの定義（推定ラベルの作り方と閾値）

すべて **生カウント**（`raw.X`）に対して判定する。閾値は HYP-20260925-02 の定義をそのまま使い、本 EXP で変えない。

- 対象 3 遺伝子: CCL19 `ENSG00000172724` / CCL21 `ENSG00000137077` / CXCL13 `ENSG00000156234`
- **goal `FIB_org`**: 3 遺伝子の**いずれか 1 つでも 1 カウント以上**（`max(counts) >= 1`）
- **start `FIB_neg_TLShigh`**: 3 遺伝子が**すべて 0 カウント**、かつ donor の `TLS_group == high`、かつ `disease ∈ {acute kidney injury, chronic kidney disease}`
- **alt `FIB_neg_TLSlow`**: 3 遺伝子がすべて 0 カウント、かつ donor の `TLS_group == low`
- どれにも当てはまらない細胞は除く（`drop_unlabeled: true`）

**事前に確認した細胞数**（本日 2026-09-27、`data/raw` を読み取りのみで集計。`work/label_counts.json`, `work/fib_labels_expected.csv`）:

| 状態 | 細胞数 | donor 数 | 疾患内訳 | 検出遺伝子数の中央値 |
|---|---|---|---|---|
| goal `FIB_org` | **777** | 65 | AKI 359 / CKD 218 / normal 200 | 2,064 |
| start `FIB_neg_TLShigh` | 5,748 | 22 | CKD・AKI のみ | 1,547 |
| alt `FIB_neg_TLSlow` | 3,768 | 25 | normal 中心 | 1,858 |
| （参考）線維芽細胞系 全体 | 21,769 | 80 | — | 1,609 |

- goal 777 細胞の内訳: CCL19⁺ 252 / CCL21⁺ 411 / CXCL13⁺ 187、細胞種は FIB 474 / M-FIB 155 / aFIB 66 / MYOF 58 / dFIB 21 / cycMYOF 3、TLS 群は high 562 / mid 141 / low 69 / 除外 donor 5
- 仮説ファイルの「AKI 349」は本日の集計では **359** だった（他の数値は一致）。集計手順の差と考えられる。ラベル定義は変えず、実行時の `work/input.summary.json` の値を正とする

### 2.3 先に洗い出した交絡（結果を見る前に把握したもの）

1. **細胞の検出された遺伝子数と総 UMI 数（最大の交絡）**: goal は検出遺伝子数の中央値 2,064、start は 1,547。3 遺伝子のうち 1 つでも 1 カウント拾えるかは、そもそも検出された遺伝子数が多い核ほど起きやすい。
   **本日、埋め込みを使わずに検出された遺伝子数と総 UMI 数だけで識別を試したところ、donor 分割 5-fold の AUC = 0.676（fold: 0.67 / 0.67 / 0.55 / 0.72 / 0.77）だった**（`work/depth_baseline_pre.json`）。つまり「AUC ≥ 0.65」は **検出された遺伝子数と総 UMI 数だけで偶然クリアしてしまう**。
   → **対策**: 主判定は **検出された遺伝子数の分布をそろえた細胞集合**で行う（§4.2）。log(検出遺伝子数)の最近傍 1:1 マッチ（許容 0.10）で goal 557・start 557 細胞になり、検出された遺伝子数と総 UMI 数だけの AUC は 0.526 に落ちることを確認済み。
2. **donor の交絡**: goal は 65 donor に広く散らばり、start は 22 donor に集中する。21 donor は両方に細胞を出す。→ 分割は必ず **donor 単位**（§4.3）。さらに donor ラベルを保存した並べ替え検定で「donor 構成では説明できない分離か」を見る（§4.5）。
3. **疾患の交絡**: goal には normal donor 由来が 200 細胞入るが、start は AKI/CKD のみ。→ 感度分析で goal を AKI/CKD に限る。
4. **部位の交絡**: goal は皮質 75%、start は皮質 86%、alt は皮質 60%。→ 感度分析で皮質（`region == C`）のみ。
5. **細胞種構成の交絡**: goal には M-FIB（髄質）・MYOF が混ざる。→ 感度分析で `subclass_l2 == FIB` のみ。
6. **循環性**: goal ラベルは CCL19/CCL21/CXCL13 の検出そのもの。埋め込みが「その 3 遺伝子が入力トークン列に含まれるか」を見ているだけでも AUC は上がる。**この事前チェックは「ラベルが埋め込み上で分離すること」しか示さない**。生物学的な意味の主張には使わない（限界として報告に必ず書く）。

## 3. モデル

| 項目 | 値 | 理由 |
|---|---|---|
| モデル | **Pretrained** `models/Geneformer/Geneformer-V2-316M`（V2, `emb_layer: -1`） | (1) 本チェックの目的は「**ISP が実際に使う埋め込み空間**で start と goal が分かれるか」を測ることで、ISP は事前学習モデルで行う計画（HYP-20260925-02）だから。別のモデルで測ると前提確認にならない。(2) この細胞系（腎線維芽細胞・TLS オーガナイザー様）に対する fine-tuning 済みモデルは存在しない。EXP-20260926-003 の FT モデルは近位尿細管の障害軸用で、細胞も分類ラベルも違う。(3) FT を新たに作るとリモート GPU が必要で、事前チェックの規模を超える（G2 相当）。FT した埋め込みはラベルに合わせて分離が強くなるため、事前チェックの指標としては甘くなる |
| 埋め込み種別 | **`emb_mode: cls`** | V2 は入力先頭が `<cls>`。guardrails（2026-09-25）で `cls` 必須。状態のセントロイドと precheck の細胞埋め込みは同じ種別で取る（`run_isp.py` で連動） |
| dtype | **bf16**（Labretary 2026-09-27 変更。当初案は fp32） | PI 方針で bf16 が標準（2026-09-26）。後の摂動の本番も bf16 で回すので、前提も同じ精度の埋め込みで測る。bf16 は細胞状態間の距離を順位や符号で比べる用途で問題になった（shift ≈1e-4）。本チェックは状態の分離（AUC）を見るので影響は小さいと見込むが、未確認であることを限界に書く |
| デバイス | **mps**（Mac）、`nproc: 1` | guardrails（2026-09-26）: MPS でモデルを使うプロセスでは `datasets` の `num_proc > 1` が落ちる |
| `forward_batch_size` | 8 | EXP-20260926-001 の MPS 実測値 |

## 4. 事前チェックの設計

### 4.1 手順（`run_precheck.sh` に固定）

1. `scripts/inspect_h5ad.py` で A1 を確認（**実施済み**: `counts_integer_like: true`, `has_raw: true`, Ensembl 比率 1.00, Geneformer 語彙との重なり 19,225 遺伝子 → `work/inspect_A1.json`）
2. `scripts/geneformer/prepare_h5ad.py --config config.yaml` — `raw.X` の整数カウントを X に置き、線維芽細胞系のみ抽出、`TLS_group` を donor から写す → `work/input_fib.h5ad`（非整数なら停止する。guardrails）
3. `analysis/label_fib_states.py --config config.yaml` — **tokenize の前に生カウントであることを再確認**（X.data がすべて整数、行和 ≈ `n_counts`、`raw` を使ったことの記録）してから 3 遺伝子の検出で `fib_tls_state` を付ける → `work/input.h5ad` と `work/label_summary.json`
4. `scripts/geneformer/run_isp.py --config work/config-seed{S}.yaml --precheck`（S = 0, 1, 2）— tokenize（初回のみ、`work/tokenized.dataset`）→ 状態ごとの細胞抽出 → CLS 埋め込み → donor 分割 AUC。**`--precheck` 以外のモードでは起動しない**（ISP は G1・G2 の後）
5. `analysis/precheck_ci.py --config config.yaml --seed {S}` — 検出された遺伝子数の分布をそろえる操作・CI・並べ替え検定・感度分析（追加の forward は不要。4 の埋め込み CSV を再利用）
6. `analysis/precheck_ci.py --config config.yaml --summary` — 3 seed をまとめ、§5 の基準を機械的に当てはめて `results-precheck/summary.json` を書く
7. `scripts/capture_env.sh` で `results-precheck/env.json`

### 4.2 細胞の選び方（事前登録）

- 状態ごと **最大 800 細胞**、**donor ごと最大 60 細胞**（`state_embs`）。seed ごとに抽出をやり直す
- seed 0 での期待値（本日算出、`work/expected_capped_counts.json`）: goal 711 細胞 / 65 donor、start 800 / 22 donor、alt 800 / 25 donor = **2,311 細胞/seed**
- **検出された遺伝子数の分布をそろえる操作（主判定の集合）**: goal の各細胞に対し、log(トークン長)が最も近い start 細胞を 1 つずつ、置換なしで対応付ける（許容 |Δlog| ≤ 0.10）。トークン長は tokenize 済みデータの `length`（入力に載った遺伝子数）を使い、埋め込みは一切見ないで決める。seed 0 の期待値は両群 557 細胞、検出遺伝子数の中央値はどちらも 1,687
- start vs alt も同じ方法でマッチさせて算出する（ISP の文脈比較の前提確認。参考値）

### 4.3 donor の扱い（学習側と評価側を分ける）

- 交差検証は **`GroupKFold(n_splits=5)`、group = `donor_id`**。ある donor の細胞は必ず 1 つの fold にだけ入るので、**同じ donor が識別器の学習側と評価側の両方に入ることはない**
- goal と start の両方に細胞を出す 21 donor も、donor 単位で 1 つの fold に入る（細胞単位で分けない）
- 片方のクラスしかない fold は AUC が定義できないので除き、使った fold 数を記録する
- 層化分割（`StratifiedKFold`）の値も `run_isp.py` が出すが、**donor が学習・評価にまたがるため参考値**であり、判定には使わない（plan で固定）
- donor ごと 60 細胞の上限で、1 donor が識別を支配しないようにする

### 4.4 AUC の算出法と信頼区間

- 特徴量: CLS 埋め込み（1,024 次元、`emb_layer: -1`）。前処理は `StandardScaler`（fold の学習側で fit）
- 識別器: ロジスティック回帰（`C = 1.0`, `max_iter = 5000`）。start = 1, goal = 0 のラベルで学習し、`roc_auc_score` を評価側で計算
- **fold 平均 AUC**（5 fold の平均 ± SD）を主指標にする
- **95% 信頼区間**: 各 fold の評価側予測を集めた out-of-fold 予測に対し、**donor 単位のブートストラップ**（donor を復元抽出、2,000 回、seed 0）で AUC を再計算し、2.5% / 97.5% 分位を CI とする。細胞単位のブートストラップは donor 内の相関を無視して CI を狭く見せるので使わない
- **seed 間**: seed 0/1/2（細胞抽出が変わる）で同じ手順を繰り返し、3 seed の平均と範囲を報告する

### 4.5 並べ替え検定（偶然に 0.65 を超える確率の見積もり）

- 帰無分布: **donor の中で** `fib_tls_state` ラベルを入れ替える（donor ごとのクラス比率は保つ）操作を 200 回行い、そのつど同じ交差検証で fold 平均 AUC を出す
- p = (#{AUC_perm ≥ AUC_obs} + 1) / 201
- この帰無仮説は donor 構成の効果を残したままなので**保守的**（donor 由来の分離は「偶然」に数えられる）。狙いは「donor の違いでは説明できない分離か」を見ること

## 5. 判定基準（事前登録。結果を見る前に固定）

主判定は §4.2 の **検出された遺伝子数の分布をそろえた後の集合**で行う。

| 記号 | 基準 | 閾値の根拠 |
|---|---|---|
| **C1（主）** | goal vs start の fold 平均 AUC（検出された遺伝子数の分布をそろえた後）が **3 seed の平均で ≥ 0.65**、かつ **3 seed のうち 2 つ以上で個別にも ≥ 0.65** | HYP-20260925-02 が事前登録した閾値。偶然（0.5）から明確に離れており、ISP のセントロイドシフトを解釈できる最低限の分離 |
| **C1b** | C1 の donor ブートストラップ **95% CI の下限 ≥ 0.55**（3 seed の中央値） | 点推定だけで通さない。donor 単位の CI 下限が 0.55 を超えることを求める |
| **C1c** | 並べ替え検定 **p < 0.05**（3 seed のうち 2 つ以上） | 200 回の並べ替えで得られる最小 p は 1/201 = 0.005。donor 構成を保った保守的な帰無分布 |
| **C2（交絡）** | 検出された遺伝子数と総 UMI 数のみ（log 総カウント・log トークン長）の AUC が、同じマッチ集合で **≤ 0.60**、かつ **C1 の AUC より 0.05 以上低い** | マッチ前の検出された遺伝子数と総 UMI 数のみの AUC は 0.676（実測）。マッチで 0.526 に落ちることを確認済み。マッチが効いていることの確認 |
| **C3（陽性対照 / パイプライン）** | 同じ埋め込みで **FIB vs M-FIB**（著者ラベル、皮質 vs 髄質の線維芽細胞）の donor 分割 AUC が **≥ 0.90** | 埋め込み抽出が正常なら確実に分かれるはずの対照。追加の計算は不要。ここが低ければ埋め込み側の不具合を疑う |

**合格（ISP へ進めてよい）** = C1 かつ C1b かつ C1c かつ C2 かつ C3 をすべて満たす。
**不合格** = C1 を満たさない。
**判定不能** = C1 は満たすが C2 または C3 を満たさない（分離が検出された遺伝子数や総 UMI 数、パイプラインの不具合で説明できてしまう）。

- 合否は本 EXP の範囲。**仮説の Status（proposed → accepted など）は変えない**。G1 の判断は PI が行う
- 参考として、マッチ前の AUC、start vs alt の AUC、感度分析（皮質のみ / goal を AKI・CKD のみ / `subclass_l2 == FIB` のみ / start と同じ donor の goal 細胞のみ）も同時に報告する。**これらは判定には使わない**（事後に基準を選び直さないため）

## 6. コントロール

| 種別 | 内容 | 目的 |
|---|---|---|
| 陽性対照 | FIB vs M-FIB の識別（C3） | 埋め込み抽出・交差検証の配管が正しいことの確認。働かなければ C1 の解釈は不能 |
| 陰性対照 | 検出された遺伝子数と総 UMI 数のみの識別（C2） | ラベルが検出された遺伝子数と総 UMI 数の代理でないことの確認 |
| 比較対照 | start vs alt（`FIB_neg_TLSlow`）の AUC | 後続 ISP の文脈比較（TLS-high vs TLS-low）の前提。**検出された遺伝子数と総 UMI 数のみでも 0.676 だった**（本日実測）ので、こちらも検出された遺伝子数の分布をそろえた後の値を主に見る |
| 帰無 | donor 内ラベル並べ替え 200 回（C1c） | 偶然合格する確率の見積もり |

本 EXP に摂動は含まれないため、ISP 用の陽性・陰性対照遺伝子（CCL19/CXCL13、ハウスキーピング、発現頻度マッチのランダム 20 遺伝子）は `config.yaml` に書くだけで**実行しない**。

## 7. 頑健性チェック

- seed 0/1/2（細胞抽出が変わる）で 3 回。3 seed の平均・範囲を報告
- 検出された遺伝子数の分布をそろえた場合／そろえない場合の両方を報告（主判定はそろえた場合）
- 感度分析 4 種（§5 の注）
- `emb_layer` は -1 のみ（層の感度分析は ISP 本番の範囲。ここでは増やさない）

## 8. 計算量見積もり

| 段階 | 細胞・forward 回数 | 見積もり時間（Mac, MPS） | 根拠 |
|---|---|---|---|
| A1 の読み取り・prepare | 304,989 細胞を 1 回走査（うち 21,769 を保持） | 5 分以内 | 本日の同等の全走査で読み取り 12.8 秒（`work/expected_capped_counts.json`）。h5ad 書き出しを含めても数分 |
| ラベル付与 | 21,769 細胞 | 1 分未満 | スパース行列の 3 列参照 |
| tokenize（1 回だけ、`nproc: 1`） | 21,769 細胞 | 1〜3 分 | EXP-20260926-001 は 17,000 細胞規模で約 10 秒（`nproc: 8`）。1 プロセスに落として余裕を見る |
| CLS 埋め込み（precheck 本体） | 2,311 細胞 × 3 seed = **6,933 forward** | 6〜8 分/seed × 3 = **18〜24 分** | EXP-20260926-001 phase 0 の precheck は 1,405 細胞で 275 秒（約 0.20 s/細胞、中央値 2,524 遺伝子）。線維芽細胞は中央値 1,609 遺伝子なのでやや速い。fp32 はモデルロードを含めて同程度 |
| 解析（CI・並べ替え 200 回・感度分析） | 追加 forward なし | 10〜30 分（CPU） | 1,114 細胞 × 1,024 次元のロジスティック回帰を 200 × 5 回 |
| **合計** | **forward 6,933 回 / 使用細胞 2,311（延べ 6,933）** | **約 0.6〜1.1 時間** | |

`run_isp.py --config config.yaml --dry-run`（2026-09-27 実行、計算なし）は `{"forward_passes": 9300, "gpu_hours_estimate": 0.02, "needs_G2_approval": false}` を返した。ただしこの値は **config の ISP セクション（摂動 6 + 対照 23 遺伝子）を A100 相当の速度で見積もったもの**で、本 EXP の事前チェックの費用でも、MPS の実測でもない。判断に使うのは上の表である。ISP の本番費用は、実測の `compute.sec_per_cell_gene` を入れて G2 申請時に出し直す。

**G2 の要否: 不要**（`Approved-by:` は空のまま、PI が必要と判断すれば記入）。
根拠は CLAUDE.md §3 の目安（GPU 2 時間超 / 1 万細胞超 / 全遺伝子スクリーニング）のいずれにも当たらないこと: 実行時間は最大 1.1 時間、埋め込みを取る細胞は 2,311（延べ 6,933、いずれも 1 万未満）、摂動・スクリーニングは行わない。`run_isp.py --precheck` の自動見積もり（`estimate()`）でも `needs_G2_approval: false` になる想定だが、実行ログの値を `results-precheck/` に残す。

**ただし本番 ISP は別**: HYP-20260925-02 の見積もりは 1.2〜2.0 時間（start = TLS-high の 1 本）、alt を start にした 2 本目を含めると 2.4〜4.0 時間 → **G2 必須**。ISP の計画は本 EXP の合格と G1 採択の後に別途 `plan.md` を作って事前登録する。

## 9. 想定される失敗と対策

| 失敗 | 兆候 | 対策（事前に決めておく） |
|---|---|---|
| **C1 不合格（AUC < 0.65）** | 主判定が閾値未満 | **ISP は行わない**。HYP-20260925-02 の事前チェック条項どおり、goal の定義をやり直す案を PI に提出する。案 (a) B4 Meylan 2022 Visium の TLS スポット由来シグネチャーでスコア化して上位を goal にする、案 (b) 閾値を上げる（3 遺伝子のいずれか ≥ 2 カウント、または 2 遺伝子以上が陽性）。**いずれも新しい EXP として事前登録し直す**（本 EXP の基準を後から書き換えない）。仮説の Status は変えず、判断は PI（G1）に返す |
| **C2 不合格（検出された遺伝子数と総 UMI 数で説明できる）** | 検出された遺伝子数と総 UMI 数のみの AUC が C1 と同水準 | 判定不能として報告。マッチの許容（0.10）や特徴量（総カウント／トークン長）を変えるのではなく、goal 定義の見直し（上記案 a/b）と、検出された遺伝子数の分布をそろえた再抽出を新 EXP で提案する |
| **C3 不合格（FIB vs M-FIB < 0.90）** | 対照が働かない | 判定不能。埋め込み抽出（`emb_mode`、層、モデルのロード、dtype）を点検してから再実行。C1 の値は解釈しない |
| tokenize 前のカウントが非整数 | `prepare_h5ad.py` が停止、または `label_fib_states.py` の再確認で不一致 | そこで止める。`raw.X` を使えているか（`has_raw: true`）を確認。A1 は本日 `counts_integer_like: true` を確認済み |
| MPS で `datasets` の子プロセスが落ちる | `One of the subprocesses has abruptly died` | `nproc: 1` を config で固定済み（guardrails 2026-09-26）。ログで `nproc` を確認 |
| goal の fold にクラスが片方しかない | `n_folds_used < 5` | 記録して報告。3 fold 未満しか使えない場合は判定不能とし、donor ごとの上限（60）と状態ごとの上限（800）を見直す案を出す（本 EXP では変えない） |
| 他の実験の tokenized を誤って参照 | — | `config.yaml` の `tokenized_dataset` は本 EXP の `work/` のみ。他 EXP の `*.dataset` は指さない（guardrails 2026-09-26） |
| 途中でエージェントが止まる | — | 各段階の出力（`work/*.json`, `results-precheck/seed*/`）をファイルに残し、次のエージェントが続きから実行できるようにする |

## 10. 承認

- 本 EXP（事前チェック）: G2 **不要**（§8）。実行は PI の指示があれば直ちに可能
- 本 EXP の合否は ISP 実行の前提であり、**仮説の採択（G1）ではない**。G1 は PI が判断する
- 本番 ISP: **G2 必須**（§8）。本 EXP 合格・G1 採択の後に別 `plan.md` を事前登録する

Approved-by:

## 逸脱の記録（実装の不具合修正、基準は不変）— 2026-09-27

**事象**: 実行時、`analysis/precheck_ci.py::load_cells` が「`cell_embs.csv` と `state-embs-seed{S}-n800-d60.dataset` の (状態, donor) の並びが一致しない」で停止した（seed 0/1/2 とも `run_isp.py --precheck` 自体は完走済み）。

**原因**: Geneformer の `EmbExtractor.extract_embs` は内部の `downsample_and_sort`（`../Geneformer/geneformer_hf/geneformer/perturber_utils.py` の `downsample_and_sort`）で、`max_ncells=None` のときシャッフルはせず、`length` 降順の**安定ソート**（同じ length のもの同士は元の順序を保つ）で行を並べ替えてから埋め込みを計算する。`cell_embs.csv` はこの並べ替え後の順序で `fib_tls_state`/`donor_id` を保持するが、`state-embs-seed{S}-n800-d60.dataset` は並べ替え前の順序のままだったため、位置で突き合わせると食い違う。

**確認した事実（結合前に検証）**:
1. `state-embs-seed{S}-n800-d60.dataset` を「length 降順・安定ソート（`np.argsort(-length, kind="stable")`）」で並べ替えると、3 seed とも全 2,311 行で `cell_embs.csv` の `(fib_tls_state, donor_id)` 列に**完全一致**した（位置ごとの一致 2,311/2,311、3 seed とも）。つまり並べ替え規則自体はこの安定ソートで説明できる。
2. ただし、同じ `(fib_tls_state, donor_id, length)` の組み合わせを持つ細胞が複数存在する「タイ」があり（seed0: 32 組・101 細胞 / seed1: 26 組・89 細胞 / seed2: 37 組・103 細胞、全体の 3.9〜4.5%）、そのうち **`n_counts` はほぼ全てのタイ組で異なり**（32/32, 26/26, 37/37 組）、`subclass_l2` も約半数のタイ組で異なった（15/32, 14/26, 10/37 組）。length・donor・状態だけでは、タイの中でどの細胞がどの埋め込み行に対応するか一意に決まらない。
3. 事前登録した §4.4 の主判定（C1〜C1c）は埋め込み `X` と `(state, donor)` だけを使うため、この位置合わせが崩れても **A の判定自体は誤らない**（`cell_embs.csv` 内で埋め込みと `(state, donor)` は元々正しく対応している）。しかし C2（検出された遺伝子数と総 UMI 数のみの AUC、`n_counts` を使う）と C3・感度分析（`subclass_l2`/`region`/`disease` を使う）は、タイのある細胞（3.9〜4.5%）で誤った補助情報を結合する恐れがあった。

**対応（実施）**: 位置ではなく **一意な行 ID（`row_id`）で結合する**方式に修正した（判定基準・細胞選択・検出された遺伝子数の分布をそろえる方法は変えていない）。
1. `scripts/geneformer/run_isp.py` に `state_emb_dataset_with_rowid()` を追加。`state-embs-seed{S}-n800-d60.dataset` に `row_id`（元の行番号、0-indexed）列を足した `state-embs-seed{S}-n800-d60-rowid.dataset` を作る（追加の forward なし。列を足すだけ）
2. `precheck()` の `EmbExtractor` 呼び出しで、この `-rowid.dataset` を入力にし、`emb_label` に `"row_id"` を追加（`[key, dkey, "row_id"]`）。これにより `cell_embs.csv` にも `row_id` 列が残る
3. `analysis/precheck_ci.py::load_cells` を、位置合わせではなく `row_id` での ID 結合に変更。結合前後で行数・`(state, donor)` ごとの件数が一致することを `assert` する
4. 埋め込みの取り直しが必要だったため、3 seed 分の `run_isp.py --config work/config-seed{S}.yaml --precheck` を再実行した（forward 再計算、約 26 分。判定基準・細胞選択・状態ラベルの定義・検出された遺伝子数の分布をそろえる際の許容は変更していない）

**変更していないもの**: C1〜C3 の閾値、検出された遺伝子数の分布をそろえる方法（log トークン長・許容 0.10）、細胞選択（状態ごと最大 800・donor ごと最大 60）、goal/start/alt の定義、seed（0/1/2）、判定の合否ロジック（`summarize()`）。仮説 HYP-20260925-02 の Status も変更していない。

# EXP-20260926-004 実験計画 — FT モデルの ISP の妥当性確認: (a) DE マッチ比較遺伝子、(c) donor 単位の解析と held-out donor の拡張

- 仮説: [HYP-20260926-03](../../hypotheses/HYP-20260926-03.md)（方法仮説, Status: proposed。PI は 2026-09-26 に (a)(b)(c)「全て実行」を指示。G1 の記録は Labretary）
- 前提の実験: [EXP-20260926-003](../EXP-20260926-003/plan.md)（FT Geneformer V2-316M 2 クラス分類器、対照のみ ISP。C1 は HAVCR1 で不合格、VCAM1 は null 上位、split-half 0.46）
- 作成: Labretary（experiment-planner）/ 2026-09-26
- 担当範囲: (a) と (c)。(b) 陽性対照の見直しは別担当（本実験では扱わない）
- **事前登録の時点**: §A（DE マッチ）は DE マッチ遺伝子の ISP を実行する前、§C（donor 単位）は donor 単位の指標を対照で計算する前（EXP-003 の細胞単位の結果は既知）に記入した。記入後の変更は §11 に日時と理由を残す
- 盲検: 主要対象 TP53 / CDKN1A / LTBR と参照 SPP1 / LCN2 は回さない・解析しない（freq_matched の頻度参照としてのみ config に残る。EXP-003 と同じ）

## 1. 問い
- (a) FT モデルで VCAM1 の delete が null 上位 10% に入ったのは、「分類器が使う aPT マーカーを消した」だけで説明できるか。VCAM1 / HAVCR1 と同程度に aPT で上昇し（log2FC・発現細胞数・入力内の順位でマッチ）、傷害・修復の駆動因子として知られていない遺伝子（DE マッチ比較遺伝子）を消したときと比べて、VCAM1 は goal（健常 PT）方向へより大きく動くか
- (c-1) EXP-003 の細胞単位の指標を donor 単位に集約すると、C1〜C3・C3-panel の判定と信頼性（split-half 0.46）はどう変わるか（GPU 不要）
- (c-2) FT に使っていない CKD の TLS-mid donor（15 人）を start に加えて donor 数を 14 → 29 に増やすと、判定と donor 単位の信頼性は再現・改善するか

方法の検証であり、生物学的な仮説の検証ではない。in silico の結果は「分類器・埋め込みの挙動」についての所見で、遺伝子の機能の実験的事実ではない。

## 2. データ
| 項目 | (a)・(c-1) | (c-2) 拡張 |
|---|---|---|
| データセット | A2 KPMP v1.5 scRNA（Lake 2023, PMID 37468583）。EXP-003 と同じ tokenized（`work/tokenized.dataset` = データ部分の APFS clone, sha256 `ba3a14c0…` 一致, `work/tokenized.SHA256SUMS`） | 同じ raw（`data/raw/kidney_tls/A2_KPMP_v1.5_scRNA_Lake2023.h5ad`, 読み取りのみ）から `prepare_h5ad.py` で新しい入力 `expanded/work/input.h5ad` を作り tokenize |
| start | aPT_CKD_TLShigh 7,533 細胞 / 14 donor（全員 FT で held-out） | aPT ∧ CKD ∧ TLS ∈ {high, mid} = 14,631 細胞 / 29 donor（TLS-mid 15 donor 7,098 細胞は EXP-001 の入力から除かれ FT に未使用。§C.5 で漏洩なしを assert） |
| goal | PT_healthy。セントロイドは held-out 健常 9 donor（EXP-003 と同じ） | 同じ 9 donor |
| alt | aPT_CKD_TLSlow（FT の train/eval → 参考値） | 同じ |
| 状態ラベル | `isp_state`（推定。TLS 群は donor pseudo-bulk の TLS スコア三分位で組織学的な TLS ではない） | 同じ定義（TLS-mid も同じ表 `donor_A2.csv` の値をそのまま使う） |

## 3. モデル
- FT 済み Geneformer V2-316M CellClassifier（`experiments/EXP-20260926-003/models/ft-main/final`, 2 クラス 0 = aPT_CKD, 1 = PT_healthy）。emb_mode cls, emb_layer −1, bf16 forward + fp32 コサイン（PI 標準）
- 実行: spark-cf5b（GB10, CUDA）。EXP-003 phase 3 と同じ機械（デバイス差を比較に持ち込まない）。`capture.embeddings`・`capture.classifier_logits` あり

## A. (a) DE マッチ比較遺伝子（事前登録）
### A.1 選定の考え方
- 陽性対照 VCAM1・HAVCR1 は aPT のマーカーでもあり、start と goal を分けるよう学習した分類器が直接使っている可能性がある（EXP-003 plan §4.5 / HYP-20260926-02 §循環性）。「aPT で強く上昇している遺伝子なら何を消しても goal 方向へ動く」なら、VCAM1 の C1 合格は駆動因子の証拠にならない
- そこで、VCAM1 / HAVCR1 と**発現の差・頻度・入力内の順位が同程度**で、**既知の駆動因子でない**遺伝子を比較群にする。選定は発現データのみ（ISP の結果は読まない）

### A.2 選定規則（`analysis/de_matched_select.py`。統計は `outputs/literature/2026-09-26-negative-controls/genome_stats.csv`: EXP-001 input.h5ad の生カウントから donor pseudo-bulk で計算済み）
| 規則 | 内容 |
|---|---|
| E1 | Geneformer V2 の辞書に収載 |
| E2 | start（TLS-high aPT）vs goal（健常 PT）の donor pseudo-bulk log2FC ≥ 2 かつ Welch BH q < 0.05（上昇） |
| E3 | start で Geneformer 入力に入る細胞 ≥ 440（HAVCR1 の 659 の 1/1.5）かつ 25 細胞以上で検出される start donor ≥ 8（200 細胞・donor ≤ 25 の抽出が可能） |
| E4（自動除外） | negctrl_candidates.py の既知リスト（傷害・aPT マーカー、TLS・免疫、解離ストレス、性染色体、細胞周期、ストレス応答）＋本実験の追加カテゴリ（傷害・老化のメディエーター、補体・IFN・炎症、修復シグナル経路 HGF/MET・EGFR・Wnt・Notch・Hedgehog・TGF-β・YAP、ECM・線維化、酸化ストレス応答、低酸素）、盲検対象、陽性対照 |
| E5（文献レビュー） | E4 を通った候補を距離順に文献知識で確認し、PT の傷害・修復・老化・TLS の駆動因子とする報告があれば除外。判断と不確実性を記録（下表） |
| 距離 | d_r(g) = max(\|Δlog2FC\| / 1.0, \|ln(細胞数比)\| / ln 1.5, \|Δ相対順位\| / 0.10)、r = VCAM1, HAVCR1 |
| 選び方 | r ごとに d の小さい順に E4・E5 を通る遺伝子を 5 つ → 計 10 遺伝子（VCAM1 帯 5 + HAVCR1 帯 5） |

- プール: E1〜E3 を満たす 1,690 遺伝子中の上昇遺伝子（`de_matched/selection.json`）

参照（陽性対照）:
| 遺伝子 | log2FC | start の入力内細胞数（検出率） | 健常の検出率 | 入力内の相対順位 | donor ≥ 25 細胞 |
|---|---|---|---|---|---|
| VCAM1 | 5.26 | 4,674（0.62） | 0.062 | 0.36 | 13 |
| HAVCR1 | 4.25 | 659（0.088） | 0.0005 | 0.25 | 10 |

選定結果（`de_matched/selected.csv`）:
| 帯 | 遺伝子 | d | log2FC | 細胞数（検出率） | 健常の検出率 | 相対順位 | donor ≥ 25 | 除外しなかった根拠（不確実性） |
|---|---|---|---|---|---|---|---|---|
| VCAM1 | PPFIBP1 | 0.98 | 4.28 | 3,761（0.50） | 0.052 | 0.42 | 13 | liprin-β1（接着斑の足場）。PT 傷害の駆動因子の報告を把握していない（中） |
| VCAM1 | AKAP12 | 1.22 | 6.48 | 3,058（0.41） | 0.008 | 0.41 | 13 | gravin。傷害 PT で上昇する記載はありうるが因果の証拠は把握していない（**高**） |
| VCAM1 | SELENOM | 1.24 | 4.02 | 3,561（0.48） | 0.068 | 0.48 | 13 | セレノプロテイン M。腎傷害での役割の報告なし（中） |
| VCAM1 | DMKN | 1.27 | 6.12 | 2,790（0.37） | 0.004 | 0.36 | 13 | dermokine。腎での機能報告なし（中） |
| VCAM1 | SCIN | 1.51 | 6.31 | 2,539（0.34） | 0.004 | 0.50 | 13 | scinderin。腎傷害での役割の報告なし（中） |
| HAVCR1 | ANKRD2 | 0.21 | 4.46 | 640（0.085） | 0.000 | 0.24 | 8 | 筋の ankyrin repeat タンパク質。腎での機能報告なし。アンビエント・ダブレット由来の可能性（中） |
| HAVCR1 | ARRDC4 | 0.37 | 3.88 | 664（0.089） | 0.004 | 0.28 | 11 | α-arrestin（糖代謝応答）。腎傷害の駆動因子の報告なし（中） |
| HAVCR1 | IDUA | 0.57 | 3.68 | 591（0.079） | 0.002 | 0.23 | 11 | リソソーム酵素（低） |
| HAVCR1 | SLITRK4 | 0.63 | 4.45 | 801（0.107） | 0.001 | 0.19 | 12 | 神経の接着分子（低） |
| HAVCR1 | ADAMTS16 | 0.65 | 4.62 | 857（0.114） | 0.001 | 0.25 | 13 | 血圧関連座位。PT 傷害の駆動因子の報告なし（中） |

除外した候補（距離順に現れたもの。`de_matched/candidates.csv`）:
- E4 自動: ITGB8, CD24, CLU, TPM1（既知の aPT/傷害マーカー）、COL4A2, COL4A1（ECM・線維化）、TNFRSF12A（Fn14, 腎傷害・炎症のメディエーター）、PON2・NFE2L3（酸化ストレス応答）、WNT2B（Wnt 経路）
- E5 文献: MACC1（MET の転写調節因子 → HGF/MET 修復経路, 経路単位で保守的に除外, 不確実性 高）、C1QTNF3（CTRP3。腎線維化・糖尿病性腎症モデルで保護的との報告がある — 記憶に基づき未検証, 高）、SCARA3（別名 CSR1 = cellular stress response 1, 酸化ストレスで誘導, 中）
- 注意: 「駆動因子の報告を把握していない」は文献データベースの系統的検索ではなく知識に基づく判断（literature-scout による確認は未実施）。AKAP12 は特に不確実。帯内の距離は VCAM1 帯で 0.98〜1.51（VCAM1 ほど高頻度かつ強く上昇し、既知マーカーでない遺伝子は少ない）、HAVCR1 帯で 0.21〜0.65

### A.3 問いと判定（事前登録）
主指標は EXP-003 と同じ **a2**（長さ 5 分位ビンごとの null 平均 shift_goal を細胞ごとに引き、遺伝子×seed で平均。null は EXP-003 results-ft の同 seed の 20 遺伝子）。DE マッチ遺伝子は EXP-003 の null と統合して同じ補正をかける（null 細胞だけで補正を推定するので、DE 遺伝子を加えても null・陽性・HK の値は変わらない）。

| # | 量 | 定義 |
|---|---|---|
| A-P（主） | VCAM1 の DE 内分位 | seed ごとに、DE マッチ 10 遺伝子のうち a2 が VCAM1 より小さい割合 |
| A-S1 | 帯内比較 | VCAM1 > VCAM1 帯の 5 遺伝子すべて、HAVCR1 > HAVCR1 帯の 5 遺伝子すべて（seed ごと） |
| A-S2 | DE マッチの C1 偽陽性率（false positive rate） | DE 遺伝子が C1 と同じ基準（null 分位 ≥ 0.90 が 3 seed 中 2 以上）を満たす数 / 10 |
| A-S3 | DE マッチ全体のずれ | DE 遺伝子の null 分位の平均（seed ごと, 帰無では ≈ 0.5）。遺伝子単位（3 seed 平均）で 0.5 に対する Wilcoxon 符号順位（片側, 参考） |
| A-R | 副指標 | A-P〜A-S3 を raw shift_goal と dP_goal_a2（ΔP(PT_healthy) = softmax(logits_pert)[1] − softmax(logits_orig)[1] に a2 と同じ長さビン補正）でも計算（判定には使わない。a2 と食い違えば報告） |

判定（a2, 事前登録）:
- **VCAM1 の DE 内分位 ≥ 0.90（DE 10 遺伝子中 9 以上を上回る）が 3 seed 中 2 以上** → 「VCAM1 の goal 方向の移動は、同程度に上昇する非ドライバー遺伝子の除去より大きい（マーカー除去だけでは説明できない）」。ただし駆動因子の証拠ではなく「分類器・埋め込みで VCAM1 が特に効く」ことの示唆にとどめる
- 満たさない → 「**VCAM1 の C1 合格はマーカー除去と区別できない**」と解釈する（EXP-003 の C1 の VCAM1 合格を駆動因子の証拠として使わない）
- 付記: A-S2 が 3/10 以上なら「aPT で上昇する非ドライバー遺伝子の除去でも C1 を満たしうる → 現在の C1 は駆動因子に特異的でない」と報告する。A-S3 の平均分位が 0.75 以上なら「上昇遺伝子の除去一般が goal 方向へ動かす（マーカー除去効果）」の示唆とする
- HAVCR1（EXP-003 で C1 不合格）も同じ量を報告のみ

### A.4 実行と統合の妥当性
- `config-de.yaml`: EXP-003 `config.yaml` と同じ（モデル・細胞抽出・セントロイド・seed 3/4/5・capture）で、摂動する遺伝子だけを変える: DE マッチ 10（role `de_matched`）＋一致確認用に VCAM1・GAPDH（全 seed）と EXP-003 seed 3 の null 1 遺伝子 ENSG00000112297（role `recheck_null`）。`negative_random_n: 0`（null は回さない）
- 細胞は遺伝子と seed だけで決まる（`random.Random(f"{seed}-{gene}")`）ので、他の遺伝子を加えても変わらない。状態セントロイドの細胞も seed と config の抽出条件だけで決まる
- **統合の条件（事前登録）**: 再実行した VCAM1・GAPDH（3 seed）と null 1 遺伝子（seed 3）で、(i) 細胞（cell_uid の並び）が完全一致、(ii) 状態セントロイドベクトル（state_emb_vectors.npz）の最大差、(iii) 細胞ごと shift_goal の最大差 ≤ 1×10⁻⁶。満たさなければ統合せず、差の大きさを報告し、null・陽性・HK をこの run の中で回し直す案（約 0.6 h）を PI に出す
- 見積もり: 13 遺伝子 × 200 細胞 × 3 seed = 7,800 細胞×遺伝子。EXP-003 の spark-cf5b 実測（20,400 / 44 min ≈ 0.13 s）から **約 17〜20 分**（G2 不要。PI 実行承認済み）

## C. (c) donor 単位の解析（事前登録）
### C.1 donor 単位の指標の定義（対照で計算する前に固定）
入力は EXP-003 と同じ細胞単位の a2（`reanalysis_cosine.corrected` の m_a2, seed ごと）。raw・dP_goal_a2 も同じ手順で副次に計算する。
| 名前 | 定義 |
|---|---|
| donor 平均 x(s,g,d) | seed s・遺伝子 g・donor d の細胞の a2 の平均。**細胞 3 未満の donor は除く** |
| **DL-med（主）** | S(s,g) = median_d x(s,g,d)。null 分布は同 seed の null 遺伝子の S |
| DL-mean（副） | S(s,g) = mean_d x(s,g,d)（donor を等重みにした平均） |
| DL-sign（報告） | 遺伝子ごとに、x を seed で平均した donor 値 x̄(g,d) > 0 の donor 割合。null 遺伝子の同じ割合の分布での分位 |
| DL-excess（報告） | 陽性対照ごとに y(g,d) = seed 平均の [x(s,g,d) − 同 seed・同 donor の null 遺伝子の x の平均]。y > 0 の donor 割合と Wilcoxon 符号順位（片側, donor が単位） |
- 判定基準は EXP-003 と同じ閾値の C1・C2・C3・C3-panel を DL-med（主）と DL-mean で計算する
- 信頼性:
  - **donor split-half**: seed ごとに donor を 2 群にランダム分割（20 回, rng seed 0。EXP-002/003 の C4i と同じ手順）し、各半分で DL-med を計算して 34 遺伝子の Spearman ρ を平均。同じ分割で細胞単位の平均（EXP-003 の C4i 相当）も計算して並べる。Spearman–Brown 補正（2ρ/(1+ρ)）も報告
  - **leave-one-donor-out（LODO）**: donor を 1 人ずつ除いて DL-med を再計算し、全 donor との Spearman ρ（seed ごと, 最小・中央値）と、除いたときに C1・C2・C3・C3-panel と VCAM1/HAVCR1 の C1 状態が変わる donor の数
- 解釈（事前登録）: DL-med で C1〜C3・C3-panel の判定が細胞単位と同じなら「donor の重みづけは結論を変えない」。HAVCR1 が DL-med で C1 を満たしても HYP-20260926-02 の判定（細胞単位の事前登録）は変えず、感度分析として報告する。donor split-half ≥ 0.7 なら「donor 単位の指標は順位の再現性が十分」、< 0.7 なら「donor 間の違いによる不安定さが残る」
- 比較: 同じ定義を EXP-003 results-pretrained-matched（Pretrained, 同じ細胞）にも適用
- 実装: `analysis/donor_level.py`（CPU のみ）。EXP-003 の gene_seed_metrics.csv の m_a2 と、再実装の細胞単位の a2 が一致すること（最大差 ≤ 1×10⁻¹²）を先に確認する

### C.2 (c-2) held-out donor の拡張（設計）
- start = CKD の aPT のうち TLS-high（14 donor, 7,533 細胞）＋ TLS-mid（15 donor, 7,098 細胞）= 29 donor・14,631 細胞（ラベル `aPT_CKD_TLShighmid`）。goal・alt は EXP-003 と同じ定義。TLS-mid は EXP-001 の入力で除外され、FT の train/eval（split_manifest.json の 21 donor）に含まれない（prepare 後に assert）
- 入力: `expanded/config-expanded.yaml` の `prepare`（EXP-001 と同じフィルタ・同じ donor 表。state_rules だけ変更）→ `expanded/work/input.h5ad` → tokenize（Mac CPU）
- 細胞抽出: **遺伝子あたり最大 350 細胞・donor あたり最大 12 細胞**（29 × 12 = 348。donor 単位の平均に各 donor の細胞を確保するため donor 上限を 25 → 12 に下げ、総数を 200 → 350 に上げる）、最小 100
- 遺伝子: EXP-003 と同じ対照 14（陽性 2・HK 3・パネル 7・RP 2）＋新しい start での発現頻度マッチ null 20（seed ごと。start が変わるので EXP-003 と別の null 遺伝子になる）
  - **DE マッチ 10 遺伝子を加えるかの規則（結果を見る前に固定）**: (a) の実測時間＋拡張の見積もり（44 遺伝子）≤ 2.1 h なら加える（(a) の再現確認になる）。超えるなら 34 遺伝子のみ
- 状態セントロイド: EXP-003 と同じ抽出条件（状態ごと 500・donor ≤ 60、goal は held-out 健常 9 donor）。start の細胞が変わるので、抽出の乱数列の都合で goal セントロイドの細胞は EXP-003 と一部異なる（同じ 9 donor）
- seed 3/4/5、delete、bf16、capture あり
- 評価（事前登録）: 細胞単位（EXP-003 と同じ C1〜C5・C3-panel, evaluate_controls.py の関数）と donor 単位（§C.1 と同じ定義, 29 donor）。追加で (i) TLS-high 14 donor の donor 値の、EXP-003 と拡張 run の間の一致（固定対照 14 遺伝子の DL-med の Spearman, 参考）、(ii) donor split-half を 14 donor（TLS-high のみ）と 29 donor で比較。TLS-high と TLS-mid の比較は記述のみ（仮説検定をしない）
- 見積もり: 34 遺伝子 × 約 340 細胞（HAVCR1・RP は少ない）× 3 seed ≈ 35,000 細胞×遺伝子 × 0.13 s ≈ **1.25 h**（44 遺伝子なら ≈ 1.6 h）＋セントロイド 3 × 1,500 細胞 ≈ 0.05 h
- (a) と合わせて **約 1.6〜2.0 h**（PI 指示の上限 ~2 h 以内）→ 実行する。(a) の実測で 2.1 h を超える見込みになったら DE 遺伝子を外し、それでも超えるなら準備だけして G2 に回す

## 5. コントロール
- ポジティブ: VCAM1、HAVCR1（EXP-003 と同じ。見直しは (b) の担当）
- ネガティブ: 発現頻度マッチ null 20/seed、HK（ACTB, GAPDH, PPIA）、陰性対照パネル 7（C3-panel）、探索枠 RP 2
- (a) の比較群: DE マッチ 10（陰性対照ではなく「マーカー除去」の対照）

## 7. 頑健性チェック
- seed 3/4/5（各 seed で細胞と null が変わる）
- (a) の統合前の一致確認（A.4）。a2 と raw・dP の一致
- (c) DL-med と DL-mean、LODO、Pretrained との比較、拡張 run での再現

## 8. 計算量
| 項目 | 場所 | 見積もり |
|---|---|---|
| (a) DE マッチ ISP | spark-cf5b | 0.3 h |
| (c-1) donor 単位の再解析 | Mac CPU | 数分 |
| (c-2) 入力作成・tokenize | Mac CPU | ~10 分 |
| (c-2) 拡張 ISP | spark-cf5b | 1.25〜1.6 h |
| 合計 GPU | | **1.6〜2.0 h**（G2 の目安 2 h 以下。PI は実行を承認済み） |
- 1 万細胞超の条件: 拡張の start は 14,631 細胞あるが、ISP に使うのは遺伝子あたり ≤ 350 細胞（ユニーク細胞は最大 ~1 万程度）。全遺伝子スクリーニングではない

## 9. 想定される失敗と対策
| 失敗 | 対策 |
|---|---|
| 再実行で EXP-003 と一致しない（GPU の非決定性など） | 差の大きさを報告。1×10⁻⁶ を超えたら統合しない（A.4） |
| DE マッチ遺伝子の抽出細胞が 100 未満 | run_isp が skip → 報告し、帯の次点は追加しない（事後の差し替えをしない） |
| TLS-mid donor の混入が FT に漏れている | prepare 後に split_manifest の train/eval donor と照合し、重なれば停止 |
| リモートの接続断 | 遺伝子単位で再開（run_isp の既存機能） |
| リモートにヒト由来データが残る | 回収・件数確認後に `isp_remote.sh clean`（他の担当が同じホストを使っていないことを確認してから） |

## 10. 承認
- PI: 2026-09-26「(a) DE マッチ比較遺伝子 (b) 陽性対照の見直し (c) donor 単位解析 を全て実行」「FT モデル前提、リモート PC を積極的に」。GPU 合計 ≤ ~2 h の範囲で実行（超える場合は G2）
Approved-by: PetadimensionLab (PI), 2026-09-26 — 実行指示（全て実行, GPU ≤ ~2 h）。G1（HYP-20260926-03）は Labretary が記録

## 11. 変更履歴
- 2026-09-26 初版（ISP 実行前）
- 2026-09-26 20:07（(a) の ISP 完了後、(c-2) の ISP 実行前）: §C.2 の事前登録した規則に従い、(c-2) に DE マッチ 10 遺伝子を加える。根拠 = (a) の実測 0.31 h（ISP 935 s + セントロイド 176 s, 44 遺伝子×seed の実測 24.2 s/200 細胞 = 0.121 s/細胞×遺伝子）＋ (c-2) の見積もり 1.49 h（44 遺伝子 ≈ 42,000 細胞×遺伝子 + セントロイド 3 × 90 s）= 1.80 h ≤ 2.1 h。(a) の結果（§A の判定）は (c-2) の設計に影響させていない（遺伝子の選定・細胞数・seed・基準はすべて (a) の実行前に確定済み）
- 2026-09-26 21:35（(c-2) の ISP 完了後、評価の実行前）: (c-2) の細胞単位の評価は `analysis/eval_expanded_cells.py` から EXP-002 `evaluate_controls.py` をそのまま呼ぶ（判定の定義・閾値・a2 の実装は変更しない）。差し替えたのは拡張 run 固有の 3 点のみ = C5 の発現細胞数を新しい start の detection に、a4/d4 の削除位置を本実験の tokenized に、ベクトル指標の start セントロイド名 `aPT_CKD_TLShigh` を `aPT_CKD_TLShighmid` に対応づけ。遺伝子集合の影響を分けるため 34 遺伝子（EXP-003 と同じ集合）と 44 遺伝子（DE マッチ込み）の両方を出す。コードとconfig の sha256 は結果を見る前に `logs/prereg-c2-eval.sha256` に記録

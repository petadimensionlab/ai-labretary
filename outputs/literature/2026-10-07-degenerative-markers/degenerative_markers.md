# Lake 2023（KPMP アトラス）退化した状態（degenerative）のマーカー遺伝子 — 一次情報からの取得記録

- 作成: 2026-10-07 / `literature-scout`
- 用途: HYP-20261007-02 の着手条件（着手前チェックの不合格条件「20 個未満」の回避）
- 出典: Lake BB, Menon R, Winfree S ら、**Nature 619:585–594 (2023)**。
  PMID 37468583、DOI 10.1038/s41586-023-05769-3（一次情報で確認済み。既存調査
  `outputs/literature/2026-09-30-independent-goal-state/report.md` で確認済みのものと同一）

## 1. 依頼された「補足表14」についての重要な訂正

**依頼文は補足表 14 を指定していたが、補足表 14 は退化した状態のリストではない。**
取得した補足資料 ZIP（下記 URL）を展開し、全 37 枚の表の表題を機械的に確認したところ、

- **補足表 12 = 「Conserved degenerative state marker genes」（退化した状態の保存的マーカー遺伝子）**
- 補足表 14 = 「Conserved adaptive stromal (aStr) state marker genes」（**退化ではなく活性化・適応した間質の状態**。aFIB・MYOF 用で、以前の調査で使った表と同じ系統）

**退化した状態の遺伝子リストは補足表 12 にある。** 以下は補足表 12 から取得した内容である。

取得元 URL（既存調査と同じ経路。`www.nature.com` は自動取得を遮断するため static-content 経由）:
`https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-05769-3/MediaObjects/41586_2023_5769_MOESM6_ESM.zip`
展開後のファイル: `2021-07-12390D-Supplementary Table 12.xlsx`（シート名 `S12`）

## 2. 表の構造と読み取り方法

補足表 12 の表題: 「List of degenerative marker genes identified separately for snCv3
(single-nucleus, 10x Chromium v3) and scCv3 (single-cell, 10x Chromium v3). Filter criteria
used to achieve minimal and expanded gene sets are provided.」

**注意: 発現差の数値（log2fc の具体的な値）そのものは、この表には列として出力されていない。**
log2fc > 0.25（最小セット）、log2fc > 0.6（拡張セット）という**しきい値（フィルター条件）**が
表題に明記されており、掲載されている遺伝子は全てこのしきい値を満たした遺伝子である
（= 退化した状態で発現が高い遺伝子と判定された遺伝子）。p 値は Wilcoxon 順位和検定、
調整済み p 値は Benjamini–Hochberg 法。「Number of Cell Types」列は、その遺伝子が
退化した状態のマーカーとして検出された腎臓の細胞型（亜分類）の数を表す
（snCv3 側は腎臓全体で最大 10 種類の細胞型中の該当数）。

4 つの列区分がある（シート内の行 7 がヘッダー、行 8 以降がデータ）:

| 列区分 | 技術 | 条件 | 本調査で取得した遺伝子数 |
|---|---|---|---|
| Minimal Conserved Gene Set（列 B） | snCv3（単一核） | log2fc > 0.25 かつ pval < 0.05 で Ref・AKI・CKD の全群一致、かつ 4 細胞型以上 | **154** |
| Expanded Gene Set（列 C・D、および D 列継続部） | snCv3（単一核） | log2fc > 0.6 かつ pval_adj < 0.05（Ref・AKI・CKD のいずれかで） | **431** |
| Minimal Conserved Gene Set（列 E・F） | scCv3（単一細胞） | 同上の条件を scCv3 データに適用 | **329** |
| Expanded Gene Set（列 G・H） | scCv3（単一細胞） | 同上 | **697** |

4 リストを合わせたユニーク遺伝子数は **882 個**。

**シート内のレイアウトの癖**: snCv3 の Expanded Gene Set は、行 161 までは列 C（遺伝子名）・
列 D（細胞型数）の組で記載され、それ以降（細胞型数が少ない遺伝子）は列 D に遺伝子名、
列 I に細胞型数が記載されている（Excel 上の見た目の列割り付けの都合によるもので、
内容としては 1 つの連続したリストである）。本調査ではこの 2 つの区間を結合して
431 個として扱った。

## 3. 本調査で採用する最有力候補 — snCv3 Minimal Conserved Gene Set（154 遺伝子）

**理由**: 「Ref・AKI・CKD の全群で有意」かつ「4 細胞型以上で検出」という、4 つの区分の中で
最も厳しい条件を満たすリストであり、かつ単一核データ（本プロジェクトが主に扱う単一核 RNA-seq
と技術が一致する）由来である。20 個以上という着手条件を大きく上回る。

頻度（該当する細胞型の数、最大 10）の内訳: 10 件中 2 個、9 件中 4 個、8 件中 17 個、
7 件中 19 個、6 件中 29 個、5 件中 33 個、4 件中 37 個、3 件中 11 個、2 件中 1 個、
不明（Expanded リストに対応する細胞型数の記載なし）1 個（CYSTM1）。

154 個全て（補足表 12 の掲載順）:

ACTB, FTL, ACTG1, APOO, EEF1A1, FTH1, HLA-A, RPL13, RPLP1, RPS9,
ATP1A1, B2M, DHFR, ITM2B, MYL6, PEBP1, PTMA, RPL10, RPL14, RPL29,
RPL7, RPS15, RPS24, RPS27A, RPS3, S100A6, SERF2, TMSB10, TMSB4X, TPT1,
ATP1B1, ATP5F1E, BAIAP2L1, CD24, CD81, DSTN, GAPDH, GSTP1, H3F3B, HSPB1,
IGFBP7, PSAP, RPL11, RPL15, RPL19, RPL3, RPL8, RPS14, RPS2, RPS23,
RPS4X, RPS6, SPP1, TMBIM6, CCNI, EEF2, EIF1, ENO1, FAU, HLA-B,
NDUFA4, PPIA, RPL10A, RPL13A, RPL28, RPL32, RPL4, RPL41, RPLP0, RPLP2,
RPS11, RPS12, RPS19, RPS8, S100A10, SRP14, STEAP1B, UBC, AMBRA1, APLP2,
ATP5F1B, ATP6V0B, CALM1, CDH16, COX4I1, COX7C, CST3, DEFB1, GPX3, HSP90AB1,
HSPA8, LDHB, MTRNR2L12, MYL12B, NORAD, OOEP, PPDPF, RPL18, RPL27A, RPL34,
RPL36, RPL37A, RPL6, RPS16, RPS20, RPS25, SLC25A6, TOMM7, TXNIP, UBA52,
UBB, UQCRB, YWHAB, ATP5F1A, ATP5MG, ATP6V1F, BSG, CANX, CD63, CD74,
COX6B1, CYSTM1, EEF1D, ERP29, FP236383.1, FXYD2, GNAS, GPX4, GRN, HINT1,
HLA-C, HSP90AA1, HSP90B1, JUND, LINGO1, MTRNR2L8, NDUFB7, PFDN5, RACK1, RPL23,
RPL30, RPL35, RPL35A, RPL7A, RPS18, RPS7, RPSA, SDC4, SELENOP, SKP1,
SLC25A3, TMEM213, TMEM59, TXNRD1

**この内容の性質について注意（推測ではなく、表題からそのまま読み取れること）**: このリストは
リボソームタンパク質（RPL/RPS 系列）・ミトコンドリア/酸化的リン酸化関連
（COX, NDUFA4, ATP5 系列）・小胞体ストレス・熱ショックタンパク質（HSP90, HSPA8, HSPB1）・
HLA クラス I（HLA-A/B/C, B2M）など、**細胞型を問わず共通して見られるストレス・分化喪失の
シグネチャーが中心**であり、線維芽細胞（dFIB）や筋線維芽細胞（dM-FIB）に固有の遺伝子という
わけではない。これは表の定義どおりで、「複数の細胞型にまたがる保存的な退化マーカー」である。

## 4. 盲検中の遺伝子についての注意

**SPP1 がこの 154 遺伝子リストに含まれている。** SPP1 は本プロジェクトの盲検対象遺伝子の
1 つである。今回の取得は公表されている補足表をそのまま転記する作業であり、SPP1 についての
摂動・個別解析は行っていない。このリストを使う担当者は、盲検解除の指示が出るまで
SPP1 を摂動・個別解析の対象から外すこと。

## 5. scCv3（単一細胞）側との重複（技術間の再現性の目安）

snCv3 Minimal（154）と scCv3 Minimal（329）の重複は **95 遺伝子**。
scCv3 側は単一細胞という異なる技術であり、単一核データでは測りにくい短い転写産物や
分泌因子を含みうる。重複する 95 個は、技術によらず再現する退化マーカーの候補として扱える。

## 6. これまでの調査との関係

以前の調査（`2026-09-30-independent-goal-state/report.md`）で確認していた dFIB（6 個）・
dM-FIB（10 個）は、**補足表 8（NSForest による細胞型識別用の最小マーカー）**から取得したもので、
「線維芽細胞サブタイプを他と区別するための最小限の組み合わせマーカー」であり、
「退化した状態で発現が高い遺伝子の網羅的なリスト」ではない。今回の補足表 12 は後者の目的に
合致し、かつ 20 個以上という条件を満たす。ただし細胞型非依存（全身共通）のリストである点が
dFIB/dM-FIB 特異的なリストとの違いであり、使用する際はこの違いを踏まえる必要がある。

## 7. 結論（確認できたこと・できなかったこと）

- **確認できたこと**: 補足表 12（Conserved degenerative state marker genes）から、
  snCv3 Minimal 154 個・scCv3 Minimal 329 個・snCv3 Expanded 431 個・scCv3 Expanded 697 個
  （合算ユニーク 882 個）を一次ファイルから取得した。20 個以上という条件は満たす。
- **確認できなかったこと / 依頼と異なっていたこと**: 依頼文が指定した補足表 14 は
  退化した状態のリストではなく、活性化・適応した間質状態（aStr）のリストだった。
  退化した状態のリストは補足表 12 である。
- **この表にない情報**: 遺伝子ごとの log2fc の具体的な数値（しきい値によるフィルターのみ）。
  dFIB・dM-FIB という線維芽細胞サブタイプに限定した発現差リストでもない
  （全細胞型共通の保存的マーカー）。
- 二次情報源は使用していない。全て一次ファイル（補足表原本）からの直接転記である。

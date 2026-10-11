# HYP-20261007-02 段階 1 の着手前: 2 つの遺伝子リスト（一次情報から取得）

- 作成: 2026-10-11 / `literature-scout`
- 用途: HYP-20261007-02「退化した線維芽細胞は細胞老化の発現特徴を持つか」の段階 1 の指標（SenMayo）と、段階 2 の陽性対照（補足表 12）を**測定の前に固定する**ための取得
- 本調査は**リストの取得のみ**。得点は計算していない。盲検中の 5 遺伝子（TP53, CDKN1A, LTBR, SPP1, LCN2）の発現・検出率は一切計算していない（リスト内の名前の有無だけを扱う）

## 0. 成果物

| ファイル | 内容 |
|---|---|
| `senmayo_125.csv` | SenMayo ヒト遺伝子 125 個（`gene_symbol` と、辞書・対象データへの収載有無） |
| `lake2023_table12.csv` | Lake 2023 補足表 12 の 4 区分をまとめた 882 個（`gene_symbol` と区分フラグ・細胞型数） |
| `README.md` | 本ファイル |

## 1. 出典（PubMed で題名・雑誌・年を照合済み。2026-10-11 確認）

### 1.1 SenMayo

- 著者: Saul D ほか（Khosla S が最終著者）
- 表題: **A new gene set identifies senescent cells and predicts senescence-associated pathways across tissues**
- 雑誌: Nature Communications 13, 4827（2022 Aug 16）
- PMID: **35974106** / DOI: 10.1038/s41467-022-32552-1 / PMCID: PMC9381717
- HYP-20261007-02 で以前誤って渡された PMID 36371402 はセレノシステインの化学の論文で無関係（仮説ファイル末尾の訂正と一致）

### 1.2 Lake 2023（KPMP アトラス）

- 著者: Lake BB, Menon R, Winfree S ほか
- 表題: **An atlas of healthy and injured cell states and niches in the human kidney**
- 雑誌: Nature 619, 585–594（2023 Jul）
- PMID: **37468583** / DOI: 10.1038/s41586-023-05769-3

## 2. 取得方法（一次情報）

いずれも出版社の補足ファイル（一次情報）を直接ダウンロードして展開し、セルを直接読んで転記した。要旨・二次情報・記憶による補完は使っていない。

### 2.1 SenMayo（`senmayo_125.csv`）

- 取得元: `https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-32552-1/MediaObjects/41467_2022_32552_MOESM4_ESM.xlsx`
- ファイル: **Supplementary Data 1**（`41467_2022_32552_MOESM4_ESM.xlsx`）
- 展開: シートは `human` と `mouse`。`human` シート（列: `Gene(human)`, `Classification`, `State`, `Reference`）の 2 行目以降が遺伝子記号。**ヘッダーを除いて 125 行、重複なし（125 個のユニーク値）を確認**
- ファイルの sha256: `812770a73c8864fc6b45e1cc7b9deabbb2e50a0fbee564aa1b4a219debb90bcd`
- **補足ファイルに Ensembl ID は無い**（遺伝子記号のみ）。`senmayo_125.csv` の `ensemble_id_geneformer_dict_secondary` 列は、Geneformer の辞書 `models/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl` から機械的に付けた**二次情報**であり、論文の記載ではない

### 2.2 Lake 2023 補足表 12（`lake2023_table12.csv`）

- 取得元: `https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-05769-3/MediaObjects/41586_2023_5769_MOESM6_ESM.zip`
- ZIP 展開後のファイル: **`2021-07-12390D-Supplementary Table 12.xlsx`**
- **シート名: `S12`**
- **表題: 「Supplementary Table 12. Conserved degenerative state marker genes.」**
- ファイルの sha256: `a0ef4883e4591bff769005036334d5e372063fe3ab5c179562389b33f8690be6`
- 補足表 14 は「Conserved adaptive stromal (aStr) state marker genes」であり、**退化ではなく活性化・適応した間質のリスト**。仮説ファイル末尾の訂正どおり、退化した状態のリストは補足表 12 である（表 14 の表題は本調査でも機械的に照合して確認）
- **補足表 12 に Ensembl ID は無い**（遺伝子記号と細胞型数のみ）

#### 補足表 12 の表構造（`S12`）

`S12` シートは 4 区分からなる。ヘッダーは 7 行目、データは 8 行目以降。

| 区分 | 列 | 件数（本調査で実測） |
|---|---|---|
| snCv3 Minimal（単一核） | 列 B | **154** |
| snCv3 Expanded（単一核） | 列 C（行 8–161）＋列 D 継続（行 162 以降） | **431** |
| scCv3 Minimal（単一細胞） | 列 E | **329** |
| scCv3 Expanded（単一細胞） | 列 G | **697** |

- 4 区分の**合算ユニーク数は 882 個**
- snCv3 の Expanded は、行 161 までは列 C（遺伝子）・列 D（細胞型数）の組、行 162 以降は列 D（遺伝子）・列 I（細胞型数）に記載が移る（Excel 上の見た目の都合）。本調査では両区間を結合して 431 個とした
- 最小セットの条件（表題に明記）: 単一核 snCv3 は log2fc > 0.25 かつ pval < 0.05 が Ref・AKI・CKD の全群で一致し、4 細胞型以上。単一細胞 scCv3 も同条件
- **遺伝子ごとの log2fc の具体値はこの表には無い**（しきい値でフィルター済みの掲載のみ）
- このリストは**線維芽細胞サブタイプに限定されない、全細胞型に共通する保存的マーカー**である（リボソームタンパク質・ミトコンドリア/酸化的リン酸化・熱ショックタンパク質など）。dFIB/dM-FIB に固有のリストではない

## 3. 各リストの個数（本調査で実測）

| リスト | 個数 |
|---|---|
| SenMayo（ヒト） | **125** |
| Lake 2023 補足表 12 / snCv3 Minimal | 154 |
| Lake 2023 補足表 12 / snCv3 Expanded | 431 |
| Lake 2023 補足表 12 / scCv3 Minimal | 329 |
| Lake 2023 補足表 12 / scCv3 Expanded | 697 |
| Lake 2023 補足表 12 / 4 区分の合算ユニーク | 882 |

## 4. Geneformer 辞書・対象データに無い遺伝子（本調査で実測）

照合先: Geneformer 辞書 `gene_name_id_dict_gc104M.pkl`（63,675 記号）、A2 の `var`（`feature_name`、30,630 記号）、A1 の `var`（`feature_name`、35,432 記号）。いずれも記号の文字列一致。

| リスト | 個数 | 辞書に無い | A2 に無い | A1 に無い | A2 で使える |
|---|---|---|---|---|---|
| SenMayo | 125 | **0** | **2**（CCL1, CCL3L1） | 1（CCL3L1） | **123** |
| Lake 表 12 / snCv3 Minimal | 154 | 2 | 4 | 4 | 150 |
| Lake 表 12 / snCv3 Expanded | 431 | 8 | 10 | 10 | 421 |
| Lake 表 12 / scCv3 Minimal | 329 | 4 | 5 | 5 | 324 |
| Lake 表 12 / scCv3 Expanded | 697 | 14 | 16 | 16 | 681 |
| Lake 表 12 / 合算ユニーク | 882 | **20** | **24** | **24** | 858 |

- SenMayo は **125 個すべてが Geneformer 辞書にある**。A2 に無いのは CCL1 と CCL3L1 の 2 個（いずれも盲検対象ではない）
- Lake 表 12 の 882 個で辞書に無い 20 個: AC010894.3, AES, ATP5MD, ATP5MPL, C19orf70, C9orf16, FAM96B, FP236383.1, FP671120.1, H1FX, H2AFV, H2AFZ, H3F3A, H3F3B, HIST1H4C, METTL7A, MINOS1, PHB, PLA2G16, RARRES3
  - 多くは**旧遺伝子記号**（例: H2AFV→MACROH2A2、H2AFZ→H2AZ1、H3F3A/B→H3-3A/3B、FAM96B→CIAO2B、METTL7A→METTL7A 改名、AES→TLE5、PHB→PHB1、PLA2G16→PLAAT3、RARRES3→PLAAT4）や、Ensembl 未登録のクローン/ORF 名（AC010894.3, FP236383.1, FP671120.1）。Geneformer 辞書は記号ベースのため旧記号と一致しない
  - Lake 側のみに無い A2（24 個）は上記 20 個に C4orf48, MRPS36, MTRNR2L12, MTRNR2L8 を加えたもの
- A1 と A2 の差は CCL1 のみ（SenMayo では）。A1 も `feature_name` に記号を持つ

## 5. 盲検中の 5 遺伝子について

- プロジェクトの盲検対象: **TP53, CDKN1A, LTBR, SPP1, LCN2**
- **SPP1 は SenMayo 125 にも Lake 表 12 の 882 個にも含まれる**（リストの掲載名としてのみ。発現・検出率は一切計算していない）
- 他の 4 個（TP53, CDKN1A, LTBR, LCN2）は**どちらのリストにも含まれない**（記号一致の確認のみ）
- これらのリストを使う解析では、盲検解除の指示があるまで 5 遺伝子を摂動・個別解析の対象から外すこと

## 6. 循環性の確認（結果を見る前の記録）

- SenMayo 125 と Lake 表 12 の 882 個の**重なりは 10 個**（記号一致。集合演算のみで、発現値は扱っていない）
- 重なる 10 個: **BEX3, CD9, CTSB, IGFBP2, IGFBP3, IGFBP5, IGFBP7, IL32, MIF, SPP1**
- このうち SPP1 は盲検対象（名字のみの記録。値は扱わない）
- この重なりは HYP-20261007-02 §7「循環性」の判断材料として、結果を見る前に記録する。陽性対照（Lake 表 12）と指標（SenMayo）が重なる 10 個の扱いは、計画担当と PI が結果を見る前に決める必要がある

## 7. 確認できたこと / 確認できなかったこと

- **確認できたこと**: SenMayo のヒト 125 個を Supplementary Data 1 から、Lake 2023 補足表 12 の 4 区分（154 / 431 / 329 / 697、合算ユニーク 882）を補足 ZIP の `S12` シートから、いずれも一次ファイルのセルから直接転記した。個数は仮説ファイルの記載（125 個、補足表 12）と一致する
- **確認できなかったこと**: 補足表 12 には遺伝子ごとの log2fc の数値が無い（しきい値のみ）。Ensembl ID はどちらの補足にも無い（Geneformer 辞書による二次的な付与のみ）
- **二次情報を使った箇所**: `senmayo_125.csv` の Ensembl 列（Geneformer 辞書から機械的に付与）。遺伝子記号そのものは一次ファイルからの直接転記

## 8. 注意（HYP-20261007-02 の設計への申し送り）

1. **本調査はリストの取得のみ**。SenMayo の得点・遺伝子ごとの順位表は計算していない
2. Lake 表 12 は**全細胞型共通の保存的マーカー**であり、dFIB/dM-FIB 固有のリストではない。段階 2 の陽性対照に使う際はこの性質を踏まえる必要がある（判断は PI と計画担当）
3. 段階 1 の判定は A2 で行う設計（HYP-20261007-02 §3.2）。SenMayo は A2 で 123 個が使える（T-5 の「25 個未満」には当たらない）
4. 数値はすべて記号の文字列一致による件数であり、発現・検出率ではない

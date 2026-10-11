# SenMayo 遺伝子セット（125 遺伝子・ヒト）— 一次情報からの取得記録

- 作成: 2026-10-07 / `literature-scout`
- 用途: HYP-20261007-02 の着手条件（細胞老化の参照遺伝子セット）

## 1. 依頼されたPMIDについての重要な訂正

依頼文に記載された **PMID 36371402 は、Saul らの SenMayo 論文ではない。**
Europe PMC で実在を確認したところ、PMID 36371402 は
「Site-selective photocatalytic functionalization of peptides and proteins at selenocysteine」
（Dowman LJ ら、Nat Commun 2022；ペプチド化学の論文で老化とは無関係）であった。

**正しい文献は以下である。一次情報（Europe PMC のメタデータ）で確認済み。**

- 著者: Saul D, Kosinsky RL, Atkinson EJ, Doolittle ML, Zhang X, LeBrasseur NK, Pignolo RJ,
  Robbins PD, Niedernhofer LJ, Ikeno Y, Jurk D, Passos JF, Hickson LJ, Xue A, Monroe DG,
  Tchkonia T, Kirkland JL, Farr JN, Khosla S.
- 表題: A new gene set identifies senescent cells and predicts senescence-associated pathways across tissues
- 誌名・巻号: Nature Communications 13, 4827 (2022)
- **正しい PMID: 35974106**（依頼文の 36371402 ではない）
- PMCID: PMC9381717
- DOI: 10.1038/s41467-022-32552-1

以降はこの正しい文献（PMID 35974106）に基づく。

## 2. 遺伝子リストの取得経路（一次情報）

- Europe PMC の全文 XML（`PMC9381717/fullTextXML`）から補足資料一覧を確認
- 本文表ではなく、補足資料 **Supplementary Data 1**
  （ファイル名 `41467_2022_32552_MOESM4_ESM.xlsx`）を取得
  - 取得元 URL:
    `https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-32552-1/MediaObjects/41467_2022_32552_MOESM4_ESM.xlsx`
- このファイルには `human` シートと `mouse` シートがある。
  `human` シート（列: `Gene(human)`, `Classification`, `State`, `Reference`）に
  **ヘッダーを除いて 125 行**あり、遺伝子記号の重複はない（`125` 個のユニーク値を確認）。
  これが SenMayo のヒト遺伝子セットである。
- 参考として `Supplementary Data 2`（`MOESM5_ESM.xlsx`）も確認した。
  こちらは SenMayo の PPI ネットワーク統計（ノード数 125・エッジ数 1,847 等）であり、
  Fig. 1e–f に対応する。ノード数 125 は本文の図と一致し、取得した遺伝子数と整合する。

## 3. 125 個の遺伝子記号（一次ファイルからそのまま転記）

ACVR1B, ANG, ANGPT1, ANGPTL4, AREG, AXL, BEX3, BMP2, BMP6, C3,
CCL1, CCL13, CCL16, CCL2, CCL20, CCL24, CCL26, CCL3, CCL3L1, CCL4,
CCL5, CCL7, CCL8, CD55, CD9, CSF1, CSF2, CSF2RB, CST4, CTNNB1,
CTSB, CXCL1, CXCL10, CXCL12, CXCL16, CXCL2, CXCL3, CXCL8, CXCR2, DKK1,
EDN1, EGF, EGFR, EREG, ESM1, ETS2, FAS, FGF1, FGF2, FGF7,
GDF15, GEM, GMFG, HGF, HMGB1, ICAM1, ICAM3, IGF1, IGFBP1, IGFBP2,
IGFBP3, IGFBP4, IGFBP5, IGFBP6, IGFBP7, IL10, IL13, IL15, IL18, IL1A,
IL1B, IL2, IL32, IL6, IL6ST, IL7, INHA, IQGAP2, ITGA2, ITPKA,
JUN, KITLG, LCP1, MIF, MMP1, MMP10, MMP12, MMP13, MMP14, MMP2,
MMP3, MMP9, NAP1L4, NRG1, PAPPA, PECAM1, PGF, PIGF, PLAT, PLAU,
PLAUR, PTBP1, PTGER2, PTGES, RPS6KA5, SCAMP4, SELPLG, SEMA3F, SERPINB4, SERPINE1,
SERPINE2, SPP1, SPX, TIMP2, TNF, TNFRSF10C, TNFRSF11B, TNFRSF1A, TNFRSF1B, TUBGCP2,
VEGFA, VEGFC, VGF, WNT16, WNT2

（計 125 個。1 列の遺伝子記号は `senmayo_125.tsv` を参照）

## 4. 盲検中の遺伝子についての注意

**SPP1 がこのリストに含まれている。** SPP1 はプロジェクトの盲検対象遺伝子の 1 つである。
今回の取得は「公表されている遺伝子セットをそのまま転記する」作業であり、SPP1 についての
摂動・解析は一切行っていない。このリストを今後の解析に使う担当者は、
盲検解除の指示が出るまで SPP1 を摂動・個別解析の対象から外す必要がある。

## 5. Ensembl ID について（二次情報源・注意が必要）

- 一次論文の補足資料には Ensembl ID は含まれていない（遺伝子記号のみ）。
- そのため Ensembl ID は **Ensembl REST API**（`rest.ensembl.org/lookup/symbol/homo_sapiens`、
  2026-10-07 時点のデータベース内容）を用いて**別途、機械的に付与した二次情報**である。
  一次論文の記載ではない点に注意。
- 125 個中 **124 個は一意の Ensembl Gene ID が得られた**。
- **CCL3L1 のみ一意に定まらなかった。** CCL3L1 はコピー数多型を持つ遺伝子で、
  Ensembl 上で CCL3L1/CCL3L3 系統の複数の遺伝子 ID
  （ENSG00000276085, ENSG00000277336, ENSG00000277768 等）が重複して存在するため。
  `senmayo_125.tsv` にはこれら候補をすべて記載し、「一意に確定せず」と注記した。
- 詳細は `senmayo_125.tsv` の 2 列目を参照。

## 6. 結論（この調査の範囲内で確認できたこと・できなかったこと）

- **確認できたこと**: SenMayo のヒト遺伝子 125 個を一次ファイル（Supplementary Data 1）から
  過不足なく取得した。重複なし。
- **確認できたこと（訂正）**: 依頼された PMID 36371402 は誤りで、正しくは PMID 35974106 である。
- **二次情報源を使った箇所**: Ensembl ID の付与（Ensembl REST API）。遺伝子記号そのものは
  一次ファイルからの直接転記であり、記憶による補完は行っていない。

# 単一細胞解析の用語照合：library size / sequencing depth / detection rate

- 作成: literature-scout / 2026-10-11
- 目的: 実験計画書・レポートで使っている内部名（`n_umi`、`n_genes`、総称としての「深さ」、「深さをそろえる（depth matching）」、「検出率」）が、
  査読論文・教科書で実際に何と呼ばれているかを一次資料で確認する。
- 方針: 事実のみ。仮説・解釈は書かない。確認できないものは「無い／未確認」と明記する。
- 用語の対応先: `docs/用語集.md` §1・§3。

---

## 0. 結論（先に要点）

| 内部名 | 一次資料で確認できた標準用語 | 判定 |
|---|---|---|
| `n_umi` | **total UMI count per cell** / **number of UMIs per cell** / 文脈により **library size** | 標準用語あり |
| `n_genes` | **number of genes detected per cell**（Seurat の `nFeature_RNA`、scater の `detected`） | 標準用語あり |
| 1 と 2 の総称としての「深さ」 | **library size**（総カウント数・総 UMI 数）と **sequencing depth / read depth**（生リード数）は**別の量**。総称する単一の標準語は無い | **「深さ」という総称は無い（別概念を混ぜている）** |
| 「深さをそろえる」`depth matching` | **downsampling** / **subsampling** / **rarefaction**（coverage equalization）。`depth matching` という定義済み用語は見つからない | **`depth matching` は標準用語ではない** |
| snRNA-seq の「検出される割合の低さ」 | **gene detection rate**（観測された割合）と **dropout rate**（未検出の条件付き確率）は別物。**detection sensitivity** は手法の感度 | 標準用語あり（ただし dropout は定義が揺れる） |

**最も重要な事実**: `docs/用語集.md` の `n_umi`＝「総 UMI 数」、`n_genes`＝「検出された遺伝子数」の対応は、文献の標準用語と一致している。
問題は **1 と 2 をまとめる総称を「深さ（depth）」としている点** と、**「depth matching / 深さマッチ」という操作名** である。前者は 2 つの別の量を 1 語で呼んでおり、後者は文献に定義済みの用語が無い。

---

## 1. `n_umi`（各細胞の総 UMI 数）

| 内部名 | 学術用語（英語） | 標準的な日本語訳 | 出典 |
|---|---|---|---|
| `n_umi` | **total UMI count per cell** / **number of UMIs per cell** / **UMI counts per cell** | **総 UMI 数**、**細胞あたりの UMI 数**（日本語論文でも「UMI カウント数」「UMI 数」と表記） | 下記 [S1][S5][S6] |
| （同義・文脈依存） | **library size**（＝その細胞の全カウントの総和。UMI でも生リードでもあり得る） | **ライブラリサイズ**（カタカナ表記が多い。訳語を付けるなら「総カウント数」） | [S2][S3][S10] |

- [S1] scater の `perCellQCMetrics()` は各細胞について `sum`（その細胞の総カウント）を出す。
- [S5] 10x Genomics は「Median UMI counts per cell」と呼ぶ。
- [S6] Hafemeister & Satija は「total UMI counts / cell」を図のラベルに使う。
- [S2][S3] OSCA は library size を「各細胞で全遺伝子にまたがるカウントの総和」と定義。
- [S10] Jiang らは「cell library size（各細胞の総リード数または総 UMI 数）」と定義。

**注意（重要）**: 総 UMI 数は library size とも呼ばれるが、**UMI 数は生リード数ではない**。UMI は同一分子の重複を 1 にまとめた数で、生リード数より少ない（リードの一部）。

---

## 2. `n_genes`（各細胞で検出された遺伝子数）

| 内部名 | 学術用語（英語） | 標準的な日本語訳 | 出典 |
|---|---|---|---|
| `n_genes` | **number of genes detected per cell** / **genes detected per cell** / **number of detected genes** | **検出遺伝子数**、**検出された遺伝子数** | [S1][S5][S9] |
| （ツール内の列名） | Seurat の **`nFeature_RNA`**、scater の **`detected`** | ― | [S1][S5] |
| （同義・QC 文脈） | **gene complexity** / **transcriptional diversity**（検出遺伝子数を QC 指標として呼ぶときの別名） | **遺伝子の複雑さ**、**転写の多様性**（訳語は定着していない） | [S9] |

- [S1] scater の `detected` 列は「検出された遺伝子の数」。
- [S5] 10x は「Median genes per cell」。「検出（detection）は少なくとも 1 UMI が存在すること」と定義。
- [S9] Subramanian ら（ddqc）は QC 指標として「gene complexity」「transcriptional diversity」を使う。

**注意**: 検出遺伝子数は「検出」の閾値（≥1 UMI か、≥3 UMI か等）に依存する。10x は ≥1 UMI と定義している [S5]。

---

## 3. 1 と 2 をまとめるときの学術的な総称（最重要）

**「深さ」という単一の総称は文献に無い。** 1 と 2 は別々の量であり、文献では次のように区別される。

| 概念（英語） | 何を指すか | 単位 | 標準的な日本語訳 | 出典 |
|---|---|---|---|---|
| **library size**（= **count depth** とも） | その細胞に割り当てられた**カウント（UMI、または UMI を使わない手法ではリード）の総和** | 分子数（counts） | **ライブラリサイズ**、**総カウント数** | [S2][S3][S4][S10] |
| **count depth** | バーコードあたりのカウント数（library size とほぼ同義） | counts | （訳語は定着せず。カウント数） | [S4] |
| **sequencing depth**（= **read depth**） | その細胞に費やした**生のシーケンスリード数**。**カウント数とは違う** | リード数（reads） | **シークエンス深度**、**シーケンシング深度**（日本語レビューでも「シークエンス深度，sequence depth」と併記） | [S7][S8][S17] |

- [S7] Haque らの用語集は **sequencing depth を「単一サンプルに費やしたシーケンス能力の量。例えば 1 細胞あたりの生リード数として報告される」** と定義する。
- [S8] Zhang らは **cell-wise sequencing depth（＝1 細胞あたりに期待されるリード数）** と定義する。
- [S2][S3] library size は「カウントの総和」＝分子数側の量。
- [S6] Hafemeister & Satija は「total UMI counts / cell」を「sequencing depth」とラベルしており、**総 UMI 数を sequencing depth と呼ぶ用法も存在する**。ただしこれは UMI（分子数）を深度（リード数）と同一視する用法で、[S7][S8] の定義とは**食い違う**。同じ尺度で混ぜてはいけない。

**したがって**:
- 内部の 1（総 UMI 数）と 2（検出遺伝子数）を 1 語でまとめる学術的な総称は**存在しない**。
- まとめたいときの候補は「**library size と gene complexity**」のように**2 つ並べる**か、両者に共通する上位概念として「**cell-specific technical bias**（細胞ごとの技術的な偏り）」[S2][S3] を使う。
- 「深さ」という語でまとめると、library size（分子数）と sequencing depth（リード数）という**別の量**に読めてしまう。

---

## 4. 「比較する群の間で 1 や 2 の分布を 1 対 1 で重ねてから集計する」操作

### 4.1 文献にある標準語

| 操作（英語） | 内容 | 標準的な日本語訳 | 出典 |
|---|---|---|---|
| **downsampling** / **subsampling** | 各細胞（または各バッチ）からカウント／リードをランダムに間引いて、共通の総数にそろえる | **ダウンサンプリング**（間引き） | [S4][S19][S20] |
| **rarefaction** | 共通の値まで間引く操作を多数回繰り返し、指標を平均する。微生物生態学由来 | **レアファクション**（希薄化。日本語訳はほぼ無い） | [S21] |
| **coverage equalization** | バッチ間の被覆（1 細胞あたりの総カウント）を等しくする。物理的には cDNA ライブラリの濃度をそろえる操作を指す | **カバレッジの均等化** | [S19]（理論面は [S22]） |
| **full-quantile normalization / distribution matching** | 群間で**分布全体の分位を一致させる**操作（バルク RNA-seq・マイクロアレイ由来） | **分位正規化**、**分布の一致** | [S22] |

- [S19] Bioconductor の scuttle の `downsampleBatches()` は「全バッチの 1 細胞あたり平均総カウントを同じにする。これは cellranger aggr のダウンサンプリング機能を模倣する」と説明する。
- [S4] Luecken & Theis は downsampling を「全細胞をあらかじめ指定したカウント数以下にそろえるため、リードまたはカウントをランダムにサンプリングする操作」と定義。
- [S21] Schloss は rarefaction を「最小 library size を選び、それを下回るサンプルを除き、残りを復元抽出なしで共通サイズに間引く操作（を多数回繰り返す）」と定義。

### 4.2 `depth matching` という用語について

- **査読論文で定義された `depth matching` という用語は見つからなかった。** 単一細胞解析の標準操作は **downsampling / subsampling / rarefaction** [S4][S19][S20][S21]。
- 形容詞的に **"depth-matched"** を使う例はある（例: 2013 年の単一細胞手法比較で「sequence-depth matched」、ある校正査読の「depth-matched analyses」）。ただし定義済みの用語ではなく、その場限りの記述。
- 内部の操作（「分布を **1 対 1 で重ねる**」）に特に近い既存語は:
  - 群全体を共通の値まで間引く **downsampling** [S19]
  - 分布全体の分位を一致させる **full-quantile normalization** [S22]（分野はバルク）
  - 共変量で個体を 1 対 1 に組にする **matching（マッチング）**（統計一般の用語で、単一細胞解析の用語ではない）

**判定**: 「1 対 1 で重ねる」という**厳密な対応付け**を指す単一細胞解析の標準用語は**無い**。内部で `depth matching` を使い続けるしかない場合は、**造語（プロジェクト内の操作名）であることと、何をそろえたか（総 UMI 数か検出遺伝子数か）を毎回明記する**必要がある。

---

## 5. snRNA-seq に特有の「ある遺伝子が検出される割合」の低さ

2 つの用語は**定義が異なる**。

| 英語 | 定義 | 標準的な日本語訳 | 出典 |
|---|---|---|---|
| **gene detection rate** | ある遺伝子が検出された細胞の割合（観測された割合）。「検出」＝カウント > 0（または閾値超え） | **遺伝子検出率** | [S15][S16][S17] |
| **detection sensitivity** | 手法・プロトコルが転写産物を検出する能力（感度）。率ではなく性能 | **検出感度** | [S12][S17] |
| **dropout rate** | ある遺伝子が、発現しているのに**ある細胞で検出されない**確率（未検出の条件付き確率）。ゼロの割合として測られることが多い | **ドロップアウト率**（訳語は定着していない） | [S12] |
| **cellular detection rate (CDR)** | **1 細胞内で検出された遺伝子の割合**（細胞単位の量。上の gene detection rate とは向きが違う） | **細胞内検出率**（訳語は定着していない） | [S11] |
| **data sparsity** | データ全体にゼロが多いこと | **データの疎性／スパース性** | [S10] |
| **zero inflation** | 特定の分布（Poisson / NB）を仮定したときに、**期待を超えるゼロの割合**（統計モデルの概念） | **ゼロ過剰**、**ゼロインフレーション** | [S10][S14] |

- [S11] Finak ら（MAST）が **CDR（cellular detection rate）＝各細胞で検出表現される遺伝子の割合** を導入した。
- [S12] Yamawaki らは **dropout rate を「ある遺伝子がある細胞で検出されない条件付き確率」** と定義し、**gene detection sensitivity** と区別している。
- [S15] Mereu ら、[S16] Selewa ら（Drop-seq vs DroNc-seq）、[S17] Salmen ら（VASA-seq）は **gene detection rate** を使う。Selewa らは **snRNA-seq の gene detection rate が scRNA-seq より有意に低い**と報告。
- [S10] Jiang らは **「dropout」「excess zeros」「zero inflation」の使い分けが文献で混乱している**と指摘し、それ自体では曖昧な「dropout」を使わず「zero inflation」を使うと明言する。
- [S13] Sarkar & Stephens も、**「dropout」「missing data」という用語自体が多くの場合役に立たず混乱を招く**と論じている。
- [S14] Svensson はドロップレット scRNA-seq はゼロ過剰ではない（zero-inflated ではない）と示した。
- [S18] Chamberlin らは **snRNA-seq のデータの疎性（sparsity）は核が細胞質 mRNA をほとんど含まないため**であり、イントロンリードの取り込み等で扱われる、と説明する。

**snRNA-seq に特有の呼び方**として文献が使う表現:
- **lower transcript yield from nuclei**（核からの転写産物収量が低い）[S18]
- **excessive data sparsity**（過剰なデータの疎性）[S18]
- **lower gene detection rate / lower detection sensitivity**（低い遺伝子検出率・検出感度）[S16]
- 単一の専用語は無く、「**低い detection rate**」「**高い dropout**」として一般的な語を流用している。

**detection rate と dropout の違い（要点）**: detection rate は**観測された「検出された割合」**（高いほど良い）、dropout rate は**「未検出の確率」**（高いほど悪い）で、**互いに補数の関係**にありながら、dropout は「本当は発現しているのに検出されない」という含意を持つ点が異なる [S10][S12][S13]。

---

## 6. 用語集への反映案（提案。判断は PI）

`docs/用語集.md` の該当 2 行について、次の修正を提案する（**適用はしていない**）。

- 現行 §3「読み取りの多さ」→ 「**library size（総カウント数・総 UMI 数）と sequencing depth（生リード数）は別の量**」と明記し、「深さ」という総称を使わない。
- 現行 §3「読み取りの多さをそろえる（depth matching, 深さマッチ）」→ 英語標準語は **downsampling / subsampling / rarefaction**（[S4][S19][S21]）。`depth matching` は**プロジェクト内の操作名（造語）**であり、そろえた量を必ず併記する、と注記。
- §1 の単一核 RNA-seq の行に、**gene detection rate（遺伝子検出率）** と **dropout rate（ドロップアウト率）は別**、と追記。

---

## 7. 出典一覧（実在確認済み）

- [S1] McCarthy DJ, Campbell KR, Lun ATL, Wills QF. Scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R. *Bioinformatics* 2017;33(8):1179–1186. DOI: 10.1093/bioinformatics/btw777（Crossref で確認）
- [S2] Amezquita RA, Lun ATL, Becht E, et al. Orchestrating single-cell analysis with Bioconductor. *Nat Methods* 2020;17:137–145. DOI: 10.1038/s41592-019-0654-x（Crossref で確認）
- [S3] OSCA（Orchestrating Single-Cell Analysis with Bioconductor）, Ch.1 Quality Control / Ch.2 Normalization. https://bioconductor.org/books/release/OSCA.basic/normalization.html / https://bioconductor.org/books/release/OSCA.basic/quality-control.html （オンライン教科書。`sum`＝総カウント、`detected`＝検出遺伝子数、library size の定義）
- [S4] Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. *Mol Syst Biol* 2019;15(6):e8746. DOI: 10.15252/msb.20188746（Crossref で確認）
- [S5] 10x Genomics. *Interpreting Cell Ranger Web Summary Files for Single Cell Gene Expression Assay*. Technical Note CG000329 RevA. https://cdn.10xgenomics.com/image/upload/v1660261286/support-documents/CG000329_TechnicalNote_InterpretingCellRangerWebSummaryFiles_RevA.pdf （**製造者テクニカルノート。査読なし**）
- [S6] Hafemeister C, Satija R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. *Genome Biol* 2019;20:296. DOI: 10.1186/s13059-019-1874-1（Crossref で確認）
- [S7] Haque A, Engel JA, Teichmann SA, Lönnberg T. A practical guide to single-cell RNA-sequencing for biomedical research and clinical applications. *Genome Med* 2017;9:75. DOI: 10.1186/s13073-017-0467-4（Crossref で確認。用語集で sequencing depth を定義）
- [S8] Zhang MJ, Ntranos V, Tse D. Determining sequencing depth in a single-cell RNA-seq experiment. *Nat Commun* 2020;11:810. DOI: 10.1038/s41467-020-14482-y（Crossref で確認）
- [S9] Subramanian A, Alperovich M, Yang Y, Li B. Biology-inspired data-driven quality control for scientific discovery in single-cell transcriptomics. *Genome Biol* 2022;23:267. DOI: 10.1186/s13059-022-02820-w. PMID: 36575523
- [S10] Jiang R, Sun T, Song D, Li JJ. Statistics or biology: the zero-inflation controversy about scRNA-seq data. *Genome Biol* 2022;23:31. DOI: 10.1186/s13059-022-02601-5. PMID: 35063006
- [S11] Finak G, McDavid A, Yajima M, et al. MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data. *Genome Biol* 2015;16:278. DOI: 10.1186/s13059-015-0844-5. PMID: 26653891（Crossref で確認。CDR の定義）
- [S12] Yamawaki TM, Lu DR, Ellwanger DC, et al. Systematic comparison of high-throughput single-cell RNA-seq methods for immune cell profiling. *BMC Genomics* 2021;22:66. DOI: 10.1186/s12864-020-07358-4（Crossref で確認。dropout rate の定義）
- [S13] Sarkar A, Stephens M. Separating measurement and expression models clarifies confusion in single-cell RNA sequencing analysis. *Nat Genet* 2021;53:760–768. DOI: 10.1038/s41588-021-00873-4（Crossref で確認）
- [S14] Svensson V. Droplet scRNA-seq is not zero-inflated. *Nat Biotechnol* 2020;38:147–150. DOI: 10.1038/s41587-019-0379-5（Crossref で確認）
- [S15] Mereu E, Lafzi A, Moutinho C, et al. Benchmarking single-cell RNA-sequencing protocols for cell atlas projects. *Nat Biotechnol* 2020;38:747–755. DOI: 10.1038/s41587-020-0469-4（Crossref で確認。gene detection rate）
- [S16] Selewa A, Dohn R, Eckart H, et al. Systematic Comparison of High-throughput Single-Cell and Single-Nucleus Transcriptomes during Cardiomyocyte Differentiation. *Sci Rep* 2020;10:1454. DOI: 10.1038/s41598-020-58327-6. PMID: 31996743（Crossref で確認。snRNA の gene detection rate が低い）
- [S17] Salmen F, De Jonghe J, Kaminski TS, et al. High-throughput total RNA sequencing in single cells using VASA-seq. *Nat Biotechnol* 2022;40:1780–1793. DOI: 10.1038/s41587-022-01361-8（DOI は出版社 URL から。**Crossref 未確認**）
- [S18] Chamberlin JT, Lee Y, Marth GT, Quinlan AR. Differences in molecular sampling and data processing explain variation among single-cell and single-nucleus RNA-seq experiments. *Genome Res* 2024;34(2):179–188. DOI: 10.1101/gr.278253.123. PMID: 38355308（URL と PubMed 由来。**DOI は Crossref 未確認**）
- [S19] scuttle（Bioconductor）. `downsampleMatrix` / `downsampleBatches` のドキュメント. https://rdrr.io/bioc/scuttle/man/downsampleBatches.html（"Downsample batches to equal coverage … mimics the downsampling functionality of cellranger aggr"）
- [S20] Luecken & Theis 2019（=[S4]）内の downsampling の定義。
- [S21] Schloss PD. Rarefaction is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses. *mSphere* 2024;9(3):e00354-23. DOI: 10.1128/msphere.00354-23（Crossref で確認。rarefaction の定義）
- [S22] Cole MB, Risso D, Wagner A, et al. Performance Assessment and Selection of Normalization Procedures for Single-Cell RNA-seq. *Cell Syst* 2019. プレプリント PDF: https://research.unipd.it/retrieve/e14fb26c-2a70-3de1-e053-1705fe0ac030/Cole2019preprint.pdf（full-quantile normalization の定義。**DOI・PMID 未確認**）

### 日本語の出典（訳語の根拠）
- 濱田（ほか）「シングルセル解析の動向と展望」*JSBi Bioinformatics Review* 2022;3(2). https://www.jstage.jst.go.jp/article/jsbibr/3/2/3_jsbibr.2022.4/_html/-char/ja （「sequence depth」「リード数」「UMI 1 カウント数」「検出感度」の併記）
- 顕微鏡誌 単一細胞トランスクリプトーム総説 2024;59(2):67. https://www.jstage.jst.go.jp/article/kenbikyo/59/2/59_67/_pdf/-char/ja （「一細胞あたりのリード数（シークエンス深度，sequence depth）」「検出遺伝子数」「検出率」）
- Basepair ブログ（**査読なし**）. https://basepairtech.jp/blog/2553/ （「Median genes/features per cell＝細胞あたり検出された遺伝子」「Median UMI count per cell＝細胞あたり検出された転写産物」）

---

## 8. 限界・確認できなかったこと

1. **`depth matching` は査読文献で定義された用語として見つからなかった。** 見つかったのは downsampling / subsampling / rarefaction / coverage equalization。形容詞の "depth-matched" は散発的に使われる（定義なし）。
2. **1 と 2 をまとめる単一の学術的総称は無い。** library size（分子数側）と sequencing depth（リード数側）は別の量で、まとめる語は確認できなかった。
3. **sequencing depth の用法は文献内で揺れる。** Hafemeister & Satija は total UMI counts を sequencing depth と呼ぶ [S6] が、Haque [S7]・Zhang [S8] は生リード数と定義する。**同じ尺度で混ぜないよう注意**。
4. [S17][S18] の DOI、[S22] の DOI・PMID は Crossref で実在確認できていない（未確認と明記）。
5. 日本語の「標準訳語」は、定着した訳が無く英語のまま使われる語が多い（library size→ライブラリサイズ、rarefaction→レアファクション等）。日本語一次資料はレビュー記事で、訳語を確定する公的規範は見つからなかった。

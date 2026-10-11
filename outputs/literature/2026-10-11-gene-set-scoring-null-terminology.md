# 遺伝子セット解析の用語照合：内部名（偽セット・較正・検出核・セット平均）と一次資料の標準用語

- 作成: literature-scout / 2026-10-11
- 目的: 単一細胞 RNA-seq のデータ解析レポートで使っている内部語（偽セット、較正、検出核、セット平均）が、
  査読付き論文・教科書・公式マニュアルで実際にどう呼ばれているかを一次資料で確認する。
- 方針: 事実のみ。仮説・解釈は書かない。確認できないものは「無い／未確認」と明記し、作り出さない。
- 関連: `docs/用語集.md` §1（遺伝子検出率）、§3.5（マーカーの検出率）、§4.5（主判定）。
- 一次資料の優先: 査読論文・学会誌・公式マニュアル。ブログ・GitHub issue は補助にとどめ、その旨を書く。

---

## 0. 結論（先に要点）

| 内部名 | 一次資料で確認できた標準用語（英語） | 定着した日本語訳 | 判定 |
|---|---|---|---|
| **偽セット** | 単独の対応語は無い。近い概念は **random gene set (size-matched)** / **gene-set permutation** / **competitive null**。復元ありで引くなら統計学では **bootstrap resampling** | 訳語は無い（英語のまま／説明的に書く） | **「sham set」は学術語として存在しない**。作成物は「ランダム遺伝子セット」等と説明的に書く |
| **較正** | **calibration**（well-calibrated p-values / p-value calibration）。手順の実体は **type I error rate control** / **uniformity of p-values** | **較正**（日本語の統計文献で定着。確率較正・riskk 較正・ブートストラップ較正） | **較正は calibration の訳として定着している**。ただし意味は「モデルの予測確率の整合」で、本件は「名目水準どおりの第一種の過誤率」 |
| **検出核** | **number of cells/nuclei in which the gene is detected** / ツール列名 `n_cells`（per gene） | 単独の定着訳は無い。「その遺伝子が検出された核の数」と説明的に書く | 割合は **gene detection rate（遺伝子検出率）**。**数（count）の専用語は無い** |
| **セット平均** | **mean gene set score** / **set-averaged score** / **average score** | 単独の定着訳は無い。「遺伝子セットの平均スコア」と書く | **「セット平均」という単独語は学術用語として確認できない** |

**最も重要な事実**:
- 内部の「偽セット」は、**英語の学術文献に「sham set」という語は無い**。近い概念は GSEA の **gene-set permutation（サイズを揃えたランダム遺伝子セット）**と、Goeman & Bühlmann の **competitive null（競合帰無仮説）**である。
- 「較正」は**日本語で定着した訳語**（校閲付き日本語統計文献で確認）。ただし本件が測っているのは「名目水準どおりの第一種の過誤率」であり、英語では **calibration** よりも **type I error rate control** と呼ばれることが多い。
- 「検出核」「セット平均」は、**単独では学術用語として確認できない**。いずれも「何を数えた・何を平均したか」を添えて説明的に書くのが適切。

---

## 1. 「偽セット」— 対照プールから復元ありで引いて帰無分布を作る

### 1.1 手順の分解

内部の手順は 3 つに分解できる。それぞれに対応する学術用語が異なる。

| 手順の要素 | 内容 | 対応する学術用語（英語） | 出典 |
|---|---|---|---|
| (a) 帰無仮説の種類 | 観測セットの平均を、**同サイズのランダムな遺伝子集合の平均**と比べる | **competitive null hypothesis（競合帰無仮説）**、「遺伝子をサンプリング単位とする **gene-sampling model**」 | Goeman & Bühlmann 2007 [S3]、Tian 2005 [S4] |
| (b) ランダム集合の作り方（GSEA の用語） | サイズを揃えたランダム遺伝子セットを作り、そのスコアを計算して帰無分布にする | **gene_set permutation**（GSEA の permutation type の一つ）。singscore はこれを **permutation test** と呼ぶ | GSEA 公式マニュアル [S2]、Subramanian 2005 [S1]、singscore 公式マニュアル [S7] |
| (c) 引くときの復元 | **復元あり**（重複を許す）で引く | **bootstrap resampling / sampling with replacement**。復元なしなら permutation、処置の割り付けをランダムに変えるなら randomization test | Efron 1979 [S5] |

### 1.2 permutation / bootstrap / resampling の区別（一次資料の定義）

| 英語 | 復元 | 定義（一次資料の記述） | 出典 |
|---|---|---|---|
| **permutation test** | **復元なし**（without replacement） | 「並べ替えの複製は復元なしで計算される」 | Efron 1979 [S5]、Good 1994 [S8] |
| **randomization test** | 復元なし（処置ラベルをランダムに割り付ける） | 「古典的な順列検定・ランダム化検定でなされるように、結果を処置群へ割り当て直す」 | [S5] の系譜、関連文献 [S9] |
| **bootstrap** | **復元あり**（with replacement） | Efron の原文: 「X\* の値は集合 {x1,…,xn} から**復元あり**で選ばれるので、これは順列分布ではない」 | Efron 1979 [S5] |
| **Monte Carlo test** | どちらもあり得る | 「順列分布を大きな Monte Carlo 標本で近似するのが普通」 | Larson & Owen 2015 [S10] |

**判定**:
- 内部の手順は**復元あり**なので、厳密な統計学の分類では **bootstrap（復元抽出による帰無分布）／Monte Carlo resampling** に当たる。**permutation test ではない**（permutation は復元なし）。
- ただし、遺伝子セット解析の**分野の慣行**では、サイズを揃えたランダム遺伝子セットで帰無分布を作る操作を **gene-set permutation / permutation test** と呼ぶ（GSEA 公式マニュアル [S2]、singscore 公式マニュアル [S7]）。厳密名と分野の通称が食い違う。
- **どちらの名前を使うにせよ、「復元ありで引いた」ことを必ず併記する**。書かないと、復元なしの permutation と読まれる。

### 1.3 「sham set」という語について

- **査読論文・教科書で「sham set」という学術用語は確認できなかった。**
  - Crossref で "sham gene set enrichment null distribution" を検索しても、ヒットするのは図のキャプションのみで、用語としての使用は無い。
  - 一般語の "sham" は医学で「シャム手術・シャム刺激（実際には効果のない対照処置）」を指す（例: デジタル治療の sham control、[S11]）。**遺伝子セット解析の用語ではない。**
- **したがって「偽セット」を英語に直すとき "sham set" を使うべきではない。** 使うと医学の対照処置と誤読される。
- 推奨する書き方（学術用語）:
  - 「サイズを揃えた**ランダム遺伝子セット**（random gene sets of matched size）」
  - 「**gene-set permutation**（遺伝子ラベルをランダム化する帰無分布）」
  - GSEA の用語を使うなら「**gene_set permutation**」（[S2]）

---

## 2. 「較正」— 名目水準どおりに第一種の過誤が起きることの確認

### 2.1 英語

| 内部の操作 | 英語 | 定義・出典 |
|---|---|---|
| 対照プール自身を検定し、偶然 99% 点を超える割合が片側 1% 付近に収まることの確認 | **calibration（較正）**／**type I error rate control（第一種の過誤率の制御）**／**uniformity of p-values（p値の一様性）** | p値が帰無仮説の下で一様分布になることが「較正されている（well-calibrated）」の意味。Wasserstein 2017 [S12]、Berger et al. 2001 [S13] |
| 同（多重検定の文脈） | **p-value calibration** | 「p値の較正は…帰無仮説の下での一様分布 Uniform(0,1) を基礎とする」[S12] |

- **calibration（較正）** は、モデルが「予測確率 p の集団で実際の割合が p になる」ことを指す語（[S13] ほか）。本件（名目水準での第一種の過誤率）も同じ「名目と実測がずれない」性質なので、calibration と呼べる。
- より直接的な英語は **type I error rate control / empirical type I error rate**。「較正」だけだと「確率予測モデルの較正」と混同されるので、**英語を併記するときは `calibration (type I error rate control)` のように補う**のが安全。

### 2.2 日本語「較正」は定着しているか → **定着している**

査読付き・学会の日本語一次資料で「較正」の使用を確認した。

| 出典 | 分野 | 「較正」の使われ方 |
|---|---|---|
| 日本計量生物学会 ニュースレター 151号 [J1] | 生物統計 | 「頻度流の**第 I 種過誤確率に対して較正**する必要がない」「第一種過誤確率の漸近的な一様制御」 |
| 計算統計シンポジウム 講演要旨 [J2] | 統計計算 | 「ロジスティック回帰モデルにおける**ブートストラップ較正**と経験尤度」 |
| 計量生物学 41(1)「イベント予測モデルの評価指標」[J3] | 医学統計 | 「**リスク較正**（risk calibration）」「較正プロット」「大域的較正」「較正スロープ」 |
| 医用画像情報学会雑誌 40(4)「機械学習における予測の不確実性定量化」[J4] | 機械学習 | 「**確率較正**（probability calibration）」「再較正（re-calibration）」 |
| 人工知能学会 全国大会 2023「決定損失の期待値と分散を用いた較正」[J5] | 機械学習 | 表題に「較正」。「較正手法（calibration）」 |

**判定**: 「較正」は**日本語の統計・機械学習文献で calibration の訳語として定着している**。内部語の「較正」はそのまま使ってよい。

**注意（同じ漢字の別語）**: 計量・測定の分野では、**「校正（calibration）」**が JIS Z 8103 で定義された正式用語（「測定標準によって提供される不確かさを伴う量の値と…指示値との関係を確立する操作」、[J6]）。**統計の「較正」と計測の「校正」は分野が違う別の語**。本件は統計の較正なので「較正」が正しい。JIS の文脈（測定器）と混ぜない。

---

## 3. 「検出核」— ある遺伝子が検出された核の数

### 3.1 英語

| 内部名 | 学術用語（英語） | 定義・出典 |
|---|---|---|
| 検出核（その遺伝子が検出された核の数） | **number of cells/nuclei in which the gene is detected** | 「遺伝子ごとに、その発現が測定された細胞の総数」＝ Scanpy の列 `n_cells`（[S14]）、Pegasus の `n_cells`「各遺伝子が測定された細胞の総数」（[S15]） |
| （同義・QC の閾値表現） | **expressed in at least N cells** / **detected in N cells** | Scanpy `filter_genes(min_cells=…)`「少なくとも min_cells 個の細胞で発現している遺伝子を残す」（[S16]） |
| （同義・査読論文） | 「遺伝子 i の発現が閾値 θ 以上である細胞の数」 | SINCERA [S17] |
| （割合。数ではなく比率） | **gene detection rate** / **fraction of cells in which a gene is detected** | 単一細胞の differential detection（DD）の定義「群間で遺伝子が検出される細胞の割合の差」（[S18]）、scBFA [S19] |

- ツール内の列名: **Scanpy は `n_cells`（per gene）**（「その発現が測定された細胞の数」[S14]）、**Pegasus は `n_cells`**（[S15]）。Seurat の `FindMarkers` は割合 `pct.1` / `pct.2`（「発現が検出された細胞の割合」）。
- **数（count）と割合（rate）は別**。内部の「検出核」は数なので、英語では **number of nuclei in which the gene is detected**。割合は別語（gene detection rate）。

### 3.2 日本語

- **数（count）に対応する単独の定着訳は確認できなかった。** 説明的に「**その遺伝子が検出された核の数**（number of nuclei in which the gene is detected）」と書くのが適切。
- 割合（rate）は `docs/用語集.md` §1 の「**遺伝子検出率（gene detection rate）**」が既にある。**「検出核」（数）と「遺伝子検出率」（割合）は別の量**なので、用語集でも区別を明記するのが望ましい。

---

## 4. 「セット平均」— 遺伝子セットのスコアの平均

### 4.1 英語

| 内部名 | 学術用語（英語） | 定義・出典 |
|---|---|---|
| セット平均（複数遺伝子の標準化スコアの平均） | **mean gene set score** / **set-averaged score** / **average score** / **mean score of the gene set** | GOAT（プレプリント）: 「遺伝子セットのスコアは、それぞれの遺伝子スコアの**平均**として定義される」（[S20]、**査読前**）。gssnng は `average_score` / `mean_z_score` 関数を実装（[S21]）。GSABenchmark は「各メソッドの**平均スコア**（mean score）」を集計（[S22]） |
| （同義・単一サンプル法） | **signature score** / **gene set score**（遺伝子セットの代表値） | singscore は遺伝子セットの**平均パーセンタイル順位**をスコアとする（[S6]） |
| （Seurat の実装） | **AddModuleScore** = セット内遺伝子の**平均発現量**（コントロール遺伝子群の平均を差し引く） | Seurat の公式チュートリアル／Tirosh 2016 由来（[S23]） |

- 「セット平均」に対応する **mean gene set score / set-averaged score** は文献にある。**ただし「平均」は何の平均か（標準化スコアか、順位か、発現量か）を必ず添える。**

### 4.2 日本語

- **「セット平均」という単独語は、学術用語として確認できなかった。**
- 日本語の説明では「**遺伝子セットの平均スコア**」「遺伝子セットに含まれる遺伝子の発現量の平均」のように書かれる（例: Seurat AddModuleScore の日本語解説 [S24]、東京大学 統計遺伝学講義資料 [S25]）。
- P-1（`docs/用語集.md` §4.5）では、本文で「**候補の 8 遺伝子の標準化スコアの平均**（mean score of the 8-gene set）」のように、**何の平均か**を添えるのが適切。

---

## 5. 英語のままで使うのが適切な語（訳すと誤解を招くもの）

| 語 | 日本語訳を付けない理由 |
|---|---|
| **gene-set permutation / gene_set permutation** | 分野固有の操作名。訳すと「遺伝子セットの並べ替え」となり、何を並べ替えるか（遺伝子ラベルか表現型ラベルか）が消える |
| **competitive null hypothesis** | 「競合帰無仮説」は訳語として確認できるが、Goeman & Bühlmann の定義（遺伝子サンプリングモデルに依存）を日本語で再現しにくい。初出時に英語を併記 |
| **bootstrap / resampling** | 「ブートストラップ」「リサンプリング」で通用。無理に訳さない |
| **detection rate**（割合） | 「検出率」でよいが、「検出核」（数）と混ざりやすい。英語を併記 |
| **calibration** | 「較正」で定着。ただし「校正（計測）」と別語であることを明記 |

---

## 6. 用語集への反映案（提案。適用はしていない。判断は PI）

`docs/用語集.md` への追記候補（**適用していない**）:

- §3 または §1 に、**「その遺伝子が検出された核の数」（数）と「遺伝子検出率」（割合）は別の量**であることを明記。英語は number of nuclei in which the gene is detected（列名 `n_cells`）。
- 本文では、**「サイズを揃えたランダム遺伝子セット（gene-set permutation）」**と書き、**復元ありで引いたこと**を併記する。**「sham set」とは書かない。**
- 「較正」は日本語で定着しているので維持。ただし**意味は「名目水準どおりの第一種の過誤率」**であり、計測の「校正（JIS）」とは別語であることを注記。
- P-1 では、本文で「**候補の 8 遺伝子の標準化スコアの平均**」のように何の平均かを添える。

---

## 7. 出典一覧（実在確認済み）

### 英語（査読論文・公式マニュアル）
- [S1] Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 2005;102(43):15545–15550. DOI: 10.1073/pnas.0506580102. PMID: 16199517（Crossref・PubMed で確認）
- [S2] Broad Institute / GSEA. *GSEA GenePattern Module Documentation v20*（公式マニュアル）. https://gsea-msigdb.github.io/gsea-gpmodule/v20/ — permutation type に **phenotype**（表現型ラベルを並べ替える）と **gene_set**（サイズを揃えたランダム遺伝子セットを作る）があると明記。gene_set permutation の説明文で確認。
- [S3] Goeman JJ, Bühlmann P. Analyzing gene expression data in terms of gene sets: methodological issues. *Bioinformatics* 2007;23(8):980–987. DOI: 10.1093/bioinformatics/btm051（Crossref で確認。**competitive** vs **self-contained**、**gene sampling** vs **subject sampling** の定義）
- [S4] Tian L, Greenberg SA, Kong SW, et al. Discovering statistically significant pathways in expression profiling studies. *PNAS* 2005;102(38):13544–13549. DOI: 10.1073/pnas.0506577102. PMID: 16174746（Crossref・PubMed で確認。帰無仮説 Q1/Q2、遺伝子セット構造をランダムにする検定）
- [S5] Efron B. Bootstrap Methods: Another Look at the Jackknife. *The Annals of Statistics* 1979;7(1):1–26. DOI: 10.1214/aos/1176344552（Crossref で確認。**bootstrap は復元ありで、順列分布ではない**と原文に明記）
- [S6] Foroutan M, Bhuva DD, Lyu R, et al. Single sample scoring of molecular phenotypes. *BMC Bioinformatics* 2018;19:406. DOI: 10.1186/s12859-018-2435-4（Crossref で確認。singscore、平均順位スコア）
- [S7] Bioconductor. *singscore 公式マニュアル / ビネット*（公式マニュアル）. https://bioconductor.posit.co/packages/3.24/bioc/vignettes/singscore/inst/doc/singscore.html — 「Hypothesis testing … is performed using a permutation test」、`generateNull()` はスコア対象と**同じ遺伝子数のランダム遺伝子セット**を作ると明記。
- [S8] Good P. *Permutation Tests: A Practical Guide to Resampling Methods for Testing Hypotheses*. Springer; 1994. DOI: 10.1007/978-1-4757-2346-5（Crossref で存在確認。順列検定の成書）
- [S9] Reallocating and Resampling: A Comparison for Inference. arXiv:1708.02102. https://arxiv.org/abs/1708.02102 — 「resampling＝復元あり（bootstrap）」「reallocating＝復元なし（permutation/randomization test）」の対比（**プレプリント**）
- [S10] Larson JL, Owen AB. Moment based gene set tests. *BMC Bioinformatics* 2015;16:132. DOI: 10.1186/s12859-015-0571-7（Crossref で確認。Monte Carlo 標本による順列分布の近似、competitive null の説明）
- [S12] Wasserstein RL. *P*-value calibration in multiple hypotheses testing. *Statistics in Medicine* 2017;36(26):4248–4266. DOI: 10.1002/sim.7330（Crossref で確認。p値較正＝帰無仮説下の Uniform(0,1) を基礎とする）
- [S13] Berger JO, Bayarri MJ, Sellke T. Calibration of ρ Values for Testing Precise Null Hypotheses. *The American Statistician* 2001;55(1):62–71. DOI: 10.1198/000313001300339950（Crossref で確認。p値の较正の原典の一つ）
- [S14] Scanpy. *scanpy.pp.calculate_qc_metrics*（公式マニュアル）. https://scanpy.scverse.org/en/stable/generated/scanpy.pp.calculate_qc_metrics.html — `n_cells_by_{expr_type}` は「Number of cells this expression is measured in」と定義。
- [S15] Pegasus. *pegasus.tools.preprocessing*（公式マニュアル）. https://pegasus.readthedocs.io/en/latest/_modules/pegasus/tools/preprocessing.html — `identify_robust_genes` で `data.var["n_cells"]`＝「Total number of cells in which each gene is measured」。
- [S16] Scanpy. *scanpy.pp.filter_genes*（公式マニュアル）. https://scanpy.scverse.org/en/stable/generated/scanpy.pp.filter_genes.html — `min_cells`＝「expressed in at least min_cells cells」、列 `n_cells` を遺伝子ごとに付与。
- [S17] Guo M, Wang H, Potter SS, Whitsett JA, Xu Y. SINCERA: A Pipeline for Single-Cell RNA-Seq Profiling Analysis. *PLoS Computational Biology* 2015;11(11):e1004575. DOI: 10.1371/journal.pcbi.1004575（Crossref で確認。遺伝子発現が閾値以上の細胞数の定義）
- [S18] Gilis J, Perin L, Malfait M, et al. Differential detection workflows for multi-sample single-cell RNA-seq data. *BMC Genomics* 2025;26:xxxx. DOI: 10.1186/s12864-025-12102-x（Crossref で確認。differential detection（DD）＝「群間で遺伝子が検出される細胞の割合の差」）
- [S19] Song F, et al. scBFA: modeling detection patterns to mitigate technical noise in large-scale single-cell genomics data. *Genome Biology* 2019;20:206. DOI: 10.1186/s13059-019-1806-0（Crossref で確認。gene detection rate（GDR）の定義）
- [S21] gssnng: Single-cell gene set scoring with nearest neighbor graph smoothed data. PMC10599965. https://pmc.ncbi.nlm.nih.gov/articles/PMC10599965/ — 遺伝子セットスコア関数として `average_score` / `mean_z_score` / `median_score` を挙げる。
- [S22] Stoica AF, et al. *GSABenchmark* リファレンスマニュアル（Bioconductor）. https://bioconductor.uib.no/packages/3.23/bioc/manuals/GSABenchmark/man/GSABenchmark.pdf — 「mean score for each method」の記述。
- [S23] Seurat. *AddModuleScore 公式チュートリアル*. https://satijalab.org/seurat/ — セット内遺伝子の平均発現量からコントロール遺伝子群の平均を引く。原著は Tirosh I, et al. *Science* 2016;352:189–196. DOI: 10.1126/science.aad0501。

### 日本語（査読誌・学会）
- [J1] 日本計量生物学会 ニュースレター 第151号. https://biometrics.gr.jp/newsletter/all/kaiho151.pdf — 「頻度流の第 I 種過誤確率に対して**較正**する必要がない」「第一種過誤確率の漸近的な一様制御」「分布的ロバスト最適化に基づく外部データ借用重みの理論的**較正**」
- [J2] 日本計算統計学会 講演要旨「ロジスティック回帰モデルにおける**ブートストラップ較正**と経験尤度」. https://www.jstage.jst.go.jp/article/jscssymo/15/0/15_27_1/_pdf/-char/ja
- [J3] 計量生物学 41(1)「イベント予測モデルの評価指標」. https://www.jstage.jst.go.jp/article/jjb/41/1/41_1/_pdf — 「**リスク較正**（risk calibration）」「較正プロット」「大域的較正」「較正スロープ」
- [J4] 医用画像情報学会雑誌 40(4)「機械学習における予測の不確実性定量化」. https://www.jstage.jst.go.jp/article/mii/40/4/40_75/_pdf/-char/ja — 「**確率較正**（probability calibration）」「再較正（re-calibration）」
- [J5] 人工知能学会 全国大会 2023「決定損失の期待値と分散を用いた較正」. https://www.jstage.jst.go.jp/article/pjsai/JSAI2023/0/JSAI2023_3E1GS204/_pdf/-char/en
- [J6] NITE「校正における測定不確かさの評価」. https://www.nite.go.jp/data/000001385.pdf — JIS Z 8103 の「**校正（calibration）**」定義（計測分野の別語）
- [S24] Seurat AddModuleScore の日本語解説（**ブログ。査読なし**）. https://zenn.dev/rchiji/books/fdd68b85675c8d/viewer/ad7bb6 — 補助的出典にとどめる
- [S25] 東京大学 統計遺伝学 講義資料（GenomeDataAnalysis7）. https://genome.m.u-tokyo.ac.jp/files/StatGenSummerSchool2025_GenomeDataAnalysis7.pdf — 「セット内に含まれる遺伝子の発現量の平均値（＝gene score）」（講義資料）

### 確認に用いた実在確認の方法
- DOI は Crossref REST API（`https://api.crossref.org/works/<DOI>`）で title・掲載誌・発行年を照合。
- PMID は NCBI E-utilities（esummary）で照合。
- 「sham set」の不存在は、Crossref の書誌検索（"sham gene set enrichment null distribution"）で用語としてヒットしないことで確認。
- 上記 [S8] の DOI は Crossref で存在確認のみ（本文は未取得）。

---

## 8. 限界・確認できなかったこと（正直に）

1. **「sham set」という学術用語は、査読論文・教科書・公式マニュアルで確認できなかった。** 医学の "sham"（シャム処置）は多数あるが、遺伝子セット解析の語ではない。**「無い」と断定できる範囲は、Crossref 書誌検索と web 検索で用語としてヒットしなかったことまで。**
2. **permutation と bootstrap の名称は、厳密な統計学と遺伝子セット解析の分野の慣行で食い違う。** 復元ありの操作を分野が "permutation" と呼ぶ例がある（[S2][S7]）。本報告では両論を併記し、断定していない。**どちらを使うかは PI の判断。**
3. **「検出核」「セット平均」に対応する単独の定着訳語は確認できなかった。** 説明的に書く案を示したが、訳語を確定する公的規範（JIS・学会用語集など）は見つからなかった。
4. [S20] GOAT の記述（遺伝子セットのスコア＝遺伝子スコアの平均）は **bioRxiv プレプリントであり査読前**。ここでは補助証拠として扱い、本文では [S21][S23] の査読・公式資料を主根拠にした。
5. 日本語の「較正」の出典 [J1]–[J5] は、いずれも**確率予測モデルの較正**（予測確率と実測頻度の一致）の文脈で「較正」を使っており、**本件（名目水準どおりの第一種の過誤率）と全く同じ用法の日本語一次資料は見つからなかった**。ただし [J1] に「第 I 種過誤確率に対して較正」という同一の用法があったため、同じ「較正」で通じると判断した。
6. [S9]（arXiv）、[S24]（ブログ）は**査読なし**。補助的出典にとどめた。
7. DOI のうち [S8] は Crossref で存在のみ確認し本文未取得。それ以外の英語一次資料 [S1][S3][S4][S5][S6][S10][S12][S13][S17][S18][S19] は Crossref または PubMed で書誌を照合済み。

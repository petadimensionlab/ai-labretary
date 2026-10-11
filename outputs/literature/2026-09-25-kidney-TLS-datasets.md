# 腎臓 TLS/TLT × 細胞老化 — 公開 scRNA/snRNA/空間データセット・カタログ

- 作成日: 2026-09-25（literature-scout）
- テーマ: 腎臓の三次リンパ組織（TLS/TLT）が腎臓病・腎癌で果たす役割を、Geneformer（human V2-316M）の in silico perturbation で解析し、細胞老化と結び付ける
- 検証方法: CELLxGENE Discover Curation API（`/curation/v1/collections` と `/collections/{id}`）の JSON を取得し、細胞数・h5ad URL・ファイルサイズ（バイト値）・疾患・組織・アッセイ・発生段階・`raw_data_location` を機械的に確認した。GEO は NCBI E-utilities（esummary）と GEO FTP の suppl ディレクトリ一覧で確認。論文の書誌情報（PMID）は Europe PMC REST API で DOI から照合。DDBJ GEA はディレクトリ一覧で確認。**データ本体のダウンロードはしていない**（h5ad の obs 列名など、ファイルを開かないと分からない項目は「未確認」と明記した）。
- サイズは API／FTP が返した値（1 GB = 10^9 byte で換算）。

---

## 0. サマリー表（優先度順）

| # | 優先度 | データセット | 種 | アッセイ | 細胞数 | サイズ | 状態ラベル（TLS/老化関連） | アクセス |
|---|---|---|---|---|---|---|---|---|
| A1 | **A** | KPMP Human Kidney Atlas v1.5 snRNA（Lake 2023 系列） | Human | 10x 3' v3 sn | 304,989 | 3.03 GB | AKI / CKD / normal、年齢 14〜70 代、B・形質細胞・T・DC | Open（CELLxGENE, CC BY 4.0） |
| A2 | **A** | KPMP Human Kidney Atlas v1.5 scRNA | Human | 10x 3' v3 sc | 225,177 | 1.82 GB | 同上 | Open |
| A3 | **A** | Li R 2022 Cancer Cell — RCC 腫瘍内/境界/正常/末梢血 scRNA | Human | 10x 5' | 270,855 | 1.62 GB | ccRCC（nonpapillary RCC）/良性/オンコサイトーマ/normal、IgA・IgG 形質細胞、class-switched memory B | Open（処理済み）/ raw は EGA（controlled） |
| A4 | **A** | Krishna 2021 Cancer Cell — ccRCC 多領域免疫 scRNA+TCR | Human | 10x 5' v2 | 167,283 | 1.18 GB | 腫瘍多領域/正常腎/リンパ節/末梢血、ICB 前後、B・Treg・naive T | Open |
| A5 | **A** | Zhang 2021 PNAS — ccRCC＋正常腎 scRNA（GSE159115） | Human | 10x 3' v2 | 20,509 + 6,044 | 0.23 + 0.06 GB | ccRCC vs normal、B・形質細胞 | Open |
| A6 | **A** | Reck 2025 Nat Commun — 尿管閉塞腎 snRNA（multiome の RNA 部分） | Human | 10x multiome | 46,957 | 0.26 GB | obstructive nephropathy vs normal、60〜84 歳、B・T | Open |
| B1 | B | KPMP Human Kidney Atlas v2.0 snRNA（Lake 2025 bioRxiv） | Human | 10x 3' v3 + multiome | 1,388,643 | 13.82 GB | AKI/CKD/normal、223 donors、10 代〜90 代 | Open（大容量・プレプリント） |
| B2 | B | KPMP Human Kidney Atlas v2.0 scRNA | Human | 10x 3' v3 | 348,984 | 4.39 GB | 同上、112 donors、細分化 T/B/DC | Open |
| B3 | B | KPMP Mouse Kidney v2.0 snRNA（Lake 2025） | Mouse | 10x 3' v2 + multiome | 315,764 | 1.82 GB | AKI/normal、2〜20 か月齢、B・形質細胞 | Open |
| B4 | B | Meylan 2022 Immunity — RCC Visium（GSE175540） | Human | Visium（空間） | 24 切片 | RAW.tar 1.12 GB | **スポット単位 TLS アノテーション CSV あり**、TLS+ / TLS− 腫瘍 | Open（GEO） |
| B5 | B | Yoshikawa 2023 JASN — 高齢マウス IRI 後 TLT 腎 snRNA（E-GEAD-595） | Mouse | 10x 3' v3/v3.1 sn | 7,485（sham）+ 15,968（IRI×3） | processed.zip 907 MB | 高齢（12 か月）IRI 30 日（TLT+）vs sham | Open（DDBJ GEA） |
| B6 | B | Tabula Muris Senis — Kidney（droplet） | Mouse | 10x 3' v2 | 21,647 | 0.34 GB | 1・3・18・20+ か月齢（加齢） | Open（CELLxGENE） |
| B7 | B | Chen 2025 Nat Genet — マウス腎 生涯 multiome | Mouse | 10x multiome | 203,139 | 1.56 GB | 胎生〜21 か月齢、性差、形質細胞・T | Open（CELLxGENE） |
| B8 | B | Wilson 2022 Nat Commun — 糖尿病性腎症 snRNA | Human | 10x 3' v3 / 5' v1 sn | 39,176 | 0.40 GB | T2DM（DKD）vs normal、50〜78 歳 | Open（CELLxGENE） |
| B9 | B | Malone 2020 — 腎移植拒絶 scRNA（GSE145927） | Human | scRNA（詳細未確認） | 未確認 | RAW.tar 247 MB | 拒絶腎、ドナー/レシピエント キメリズム | Open（GEO） |
| B10 | B | PICseq 腎移植拒絶（GSE189536） | Human | scRNA（相互作用細胞） | 未確認 | matrix 322 MB | 拒絶 vs 非拒絶（詳細未確認） | Open（GEO） |
| B11 | B | Kirita 2020 PNAS — マウス IRI 時系列 snRNA（GSE139107） | Mouse | snRNA | 未確認 | dge.txt 計約 126 MB | 4h〜6 週（AKI→修復不全） | Open（GEO） |
| C1 | C | Arazi 2019 Nat Immunol — ループス腎炎 AMP | Human | CEL-seq2 | — | — | LN vs control、age-associated B cells | **Controlled**（dbGaP phs001457.v1.p1 / ImmPort SDY997） |
| C2 | C | Bi 2021 Cancer Cell — 進行 RCC ICB 前後 | Human | scRNA | — | — | ICB 曝露前後 | raw **Controlled**（dbGaP phs002065.v1.p1）、処理済みは Broad SCP1288 |
| C3 | C | Braun 2021 Cancer Cell — ccRCC 病期別 | Human | scRNA+TCR | 164,722（論文記載） | — | 早期/局所進行/転移 | **アクセッション未確認** |
| C4 | C | Obradovic 2021 Cell — ccRCC TME | Human | scRNA | — | — | TREM2+ マクロファージ、再発 | **アクセッション未確認** |
| C5 | C | Stewart 2019 / Young 2018（成人正常腎・胎児腎） | Human | 10x 3' v2 | 40,268 など | 0.19〜0.23 GB | normal のみ（参照） | Open（CELLxGENE） |
| C6 | C | McEvoy 2022（生体ドナー腎）、Muto 2021、Tabula Sapiens Kidney | Human | 10x 3' v3 等 | 27,675 / 19,985 / 11,376 | 0.48 / 0.16 / 0.45 GB | normal のみ（参照・年齢分布） | Open（CELLxGENE） |
| C7 | C | 他の GEO ccRCC/IgAN/マウス加齢（§2 参照） | 混在 | — | — | 11 MB〜1.4 GB | 小規模・ラベル不十分 | Open |

**優先度 A の合計: 約 8.20 GB**（3.03 + 1.82 + 1.62 + 1.18 + 0.23 + 0.06 + 0.26 GB）。
注: KPMP v1.5 sn と sc は同一アトラスの別モダリティで、v1.0／v1.5／v2.0 は細胞が重複する（`is_primary_data=False` のものがある）。複数バージョンを混ぜて使ってはいけない。

---

## 1. データセット詳細

### A1/A2. KPMP Human Kidney Atlas v1.5（snRNA / scRNA）
- 論文: Lake BB et al., *Nature* 2023, "An atlas of healthy and injured cell states and niches in the human kidney." PMID 37468583, DOI 10.1038/s41586-023-05769-3（v1.0 のコレクションは bcb61471-2a44-4d00-a0af-ff085512674c。v1.5 は別コレクション 0f528c8a-a25c-4840-8fa3-d156fa11086f で、API 上は DOI 未登録）
- snRNA v1.5: dataset a12ccb9b-4fbe-457d-8590-ac78053259ef, 304,989 nuclei, 10x 3' v3, 3,034,341,738 byte
  - URL: https://datasets.cellxgene.cziscience.com/02d5db96-d67e-45f2-b556-b6bd99753665.h5ad
- scRNA v1.5: dataset dea717d4-7bc0-4e46-950f-fd7e1cc8df7d, 225,177 cells, 10x 3' v3, 1,821,593,028 byte
  - URL: https://datasets.cellxgene.cziscience.com/08648660-9898-4644-a9a1-ff3a159a0041.h5ad
- 組織: kidney / cortex / medulla / papilla。疾患: acute kidney injury, chronic kidney disease, normal。donor 数: sn 80、sc 77（API の donor_id）
- 細胞型: B cell, plasma cell, T cell, cytotoxic T, Treg, cDC, pDC など（CL 語彙）。年齢は「14-year-old」〜「seventh decade」（10 年刻みのものがある）
- 定義できる状態: **injured（AKI/CKD）vs healthy**。上皮の adaptive/maladaptive（failed-repair）状態を使えば、老化様の尿細管状態を start、健常状態を goal に設定できる。年齢（decade）で若年 vs 高齢の比較も可能。**TLS そのもののラベルはない**（B・形質細胞・T の存在比で代替指標を作ることはできる）
- アクセス: Open。CELLxGENE Discover 掲載データは CC BY 4.0（データセットごとの表示は未確認）。KPMP の利用規約（謝辞表記）にも従う
- Geneformer 準備度: `raw_data_location = raw.X`（生カウントあり、API で確認）。CELLxGENE スキーマ 7.1.0 では var のインデックスが Ensembl ID。`n_counts` は標準では含まれない → raw.X の行和から計算する必要がある
- 関連性: TLS — CKD では間質の B・形質細胞・T 浸潤があり、腎 TLT が形成される背景となる。老化 — failed-repair PT は老化・SASP 様の状態であり、年齢レンジも広い
- 優先度 A（最初に使う主データ）

### A3. Li R 2022 Cancer Cell — RCC の腫瘍内・境界・正常組織
- Li R, Ferdinand JR, Loudon KW et al., *Cancer Cell* 2022. PMID 36423636, DOI 10.1016/j.ccell.2022.11.001。CELLxGENE コレクション f7cecffa-00b4-4560-a29a-8ad626b8ee08（19 datasets）
- scRNA: dataset 5af90777-6760-4003-9dba-8f945fec6fdf, 270,855 cells, "10x 5' transcription profiling", 1,624,330,103 byte, 12 donors
  - URL: https://datasets.cellxgene.cziscience.com/8784de63-bec2-49ae-b147-8f16108cecce.h5ad
- 疾患ラベル: nonpapillary renal cell carcinoma（ccRCC を含む）, kidney benign neoplasm, kidney oncocytoma, normal。論文では腫瘍コア・腫瘍-正常境界・周囲正常・末梢血をサンプリング（領域が obs 列にあるかは未確認）
- 細胞型: IgA/IgG plasma cell, class switched memory B, naive B, Treg, 各種 T・DC（45 種）
- 付随する Visium（同コレクション）: 腫瘍-正常境界 eaf0c852-1dff-47ab-8de2-a33a59969a40（27,912 spots, 95.7 MB）、腫瘍コア 318cb2a6-743a-4404-9515-413d5f64b9ff（10,924 spots, 53.5 MB）ほか
- アクセス: 処理済み h5ad は Open。raw は EGA（EGAD00001008030, EGAD00001008781; **controlled**、コレクションのリンクで確認）
- Geneformer 準備度: raw.X あり、Ensembl ID。n_counts は計算する
- 関連性: TLS — 形質細胞の IgA/IgG 別アノテーションがあり、Meylan 2022 の TLS→形質細胞の流れを腫瘍/境界/正常の軸で検証できる。老化 — 高齢患者（50〜80 代）のみ。老化は二次的な軸になる
- 優先度 A

### A4. Krishna 2021 Cancer Cell — ccRCC 多領域の免疫景観
- Krishna C et al., *Cancer Cell* 2021. PMID 33861994, DOI 10.1016/j.ccell.2021.03.007。CELLxGENE コレクション 3f50314f-bdc9-40c6-8e4a-b0901ebfbe4c
- dataset bd65a70f-b274-4133-b9dd-0d1431b6af34, 167,283 cells, 10x 5' v2, 1,177,999,600 byte, 6 donors（ICB 未治療 2、治療後 4）
  - URL: https://datasets.cellxgene.cziscience.com/6602fa5a-d663-4fe0-8dc3-90191c3a013b.h5ad
- 論文記載の構成: 腫瘍多領域、リンパ節、正常腎、末梢血。paired TCR。raw は SRA PRJNA705464（公開、コレクションのリンクで確認）。HCA project 12f32054-8f18-4dae-8959-bfce7e3108e7
- 細胞型: B cell, naive T, Treg, effector CD4, cDC など 18 種
- 定義できる状態: 腫瘍 vs 正常 vs リンパ節（**二次リンパ組織を参照にできる**）、ICB 前後。TLS+/− のラベルはない
- Geneformer 準備度: raw.X あり、Ensembl ID
- 関連性: TLS — リンパ節の B/T 状態を「TLS 様状態」の goal として参照できる。老化 — 限定的（6 名、50〜69 歳）
- 優先度 A

### A5. Zhang 2021 PNAS — ccRCC・chRCC・正常腎
- Zhang Y et al., *PNAS* 2021, "Single-cell analyses of renal cell cancers reveal insights into tumor microenvironment, cell of origin, and therapy response." PMID 34099557, GEO GSE159115（14 samples、E-utilities で確認）
- CELLxGENE コレクション 1df8c90d-d299-4b2e-a54d-a5a80f36e780
  - ccRCC: be39785b-67cb-4177-be19-a40ee3747e45, 20,509 cells, 232,252,366 byte — https://datasets.cellxgene.cziscience.com/f77a9049-7cbe-477d-8d3f-b417bba03147.h5ad
  - normal: f801b7a9-80a6-4d09-9161-71474deb58ae, 6,044 cells, 58,231,471 byte — https://datasets.cellxgene.cziscience.com/5f0b9675-aed8-4870-8ad8-7a3f38037486.h5ad
  - chRCC: 4c6f9f26（2,576 cells, 0.04 GB）も同コレクションにある
- 細胞型: B cell, plasma cell, T cell など。10x 3' v2、65〜76 歳
- 関連性: 小さく軽いので、パイプラインの smoke テストや腫瘍 vs 正常の再現確認に向く
- 優先度 A（小規模の検証用）

### A6. Reck 2025 Nat Commun — 尿管閉塞によるヒト腎障害 snRNA
- Reck M, Baird DP et al., *Nat Commun* 2025, "Multiomic analysis of human kidney disease identifies a tractable inflammatory and pro-fibrotic tubular cell phenotype." PMID 40399382, DOI 10.1038/s41467-025-59997-4。CELLxGENE コレクション 307ce143-7cc8-4813-99b1-3797834149c9
- dataset 867757c1-3b1a-49d9-a0cd-17767eb160cc, 46,957 nuclei, 10x multiome（RNA 部分）, 258,738,667 byte, 12 donors, 60〜84 歳
  - URL: https://datasets.cellxgene.cziscience.com/0fe5eee4-380d-4bd9-8735-ada5e03021d9.h5ad
- 疾患: obstructive nephropathy vs normal。B cell, T cell を含む
- 関連性: TLS — 炎症・線維化を起こす尿細管表現型は、マウスの Vcam1+ PT（TLT 周囲に局在、Yoshikawa 2023）に相当しうる（相同性は未検証）。老化 — 全例高齢
- Geneformer 準備度: raw.X あり
- 優先度 A

### B1/B2/B3. KPMP Atlas v2.0（Lake BB et al., bioRxiv 2025）
- "Cellular and Spatial Drivers of Unresolved Injury and Functional Decline in the Human Kidney", DOI 10.1101/2025.09.26.678707（**査読前プレプリント**。Europe PMC で PMID なし）。コレクション 9c9d04c4-8899-417f-bb6f-6107dcadf14f。GitHub KPMP-Scientific/KPMP-Atlas-v2、KPMP DOI 10.48698/16dd-vj20
- B1 snRNA human: 7ff0197b-d175-49bf-b4fa-150fe0995d93, 1,388,643 nuclei, 13,818,560,828 byte, 223 donors, 年齢は 14 歳〜「tenth decade」 — https://datasets.cellxgene.cziscience.com/06b8d9c2-0bde-48bd-88d9-730993fc95e2.h5ad
- B2 scRNA human: 91f31e05-56d8-46fc-b408-d90c9228a81b, 348,984 cells, 4,392,345,855 byte, 112 donors — https://datasets.cellxgene.cziscience.com/255ce81a-eb4e-48ff-9f18-e50455ccb24d.h5ad
- B3 snRNA mouse: 9ac3d74d-1d94-4ff9-abf9-fe5abe4bd2a6, 315,764 nuclei, 1,821,316,660 byte, 70 mice, 2〜20 か月齢, AKI/normal — https://datasets.cellxgene.cziscience.com/19aab2ef-654a-4047-8f60-3529c041851d.h5ad
- 特長: 細胞型が細かい（human sn で 63 種。migratory DC、naive Treg、TEMRA、cycling T など）。**"unresolved injury"（修復不全）状態**が定義されている。年齢レンジが最も広く、ヒトの加齢解析に最適
- 注意: 1.39M 細胞のトークン化・推論は §3 の「1 万細胞超」に当たる → サブサンプリング計画と PI 承認（G2）が必要
- 優先度 B（本命の拡張データ。サイズとプレプリントであることから B とした）

### B4. Meylan 2022 Immunity — RCC Visium（TLS アノテーション付き）
- Meylan M, Petitprez F, Becht E et al., *Immunity* 2022. PMID 35231421, DOI 10.1016/j.immuni.2022.02.001。GEO **GSE175540**（Homo sapiens、24 samples、E-utilities で PMID 35231421 に紐付くことを確認）
- ファイル: GSE175540_RAW.tar 1,122,877,440 byte。サンプルごとに `*_filtered_feature_bc_matrix.h5`（約 17〜40 MB）、`*_TLS_annotation.csv.gz`、画像、`tissue_positions_list`。別途 bulk RNA-seq TPM（10 MB）
  - https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175540/suppl/
  - コード: https://github.com/maximemeylan/Meylan_et_al_2022
- 論文記載: FFPE 12 + 凍結 12 の ccRCC 切片（GEO の件数 24 と一致）
- 定義できる状態: **スポット単位の TLS / 非 TLS**、TLS+ / TLS− 腫瘍。Visium スポットは多細胞の混合なので Geneformer の単一細胞仮定に合わない → TLS シグネチャーの抽出と検証に使い、perturbation の主対象にはしない
- Geneformer 準備度: 10x h5（生カウント、gene_ids 列に Ensembl ID があるのが通例。ファイル未開封のため未確認）
- 優先度 B（TLS の「goal 状態」遺伝子セットを定義する根拠データとして重要）

### B5. Yoshikawa 2023 JASN — 高齢マウス IRI 後の TLT 腎 snRNA
- Yoshikawa T et al., *J Am Soc Nephrol* 2023, "Tertiary Lymphoid Tissues Are Microenvironments with Intensive Interactions between Immune Cells and Proinflammatory Parenchymal Cells in Aged Kidneys." PMID 37548710, DOI 10.1681/ASN.0000000000000202（柳田研）
- アクセッション: raw FASTQ DRA015783、Cell Ranger の raw/filtered 行列は **DDBJ GEA E-GEAD-595**（ディレクトリで確認: `E-GEAD-595.processed.zip` 907 MB, 2023-08-04）
  - https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-595/E-GEAD-595.processed.zip
- デザイン（本文）: 12 か月齢の雄 C57BL/6J。片側 IRI 30 日後 3 匹（TLT+）、sham 1 匹。10x 3' v3/v3.1。sham 7,485 nuclei、IRI 統合 15,968 nuclei
- 所見: TLT 内の炎症性線維芽細胞と TLT 外の線維化線維芽細胞を区別。NF-κB/IFN 活性化を伴う Vcam1+ 障害 PT が TLT 周囲に局在
- 状態定義: aged-injured（TLT+）vs aged-sham。**若齢対照はない**。sham が n=1 のため統計は弱い
- Geneformer: マウスなので、ヒトオルソログへの変換（1:1 対応のみ）またはマウス Geneformer が必要。ゲノムは mm10 と推定（本文の他解析の記載から。行列の gene ID は未確認）
- 優先度 B（生物学的な関連は最も直接的。種の問題で B とした）

### B6. Tabula Muris Senis — Kidney（加齢）
- Tabula Muris Consortium, *Nature* 2020. PMID 32669714, DOI 10.1038/s41586-020-2496-1
- droplet: b8c618e5-4b3d-4566-8a3f-7e40047f5c54, 21,647 cells, 10x 3' v2, 0.34 GB — https://datasets.cellxgene.cziscience.com/cfd31d35-ceb2-4990-b0b5-3434ecbb7967.h5ad
- FACS/Smart-seq2: 524179b0-b406-4723-9c46-293ffa77ca81, 1,833 cells, 0.03 GB — https://datasets.cellxgene.cziscience.com/7769f88b-96e2-46e8-a2b2-2b7282804da7.h5ad
- 年齢: 1 か月、3 か月、18 か月、20 か月以上。障害モデルなし（自然加齢による TLT の有無は不明）
- 優先度 B（マウスの加齢軸の参照）

### B7. Chen 2025 Nat Genet — マウス腎の生涯 multiome
- Chen S et al., *Nat Genet* 2025. PMID 40259083, DOI 10.1038/s41588-025-02161-x。dataset 25818bf7-e2a7-41ec-8ff2-bc369c0ff4f5, 203,139 nuclei, 1.56 GB, 32 mice, 胎生〜21 か月, plasma cell/T/cDC — https://datasets.cellxgene.cziscience.com/29eec938-c857-4dc0-a3c1-db864bdea6dd.h5ad
- 優先度 B（加齢・性差の軸）

### B8. Wilson 2022 Nat Commun — DKD snRNA
- Wilson PC, Muto Y et al., *Nat Commun* 2022. PMID 36068241, DOI 10.1038/s41467-022-32972-z。dataset e067e5ca-e53e-485f-aa8e-efd5435229c8, 39,176 nuclei（primary 20,215）, 403,168,674 byte — https://datasets.cellxgene.cziscience.com/473d2e1c-a579-4682-96aa-b5b166d7b23e.h5ad
- T2DM vs normal。API 上は B/T の細胞型ラベルなし（13 種）→ 免疫細胞は少ないと推定
- 優先度 B（CKD 上皮の老化軸の補助）

### B9–B11. GEO（腎移植・マウス AKI）
- **GSE145927**: Malone AF et al. 2020（PMID 32669324）、腎移植拒絶のドナー/レシピエント免疫キメリズム、17 samples、`GSE145927_RAW.tar` 247 MB ＋ metadata.csv.gz 2.1 MB。形式は CSV（カウント行列か正規化値かは未確認）。移植腎の TLO との関連を探るための候補
- **GSE189536**: 腎移植拒絶の Physically Interacting Cells（PICseq）、PMID 37638864、32 samples、`matrix.mtx.gz` 322 MB（MTX 形式、生カウントの可能性が高いが未確認）
- **GSE139107**: Kirita Y et al., *PNAS* 2020（PMID 32571916）、マウス IRI 時系列 snRNA（control, 4h, 12h, 2d, 14d, 6wk）、dge.txt.gz 6 ファイル計約 126 MB＋metadata。若齢マウスのみ。failed-repair PT（Vcam1+）の元の定義
- 優先度 B

### C1. Arazi 2019 Nat Immunol — ループス腎炎（AMP）
- PMID 31209404, DOI 10.1038/s41590-019-0398-x。ImmPort SDY997、raw は **dbGaP phs001457.v1.p1（controlled）**（検索結果で確認。dbGaP ページは直接未確認）。HCA project 2d559a6e-7cd9-432f-9f6e-0e4df03b0888
- age-associated B cell（ABC）シグネチャーの原典として重要。データは PI の承認と DAR 手続きなしには使用しない
- 補足: Kondo M, Sato Y et al., *Kidney Int Rep* 2026（PMID 42707204）はループス腎炎の TLS が加齢・慢性障害と関連すると報告（組織学的研究で、scRNA データなし・未確認）

### C2–C4. 進行 RCC（controlled または未確認）
- **Bi K et al., Cancer Cell 2021**（PMID 33711272, DOI 10.1016/j.ccell.2021.02.015）: raw は dbGaP phs002065.v1.p1（controlled）。処理済みは Broad Single Cell Portal SCP1288（本文で確認。SCP はダウンロードに Google ログインが必要な場合がある。ライセンスは未確認）
- **Braun DA et al., Cancer Cell 2021**（PMID 33711273, DOI 10.1016/j.ccell.2021.02.013）: 164,722 cells（抄録）。Europe PMC で全文を取得できず、**アクセッション未確認**
- **Obradovic A et al., Cell 2021**（PMID 34019793）: 全文を取得できず、**アクセッション未確認**（GitHub Aleksobrad/single-cell-rcc-pipeline の存在のみ検索結果に出た）

### C5–C7. 参照・小規模データ
- Stewart BJ 2019 Science（PMID 31604275）成人正常腎 40,268 cells（9ea768a2, 0.19 GB）、Young MD 2018 Science（PMID 30093597）mature 40,268 cells（105c7dad, 0.23 GB）— 同一データの別掲載の可能性あり
- McEvoy CM 2022 Nat Commun（PMID 36496458）生体ドナー腎 27,675 cells（0bae7ebf, 0.48 GB, 33〜65 歳, 19 donors）— 健常の年齢軸として有用
- Muto 2021 Nat Commun（DOI 10.1038/s41467-021-22368-w）snRNA 19,985（9df60c57）、Tabula Sapiens Kidney 11,376 cells（2423ce2c, 0.45 GB, 61 歳 1 名）、HCA Michigan 正常腎 48,783（0b4a15a7, 0.66 GB）、mBDRC 97,125（60a29d0b, 1.41 GB, bioRxiv）
- GEO（すべて E-utilities で存在を確認）:
  - GSE121636（ccRCC 免疫 + PBMC, PMID 33504936, RAW 131 MB）
  - GSE156632（ccRCC, PMID 35605062, 86 MB）
  - GSE152938（RCC 多病型, PMID 34722263, 367 MB）
  - GSE171306（両側性 RCC, PMID 34168986, 76 MB）
  - GSE207493（ccRCC scRNA 1.4 GB ＋ scATAC 36 GB, PMID 36607615）
  - GSE233236（汎癌 B 細胞統合。腎に限定されない, PMID 39047727, count.mtx 194 MB）
  - GSE171314（IgA 腎症 5 samples, PMID 33936064, TXT 37 MB）
  - GSE131685（正常腎 3 名, PMID 31896769, 85 MB）
  - マウス加齢: GSE140010（若齢 3〜5 か月 vs 高齢 18 か月 IRI 間質, PMID 33951342, 11 MB。GEO 上 1 sample）、GSE262959（C57BL/6 加齢腎アトラス, BMC Nephrol 2024 PMID 38641839。suppl は `result.txt.gz` 5.1 MB のみで、生カウントかどうか不明）
  - Sato Y 2022 JCI（CD153/CD30、PMID 34813503）の CD45+ scRNA は DDBJ DRA011230（FASTQ のみ。2,344 cells。再処理が必要）

---

## 2. ポジティブ/ネガティブコントロール用 遺伝子セット（文献）

| セット | 遺伝子 | 用途 | 出典 |
|---|---|---|---|
| TLS 12-chemokine | CCL2, CCL3, CCL4, CCL5, CCL8, CCL18, CCL19, CCL21, CXCL9, CXCL10, CXCL11, CXCL13 | TLS 状態スコア／goal 定義 | Messina JL et al., Sci Rep 2012, PMID 23097687 |
| TLS コア（リンパ組織形成） | CXCL13, CCL19, CCL21, LTB, LTA, CXCR5, CCR7, MS4A1, CD79A/B, CR2 (CD21), LAMP3 | TLS 誘導のポジコン（過剰発現/活性化） | Sautès-Fridman C et al., Nat Rev Cancer 2019, PMID 31092904；腎 TLT: Sato Y et al., JCI Insight 2016, PMID 27699223 |
| 腎 TLT 線維芽細胞 | CXCL13, CCL19, NGFR (p75NTR), PDPN, ALDH1A2（レチノイン酸産生） | 高齢腎 TLT のストローマ状態 | Sato 2016（PMID 27699223）；Arai H et al., Immunol Rev 2021, PMID 33951198（遺伝子は総説・原著で要確認） |
| 老化関連 T 細胞（SA-T） | TNFSF8 (CD153), TNFRSF8 (CD30), PDCD1 | TLT 拡大と老化をつなぐ軸（CD153/CD30 阻害で TLT が縮小） | Sato Y et al., J Clin Invest 2022, PMID 34813503 |
| 加齢関連 B 細胞（ABC） | ITGAX (CD11c), TBX21 (T-bet), FCRL5, ZEB2 | LN・加齢 TLT の B 細胞状態 | Arazi A et al., Nat Immunol 2019, PMID 31209404（ITGAX/TBX21 以外は一般知見。原著の遺伝子リストは要確認） |
| TLS 形質細胞（RCC） | IGHG1–4, IGHA1, JCHAIN, MZB1, XBP1, PRDM1 | TLS→形質細胞の分化軸 | Meylan M et al., Immunity 2022, PMID 35231421（原著の "TLS signature" の正確な遺伝子リストは Table で要確認） |
| TLS 9 遺伝子（メラノーマ） | CD79B, CD1D, CCR6, LAT, SKAP1, CETP, EIF1AY, RBP5, PTGDS | 腎外の外部シグネチャーとして | Cabrita R et al., Nature 2020, PMID 31942071（遺伝子は記憶に基づく。原著で要照合） |
| 細胞老化コア | CDKN1A (p21), CDKN2A (p16INK4a), TP53, LMNB1（低下）, GLB1 | 老化状態のポジコン | Gorgoulis V et al., Cell 2019, PMID 31675495 |
| SASP | IL6, CXCL8, IL1A, IL1B, CCL2, SERPINE1, MMP1, MMP3, IGFBP3, IGFBP7, GDF15 | SASP スコア | Coppé JP et al., PLoS Biol 2008, PMID 19053174 |
| SenMayo（125 遺伝子） | 原著 Supplementary の全リスト | 老化スコア（組織横断で検証済み） | Saul D et al., Nat Commun 2022, PMID 35974106 |
| 腎の障害・修復不全 PT | HAVCR1, LCN2, VCAM1, CCL2, CDKN1A, SPP1 | 修復不全＝老化様上皮の start 状態 | Kirita Y et al., PNAS 2020, PMID 32571916；Lake 2023, PMID 37468583；Yoshikawa 2023, PMID 37548710 |
| ネガティブコントロール | ハウスキーピング（ACTB, GAPDH, RPLP0, PPIA, TBP）＋ 発現量をマッチさせたランダム遺伝子（≥3 セット） | ヌル分布 | CLAUDE.md §6 の規定 |

---

## 3. 未確認事項・ギャップ

1. **ヒトの単一細胞データに TLS+/TLS− のラベルが付いたものは見つからなかった。** ヒト腎で TLS を直接アノテーションしたデータは、Meylan 2022 の Visium（空間、RCC）だけ確認できた。単一細胞で TLS 状態を定義するには、(a) Visium の TLS スポットから遺伝子シグネチャーを作って sc データに移す、(b) 腫瘍/サンプル単位で CXCL13/CCL19/B 細胞比率から TLS-high/low を推定する、のどちらかが必要（いずれも推定であり、正解ラベルではない）。
2. **腎疾患（非腫瘍）のヒト TLT scRNA** — IgAN・LN・移植腎で TLT をアノテーションした公開 scRNA は確認できなかった。AJP 2024 の論文（"Single-Cell Analysis ... TLS ... Kidney Injury and CKD"）は 403 エラーで本文を確認できず（公開データの再解析と推定）。
3. CELLxGENE の h5ad の **obs 列名**（例: KPMP の condition/region/年齢区分、Li R の腫瘍領域）は開かないと分からず、未確認。`n_counts` は含まれないため計算が必要。Ensembl ID の版（CELLxGENE スキーマ 7 系の GENCODE）と Geneformer V2 のトークン辞書の対応率は未確認。
4. CELLxGENE の CC BY 4.0 はポータルのポリシーに基づく記載で、データセット個別のライセンス表示は API から取得していない。
5. Braun 2021・Obradovic 2021 のアクセッションは未確認。Bi 2021 SCP1288 のダウンロード条件は未確認。LN（Arazi）と RCC（Bi）の raw は controlled。
6. GSE145927／GSE189536／GSE262959 は行列の形式（生カウントか）、細胞数、ラベルが未確認。
7. マウスの TLT データ（E-GEAD-595, DRA011230）は小規模（n=1〜3）で、若齢対照がない。KPMP マウス v2.0（2〜20 か月、AKI）に「高齢×AKI」の組合せが含まれるかは未確認（年齢と疾患のクロス集計は h5ad を見ないと分からない）。
8. KPMP v2.0 は査読前（bioRxiv 2025）。
9. 空間データは今回、Li R 2022 と Meylan 2022 の Visium と Marshall 2022 Slide-seqV2（CELLxGENE 8e880741、正常腎と DKD マウス）だけ確認した。KPMP の空間データ（Visium/CODEX）は未調査。

## 4. 推奨する次のアクション
1. PI 承認後、A1〜A6（約 8.2 GB）を `data/raw/` にダウンロード（別途 PI の許可が必要）。まず A5（0.29 GB）で Geneformer のトークン化 smoke テストを行う
2. B4（Meylan Visium, 1.1 GB）から TLS スポットのシグネチャーを抽出し、goal 状態の定義に使う
3. KPMP v2.0 はサブサンプリング計画（例: 免疫＋PT で 1 万細胞/条件以下）を立て、G2 承認を得る
4. マウス TLT（B5）はオルソログ変換の方針を決めてから使う

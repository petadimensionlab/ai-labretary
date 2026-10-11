# A6 管状系（正常 PT → PT Injured/Inflammatory）の「退化・老化・傷害・修復」陽性対照 — 候補材料

- 作成: 2026-10-11 / `literature-scout`
- 用途: EXP-20261011-003 §5.1・§12-3 の「管状系の退化軸の陽性対照を A6 用に組み直す」ための**候補材料**。
  **最終的な固定・採用・閾値は planner / PI が行う。本調査は提案までで、判定基準も計画も作っていない。**
- 対象データ: `data/raw/kidney_tls/A6_UUO_Reck2025.h5ad`（読み取り専用・sha256 `8d9cca13…`、46,957 核）
- 本調査はリストと検出率の提示のみ。埋め込み・摂動・得点・順位表は計算していない。
- **盲検中の 5 遺伝子（TP53, CDKN1A, LTBR, SPP1, LCN2）は候補に含めない**（本ファイルにも値も出さない）。

## 0. 結論（何が得られたか）

- Reck ら 2025 の**補足データ 5（Supplementary Data 5）**に、A6 自身の
  **「PT Injured+PT Inflammatory（上昇）対 正常 PT（下降）」**と
  **「PT Inflammatory（上昇）対 PT Injured（下降）」**の遺伝子ごとの差次発現（L2FC・%Expressed）が
  一次情報として公開されている（第 I 節・第 II 節）。これが A6 の PT Injured/Inflammatory の**直接の一次出典**。
- そこから **40 個**を候補材料として並べた（`candidate_markers.csv`）。すべて A6 の var と Geneformer 辞書（Ensembl 一致）にある。
- A6 の生カウント（raw.X ≥ 1）で測った検出率は、**Reck 論文の `%Expressed - Celltype` とほぼ一致**した
  （例: CDH6 72.9% 対 72.9%、TPM1 50.5% 対 50.6%、DCDC2 62.2% 対 61.9%）。測り方は §3 のとおり。
- 活性化側＝間質系の陽性対照 **71 と重なる遺伝子が 28 個**ある（`overlap_with_activation71.csv`）。
  **これらは 2 本の軸で共通に上がるため、管状系の陽性対照としては「軸に固有でない」候補**になる（§5）。

## 1. 出典（一次資料。PMID / DOI は PubMed と出版社ページで確認）

| 記号 | 文献 | PMID / DOI | 本調査で使った箇所 |
|---|---|---|---|
| **Reck2025** | Reck M, Baird DP ら, Nat Commun 16:4745 (2025) 「Multiomic analysis of human kidney disease identifies a tractable inflammatory and pro-fibrotic tubular cell phenotype」 | **PMID 40399382** / **DOI 10.1038/s41467-025-59997-4** | **補足データ 5 第 I・II 節**（PT の差次発現）、本文 Fig 1d・Fig 2c・Fig 7。A6 データの出典 |
| **Lake2023** | Lake BB, Menon R, Winfree S ら, Nature 619:585–594 (2023) 「An atlas of healthy and injured cell states and niches in the human kidney」（KPMP アトラス） | **PMID 37468583** / **DOI 10.1038/s41586-023-05769-3** | aPT マーカー（Extended Data Fig 7c: VCAM1, DCDC2, HAVCR1。Extended Data Fig 7f・補足表 11: CDH6。補足表 13: PROM1, DCDC2）。退化マーカー（本文・補足図 2: 傷害マーカー） |
| **Gisch2024** | Gisch DL, Brennan M, Lake BB ら, Nat Commun 15:433 (2024) 「The chromatin landscape of healthy and injured cell types in the human kidney」 | **PMID 38199997** / **DOI 10.1038/s41467-023-44467-6** | aPT マーカー（Fig 4c: ITGB3, PROM1, TPM1）、aPT の転写因子（Fig 6a: KLF6, ELF3, SOX4, KLF10, NFKB1） |
| （参考）Saul2022 | Saul D ら（SenMayo）, Nat Commun 13:4827 (2022) | PMID 35974106 / DOI 10.1038/s41467-022-32552-1 | 老化マーカー。ただし主要な SenMayo 遺伝子は A6 の PT で検出されない（§6-3）。既存取得 `outputs/literature/2026-10-11-senescence-and-degenerated-lists/` |

- Reck2025 補足データ 5 の取得元（一次ファイル）:
  `https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-59997-4/MediaObjects/41467_2025_59997_MOESM3_ESM.xlsx`
  （sha256 `4217d777972ea3b93804b824f2ce1ac037c26021630a585c563cb17d092e3912`）。
  シート `Supplementary Data 5` の第 I 節（列 B–G）と第 II 節（列 I–N）を直接読んだ。
  第 I 節の見出しは「PT Injured+PT Inflammatory (up) vs PT Healthy (PT S1-S3) (down)」、
  第 II 節は「PT Inflammatory (up) vs PT Injured (down)」。フィルターは「3% 以上の細胞で発現、Wilcoxon、BH 補正」。

### 1.1 指定との食い違い（確認できたこと）

- 依頼の「**Lake ら Nat Commun 2023（KPMP）の近位尿細管（S1/S2/S3・Injured/Inflammatory）マーカー**」は、
  一次資料では次のように分かれる。
  - **KPMP アトラスは Lake ら Nature 2023**（Nat Commun ではない）。Lake2023 が持つのは **aPT（adaptive PT）**で、
    見出しの「PT Injured / PT Inflammatory」という**ラベル自体は A6（Reck2025）のもの**である（`celltype_l2` で確認）。
  - 「Lake ら Nat Commun 2023」に当たる一次資料は **Gisch ら Nat Commun 2024**（KPMP のクロマチンアトラス）と読める。
- したがって本調査は、**A6 の PT Injured/Inflammatory には Reck2025 補足データ 5** を、
  **aPT の裏づけには Lake2023 と Gisch2024** を一次出典として使った。この読み替えが PI の意図と違う場合は、材料を入れ替える。

## 2. A6 の状態（`precheck_a6.json` と一致。本調査で数え直した）

| 状態 | ラベル | 核数 |
|---|---|---|
| 目標（退化・老化） | PT Injured ∪ PT Inflammatory | **2,632**（Injured 2,019 ＋ Inflammatory 613） |
| 出発（正常 PT） | PT S1 ∪ PT S2 ∪ PT S3 | **3,977** |

## 3. 検出率の測り方（結果を見る前の固定ではなく、材料の併記）

- `data/raw/kidney_tls/A6_UUO_Reck2025.h5ad` の **raw.X（整数カウント）** を読み、遺伝子 g の
  検出率 = 「生カウント ≥ 1 の核の割合」。`precheck_a6.py` と同じ定義。
- 読み取りのみ。書き込みなし。得点・順位・埋め込みは計算していない。
- **照合**: 本調査の A6 検出率は Reck2025 補足データ 5 の `%Expressed - Celltype`（= PT Injured+Inflammatory での
  発現細胞割合）と数ポイント以内で一致した（CDH6・TPM1・DCDC2・ITGB8・PROM1 などで確認）。
  これは「A6 h5ad の生カウント」と「Reck 論文の集計」が同じデータを指していることの確認である。

## 4. 候補リスト（40 個）

`candidate_markers.csv` の要約。**検出率は A6 の生カウント ≥1。**「Reck I / II」は補足データ 5 の第 I 節・第 II 節の
L2FC（空欄はその節に載っていない＝その比較では差次発現の基準を満たしていない）。**向きはすべて上昇側**。

| 遺伝子 | 群 | 一次出典（図表） | Reck I L2FC | Reck II L2FC | A6 PT 目標 検出率 (n=2,632) | A6 正常 PT 検出率 (n=3,977) | 根拠（1 行） |
|---|---|---|---|---|---|---|---|
| HAVCR1 | T1 古典的 | Reck SD5-I/II・Fig 1d; Lake2023 ED Fig 7c | 0.50 | 0.42 | **32.0%** | 4.6% | KIM-1。A6 の PT Injured で最も標準的な傷害マーカー |
| VCAM1 | T1 古典的 | Reck SD5-I・Fig 2; Lake2023 ED Fig 7c | 0.77 | 0.27 | **31.3%** | 2.4% | aPT/failed-repair の代表。正常 PT でほぼ消えている |
| DCDC2 | T1 古典的 | Reck SD5-I/II・Fig 1d; Lake2023 ED Fig 7c, ST13 | 1.24 | 0.82 | **62.2%** | 14.0% | aPT の代表マーカー。検出率が高く目標で強い |
| CDH6 | T1 古典的 | Reck SD5-I; Lake2023 ED Fig 7f, ST11 | 1.51 | 0.32 | **72.9%** | 43.7% | 委縮・拡張尿細管の aPT マーカー。検出率が最も高い部類 |
| ITGB8 | T1 古典的 | Reck SD5-I/II・Fig 1d | 1.61 | 0.88 | **76.9%** | 25.9% | Reck の核心傷害遺伝子。L2FC が大きく検出率も高い |
| ITGB6 | T1 古典的 | Reck SD5-I/II・Fig 1d; Lake2023（aTAL） | 0.40 | 0.44 | 20.2% | 1.0% | 傷害・修復上皮マーカー。正常 PT でほぼ 0 |
| ITGB3 | T1 古典的 | Reck SD5-I; Gisch2024 Fig 4c | 0.40 | 0.24 | 22.3% | 1.0% | aPT マーカー。正常 PT でほぼ 0 |
| PROM1 | T1 古典的 | Reck SD5-I/II・Fig 1d; Lake2023 ST13; Gisch2024 Fig 4c | 0.73 | 0.49 | **39.7%** | 1.7% | CD133。aPT/aEpi の代表。正常 PT でほぼ 0 |
| TPM1 | T1 古典的 | Reck SD5-I/II; Gisch2024 Fig 4c | 0.92 | 0.96 | **50.5%** | 9.4% | aPT マーカー。目標で 5 割、正常で 1 割 |
| IGFBP7 | T1 古典的 | Reck SD5-I; Lake2023 本文・補足図 2 | 0.65 | 0.30 | **49.2%** | 20.5% | Lake が挙げる既知の傷害・退化マーカー |
| KLF6 | T2 転写因子 | Reck SD5-I/II; Gisch2024 Fig 6a | 0.65 | 0.45 | **43.0%** | 14.7% | aPT の転写因子ネットワークの中心 |
| ELF3 | T2 転写因子 | Reck SD5-I/II; Gisch2024 Fig 6a | 0.39 | 0.21 | 20.9% | 2.9% | aPT の転写因子。正常 PT で低い |
| SOX4 | T2 転写因子 | Reck SD5-I/II; Gisch2024 Fig 6a | 0.29 | 0.26 | 16.9% | 2.3% | 傷害・修復を促す転写因子として報告 |
| NFKB1 | T2 転写因子 | Reck SD5-I/II・Fig 7; Gisch2024 Fig 6a | 0.31 | 0.33 | 28.8% | 11.8% | Reck が炎症 PT の駆動因子とする NF-κB |
| NFKBIZ | T2 転写因子 | Reck SD5-I/II | 0.58 | 0.74 | 28.8% | 2.1% | NF-κB 標的。**PT Inflammatory で 57%** と目標特異的 |
| RELB | T2 転写因子 | Reck SD5-I | 0.11 | 0.11 | 8.1% | 2.0% | NF-κB 経路。L2FC は小さいが正常 PT と分かれる |
| JUN | T2 転写因子 | Reck SD5-I・Fig 7 | 0.16 | — | 12.7% | 2.8% | AP-1。Reck が炎症 PT の駆動因子とする |
| ICAM1 | T3 炎症 | Reck SD5-II・Fig 2 | — | 0.19 | 3.7% | 0.6% | Reck の炎症 PT の代表（VCAM1+ICAM1+）。**検出率は低い** |
| CCL2 | T3 炎症 | Reck SD5-II・Fig 2 | — | 0.15 | 3.7% | 0.1% | 炎症 PT のケモカイン。**検出率は低い** |
| CXCL1 | T3 炎症 | Reck SD5-I/II・Fig 2 | 0.11 | 0.34 | 2.3% | 0.0% | 炎症 PT のケモカイン。**検出率は低い** |
| TNC | T3 炎症 | Reck SD5-I/II | 0.20 | 0.54 | 6.9% | 0.6% | テネイシン C。炎症 PT 側で上がる |
| ANXA2 | T3 炎症 | Reck SD5-I/II | 0.27 | 0.35 | 16.2% | 0.9% | アネキシン A2。正常 PT でほぼ 0 |
| MMP7 | T3 炎症 | Reck SD5-I/II | 0.22 | 0.47 | 9.6% | 0.6% | マトリックスメタロプロテアーゼ。正常 PT でほぼ 0 |
| C3 | T3 炎症 | Reck SD5-I/II | 0.21 | 0.26 | 9.0% | 0.4% | 補体 C3。傷害尿細管の報告がある |
| SOX9 | T3 炎症 | Reck SD5-II | — | 0.12 | 6.4% | 1.1% | 修復・発生の転写因子。Reck の第 II 節で上昇 |
| SOD2 | T4 ストレス・老化 | Reck SD5-I/II | 0.53 | 0.72 | **32.8%** | 14.3% | 酸化ストレス。炎症 PT でさらに上がる |
| HDAC9 | T4 ストレス・老化 | Reck SD5-I/II | 0.24 | 0.48 | 12.3% | 1.1% | **Reck が老化関連・G2/M 停止として挙げる**。炎症 PT で 31% |
| GDF15 | T4 ストレス・老化 | Reck SD5-I | 0.15 | — | 8.2% | 1.2% | ストレス・老化関連。正常 PT で低い |
| BCL2 | T4 ストレス・老化 | Reck SD5-I | 0.67 | 0.24 | **69.8%** | 45.9% | 抗アポトーシス（Reck の老化プログラム）。検出率は高い |
| HSPB1 | T4 ストレス・老化 | Reck SD5-I/II; Lake2023 補足表 12 | 0.17 | 0.12 | 12.0% | 4.7% | 熱ショック・ストレス。Lake の保存的退化遺伝子にもある |
| HSPA1A | T4 ストレス・老化 | Reck SD5-I | 0.21 | — | 11.9% | 2.5% | 熱ショックタンパク質。正常 PT で低い |
| HSPA1B | T4 ストレス・老化 | Reck SD5-I | 0.17 | — | 11.0% | 2.3% | 熱ショックタンパク質。HSPA1A と同挙動 |
| BIRC3 | T4 ストレス・老化 | Reck SD5-I/II | 0.20 | 0.33 | 8.3% | 0.3% | NF-κB 標的・抗アポトーシス。正常 PT でほぼ 0 |
| BICC1 | T5 高検出 | Reck SD5-I | 1.04 | 0.30 | **92.0%** | 81.0% | L2FC 大だが正常 PT でも高く、状態の分離力は弱い |
| RHEX | T5 高検出 | Reck SD5-I | 1.48 | 0.39 | **80.5%** | 38.5% | L2FC 大・検出率高。機能の同定は薄い |
| GLIS3 | T5 高検出 | Reck SD5-I | 0.81 | 0.56 | **79.7%** | 57.0% | 転写因子。検出率が高い |
| UGT2B7 | T5 高検出 | Reck SD5-I | 0.86 | — | **79.5%** | 48.2% | 代謝酵素。検出率が高い |
| PPFIBP1 | T5 高検出 | Reck SD5-I | 1.05 | 0.21 | **71.7%** | 33.4% | L2FC 大・検出率高。機能の同定は薄い |
| FNDC3B | T5 高検出 | Reck SD5-I | 0.71 | 0.58 | **71.5%** | 41.8% | 検出率高。機能の同定は薄い |
| SEMA5A | T5 高検出 | Reck SD5-I | 0.80 | −0.13 | **68.7%** | 43.7% | 検出率高。第 II 節では逆符号（炎症＞傷害でない） |

- 「Reck I」= 正常 PT に対する PT Injured+Inflammatory の上昇（L2FC）。「Reck II」= PT Injured に対する PT Inflammatory の上昇。
- 群 T1・T2 は**目標に固有でよく知られた**マーカー。T3 は炎症 PT 寄り。T4 はストレス・老化寄り。
  T5 は検出率が高いだけで**状態への固有性が弱い**（正常 PT でも高く出るものが多い）。

## 5. 活性化側＝間質系の陽性対照 71 と重なるもの（28 個）

「71」の各遺伝子が、Reck 補足データ 5 の**第 I 節または第 II 節で上昇**していたもの。`overlap_with_activation71.csv`。

| 遺伝子 | Reck 節 | Reck L2FC / %exp | A6 PT 目標 検出率 | A6 正常 PT 検出率 |
|---|---|---|---|---|
| RAD51B | I | 0.35 / 0.419 | 43.8% | 25.6% |
| SDC2 | I | 0.14 / 0.144 | 14.0% | 6.2% |
| FRY | I | 0.26 / 0.179 | 17.6% | 2.3% |
| SIPA1L1 | I | 0.15 / 0.660 | 66.9% | 63.3% |
| MICAL3 | I+II | 0.15 / 0.464; 0.46 / 0.615 | 46.4% | 41.3% |
| ARMC9 | I | 0.12 / 0.204 | 20.1% | 12.5% |
| ANK2 | I | 0.34 / 0.430 | 43.1% | 41.2% |
| MAML3 | I | 0.28 / 0.318 | 32.0% | 17.6% |
| PAWR | I | 0.43 / 0.466 | 46.7% | 25.4% |
| DCLK1 | I+II | 0.30 / 0.128; 0.13 / 0.170 | 11.7% | 1.3% |
| ESYT2 | I | 0.15 / 0.271 | 27.1% | 18.3% |
| FOXP1 | I | 0.41 / 0.556 | 58.9% | 39.7% |
| PRICKLE1 | I+II | 0.70 / 0.372; 0.86 / 0.630 | 36.4% | 9.3% |
| ANTXR1 | I | 0.19 / 0.126 | 11.4% | 0.4% |
| JAZF1 | I+II | 0.21 / 0.156; 0.14 / 0.228 | 15.5% | 3.3% |
| SESTD1 | I+II | 0.19 / 0.333; 0.21 / 0.439 | 33.3% | 23.5% |
| COL8A1 | I+II | 0.16 / 0.079; 0.15 / 0.119 | 7.5% | 1.3% |
| PARD3B | I | 0.59 / 0.780 | 79.5% | 61.1% |
| NFIB | I | 0.10 / 0.839 | 84.0% | 82.5% |
| PDE5A | I | 0.17 / 0.151 | 14.9% | 4.8% |
| ZEB1 | I+II | 0.16 / 0.219; 0.31 / 0.336 | 21.9% | 14.6% |
| PTPRD | I | 0.14 / 0.839 | 83.9% | 80.6% |
| PLCB1 | I+II | 0.37 / 0.241; 0.19 / 0.292 | 25.8% | 11.1% |
| ATP8B1 | II | 0.10 / 0.131 | 8.5% | 2.5% |
| SLC20A2 | II | 0.21 / 0.282 | 20.7% | 24.3% |
| AFDN | II | 0.15 / 0.597 | 52.6% | 51.9% |
| NAV3 | II | 0.17 / 0.098 | 2.9% | 1.4% |
| PDGFC | II | 0.31 / 0.581 | 51.9% | 53.1% |

- **本調査の 40 候補のうち、71 と名前が重なるものは 0 個**（`candidate_markers.csv` の `in_activation_71` は全て no）。
  一方で、上の 28 個は「71 の側の遺伝子が A6 の管状系の目標でも上がる」もので、**2 本の軸に共通**である。
- 使い方の含意（判断はしない）: 管状系の陽性対照に**軸固有性**を求めるなら、上の 28 個は除く候補になる
  （とくに NFIB・PTPRD・SIPA1L1・AFDN・PDGFC・PARD3B は正常 PT でも高く、分離力が弱い）。
  逆に「共通に上がる遺伝子」も材料として残すなら、この一覧をそのまま使う。

## 6. 確認できたこと / 確認できず（一次資料で確かめられなかったもの）

1. **確認できたこと**
   - Reck2025 補足データ 5 の第 I 節・第 II 節から、PT Injured/Inflammatory の差次発現（L2FC・%Expressed）を
     一次ファイルのセルから直接読んだ（`41467_2025_59997_MOESM3_ESM.xlsx`、シート `Supplementary Data 5`）。
   - A6 の `celltype_l2` に「PT Injured」2,019 核・「PT Inflammatory」613 核が実在することを確認した。
   - 40 候補・28 重なり遺伝子が、すべて A6 の var と Geneformer の辞書（Ensembl ID 一致）にあることを確認した。
2. **確認できず / 材料として弱いと分かったもの**
   - **TGM2・CXCL2・CXCL8・IL6・TIMP1** は Reck 補足データ 5 の PT の節に見当たらない（基準を満たしていない）か、
     A6 での検出がほぼ 0（例: TIMP1 は PT 目標で 1.2% だが Myofibroblast 側 18.3%）。管状系の陽性対照には使えない。
   - **老化の主要マーカー CDKN2A・SERPINE1・IGFBP3・IGFBP2** は A6 の PT 目標で検出率が 0〜0.8% で使えない。
     SenMayo（Saul2022）の遺伝子も多くが A6 の PT で検出されない（既存取得の README 参照）。
   - **CST3・CLU**（Lake2023 の退化マーカー）は A6 の PT 目標で 6.5%・4.0% と低く、陽性対照としては弱い。
   - Reck2025 の **CosMx（6,000 遺伝子パネル）の差次発現（補足データ 9〜11）**は本調査では展開していない
     （補足データ 5 の snRNA 側を優先）。追加の炎症遺伝子候補が要る場合は別途取得できる。
   - Reck2025 の**図の番号**（Fig 1d・Fig 2c・Fig 7）は本文の記述に基づく。補足データ 5 の節見出しは一次ファイルで確認済み。
3. **盲検について**: TP53・CDKN1A・LTBR・SPP1・LCN2 は候補に含めず、本ファイルに値も書いていない。

## 7. 限界・しないこと

- 本調査は**文献と検出率の一覧**である。陽性対照の採否・個数・検出率の床・閾値は決めていない（planner / PI の仕事）。
- A6 は snRNA（10x multiome）。**検出率は単一核での値**であり、sscRNA での値ではない。
- 検出率が高いことは「陽性対照として働く」ことの**必要条件であって十分条件ではない**。軸に沿うかは埋め込み後にしか測れない
  （EXP-20261011-003 §4.4 の P-4 は相関・記述）。
- 「上がる」は Reck2025 補足データ 5 の集計（L2FC）と A6 の検出率に基づく。**因果や機構は主張しない。**

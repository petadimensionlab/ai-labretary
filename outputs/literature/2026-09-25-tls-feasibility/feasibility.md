# 腎 TLS × 細胞老化 — Geneformer ISP の実現可能性評価（実データ）

- 作成: 2026-09-26（hypothesis-builder / Labretary）。対象: `data/raw/kidney_tls/` の A1〜A6（**読み取り専用**。raw.X の整数カウントを h5py でチャンク読みし、元ファイルは変更していない）
- スクリプト: `scripts/analysis/tls_feasibility.py`（サンプル単位・状態単位の集計）、`scripts/analysis/tls_feasibility_extra.py`（推定サブセット・仮説用の細胞数）
- 出力 CSV（このディレクトリ）: `donor_<DS>.csv`, `states_<DS>.csv`, `assoc_summary.csv`, `special_cells.csv`, `donor_assoc_extra.csv`, `start_state_genes.csv`, `cxcl13_source.csv`, `tumour_cd8_cxcl13.csv`, `hyp_state_counts.csv`, `token_dict.csv`
- **注意: ここで出てくる「TLS-high / TLS-low」はすべて推定ラベルです。** 組織学的な TLS の有無はどのデータにも付いていません（カタログ §3-1）

---

## 0. 結論（先に）

1. **腎疾患（A1/A2/A6）では、サンプル単位の TLS スコア（推定）を定義でき、疾患・障害 PT 比率と関連した。** TLS スコアは B＋形質細胞比率と強く相関した（Spearman ρ = 0.75 / 0.80 / 0.81、n = 75 / 77 / 12 donor）。CKD・AKI は normal より高かった（A1: Kruskal p = 3.6×10⁻⁴、A2: p = 9.9×10⁻⁵、A6: 閉塞腎と正常の Mann-Whitney p = 0.0025）。
2. **年齢との関連は弱く、疾患で交絡している。** A1 では疾患で調整した OLS の年齢係数 = +0.012 /年（p = 0.019）、normal だけに絞ると ρ = 0.66（n = 23、p = 6.8×10⁻⁴）。一方 A2（調整後 p = 0.22）と A6（p = 0.89）では再現しなかった。したがって「加齢で TLS が増える」は **示唆にとどまる**（A1 の 1 データセットのみ）。
3. **障害 PT（aPT / PT Injured）の比率は TLS スコアと相関した**（A1 ρ = 0.54、n = 75、p < 10⁻⁵ / A2 ρ = 0.43、n = 74、p = 1×10⁻⁴ / A6 ρ = 0.84、n = 12、p = 6×10⁻⁴）。ただし疾患の交絡を除いていない。
4. **ISP の state として使える細胞数が確保できるのは、aPT（A2: CKD 内で TLS-high 7,533 細胞 / 14 donor、TLS-low 4,002 / 8 donor）、CCL19/CCL21/CXCL13 陽性の線維芽細胞（A1: 777 細胞 / 65 donor、推定ラベル）、ccRCC の B 細胞→形質細胞（A3: 腫瘍＋境界部で活性化 B 1,498、IgG 形質細胞 283）である。**
5. **老化関連 T 細胞（SA-T 様、TNFSF8⁺PDCD1⁺ T）と加齢関連 B 細胞（ABC 様、ITGAX⁺ かつ TBX21⁺ または FCRL5⁺ の B）は非腫瘍データにほとんどいない**（SA-T 様: A1 14、A2 68、A6 1 細胞 / ABC 様: A1 6、A2 23、A6 19 細胞）。→ **これらは ISP の state にはできない。** 摂動する遺伝子プログラムとしてのみ使う。
6. **ccRCC では、CXCL13 のほとんどが CD8 T 細胞から出ている**（A3: CXCL13 カウントの 72.5% が CD8⁺T_RM-CXCL13 / A4: 74% が exhausted＋proliferating CD8 / A5: 75% が CD8 T）。→ 腫瘍の pseudo-bulk TLS スコアは、TLS より **CD8 T の疲弊状態を反映している可能性が高い**。腫瘍内の B＋形質細胞は少なく（リンパ球の 0.3〜4.5%）、CXCL13⁺ CD8 比率との相関は有意でなかった（A3: ρ = 0.31、n = 11 腫瘍 / 11 donor、p = 0.36）。
7. **老化マーカーの検出率が低い。** CDKN2A（p16）は PT 系で 1% 未満（A1 aPT 0.3%、A2 aPT 0.6%）で、delete 摂動の対象にはできない。CDKN1A（p21）は A2 aPT で 7.8%（PT-S1/S2 は 3.8%）。SenMayo は検証済みの遺伝子リストが手元にないため**使っていない**（CDKN1A/CDKN2A と SASP 10 遺伝子で代用）。
8. **Geneformer トークン辞書**: 解析した 76 遺伝子のうち、辞書にないのは IGHG1 と IGHA1 だけだった（Ig 定常領域は辞書に入っていない）。CXCL13, CCL19, CCL21, LTB, LTA, CXCR5, CCR7, MS4A1, CDKN1A, CDKN2A, VCAM1, HAVCR1, TNFSF8, PDCD1, ITGAX, TBX21, PRDM1, XBP1, LTBR, RELB, NFKB2, ハウスキーピング 5 遺伝子はすべて収載されていた。A4 は var が 15,934 遺伝子しかなく、CCL3/4/5/18、RPLP0、IGHG1 が欠けている。

---

## 1. 方法（要点）

- **B / 形質細胞**: CELLxGENE の `cell_type`（A3 は著者の `broad_type`）で判定した。リンパ球 = T・NK・B・形質細胞・ILC・MAIT
- **pseudo-bulk CPM**: サンプル内の対象遺伝子カウントの合計 ÷ 全カウントの合計 × 10⁶（raw.X から計算）
- **TLS スコア（推定）**: CXCL13, CCL19, CCL21, LTB, LTA, CXCR5, CCR7, MS4A1 の log1p(CPM) をそれぞれ z 化し、その平均をとった。z 化の基準は腎組織サンプルだけ（A3 の血液・脂肪・副腎・転移・血栓、A4 の PBMC・リンパ節は基準から外した）。**細胞組成に強く依存する指標**である（CCR7/LTB は naive T にも多く、CXCL13 は腫瘍では CD8 T 由来）
- **TLS 群**: 基準サンプルのうち 300 細胞以上のものを三分位に分けた（high / mid / low）。サンプル単位は、A1/A2/A6 が donor、A3 が donor × 部位、A4 が部位（author_sample、6 donor）、A5 が donor × 腫瘍/正常
- **老化スコア**: 細胞ごとの log1p(CP10k) で、sen = mean(CDKN1A, CDKN2A)、SASP = 10 遺伝子（IL6, CXCL8, IL1A, IL1B, SERPINE1, MMP1, MMP3, IGFBP3, IGFBP7, GDF15; Coppé 2008 PMID 19053174）の平均
- 統計: Spearman、Mann-Whitney、Kruskal、腫瘍 vs 正常の donor 内対応 Wilcoxon、OLS（TLS ~ 年齢 + 疾患）。多重比較の補正はしていない（探索的）

## 2. サンプル単位の結果（`assoc_summary.csv`, `donor_*.csv`）

| DS | 基準サンプル数（donor） | B+PL / 全細胞 中央値 [IQR] | TLS vs B+PL 比率 ρ | 疾患/部位との関連 | TLS vs 年齢 |
|---|---|---|---|---|---|
| A1 KPMP sn | 75 (75) | 0.19% [0.08–1.06] | 0.75 (p = 9×10⁻¹⁵) | normal −0.45 < CKD −0.12 < AKI 0.03（中央値, Kruskal p = 3.6×10⁻⁴） | ρ = 0.33 (p = 0.004)。疾患で調整した β = +0.012/年 (p = 0.019)。normal のみ ρ = 0.66 (n = 23) |
| A2 KPMP sc | 77 (77) | 1.04% [0.17–2.52] | 0.80 (p = 4×10⁻¹⁸) | normal −0.66 < CKD 0.18 < AKI 0.45（p = 9.9×10⁻⁵） | ρ = 0.25 (p = 0.027)。調整後 p = 0.22。normal のみ ρ = 0.32 (n = 26, p = 0.11) |
| A6 UUO sn | 12 (12) | 1.10% [0.38–4.33] | 0.81 (p = 0.0014) | 閉塞腎 0.46 vs 正常 −0.43（p = 0.0025） | 有意な関連なし（60〜84 歳のみ） |
| A3 Li RCC | 25 (11) | 1.54% | 0.59 (p = 0.0017) | 腫瘍と正常の donor 内比較で 8/9 donor が腫瘍 > 正常（Wilcoxon p = 0.012）。独立群の比較では有意差なし | ρ = −0.58（donor 内の複数サンプルによる擬似反復があり、解釈しない） |
| A4 Krishna | 24 (6) | 0.86% | 0.43 (p = 0.037) | 腫瘍と正常に差なし（p = 0.77） | n = 6 donor で評価不能 |
| A5 Zhang | 12 (8) | 0.37% | 0.69 (p = 0.012) | 腫瘍と正常に差なし（p = 0.27）。対応あり 3/4 donor | n が小さい |

**TLS 群と疾患のクロス表（交絡の大きさ）**
- A2 TLS-high 26 donor: AKI 8、CKD 14、normal 4 / TLS-low 26 donor: AKI 1、CKD 8、normal 17
- A1 TLS-high: AKI 8、CKD 14、normal 3 / TLS-low: AKI 2、CKD 8、normal 15
- A6: TLS-high 4 donor はすべて閉塞腎、TLS-low 4 donor はすべて正常 → **A6 では TLS 群と疾患を分離できない**
- → 仮説では **CKD donor の中だけで TLS-high と TLS-low を比べる**（A2 CKD: high 14 donor、low 8 donor）

## 3. 細胞状態ごとの細胞数と老化スコア（`states_*.csv`, `special_cells.csv`, `hyp_state_counts.csv`）

### 3.1 障害 PT（start 候補）と健常 PT（goal 候補）
| DS | 状態 | 細胞数 | TLS-high / TLS-low | CDKN1A⁺ | CDKN2A⁺ | VCAM1⁺ | sen（平均 ± SD） | SASP（平均 ± SD） | 中央値遺伝子数 |
|---|---|---|---|---|---|---|---|---|---|
| A2 | aPT | 29,397 | 12,274 / 6,185 | 7.8% | 0.6% | 60.7% | 0.047 ± 0.176 | 0.283 | 1,660 |
| A2 | dPT | 6,721 | 3,132 / 1,274 | 4.6% | 0.1% | 15.7% | 0.032 ± 0.162 | 0.092 | 926 |
| A2 | PT-S1/S2 | 13,019 | 2,377 / 5,941 | 3.8% | 0.0% | 9.4% | 0.025 ± 0.134 | 0.085 | 988 |
| A1 | aPT | 28,396 | 13,827 / 6,204 | 4.7% | 0.3% | 35.7% | 0.025 ± 0.123 | 0.184 | 2,524 |
| A1 | PT-S1 | 13,126 | 4,732 / 5,546 | 1.6% | 0.1% | 8.7% | 0.011 ± 0.095 | 0.105 | 1,390 |
| A6 | PT Injured | 2,019 | 1,180 / 128 | 2.2% | 0% | 29.6% | 0.020 ± 0.136 | 0.122 | 1,312 |
| A6 | PT Inflammatory | 613 | 460 / 22 | 4.2% | 0% | 37.0% | 0.040 ± 0.192 | 0.179 | 1,341 |
| A6 | PT S1 | 1,815 | 98 / 297 | 1.4% | 0% | 2.7% | 0.012 ± 0.102 | 0.070 | 1,263 |

- **A2 CKD donor に限った場合**: aPT は TLS-high 7,533（14 donor）、TLS-low 4,002（8 donor）。健常 PT（PT-S1/S2＋S3、normal donor）は 6,149
- A2 aPT の中で TLS-high と TLS-low の donor（22 対 22）の検出率を比べると、TLS-high で高かったのは SPP1（79% vs 45%, p < 0.001）、LCN2（9.6% vs 1.5%, p < 0.001）、HAVCR1（14% vs 7%, p = 0.007）、TP53（16% vs 12%, p = 0.022）、CXCL8（p = 0.004）だった。CDKN1A は差がなく（9.9% vs 11.7%, p = 0.43）、NFKB2 は TLS-high でむしろ低かった（12% vs 17%, p = 0.02）。IL7R・LTB・CCL19 の陽性率も TLS-high で高いが、周囲の免疫細胞からの ambient RNA の可能性がある（`start_state_genes.csv`）

### 3.2 線維芽細胞
- A1 の線維芽細胞系（FIB, aFIB, dFIB, M-FIB, dM-FIB, MYOF, cycMYOF）は 21,769 細胞。**CCL19, CCL21, CXCL13 のいずれかが 1 カウント以上**ある細胞（推定の「TLS オーガナイザー様」）は 777 細胞（65 donor）で、TLS-high 562 / mid 141 / low 69 だった。内訳は CCL19⁺ 252、CCL21⁺ 411
- donor ごとのオーガナイザー様比率は TLS スコアと相関した（ρ = 0.59、n = 71、p < 10⁻⁶）。ただし TLS スコア自体に CCL19/CCL21/CXCL13 が含まれるため、**一部は循環的な相関**である
- A2 では同じ条件の細胞が 123（30 donor）、A6 では 20、A3 では 11（A3 は CCL19 の検出が全細胞で 0）、A4 では 65。→ 線維芽細胞の仮説は **A1 でしか組めない**
- 老化: A1 FIB の CDKN1A⁺ は 1.0%、aFIB は 2.9%。snRNA なので検出率が低い。ccRCC の線維芽細胞は高く、A4 の myofibroblast で CDKN1A⁺ 65%、SASP 0.63、A3 の Stress Fibro で 72%

### 3.3 B 細胞・形質細胞・T 細胞
| DS | 状態 | 細胞数 | TLS-high / TLS-low | 備考 |
|---|---|---|---|---|
| A2 | B | 3,261 | 2,392 / 105 | TLS-low 側が少ない（TLS スコアの定義上、当然そうなる） |
| A2 | PL | 1,051 | 678 / 53 | FCRL5⁺ 65% |
| A6 | Memory B | 3,356 | 3,076 / 14 | 3,315 が閉塞腎 |
| A3 | 活性化 B（AREG＋RHOB） | 2,017（腫瘍＋境界 1,498） | 780 / 35（腫瘍＋境界） | CDKN1A⁺ 27〜31%（naive B は 2.5%） |
| A3 | Plasma-IgG | 786（腫瘍＋境界 283） | — | 中央値遺伝子数 417。IGHG1 は辞書にない |
| A3 | Plasma-IgA | 563（腫瘍＋境界 253） | — | |
| A3 | CD8⁺T_RM-CXCL13 | 2,449（腫瘍 2,059） | — | CXCL13⁺ 99%、CDKN2A⁺ 21%、TNFSF8⁺ 18% |
| A4 | CD8A⁺ Exhausted | 17,893 | — | CXCL13⁺ 51% |

- SA-T 様・ABC 様の細胞数は §0-5 のとおり。ccRCC には SA-T 様が A3 で 4,646、A4 で 2,056 細胞あるが、腫瘍内の PDCD1 は疲弊マーカーでもあるため、**老化 T 細胞と同一視できない**

## 4. 計算コストの見積もり方法
- 下限: 1 forward の秒数 = 中央値遺伝子数 ÷ 7,200 tok/s（MPS bf16、batch 8、PI 周知の実測値）。1 細胞 × 1 遺伝子あたり元の入力と摂動後の 2 forward とした
- 上限: EXP-20260925-001 の fp32 smoke（macrophage 50 細胞、1 遺伝子あたり約 40 秒）から、1 細胞 × 1 遺伝子 ≈ 0.8 秒
- delete / inhibit は**その遺伝子を発現している細胞にしか効かない**（Geneformer は発現していない細胞を飛ばす）。発現率の低い遺伝子（CDKN1A 7.8% など）は、発現細胞を濃縮して抽出する必要がある
- どの仮説でも、まず smoke（10 細胞 × 3 遺伝子）で bf16 の実測値を取ってから本実行する

## 5. できないこと・限界
- TLS 状態は推定。組織学（Meylan 2022 Visium, B4）による検証はまだ。TLS スコアは組成指標であり、TLS 構造（B 濾胞・HEV・FDC）の有無は判定できない
- A1/A2 は `is_primary_data = False`（KPMP の別版と細胞が重複する）。KPMP v2.0 と混ぜてはいけない
- A3/A4 は 10x 5'、A5 は 3' v2、A1/A6 は snRNA。Geneformer の事前学習の分布とずれている
- donor 数: A6 12、A3 11、A4 6、A5 8。ccRCC で TLS-high と TLS-low を比べる統計的な検出力は低い
- 多重比較は補正していない。これらの数値は**仮説を作るための探索的な記述**である

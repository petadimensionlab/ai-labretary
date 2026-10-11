# 三次リンパ組織オーガナイザー転換の陽性対照候補（文献調査）

- 作成: literature-scout（Labretary 委託）/ 2026-09-27
- 目的: EXP-20260927-004 の事前登録のための陽性対照候補の一次情報確認
- 方法: PubMed / PMC / Cell / Nature 誌の一次情報を検索し、PMID・DOI を確認できたもののみ記載。仮説の提案・判定基準の設計は行っていない。
- 用語は `docs/用語集.md` に合わせた（例: 「三次リンパ組織」「腎線維芽細胞」「三次リンパ組織オーガナイザー様の状態」）。

---

## Q1（最重要）: 線維芽細胞／間葉系細胞を三次リンパ組織オーガナイザー様（LTo 様）へ転換させることが実験的に示された遺伝子

### 結論を先に書く
- **「線維芽細胞に単一遺伝子を過剰発現させて（内因性に）LTo 様へ転換させた」という、遺伝子導入型の直接的な実験は見つからなかった**。
- 見つかったのは主に **(a) LTBR への外部からのリガンド刺激（アゴニスト抗体または LTA/LTB の異所性発現）で間葉系前駆細胞が LTo 様（CCL19/CCL21/CXCL13 産生）に転換した**という報告と、**(b) LTBR を線維芽細胞系譜で欠失させると LTo 様への成熟が止まる（喪失）**という報告である。
- 「NIK を線維芽細胞に過剰発現させたら CCL19/CCL21 が誘導された」という趣旨のご指摘の候補は、検索した範囲では**確認できなかった**（未確認）。非古典的 NF-κB（RELB, NFKB2/p52, NIK）が下流で必要という報告はあるが、NIK 自体の異所性過剰発現実験は見当たらない。

### (a) 受容体シグナルを外から与えて LTo 様へ転換した報告

| # | 遺伝子・経路 | 与え方 | 対象細胞 | 誘導されたもの | 出典 |
|---|---|---|---|---|---|
| 1 | LTBR（アゴニスト抗体で刺激） | **受容体シグナルを外から与えた**（培養した脂肪組織および FACS で分離した前駆脂肪細胞に、アゴニスト性抗 LTβR 抗体を投与） | マウス胎児性の前駆脂肪細胞（adipocyte precursor） | Cxcl13, Ccl21, Ccl19, Icam1, Vcam1 の誘導。RelB / NF-κB2（非古典的 NF-κB）経路を介して Pparg・Cebpa（脂肪分化因子）を抑制し、脂肪分化から リンパ節間質細胞様分化へ運命を転換 | Bénézech C, et al. *Immunity*. 2012;37(4):721-734. PMID: 22940098 |
| 2 | LTA（リンフォトキシンα）の異所性発現 | **受容体シグナルを外から与えた**（ラットインスリンプロモーター下で膵島β細胞に LTA を異所性発現させるトランスジェニックマウス, RIP-LTα） | 膵臓・腎臓の非リンパ組織（周囲の間質・血管内皮を含む） | リンパ節様の慢性炎症巣（T/B 細胞域の分離、一次・二次濾胞、高内皮細静脈=HEV）を新生（lymphoid neogenesis）。この論文自体は CCL19/CCL21/CXCL13 の直接測定ではなく、組織学的なリンパ組織様構造の記載が中心 | Kratz A, Campos-Neto A, Hanson MS, Ruddle NH. *J Exp Med*. 1996;183(4):1461-1472. PMID: 8666904 |
| 3 | LTA+LTB（LTαβ）の異所性発現 | **受容体シグナルを外から与えた**（RIP-LTαβ トランスジェニックマウス、LTβ欠損との比較） | 膵臓の非リンパ組織 | 末梢リンパ節アドレシン（PNAd）と HEV 特異的スルホトランスフェラーゼ（HEC-6ST）の発現を伴うリンパ組織新生を誘導。CCL19/CCL21/CXCL13 の直接測定は本論文の主眼ではない | Drayton DL, Ying X, Lee J, Lesslauer W, Ruddle NH. *J Exp Med*. 2003;197(9):1153-1163. PMID: 12732657 |

**解釈上の注意**: #2, #3 は「線維芽細胞という細胞自律的な遺伝子過剰発現」ではなく、「上皮細胞（膵島β細胞）にリガンドを異所性発現させ、周囲の間質にパラクラインで作用させた」モデルである。CCL19/CCL21/CXCL13 の直接的な誘導は #1（Bénézech 2012）が最も明確。

### (b) 遺伝子欠失（ノックアウト）で LTo 様への成熟が失われた報告

| # | 遺伝子 | 欠失の方法 | 対象細胞 | 失われたもの | 出典 |
|---|---|---|---|---|---|
| 4 | Ltbr（LTBR） | Ccl19-Cre による条件付き欠失（BAC トランスジェニックの Ccl19-Cre で間葉系間質細胞・FRC 系譜に限定） | マウスリンパ節の線維芽性細網細胞（FRC）前駆細胞 | LTBR 欠失により FRC が「筋線維芽細胞様前駆細胞」の段階で成熟が止まり、恒常性ケモカイン（CCL19, CCL21 等）と IL-7 の産生を欠く。T/B 域の区画化自体は保たれるが、免疫応答に必要な成熟 FRC 表現型を欠く | Chai Q, Onder L, Scandella E, et al. (Ludewig B 研究室). *Immunity*. 2013;38(5):1013-1024. PMID: 23623380 |

**判定への意味**: #4 は「必要性（necessity）」の証拠であり、「十分性（sufficiency、過剰発現で転換できる）」の証拠ではない。陽性対照として使う場合は、#1 の LTBR 経路刺激（十分性）と組み合わせて解釈する必要がある。

### ご指定の具体的候補についての確認結果

| 候補 | 確認結果 |
|---|---|
| LTBR / LTA / LTB / TNFSF14(LIGHT) の系 | LTBR・LTA/LTB は上記 #1〜#3 で十分性・#4 で必要性を確認。**TNFSF14(LIGHT) 単独でのfibroblast→LTo転換を示す一次文献は検索範囲では見つからなかった（未確認）** |
| 非古典的 NF-κB: RELB, NFKB2, MAP3K14(NIK), CHUK | RELB・NFKB2 は #1（Bénézech 2012）で LTBR 下流の必須経路として関与が確認された（ただし RELB/NFKB2 自体を単独で過剰発現させた実験ではなく、LTBR 刺激後の下流経路としての関与）。**MAP3K14(NIK) 単独の異所性過剰発現で CCL19/CCL21 が誘導されたという直接の一次文献は見つからなかった（未確認）**。CHUK(IKKα) も同様に未確認 |
| 転写因子: NR4A1, FOXL2, TBX5, TCF21, TWIST1, PRRX1, RUNX1 | いずれも「線維芽細胞を LTo 様に転換させた」という一次文献は見つからなかった（未確認・該当なし） |
| ID2 / RORC | ご指摘の通り ILC3（リンパ組織誘導細胞, LTi）側の転写因子であり、線維芽細胞（間質側）の転換因子としての報告はない（該当なし） |
| TNFAIP3 の喪失 | 線維芽細胞での TNFAIP3(A20) 欠失が LTo 様転換を誘導したという一次文献は見つからなかった（未確認）。A20 欠失患者由来線維芽細胞は TNF-α 応答が増強するという報告はあるが、CCL19/21/13 誘導やオーガナイザー転換の記載はない |
| IL17/IL22 シグナル | Nayar S, et al. *PNAS*. 2019;116(27):13490-13497. PMID: 31213547 — 「免疫線維芽細胞（immunofibroblast）」が TLS 形成を駆動し、既存の組織常在線維芽細胞が **IL-13 および IL-22 依存的に拡大・免疫学的形質を獲得する**ことを示した。ただし単一遺伝子の過剰発現による細胞状態の転換ではなく、サイトカイン刺激による細胞集団の拡大・活性化のモデルである点に注意 |
| TGFB, PDGFRB 系 | 線維芽細胞から LTo 様への転換を直接示す一次文献は見つからなかった（未確認） |
| 構成的活性型 NIK 過剰発現で CCL19/CCL21 が出るという報告 | **見つからなかった**。検索した範囲では、NIK は非古典的 NF-κB の下流エフェクターとして間接的に言及されるのみで、線維芽細胞への NIK 異所性発現実験そのものは確認できていない |

---

## Q2: 線維芽細胞／間葉系細胞で単一遺伝子過剰発現による明確な細胞状態転換が確立している例

Q1 で「線維芽細胞を単一遺伝子でLTo様に転換した」直接の実験系が見つからなかったため、代替として、間葉系細胞の単一遺伝子転換（transdifferentiation）の確立例を挙げる。これらは「陽性対照が働くかどうか」を確認するための手法的な対照（別の細胞状態への転換が起きるかを見る参照点）として使える可能性がある。

| 遺伝子 | 転換元 → 転換先 | 効率・根拠 | 出典 |
|---|---|---|---|
| MYOD1（MyoD） | 線維芽細胞（マウス C3H10T1/2 細胞株など） → 筋芽細胞 | 単一の cDNA のトランスフェクションのみで、筋特異的遺伝子群の活性化と筋芽細胞様の形態・分化を誘導。分化転換の古典的な最初の実証例 | Davis RL, Weintraub H, Lassar AB. *Cell*. 1987;51(6):987-1000. PMID: 3690668 |
| PPARG2（PPARγ2） | 線維芽細胞（NIH 3T3 など） → 脂肪細胞様 | レトロウイルスによる PPARγ2 の強制発現のみで脂肪分化を刺激。リガンド添加で用量依存的に分化が促進 | Tontonoz P, Hu E, Spiegelman BM. *Cell*. 1994;79(7):1147-1156. PMID: 8001151 |
| CBFA1/OSF2（RUNX2） | 非骨芽細胞（線維芽細胞を含む非骨芽細胞系） → 骨芽細胞様 | 非骨芽細胞への強制発現で、骨芽細胞特異的遺伝子群（オステオカルシン等）の主要なものを誘導 | Ducy P, Zhang R, Geoffroy V, Ridall AL, Karsenty G. *Cell*. 1997;89(5):747-754. PMID: 9182762 |

補足: ご提示の PPARG+CEBPA（2遺伝子）や SPI1 単独での骨髄系転換は、時間の制約により一次文献の確認まで至っていない（未確認）。単一遺伝子で確立している例としては上記 3 件が最も一次情報で裏付けが強い。

---

## Q3: 三次リンパ組織オーガナイザー様線維芽細胞（FRC 系列を含む）のマーカー遺伝子（CCL19/CCL21/CXCL13 以外）

| 遺伝子 | (i) FRC/LTo での発現の根拠 | (ii) マーカーか、分化を駆動する証拠があるか | 出典 |
|---|---|---|---|
| PDPN（ポドプラニン） | FRC の代表的な陽性マーカー（CD45陰性・CD31陰性・PDPN陽性で FRC を定義） | マーカーとしての記載が中心。分化を駆動する直接証拠は見当たらない | Fletcher AL, Acton SE, Knoblich K. *Nat Rev Immunol*. 2015;15(6):350-361. PMID: 25998961 |
| ICAM1, VCAM1 | LTBR 刺激後の前駆脂肪細胞で CCL19/21/13 とともに誘導される接着分子 | LTBR 刺激下流の誘導遺伝子として記載。分化そのものを駆動するという証拠はなし（マーカー的記載） | Bénézech C, et al. *Immunity*. 2012;37(4):721-734. PMID: 22940098 |
| MADCAM1 | 辺縁部細網細胞（marginal reticular cell, MRC）が VCAM1+ICAM1+MADCAM1+ で、CXCL13 産生・LTo 由来と考えられる | 「LTo から直接連続する子孫」という系譜的な記載があり、マーカーというより **LTo 系譜の同定基準**として使われている。分化を駆動する遺伝子としての機能実証ではない | Katakai T, et al. *J Immunol*. 2008;181(9):6189-6200. PMID: 18941209 / Katakai T. *Front Immunol*. 2012;3:200. PMID: 22807928 |
| RANKL（TNFSF11） | MRC・胎生期 LTo 細胞で発現し、LTi 細胞との相互作用に関わるとされる | 上記と同じ文献群でマーカー的に記載。線維芽細胞側の分化を駆動する証拠ではなく、免疫細胞側（LTi, FDC 前駆細胞）への作用が主眼 | 同上（PMID: 18941209, 22807928） |
| CR2（補体受容体2） | **注意**: CR2 は古典的には**濾胞樹状細胞（FDC）**のマーカーであり、FRC/LTo 系譜のマーカーとして一次文献で確認できなかった（未確認・別系譜の可能性が高い） | — | 未確認 |
| PI16, CXCL12 | 組織横断的な線維芽細胞アトラス（17 組織・約 23 万細胞、ヒト・マウス）で、Pi16 と Col15a1 が「普遍的（universal）」な 2 大線維芽細胞サブタイプのマーカーとして同定された。CXCL12 は多くの間質・FRC 系文献で恒常性ケモカインとして併記される | この論文自体は横断的な分類・記載研究であり、Pi16 陽性細胞が LTo/FRC へ分化することを駆動的に示した実験ではない（マーカー・分類としての記載） | Buechler MB, Pradhan RN, Krishnamurty AT, et al. *Nature*. 2021;593(7860):575-579. PMID: 33981032 |
| TNFSF13B(BAFF), GREM1 | 検索範囲では FRC/LTo との関連を直接確認できる一次文献に到達できなかった（未確認）。時間の制約により確認を打ち切った | 未確認 | 未確認 |

---

## まとめ（確認できたことと推測の区別）

**確認できたこと（一次文献あり）**
- LTBR へのアゴニスト刺激により、間葉系前駆細胞（前駆脂肪細胞）が CCL19/CCL21/CXCL13/ICAM1/VCAM1 陽性の LTo 様間質細胞へ転換する（PMID 22940098）。
- LTBR を線維芽細胞系譜（Ccl19-Cre）で欠失させると、FRC の成熟が筋線維芽細胞様前駆段階で止まり、CCL19/CCL21/IL-7 産生を失う（PMID 23623380）。
- 異所性 LTA/LTAB 発現（上皮細胞発現、パラクライン作用）で非リンパ組織にリンパ組織様構造が新生する（PMID 8666904, 12732657）。
- 単一遺伝子過剰発現による間葉系細胞の明確な転換は、筋（MYOD1, PMID 3690668）、脂肪（PPARG2, PMID 8001151）、骨芽細胞（RUNX2/CBFA1, PMID 9182762）で確立している。

**推測・未確認（本文中に明記済み）**
- 「NIK（MAP3K14）の構成的活性型を線維芽細胞に過剰発現させると CCL19/CCL21 が誘導される」という趣旨の一次文献は見つからなかった。
- RELB, NFKB2, CHUK, TNFAIP3, NR4A1, FOXL2, TBX5, TCF21, TWIST1, PRRX1, RUNX1, TNFSF14(LIGHT) の単独過剰発現／欠失による線維芽細胞→LTo様転換は、検索範囲では確認できなかった。
- CR2, TNFSF13B(BAFF), GREM1 の FRC/LTo での役割は未確認、または別系譜（FDC）である可能性が高い。

**陽性対照としての示唆（判定は PI に委ねる）**
- 「十分性」の陽性対照候補としては、現時点で一次文献の裏付けが最も強いのは **LTBR 経路への刺激（アゴニスト、または LTA/LTB 発現)** であり、腎線維芽細胞に対する in silico perturbation の陽性対照として使う場合は「単一遺伝子の過剰発現型摂動」というより「受容体シグナル経路の活性化」を模す形になる点に注意が必要。
- 「必要性」の対照としては LTBR 自体（削除型摂動）が使える可能性がある。
- NIK 単独過剰発現を陽性対照に使う計画がある場合、現時点でその十分性を直接支持する一次文献は見つかっていない。計画段階でこの点を踏まえた判断が必要。

---

## 出典一覧（PMID 確認済み）

1. Bénézech C, et al. Lymphotoxin-β receptor signaling through NF-κB2-RelB pathway reprograms adipocyte precursors as lymph node stromal cells. Immunity. 2012;37(4):721-734. PMID: 22940098
2. Kratz A, Campos-Neto A, Hanson MS, Ruddle NH. Chronic inflammation caused by lymphotoxin is lymphoid neogenesis. J Exp Med. 1996;183(4):1461-1472. PMID: 8666904
3. Drayton DL, Ying X, Lee J, Lesslauer W, Ruddle NH. Ectopic LTαβ directs lymphoid organ neogenesis with concomitant expression of peripheral node addressin and a HEV-restricted sulfotransferase. J Exp Med. 2003;197(9):1153-1163. PMID: 12732657
4. Chai Q, Onder L, Scandella E, et al. Maturation of lymph node fibroblastic reticular cells from myofibroblastic precursors is critical for antiviral immunity. Immunity. 2013;38(5):1013-1024. PMID: 23623380
5. Davis RL, Weintraub H, Lassar AB. Expression of a single transfected cDNA converts fibroblasts to myoblasts. Cell. 1987;51(6):987-1000. PMID: 3690668
6. Tontonoz P, Hu E, Spiegelman BM. Stimulation of adipogenesis in fibroblasts by PPARγ2, a lipid-activated transcription factor. Cell. 1994;79(7):1147-1156. PMID: 8001151
7. Ducy P, Zhang R, Geoffroy V, Ridall AL, Karsenty G. Osf2/Cbfa1: a transcriptional activator of osteoblast differentiation. Cell. 1997;89(5):747-754. PMID: 9182762
8. Fletcher AL, Acton SE, Knoblich K. Lymph node fibroblastic reticular cells in health and disease. Nat Rev Immunol. 2015;15(6):350-361. PMID: 25998961
9. Katakai T, et al. Organizer-like reticular stromal cell layer common to adult secondary lymphoid organs. J Immunol. 2008;181(9):6189-6200. PMID: 18941209
10. Katakai T. Marginal reticular cells: a stromal subset directly descended from the lymphoid tissue organizer. Front Immunol. 2012;3:200. PMID: 22807928
11. Nayar S, et al. Immunofibroblasts are pivotal drivers of tertiary lymphoid structure formation and local pathology. PNAS. 2019;116(27):13490-13497. PMID: 31213547
12. Buechler MB, Pradhan RN, Krishnamurty AT, et al. Cross-tissue organization of the fibroblast lineage. Nature. 2021;593(7860):575-579. PMID: 33981032

（この一覧は WebSearch で得たタイトル・巻号・PMID をもとにしており、フルテキストの詳細な方法論までは一部の論文で未確認。特に #2, #3 は CCL19/CCL21/CXCL13 の直接測定について要旨レベルの確認にとどまる。）

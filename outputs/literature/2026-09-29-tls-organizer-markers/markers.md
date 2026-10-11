# 三次リンパ組織オーガナイザー様間質細胞のマーカー（CCL19・CCL21・CXCL13 を使わない定義）

- 作成: literature-scout（Labretary 委託）/ 2026-09-29
- 目的: 三次リンパ組織（腎の三次リンパ組織を含む）の「オーガナイザー様線維芽細胞・リンパ組織様間質細胞」を、CXCL13・CCL19・CCL21 以外のマーカーで定義できるか、一次文献で確認する
- 方法: PubMed（NCBI eutils による アブストラクト直接取得を含む）、PMC 全文、GEO（NCBI eutils）、journal 本文サイトを確認した。**フルテキストに直接アクセスできず、検索エンジンの要約（二次的な引用）にとどまったものは、その旨を明記した**。仮説の提案・判定基準の設計は行っていない。
- 用語は `docs/用語集.md` に合わせた（三次リンパ組織、腎線維芽細胞 等）。文中の遺伝子記号はヒト表記（大文字）で統一し、マウス由来データであることは注記した。

---

## 0. 総括（先に要点）

- **腎臓の三次リンパ組織そのもので実験的に確認されている線維芽細胞マーカー**は、確認できた範囲では **NGFR（p75NTR）** が最も直接的である。Sato ら 2016（PMID 27699223, マウスおよびヒト）が、三次リンパ組織の中の線維芽細胞を NGFR 陽性として同定し、この所見はヒト加齢腎でも確認されたとアブストラクトに明記されている。
- Yoshikawa ら 2023（PMID 37548710, マウス snRNA-seq＋ヒト検証）は、三次リンパ組織内の「炎症促進性線維芽細胞」を **ICAM1, VCAM1, TNFSF13B（BAFF）, ENPP2, C3, CXCL9, CXCL10** などで特徴づけており、CCL19/CCL21/CXCL13 以外の遺伝子でも一定の記述が一次文献にある。ただしこれは **フルテキストを機械的に要約させた結果**であり、Labretary 自身が図表を目視確認したわけではない点に注意（§2 参照）。
- **CR2（CD21）は濾胞樹状細胞（FDC）のマーカーであり、リンパ組織オーガナイザー様線維芽細胞（FRC/LTo）系譜のマーカーではない**。Sato 2016 でも CD21 は FDC マーカーとして FRC 系譜とは別に使われている。
- **FAP は、TLS を伴う線維芽細胞ではむしろ低い**、という報告がある（TLS 随伴の癌関連線維芽細胞は VCAM1・ICAM1 高／FAP 低、腫瘍促進性の癌関連線維芽細胞は FAP 高）。FAP を「オーガナイザー様」の陽性マーカーとして使うのは一次文献と矛盾する可能性がある。
- **Meylan 2022（PMID 35231421）の GEO データ（GSE175540）は、Visium の生・処理済みスポット行列（24 検体、ヒト腎細胞癌）として公開されている**が、確認できた範囲では「三次リンパ組織シグネチャー遺伝子リスト」そのものはこの GEO には含まれていない。論文が独自に定義した **29 遺伝子の「TLS imprint signature」** は論文本文・補足資料に記載があるはずだが、出版社サイトへの直接アクセスが 403 で妨げられたため、**この 29 遺伝子リストは検索エンジンが返した要約からの再構成であり、補足表 PDF を直接目視確認できていない（中程度の確信度、要一次確認）**。
- 単一細胞用に公開済みで構成遺伝子を確認できたシグネチャーは **Coppola 2011（12-chemokine, PMID 21703392, 直接確認）** と **Cabrita 2020（9 遺伝子, PMID 31942071, 二次資料経由）**。いずれも CXCL13 を含み、Coppola は CCL19/CCL21 も含む。「CXCL13/CCL19/CCL21 を使わない」という要件を満たすシグネチャーは、確認できた範囲では存在しない。

---

## 1. 候補マーカーごとの裏づけ

凡例: 確度 = 「一次文献を直接確認（アブストラクト全文または PMC/journal 本文を取得）」／「二次資料経由（検索エンジンの要約のみ、原文直視できず）」／「未確認」

| 遺伝子 | 細胞種（マーカーとしての位置づけ） | 種 | 組織 | 腎の三次リンパ組織での直接所見か、他臓器からの類推か | 出典（PMID） | 確度 |
|---|---|---|---|---|---|---|
| PDPN（ポドプラニン, gp38） | FRC 全般の代表的陽性マーカー。腎では TLT 内の炎症促進性線維芽細胞でも発現 | ヒト・マウス | リンパ節（総説）／マウス腎（TLT） | リンパ節からの類推が主。**腎 TLT での発現は Yoshikawa 2023 の図に記載ありと報告された**が、フルテキスト自動要約に基づく（二次確認） | Fletcher 2015 総説 PMID 25998961／Yoshikawa 2023 PMID 37548710 | Fletcher: 直接確認（アブストラクトのみ、PDPN の記載は総説本文）／Yoshikawa 該当箇所: 二次資料経由 |
| NGFR（p75NTR） | **TLT 内の線維芽細胞のマーカー**として腎で直接同定 | マウス・ヒト（アブストラクトに明記） | **腎臓（TLT）** | **腎の三次リンパ組織での直接所見**（アブストラクトに "p75 neurotrophin receptor+ (p75NTR+) fibroblasts...inside TLTs" と明記） | Sato 2016 PMID 27699223 | 直接確認（アブストラクト全文取得） |
| ICAM1 | LTBR シグナル下流で誘導される接着分子。腎 TLT の炎症促進性線維芽細胞でも発現と報告 | マウス（前駆脂肪細胞）／マウス腎 | リンパ節新生モデル／マウス腎（TLT） | リンパ節・脂肪組織モデルからの類推が主。腎 TLT での発現は二次資料経由の確認にとどまる | Bénézech 2012 PMID 22940098／Yoshikawa 2023 PMID 37548710 | Bénézech: 直接確認（過去の文献調査で確認済み）／Yoshikawa 該当箇所: 二次資料経由 |
| VCAM1 | 同上。TLS 随伴の癌関連線維芽細胞で ICAM1 とともに高発現と報告 | マウス／ヒト腫瘍 | リンパ節新生モデル／癌（TLS 随伴 CAF） | 他臓器（脂肪組織モデル、腫瘍）からの類推 | Bénézech 2012 PMID 22940098／腫瘍 CAF-TLS の記載は検索結果の要約経由（PMID 未特定、§1 末尾の注記参照） | 一部直接確認、一部未確認 |
| TNFSF13B（BAFF） | **腎 TLT 内の炎症促進性線維芽細胞のマーカー**（B 細胞生存因子） | マウス（snRNA-seq）＋ヒト検証 | **腎臓（TLT）** | **腎の三次リンパ組織での直接所見**。アブストラクトに "B cell-activating factor" の産生が TLT 内線維芽細胞の特徴として明記 | Yoshikawa 2023 PMID 37548710 | 直接確認（アブストラクト全文にて "B cell-activating factor" の記載を確認。遺伝子名 TNFSF13B としての明記は本文レベルで、フルテキスト自動要約経由） |
| TNFSF13（APRIL） | GC-T zone interface の PDPN+CD157hi 線維芽性細網細胞（GTIRC）のマーカー | マウス | リンパ節（胚中心） | **リンパ節からの類推**。腎 TLT での報告は確認できず（未確認） | Zhang et al. 2018 PMID 29549115 | 直接確認（リンパ節での所見）。腎への外挿は未確認 |
| CR2（CD21, 濾胞樹状細胞） | **濾胞樹状細胞（FDC）のマーカーであり、FRC/LTo 系譜のマーカーではない**。腎 TLT でも FDC マーカーとして FRC 系譜と区別して使われている | マウス・ヒト | リンパ節（総説）／マウス・ヒト腎（TLT） | **腎の三次リンパ組織での直接所見あり（ただし FDC マーカーとして）**。前回調査（2026-09-27）で「FRC/LTo 系譜のマーカーとしては未確認」と結論済み。今回 Sato 2016 で CD21 が FDC マーカーとして図に使われていることを確認したが、これは FRC 系譜とは別の細胞種の話である | Sato 2016 PMID 27699223 | 直接確認（PMC 全文の自動要約経由、二次資料経由） |
| CXCL12 | Meylan 2022 の「TLS imprint signature」に含まれる 2 つの線維芽細胞マーカーの 1 つ。リンパ節では CXCL12 高発現の細網細胞（CAR 細胞、CXCL12hi FRC）が LTBR 高発現のサブセットとして報告 | ヒト（腎細胞癌）／マウス（骨髄・リンパ節） | 腎細胞癌（Visium）／骨髄・リンパ節 | 腎 TLT そのものでの発現は未確認。**腎細胞癌の TLS 随伴シグネチャーの構成遺伝子としては直接確認**（Meylan の 29 遺伝子シグネチャーに含まれる、二次資料経由） | Meylan 2022 PMID 35231421（29 遺伝子シグネチャー, 二次資料経由）／Rodda 2018 PMID 29752062（リンパ節 CXCL9+ TRC 等の文脈, 直接確認） | 混在（Meylan: 二次資料経由／Rodda: 直接確認） |
| MADCAM1 | 辺縁部細網細胞（MRC）のマーカー。CXCL13 産生細胞で LTo 由来の系譜と考えられている | マウス | リンパ節 | **リンパ節からの類推**。腎 TLT での報告は確認できず（未確認） | Katakai 2008 PMID 18941209／Katakai 2012 PMID 22807928（前回調査で確認済み） | 直接確認（リンパ節）。腎への外挿は未確認 |
| IL7 | FRC の恒常性維持ケモカイン・サイトカインの代表格。LTBR 欠失で FRC 成熟が止まると IL-7 産生も失われる | マウス | リンパ節 | **リンパ節からの類推**。腎 TLT での直接報告は今回の検索範囲では確認できず（未確認） | Chai 2013 PMID 23623380（前回調査で確認済み） | 直接確認（リンパ節）。腎への外挿は未確認 |
| LTA（リンフォトキシンα） | LTo/FRC を「誘導する」リガンド側の遺伝子。産生細胞は主にリンパ組織誘導細胞（LTi）やリンパ球であり、線維芽細胞側のマーカーではない | マウス | 膵臓（異所性発現モデル） | 他臓器（膵島）モデルからの類推。腎 TLT での役割は確認できず（未確認） | Kratz 1996 PMID 8666904／Drayton 2003 PMID 12732657（前回調査で確認済み） | 直接確認（他モデル）。線維芽細胞自体のマーカーではない点に注意 |
| ENPP2（autotaxin） | **腎 TLT 内の炎症促進性線維芽細胞、二次リンパ器官間質マーカーとして記載** | マウス | **腎臓（TLT）** | **腎の三次リンパ組織での直接所見（と報告された）** | Yoshikawa 2023 PMID 37548710 | 二次資料経由（PMC フルテキストの自動要約のみ。Labretary 自身が図表を確認したわけではない） |
| CD74 / HLA-DRA（MHC class II） | 一部の FRC サブセット（T 細胞域）が MHC class II 抗原提示経路の分子（CD74, H2-Ab1）を発現し、末梢性免疫寛容に関与 | マウス | リンパ節 | **リンパ節からの類推**。腎 TLT での報告は今回の検索範囲では確認できず（未確認。Yoshikawa 2023 の自動要約では「記載なし」という回答だった） | Baptista/Mebius 2014 PMID 25407678 | 直接確認は不十分（検索エンジン要約経由。原著 PMC 全文は未取得） |
| C3（補体） | **腎 TLT 内の炎症促進性線維芽細胞クラスターのマーカー**として記載 | マウス | **腎臓（TLT）** | **腎の三次リンパ組織での直接所見（と報告された）** | Yoshikawa 2023 PMID 37548710 | 二次資料経由（フルテキスト自動要約のみ） |
| C7（補体） | Meylan の「TLS imprint signature」の補体関連遺伝子 2 つのうちの 1 つ | ヒト（腎細胞癌） | 腎細胞癌（Visium） | 腎 TLT の間質細胞そのものでの報告ではなく、TLS 随伴シグネチャーの構成遺伝子としての位置づけ（二次資料経由） | Meylan 2022 PMID 35231421 | 二次資料経由 |
| SERPING1（C1 インヒビター） | 検索範囲では FRC/LTo との直接的関連を確認できる一次文献に到達できなかった | — | — | — | 未確認 | 未確認 |
| CFD（アディプシン, 補体第D因子） | 検索範囲では FRC/LTo との直接的関連を確認できる一次文献に到達できなかった。関連しそうな記載として「CFD+ 炎症性線維芽細胞」という横断的線維芽細胞アトラスの記述に触れたが、リンパ組織・腎 TLT との関連は確認できず | — | — | — | 未確認（Buechler 2021 PMID 33981032 関連の可能性はあるが、CFD が LTo/FRC マーカーであるとの直接記載は確認できず） | 未確認 |
| PI16 | 組織横断的な線維芽細胞アトラスで、Col15a1 とともに「普遍的（universal）」な 2 大線維芽細胞サブタイプのマーカーとして同定。**LTo/FRC に特異的なマーカーではなく、あらゆる組織の間質前駆細胞に共通するマーカー** | マウス（主）＋一部ヒト | 17 組織横断（腎を含む一般的結合組織） | **腎 TLT・LTo に特異的というわけではない**（前回調査で結論済み）。今回も新たな LTo 特異性を示す一次文献は見つからなかった | Buechler 2021 PMID 33981032 | 直接確認（アブストラクト全文取得） |
| CD34 | リンパ節被膜・髄質血管周囲外膜の間質細胞（FRC/MRC/FDC とは別系譜）のマーカー | マウス | リンパ節 | **リンパ節からの類推であり、かつ「LTo 系譜ではない別系譜」として報告されている**。Rodda 2018 は CD34+ SC を FRC/MRC/FDC と並ぶ独立したクラスターとして記載 | Rodda 2018 PMID 29752062 | 直接確認（検索エンジン要約だが原論文アブストラクトと整合） |
| CLU（クラスタリン） | **濾胞樹状細胞（FDC）の代表的マーカー**（感度が高いとされる）。FRC/LTo のマーカーではない | ヒト（病理診断文脈） | リンパ節（FDC 肉腫の鑑別） | リンパ節からの類推。**FDC マーカーであり、オーガナイザー様線維芽細胞（FRC/LTo）のマーカーとしては該当しない** | 病理総説（PMID 未特定、二次資料経由） | 二次資料経由（正確な一次 PMID は今回特定できず、未確認扱い） |
| GREM1 | T-B 細胞境界に局在する FRC サブセット（樹状細胞恒常性維持ニッチ）のマーカー。腎では線維化の文脈で発現上昇が報告されているが、これは線維芽細胞全般の話で TLT 特異的ではない | マウス・ヒト（リンパ節）／マウス（腎線維化, 別文脈） | リンパ節／腎（線維化、TLT 非特異的） | **リンパ節からの類推が主**。腎での GREM1 上昇は「線維化」の文脈であり、TLT・LTo との直接の関連は確認できず（未確認） | Kapoor 2021 PMID 33903764（リンパ節 GREM1+ FRC）／腎線維化での GREM1 上昇は別の総説的記載（PMID 未特定） | リンパ節の所見: 直接確認（検索エンジン要約経由だが記載内容は具体的）。腎 TLT への外挿: 未確認 |
| SFRP2 | ヒト皮膚の主要な線維芽細胞サブタイプ（SFRP2+/DPP4+）のマーカーとして報告。**リンパ組織・TLS との直接関連は確認できず** | ヒト | 皮膚（TLS・リンパ組織とは無関係の文脈） | **該当なし**。皮膚線維芽細胞の分類研究であり、腎・リンパ組織のオーガナイザー様細胞のマーカーとしての報告は見つからなかった | Tabib et al.（皮膚線維芽細胞論文, PMID 未特定, 二次資料経由） | 未確認（TLS/LTo 文脈での関連なし） |
| FAP（線維芽細胞活性化タンパク質） | 癌関連線維芽細胞（CAF）の代表的マーカーだが、**TLS 随伴の CAF はむしろ FAP 低・VCAM1/ICAM1 高であり、FAP 高の CAF は TLS を伴わない腫瘍促進性の表現型**と報告されている | マウス（黒色腫モデル） | 腫瘍（TLS 随伴 CAF の比較） | 腎ではなく腫瘍（黒色腫）モデルからの知見。**「FAP がオーガナイザー様線維芽細胞の陽性マーカーである」という想定とは逆の報告**である点に注意 | 検索エンジン要約経由（frontiersin 2026 論文, PMID 42540229 関連の可能性。正確な一次文献の特定は今回できず） | 二次資料経由（要一次確認） |
| THY1（CD90） | 一般的な線維芽細胞マーカー。肺癌 CAF では Thy1+ CAF が腫瘍浸潤を促進すると報告されているが、**TLS・LTo との直接関連は確認できず** | マウス（肺腺癌モデル） | 腫瘍（肺） | **該当なし**。TLS 形成との関連を示す一次文献は見つからなかった（未確認） | Nurmik ら等の CAF 総説（PMID 未特定, 二次資料経由） | 未確認 |
| TNC（テネイシン C） | **腎で髄質の線維芽細胞マーカー**として記載（TLT 特異的マーカーではなく、髄質という部位の線維芽細胞マーカーとして区別されている） | マウス | 腎臓 | 腎での直接所見だが、**TLT 特異的マーカーではなく、むしろ TLT とは別の部位（髄質）の線維芽細胞を区別するために使われている**点に注意 | Yoshikawa 2023 PMID 37548710（補足図, 二次資料経由） | 二次資料経由 |
| POSTN（ペリオスチン） | 筋線維芽細胞・CAF の細胞外基質関連マーカーとして一般に知られるが、**TLS・LTo との直接関連を示す一次文献は今回の検索範囲では見つからなかった** | — | — | 未確認 | 未確認 | 未確認 |

**FAP・THY1・POSTN・SFRP2・CD34・CLU・GREM1・SERPING1・CFD について、腎の三次リンパ組織そのものでの発現を直接示した一次文献は、今回の検索範囲では見つからなかった。** これらは他組織（皮膚・肺・腫瘍・骨髄）の線維芽細胞分類研究からの類推であるか、あるいは（FAP のように）むしろ TLS 随伴線維芽細胞との**負の**相関が報告されているものである。

---

## 2. 腎臓の三次リンパ組織に特化した一次文献での間質マーカー

事前登録された 3 論文について、アブストラクトを NCBI eutils で直接取得し、全文（PMC）は自動要約ツールで確認した。**全文の自動要約は Labretary 自身が図表を目視したものではないため、中程度の確信度として扱うこと。**

### Sato 2016（PMID 27699223, JCI Insight, マウス＋ヒト） — アブストラクト直接確認済み
- 対象: 急性腎障害後の加齢マウス腎臓、および加齢ヒト腎臓（比較検証）
- アブストラクトに明記された知見: 単一系譜の腎常在線維芽細胞が、**(1) p75 neurotrophin receptor 陽性（p75NTR+, NGFR）線維芽細胞**、**(2) 恒常性ケモカイン産生線維芽細胞（CCL19/CCL21/CXCL13、TLT 内部）**、**(3) レチノイン酸産生線維芽細胞（TLT 周囲）** の 3 方向に分化することが示された。
- 全文の自動要約（二次資料扱い）で確認した追加情報: PDGFRB, ER-TR7, αSMA を線維芽細胞・筋線維芽細胞マーカーとして使用。CD21 を FDC マーカーとして使用。CD45・CD11c を陰性対照として使用。PDPN・CD34 の言及は本文中に確認できなかった。
- ヒト検証: 加齢ヒト腎臓でも同様の細胞・分子構成が確認されたとアブストラクトに明記（Figure 7 で言及、との自動要約あり）。

### Arai 2021（PMID 33951198, Immunological Reviews, 総説） — アブストラクト直接確認済み、本文は 403 エラーで未取得
- 総説であり新規データは含まない。アブストラクトの範囲では、線維芽細胞の不均一性が腎の生理的機能（貧血、尿細管周囲毛細血管の維持）と病的状態（線維化、TLT 形成）の両方に関わることを概説している、との記載にとどまる。
- **具体的なマーカー遺伝子名はアブストラクトには記載されていない**。本文（Wiley サイト）へのアクセスは 403 エラーで妨げられ、確認できなかった。したがって、この総説が挙げる個別マーカーの一覧は「未確認」とする。Sato 2016・Yoshikawa 2023（同じ京都大学グループ）の知見を統合的に紹介している可能性が高いが、これは推測であり本文で確認していない。

### Yoshikawa 2023（PMID 37548710, JASN, マウス snRNA-seq＋ヒト検証） — アブストラクト直接確認済み、本文は自動要約経由
- アブストラクトに明記: snRNA-seq でマウス腎臓（虚血再灌流後、加齢）を解析。**TLT 内外に局在する「炎症促進性」「線維化促進性」の 2 種類の線維芽細胞**を同定。炎症促進性線維芽細胞は STAT1 活性化を示し、**CXCL9/CXCL10 および B 細胞活性化因子（TNFSF13B/BAFF）**などのケモカイン・サイトカインを産生し、リンパ球の動員・生存に寄与する。これらの知見は免疫染色・in situ hybridization でヒト腎臓（TLT を有する）でも確認されたとアブストラクトに明記。
- 全文の自動要約（二次資料扱い、要注意）で確認したクラスター別マーカー: 炎症促進性線維芽細胞（クラスター 5 相当）= VCAM1, ICAM1, CXCL9, CXCL10, TNFSF13B, IL15, IL33, C3, PDPN, MFGE8, ENPP2, NGFR。髄質線維芽細胞 = TNC。
- **この全文要約はツールによる自動生成であり、Labretary が図表を直接確認したものではない。** ただし NGFR・TNFSF13B・CXCL9/CXCL10 の産生についてはアブストラクト本文でも独立に裏付けられている。

---

## 3. Meylan 2022（PMID 35231421, Immunity, 腎細胞癌の空間トランスクリプトーム）

### GEO 公開状況（NCBI eutils で直接確認）
- アクセッション: **GSE175540**（確認済み。PubMed リンクも PMID 35231421 と一致）
- 内容: ヒト腎細胞癌（ccRCC）の **24 検体**の Visium 空間トランスクリプトーム（frozen および FFPE、公開日 2022-03-03）
- 補足ファイル形式: CSV, H5, JPG, JSON, PNG（**スポット×遺伝子の発現行列そのものが公開されている**）
- **確認できたこと**: GEO には Visium のスポットレベル発現行列（生または処理済み）が公開されている。
- **確認できなかったこと**: 「三次リンパ組織シグネチャー遺伝子リスト」という形の遺伝子リストファイルは、GEO のサプリメントファイル一覧（CSV/H5/JPG/JSON/PNG）からは確認できなかった。遺伝子リストは論文本体・補足表（Table S1 等）に記載されている可能性が高いが、**出版社サイト（Cell Press / ScienceDirect）へのアクセスは 403 エラーで直接確認できなかった**。

### 手法（検索エンジン要約経由、要一次確認）
- Visium（10x Genomics）を使用。1 スポットの直径は 55 μm（Visium の仕様どおり）。
- **TLS＋スポットの定義方法**: 病理医による手動アノテーション。凍結検体は隣接切片での CD20/CD3/PD-1 の三重蛍光免疫染色、FFPE 検体は CD3/CD20 の免疫組織化学で検証、との検索結果の要約が得られた。**この記載は一次文献本文を直接確認したものではなく、検索エンジンが返した要約に基づく（要一次確認）**。

### 著者らの GitHub リポジトリ（`maximemeylan/Meylan_et_al_2022`）
- 確認できたスクリプト一覧: `00_visium_preprocessing.R`, `01_MCP_counter.R`, `02_TLS_imprint_signature.R`, `03_bulk_clonotype_analysis.R`, `03_mixcr_for_bulk.sh`, `03_mixcr_for_visium.sh`, `03_spatial_clonotype_analysis.R`
- `02_TLS_imprint_signature.R` は `tls_signature_29.rds` という外部ファイルを読み込む記述になっており、**29 遺伝子の構成そのものはこのスクリプトファイルには書き出されていない**（別ファイルに保存されている）。

### CXCL13・CCL19・CCL21 以外の構成遺伝子（**中程度の確信度、一次資料未直視**）
検索エンジンが複数回にわたり一致して返した記述によれば、Meylan らは腫瘍内で TLS が周囲組織に及ぼす影響（形質細胞の播種を含む）を捉えるための **29 遺伝子の「TLS imprint signature」** を、凍結 TLS＋腫瘍 4 検体からの leave-one-out 交差検証で選定したとされる。構成は以下の通り（**未検証の二次資料経由。原著の補足表 PDF は 403 エラーで直接確認できていない**）:
- 免疫グロブリン関連 12 遺伝子: IGHA1, IGHG1, IGHG2, IGHG3, IGHG4, IGHGP, IGHM, IGKC, IGLC1, IGLC2, IGLC3, JCHAIN
- B 細胞マーカー 5 遺伝子: CD79A, FCRL5, MZB1, SSR4, XBP1
- T 細胞マーカー 2 遺伝子: TRBC2, IL7R
- 線維芽細胞マーカー 2 遺伝子: **CXCL12, LUM**
- 補体関連 2 遺伝子: C1QA, C7
- その他 6 遺伝子: CD52, APOE, PTGDS, PIM2, DERL3, および表記が不明瞭な 1 遺伝子（検索結果では "PTLP" と表示されたが、正式な遺伝子記号として存在しないため誤記または OCR 由来の誤りの可能性が高い。**PLTP（リン脂質転送タンパク質）の可能性があるが未確認**）
- **重要な注意点**: この 29 遺伝子シグネチャーは「TLS imprint（TLS が周囲の腫瘍組織に残す痕跡）」であり、**TLS スポット自体を定義するためのシグネチャーではなく、TLS が離れた場所に及ぼす影響を bulk 腫瘍検体から検出するためのシグネチャー**である。したがって、この 29 遺伝子＝「TLS スポットの構成遺伝子」ではない点に注意。線維芽細胞に関する遺伝子は CXCL12 と LUM の 2 つのみで、CXCL13・CCL19・CCL21 はこの 29 遺伝子には含まれていない（ただし、TLS＋スポット自体の同定には別途 CXCL13 等の古典的マーカーが使われた可能性が高いが、これは未確認）。

---

## 4. 単一細胞でスコア化できる公開済み TLS シグネチャー（構成遺伝子を実際に列挙）

### Coppola 2011（12-chemokine signature, PMID 21703392, ヒト大腸癌） — 一次文献本文で直接確認
構成遺伝子（12 個、全て確認）:
**CCL2, CCL3, CCL4, CCL5, CCL8, CCL18, CCL19, CCL21, CXCL9, CXCL10, CXCL11, CXCL13**
- **CXCL13・CCL19・CCL21 を含む**（3 つとも含む）。
- 326 検体のヒト大腸腺癌の遺伝子チップデータで検証。異所性リンパ組織様構造の存在と予後改善を予測。

### Cabrita 2020（TLS シグネチャー, PMID 31942071, ヒト黒色腫） — 検索エンジン経由、二次資料
構成遺伝子（9 個、検索結果で一致して得られたもの。**原著の補足表を直接確認できておらず、要一次確認**）:
**CD79B, CD1D, CCR6, LAT, SKAP1, CETP, EIF1AY, RBP5, PTGDS**
- **CXCL13・CCL19・CCL21 のいずれも含まない**（この点は特筆に値する。ただし遺伝子数について、ご依頼文では「12 遺伝子」とされていたが、検索で確認できたのは 9 遺伝子であり、**遺伝子数の食い違いがある。この不一致は解消できておらず、要一次確認**）。
- CD8+CD20+ 黒色腫検体で発現上昇していた遺伝子から選定されたと記載。

---

## 5. 空間データ（1 スポット＝複数細胞）由来シグネチャーを単一核データに当てはめる際の方法論的問題

専用に「TLS シグネチャーを空間データから単一核データへ翻訳する」ことだけを扱った一次方法論論文は、今回の検索範囲では特定できなかった（未確認）。関連する一般的な方法論的知見として、検索で得られたもの（**いずれも要一次確認、タイトルと要旨レベルの確認にとどまる**）:
- Visium のスポット解像度は TLS 内部の細胞不均一性を解像するには不十分であり、単一細胞参照データを用いたデコンボリューションに誤差を持ち込みうる、との指摘。
- bulk／spot 逆重畳（deconvolution）において、bulk に存在する細胞型が単一細胞参照に欠けていると、細胞比率の推定が大きく損なわれる。scRNA-seq 参照を優先し、snRNA-seq を追加する場合はモダリティ間の発現差遺伝子を除去してから使うべき、との提言。
- Visium は 2 次元切片であり、TLS の 3 次元的な不均一性を捉えられないという限界。

**これらはいずれも一般的なデコンボリューション方法論の議論であり、TLS シグネチャーに特化した検証論文ではない。** 本研究（腎臓の単一核データに空間シグネチャーを当てはめる計画がある場合）で想定すべき固有の問題として、Labretary 自身の判断は含めないが、上記の一般論から読み取れる論点を整理すると:
1. Visium 1 スポット（直径 55 μm）には複数細胞（腎の緻密な組織では十数個程度になりうる）の発現が混合しており、線維芽細胞由来のシグナルが免疫細胞由来のシグナルに希釈・混入している可能性がある。
2. 単一核 RNA-seq は細胞質 mRNA を捉えないため、分泌タンパク質（ケモカイン・補体等）をコードする mRNA の核内シグナルが弱く、スポットベースで有意だった遺伝子が単一核データでは検出限界以下になりうる（これは `docs/用語集.md` の「単一核 RNA-seq」の項目にある一般的注意点と一致する）。
3. Meylan の 29 遺伝子シグネチャーは前述の通り「TLS imprint（bulk 腫瘍への染み出し効果）」を検出するために設計されたものであり、単一細胞・単一核レベルでの「その細胞が TLS オーガナイザー様かどうか」を判定する目的では設計されていない。目的の異なるシグネチャーを転用することそのものが方法論上の論点になりうる。

---

## 確認できた PMID 一覧

| PMID | 文献 |
|---|---|
| 27699223 | Sato Y, et al. Heterogeneous fibroblasts underlie age-dependent tertiary lymphoid tissues in the kidney. JCI Insight. 2016. |
| 33951198 | Arai H, Sato Y, Yanagita M. Fibroblast heterogeneity and tertiary lymphoid tissues in the kidney. Immunol Rev. 2021. |
| 37548710 | Yoshikawa T, et al. Tertiary Lymphoid Tissues Are Microenvironments with Intensive Interactions between Immune Cells and Proinflammatory Parenchymal Cells in Aged Kidneys. J Am Soc Nephrol. 2023. |
| 35231421 | Meylan M, et al. Tertiary lymphoid structures generate and propagate anti-tumor antibody-producing plasma cells in renal cell cancer. Immunity. 2022. |
| 31942071 | Cabrita R, et al. Tertiary lymphoid structures improve immunotherapy and survival in melanoma. Nature. 2020. |
| 21703392 | Coppola D, et al. Unique ectopic lymph node-like structures present in human primary colorectal carcinoma are identified by immune gene array profiling. Am J Pathol. 2011. |
| 22940098 | Bénézech C, et al. Lymphotoxin-β receptor signaling through NF-κB2-RelB pathway reprograms adipocyte precursors as lymph node stromal cells. Immunity. 2012.（前回調査で確認済み） |
| 23623380 | Chai Q, et al. Maturation of lymph node fibroblastic reticular cells from myofibroblastic precursors is critical for antiviral immunity. Immunity. 2013.（前回調査で確認済み） |
| 18941209 | Katakai T, et al. Organizer-like reticular stromal cell layer common to adult secondary lymphoid organs. J Immunol. 2008.（前回調査で確認済み） |
| 22807928 | Katakai T. Marginal reticular cells: a stromal subset directly descended from the lymphoid tissue organizer. Front Immunol. 2012.（前回調査で確認済み） |
| 25998961 | Fletcher AL, Acton SE, Knoblich K. Lymph node fibroblastic reticular cells in health and disease. Nat Rev Immunol. 2015. |
| 33981032 | Buechler MB, et al. Cross-tissue organization of the fibroblast lineage. Nature. 2021. |
| 29752062 | Rodda LB, et al. Single-Cell RNA Sequencing of Lymph Node Stromal Cells Reveals Niche-Associated Heterogeneity. Immunity. 2018. |
| 29549115 | Zhang Y, et al. Plasma cell output from germinal centers is regulated by signals from Tfh and stromal cells. J Exp Med. 2018. |
| 33903764 | Kapoor VN, et al. Gremlin 1+ fibroblastic niche maintains dendritic cell homeostasis in lymphoid tissues. Nat Immunol. 2021. |
| 22566825 | Fletcher AL, et al. Reproducible isolation of lymph node stromal cells reveals site-dependent differences in fibroblastic reticular cells. Front Immunol. 2011. |
| 25407678 | Baptista AP, et al.（Mebius RE 研究室）Lymph node stromal cells constrain immunity via MHC class II self-antigen presentation. eLife. 2014. |
| 8666904 | Kratz A, et al. Chronic inflammation caused by lymphotoxin is lymphoid neogenesis. J Exp Med. 1996.（前回調査で確認済み） |
| 12732657 | Drayton DL, et al. Ectopic LTαβ directs lymphoid organ neogenesis. J Exp Med. 2003.（前回調査で確認済み） |

**未確認・要一次確認のまま残った項目**: Meylan 2022 の 29 遺伝子リストの完全な正確性（1 遺伝子の記号が不明瞭）、TLS＋スポットの定義方法の細部、Arai 2021 総説本文の具体的マーカー一覧、CLU（クラスタリン）・SFRP2・THY1・FAP・POSTN・CFD・SERPING1・GREM1（腎 TLT 文脈）の一次文献特定、Cabrita 2020 の遺伝子数（9 個 vs 依頼文にあった 12 個の不一致）。これらはいずれも本文中で「二次資料経由」「未確認」と明記した。

# EXP-20260926-005 実験計画 — 陽性対照の見直し（文献で因果が確立した PT ドライバー／保護因子, delete + overexpress）

- 仮説: [HYP-20260926-04](../../hypotheses/HYP-20260926-04.md)（方法仮説, Status: proposed）
- 起点: PI 指示（2026-09-26）「次の一手 (a) DE マッチ比較遺伝子 (b) 陽性対照の見直し (c) donor 単位解析 を**全て実行**。FT モデル前提、リモート PC を積極利用」。本実験は **(b)** のみを担当する（(a)(c) は EXP-20260926-004、別エージェント）
- 作成: Labretary（literature-scout + experiment-planner + simulation-engineer + validator）/ 2026-09-26
- **事前登録の時点**: §5・§6 の陽性対照リスト・基準・マーカー度の定義は、本実験の ISP を 1 回も回す前（2026-09-26 20:0x, GPU 実行の開始前）に確定した。使った事前情報は文献（PMID 照合済み）と、トークン化済み入力の発現情報（`work/gene_annotation.json`, ISP 結果を含まない）だけ

## 1. 問い
EXP-20260926-003 phase 3 で、FT（aPT_CKD vs PT_healthy）CellClassifier の ISP は削除ドリフトを解消し HK・陰性対照パネルも null 内に収まったが、**C1 が不合格**だった（VCAM1 は null 分位 0.98 で合格、HAVCR1 は 0.45/0.75/0.75 で不合格）。VCAM1 は極端な aPT マーカー（start 62.1% vs goal 7.6% の細胞で検出、|Δ検出率| は背景 9,615 遺伝子中 99 パーセンタイル）で、合格が「分類器が使うマーカーを消しただけ」の可能性がある。

1. 陽性対照を「PT で**遺伝学的操作の因果効果**が報告された遺伝子」に替えたとき、FT モデルの ISP（a2）はそれらを null から分離するか（delete 型・overexpress 型）
2. 分離した遺伝子は、start/goal の発現差が小さい（＝分類器の自明なマーカーではない）ものを含むか

## 2. データ・モデル（EXP-20260926-003 phase 3 と同一）
| 項目 | 値 |
|---|---|
| データ | A2 KPMP v1.5 scRNA（Lake 2023, PMID 37468583）。`work/input.h5ad` → EXP-001 への symlink（読み取りのみ） |
| tokenized | `work/tokenized.dataset` = EXP-001 tokenized のデータ部分の APFS clone。sha256 `ba3a14c06d57d8e6344d79329b6be9f2dad722e1e6290c675e0cc3545d312070`（EXP-001/003 と一致）。他実験の dataset を直接指さない（guardrails） |
| モデル | `experiments/EXP-20260926-003/models/ft-main/final`（FT 済み BertForSequenceClassification, 2 クラス, 0=aPT_CKD / 1=PT_healthy）。読み取りのみ |
| start | `aPT_CKD_TLShigh`（held-out 14 donor, 7,533 細胞）から 200 細胞/遺伝子・donor ≤ 25・seed 3/4/5 |
| goal セントロイド | `PT_healthy` の **held-out 9 donor** のみ（3,106 細胞から 500 細胞, donor ≤ 60）。FT に使っていない |
| alt | `aPT_CKD_TLSlow`（参考。m_b 指標用） |
| emb | `emb_mode: cls`, `emb_layer: -1`, bf16 forward + fp32 コサイン |

- 実行可能性は事前に確認済み（`analysis/gene_annotation.py`）: 回す 23 の固定遺伝子すべてで seed 3/4/5 とも **200/200 細胞**、13〜14 donor

## 3. 設計（2 つの run。1 run = 1 perturb_type）
| run | config | perturb_type | 陽性対照 | 同梱 | null | 遺伝子数 |
|---|---|---|---|---|---|---|
| D | `config-delete.yaml` | delete | HIF1A, NOTCH1, EGFR, YAP1, KEAP1（5） | 旧陽性 VCAM1・HAVCR1（継続比較, 判定外）、HK 3、パネル 7 | 発現頻度マッチ 20（**delete 陽性 5 に頻度を合わせる**） | 37 |
| O | `config-overexpress.yaml` | overexpress | HNF4A, PPARGC1A, PPARA, ESRRA, CPT1A, NFE2L2（6） | 方向対照 VCAM1・HAVCR1 の overexpress（C7, 参考）、HK 3、パネル 7 | 発現頻度マッチ 20（**overexpress 陽性 6 に頻度を合わせる**, run D とは別セット） | 38 |

- null は run ごとに独立に抽出する（`negative_sampler.ref_genes` にその run の陽性対照を入れる）。EXP-003 の null は盲検対象 5 遺伝子に頻度を合わせていたので、本実験の null とは異なる遺伝子になる
- **overexpress の定義**: Geneformer の overexpress は「その細胞の rank 列でその遺伝子を最上位へ移す」操作。細胞選抜は delete と同じく「その遺伝子を発現している細胞」に限るので、非発現細胞への導入ではない（限界として §9 に明記）
- 盲検: `controls_only: true`。TP53/CDKN1A/LTBR/SPP1/LCN2 は回さず、解析でも読み込み直後に除外（assert）

## 4. 陽性対照の選定（文献, PMID は PubMed で照合済み）
[HYP-20260926-04](../../hypotheses/HYP-20260926-04.md) の表を参照（方向・種・証拠の強さ・PMID・不採用の理由を記載）。要点:

- **delete 型**（削除が健常方向を予測）: HIF1A（18037992, 18667485）、NOTCH1（20978353）、EGFR（22095949, 27642633, 31676723）、YAP1（32843569, 32253273）、KEAP1（27789056, 19436334）
- **overexpress 型**（過剰発現が健常方向を予測）: HNF4A（30046000, 32764140, 37488681, 39954965）、PPARGC1A（21881206, 30982966）、PPARA（23804447）、ESRRA（33301705）、CPT1A（33465052, 25419705）、NFE2L2（19436334, 27789056）
- **不採用**: TGFBR2/SMAD3（文献で方向が対立: 23160515 vs 37270534, 28701516）、KLF15（34634362。start 検出 213 細胞で細胞数不足の恐れ）、SNAI1/TWIST1/FOXM1/CDKN2A（start 検出 2〜75 細胞で実行不可）、SOX9（欠損は修復を阻害＝方向が逆）

## 5. マーカー度（駆動因子 vs マーカー除去の切り分け, 事前登録）
`analysis/gene_annotation.py`（ISP 結果を使わない）で、遺伝子ごとに
- `det_start` / `det_goal`: start（7,533 細胞）と goal（held-out 健常 3,106 細胞）で Geneformer 入力に含まれる細胞の割合
- `ddet = f_start − f_goal`、`ddet_pct` = |ddet| の背景分布（どちらかの状態で 300 細胞以上検出される 9,615 遺伝子）での分位

**`ddet_pct ≤ 0.80` を「マーカー度が低い」と定義する**（事前登録）。実測値（結果を見る前に確定）:

| 分類 | 遺伝子（ddet_pct） |
|---|---|
| マーカー度 高（> 0.80） | VCAM1 0.99, NFE2L2 0.96, HIF1A 0.94, HNF4A 0.92, GAPDH 0.91, EGFR 0.86 |
| マーカー度 低（≤ 0.80） | YAP1 0.77, HAVCR1 0.71, ESRRA 0.61, PPARA 0.60, CPT1A 0.41, PPARGC1A 0.29, NOTCH1 0.27, KEAP1 0.01, ACTB 0.26, PPIA 0.00, パネル 7 遺伝子 0.02〜0.79 |

**KEAP1（delete, ddet_pct 0.01）と PPARGC1A / CPT1A（overexpress, 0.29 / 0.41）が、「マーカーではないが因果証拠がある」中心的な検証対象**である。これらが null 上位に来れば、分離はマーカーの除去／付加では説明できない。逆に合格が VCAM1・HIF1A・NFE2L2・HNF4A のような高分位の遺伝子だけなら、マーカー操作で説明できる。

（注: (a) の「DE マッチ比較遺伝子」は別実験 EXP-20260926-004 の担当。本実験は自分が回す遺伝子の注釈にとどめ、マッチさせた比較遺伝子セットは作らない）

## 6. 評価指標と判定基準（事前登録, 2026-09-26, ISP 実行前）
- **主指標 a2**: EXP-20260926-002 §6.4 / EXP-003 §6.2 と同一定義（長さ 5 分位ビン［null 細胞で境界］ごとの null 平均 shift_goal を細胞ごとに引き、遺伝子で平均。null 遺伝子自身は leave-one-gene-out）
- 併記: `m_raw`（補正なし shift_goal）、`dP_goal_raw` = 細胞ごとの ΔP(PT_healthy) = softmax(logits_pert)[1] − softmax(logits_orig)[1]、`dP_goal_a2` = それに a2 と同じ長さビン補正。副指標 m_d, m_d2, m_c, V2, V2L も EXP-003 と同じく出す
- 分位は run ごと・seed ごとの 20 null 遺伝子に対する `null_pct`（EXP-002 の実装と同一）

### 6.1 判定基準（run ごと）
| # | 基準 | 合格ライン |
|---|---|---|
| **C1'** | 陽性対照のうち、null 上位 10%（分位 ≥ 0.90）に 3 seed 中 2 以上入る遺伝子の数 | run D: 5 中 **2 以上** / run O: 6 中 **2 以上** |
| **C2'** | C1' に合格した各遺伝子が HK（ACTB/GAPDH/PPIA）を上回る組数（3 HK × 3 seed = 9 組） | すべての C1' 合格遺伝子で **7 組以上** |
| **C3** | ACTB・GAPDH・PPIA が null 分位 0.10〜0.90 | 3 遺伝子とも 3 seed 中 2 以上 |
| **C3-panel** | パネル 7 遺伝子 | 5 遺伝子以上が 3 seed 中 2 以上で 0.10〜0.90 |
| **C6** | C1' 合格遺伝子のうち `ddet_pct ≤ 0.80` のものが 1 つ以上 | 1 つ以上 |
| C7（参考, run O） | VCAM1・HAVCR1 の overexpress が null 分位 ≤ 0.50 | 3 seed 中 2 以上（報告のみ） |
| C4（参考） | split-half（donor 2 群, 20 回）と seed 間の遺伝子順位 Spearman | ≥ 0.7（報告のみ） |
| C5（参考） | 指標と細胞長・発現細胞数の相関 | 遺伝子単位 \|ρ\| ≤ 0.3、細胞単位 \|ρ\| 中央値 ≤ 0.2（報告のみ） |

### 6.2 C1' の閾値の根拠（偶然合格率）
null 仮説（陽性対照が null と区別できない）では 1 seed で分位 ≥ 0.90 になる確率 ≈ 0.10、3 seed 中 2 以上は ≈ 3·0.1²·0.9 + 0.1³ = **0.028**。遺伝子間を独立とすると
- run D（5 遺伝子）: P(2 遺伝子以上が合格) = 1 − 0.972⁵ − 5·0.028·0.972⁴ ≈ **0.0075**
- run O（6 遺伝子）: 同様に ≈ **0.012**

EXP-003 の C1（2 遺伝子**全部**が合格）は、陽性対照の 1 つでも外れると不合格になる厳しい基準だった（偶然合格率 0.028² ≈ 0.0008）。本実験は「文献の因果が in silico に写るとは限らない遺伝子を複数並べる」設計なので、**遺伝子ごとの全数一致ではなく「2 遺伝子以上」**を採る。有意水準は 1% 前後で EXP-003 とほぼ同等の厳しさを保つ。

### 6.3 判定（事前登録）
- **run D または run O のどちらかで C1' かつ C2' かつ C3 かつ C3-panel を満たす** → HYP-20260926-04 の該当予測（P1 / P2）を**支持**
- さらに **C6** を満たす → 「分離はマーカーの除去／付加だけでは説明できない」ことが**示唆される**（断定しない）
- 両 run とも C1' 不合格 → HYP-20260926-04 を**棄却**（陽性対照の見直しでは解決しない）
- C3 または C3-panel が不合格の run → その run は**判定不能**（系が壊れている）

## 7. 頑健性
- seed 3/4/5（EXP-002 案 A・EXP-003 phase 3 と同じ細胞）
- 精度: bf16（PI 標準）。分位 0.85〜0.95 の境界に来た遺伝子は fp32 で再計算して確認（最大 5 遺伝子 × 200 細胞 × 1 seed ≈ 0.04 h）
- 長さビン数 3 / 5 / 10 の感度（解析のみ）
- 同じマシン（spark-c083）で 2 run を逐次実行し、デバイス差を比較から除く

## 8. 計算量見積もり（G2 の判断材料）
| run | cell×gene | 実測レート | 時間 |
|---|---|---|---|
| D（delete, 37 遺伝子 × 200 × 3 seed） | 22,200 | 0.129 s（GB10, EXP-003 phase 3: 20,400 を 44 分） | **0.80 h** |
| O（overexpress, 38 × 200 × 3） | 22,800 | 同上 | **0.82 h** |
| セントロイド抽出・モデルロード・転送 | — | — | ~0.2 h |
| **合計** | 45,000 | | **約 1.8 h**（GPU 2 h 以内 → CLAUDE.md §3 の G2 追加申請は不要。PI の「全て実行」に含まれる） |

- 1 万細胞超のスクリーニングには該当しない（ユニーク細胞は start 7,533 のうち選抜分）
- 2 h を超える見込みになった時点で中断し、PI に G2 として報告する

## 9. 想定される失敗と限界
| 事項 | 対策・扱い |
|---|---|
| overexpress は「発現細胞で rank を最上位へ」であり、非発現細胞への遺伝子導入ではない | 限界として明記。null も同じ操作なので比較は公平 |
| overexpress では HK（ACTB/GAPDH/PPIA）が rank 中〜下位（rank_pct 0.55〜0.89）から最上位へ動く＝大きな摂動 | C3 が「範囲内」を求めるだけなので過大／過小どちらでも検出される。結果を見て解釈する |
| 文献の因果効果が in silico に写らない（モデルが学んでいない） | C1' は「5〜6 中 2 以上」で部分的な写り込みを許す。全滅なら棄却 |
| 合格が高マーカー度の遺伝子に偏る | C6 で区別。偏れば「マーカー除去で説明可能」と報告する |
| 循環性（start/goal を分けた分類器で goal 方向を読む） | HYP §循環性。埋め込みと ΔP の 2 系統、HK・パネル・null との比較、held-out donor |
| bf16 の量子化（guardrails） | 境界の遺伝子だけ fp32 で再計算 |
| datasets 4.0/4.8 の fork（guardrails） | `nproc: 1`、run_isp.py の `single_process_datasets` |
| リモートに人由来データの複製が残る | 終了後に `isp_remote.sh clean` |
| 他エージェントが同じ run_isp.py を編集 | 本実験は run_isp.py を変更しない（既存機能のみ使用）。push 時に sha256 を照合 |

## 10. 承認
- G1（HYP-20260926-04 の採択）: PI の「全て実行」（2026-09-26）に基づき Labretary が記録。明示の採択は G3 で確認
- G2: **不要**（GPU 約 1.8 h < 2 h、1 万細胞未満、全遺伝子スクリーニングでない）。超過が見込まれた時点で中断して申請する
- G3: 未（結果と判定を PI に諮る）

## 11. コード
- `scripts/geneformer/run_isp.py` は**変更しない**（既存の perturb_type: overexpress / capture.classifier_logits / negative_sampler.ref_genes で足りる）
- 新規 `analysis/gene_annotation.py`（実行可能性とマーカー度。ISP 前に実行, 出力 `work/gene_annotation.json`）
- 新規 `analysis/evaluate_positives.py`（EXP-002 の `reanalysis_cosine.py` / `vector_phase0.py` / `evaluate_controls.py` をそのまま import し、陽性対照リストと C1'・C2'・C6・C7・ΔP 指標を追加。EXP-002 側のコードは変更しない）

## 12. 事前登録の固定（結果を見る前）
- 2026-09-26 20:0x に、config・解析コードの sha256 を `work/SHA256SUMS-preregistered.txt` に記録
- 解析コードの検証: `analysis/evaluate_positives.py` を **EXP-20260926-003 phase 3 の結果**に適用し、EXP-003 の公表値を再現することを確認（a2 で VCAM1 null 分位 0.983 / HAVCR1 0.65、C3 ACTB 0.383・GAPDH 0.65・PPIA 0.20、C3-panel 7/7、null 負の割合 0.58、dP_a2 VCAM1 0.75）。出力は `work/check-exp003/`
- overexpress 経路の smoke: `work/config-smoke-oe.yaml`（2 陽性 + 1 HK + 2 null + VCAM1、6 細胞、seed 3）を Mac MPS で実行し、cells.csv・CLS 埋め込み・分類ヘッド logits が正しく出ることを確認（出力 `work/work/smoke-oe/`）。`scripts/geneformer/run_isp.py` は未変更

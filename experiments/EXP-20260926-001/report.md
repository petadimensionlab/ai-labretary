# 実験レポート: EXP-20260926-001 — CKD 腎の修復不全型近位尿細管（aPT）における p53–p21 軸・LTβR 削除の in silico perturbation

- Status: final
- 仮説: [HYP-20260925-01](../../hypotheses/HYP-20260925-01.md) / 検証判定: **判定不能（inconclusive）**（テンプレートの Strong/Moderate/Weak/Refuted のいずれにも当てはめない）
- G3 承認: PetadimensionLab (PI), 2026-09-26 — inconclusive
- 作成: Labretary / 2026-09-26

## 要旨（5 行以内）
ヒト CKD 腎（KPMP scRNA）の aPT で TP53・CDKN1A・LTBR を Geneformer 上で delete しても、健常 PT 方向への移動は事前登録した 2 指標のいずれでも支持されなかった（判定条件を満たした seed は 3 遺伝子とも 0/3）。ただし陽性対照がハウスキーピング遺伝子と区別できず、パイプラインの検出力が示されなかったため、結論は「判定不能」とした。主な原因は、どの遺伝子を削除しても埋め込みが全状態から一様に離れる「削除ドリフト」である。TLS-high と TLS-low の対比も、状態の分離が弱く意味をなさなかった。次の実験では削除ドリフトを差し引いた指標を開発する。

## 1. 背景と仮説
腎の三次リンパ組織（TLS/TLT）は加齢・障害腎で形成され、修復不全の尿細管と相互作用すると報告されている（Sato 2016, PMID 27699223; Yoshikawa 2023, PMID 37548710）。本実験は、細胞老化（p53–p21 軸）と TLS 由来シグナル（LTα1β2 の受容体 LTβR）が修復不全型 PT（aPT）の状態維持に関わるという仮説を in silico で検証した。

仮説: CKD 腎の aPT において TP53・CDKN1A・LTBR を delete すると、細胞状態が aPT_CKD_TLShigh から PT_healthy へ近づき、その効果は TLS-low の CKD 腎の aPT より大きい（新規予測）。

## 2. 方法
### 2.1 データ
- A2 KPMP v1.5 scRNA（Lake BB et al., *Nature* 2023, PMID 37468583, DOI 10.1038/s41586-023-05769-3）。CELLxGENE から取得（`data/raw/kidney_tls/SOURCES.md`）。CELLxGENE の方針により CC BY 4.0（データセット個別の表示は未確認）。KPMP の利用規約に従い謝辞を記載する
- raw.X の整数カウントを使用し、`prepare_h5ad.py` で 17,684 細胞 × 30,635 遺伝子の入力を作成
- 状態（**TLS 群は推定ラベル**。donor 単位 pseudo-bulk の TLS 8 遺伝子スコアの三分位。組織学的な TLS ではない）

| 状態 | 定義 | 細胞数 / donor 数 |
|---|---|---|
| start: aPT_CKD_TLShigh | aPT ∧ CKD ∧ TLS high | 7,533 / 14 |
| goal: PT_healthy | (PT-S1/S2 ∨ PT-S3) ∧ normal | 6,149 / 22 |
| alt: aPT_CKD_TLSlow | aPT ∧ CKD ∧ TLS low | 4,002 / 8 |

### 2.2 モデルと perturbation 設定（config.yaml の要点）
- Geneformer V2-316M（Pretrained）、`emb_mode: cls`、最終層、Apple M4 Pro の MPS
- 精度: bf16（forward は bf16、コサイン比較は fp32）。**PI 決定**。phase 0 の bf16 と fp32 の比較では Spearman 0.879 で、../Geneformer の基準 0.98 に届いていない
- `perturb_type: delete`、1 遺伝子ずつ。遺伝子ごと・seed ごとに、その遺伝子を発現している start 細胞から donor あたり ≤ 15、計 ≤ 100 細胞を抽出
- run 1: start = TLS-high（alt = TLS-low）、run 2: start = TLS-low（文脈依存性の比較用）。seed 0/1/2
- 対象 5 遺伝子: TP53, CDKN1A, LTBR（新規予測）、SPP1, LCN2（参照）

### 2.3 コントロールと統計
- 陽性対照: VCAM1（Kirita 2020, PMID 32571916）、HAVCR1（Lake 2023）
- 陰性対照: 発現頻度をマッチさせたランダム 20 遺伝子（seed ごと）、ハウスキーピング ACTB, GAPDH, PPIA
- 指標（事前登録）: `shift_goal`（goal セントロイドへのコサイン類似度の増加）と `shift_goal_minus_start`（goal 方向 − start 方向。phase 0 を見て追加したが、本番データを見る前の G2 で登録）
- 判定: null 上位 5% かつ FDR < 0.05（null の正規近似 p に対象 5 遺伝子内で BH）を 3 seed 中 2 以上で満たせば支持。陽性対照は null 上位 10% に 3 seed 中 2 以上。文脈依存性は shift(TLS-high) / shift(TLS-low) ≥ 1.2 かつ donor ブートストラップ 95% CI の下限 > 1（両 shift が正であることを必須とする）

## 3. 結果
### 3.1 コントロールの挙動
- 陽性対照: `shift_goal` では HAVCR1 が上位 10% に 1/3 seed のみ（fp32 再計算では VCAM1 も 1/3）→ 不合格。`shift_goal_minus_start` では VCAM1・HAVCR1 とも 3/3 seed で合格
- **ハウスキーピング遺伝子も同じように上位に入った**: GAPDH・PPIA は `shift_goal_minus_start` で 3/3 seed とも null 分位 1.00、`shift_goal` でも 2〜3/3 seed で上位 10%。陽性対照の合格は特異的な信号を示していない
- 削除ドリフト: run 1 の遺伝子 × seed 90/90 が goal・start・alt の 3 セントロイドすべてから離れた（`figures/deletion_drift.png`）

### 3.2 主要結果（図表）
run 1（tls_high）、各遺伝子 100 細胞 × 3 seed。平均 ± SD は seed 間。

| 遺伝子 | shift_goal 平均 ± SD | null 分位（seed 0/1/2） | shift_goal_minus_start 平均 ± SD | null 分位 | 判定条件を満たした seed |
|---|---|---|---|---|---|
| TP53 | −2.45×10⁻⁴ ± 3.3×10⁻⁵ | 0.55 / 0.45 / 0.35 | −2.65×10⁻⁵ ± 1.3×10⁻⁵ | 0.85 / 0.85 / 0.60 | 0/3（両指標） |
| CDKN1A | −2.07×10⁻⁴ ± 4.1×10⁻⁵ | 0.55 / 0.50 / 0.95 | −3.98×10⁻⁵ ± 9.7×10⁻⁶ | 0.50 / 0.30 / 0.80 | 0/3（両指標） |
| LTBR | −3.09×10⁻⁴ ± 1.0×10⁻⁴ | 0.10 / 0.50 / 0.35 | −4.71×10⁻⁵ ± 1.0×10⁻⁵ | 0.45 / 0.30 / 0.60 | 0/3（両指標） |

- FDR（p_norm, BH）の最小は 0.23（CDKN1A, shift_goal, seed 2）
- 事前登録規則の機械的な適用: `shift_goal` は陽性対照不合格のため inconclusive、`shift_goal_minus_start` は refuted。ただし後者の陽性対照の合格はハウスキーピング遺伝子と区別できないため、総合判定は inconclusive（G3）
- 文脈依存性: 基準を満たした遺伝子なし。LTBR の比 1.96（95% CI 1.02–4.27）は分子・分母ともに負（−3.25×10⁻⁴ / −1.65×10⁻⁴）で、TLS-high の aPT の方が健常から大きく離れることを意味し、仮説と逆向き
- 図: `figures/null_vs_genes_tls_high.png`, `figures/null_vs_genes_tls_low.png`, `figures/fp32_recheck.png`

### 3.3 頑健性
- seed 間の Spearman（run 1, shift_goal）: 0.745 / 0.718 / 0.673（最小 0.673 < 0.7）。差分指標では 0.89–0.94。run 2 は最小 0.736
- donor を 1 人ずつ除いた解析: 主要予測が成立する seed は最大 0 のままで結論は変わらない。CDKN1A の分位による棄却のみ donor の抜き方に依存
- fp32 再計算（事前登録、run 1 の境界にある 10 組、32.5 分）: 主要予測の判定は変わらなかった。変化は `shift_goal` の VCAM1（合格 → 不合格）と差分指標の SPP1（1/3 → 3/3）。bf16 と fp32 の差は `shift_goal` で mean|Δ| 1.7×10⁻⁵（|shift| の約 15%）、差分指標で 8.2×10⁻⁶（約 59%）、HAVCR1 seed 0 で差分指標の符号が反転

## 4. 検証（validation.md の要約）
| 項目 | 判定 |
|---|---|
| 陽性対照（shift_goal / 差分指標） | Fail / Pass（非特異的） |
| 陰性対照（HK が null 内に収まるか） | Fail |
| seed 間の頑健性（run 1） | Fail（最小 0.673） |
| 交絡: 細胞の長さ | Fail（遺伝子単位の Spearman −0.5〜−0.6） |
| 交絡: 発現頻度 | run 1 Pass（0.20）/ run 2 Fail（0.32） |
| donor 依存性 | 結論は不変 |
| 再現性 | Pass（config hash 一致、実行後の編集なし） |

詳細は [validation.md](validation.md)。

## 5. 考察
- この設計では、p53–p21 軸と LTβR の削除が aPT を健常 PT 方向へ特異的に動かすという予測は**支持されなかった**。ただし、陽性対照がハウスキーピング遺伝子と区別できなかったため、「これらの遺伝子に効果がない」ことが**示されたわけではない**
- 主因は削除ドリフトと考えられる。1 遺伝子の delete は入力トークン列を変え、CLS 埋め込みを全状態から一様に遠ざける。この一様成分が遺伝子特異的な方向の変化（|shift| ~10⁻⁵）より大きく、null との比較が「goal に近づくか」ではなく「ドリフトが小さいか」を測ってしまう。細胞の長さ（トークン数）と shift の相関もこの解釈と整合する
- TLS-high と TLS-low のセントロイド間距離（1 − cos = 0.0006）は start–goal（0.0073）の約 1/12 で、細胞ごとの start 方向と alt 方向のシフトの相関は 0.998。事前チェックの AUC 0.627 ± 0.206 と合わせ、推定 TLS 群による aPT の対比は埋め込み上ほとんど区別できない。文脈依存性の検証には、分離がはっきりした状態定義が必要と示唆される

## 6. 限界
- in silico perturbation はモデルが学習した発現共変動に基づく予測であり、因果関係を直接示すものではない
- TLS 群は pseudo-bulk の組成から推定したラベルで、組織学的な TLS ではない
- 精度は PI 決定により bf16。差分指標では判定閾値との差（~1×10⁻⁵）が bf16 の誤差と同じ桁。fp32 再計算は境界の遺伝子のみで、null は bf16 のまま比較した
- 差分指標は phase 0 の結果を見てから追加した（本番データを見る前に G2 で登録）
- run 2 の CDKN1A（69 細胞・7 donor）、HAVCR1（62・7）、LCN2（48・6）は細胞数が異なり、同列比較できない
- 多重比較は対象 5 遺伝子内のみで補正。seed 間の頑健性の比較は約 10 遺伝子に限られる
- 仮説設計からの逸脱: dtype を bf16 に固定（PI 決定）、`nproc` を 8 → 1（MPS 使用後の fork で停止するため）、FDR を経験的 p ではなく正規近似 p で計算（null 20 遺伝子では経験的 p の下限が 0.048 のため、G2 で決定）

## 7. 次のアクション（in silico / wet）
- **in silico（PI が指定した次の実験, EXP-20260926-002）**: 削除による一様な効果を差し引いた指標を開発する。既存の細胞ごとのシフトによる補正（null 平均の差し引き、細胞長ごとの標準化など）と、摂動前後の埋め込みベクトルから共通のドリフト方向を除いて goal 方向に射影する方法を、コントロールだけで評価して選ぶ（主要遺伝子は指標を決めるまで見ない）。結果はまだない
- in silico（その後）: 細胞長をマッチさせた null、HK を陰性対照の合否判定に含める、TLS 状態を空間データ（Meylan 2022, PMID 35231421）由来のシグネチャーで定義し直す
- wet（in silico で支持が得られた場合に検討）: ヒト近位尿細管細胞で修復不全・老化様の状態を誘導し、TP53/CDKN1A/LTBR の siRNA で HAVCR1・VCAM1・SA-β-gal の変化を見る（HYP-20260925-01 §「支持された場合の次の一手」）

## 補足: 全 run の一覧（失敗 run を含む）
| run | 内容 | 結果 |
|---|---|---|
| phase 0 a | 入力作成・dry-run | 成功 |
| phase 0 b | 事前チェック（CLS 埋め込みの AUC） | 成功: start vs goal 1.00 ± 0.00、start vs alt 0.627 ± 0.206 |
| phase 0 c（初回） | bf16 精度ペア（nproc 8） | **失敗**: 1 遺伝子目で datasets.map の子プロセスが停止（guardrails に記録。ログは残っていない） |
| phase 0 c | bf16 と fp32 の精度ペア（nproc 1） | 成功: Spearman 0.879 |
| 本番（初回） | run 1 | **失敗**: 実験ディレクトリに `logs/` がなく起動直後に停止（計算は行われていない） |
| 本番 run 1 | tls_high, 3 seed | 成功（1 h 55 min） |
| 本番 run 2 | tls_low, 3 seed | 成功（1 h 39 min） |
| fp32 再計算 | run 1 の境界 10 組 | 成功（32.5 min） |

## 再現情報
- config: `config.yaml`（sha256: c5e72b5ff4e4e89b4e609e40518d775a9595810be5e82cc48990dcd6f7a58863。env.json・results・現行 config で一致）
- env: `results/env.json`（2026-09-26 01:53:39 記録。macOS 27.0 arm64、Apple M4 Pro、Python 3.12.13、torch 2.14.0、transformers 4.46.3、datasets 4.0.0、geneformer 0.1.0）
- git commit: なし（プロジェクトは git 管理外）
- 解析コード: `scripts/geneformer/prepare_h5ad.py`, `scripts/geneformer/run_isp.py`, `scripts/analysis/validate_isp.py`, `validation-analysis/*.py`

KPMP 謝辞: KPMP 所定の謝辞文を掲載すること（**文言は KPMP の利用規約で要確認。本レポートでは未確認**）。

# EXP-20260926-005 検証（validator）

対象: `results-delete/`（run D）と `results-oe/`（run O）、判定は `*/eval/criteria.json`（事前登録 `plan.md` §6）。

## 1. 事前登録の遵守
| 項目 | 確認 | 結果 |
|---|---|---|
| config・解析コードが事前登録時から未変更 | `shasum -a 256 -c work/SHA256SUMS-preregistered.txt` | 4 ファイルとも OK |
| 判定基準を結果を見た後に変えていない | `criteria.json` は config の `analysis.criteria` をそのまま読む。config は未変更 | OK |
| 盲検（TP53/CDKN1A/LTBR/SPP1/LCN2） | `controls_only: true` で未実行、`load_run` が読み込み直後に除外して assert | OK |
| 入力データの同一性 | tokenized の `data-*.arrow` sha256 = `ba3a14c0…`（EXP-001/003 と一致）、FT 重み `e1bc1fdf…` がローカル＝リモートで一致 | OK |
| 実行環境 | 2 run とも同一ホスト spark-c083（CUDA, bf16）。`device.json`・`env-remote.json` に記録 | OK |
| 事後に追加した解析 | `analysis/report_tables.py`（表の整形＋参考 C4/C5）、`analysis/explore_rank.py`（rank 交絡）、`dM_raw`（logit 差）。いずれも**事前登録の判定には使っていない**。README に「探索的（事前登録外）」と明示 | OK |
| 成果物の実在 | `ls` で確認（run D 795 ファイル、run O 816 ファイル。いずれもリモート側と一致） | OK |

## 2. 判定の再確認
- **run D（delete）**: C1' 0/5 で不合格。C3・C3-panel は合格（系は壊れていない）ので、これは**有効な否定的結果**。
- **run O（overexpress）**: C1' 3/6 は合格だが、**C3 が不合格**（GAPDH が 3 seed とも分位 0.95）。`plan.md` §6.3 の「C3 または C3-panel が不合格の run は判定不能」に従い、**run O は判定不能**。C1'・C6 が合格していることをもって P2・P3 を支持と書いてはいけない。
- 上記より `plan.md` §6.3 の「両 run とも C1' 不合格 → 棄却」にも「どちらかが C1'+C2'+C3+C3-panel → 支持」にも当たらない。**HYP-20260926-04 は inconclusive（判定不能）**が事前登録に忠実な結論。

## 3. 符号・向きの確認（guardrails: 比や分位を報告する前に向きを見る）
- run O の null は補正前 m_raw が正（負の割合 0.02）＝ overexpress では**全遺伝子が goal 方向へ一様に動く**。a2（null 平均を引く）はこの一様分を除いており、分位上位 = 「一様分より goal 側」。run D の null 負の割合は 0.12（m_raw）で、EXP-001 で問題になった削除ドリフトは FT モデルでは解消している。
- 方向対照（C7）は VCAM1・HAVCR1 とも 3 seed で分位 0.05〜0.10、a2 は負（−22.6e-4 / −16.4e-4）＝障害マーカーの過剰発現は goal から**遠ざかる**。向きとしては期待どおりだが、下記 §4 の rank 交絡でも同じ順序になるため、これ単独では特異性の証拠にならない。

## 4. 最大の弱点 — run O の指標は「摂動の大きさ」を測っている疑い
`eval/explore_rank.json`（探索的）:
- 元のトークン位置（pos/length。大きいほど元は下位）と a2 の遺伝子単位 Spearman は run O で 0.65 / 0.47 / 0.32（seed 3/4/5）、**null 20 遺伝子だけでも seed 3 で 0.51（p=0.02）**。run D では −0.11 / −0.06 / −0.23 と無関係。
- 固定 18 遺伝子では `rank_pct_start` と a2 の ρ = 0.71（p=0.001）。GAPDH（rank_pct 0.71）・EEF1D（0.72）が上位、元から高発現の PPIA（0.55）・ESRRA（0.24）・VCAM1（0.36）が下位。
- つまり run O の分位は「元の順位が低い遺伝子を最上位に上げると埋め込みが大きく動く」という**操作の大きさの効果**で相当部分が説明できる。C3 不合格（GAPDH）はその症状であり、独立した不具合ではない。
- **設計上の穴**: null は発現**頻度**でマッチさせたが、overexpress では効きしろ（元の rank）が効果量を決めるのに、そこはマッチしていない。次に overexpress で ISP をするときは **rank マッチ null**（または rank を共変量として補正した指標）が必要。

## 5. 副指標 ΔP の信頼性
- FT 分類器は start 細胞で飽和（P(PT_healthy) 中央値 0、99 パーセンタイル 1.2e-3）。さらに logits は bf16 格子（0.031 刻み）で、Δ(logit 差) が厳密に 0 の細胞が 53%。
- そのため dP 系は seed 間一致が低く（run D 0.27 / run O 0.055）、run D では細胞単位の長さ相関が −0.82 という見かけの構造も出る。**事前登録では副指標だったが、実質的に情報を持たない**。結論には使わず、限界として明記した。guardrails に追記済み。

## 6. 再現性・頑健性
- run D: C4 split-half 0.37（参考ライン 0.7 未達）、seed 間（固定 17 遺伝子）0.70。効果が null SD の中に収まる遺伝子ばかりなので、順位が安定しないのは整合的。
- run O: C4 split-half 0.871、seed 間 0.907 と高い。ただし §4 のとおり、安定して測れているのが生物学ではなく rank の可能性がある。
- bf16 の確認（plan §7）: run D は分位 0.85〜0.95 の陽性が 0 のため fp32 再計算は不要。run O は 9 組（陽性 × seed）が境界域に入ったが、**run O 自体が C3 不合格で判定不能**のため、fp32 再計算をしても結論は変わらない（rank 交絡を解かない限り分位の解釈が定まらない）。再計算は次の実験で rank マッチ null を入れてから行うのが妥当。
- 細胞数: run D は null 1 遺伝子（seed 4）が 179 細胞、他はすべて 200。run O は全 114 遺伝子×seed が 200。事前の実行可能性見積もり（plan §2）どおり。

## 7. データ・安全性
- `data/raw/` と EXP-001/003 の資産は読み取りのみ（tokenized は APFS clone を使用、他実験の dataset を直接指していない）。
- リモート（spark-c083）の複製 3.6 GB は `isp_remote.sh clean` で削除し、`ls` で消失を確認。資格情報は記録していない。

## 8. 検証の結論
- run D の否定的結果は**有効**（系は正常で、文献ドライバー 5 遺伝子は null から分離しなかった）。
- run O は事前登録どおり**判定不能**。C1' 合格は rank 交絡で説明できる可能性が高く、支持と書いてはならない。
- したがって **HYP-20260926-04 は inconclusive**。「陽性対照を因果の確立した遺伝子に替えれば解決する」とは言えず、同時に「棄却」も断定できない（run O が判定不能のため）。

## 9. 事後検証（2026-09-27 追加・探索的／事前登録外。§2 の判定と criteria.json は変更していない）
- `experiments/EXP-20260926-005/analysis/verify_rank_null.py` で §4 の穴への事後検証を実施。出力は `results-oe/eval/rank_null/` と `results-delete/eval/rank_null/`（7 variant: 補正指標 m_a2 / m_a2_rank（relpos 5 分位ビン）/ m_a2_rank_f（relpos 10 ビン）/ m_a2_dist（pos 5 ビン）× 分位 full / rank マッチ null（relpos 最近隣 k=10、自身除外））。
- **上位の維持**: rank マッチ null でも run O の C1' 合格 3 遺伝子（HNF4A・PPARA・NFE2L2、pct +0.87〜+0.93）は全 variant で維持された。ESRRA・PPARGC1A は下位側（+0.03〜+0.20）のまま。→ 上位側は rank のみでは説明しにくい（依然として探索的）。
- **GAPDH（C3）**: rank/rank_f マッチ null では範囲内 2/3 seed でこの指標でのみ C3 合格。dist 補正では GAPDH の matched pct が +0.90〜+1.00 と下位側に寄る向きも観測され、variant 依存的。
- **交絡の除去度**: gene-level Spearman(relpos, 指標) は 5/10/pos ビンいずれでも大きくは縮まない（seed4 +0.47→+0.45〜+0.50、seed3 +0.65→+0.49〜+0.62、seed5 +0.32→+0.27〜+0.31）。＝ビン補正では残り交絡を除去しきれない（ビン内の相関が残る）。
- **run D**: 7 variant すべてで C1' 0/5（pct 上限 +0.63）。§3 の有効な否定的結果は rank 操作に対して頑健。
- **次にすべきこと（提案・PI 判断待ち）**: run O 系を継続するなら、ビン内パーセントイル検定または移動距離（pos）マッチ permutation を**事前登録**する。事前登録の run D/run O 判定（§2）と HYP status（inconclusive）はこの節でも変更しない。

# EXP-20260927-001 作業経過（experiment-planner）
（時刻は記録していない。項目は作業の順。事前登録の固定時刻は work/SHA256SUMS-preregistered.txt）
- 2026-09-27 12:46 ディレクトリ作成。templates/experiment をコピー。tokenized のデータ部分（data-*.arrow, dataset_info.json, state.json）を EXP-005 から `cp -c`（APFS clone）。sha256 ba3a14c0…（EXP-001 と一致）。input.h5ad は EXP-001 への symlink（読み取りのみ）
- 次: 事前確認（ISP 結果を使わず tokenized だけで、start 細胞での各遺伝子の元の位置の分布）→ null プールの選定 → plan.md / config.yaml
- precheck_pos.py 実行（tokenized のみ）→ work/pos_by_gene.parquet（start 7,533 細胞, 16,547 遺伝子）。遺伝子内の元の位置のばらつきが大きい（例 HNF4A q10–q90 = 202–2860）→ 細胞単位でマッチする並べ替え検定を主検定に選ぶ方針
- select_null_pool.py 初版（5×3 層の層別無作為 60 遺伝子）。細胞数は全 78 遺伝子×seed 11/12/13 で 100。マッチ許容 0.5 で使える null 数: 陽性 37〜54、HK 14〜32（PPIA が 14〜22 と不足）、VCAM1 2〜6 → 層別無作為では短い細胞・高頻度遺伝子の被覆が不足
- Labretary から追加指示: 研究報告 §7 推奨 1〜3（null を頻度＋元の順位でそろえる／陽性対照ごとに頻度・順位・発現差をそろえた比較遺伝子群を事前に決め、上回ることを合格条件に加える／提供者単位の副指標と固定遺伝子集合の折半信頼性）
- 計算時間の実測: EXP-005 の遺伝子ごと時間 = 約 3 s（固定）+ 4.5 s / 10 万トークン（r=0.53）。100 細胞なら 1 遺伝子×seed 約 14〜15 s
- 次: 比較遺伝子（陽性 1 つにつき 5）→ null プール（固定 18 遺伝子それぞれに頻度・元の位置・長さの近い 3 遺伝子 = 54）を選び直す
- select_genes.py 確定（select_null_pool.py は廃止・削除）。比較遺伝子 30（陽性 6 × 5, work/comparators.csv, 候補と除外理由は work/comparator_candidates.csv）、null 54（固定 18 × 3, work/null_pool.json）。細胞数は全 102 遺伝子 × seed 11/12/13 で 100。使える null 数（許容 0.5）: 固定遺伝子は VCAM1（7〜15, 参考のみ）以外 22 以上、比較遺伝子 27 以上。null のうち AQP1・TMEM176A・RNASEK は一部 seed で 20 未満（較正（calibration。計測分野の「校正」(JIS Z 8103) とは別語）から除外される）
- 見積もり: 306 遺伝子×seed、7.41e7 トークン → ISP 1.19 h（固定 3.18 s + 4.49 s/1e5 トークン）＋セントロイド・起動・転送 ≈ 0.2 h → 約 1.4 h（上振れ 1.6 h）
- 次: config.yaml に遺伝子リストを反映 → HYP-20260927-01 → plan.md → sha256 固定
- evaluate_matched.py（評価・判定コード）を作成。合成データ検査 synthetic_check.py（work/synthetic/）で主検定を 2 回修正:
  (1) 許容幅内の単純マッチ → 分布の端で中央寄りに引かれる偏りで null の残差相関 0.16〜0.45 → 2 次式の偏り補正を追加
  (2) 1 つの null 遺伝子から細胞を重複ありで引く疑似遺伝子 → null 間のばらつきが過大で保守的（上位 10% の割合 0.02）→ 遺伝子ラベルの並べ替え（重なる null 遺伝子の中での順位）に変更
  最終版: 交絡のみ → 合格 0（否定）、効果 2 遺伝子 → 合格 2（支持）、較正 0.104 / 0.097、残差 |ρ| ≤ 0.06。長さだけの補正（従来 a2）は同じデータで残差 ρ 0.81〜0.85
- 次: HYP-20260927-01, tracker, plan.md, sha256 固定
- 事前登録を固定（work/SHA256SUMS-preregistered.txt）。ISP は未実行

## 実行（simulation-engineer, 2026-09-27）
- PI 承認済み方針（2026-09-27）を受けて実行開始。承認済み見積もり 1.4 h（上振れ 1.6 h）、G2 の目安（2 h）未満
- 13:08 事前登録の固定を確認: `shasum -a 256 -c work/SHA256SUMS-preregistered.txt` 全 OK（11 ファイル）
- 13:10 `scripts/capture_env.sh` で `results/env.json` 保存（Mac, config_sha256 が事前登録と一致）
- 13:10 手元で `run_isp.py --dry-run` 実行 → OK（内蔵見積もり 0.58 h は既知の過小評価, guardrails: negative_sampler.method=explicit の 54 遺伝子を数えない）
- 13:11 リモート計算機 2 号（spark-c083, coi-next1@100.94.57.17, REMOTE_GF_HOME=~/workspace/research/Geneformer）へ push（サンドボックス外）
  - run_isp.py・capture_env.sh・config.yaml は cat/mv 方式で送り sha256 一致確認
  - tokenized.dataset（96M）・ft-main/final モデル（1.2G, models/Geneformer/* ではないため rsync 本体転送）・geneformer_pkg を openrsync -v で送信
  - 転送後にリモートで ls・sha256 確認: tokenized の arrow ファイル sha256 = ba3a14c0…（Mac・事前登録と一致）、model.safetensors sha256 = e1bc1fdf…（Mac と一致）
- 13:12 `isp_remote.sh run` でリモート起動（tmux セッション isp-config_）
- 13:37 seed 11（102 遺伝子）完了、経過 23 分（0.38 h）。中断規則の基準（0.6 h）未満のため続行（seed 12/13 は同一プロセスが自動継続）
- 14:26:40 全 seed（11/12/13）完了。exit=0。ログにエラーなし（WARN は compute.device=mps→cuda への想定内の切替のみ、他はライブラリの非致命的な警告）
- 実測時間: 起動 13:12:38 → 終了 14:26:40 = **74 分（1.23 h）**。事前登録の見積もり 1.4 h（上振れ 1.6 h）を下回った
- 回収（pull）: results/{run-11,12,13}-L-1, cells.csv, summary.csv, timing.jsonl, device.json, env-remote.json 一式。ls で各 seed の遺伝子ディレクトリ数（102/102/102）・cells.csv 数（102/102/102）を確認。config.sha256（回収結果側）= 事前登録 config.yaml の sha256 と一致（8c8e7ef5…）
- `isp_remote.sh clean` 実行 → リモート ~/labretary-remote（3.3G, ヒト由来データの複製含む）を削除。ssh で test -d による消失確認 = CONFIRMED_REMOVED
- 評価: `.venv/bin/python analysis/evaluate_matched.py --config config.yaml` 実行。出力 4 ファイル（criteria.json, gene_seed.csv, donor_level.csv, split_half.json）を ls で確認。MixedLM の ConvergenceWarning/SingularMatrixWarning は donor 数が少ない遺伝子×seed で発生する既知の非致命的警告（計算は完了）
- 主検定の結果: 陽性対照 6 遺伝子とも q≥0.90 が 2/3 seed 以上に達せず（(a) 不合格）。較正（C8a）・交絡の除去（C8b）・恒常発現遺伝子の中立性（C3）・中立対照パネルの中立性（C3_panel）はすべて合格。evaluate_matched.py の verdict_proposal = 「否定（移動距離をそろえた null から陽性対照が分離しない）」（提案。確定は validator → PI, G3）
- 完了。simulation-engineer としての作業終了

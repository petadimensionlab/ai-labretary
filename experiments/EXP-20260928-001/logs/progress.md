# EXP-20260928-001 実行ログ（Labretary/simulation-engineer）

- 04:02 事前確認: `shasum -a 256 -c work/SHA256SUMS-preregistered.txt` 24/24 OK（自分で実施）。
  `run_isp.py --config config.yaml --dry-run` OK（見積もり gpu_hours_estimate 表示は 0.17 とやや低いが
  measured_hours_total=3.38 とラベルされた値が出ており、config.compute.gpu_hours_estimate=6.0 と
  実測ベースの見積もりのどちらでも上限 24h 未満。needs_G2_approval=false）。
- 04:02 `scripts/capture_env.sh` で `results/env.json` 保存（host=Mac, config_sha256 一致確認）。
- 04:03 push（リモート2号 coi-next1@100.94.57.17, REMOTE_GF_HOME=~/workspace/research/Geneformer）。
  run_isp.py / capture_env.sh / config.yaml / tokenized.dataset / モデル重み(既存と同一) を転送、
  sha256 で config.yaml と tokenized.dataset の3ファイルをリモートで再計算し一致確認。
- 04:03:45 run 開始 → **即時失敗**（exit=1）。
  原因: `state_embs.rows_file` が参照する `work/state_rows-seed{41,42,43}.json` を push スクリプトが
  転送していなかった（push は tokenized_dataset と model.directory のみを rsync する設計で、
  config 中の `rows_file` 参照は対象外だった）。FileNotFoundError で seed41 の state_embs 抽出前に停止。
  GPU 消費はほぼゼロ（数秒で失敗）。
- 04:05 state_rows-seed{41,42,43}.json を cat/mv で個別送付、sha256 一致確認（3/3 OK）。
- 04:05:28 再実行開始。
- 04:07 確認: RUNNING。ログは seed41 (`=== run-41-L-1 ===`) を実行中、device=cuda/dtype=bfloat16 で問題なし。
  （注: リモートのログファイル実名は `remote-config_.log`（アンダースコア付き）。最初に設定した監視の
  grepパターンが `remote-config.log` を見ていて空振りしていたため、`isp_remote.sh status` 経由の
  監視に切り替えた。）
- 監視は Monitor（60秒間隔ポーリング, isp_remote.sh status 呼び出し）で継続。中断規則:
  seed41完了時点で開始(04:05:28)から3.0h超（07:05:28以降）なら停止・報告。全体9.0h超見込み
  （13:05:28以降に終わりそうなら）も同様。

## 完了（07:05:31, exit=0）
- 07:xx 独立に `isp_remote.sh status` を確認し NOT_RUNNING / exit=0 / 完了ログを直接確認（Labretary/PIからの連絡と一致）。
  実測 04:05:28〜07:05:31 = **3.00 h**（見積もり6.0hを大きく下回った）。中断規則（seed41完了3.0h超／全体9h超）には抵触せず。
- `isp_remote.sh pull` でresults一式を回収。
  - summary.csv 469行（ヘッダ+468 = 156遺伝子×3seed）
  - run-{41,42,43}-L-1/genes 各156ディレクトリ、cells.csv/stats/emb/pert.npy が全て揃っていることを確認（欠落0）
  - timing.jsonl の gene エントリも各seed156件（skip 0件）
  - negatives.json 各seed96件（negative_random_n=96と一致）
  - blocked_genes（LTBR/RELB/NFKB2/CDKN1A/TP53/SPP1/LCN2/CCL19/CCL21/CXCL13）が genes/ に1件も無いことを確認
  - results/config.sha256 とローカル config.yaml の sha256 が一致（実行に使われたconfigが事前登録と同一であることを確認）
- `isp_remote.sh clean` 実行 → `ssh ... ls ~/labretary-remote` がエラー（No such file）となり、複製の削除を確認。
- 評価には `work/extra_centroids-seed{41,42,43}.npz` が必要だが未生成だったため、
  `analysis/extra_centroids.py --config config.yaml --seed {41,42,43}` を Mac(MPS)でそれぞれ実行
  （seed41: 10:05完了/約4:47, seed42: 10:10完了, seed43: 10:15完了。摂動は行わないセントロイド抽出のみ、
  各seed MYOF_neg 565〜566細胞・FIB_org_unmatched 665細胞、start セントロイドとの重なり0件を確認）。
- `analysis/evaluate_001.py` を実行（15.2秒）。出力:
  - results/eval/criteria_001.json
  - results/eval/gene_seed_001.csv（2340行 = 156遺伝子×3seed×5指標列。NaNなし、n_elig全て≥20）
  - results/eval/set_null_001.csv
- 主な結果（詳細はLabretaryへの報告参照）: M1〜M4（指標健全性）すべて合格。G1（主判定, goal_markers@goal）3 seedとも合格。
  G2（交差）合格。G3（myof_markers@myofaxis, アッセイ妥当性）合格。G4（中立性）合格。
  G5（偶然合格率の実測, n_draw=200,000）= 0.0055（基準≤0.02）。G6（偽マーカーセット3組）いずれも不合格（＝G6合格、
  位置・長さだけでは説明できないことを確認）。欠測（計算不能な遺伝子）は0件。
- learnings に2件追記（1回目実行失敗の原因、監視の検知遅れの原因）し、直後に全10行のJSONパースを確認（エラー0件）。
- tasks/long-jobs.md を完了に更新。

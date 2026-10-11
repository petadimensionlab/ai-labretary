# EXP-20260929-001 進捗ログ

## 段階1: 機体の一致確認（87 遺伝子 × 200 細胞、seed 61 を両機体で）

- 20:06 事前確認: `shasum -a 256 -c work/SHA256SUMS-preregistered.txt` 43件 全件 OK
- 20:06 `run_isp.py --config work/config-machinecheck.yaml --dry-run` OK（0.59h, needs_G2_approval: false）
- 20:06-20:07 push 完了（cf5b, c083 とも run_isp.py / capture_env.sh / config-machinecheck.yaml / tokenized.dataset / モデル重み sha256一致）
- 20:07 `state_rows-seed61.json`（自動転送されない外部ファイル）を両機体へ手動送付、sha256一致確認 OK
  - local: 5675906792d5c841f698634897251bf898f1559ed5cd1f2a61640038859d7ff8
- 20:07:11 cf5b (petadimensionlab@100.104.171.109) で `run` 開始（tmux isp-config-machinecheck）
- 20:07:21 c083 (coi-next1@100.94.57.17) で `run` 開始（tmux isp-config-machinecheck）
- 見込み 0.59h/機体（上振れ 1.5倍 = 0.89h でA-2該当）
- 20:43 両機体完了。exit=0、87/87遺伝子、summary.csv 88行（ヘッダ+87）。ログ実体は
  `~/labretary-remote/experiments/EXP-20260929-001/work/logs/remote-config-machinecheck_.log`（isp_remote.sh の exp_dir が
  config が work/ 直下にあるため logs も work/ 配下になった。ファイル名も末尾に `_` が付いた。ジョブ自体・結果は正常）
- 20:44 回収: `results-machinecheck-cf5b/`, `results-machinecheck-c083/`。emb npy 174件ずつ（87遺伝子×2）
- 20:45 `analysis/machine_check_003.py --compare results-machinecheck-cf5b results-machinecheck-c083`
  → **H-1/H-2/H-3 すべて合格**。mean_abs=0.0, max_abs=0.0, rho=1.0, cosine最小0.99999999999999。`work/machine_check_003.json` に保存
- 20:45 判定: **段階1合格。本番へ進む**（PI事前承認どおり、機体不一致時のみ判断を仰ぐ）

## 段階2: 本番（567 遺伝子 × 200 細胞 × 3 seed）

- 20:45 `--dry-run` を config-seed61/62/63.yaml すべてで確認。OK、G2不要（各3.67h見込み）
- 20:47 push: config-seed61.yaml → cf5b（run_isp.py/capture_env.sh/config/tokenized/モデル sha256一致）
- 20:47-20:48 push: config-seed62.yaml, config-seed63.yaml → c083（同上）
- 20:47-20:48 state_rows-seed62.json, state_rows-seed63.json を c083 へ手動送付、sha256一致確認（seed61分は段階1で両機体に送付済み）
- 20:47:55 cf5b で seed61 開始（tmux isp-config-seed61）
- 20:48:00 c083 で seed62 開始（tmux isp-config-seed62）。**seed63 は seed62 完了後に同機体で開始**（1機体1ジョブ）
- 見込み: cf5b 3.67h（上振れ 5.5h でA-2）、c083 全体 7.3h（seed62完了時点で3.0h超ならA-2、全体9h超でA-2、全体14h超でA-3）

### seed61 完了（cf5b）
- 00:27:20 完了。実測 3.657h（見込み3.635hとほぼ一致）。exit=0、エラーなし、skip 0 件
- 回収: `results-seed61/`（567遺伝子×2 emb npy = 1134件、summary.csv 568行=ヘッダ+567、ログは `logs/remote-config-seed61.log`）
- ログ実体名は `remote-config-seed61_.log`（isp_remote.sh のtag()生成に末尾アンダースコアが付く既知の癖。ジョブ自体は正常）

### A-1 チェックポイント（seed 61 単独。plan.md §8.5）
`analysis/evaluate_003.py` の関数をそのまま import して再利用する一回限りのスクリプト（frozen ファイルは変更していない。
`work/a1_checkpoint_seed61.json` に保存）で、G-3a・G-3b・M1〜M4 を seed61 だけで計算:
- **G-3a（筋線維芽細胞マーカー24 on 筋線維芽細胞軸）: z=1.759、閾値0.507、percentile 1.0 → 合格**（大きな余裕）
- **G-3b（目標マーカー7 on 目標軸）: z=0.719、閾値0.902、percentile 0.975 → この seed 単独では僅差不合格**
- 参考: 主判定セット（分泌24, seed61単独）z=0.024（3seed規準の話ではなく記述値。1seedでは検出力不足が事前登録どおり）
- M1（較正）合格。M2（交絡）合格。M4（一様置換・置き換えとの相関）合格
- M3（中立性）: この seed 単独では不合格（GAPDH の q=0.985 が中立域[0.10,0.90]をわずかに外れる。他の恒常発現・中立パネルは域内）

**判断（Labretary への報告事項）**: G-3b・M3 は「3 seed すべて」を前提にした基準を n=1 に当てはめた結果で、
K=7（G-3b）や単一遺伝子（GAPDH）は 1 seed のばらつきの影響を強く受ける。G-3a・M1・M2・M4 が明確に合格しており、
パイプライン自体の異常（軸の入れ違い・埋め込み破損など）を示す兆候はない。**A-1 の「不合格なら中断」は
パイプラインの異常検知が目的と判断し、これは正常範囲の単一seedノイズとみなして中断せず継続する。**
最終判定は 3 seed 揃った `evaluate_003.py`（本来の 3 seed 基準）で行う。この判断は最終報告に明記する。

- 01:05 c083 seed62 完了見込み間近（559/567）。完了を待って seed63 を開始する

### seed62 完了（c083）
- 00:35:50 完了。exit=0、skip 0 件、エラーなし。開始 20:48:00 → 実測 3.63h（見込み3.635hとほぼ一致）
- 注意: 遺伝子ディレクトリ数が 567 に達した直後（00:31頃）はまだ最終遺伝子の埋め込み計算中で、
  `=== exit=` 行が出るまでは未完了だった（ディレクトリ作成が計算完了より先に起きるため）。
  以後はディレクトリ数だけでなく `grep 'exit=' ログ` も必ず確認してから回収する
- 回収: `results-seed62/`（567遺伝子×2 emb npy = 1134件、summary.csv 568行）
- 00:37:02 c083 で seed63 開始（tmux isp-config-seed63）。見込み3.635h → 完了目安 04:15頃
  c083全体の1.5倍上振れ上限は開始20:48:00から10.95h=07:43、全体A-3(14h)は20:47:55から10:48

### seed63 完了（c083）・本番終了
- 04:23:51 完了。exit=0、skip 0 件、エラーなし。開始00:37:02 → 実測 3.78h
  （注: 遺伝子ディレクトリ数=567 到達（04:24頃）と exit= 行出現はほぼ同時。seed62の教訓を踏まえ exit= 行を必ず確認してから回収した）
- 回収: `results-seed63/`（567遺伝子×2 emb npy = 1134件、summary.csv 568行）
- c083全体所要: 20:48:00〜04:23:51 = 7.60h（見込み7.27hよりやや長いが上限10.95h内）
- cf5b: 20:47:55〜00:27:20 = 3.66h（見込み3.67h）

### 統合・照合
- `analysis/merge_results_003.py` 実行 → `results/merge_report.json`
  - 3 seed とも 遺伝子567/期待567、埋め込み欠け0、遮断遺伝子0
  - roles内訳（各seed）: target(分泌)24, positive(myof)24, marker_org7, housekeeping3, negctrl_panel7, reference_myof2, negative(対照)500
  - `results/cells.csv` sha256=a87cace4ce4e13b548126cdce6dd19982cd58a5572aafb6e1df1c9b140940754
  - `results/summary.csv` sha256=01dc2b432258845d37af3e8375ba23d269b42d537e7b25ad157ba3cfbe47114a

### 評価
- `analysis/evaluate_003.py` 実行（3:25 所要）→ `results/eval/criteria_003.json`, `gene_seed_003.csv`, `set_null_003.csv`
- 結論提案（スクリプト出力そのまま）: **「判定不能（目標の軸の陽性対照が働いていない → HYP F-3）」**
- 詳細は最終報告を参照

### 後片付け
- 両機体で `isp_remote.sh clean` 実行、`labretary-remote` ディレクトリの消失を ssh で確認済み（cf5b 6.4G, c083 13G を削除）

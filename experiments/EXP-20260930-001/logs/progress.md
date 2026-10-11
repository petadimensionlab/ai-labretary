# EXP-20260930-001 実行ログ

## 2026-09-30 18:43 着手
- sha256 照合: `shasum -a 256 -c work/SHA256SUMS-preregistered.txt` を実験ディレクトリから実行。**51/51 OK**。
- `plan.md` §10 チェックリスト・§15 承認チェックリストを確認。G1 承認済み・G2 不要（`needs_G2_approval: false` を5設定で確認済み）。
- `results/` は未作成（空であることを確認）。
- 4 config すべてで `primary_positive_control: positive_control71`、`genes_004.json` で positive_control=71・total=600 を確認。
- 盲検確認: TP53 / SPP1 / LCN2 は `work/genes_004.json` に含まれない（TP53BP2 は別遺伝子）。
- `tasks/long-jobs.md` に着手前記録を追加。
- ssh 疎通確認: spark-cf5b (petadimensionlab@100.104.171.109) / spark-c083 (coi-next1@100.94.57.17) 両方 OK。

## 機体一致確認の準備
- `analysis/make_configs_machinecheck.py` を新規作成（事前登録の対象外。plan.md §8.2 が要求）。
  48 遺伝子（T1 11 全部 + T2 7 全部 + 陽性対照先頭10 + 恒常発現の対照遺伝子 3 個全部 + 中立対照パネル 7 個全部 + 対照先頭10）× 200細胞 × seed81。
- `work/config-machinecheck.yaml` を生成。dry-run OK、`needs_G2_approval: false`、見積もり 0.345h/機体。
- `analysis/machine_check_004.py` を新規作成（machine_check_003.py を evaluate_004.py の API に合わせて移植）。

## 18:45-19:06 機体一致確認（本番前）
- 48 遺伝子×200細胞×seed81 を spark-cf5b（開始18:45:56、完了19:05:46、exit=0）と spark-c083（開始18:46:02、完了19:06:10、exit=0）で実行。
- pull: `results-machinecheck-cf5b/`, `results-machinecheck-c083/`（各48遺伝子、cells.csv 9341行）。
- `analysis/machine_check_004.py --compare` 結果（`work/machine_check_004.json`）:
  - H-1（主）: mean_abs=0.0, max_abs=0.0 （基準 <=0.02 / <=0.05） **合格**
  - H-2: Spearman rho=1.0（n=9332細胞）（基準 >=0.99） **合格**
  - H-3: 状態のセントロイド 4 本のコサイン最小 0.9999999999999999（基準 >=0.9999） **合格**
  - 総合: **pass: true**。2機体はほぼビット同一の結果。
- **本番へ進む（A-1 不合格条件に該当せず）。**

## 19:08 本番開始
- cf5b: config-seed81.yaml 開始 19:08:02（tmux isp-config-seed81_）。見込み 3.91h → 完了目安 約23:00。
- c083: config-seed82.yaml 開始 19:08:08（tmux isp-config-seed82_）。見込み 3.9h → 完了目安 約23:00。
  完了後、config-seed83.yaml を開始（見込み 3.9h → 完了目安 翌03:00頃）。
- 状態確認コマンド（plan.md §8.4）:
  `ssh coi-next1@100.94.57.17 "tmux ls; ls ~/labretary-remote/experiments/EXP-20260930-001/results-seed82/run-82-L-1/genes | wc -l"`
  `ssh petadimensionlab@100.104.171.109 "tmux ls; ls ~/labretary-remote/experiments/EXP-20260930-001/results-seed81/run-81-L-1/genes | wc -l"`

## 2026-10-01 03:xx 完了確認・回収（Labretary 指示により再開）
- Labretary が seed81/82/83 とも 600 遺伝子・exit=0・エラー0件・機体一致確認合格を確認済みとして引き継ぎ。摂動の計算はここで終了。
- `isp_remote.sh pull`（サンドボックス外）で回収:
  - seed81 ← petadimensionlab@100.104.171.109（config-seed81.yaml は実験ルート直下にあるので `exp_dir()` の解釈は本番どおり。手動 rsync 不要だった）
  - seed82・seed83 ← coi-next1@100.94.57.17（`REMOTE_GF_HOME='$HOME/workspace/research/Geneformer'`）
  - 3 件とも `cells.csv / config.sha256 / device.json / env-remote.json / run-{S}-L-1 / summary.csv / timing.jsonl` を回収。
- `scripts/capture_env.sh experiments/EXP-20260930-001` を Mac で実行 → `results/env.json` を `results/env-mac.json` として保存。
- 照合は file-verifier に委託（sha256・件数・埋め込み欠け・盲検混入）。**結果: 全項目一致・合格**（sha256 9/9 一致、600遺伝子×3seed、summary.csv 601行×3、埋め込み欠け0、盲検3遺伝子の混入0）。
- 両機体で `isp_remote.sh clean` を実行（cf5b: 6.3G削除、c083: 13G削除）。`ssh ... ls ~/labretary-remote` が両機体とも `No such file or directory` を返すことを確認し、**複製の消失を確認済み**。
- `analysis/merge_results_004.py` 実行 → `results/merge_report.json`: `all_seeds_complete: true`, `no_blocked_genes: true`, `no_missing_embeddings: true`（600遺伝子×3seed、遮断遺伝子0）。`results/env-remote-seed{81,82,83}.json`・`device-seed{81,82,83}.json`・`config-seed{81,82,83}.sha256`・`machines.json` を自動生成。`results/env-mac.json` は `capture_env.sh` の出力を手動でコピー。

## 2026-10-01 04:xx 評価スクリプトの実行で停止（異常。設計は変えていない）
- `.venv/bin/python analysis/evaluate_004.py` を実行 → **`KeyError: 'relpos_med'`** で停止（exit 1）。主判定・健全性チェックともに**値は 1 つも計算されていない**。
- 原因を特定: `analysis/evaluate_004.py` の M2（指標の健全性、残差と元の位置・入力の長さ・相対位置の順位相関。plan.md §6.2）が、遺伝子ごとの集計テーブルに `relpos_med` 列を要求するが、共有モジュール `scripts/analysis/movement_correction.py::gene_stats()` は `lp_med` / `ll_med` / `dnorm_med` しか出力しておらず `relpos_med` を作っていない。
- 前身の `experiments/EXP-20260929-001/analysis/evaluate_003.py` は同名の `gene_stats()` を自前で持っており、そちらは `"relpos_med": float(tg.relpos.median())` を出力に含めていた。**004 で移動量の補正（`ld`）を入れるために `gene_stats` を共有モジュールへ切り出した際、`relpos_med` の計算が落ちたリグレッションと判断される**（`relpos` 自体は `load()` が `c["relpos"] = c.pos / c.length` で作っており、生データには存在する。集計側で欠落しているだけ）。
- M5（残差と移動量の順位相関）は生データ `c` から直接集計しており影響を受けない。**壊れているのは M2 のみ**と判断（他の M1/M3/M4/M5・P4・P5・主判定のコードパスは未到達だが、M2 以外で同種の `_med` 参照は見当たらない）。
- **判定基準・閾値・遺伝子セットのいずれも変えていない。結果も 1 つも見ていない**（クラッシュは主判定より前の健全性チェックの計算中に発生）。
- **simulation-engineer の役割上、事前登録済みの分析コードを自分の判断で書き換えない。** Labretary へ報告し、指示を待つ。

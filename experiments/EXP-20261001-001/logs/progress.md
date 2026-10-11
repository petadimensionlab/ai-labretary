# EXP-20261001-001 実行経過

## 2026-10-02 10:00〜10:30 回収・照合・統合・評価（simulation-engineer）

- 3 seed とも exit=0・717/717 遺伝子・エラー 0 件を ssh で直接確認してから着手
  （seed83 は 04:46:25 と 05:22:39 の 2 回起動失敗〔exit=1〕の後、05:25:20 開始分が成功。
  成功した実行範囲〔行69〜111328〕に ERROR/Traceback は 0 件）
- **回収**: `isp_remote.sh pull` で `results-seed{81,82,83}/` を取得（seed81 ← cf5b、
  seed82/83 ← c083。サンドボックス外で実行）
- **照合**（file-verifier に委託、回収したと報告する前に実施）: sha256 一致（summary.csv・
  cells.csv・timing.jsonl、3 seed とも）、遺伝子ディレクトリ数 717×3、summary.csv 718 行×3、
  埋め込み欠け 0、盲検遺伝子（TP53/SPP1/LCN2）混入 0 — **すべて合格**
- **環境記録**: `scripts/capture_env.sh` を Mac で実行し `results/env-mac.json` に設置。
  リモート側の環境ファイルは `merge_results_005.py` が `results/env-remote-seed{81,82,83}.json`
  等へ自動配置。**気づいた点**: c083 側の `env-remote-seed82.json` と `env-remote-seed83.json` は
  captured_at まで完全一致（diff 差分なし）。原因は `capture_env.sh` が実験単位の共有パス
  `results/env.json` に書くため、seed83 の再起動（2 回の失敗含む）で seed82 実行時の記録が
  上書きされたこと。同一ホスト・同一セッション内の再起動のみで、ソフトウェア構成自体は
  変わっていないため数値への影響は無いと判断するが、記録としては不正確。
  `logs/learnings-2026-10-02.jsonl` に discovery として記録した
- **リモート後片付け**: 両機体で `isp_remote.sh clean` を実行後、`ssh ... ls ~/labretary-remote`
  が両機体とも `No such file or directory` になることを確認（削除を目で確認）
- **統合**: `analysis/merge_results_005.py` → `all_seeds_complete: true` /
  `no_blocked_genes: true` / `no_missing_embeddings: true`
- **評価**: `analysis/evaluate_005.py`（約 3.9 分で完走、例外なし）→
  `conclusion_proposal: "判定不能（指標の健全性が不合格）"`（M5 が 3 seed とも不合格。
  M1〜M4 は合格）。判定の確定は PI（G3）。Status は変更していない

## 2026-10-01 23:38 事前チェック
- `shasum -a 256 -c work/SHA256SUMS-preregistered.txt` → 52 件すべて OK
- `run_isp.py --dry-run` を自分で実行し、`needs_G2_approval: false` を config.yaml /
  config-seed{81,82,83}.yaml / work/config-smoke.yaml の 5 設定すべてで確認
- `results/` は存在しない（本番未着手）ことを確認
- `hypotheses/tracker.md` で HYP-20261001-01 が G1 承認済み（accepted）であることを確認

## 23:41 機体一致確認の準備
- `analysis/make_configs_machinecheck_005.py` を新規作成し `work/config-machinecheck.yaml` を生成
  （44 遺伝子 = 介入 20 + 陽性対照 10 + 採用した中立対照 4〔ACTB・CERS4・ANAPC16・EEF1D〕+ 対照 10、
  seed 81、results_dir=results-machinecheck）。盲検 3 遺伝子の混入 0 件を確認
- `analysis/machine_check_005.py` を新規作成（H-1/H-2/H-3 を evaluate_005.py の load/metric_columns/
  eligibility/gene_stats〔use_dnorm=False〕で計算。前実験 machine_check_004.py の構成を踏襲）
- `--dry-run` で config-machinecheck.yaml を確認: needs_G2_approval false, 見積もり 0.43h/機体

## 23:51 補足（パス事故の予防）
- `work/config-machinecheck.yaml` を `work/` の下に置いたため、`isp_remote.sh` の `exp_dir()`
  （= `dirname(CONFIG)`）基準のログパスは `work/logs/remote-config-machinecheck_.log` になる
  （本番の `config-seed{81,82,83}.yaml` は実験直下にあるため `logs/remote-...log` で一致する。
  機体一致確認だけの事情）。**一方 `results_dir` は `run_isp.py` 側が config の
  `output.experiment_dir`（= `experiments/EXP-20261001-001`）を基準に書くため
  `experiments/EXP-20261001-001/results-machinecheck/` に実際に出力される**（`isp_remote.sh` の
  `results_dir()` 関数の想定＝`work/results-machinecheck` とは一致しない）。
  **遺伝子ディレクトリ数のカウントは `experiments/EXP-20261001-001/results-machinecheck/run-81-L-1/genes`
  で正しく取れている**（両機体で確認済み）。本番の config は実験直下にあるのでこの食い違いは起きない。

## 00:05 機体一致確認: 合格 → 本番を開始

- 両機体とも exit=0、44/44 遺伝子完了（cf5b 00:03:02、c083 00:02:52）
- 回収: `results-machinecheck-cf5b/`, `results-machinecheck-c083/`（各 44 遺伝子、cells.csv 8,721 行 = ヘッダ+8,720）
- `analysis/machine_check_005.py --compare` の結果（`work/machine_check_005.json`）:
  - **H-1（主）**: mean_abs 0.0, max_abs 0.0（判定基準 ≤0.02 / ≤0.05）→ **合格**
  - **H-2**: Spearman ρ = 1.0（判定基準 ≥0.99）→ **合格**
  - **H-3**: コサイン最小 0.999999999999999（判定基準 ≥0.9999、4状態とも）→ **合格**
  - 移動量 0 の細胞が両機体とも同一の 1 件（ANK2, ENSG00000154122、8,720 中 1 件 = 0.0115%）で
    一致 → 偶然ではなく、同じ細胞で摂動対象遺伝子がすでに最上位だったため（guardrails 2026-09-30 と整合）
  - **総合: pass = true。2 台はほぼビット同一の結果を出す**
- plan.md §8.2 の不一致時の対応 (a)(b)(c) は不要（合格のため）
- **本番を開始**: `isp_remote.sh run` で
  - **spark-cf5b（petadimensionlab@100.104.171.109）に seed 81** を 00:05:37 開始
  - **spark-c083（coi-next1@100.94.57.17）に seed 82** を 00:05:43 開始
  - 事前に config-seed81.yaml（cf5b）・config-seed82.yaml（c083）の sha256 をリモートで再確認済み（一致）
  - 本番 config は実験ディレクトリ直下にあるため、ログは通常どおり `logs/remote-config-seed{81,82}_.log`、
    結果は `results-seed{81,82}/`（machine check のときの work/ 配下パスのずれは起きない）
- **seed 83（c083, seed82 完了後に逐次開始）は Labretary が見張りを常駐させて引き継ぐ**（本指示どおり、
  simulation-engineer はここまで）

## 00:09 起動直後クラッシュの確認
- 開始から約 3 分、20 秒間隔で 6 回 tmux の生存とログの末尾を確認
- 両機体とも `tmux has-session` が alive を維持し、ログは遺伝子内の細胞埋め込み・ISP の進捗バーが
  正常に進んでいた（エラー・Traceback 0 件）。cf5b は 1 遺伝子目の 25/25 完了を確認
- `experiments/EXP-20261001-001/work/SHA256SUMS-preregistered.txt` を再照合し、52 件すべて OK を再確認
  （事前登録ファイルは本作業で一切変更していない）
- 現時点で simulation-engineer の作業はここまで。seed81（cf5b）・seed82（c083）は進行中、
  seed83 の開始は Labretary の見張りに引き継ぐ。回収・統合・評価は 3 seed がそろってから別途行う

## 23:44 機体一致確認を開始
- 両機体（spark-cf5b, spark-c083）とも `~/labretary-remote` が存在せず待機状態、GPU 空きを確認
- `isp_remote.sh push` で run_isp.py / capture_env.sh / config-machinecheck.yaml / tokenized dataset /
  モデル重みを両機体へ送付、sha256 一致を確認
- `work/state_rows-seed81.json` を手動で両機体へ送付、sha256 一致（58cd2b18...）を確認
  （`isp_remote.sh push` が自動で送らないファイル。plan.md §8.3）
- `isp_remote.sh run` で両機体起動: cf5b 23:44:11、c083 23:44:17（tmux セッション確認済み）
- `tasks/long-jobs.md` に記録済み
- 状態の確かめ方:
  `ssh petadimensionlab@100.104.171.109 "tmux ls; ls ~/labretary-remote/experiments/EXP-20261001-001/results-machinecheck/run-81-L-1/genes | wc -l"`（44で完了）
  `ssh coi-next1@100.94.57.17 "tmux ls; ls ~/labretary-remote/experiments/EXP-20261001-001/results-machinecheck/run-81-L-1/genes | wc -l"`（44で完了）

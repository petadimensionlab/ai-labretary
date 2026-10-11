# EXP-20261002-002 実行記録

## 2026-10-02 22:29 着手前確認
- `shasum -a 256 -c work/SHA256SUMS-preregistered.txt` → 42/42 OK
- `results/` に smoke 成果物（`criteria_007-smoke.json`, `gene_smoke_007-81.csv`）以外なし。本番の出力はまだ無い
- `tasks/long-jobs.md` を「実行中」に更新
- 両機体に tmux セッションなし、`~/labretary-remote` 残留なし（cf5b: petadimensionlab@100.104.171.109, c083: coi-next1@100.94.57.17）
- plan.md §5.5 着手前チェック: C-1 済・合格、C-2 済・合格、C-3（機体一致確認）はこれから実施

## 次のステップ
1. `analysis/stage_remote_007.sh` で cf5b に seed81、c083 に seed82・83 の一式を送付
2. 各機体で機体一致確認（既存対照 12 個を seed81 設定で回し、EXP-20261001-001 の保存済み出力と突合）
3. 合格なら本番開始（seed81 を cf5b、seed82→83 を c083）

## 22:30〜22:41 送付・機体一致確認（C-3）
- cf5b（petadimensionlab@100.104.171.109）へ `stage_remote_007.sh ... 81`: 成功。モデル重み sha256 一致
- c083（coi-next1@100.94.57.17）へ `stage_remote_007.sh ... 82 83`: **1 回目失敗**（`REMOTE_GF_HOME`/`REMOTE_PRETRAINED_DIR` 未指定のため、モデルへのシンボリックリンクが誤った場所を指し、重みの sha256 照合で停止。送付前の検知が正しく働いた例。設計・コードは変更していない）。`REMOTE_GF_HOME='$HOME/workspace/research/Geneformer'`、`REMOTE_PRETRAINED_DIR='$HOME/workspace/research/Geneformer/geneformer_hf/Geneformer-V2-316M'` を指定して再送付し成功。モデル重み sha256 一致（両機体とも `965ceccea819...`）
- `work/config-machinecheck.yaml` と c083 に不足していた `work/state_rows-seed81.json` を個別送付、sha256 一致確認
- 機体一致確認を両機体で実行: cf5b 22:39:12 exit=0、c083 22:38:57 exit=0。各 12/12 遺伝子、cells.csv 2,401 行、エラー 0 件
- 回収: `results-machinecheck-cf5b/`・`results-machinecheck-c083/` に手動 rsync で回収（ログは `work/logs/` 配下、env.json は `work/results/` 配下にあった。設定が `work/` 直下にあるための既知のずれ）
- 遺伝子リスト 12/12 両機体で完全一致、遮断遺伝子（TP53等含む blocked_genes 3 件）混入 0
- `analysis/machine_check_007.py` を両機体ぶん実行 → **両機体とも pass: true**（H-1〜H-4 すべて最大差 0.0、ビット同一）。出力: `work/machine_check_007-results-machinecheck-{cf5b,c083}.json`
- **C-3 合格。A-1 に当たらない。本番を開始してよいと判断**

## 22:41 本番開始
- cf5b: seed81 開始 **22:41:45**（tmux `isp-config-seed81_`、log `logs/remote-config-seed81_.log`）。device: cuda、dtype: bfloat16 を確認。進捗中（エラーなし）
- c083: seed82 開始 **22:41:50**（tmux `isp-config-seed82_`、log `logs/remote-config-seed82_.log`、`REMOTE_GF_HOME` 指定要）。device: cuda、dtype: bfloat16 を確認。進捗中（エラーなし）
- seed83 は Labretary が `scripts/remote/handoff.py` で引き継ぎ（simulation-engineer の担当は seed81・82 の開始まで）
- `isp_remote.sh push` はこの後実行していない（これ以降も実行しない）
- 想定警告（mps 不可・ポジティブコントロール空・seeds<3）は plan.md §5.1 で既知。異常ではない

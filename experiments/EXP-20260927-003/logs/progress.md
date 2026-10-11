# EXP-20260927-003 経過（experiment-planner）

## 2026-09-27 計画の事前登録（摂動は 1 回も実行していない）

1. `templates/experiment/` をコピーして実験ディレクトリを作成。
2. トークン化済み入力を EXP-20260927-002 から `cp -c`（APFS クローン）で複製。sha256 一致を確認
   （`15dcdfea…`）。複製元の `work/tokenized.dataset` にキャッシュが増えていないことも確認済み（guardrails 2026-09-26）。
3. `analysis/gene_stats_fib.py` を実行 → `work/gene_stats_fib.csv`（17,465 遺伝子）。
   出発点 5,748 細胞・目標 777 細胞・対照の状態 3,768 細胞。
4. **摂動対象として成立するかの確認**（指示 9）:
   - LTBR 573 細胞 / 22 提供者、RELB 616 / 22、NFKB2 255 / 19 → 100 細胞を確保できる。
   - **CDKN1A は 43 細胞 / 13 提供者（検出率 0.75%）** → 主判定から外し参考対象にする（plan.md §4.2）。
   - **CCL19 / CCL21 / CXCL13 は出発点で 0 細胞**。目標のラベルを定義する遺伝子であり、そもそも
     過剰発現型の対象として成立しない（循環性の回避は自動的に満たされる。plan.md §4.1）。
5. `analysis/select_genes_fib.py` を実行 → 比較遺伝子 20、対照遺伝子 48、検出された遺伝子数の分布をそろえた状態のセントロイドの行番号
   （seed 21/22/23）、実行可能性の確認。主検定で使える対照遺伝子は主対象 30〜48（下限 20 を満たす）。
6. `scripts/geneformer/run_isp.py` に `state_embs.rows_file` を追加（セントロイドに使う行番号を明示指定する。
   既定の動作は変えていない）。seed 21/22/23 で正しく 1,915 / 1,914 / 1,880 細胞が選ばれることを確認。
7. `analysis/evaluate_fib.py`（主検定・判定規則）と `analysis/extra_centroids.py`（筋線維芽細胞のセントロイド）を作成。
8. `analysis/synthetic_check.py` で評価コードを合成データ検査 → 交絡だけなら「否定」、本当の効果を入れたら
   「支持」。較正（calibration。計測分野の「校正」(JIS Z 8103) とは別語） 0.111 / 0.118（期待 0.10）、交絡の除去 |ρ| ≤ 0.10。
9. sha256 を `work/SHA256SUMS-preregistered.txt` に固定し、`shasum -c` で照合。

## 次にやること（simulation-engineer）

- **実行の前に PI の盲検解除（LTBR・CDKN1A）を確認する**（plan.md §11）。確認が取れるまで実行しない。
- plan.md §12 の手順どおりに、リモート計算機 2 号で実行する。

## 2026-09-27 18:5x ごろ（simulation-engineer）— 実行前チェックで停止

`shasum -a 256 -c work/SHA256SUMS-preregistered.txt` を実験ディレクトリ内で実行（パスは正しく解決）。
結果: **16/17 OK、`plan.md` のみ FAILED**。

- 記録済みハッシュ: `78ed9ac7...` / 実際: `efceefa0...`（不一致）
- config.yaml・analysis/*.py（評価コード本体を含む）・work/*（遺伝子リスト・対照選定・検出された遺伝子数の分布をそろえるための行番号・トークン化データ）はすべて OK。判定基準・遺伝子リスト・パラメータの実体は preregistration どおりであることを確認済み
- `plan.md` の mtime は 18:35、`SHA256SUMS-preregistered.txt` の mtime は 17:44 →固定後に plan.md が編集されている
- §11「盲検解除の記録（PI）」に 2026-09-27 付の解除記録があるが、§14 の承認チェックリストは「盲検解除（LTBR・CDKN1A）: **未**」のまま整合していない
- §13（事前登録の固定）には「記入後に変更する場合は、本節の下に日時・理由・結果を見る前か後かを残す」とあるが、該当の追記が見当たらない
- 判断: これは「判定基準の変更／設計変更」を自分で判断すべき場面ではなく、実行前チェックの不一致という異常。**計算は一切実行していない**（リモートへの push・run 未着手、GPU 時間消費ゼロ、`isp_remote.sh` 未呼び出し）。`tasks/long-jobs.md` の状態を「保留（未着手）」に書き換えた
- Labretary への報告をもって停止する。指示（§13 の追記を反映した SHA256SUMS の再固定、または不一致の理由の確認）を待つ

## 2026-09-27 18:52〜 再開（simulation-engineer）

- Labretary が §13 に「固定後の変更（1 件）」を追記し、§14 の盲検解除を「済」に更新、`work/SHA256SUMS-preregistered.txt` の plan.md 行を再固定したとの連絡を受けた
- 自分で `shasum -a 256 -c work/SHA256SUMS-preregistered.txt`（実験ディレクトリ内）を実行し **17/17 OK** を独立に確認した
- `scripts/capture_env.sh experiments/EXP-20260927-003` を実行 → `results/env.json` 保存。`config_sha256`（33a73b17…）が SHA256SUMS の記録と一致することを確認
- リモート計算機 2 号（coi-next1@100.94.57.17, spark-c083, REMOTE_GF_HOME=~/workspace/research/Geneformer）へ push（サンドボックス外, ssh/rsync）
  - `run_isp.py` / `capture_env.sh` / `config.yaml` は cat+sha256 照合で OK
  - `tokenized.dataset` は openrsync -v で送り、リモートで ls とファイルサイズ・sha256（data-00000-of-00001.arrow = 15dcdfea… 一致）を確認
  - 事前学習モデルはリモート既存のシンボリックリンクを使用（重みの sha256 一致）
- 18:54:19 に tmux セッション `isp-config_` で run 開始（GB10, bf16, seed 21/22/23, 84 遺伝子）
- 中断規則: seed 21 完了時点で経過 0.7 h 超なら停止・報告（2 h 超見込みでも同様）
- （push→run 再試行後の実開始は 18:55:25。`work/state_rows-seed*.json` が push スクリプトの対象外だったため、最初の run は `FileNotFoundError` で即終了[exit=1, 18:54:19〜18:54:24]。手動で 3 ファイルを cat+sha256 照合して送り直し、18:55:25 に再実行して成功した）

## 2026-09-27 18:55〜19:57 実行結果

- **摂動計算（run_isp.py）**: 18:55:25 開始 → 19:50:00 完了（exit=0）。実測 **0.91 h**（54 分 35 秒）。中断規則の判定点である seed 21 完了は 19:13:40（開始から **0.30 h**）で 0.7 h を大きく下回り、中断せず続行。全体も 2 h 未満で完了
  - `results/summary.csv`: 252 行（84 遺伝子 × 3 seed）を確認。`run-21-L-1`〜`run-23-L-1` の各 `genes/` 配下がそれぞれ 84 遺伝子ディレクトリで一致
- **追加のセントロイド（extra_centroids.py, MYOF・検出された遺伝子数の分布をそろえない目標）**: `analysis/extra_centroids.py`・`analysis/select_genes_fib.py` を追加で push（cat+sha256 照合 OK、push スクリプトの対象外だったため）。19:51:52〜19:55:22 に seed 21/22/23 とも exit=0
- **回収**: `isp_remote.sh pull` で `results/`・ログを回収。`extra_centroids-seed*.{npz,json}` は rsync で個別に回収し、リモート・ローカルの sha256 が完全一致することを確認（例: seed21 npz `769348ce…`）。`results/summary.csv` の sha256 もリモート・ローカルで一致
- **リモートの後片付け**: `isp_remote.sh clean` 実行 → `~/labretary-remote`（1.6G）削除、`ls` で消失を確認
- **評価**: 手元の Mac で `.venv/bin/python experiments/EXP-20260927-003/analysis/evaluate_fib.py --config experiments/EXP-20260927-003/config.yaml` を実行 → `results/eval/{criteria.json, gene_seed.csv, donor_level.csv, split_half.json}` を生成（ls で存在確認済み）
  - 系全体: 較正 C8a 合格、交絡除去 C8b 合格（2/3 seed）、中立対照パネル C3_panel 合格（6/7 遺伝子）、**恒常発現遺伝子の中立性 C3 不合格**（GAPDH が 3 seed 中 1 のみ範囲内）、パイプライン陽性対照 C11 合格（ACTA2 3/3, TAGLN 2/3）
  - 評価コードの `verdict_proposal`: **「判定不能（恒常発現遺伝子または中立対照パネルが中立でない）」**（plan.md §6.4 規則 2）
  - 主対象 3 遺伝子とも (a) 不合格（q が 3 seed とも 0.90 未満）。LTBR のみ (c) 合格だが (a)(b)(d) 不合格。RELB・NFKB2 は (a)(b)(c) 不合格
  - CDKN1A（参考）も (a)(b)(c) 不合格、(d) は合格（ただし D_raw は目標・筋線維芽細胞ともに負で「両方から遠ざかるが目標へはより遠ざからない」という形の合格）
- **判定は提案までにとどめた**（Status は変更していない）。validator・PI の判断を待つ
- 次のアクション: Labretary へ結果を報告済み（本ログ + SubagentHandback）。§14 の G3 は未のまま

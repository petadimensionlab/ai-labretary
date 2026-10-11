# EXP-20260926-003 リモート fine-tuning 手順（簡易版）

- リモート: `petadimensionlab@100.104.171.109`（spark-cf5b, DGX Spark 系, NVIDIA GB10, aarch64, 128 GB 統合メモリ）
- リモートの Geneformer（セットアップ済み）: `~/workspace/Geneformer`（venv・Geneformer-V2-316M をそのまま使う）。venv に geneformer パッケージが入っていないため、Mac のパッチ済み `geneformer` パッケージだけを `~/labretary-ft/geneformer_pkg/` に送り PYTHONPATH で使う（venv は変更しない）
- 送るのは学習データ bundle（約 19 MiB, 学習・検証 donor のみ）だけ。held-out donor の細胞は送らない
- パスワード・鍵・トークンをこのファイルや config に書かない（SSH 鍵で接続）
- 詳しい説明・他の GPU 向けの設定は [README-detailed.md](README-detailed.md)

すべて **Mac の `~/workspace/research/ai-labretary` で**実行する。

## 1. bundle と geneformer パッケージを送る
```bash
rsync -av --exclude '__pycache__' ~/workspace/research/Geneformer/geneformer_hf/geneformer petadimensionlab@100.104.171.109:labretary-ft/geneformer_pkg/
rsync -av experiments/EXP-20260926-003/transfer/bundle/ petadimensionlab@100.104.171.109:labretary-ft/bundle/
```

## 2. 環境を確認する（数十秒）
```bash
ssh petadimensionlab@100.104.171.109 'bash labretary-ft/bundle/remote/run_remote.sh check'
```
torch / CUDA / bf16 / geneformer の場所 / モデル重み（Mac と同一か）/ bundle の sha256 を表示して、学習の設定を確認する（学習はしない）。
`NG` が出たら止めて Labretary に連絡。

## 3. 学習する（G2 承認後）
```bash
ssh petadimensionlab@100.104.171.109 'bash labretary-ft/bundle/remote/run_remote.sh probe'                 # 20 step で速度とメモリを測る（数分）
ssh petadimensionlab@100.104.171.109 'tmux new -d -s ft "bash labretary-ft/bundle/remote/run_remote.sh train"'   # 本番（切断しても続く）
ssh petadimensionlab@100.104.171.109 'tail -n 5 labretary-ft/out/train.log'                                # 進み具合
```
中断したら同じ `train` のコマンドをもう一度実行すると、最新の checkpoint から再開する。

## 4. 結果を戻す
```bash
rsync -av --exclude 'checkpoints/' petadimensionlab@100.104.171.109:labretary-ft/out/ft-main/ experiments/EXP-20260926-003/models/ft-main/
.venv/bin/python experiments/EXP-20260926-003/remote/verify_sha256.py check experiments/EXP-20260926-003/models/ft-main/final
```
確認できたら、リモートの `~/labretary-ft/` は削除してよい（ヒト由来データの複製を残さない）。
その後の評価（ゲート）と ISP は Mac で行う（plan.md §6.1）。

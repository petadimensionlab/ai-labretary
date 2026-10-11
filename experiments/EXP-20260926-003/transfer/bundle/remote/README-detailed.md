# EXP-20260926-003 リモート GPU での fine-tuning 手順（Tailscale 経由）

PI 指示（2026-09-26）: **fine-tuning だけ**を別の NVIDIA GPU マシンで行う。ISP（phase 3）とゲート評価（phase 2）は Mac（MPS）で行う。
**G2 承認前に手順 4 以降を実行しないこと**（リモートの GPU 時間も G2 の対象。plan.md §8, §10）。

- 送るもの: `transfer/bundle/`（学習・検証 donor の tokenized だけ, 約 19 MiB）。Python 環境・Geneformer・重みはリモートにセットアップ済みのものを使う
- 送らないもの: held-out donor の細胞（`work/ft/test.dataset`, `work/tokenized.dataset`, `work/input.h5ad`）。prepare が bundle を読み直して漏洩がないことを確認済み（`split_manifest.json` の `heldout` は件数のみで ID もない）
- 返すもの: `out/ft-main/`（`final/` = 最良モデル + SHA256SUMS、`env.json`, `run_record.json`, `train_log.jsonl`, `training_summary.json`, `train.log`）。checkpoint は返さなくてよい

以下の `<REMOTE_USER>`, `<REMOTE_HOST>`（Tailscale の MagicDNS 名または 100.x アドレス）は置き換える。作業ディレクトリはリモートの `~/labretary-ft`。
パスワード・トークン・鍵をこのファイルや config に書かない（SSH 鍵 / Tailscale SSH を使う）。

## 0. リモートの確認（最初に 1 回）
```bash
ssh <REMOTE_USER>@<REMOTE_HOST> 'uname -m; cat /etc/os-release | head -2; nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv; df -h ~ | tail -1'
```
- `uname -m` が `aarch64` でも手順は同じ（torch の CUDA wheel が aarch64 用になるだけ）。WSL2 は Linux と同じ手順（Windows 側に NVIDIA の WSL 対応ドライバが必要）
- 空き容量 ≥ 12 GB（checkpoint ~3.8 GB × 3）
- 結果を plan.md §8 の表（GPU 名・VRAM・driver）に記録する

## 1. ローカル（Mac）: bundle を作る・照合する
```bash
cd ~/workspace/research/ai-labretary
.venv/bin/python scripts/geneformer/finetune_classifier.py --config experiments/EXP-20260926-003/config-finetune.yaml --stage prepare --force
.venv/bin/python experiments/EXP-20260926-003/remote/verify_sha256.py check experiments/EXP-20260926-003/transfer/bundle --no-extra
```
（2026-09-26 に作成済み。config-finetune.yaml を変えたら作り直す）

## 2. ローカル → リモート: bundle
リモートの Python 環境（uv）と Geneformer（`~/workspace/research/Geneformer`, venv・Geneformer-V2-316M を含む）は**セットアップ済み**なので、送るのは bundle だけ。
```bash
rsync -av experiments/EXP-20260926-003/transfer/bundle/ <REMOTE_USER>@<REMOTE_HOST>:labretary-ft/bundle/
```

## 3. リモート: 照合と dry-run
```bash
cd ~/labretary-ft
bash bundle/remote/run_remote.sh check          # 環境（torch / CUDA / bf16 / geneformer）と重みの sha256（MODEL_SHA256）も照合                              # 重みの sha256（MODEL_SHA256）も照合
```
`env.json` と同じ情報（GPU 名, driver, CUDA, 各パッケージ版, geneformer の .py digest, config sha256, seed）は train 開始時に出力先へ書かれる。
geneformer の digest（`env.json` → `geneformer.py_digest`）が Mac の値 `c050bc47…` と同じならパッチ込みで同一。

## 4. リモート: probe（G2 承認後）→ 本番
VRAM に合わせて profile を選ぶ（実効バッチ 16 = micro-batch × 勾配累積 はどれも同じ。結果が同等になるようにしている）:

| VRAM | 例 | profile | micro-batch × 累積 | 備考 |
|---|---|---|---|---|
| 16 GB | RTX 4080 / A4000 / T4 級 | `vram16` | 2 × 8 | T4 は bf16 非対応 → 停止（torch.cuda.is_bf16_supported() が False なら停止） |
| 24 GB | RTX 4090 / 3090 / L4 / A10 | `vram24` | 4 × 4 | |
| 40–48 GB | A100-40GB / L40S / RTX 6000 Ada / A6000 | `vram48` | 8 × 2 | |
| 80 GB 以上 | A100-80GB / H100 / GH200 / 128 GB 統合メモリ機 | `vram80` | 16 × 1 | |

- いずれも gradient checkpointing・bf16。最長 4,096 トークンの細胞が同じバッチに入った場合を想定した保守的な値
- **まず probe**（20 optimizer step, 数分）: ピークメモリ（`probe_summary.json` の `cuda_peak_mem_gb`）が VRAM の 85% を超えたら 1 段小さい profile、40% 未満なら 1 段大きい profile にしてよい（実効バッチは変わらない）
```bash
bash bundle/remote/run_remote.sh probe [PROFILE]  # 既定 vram80。→ ~/labretary-ft/out/ft-main-probe/probe_summary.json（s/step と 1 epoch の見積もり）
bash bundle/remote/run_remote.sh train [PROFILE]  # 本番。tmux / screen の中で実行する
```
- 中断（接続切れ・OOM・再起動）したら**同じコマンドをもう一度**実行 → 最新の step checkpoint から再開（config の sha256 が変わっていたら停止）
- 途中で profile を変えても再開できる（実効バッチが同じなので学習率スケジュールは変わらない。変更は `run_record.json` に警告として残る）
- 終了時に `out/ft-main/final/SHA256SUMS` が作られ、照合まで自動で行う

## 5. リモート → ローカル: 結果を戻す
```bash
# ローカル（Mac）で実行。checkpoint（数 GB）は除く
rsync -av --exclude 'checkpoints/' \
    <REMOTE_USER>@<REMOTE_HOST>:labretary-ft/out/ft-main/ experiments/EXP-20260926-003/models/ft-main/
rsync -av <REMOTE_USER>@<REMOTE_HOST>:labretary-ft/out/ft-main-probe/probe_summary.json experiments/EXP-20260926-003/models/ft-main/ 2>/dev/null || true
.venv/bin/python experiments/EXP-20260926-003/remote/verify_sha256.py check experiments/EXP-20260926-003/models/ft-main/final
ls -la experiments/EXP-20260926-003/models/ft-main experiments/EXP-20260926-003/models/ft-main/final
```
- `final/split_manifest.json` の `config_sha256` と `final/config-finetune.yaml` がローカルの config と同じことを確認
- リモートの `bundle/` と `out/` は結果の確認後に削除してよい（ヒト由来データの複製を残さない。KPMP の利用条件）

## 6. ローカル: ゲート評価（phase 2, MPS）
MPS で他のジョブが動いていないことを確認してから:
```bash
.venv/bin/python scripts/geneformer/finetune_classifier.py --config experiments/EXP-20260926-003/config-finetune.yaml --stage evaluate
# → results/finetune-eval/gate.json（go / no-go）, test_predictions.csv
```
`go: false` なら phase 3（FT モデルでの ISP）は行わない（plan.md §6.1）。

## リモートが使えない場合（MPS フォールバック）
```bash
.venv/bin/python scripts/geneformer/finetune_classifier.py --config experiments/EXP-20260926-003/config-finetune.yaml --stage train --profile mps_fallback
```
学習細胞を 2,091 に減らし、下位 12 層と埋め込みを凍結、1 epoch（見積もり 0.8–1.6 h。plan.md §8）。精度が下がる可能性があるため、ゲートの結果と合わせて報告する。

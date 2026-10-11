# ai-labretary — PetadimensionLab の「成長する AI 研究秘書」

[Claude Code で「成長する AI 秘書」を作ってみた（DevelopersIO）](https://dev.classmethod.jp/articles/ai-secretary/) の構成
（CLAUDE.md で人格定義 → サブエージェント委託 → 日中は指示・夜間は自律 → 5 トリガーの学習ループ）を、
**in silico 研究（Geneformer による perturbation 実験）** 向けに作り直したものです。

| 元記事 | このリポジトリ |
|---|---|
| Nina | **Labretary** |
| kamedas（事業企画） | **PetadimensionLab**（研究） |
| 仮説 → 顧客マッチング → 検証 → サービス化 | 仮説 → 実験計画 → 数値実験 → 検証 → 実験レポート |
| researcher / writer / reporter | literature-scout / hypothesis-builder / experiment-planner / simulation-engineer / validator / report-writer / lab-reporter |
| /morning-news, /nina-daily-report | /morning-papers, /hypothesize, /plan-experiment, /run-perturbation, /validate-experiment, /experiment-report, /lab-status, /labretary-daily-report |
| nina-level.md / boss-patterns.md / guardrails.md | labretary-level.md / pi-patterns.md / guardrails.md |
| nina-tmux.sh | **labretary-ghostty.sh**（Ghostty） |
| nina-nightshift.sh | labretary-nightshift.sh |

## 構成
> **2026-09-27 訂正（実態との違い）**: 下の図の `.claude/settings.json`・`.claude/skills/`・`.claude/hooks/`・`.mcp.json` は、このプロジェクトには**存在しない**。スラッシュコマンドは使えないので、各担当（`.claude/agents/*.md`）に委託する。data/raw の書き込み禁止や危険コマンドのブロックは hook では強制されておらず、CLAUDE.md のルールで守っている。夜間の大規模ジョブの停止は `run_isp.py` の検査だけで強制している。

```
ai-labretary/
├── CLAUDE.md                        ← Labretary の人格・判断基準・ゲート・Learning Protocol
├── .claude/
│   ├── settings.json                ← 権限（data/raw は書き込み禁止）と hooks 登録
│   ├── agents/                      ← 7 つのサブエージェント（計画と検証は opus、日報は haiku）
│   ├── skills/                      ← 8 つのワークフロー（/コマンドで呼べる）
│   ├── hooks/                       ← 危険コマンド・資格情報・夜間大規模ジョブのブロック / 行動ログ / 終了マーカー
│   ├── memory/                      ← 成長記録・PI の判断パターン・ガードレール
│   └── rules/                       ← Geneformer のドメイン知識・実験標準
├── .mcp.json                        ← PubMed / bioRxiv / ChEMBL / Clinical Trials MCP（任意）
├── hypotheses/tracker.md            ← 仮説ストックと結果（支持率＝成長指標）
├── experiments/EXP-*/               ← 1 実験 1 ディレクトリ（plan → config → results → validation → report）
├── templates/                       ← 仮説・計画・config・レポートの雛形
├── scripts/
│   ├── geneformer/run_isp.py        ← config 駆動の perturbation 実行（dry-run / smoke / mock）
│   ├── analysis/validate_isp.py     ← コントロール・null・seed 間再現性・発現量交絡の定量チェック
│   ├── analysis/plot_isp.py         ← 標準図
│   ├── capture_env.sh               ← 再現情報（git, パッケージ版, GPU, config ハッシュ）
│   └── lab_status.sh, watch_*.sh    ← ダッシュボード用
├── tasks/night-queue/               ← 夜間キュー（1 タスク 1 ファイル）
├── labretary-nightshift.sh          ← 夜間自律シフト
├── labretary-ghostty.sh             ← Ghostty ダッシュボード
├── ghostty/labretary.conf           ← Ghostty 設定（分割 keybind 等）
└── scheduler/                       ← launchd plist / crontab
```

## セットアップ
```bash
# 1. 依存
brew install jq git-lfs          # Linux は apt 等で
# Ghostty: https://ghostty.org/  / Claude Code: https://docs.claude.com/en/docs/claude-code

# 2. 既存の Geneformer（既定 ../Geneformer、GENEFORMER_HOME で変更可）を .venv と models/Geneformer にリンク
make setup

# 3. GPU なしでパイプライン確認（合成データ）
make demo
# ※ Geneformer は Apple Silicon の MPS でも動作する（config.yaml の compute.device: mps）

# 4. ダッシュボードを開く
make ghostty           # 役割ごとのウィンドウ（macOS / Linux）
make ghostty-splits    # 1 ウィンドウ 4 分割（macOS、システム設定でアクセシビリティ許可が必要）
```
最初に `.claude/rules/domain-knowledge.md` の「現在の研究テーマ」を埋めてください。Labretary はここを起点に文献・仮説を探します。

## 日中の使い方
Ghostty の「Labretary」ウィンドウ（`claude`）で話しかけます。
```
/morning-papers                          朝の新着文献 → 仮説のタネ
/hypothesize 拡張型心筋症の心筋細胞を非疾患状態へ戻す遺伝子
/plan-experiment HYP-20260925-01
/run-perturbation EXP-20260925-001 --smoke
/validate-experiment EXP-20260925-001
/experiment-report EXP-20260925-001
/lab-status
```
PI の判断が必要な地点は 3 つだけ（G1 仮説採択 / G2 大規模ジョブ承認 / G3 結論確定）。それ以外は Labretary が自明な次ステップとして進めます。

## 夜間シフト
`tasks/night-queue/` にタスクを置くと、01:00 起動・05:45 までに 1 件ずつ消化し、最後に学習セッションで learnings をメモリへ転記し、朝読む要約を `outputs/reports/daily/*-night.md` に書きます。
- 大規模な Geneformer ジョブは、タスクの frontmatter で PI が `approved: true` にしたものだけ実行（`run_isp.py` の検査で強制。hook は未設定）
- `SOURCE=github` で元記事と同じく GitHub Issue（label `labretary-night`）を消化
- スケジュールは `scheduler/` を参照。macOS は `caffeinate`、Linux は `systemd-inhibit` でスリープ抑止

## 成長の仕組み
元記事と同じく「作業の正確性」ではなく「判断の質」で測ります（`.claude/memory/labretary-level.md`）。
業務別（文献・仮説・計画・実行・検証・レポート）に、PI の「修正不要」3 連続で昇格。
Lv3 になった業務だけ夜間に自律実行を任せる、といった運用に使えます。指標は仮説の採択率・支持率、計画の差し戻し率、レポートの修正往復数、再実行率。

## 元記事からの発展点
- **ゲート制**: 研究では「勝手に GPU を 20 時間回す」ことが最大の事故なので、承認点を明示
- **検証専任エージェント（懐疑的査読者）**: ポジティブ／ネガティブコントロール、発現量マッチ null、seed・サブサンプル・層の頑健性、発現量交絡、データリークを必須チェック化
- **再現性の自動記録**: config ハッシュ、環境、全 run（失敗 run 含む）を保存
- **事前登録**: 判定基準を plan.md に先に書き、結果を見てから基準を変えない
- **版差に強いランナー**: Geneformer の各クラスの signature を見て受け付けられる引数だけ渡す

## 今後の拡張アイデア
- scGPT / GEARS / UCE など他モデルのランナーを `scripts/<model>/` に追加し、同じ summary.csv 形式で横比較
- 公開 Perturb-seq（Replogle 2022 など）を使った「予測 vs 実測」のベンチマーク実験をテンプレート化
- 組み合わせ perturbation（`combos: 1`）と anchor gene を使った相乗効果探索
- wet ラボ向けに、上位候補から sgRNA 設計リストや実験依頼書を作る `wet-handoff` エージェント
- 週報に仮説の支持率の推移グラフを自動添付

## 注意
- `run_isp.py` の Geneformer 呼び出しは公開 API に沿っていますが、手元の版で `make setup` の signature 出力を確認し、必要なら調整してください（最初の run は必ず `--smoke`）
- demo config の Ensembl ID は例です。コントロール遺伝子は literature-scout の根拠付きで選び直してください
- `--splits` はキーストローク送信で分割を作るため、Ghostty の keybind 設定やタイミングに依存します。うまくいかない場合はウィンドウ方式を使ってください
- in silico perturbation の結果はモデルが学習した共変動に基づく予測で、因果の証明ではありません

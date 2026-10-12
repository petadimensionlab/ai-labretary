# AGENTS.md — OpenCode 用の適応メモ

このプロジェクト（PetadimensionLab の AI 研究秘書 Labretary）の研究ルールは
`CLAUDE.md` にある。OpenCode では `opencode.json` の `instructions` で
`CLAUDE.md` を読み込むので、内容は同一。**規則の本体は `CLAUDE.md`。**
以下は Claude Code と OpenCode の差だけを書く。

## サブエージェントの場所

- OpenCode では `.opencode/agent/<name>.md` を読む。**`.claude/agents/*.md` は読まない。**
- 10 体を同じ役割・同じ本文で `.opencode/agent/` に置いた。名前は Claude Code と同一。
  `Task` ツールの `subagent_type` にこの名前を渡す。
- Claude Code の `tools:` は OpenCode の `permission:` に読み替えた
  （読み取り専用の担当は `edit: deny`）。

## モデルの対応

Claude Code の 3 階層（当時）を `opencode-go` のモデルに対応付けた。

すべての担当を同じモデルに統一した（2026-10-11 変更）。

| 担当 | Claude | OpenCode |
|---|---|---|
| 全 10 体 | opus / sonnet / haiku | `opencode-go/deepseek-v4-flash`（variant `high`） |

Claude Code の 3 階層（opus / sonnet / haiku）は使わず、
`opencode-go/deepseek-v4-flash` を全員に割り当てている。
variant は推論強度。

モデルを変えるときは `.opencode/agent/<name>.md` の `model:` を書き換え、
OpenCode を再起動する（設定は起動時に一度だけ読まれる）。

## ツールの差（重要）

- **`Artifact` ツールは OpenCode に無い。** `artifact-refresher` の
  「同じ URL に再公開する」手順は実行できない。OpenCode では
  `scripts/build_results_page.py` による `outputs/results-page.html` の生成までで止め、
  再公開は PI に確認する。
- `WebSearch` / `WebFetch` は OpenCode の `websearch` / `webfetch` に対応する。
- `Task` は OpenCode の同名ツールに対応する。

## 迷ったとき

- ゲート（G1 仮説採択 / G2 計画承認 / G3 結論確定）は PI だけが判断する。
- 盲検中の遺伝子（TP53, CDKN1A, LTBR, SPP1, LCN2）は解除まで触らない。
- `data/raw/` と `models/` は読み取り専用（`opencode.json` の `permission.edit` でも拒否）。

## 技能（skill）の場所（2026-10-12 追加）

- OpenCode では `.opencode/skills/<name>/SKILL.md` を読む。公式の Skill ツールで
  ロードできる（`skill` ツールの `name` に `honest-reporting` / `evidence-before-revert`）。
- Claude Code では `.claude/skills/` 相当だが、このプロジェクトは **`.opencode/skills/` のみを正とする**
  （CLAUDE.md §10 の要点がそれぞれの SKILL.md にある）。
- 2 つの技能 = 正直な報告（honest-reporting）／証拠を見るまで引き返さない（evidence-before-revert）。
  報告や検証のときはこの 2 つを効かせる。

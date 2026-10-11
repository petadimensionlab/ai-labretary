# 夜間キュー
1 タスク = 1 ファイル（`NNN-<slug>.md`）。先頭の番号順に処理され、完了後は `tasks/done/` へ移動する。

書式:
```
---
type: literature | hypothesize | plan | run | validate | report | custom
target: EXP-... / HYP-... / テーマ
approved: false   # run で大規模ジョブを許可する場合のみ true（PI が書く）
max_turns: 40
---
やってほしいことを自然文で
```

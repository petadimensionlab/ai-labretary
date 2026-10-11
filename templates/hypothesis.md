# HYP-YYYYMMDD-NN: <一文の仮説>

- Status: proposed | accepted | planned | ran | supported | refuted | inconclusive | rejected
- 起点: <文献メモ / 過去 EXP / PI の着想>
- 新規性: 既知の再確認 | 新規予測

## 仮説
<細胞系> において <遺伝子/遺伝子セット> を <delete/overexpress/inhibit/activate> すると、細胞状態が <start> から <goal> へ近づく。

## 予測（定量）
- 例: 対象遺伝子の goal_state_shift が発現量マッチしたランダム遺伝子 null の上位 5% に入る（FDR < 0.05）

## 反証条件
- 例: 3 seed すべてで null と区別できない / シフト方向が逆

## コントロール
- ポジティブ: <既知ドライバーと根拠>
- ネガティブ: <ランダム遺伝子 n, ハウスキーピング遺伝子>

## 実験コスト概算
- 細胞数 × 遺伝子数 × 条件数 → GPU 時間 ~ <h>

## 支持された場合の次の一手（wet 実験）
-

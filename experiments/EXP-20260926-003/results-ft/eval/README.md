# EXP-20260926-003 phase 3–4 の結果（FT モデル vs Pretrained, 同じ細胞・同じ遺伝子・同じ goal セントロイド・同じマシン）

- 実行: spark-cf5b（NVIDIA GB10, CUDA, bf16 forward + fp32 コサイン）。FT 約 44 min、Pretrained（config-pretrained-matched）約 46 min
- 対照のみ 34 遺伝子 × 200 細胞 × seed 3/4/5（EXP-002 案 A と同じ細胞・同じ null）。goal セントロイドは FT に使っていない健常 9 donor
- 評価: EXP-002 `analysis/evaluate_controls.py`（criteria.csv）。判定規則は plan.md §6.2 / HYP-20260926-02（事前登録）

## P1（主, a2）: C1・C2・C3・C3-panel をすべて満たせば支持
| 基準 | FT | Pretrained（同一条件） |
|---|---|---|
| C1 陽性が null 上位 10%（3 seed 中 2） | VCAM1 1.00/1.00/0.95 ✓、**HAVCR1 0.45/0.75/0.75 ✗** → **Fail** | VCAM1 0.617 平均 / HAVCR1 0.80 → Fail |
| C2 陽性 > HK（18 組中 15） | **16/18 ✓** | 10/18 ✗ |
| C3 HK が null 分位 0.10–0.90 | ACTB 0.38 / GAPDH 0.65 / PPIA 0.20 ✓ | ACTB 0.17 / GAPDH 0.80 / PPIA 0.83 ✓ |
| C3-panel（7 中 5 以上） | 7/7 ✓ | 7/7 ✓ |
| C4（参考）split-half | 0.46 ✗ | 0.32 ✗ |
| C5（参考）長さとの相関 | 0.11 ✓ | −0.03 ✓ |
→ **P1 は事前登録の規則で棄却**（C1 の HAVCR1 が不合格）。fp32 再計算の対象（分位 0.85–0.95）は VCAM1 seed 5（0.95）のみで、VCAM1 は seed 3・4 で合格、HAVCR1 は境界外のため、判定は変わらない（再計算は未実施）

## P2（副）: 削除ドリフト
| | FT | Pretrained |
|---|---|---|
| null 遺伝子の補正前 shift_goal が負の割合（seed 3/4/5） | **0.10 / 0.15 / 0.15** | 1.00 / 1.00 / 1.00 |
| null 平均 shift_goal ÷ (1 − cos(start, goal)) | **≈ 0.000** | −0.036 / −0.038 / −0.036 |
→ **支持**: FT で削除ドリフトはほぼ消えた

## P3（副）: C2 の勝ち数 FT ≥ Pretrained
a2 で 16 vs 10、補正前で 17 vs 6 → **支持**

## 参考: 分類ヘッドの ΔP(PT_healthy)（a2 と同じ長さビン補正）
- VCAM1 分位 0.75/0.90/0.60、HAVCR1 0.35/0.50/0.15、HK 0.0–0.8、パネル 6–7/7 が範囲内、null 平均 ΔP ≈ 0
- 陽性対照は分類確率でも null から際立たない（埋め込み指標と一致）

## 解釈（示唆）
- fine-tuning は削除ドリフトと HK の偽陽性（Pretrained で GAPDH/PPIA が null 最上位）を解消した。方法上の大きな改善
- それでも HAVCR1 の delete は健常方向への特異的な移動を示さず、VCAM1 だけが一貫して上位。VCAM1 は分類器が直接使うマーカーの可能性があり、C1 合格でも駆動因子の証拠にはならない（HYP §循環性）
- split-half は 0.46 で、donor 間の違いによる不安定さは残る

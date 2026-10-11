# (c-1) donor 単位の再解析（EXP-003 の結果を CPU で再集計。GPU 不使用）

- 事前登録: [../plan.md](../plan.md) §C.1（定義は対照で計算する前に固定）。実装 `../analysis/{common.py,donor_level.py}`
- 対象: `experiments/EXP-20260926-003/results-ft`（FT）と `results-pretrained-matched`（Pretrained, 同じ細胞）。TLS-high aPT 14 donor（うち a2 の集計に使える donor は遺伝子×seed あたり中央値 13、最小 11。細胞 3 未満の donor は除く）
- 妥当性の確認: 本スクリプトの a2 は EXP-003 の `gene_seed_metrics.csv` と最大差 9.6×10⁻¹⁷（= 一致）

## 判定（主指標 a2, FT）
| 集約 | C1（VCAM1 / HAVCR1 の null 分位, seed 3/4/5） | C2（陽性>HK） | C3（HK が 0.10–0.90） | C3-panel |
|---|---|---|---|---|
| 細胞平均（EXP-003 と同じ） | ✗ VCAM1 1.00/1.00/0.95、HAVCR1 0.45/0.75/0.75 | 16/18 ✓ | ✓ | 7/7 ✓ |
| **DL-med（donor 中央値, 主）** | ✗ VCAM1 1.00/1.00/0.95、HAVCR1 0.10/0.85/0.70 | 14/18 ✗ | ✓ | 7/7 ✓ |
| DL-mean（donor 等重み平均） | ✗ VCAM1 1.00/1.00/1.00、HAVCR1 0.25/0.65/0.85 | 16/18 ✓ | ✓ | 7/7 ✓ |

- **結論は変わらない**: どの集約でも C1 は HAVCR1 で不合格（EXP-003 の P1 棄却は donor の重みづけによるものではない）。VCAM1 は donor 単位でも一貫して null 最上位
- DL-med では C2 が 14/18（基準 15 未満）に下がる。内訳は HAVCR1 vs HK の組で、HAVCR1 の弱さが donor 中央値でより明確になったため
- Pretrained（同じ細胞）: DL-med で C1 ✗（VCAM1 0.55/0.75/0.80、HAVCR1 0.95/0.85/0.75）、C2 13/18、C3 ✗（PPIA 0.95/0.95/0.60 が範囲外）、panel 7/7 → FT が HK の偽陽性を解消した所見は donor 単位でも保たれる

## donor 単位の符号の一貫性（a2, FT, 13 donor）
| 遺伝子 | 値が正の donor の割合 | null 遺伝子の同割合の分布での分位 | null 超過分（donor ごとに同 seed・同 donor の null 平均を引いた値）が正の donor | Wilcoxon（片側） |
|---|---|---|---|---|
| VCAM1 | **13/13 (1.00)** | 0.98 | 13/13、中央値 +1.6×10⁻³ | p = 1.2×10⁻⁴ |
| HAVCR1 | 5/13 (0.385) | 0.30 | 7/13、中央値 +8.3×10⁻⁵ | p = 0.58 |
| GAPDH | 0.667 | 0.88 | 8/12 | p = 0.17 |
| ACTB / PPIA | 0.46 / 0.23 | 0.48 / 0.08 | 4/13 / 2/13 | — |
（null 遺伝子の「正の donor 割合」の中央値 0.46）
→ VCAM1 の効果は 13 donor すべてで同じ向き。HAVCR1 は donor 間で符号が割れる（**donor 間の不均一が HAVCR1 の不合格の一因**）

## 信頼性
| 指標 | 細胞平均 | DL-med | DL-mean |
|---|---|---|---|
| donor split-half の Spearman ρ（seed ごと 20 分割 × 3 seed の平均 ± SD, 34 遺伝子） | 0.464 ± 0.138 | 0.423 ± 0.131 | 0.442 ± 0.118 |
| Spearman–Brown 補正（全 donor 相当） | 0.634 | 0.595 | 0.613 |
| Pretrained の同じ値（ρ） | 0.320 | 0.272 | 0.216 |
| seed 間の順位一致（固定対照 5 遺伝子, 参考 C4ii） | 0.70（FT）/ 0.40（Pretrained） | — | — |

- **donor 単位に集約しても再現性は上がらない**（0.464 → 0.423）。split-half が低いのは donor の重みづけではなく、34 遺伝子の多く（null と陰性対照）が実質ゼロ効果で順位がノイズだからと**示唆される**（Spearman–Brown 補正後でも 0.6 前後で基準 0.7 未満）
- leave-one-donor-out（DL-med, a2）: 全 donor との順位相関は ρ 最小 0.943・中央値 0.969。C1（VCAM1 3 seed 合格・HAVCR1 0 seed）と C3-panel はどの donor を除いても変わらない。C3 は 14 回中 3 回で反転（HK が範囲外に出る境界例）。C2 も 14/18 → 13〜14/18 の範囲で基準未満のまま
→ **判定は特定の donor に依存していない**

## 副指標
- raw shift_goal: DL-med でも VCAM1 1.00/1.00/1.00、HAVCR1 0.10/0.85/0.65（a2 と同じ傾向）
- ΔP(PT_healthy): DL-med で C3 不合格（ACTB 1.00/0.60/0.95、GAPDH 0.80/0.45/0.95）、null の負の割合 0.95 → 分類確率は delete で一様に下がる（埋め込みの削除ドリフトに相当する偏りが確率側に残る）。ΔP を主指標にはできない

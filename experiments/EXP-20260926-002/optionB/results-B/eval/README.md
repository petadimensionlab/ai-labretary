# 案 B の結果（対照のみ, 100 細胞 × seed 1/2 + phase 0 の seed 0）

- 評価: `analysis/evaluate_controls.py --results phase0/results-phase0 optionB/results-B`（criteria.csv）
- **再現性**: 案 B は EXP-001 run 1 と同じ細胞・同じ null・同じ seed のため、shift_goal は EXP-001 と 5,000/5,000 行で完全一致（最大差 0.0）。したがってコサイン系指標（raw, d, a2 など）の判定は EXP-001 の再解析と**同一**で、独立な追加証拠にはならない
- 新しい情報はベクトル系指標（V0–V5）の 3 seed 評価のみ: どれも C2（5–7/18）と C3（GAPDH/PPIA が null 分位 0.97–1.00）を満たさない
- a2（推奨指標）: C3・C5 合格、C1（VCAM1 分位 0.70）・C2（13/18 < 15）・C4（split-half ρ 0.377）不合格
- seed 1 の null が `screen` と記録されていた（seed を別プロセスで回したため）。評価では negatives.json で戻した。run_isp.py の summarize も修正済み

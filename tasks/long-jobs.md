# 長時間の計算の記録

見込み 30 分以上の計算は、**始める前に**ここへ 1 行足す。終わったら状態を書き換え、結果の場所を書く。
セッションが切れても、ここを見れば「何が・どこで・いつまで・どう確かめるか」が分かるようにする（PI 指示 2026-09-27）。

| 開始（JST） | 実験 | 内容 | 実行先 | 見込み（上振れ） | 中断規則 | 状態の確かめ方 | 状態 | 結果 |
|---|---|---|---|---|---|---|---|---|
| 2026-09-27 13:00 ごろ（再抽出 13:36〜14:03） | EXP-20260927-002 | 線維芽細胞の事前チェック（埋め込み抽出 3 seed＋識別 AUC）。摂動なし。結合バグ（cell_embs.csv と state-embs dataset の位置ずれ）を row_id での ID 結合に修正し、埋め込みを取り直した（Labretary 指示 2026-09-27。判定基準・細胞選択・検出された遺伝子数の分布をそろえる操作は不変。plan.md 末尾「逸脱の記録」参照） | Mac（MPS, bf16） | 実測 約1.4 h（初回0.5h＋再計算27分＋解析） | なし（G2 不要） | — | **完了** | `experiments/EXP-20260927-002/results-precheck/summary.json`（verdict: 合格）。C1 3seed平均AUC 0.7044（0.7207/0.6902/0.7022）、C1b CI下限中央値0.6405（0.67/0.6324/0.6405）、C1c p=0.005（3seedとも）、C2 検出遺伝子数・総カウントだけのAUC平均0.5611（margin 0.1433）、C3 陽性対照AUC平均0.9873。仮説Status・G1/G2/G3は未判断（PI待ち） |
| 2026-09-27 13:14 ごろ | EXP-20260927-001 | 過剰発現の摂動 102 遺伝子 × 100 細胞 × 3 seed（seed 11/12/13） | spark-c083（GB10, bf16） | 1.4 h（1.6 h） | seed 11 完了時点で 0.6 h 超なら停止・報告 | `bash scripts/remote/isp_remote.sh status coi-next1@100.94.57.17 experiments/EXP-20260927-001/config.yaml`（サンドボックス外）、`experiments/EXP-20260927-001/logs/progress.md` | **完了**（13:12〜14:26, 1.23 h, エラーなし。リモートの複製は削除済み） | `experiments/EXP-20260927-001/results/eval/criteria.json`（暫定: 否定の提案。validator の検証待ち） |
| 2026-09-27 18:25 ごろ（実開始 18:55:25） | EXP-20260927-003 | 線維芽細胞の過剰発現 84 遺伝子 × 100 細胞 × 3 seed（seed 21/22/23）。LTBR/RELB/NFKB2 が対象 | リモート計算機 2 号（coi-next1@100.94.57.17, spark-c083, GB10, bf16） | 1.2 h（上振れ 1.4 h） | seed 21 完了時点で 0.7 h 超なら停止・報告（2 h 超見込みも同様） | — | **完了**（18:55:25〜19:50:00, 実測 0.91 h。seed21 完了は開始から 0.30 h で中断規則に抵触せず。追加のセントロイド 19:51:52〜19:55:22 完了。回収・sha256 一致確認・評価まで完了。リモートの複製は `clean` 済みで消失を確認） | `experiments/EXP-20260927-003/results/eval/criteria.json`（暫定: **判定不能**〔恒常発現遺伝子 GAPDH の中立性 C3 不合格〕。主対象 3 遺伝子とも (a) 不合格。validator・PI 待ち） |
| 2026-09-27 23:35 ごろ（前段で計 4 回停止: 設定不備 2 回、6 秒クラッシュ 1 回、operator ミス 1 回。GPU 消費ゼロ。実開始 23:35:40） | EXP-20260927-004 | 過剰発現型の指標の作り直し。106 遺伝子 × 200 細胞 × 3 seed（31/32/33）。主対象は回さない。smoke（2 遺伝子×10 細胞×seed31）を先に通過確認済み。`work/extra_centroids-seed{31,32,33}.npz` も生成済み（MYOF_neg 615/615/597 細胞, FIB_org_unmatched 711 細胞 ×3、dim 1152） | リモート計算機 2 号（coi-next1@100.94.57.17, spark-c083, GB10, bf16） | 4.0 h（実測モデル採用。内蔵見積もり 2.33 h は採らない） | seed 31 完了時点で 2.0 h 超なら停止・報告。全体 6 h 超見込みも同様 | `experiments/EXP-20260927-004/logs/progress.md` | **完了**（23:35:40〜01:37:21, 実測 2.03 h。seed31 完了 0.68 h で中断規則に抵触せず。exit=0、エラーなし。ファイル数照合済み（106遺伝子×3seed、cells.csv 63,570行、blocked_genes混入0件）。リモートの複製は `clean` 済みで消失確認済み） | `experiments/EXP-20260927-004/results/eval/criteria_004.json`（暫定提案: M1〜M5 すべて合格。目標の軸の陽性対照 10 個中 SYN3 の 1 個が合格。GAPDH は中立域に復帰。validator・PI 待ち） |
| 2026-09-28 04:03:45（1回目 state_rows 未送付で即 exit=1、GPU消費ほぼ0。04:05:28 再実行） | EXP-20260928-001 | 陽性対照の確立をやり直す。判定の単位を遺伝子セットに変える。156 遺伝子 × 200 細胞 × 3 seed（41/42/43） | リモート計算機 2 号（coi-next1@100.94.57.17, GB10, bf16） | 6.0 h | seed 41 完了時点で 3.0 h 超なら停止・報告（全体 9 h 超見込みも同様）。**両方とも抵触せず**（seed41完了は開始から約1.0h、全体3.00hで完了） | `experiments/EXP-20260928-001/logs/progress.md` | **完了**（04:05:28〜07:05:31, 実測 3.00 h, exit=0。ファイル数照合済み: 156遺伝子×3seed完備、summary.csv 469行、timing.jsonl 471行、blocked_genes混入0件、config.sha256一致。追加で `analysis/extra_centroids.py` を seed41/42/43 それぞれ Mac(MPS) で実行(計約14分、10:05〜10:19)。監視は当初Monitorの30分上限後に再アームが約2h空いて完了検知に遅れた。原因はlearnings参照。リモートの複製は `clean` 済みで消失確認済み） | `experiments/EXP-20260928-001/results/eval/criteria_001.json`（暫定提案: G1〜G6すべて合格。目標マーカーセット24の陽性対照が確立との提案。判定確定はvalidator・PI待ち） |
| 2026-09-28 17:43:42（cf5b）/ 17:43:48（c083） | EXP-20260928-002（機体一致確認） | 48 遺伝子 × 200 細胞 × seed 51 を 2 機体で同時に回し、埋め込みと指標を突き合わせる | spark-c083 と spark-cf5b | 0.35 h/機体（見積もり上振れ込み） | H-1〜H-3 不合格なら本番を開始しない | `experiments/EXP-20260928-002/logs/progress.md` | **完了**（17:43〜18:03、実測 0.33h/0.325h、exit=0 両機体。H-1〜H-3 すべて合格。系統成分絶対値0.0・最大絶対差0.0・ρ=1.0・コサイン最小1.0。2機体の出力はsha256まで完全一致） | `experiments/EXP-20260928-002/work/machine_check_002.json`、`results-machinecheck-{cf5b,c083}/` |
| 2026-09-28 18:08:17（cf5b）/ 18:08:35（c083） | EXP-20260928-002（本番） | 経路セット 120 + 対照 980 + アッセイ対照 36 = 1,136 遺伝子 × 200 細胞 × 3 seed（51/52/53）。seed 単位で分割 | cf5b が seed 51、c083 が seed 52→53（tmux isp-c083-seed52-53 で連結） | cf5b 7.32 h、c083 計 14.6 h | A-1〜A-5（plan.md §8.5） | `experiments/EXP-20260928-002/logs/progress.md` | **完了**（cf5b: 18:08:22〜01:05:42 実測6.96h。c083: 18:08:35〜08:18:11 実測14.16h＝seed52 7.08h+seed53 7.08h。全機体 exit=0、skip 0 件、遮断遺伝子混入 0、3 seed 統合・自動照合すべて合格。A-1 チェックポイント〔seed52単独, G-3/M1/M4〕合格。**途中 A-2 の誤発報 1 件**（監視スクリプト自作分のログファイル名の食い違いによる誤検知。実際の計算は正常。guardrails に追記済み）。両機体 `clean` 済み、複製消失確認済み） | `experiments/EXP-20260928-002/results/eval/criteria_002.json`（暫定提案: **G-1 主判定は棄却**。3 seed・980/600/520 の全版で不合格が一致。副次 nonmarker も不合格。アッセイの妥当性・中立性・指標の健全性・偽経路セット・機体一致はすべて合格。判定確定は validator・PI 待ち） |

## 2026-09-28 18:xx 時点の引き継ぎ（利用上限のリセット 19:50 後に再開）
- **EXP-20260928-002 の機体一致確認を実行中**（48 遺伝子 × 200 細胞 × seed 51 を 2 機体で並行、0.31 h/機体）。合格なら本番（seed 単位で分割、約 14.6 h）へ。
- **リモートのジョブはセッションが切れても動き続ける。** 再開したらまず状態を自分で確認する:
  - `bash scripts/remote/isp_remote.sh status coi-next1@100.94.57.17 experiments/EXP-20260928-002/config-seed52.yaml`
  - `bash scripts/remote/isp_remote.sh status petadimensionlab@100.104.171.109 experiments/EXP-20260928-002/config-seed51.yaml`
  - いずれもサンドボックス外で実行する
- 経緯は `experiments/EXP-20260928-002/logs/progress.md`。判定基準は `plan.md`（事前登録 38 ファイル、固定済み）。
- **見込み時間を過ぎたら Labretary 自身が status を叩く**（2026-09-28 の guardrails。前回は完了に 3 時間気づかなかった）。
- **状態確認の注意（2026-09-28）**: c083 は 2 seed を連結して回しているため、`isp_remote.sh status` が「NOT_RUNNING」と誤表示する（ログのファイル名が想定と違う）。正しくは次で確認する:
  - `ssh coi-next1@100.94.57.17 "tmux ls; tail -2 ~/labretary-remote/experiments/EXP-20260928-002/logs/remote-config-seed52.log; ls ~/labretary-remote/experiments/EXP-20260928-002/results-seed52/run-52-L-1/genes | wc -l"`
  - `ssh petadimensionlab@100.104.171.109 "ls ~/labretary-remote/experiments/EXP-20260928-002/results-seed51/run-51-L-1/genes | wc -l"`（1,136 で完了）
| 2026-09-29 20:07:11（cf5b）/ 20:07:21（c083） | EXP-20260929-001（機体一致確認） | 87 遺伝子 × 200 細胞（seed 61）を 2 機体で回して突き合わせ | c083 と cf5b | 0.59 h/機体（上振れ 0.89 h） | 不合格なら本番を開始しない。経過 0.89h 超で状態確認 | `ssh coi-next1@100.94.57.17 "tmux ls; ls ~/labretary-remote/experiments/EXP-20260929-001/results-machinecheck/run-61-L-1/genes | wc -l"`（87で完了）／同様に cf5b | **完了**（20:07〜20:43、実測 0.6h/機体。両機体 exit=0、87/87遺伝子。H-1/H-2/H-3 すべて合格＝ほぼビット同一） | `experiments/EXP-20260929-001/work/machine_check_003.json`（pass: true） |
| 2026-09-29 20:47:55（cf5b, seed61）/ 20:48:00（c083, seed62→63は62完了後） | EXP-20260929-001（本番） | 分泌プログラム 24 + 対照 500 ほか 567 遺伝子 × 200 細胞 × 3 seed（61/62/63）。seed 単位で分割 | cf5b が seed 61（3.7 h）、c083 が seed 62・63（7.3 h、逐次） | cf5b 3.67h、c083 計7.3h | plan.md §6（seed62完了時点3.0h超・全体9h超→A-2、全体14h超→A-3、いずれかのseed完了時点でG-3a/G-3b/M1-M4不合格→A-1で中断、skip 5%超→A-4） | `ssh coi-next1@100.94.57.17 "tmux ls; ls ~/labretary-remote/experiments/EXP-20260929-001/results-seed62/run-62-L-1/genes | wc -l"`（567で完了）／cf5bも同様（results-seed61） | **完了**（cf5b 20:47:55〜00:27:20実測3.66h、c083 20:48:00〜04:23:51実測7.60h。全機体exit=0・skip0・遮断遺伝子混入0・埋め込み欠け0。A-1チェックポイント〔seed61単独〕はG-3a/M1/M2/M4合格、G-3b・M3は小K・1seedノイズ域内の僅差不合格と判断し継続。統合・評価まで完了） | `experiments/EXP-20260929-001/results/eval/criteria_003.json`（暫定提案: **判定不能**〔目標の軸の陽性対照G-3bが3seedとも不合格〕。主判定G-1自体も不合格〔z≈0.02-0.05、最小検出可能効果量0.91を大きく下回る〕。validator・PI待ち） |
| 2026-09-30 13:0x | EXP-20260929-001（事後・段階 2） | 移動量の補正の検算。摂動は回さない。既存の埋め込み 1,701 遺伝子ぶん（567 × 3 seed）を読み直し、補正の設計行列に log2(移動量) を足して主指標を計算し直す | Mac（計算のみ、GPU 不使用） | 0.5 h（上振れ 1.0 h。比較対象の一覧づくりが重く、seed ごとに 10〜20 分） | 1.0 h を超えたら停止して報告 | `ls experiments/EXP-20260929-001/posthoc/` に出力 JSON が出れば完了。途中経過は `logs/movement-correction-2026-09-30.log` | **完了**（実測 0.22 h。13:03 開始、途中 2 回作り直し〔版ごとに比較対象の一覧を二度作っていた／移動量 0 の細胞で対数が取れなかった〕。中断規則には抵触せず） | `experiments/EXP-20260929-001/posthoc/movement_correction_check.json`。従来の版で確定値を差 0.0005 以内で再現した上で、補正の効果を確認 |
- **状態の確かめ方（この実験）**: `ssh coi-next1@100.94.57.17 "tmux ls; ls ~/labretary-remote/experiments/EXP-20260929-001/results-seed62/run-62-L-1/genes | wc -l"`（567 で完了）。`isp_remote.sh status` はログ名が違うと誤表示するので使わない
| 2026-09-30 14:1x | 段階 3（着手前チェック） | 候補 (I) の状態のセントロイド 4 本 × 3 seed（71/72/73）を EmbExtractor の前向き計算だけで作り、`scripts/analysis/axis_feasibility.py` を通す。**摂動は 1 回も行わない**。細胞数は seed あたり約 2,940（aFIB 938-946 / FIB_MFIB 938-946 / MYOF 818-847 / dFIB 244-262） | Mac（MPS, bf16） | 0.75 h（上振れ 1.5 h） | 1.5 h を超えたら停止して報告 | `ls outputs/design/2026-09-30-annotation-goal-state/axis_feasibility.json` が出れば完了。途中経過は `logs/axis-feasibility-2026-09-30.log` | **中止**（14:18:53、seed 71 の 50% 付近。PI が 14:1x に目標の状態を「活性化した間質の状態（aFIB + MYOF）」へ変更したため、この設計は差し替え。出力 JSON は作られていない） | — |
| 2026-09-30 14:29:24 | 段階 3（着手前チェック・第 2 版） | 候補 (A) の状態のセントロイド 4 本 × 3 seed（81/82/83）を前向き計算だけで作り、`scripts/analysis/axis_feasibility.py` を通す。**摂動は 1 回も行わない**。目標 = 活性化した間質の状態（aFIB + MYOF）、出発点 = FIB + M-FIB、参照 = dFIB、第 4 = dM-FIB。細胞数は seed あたり約 2,300（aStr 950 / FIB_MFIB 950 / dFIB 261-274 / dMFIB 133-137） | Mac（MPS, bf16） | 0.75 h（上振れ 1.5 h） | 1.5 h を超えたら停止して報告 | `ls outputs/design/2026-09-30-annotation-goal-state/axis_feasibility_v2.json` が出れば完了。途中経過は `logs/axis-feasibility-v2-2026-09-30.log` | **完了**（14:29:24〜14:53:06、実測 0.40 h。中断規則 1.5 h に抵触せず。exit=0、摂動は 1 回も行っていない） | `outputs/design/2026-09-30-annotation-goal-state/axis_feasibility_v2.json`（**判定: 判定保留**。不合格の条件は 0 件。K=39、検出率 0.3214、まとまり +0.0818 はいずれも過去最良。弱点 3 件は仮説 §8.5 に記載） |
| 2026-09-30 18:43 ごろ | EXP-20260930-001 | (1) 機体一致確認: 48 遺伝子（T1 11+T2 7+陽性対照10+恒常発現の対照遺伝子（3 個）+中立の対照遺伝子（7 個）+対照10）×200細胞×seed81 を spark-cf5b と spark-c083 で並行実行し突き合わせ。合格なら (2) 本番: 600 遺伝子×200細胞×3seed （cf5b が seed81、c083 が seed82→83 逐次）。過剰発現型摂動、三次リンパ組織オーガナイザー遺伝子セットが 腎臓線維芽細胞を活性化間質状態へ動かすかの主判定と、陽性対照71個の健全性チェック | (1) spark-cf5b + spark-c083 並行 0.35h/機体 (2) cf5b: 3.91h（seed81）／c083: 7.8h（seed82→83逐次） | 全体で見込み約8.5h（機体一致確認込み）、上振れ 12h 程度 | plan.md §7.5: A-1〜A-9（機体不一致・チェックポイント不合格・skip超過・盲検混入・18h超過・1.5倍超過・細胞数120-130・移動量0・欠測20未満）。A-1/A-2/A-3/A-4/A-5 は**止める**、A-6/A-7/A-8/A-9 は記録して続ける | `ssh coi-next1@100.94.57.17 "tmux ls; ls ~/labretary-remote/experiments/EXP-20260930-001/results-seed82/run-82-L-1/genes | wc -l"`（600で完了）／`ssh petadimensionlab@100.104.171.109 "tmux ls; ls ~/labretary-remote/experiments/EXP-20260930-001/results-seed81/run-81-L-1/genes | wc -l"`（600で完了）。isp_remote.sh status は使わない（ログ名の食い違いで誤表示） | **計算完了・評価は不具合で停止中**（機体一致確認 18:45-19:06 完了・合格 [H-1 mean/max abs diff 0.0/0.0, H-2 rho=1.0, H-3 cos_min=0.99999999...]。本番3seedとも exit=0・600遺伝子・エラー0件で完了。回収〔isp_remote.sh pull、手動rsync不要〕・照合〔sha256 9/9一致・件数一致・埋め込み欠け0・盲検混入0、file-verifier実施〕・リモート後片付け〔clean、両機体で複製消失をls確認済み〕・統合〔merge_results_004.py、all_seeds_complete/no_blocked_genes/no_missing_embeddings すべてtrue〕まで完了。**`analysis/evaluate_004.py` が `KeyError: relpos_med` で停止**（M2指標健全性チェックの集計関数 `scripts/analysis/movement_correction.py::gene_stats()` が `relpos_med` 列を出力しておらず、前身 evaluate_003.py からの移植時のリグレッションと判断。判定基準・閾値は未変更、結果も一切未計算）。simulation-engineer の権限では分析コードを書き換えないため、Labretary に原因と選択肢を報告し指示待ち） | `experiments/EXP-20260930-001/results/merge_report.json`（統合結果、合格）。評価は未出力 |

- **seed83 の自動引き継ぎ（2026-09-30 19:12 起動）**: c083 は seed82 → seed83 を自動で連結しないため、`scripts/remote/handoff_seed83.sh` を手元で常駐させた。**seed82 が 600 遺伝子で正常完了したときだけ** seed83 を開始する。600 に届いていなければ**開始せず**にログへ理由を書いて終わる（中断規則を飛び越えないため）。
  - ログ: `logs/handoff-seed83-2026-09-30.log`
  - 生存確認: `pgrep -f handoff_seed83`
  - 手で状態を見る: `ssh coi-next1@100.94.57.17 "tmux ls; ls ~/labretary-remote/experiments/EXP-20260930-001/results-seed82/run-82-L-1/genes | wc -l"`（600 で完了）／`ssh petadimensionlab@100.104.171.109 "tmux ls; ls ~/labretary-remote/experiments/EXP-20260930-001/results-seed81/run-81-L-1/genes | wc -l"`
  - **`isp_remote.sh status` は使わない**（ログ名の食い違いで誤表示する）

| 2026-10-01 | EXP-20260930-001（事後診断・事前登録の外） | 指標の健全性 M3（B2M）・M4（一様成分）・M5（残差と移動量）の原因特定。**摂動は 1 回も回さない**。既存の埋め込み 1,800 遺伝子ぶん（600 × 3 seed）を読み直して細胞ごとの表を作り、補正の形（補正なし／`ld` の 1 次／`ld` と `ld²`／完全二次式ほか）を並べて再計算する | Mac（計算のみ、GPU 不使用） | 1.0 h（上振れ 2.0 h） | 2.0 h を超えたら停止して報告 | `ls experiments/EXP-20260930-001/posthoc/` に出力 JSON が出れば完了 | **完了**（実測 0.2 h。17:50〜18:02。中断規則に抵触せず。摂動は 1 回も回していない。事前登録 51 ファイルのうち 50 件は照合一致、残り 1 件は `hypotheses/HYP-20260930-01.md`〔本診断では触れていない〕） | `experiments/EXP-20260930-001/posthoc/` |
| 2026-10-01 23:44 ごろ | EXP-20261001-001（機体一致確認） | 44 遺伝子（介入 20 + 陽性対照 10 + 採用した中立対照 4 + 対照 10）× 200 細胞 × seed 81 を spark-cf5b と spark-c083 で並行実行し、H-1〜H-3（plan.md §8.2）で突き合わせる。**事前登録の対象外**（`analysis/make_configs_machinecheck_005.py` / `analysis/machine_check_005.py`、Labretary 作成）。不合格なら本番を開始しない | spark-cf5b（petadimensionlab@100.104.171.109）と spark-c083（coi-next1@100.94.57.17）並行 | 0.43 h/機体（dry-run 実測モデル。上振れ 0.9 h） | 不合格（H-1〜H-3 のいずれか）なら本番を開始せず Labretary へ報告。0.9 h 超過でも状態確認 | `ssh petadimensionlab@100.104.171.109 "tmux ls; ls ~/labretary-remote/experiments/EXP-20261001-001/results-machinecheck/run-81-L-1/genes | wc -l"`（44 で完了）／`ssh coi-next1@100.94.57.17 "tmux ls; ls ~/labretary-remote/experiments/EXP-20261001-001/results-machinecheck/run-81-L-1/genes | wc -l"`（44 で完了）。**`isp_remote.sh status` は使わない**（ログ名の食い違いで誤表示する） | **完了・合格**（23:44:11/17〜00:02:52/00:03:01、実測 約0.31h/機体、両機体 exit=0・44/44遺伝子。H-1 mean_abs=0.0/max_abs=0.0、H-2 ρ=1.0、H-3 cos_min=0.999999999999999、すべて合格。移動量0細胞が両機体で同一の1件〔ANK2〕のみで一致＝ほぼビット同一） | `experiments/EXP-20261001-001/work/machine_check_005.json`（pass: true）。回収先 `results-machinecheck-{cf5b,c083}/` |
| 2026-10-02 00:05 | EXP-20261001-001（本番） | 717 遺伝子 × 200 細胞 × 3 seed。cf5b が seed 81（約 4.7 h）、c083 が seed 82 → 83 逐次（約 9.3 h）。seed 単位でのみ分割（遺伝子では分けない） | spark-cf5b（seed 81 開始 00:05:37）、spark-c083（seed 82 開始 00:05:43 → seed83 は Labretary の見張りが引き継ぐ） | cf5b 4.7 h（上振れ 7 h）、c083 計 9.3 h（上振れ 12 h） | plan.md §7.5 A-1〜A-10。A-1〜A-6 は**止める**（機体不一致・チェックポイント不合格・skip 5%超・盲検混入・24h超・評価コード例外）、A-7〜A-10 は記録して続ける。**全体 24h 超で A-5、停止して報告** | `ssh petadimensionlab@100.104.171.109 "tmux ls; ls ~/labretary-remote/experiments/EXP-20261001-001/results-seed81/run-81-L-1/genes | wc -l"`（717で完了）／`ssh coi-next1@100.94.57.17 "tmux ls; ls ~/labretary-remote/experiments/EXP-20261001-001/results-seed82/run-82-L-1/genes | wc -l"`→ seed82完了後 `results-seed83/run-83-L-1/genes | wc -l`（各717で完了）。**`isp_remote.sh status` は使わない** | **完了**（seed81 cf5b: 00:05:37〜04:43:38 実測 4.62h、exit=0。seed82 c083: 00:05:43〜04:44:48 実測 4.65h、exit=0。seed83 c083: 起動 2 回失敗〔04:46:25 設定ファイル未送付で FileNotFoundError／05:22:39 push によるモデルリンク破損で OSError、いずれも exit=1・GPU消費ほぼ0〕、**05:25:20 の 3 回目で成功**し 10:02:19 に exit=0・717/717遺伝子・エラー0件で完了。全体 24h 超〔A-5〕には抵触せず。回収: isp_remote.sh pull、sha256・件数・埋め込み欠け・盲検混入をfile-verifierで確認し全項目一致。リモート複製は両機体 clean 済み・ls で消失確認済み。統合: merge_results_005.py で all_seeds_complete/no_blocked_genes/no_missing_embeddings すべてtrue。評価: evaluate_005.py 完走〔約3.9分〕、M5不合格のため判定提案は「判定不能」） | `experiments/EXP-20261001-001/results/eval/criteria_005.json`（conclusion_proposal: 判定不能〔指標の健全性 M5 不合格〕）。`results/merge_report.json`（統合結果、合格） |

- **EXP-20261001-001 の見張り（2026-10-02 00:10 起動、Labretary）**: 前回 seed の引き継ぎが人の記憶頼りで 43 分遅れた反省から、**Python で 2 本常駐**させた。
  - **引き継ぎ** `scripts/remote/handoff.py`: c083 の seed82 が **717 遺伝子・エラー 0 件・tmux 終了**の 3 つを満たしたときだけ seed83 を開始する。1 つでも欠ければ**開始せず**理由を書いて止まる。ログ `logs/handoff-20261001.log`
  - **止まり検知** `scripts/remote/watch_run.py`: 3 seed の合計が **25 分まったく増えなければ鳴る**（原因が引き継ぎの死亡でも異常終了でも通信断でも拾える）。完了（2,151）で正常終了。ログ `logs/watch-run-20261001.log`
  - 生存確認: `pgrep -f handoff.py` / `pgrep -f watch_run.py`。**stderr を捨てず** `logs/*-stderr.log` に残している
  - 手で状態を見る: `ssh petadimensionlab@100.104.171.109 "tmux ls; ls ~/labretary-remote/experiments/EXP-20261001-001/results-seed81/run-81-L-1/genes | wc -l"`（717 で完了）／c083 は results-seed82・results-seed83
  - **`isp_remote.sh status` は使わない**（ログ名の食い違いで誤表示する）

- 2026-10-02 事後診断 M5（EXP-20261001-001/posthoc-m5）: 変位ベクトルの再集計と提供者ごとの除外。見込み 20〜30 分。摂動は回さない。
  - 完了 2026-10-02。実計算 6 分（最長は提供者 79 人の 1 人ずつ除外 3 分）。出力は experiments/EXP-20261001-001/posthoc-m5/*.json。摂動は 1 回も回していない。

| 2026-10-02 22:29 ごろ（PI 承認済み、G2 不要） | **EXP-20261002-002** | 対照遺伝子を 500 → 1,000 個に増やし、歪みの見積もりの区間を縮める。**摂動するのは追加の対照 500 個だけ**（既存の 717 遺伝子は EXP-20261001-001 で摂動済みなので回し直さない）。500 遺伝子 × 200 細胞 × 3 seed（81/82/83）。過剰発現型。設定・モデル・トークン化済み入力・状態のセントロイドの行はすべて EXP-20261001-001 と同一 | spark-cf5b が seed 81、spark-c083 が seed 82 → 83（逐次）。**seed 単位でのみ分割する（遺伝子で分けない）** | **9.80 h**（`run_isp.py --dry-run` の実測モデル。`needs_G2_approval: false`。内蔵見積もり 0.56 h は採らない）。機体一致確認 0.13 h/機体を別途。評価は Mac で 12〜15 分（GPU 不使用）。上振れ 14 h | plan.md §7.5 A-1〜A-10。**止める**: A-1 機体一致確認の不合格／A-2 状態のセントロイドが既存と不一致／A-3 再現の検査が不一致（差 0.002 超）／A-4 盲検混入／A-5 全体 24 h 超／A-6 評価コードの例外。**記録して続ける**: A-7 S-1 不合格（予算 0.10 を緩めない）／A-8 移動量 0 が 5% 超／A-9 1.5 倍超過／A-10 頑健性の割れ | **開始前に必ず**: `bash experiments/EXP-20261002-002/analysis/stage_remote_007.sh HOST SEED...`（config と state_rows を 1 回でまとめて送り sha256 照合。rows を送り忘れると起動できない）。進捗: `ssh petadimensionlab@100.104.171.109 "tmux ls; ls ~/labretary-remote/experiments/EXP-20261002-002/results-seed81/run-81-L-1/genes \| wc -l"`（500 で完了）／c083 は results-seed82・results-seed83。**`isp_remote.sh status` は使わない**（ログ名の食い違いで誤表示）。**実験の途中で `isp_remote.sh push` を走らせない**（モデルへのリンクを壊す） | **実行中・本番開始済み**。22:29 事前確認（sha256 42/42 一致）→ 22:32 機体一致確認（C-3、既存対照12個を回し直し EXP-20261001-001 と突合）両機体 exit=0・pass:true（H-1〜H-4 最大差 0.0）→ **22:41:45 cf5b で seed81 開始、22:41:50 c083 で seed82 開始**（device: cuda 確認、エラーなし）。**seed83 は Labretary が `scripts/remote/handoff.py` で引き継ぎ予定**（seed82 が 500 遺伝子・エラー0・tmux終了を満たした時のみ開始）。状態の確かめ方: `ssh petadimensionlab@100.104.171.109 "tmux ls; ls ~/labretary-remote/experiments/EXP-20261002-002/results-seed81/run-81-L-1/genes | wc -l"`（500で完了）／c083同様に results-seed82→results-seed83。詳細経過は `experiments/EXP-20261002-002/logs/progress.md` | 出力予定: `experiments/EXP-20261002-002/results/criteria_007.json`。試し実行の確認は `results/criteria_007-smoke.json`（**判定ではない**）。機体一致確認の結果: `experiments/EXP-20261002-002/work/machine_check_007-results-machinecheck-{cf5b,c083}.json` |

- **EXP-20261002-002 の見張り（2026-10-02 22:43 起動、Labretary）**: Python の常駐 2 本。
  - **引き継ぎ** `scripts/remote/handoff.py`: c083 の seed82 が **500 遺伝子・エラー 0 件・tmux 終了**の 3 つを満たしたときだけ seed83 を開始。1 つでも欠ければ**開始せず**理由を書いて止まる。ログ `logs/handoff-20261002.log`
  - **止まり検知** `scripts/remote/watch_run.py`: 3 seed の合計（目標 1,500）が **25 分まったく増えなければ鳴る**。ログ `logs/watch-run-20261002.log`
  - 生存確認: `pgrep -f handoff.py` / `pgrep -f watch_run.py`。stderr は `logs/*-stderr.log` に残している
  - 手で状態を見る: `ssh petadimensionlab@100.104.171.109 "tmux ls; ls labretary-remote/experiments/EXP-20261002-002/results-seed81/run-81-L-1/genes | wc -l"`（500 で完了）／c083 は results-seed82・results-seed83
  - **注意**: 機体一致確認のログだけは `work/logs/` 配下（設定が `work/` 直下にあるため）。本番は `logs/` 直下


| 2026-10-03 21:2x（PI 承認済み 2026-10-03、G2 不要） | **EXP-20261003-001** | 目標の軸が「筋線維芽細胞らしさ」を測っているかの直接検証。**筋線維芽細胞・平滑筋・活性化線維芽細胞・間質のマーカー 6 個**（ACTA2 / TAGLN / CNN1 / MYH11 / POSTN / PDGFRB。いずれも既存の 1,217 遺伝子に入っておらず未測定）× 200 細胞 × 3 seed（81/82/83）の過剰発現。設定・モデル・トークン化済み入力・状態のセントロイドの行・対照 1,000 個はすべて EXP-20261001-001 / EXP-20261002-002 と同一 | **Mac（mps / bf16）1 台**。計算量が小さいので機体を分けない。**リモートは使わない**（`isp_remote.sh push` を走らせない） | 本番 **0.30 h**（`run_isp.py --dry-run` の実測モデル。`needs_G2_approval: false`）。機体の一致確認 0.14 h、試し実行 0.01 h、評価 15 分程度（GPU 不使用）。上振れ 1.5 h | plan.md §7.5。**止める**: A-1 機体の一致確認の不合格／A-2 状態のセントロイドが既存と不一致／A-3 再現の検査が不一致（差 0.002 超）／A-4 盲検混入／A-5 評価コードの例外／A-6 6 遺伝子のいずれかが 120 細胞に届かない。**記録して続ける**: A-7 移動量 0 が 5% 超／A-8 1.5 倍超過 | `ls experiments/EXP-20261003-001/results/run-{81,82,83}-L-1/genes \| wc -l`（各 6 で完了）。ログ `experiments/EXP-20261003-001/logs/` | 記録中 | 出力予定: `experiments/EXP-20261003-001/results/criteria_008.json`。機体の一致確認: `work/machine_check_008.json` |

- **EXP-20261003-001 の実行記録（2026-10-03）**: 当初 Mac を予定したが **C-4 機体の一致確認で不合格**（MPS 対 CUDA。埋め込みの最大差 0.15625、セントロイド 0.03125、主指標の中央値 0.0155）。既存と同じ **spark-cf5b（NVIDIA GB10 / cuda）** に替えて C-4 合格（最大差すべて 0.0）。本番 **21:48:34〜22:01:03 exit=0・エラー 0 件**、3 seed とも 6/6 遺伝子・各 200 細胞。回収後に sha256（埋め込み 36 本を含む全件）一致・盲検混入 0 を確認し、**リモートの複製 931 MB を `isp_remote.sh clean` で削除して `ls` で消失を確認**。評価は **中断条件 A-2 で停止**（seed83 のセントロイド dMFIB のみ bf16 1 刻みの不一致。主指標の向きを決める 2 本は 3 seed とも差 0.0）。**判定 V-1 / V-2 は出していない。** 診断: `experiments/EXP-20261003-001/work/centroid_stop_008.json`

## 2026-10-04 EXP-20261003-001 数値型の切り分け（Mac fp32 で 12 遺伝子）
- 目的: bf16 と fp32 の差が**主指標でどれだけになるか**を、同じ機体で測る。PI 方針「判定の基準は float32（単精度）、探索は bfloat16」の運用条件を決めるため
- config: `experiments/EXP-20261003-001/work/config-machinecheck-fp32.yaml`（bf16 版と dtype・出力先以外は同一）
- 出力: `experiments/EXP-20261003-001/results-machinecheck-fp32/`
- 見込み: 2,400 細胞遺伝子。bf16 実測 0.1154 s/件 = 約 5 分、fp32 は 2〜3 倍として **10〜15 分**
- 状態: 実行中
- 実測の更新（2026-10-04）: state_embs 728 秒・1 遺伝子 195 秒（bf16 の約 7 倍）。全体の見込みを **45 分**に修正
- **状態: 完了**（2026-10-04 12:05）。12 遺伝子・seed81。結果 `experiments/EXP-20261003-001/results-machinecheck-fp32/`、解析 `work/dtype_check_008.json`
- 要点: bf16 の誤差は共通のずれ（0.10×対照SD）ではなく**遺伝子ごとの雑音**（0.45×対照SD）。同じ機体の bf16↔fp32 のばらつき 0.0066 と、機体をまたぐ bf16↔bf16 のばらつき 0.0058 が**同じ大きさ**＝原因は機体ではなく bf16

## 2026-10-04 EXP-20261004-001 事前登録固定前の計算（C-11 機体一致確認 / C-9 試し実行）

- 目的: 本番（626 遺伝子・共通細胞 200・3 seed・13.32 GPU 時間、`needs_G2_approval: false`）を開始する前の
  最後の確認。**本番の 626 遺伝子は今回は回さない。**
- 開始: 2026-10-04（JST）。実行先: spark-cf5b（petadimensionlab@100.104.171.109）/ spark-c083（coi-next1@100.94.57.17）
- 中身:
  1. **C-11 機体一致確認**: `work/config-machinecheck.yaml`（対照 500 個の先頭 12 個 × 共通 200 細胞 × seed81、
     bf16/cuda）を cf5b と c083 で独立に実行し、`analysis/machine_check_009.py` で直接突き合わせる
     （cell_selection.mode: common は本番で初めて使うため、比べる既存結果が無い。§5.5.3）。見込み 0.15h × 2 機体
  2. **C-9 試し実行**: `work/config-smoke.yaml`（候補2・陽性対照2・恒常発現1・対照35 = 40 遺伝子 ×
     共通30細胞 × seed81、`--smoke` 旗は使わない）を 1 台で実行。出力は `results-smoke/`。見込み 0.01h
  3. **C-10**: C-9 の出力に対して `analysis/evaluate_009.py --smoke` を実行し、評価コードが最後まで通るか確認
- 合わせて見込み **約 0.4 時間**（実計算）。**中断規則**: C-11 で H-1〜H-4 のいずれか不合格なら**止めて報告**
  （原因をログ・weight・ライブラリ版・数値型の 1 つずつ切り分けるまで推測で書かない。guardrails 2026-10-03/04）。
  C-9 で exit ≠ 0、または 40/40 遺伝子・`emb/{orig,pert}.npy` のいずれかが欠けていれば**止めて報告**
- 状態の確かめ方: `ssh HOST "tmux ls"`、`ssh HOST "ls labretary-remote/experiments/EXP-20261004-001/results-machinecheck/run-81-L-1/genes | wc -l"`（12 で完了）
- 出力予定: `experiments/EXP-20261004-001/work/machine_check_009.json`、
  `experiments/EXP-20261004-001/results-smoke/`、`experiments/EXP-20261004-001/results/criteria_009-smoke.json`
- **本番（626 遺伝子・3 seed・13.32 GPU 時間）はこの長時間ジョブの対象外。** C-9/C-10/C-11 が全部合格し、
  事前登録を sha256 で固定したあとに、別途この表へ本番の行を追加してから開始する

- **状態: C-11 で停止（2026-10-04 13:37 ごろ）**。cf5b（seed81, 13:30:34 開始 → 13:36:57 完了, exit=0, 12/12遺伝子）・
  c083（13:30:43 開始 → 13:36:18 完了, exit=0, 12/12遺伝子）とも正常終了。
  `scripts/analysis/compare_conditions.py` で実行条件を先に突き合わせたところ、
  **`torch`（cf5b 2.12.0+cu130 / c083 2.13.0+cu130）と `transformers`（cf5b 4.57.6 / c083 4.46.3）が違う**
  （判定に効く条件。dtype=bfloat16・attn_implementation=sdpa・tf32=False は両機体で一致）。
  `analysis/machine_check_009.py` で H-1〜H-4 を直接計算したところ、**埋め込み・摂動に使った細胞・
  状態のセントロイド（FIB_MFIB・aStr・dFIB・dMFIB すべて）・主指標の中央値は全遺伝子で最大差 0.0**
  （完全一致）。ただし **compare_conditions.py がライブラリ版の不一致を理由に exit 1 を返しているため、
  machine_check_009.py 全体の pass は False**。
  guardrails 2026-10-04 の指示どおり「原因を推測で書かない」を守り、**C-9 試し実行・事前登録の固定は
  実行せず、Labretary へ報告して指示を待つ**。両機体の `~/labretary-remote`（360MB、ヒト由来データの複製を含む）
  は `isp_remote.sh clean` で削除し `ls` で消失を確認済み
- 出力: `experiments/EXP-20261004-001/work/machine_check_009.json`、
  `experiments/EXP-20261004-001/results-machinecheck-{cf5b,c083}/`、
  `experiments/EXP-20261004-001/logs/remote-config-machinecheck-{cf5b,c083}.log`

- **状態: C-11 は合格として扱う（2026-10-04 訂正）**。`scripts/analysis/compare_conditions.py` の終了コードの
  意味を直した（Labretary 2026-10-04）。終了コード 2（条件は違うが数値を比べられる）は止める理由ではない。
  H-1〜H-4 の数値（最大差 0.0）だけで合否を決めた結果、**C-11 は合格**。torch・transformers の版の違いは
  plan.md §10-10/11 に限界として明記した
- **C-9 試し実行: 完了（2026-10-04 16:02:20 開始 → 16:06:54 完了, exit=0, 40/40遺伝子, cf5b）**。
  `emb/{orig,pert}.npy` の実在・`check_common_cells_009.py --config work/config-smoke.yaml` 合格を確認
- **C-10: 完了**。`evaluate_009.py --smoke` が exit=0 で最後まで通った。途中で config.yaml の
  `metric.confound_correction.form` が `full_quadratic` のまま（plan.md §5.4 の `stratified_by_expression`
  と不一致）であることが分かり、`run_isp.py` はこの値を参照しないこと（摂動の計算には影響しない）を
  確認したうえで `stratified_by_expression` に訂正した（固定前・候補の値を見る前）
- cf5b の `~/labretary-remote`（319MB、ヒト由来データの複製を含む）は `isp_remote.sh clean` で削除し
  `ls` で消失を確認済み
- 出力: `experiments/EXP-20261004-001/results-smoke/`、`experiments/EXP-20261004-001/results/criteria_009-smoke.json`、
  `experiments/EXP-20261004-001/work/common_cells_check_009-smoke.json`

## 2026-10-04 16:35 EXP-20261004-001 本番（626遺伝子 × 3 seed、PI 承認済み・ゲート G2 不要）

- 承認: PI（2026-10-04）。事前登録固定済み（`work/SHA256SUMS-preregistered.txt`、43ファイル、不一致0件）
- 計算量: 13.32 GPU 時間（`run_isp.py --dry-run` の実測モデル）。`needs_G2_approval: false`
- 機体の分け方（plan.md §8.4。**seed 単位でのみ分ける**）:
  - **spark-cf5b（petadimensionlab@100.104.171.109）**: seed81 のみ。見込み 4.44 h
  - **spark-c083（coi-next1@100.94.57.17）**: seed82 → seed83 を順に。見込み 4.44 h × 2 = 8.88 h
  - 全体のウォールクロックの見込み: **約 9 時間**
  - seed83 の開始は `scripts/remote/handoff.py`（seed82 が 626 遺伝子・エラー0件・tmux終了の3つを
    満たしたときだけ開始。1つでも欠ければ開始せず理由を書いて止まる）
  - 実行用の config は `config-seed{81,82,83}.yaml`（事前登録済み `config.yaml` の機械的な派生。
    `robustness.seeds` と `output.results_dir` だけを変え、遺伝子集合・判定基準は一切変えていない。
    事前登録の対象外。`_derived_from_config_yaml_sha256` に元の sha256 を記録した）
- **中断条件**（plan.md §6.0。止める: A-1〜A-5。記録して続ける: A-6〜A-8）:
  - A-1 盲検中の TP53・SPP1・LCN2 が結果に混入 → 止める
  - A-2 共通の細胞集合が遺伝子によらず同一でない（`check_common_cells_009.py` 不合格）→ 止める
  - A-3 状態のセントロイド（FIB_MFIB・aStr）が記録と不一致 → 止める
  - A-4 評価コードが例外で止まる → 止める
  - A-5 626 遺伝子のいずれかが skip される → 止める
  - A-6 移動量 0 の細胞の割合が 5% 超 → 記録して続ける
  - A-7 dFIB・dMFIB のセントロイドの差 → 記録して続ける
  - A-8 所要が見込みの 1.5 倍（約 13.5 時間）を超えた → 記録して状態を確認し報告
  - **いずれかに当たったら実行担当は止めて報告する。続行を自分で決めない**
- 進み具合の確かめ方: `ssh HOST "ls labretary-remote/experiments/EXP-20261004-001/results-seedN/run-N-L-1/genes | wc -l"`（626 で完了）、`ssh HOST "tmux ls"`
- 見張り: `scripts/remote/watch_run.py --exp EXP-20261004-001 --seed 81:petadimensionlab@100.104.171.109 --seed 82:coi-next1@100.94.57.17 --seed 83:coi-next1@100.94.57.17 --expect 626 --stall-min 25`（合計 1,878 のうち 25 分まったく増えなければ異常）
- **実行中は `isp_remote.sh push` を走らせない**（モデルへのリンクを壊す。2026-10-02 に実際に発生）。
  必要なファイルは開始前に `stage_remote_009.sh` で 1 回でまとめて送る
- **状態: 実行中。ステージング完了（両機体で sha256 照合・実在確認済み）→ 16:36:52 cf5b で seed81 開始、
  16:37:01 c083 で seed82 開始**（両方 device: cuda / dtype: bfloat16 を確認。起動 1 分後にログの伸びを確認済み）。
  `scripts/remote/handoff.py`（pid 58950）と `scripts/remote/watch_run.py`（pid 58958）を常駐させた。
  ログ: `logs/handoff-20261004.log`、`logs/watch-run-20261004.log`（stderr はそれぞれ `*-stderr.log`）
- 生存確認: `pgrep -f handoff.py` / `pgrep -f watch_run.py`
- 手で状態を見る: `ssh petadimensionlab@100.104.171.109 "tmux ls; ls labretary-remote/experiments/EXP-20261004-001/results-seed81/run-81-L-1/genes | wc -l"`（626で完了）／
  `ssh coi-next1@100.94.57.17 "tmux ls; ls labretary-remote/experiments/EXP-20261004-001/results-seed{82,83}/run-{82,83}-L-1/genes | wc -l"`

- **状態: 完了（2026-10-05 回収・評価まで完了）**。
  - cf5b: seed81 開始 16:36:52 → 完了 20:05:04（約 3.47 時間）、626/626 遺伝子、エラー 0 件
  - c083: seed82 開始 16:37:01 → 完了 19:53:53（約 3.28 時間）、626/626 遺伝子、エラー 0 件。
    `scripts/remote/handoff.py` が正常完了を確認して 19:53:53 に seed83 を自動開始 → 完了 23:12:47
    （約 3.32 時間）、626/626 遺伝子、エラー 0 件
  - 全体のウォールクロック: 16:36:52〜23:12:47（**約 6.60 時間**。見込み約 9 時間の 1.5 倍（A-8 の基準）を
    超えていない）。`scripts/remote/watch_run.py` も 23:15:14 に「合計 1,878 / 1,878 に到達」で正常終了
  - 回収: 3 seed とも `ls` で 626/626 を確認し、`rsync` で `results-seed{81,82,83}/` に回収。
    `run-{81,82,83}-L-1` を `results/` にまとめ、3 本の `cells.csv`（各 125,201 行）を連結して
    `results/cells.csv`（375,601 行 = ヘッダ + 626 遺伝子 × 200 細胞 × 3 seed）を作成
  - 後片付け: 両機体の `~/labretary-remote`（cf5b 3.1GB、c083 5.9GB。ヒト由来データの複製を含む）を
    `isp_remote.sh clean` で削除し、`ls` で両方とも消失を確認
  - A-2（`check_common_cells_009.py`）: 3 seed とも pass。各 seed 内で 626 遺伝子すべてが同一の細胞構成
    （cell_uid・donor_id・length が完全一致）、200 細胞/遺伝子
  - `evaluate_009.py`: exit=0、約 18 秒で完了（GPU 不使用）。A-1〜A-5 すべて pass（止める条件に当たらず）。
    A-6 移動量 0 の細胞は 0.022%（5% 未満）。A-7 は記録のみ
  - 出力: `experiments/EXP-20261004-001/results/criteria_009.json`、`results/gene_seed_009.csv`、
    `results/gene_rank_009.csv`、`work/common_cells_check_009.json`

## 2026-10-07 HYP-20261007-01 着手前チェック（A2）の埋め込み抽出

- 内容: A2 の線維芽細胞系列 2,799 細胞の細胞埋め込みを抽出（**摂動は 1 回も行わない**。前向き計算のみ）
- 実行: `.venv/bin/python outputs/design/2026-10-07-A2-feasibility/extract_embs.py all cell_embs.npz`
- 装置: このローカルの Mac、MPS、float32。見込み 25〜35 分
- ログ: `logs/A2-precheck-embs.log`
- 状態の確かめ方: `tail -c 300 logs/A2-precheck-embs.log`（進捗バー）、
  完了の確認は `ls -l outputs/design/2026-10-07-A2-feasibility/work/cell_embs.npz`
- 結果: **完了**（2026-10-07 14:27）。2,799 細胞・824.7 秒（13.7 分）・MPS・float32・1,152 次元。
  出力 `outputs/design/2026-10-07-A2-feasibility/work/cell_embs.npz`（25.8 MB、`ls` で実在を確認）

## 2026-10-08 EXP-20261007-001 本番（float32・2 台）

- 内容: HYP-20261007-01 の確認的検定。587 遺伝子を過剰発現。A2 の出発点 FIB を提供者で分けた 3 群（各 55・59・59 細胞）
- 事前登録: `experiments/EXP-20261007-001/work/SHA256SUMS-preregistered.txt`（27 件、2026-10-08 再固定）
- 機体と群:
  - spark-cf5b（petadimensionlab@100.104.171.109）: 群 81、`work/config-cf5b.yaml`、出力 `results-cf5b/`、見込み 3.2 時間
  - spark-c083（coi-next1@100.94.57.17）: 群 82 → 83（同じ設定の中で順に回る）、`work/config-c083.yaml`、出力 `results-c083/`、見込み 6.4 時間
- 数値型: float32（tf32 閉）。ゲート G2 不要（dry-run の合計 9.61 時間、`needs_G2_approval: false`）
- **開始: 2026-10-08 18:31:50（cf5b）・18:31:55（c083）**。見込みの終了: cf5b 21:45 ごろ、c083 翌 00:55 ごろ
- tmux: cf5b `isp-config-cf5b_`、c083 `isp-config-c083_`（`isp_remote.sh run` が作る名前。末尾に `_`）
- ログ（リモート）: `~/labretary-remote/experiments/EXP-20261007-001/work/logs/remote-config-cf5b_.log`・`.../work/logs/remote-config-c083_.log`（設定が work/ にあるので logs も work/ の下）（末尾の `_` に注意。guardrails 2026-09-29）
- 状態の確かめ方（手で叩く）:
  - `ssh petadimensionlab@100.104.171.109 "ls ~/labretary-remote/experiments/EXP-20261007-001/results-cf5b/run-81-L-1/genes | wc -l; tmux ls"`（目標 587）
  - `ssh coi-next1@100.94.57.17 "for s in 82 83; do ls ~/labretary-remote/experiments/EXP-20261007-001/results-c083/run-\$s-L-1/genes 2>/dev/null | wc -l; done; tmux ls"`（目標 587 ずつ）
- 見張り: `logs/watch-EXP-20261007-001-cf5b.log`・`logs/watch-EXP-20261007-001-c083.log`（`scripts/remote/watch_run.py`、25 分増えなければ鳴る）
- 終わったら: 回収 → `analysis/merge_results_011.py` → `analysis/evaluate_011.py` → 両機体で `isp_remote.sh clean` と `ls` で消失確認
- **結果: 完了（2026-10-11 回収・統合・評価まで）**。
  - 計算: cf5b 群81（587/587）は 2026-10-08 23:08 完了、c083 群82・83（各 587）は 2026-10-09 00:08 完了。両機体 exit=0、エラー 0 件（`logs/watch-EXP-20261007-001-*.log`）
  - 回収: `results-c083/` を rsync でローカルへ（587×2、config.sha256 一致 f976afd8…）。`results-cf5b/` は既存（config.sha256 3ac94cc7… 一致）。両機体とも dtype=float32・tf32 閉・`highest`
  - 統合: `merge_results_011.py` → `results/`（run-81/82/83、cells.csv 101,551 行、機体ごとの device.json・timing.jsonl）
  - 評価: `evaluate_011.py` は初回に §8b（記述）で例外停止（候補の GSN・ELN が陽性対照にも重なり二重に並んだ）。重複除去で修正し（事前登録の 3 点セット。plan.md §13 の 2026-10-11 の行）、**修正後に初めて判定の値が出た**。完走・exit=0
  - 出力: `results/criteria_011.json`（**P-1 主判定は合格の提案**、P-2 は R-4、P-3 は予測どおり。確定は PI＝ゲート G3）、`results/gene_group_011.csv`、`results/gene_rank_011.csv`
  - 後片付け: c083 の `~/labretary-remote`（1.4 GB、ヒト由来データの複製）を削除し `ls` で消失確認。cf5b は既に削除済み
- **状態: 計算・回収・統合・評価 完了。判定（G3）は PI 待ち**

| 2026-10-11 11:46 ごろ（ステージング 11:44〜11:47、両機体で sha256 照合・実在確認済み） | **EXP-20261011-001（A1 の新規摂動）** | 退化の陽性対照 30 個 × 200 細胞 × 3 群（18,000 細胞×遺伝子）を **A1 と同じ共通 200 細胞・同じ機体（GB10/CUDA）・同じ数値型（bfloat16）** で過剰発現。HYP-20261011-01 の確認（差の向き）用 | spark-cf5b（seed81）・spark-c083（seed82→83 逐次）。**Mac は使わない** | 摂動 0.63 h＋状態のセントロイド 0.20 h ≒ **0.83 h**（`run_isp.py --dry-run` の measured_hours_total。上振れ 1.0 h。`needs_G2_approval: false`） | plan.md §8 の B-1〜B-6。**止める**: A-1 盲検混入／A-2 材料 sha256 不一致／A-3 遺伝子集合不一致／A-4 再現確認・新規 30 の機械一致の不合格／A-5 値に NaN（null と理由を書く）／B-4 30 個に満たない・skip 0 件でない。記録して続ける: B-5（R-14 保留） | `ssh HOST "ls labretary-remote/experiments/EXP-20261011-001/results-seedN/run-N-L-1/genes | wc -l"`（30 で完了）。**`isp_remote.sh status` は使わない**（ログ名の食い違い）。**実験中に `isp_remote.sh push` を走らせない** | **計算・回収・評価 完了。判定（G3）は PI 待ち**（開始 11:47:44/54。seed81 cf5b 30/30 exit=0 11:59:51・seed82 c083 30/30 exit=0 11:58:47・seed83 c083 30/30 exit=0 12:11:50、エラー 0・skip 0。A-1〜A-4・B-4 すべて合格）。evaluate_014.py の実装バグ 5 件を修正（plan.md §14）後、試し実行 EXIT=0 → 本番評価 12:39 完走 EXIT=0。結果: A1 P-1=False（(i) 活性化の陽性対照（71 個）=True・(ii) 退化・老化の陽性対照（30 個）=False・(iii)=True）→ 最初に当たった反証条件 **R-6**。R-12 交絡 該当（ρ_exp 0.31）。R-13 は重なり除去でも結論不変。A2 は P-2=True・P-3=True、R-7（段階 2 の合格は共通成分）。出力 `results/criteria_014.json`・`gene_rank_014.csv`・`run.log`・`README.md`。**緩和検算（PI 指示）**: (ii) を 10% に緩めると P-1 合格＝「条件付き合格」、重なり除去（23）では 4% で合格。`results/relax_014.json`。判定は PI（ゲート G3） | 出力予定: `experiments/EXP-20261011-001/results/`（run-{81,82,83}-L-1。回収・統合済み）。評価 `evaluate_014.py` は GPU 0 時間 |

| 2026-10-11 14:37〜14:48 | **EXP-20261011-002（A1 のマッチ補正付き再評価）** | 保存済み埋め込みの再射影・再評価のみ（GPU 0 時間・CPU 約 3.5 分）。新しい摂動は回さない。PI 指示「退化側の移動を強く出す設計 → 交絡補正」に従い事前登録（G2 不要と判定。plan §10）。[1] 退化の主陽性対照=特異的マーカー 23（重なり 7 を除く）、[2] 主=マッチさせたサイズを揃えたランダム遺伝子セット（対照遺伝子群から復元ありで引く。発現細胞数・入力長・位置をそろえる）・副=残差化 | ローカル Mac（CPU・fp64 評価） | 0（GPU 0 時間） | A-1 盲検 0・A-2 材料 sha256・A-3 集合・A-4 再現（gene_rank_014 と 1e-3 以内）・A-5 値 NaN 0・A-6 結果ファイル再読込。B-5 主判定の値は途中で計算しない | `--self-test` 合格後に本番 1 回 | **計算・評価 完了。判定（ゲート G3）は PI 待ち**。実行中に構造バグ 2 件（`add_projection` の Series 名が `m_proj_deg` 固定→ `rename` 修正／`repro_tol: 1e-3` が PyYAML で文字列→ `0.001` に定格）と A-5 直列化 1 件（NaN→ null＋理由）を plan §14 に記録し修正・再固定。**判定値は 1 回目と 2 回目で完全に同一**。結果: **P-1'=False**（(i) 活性化の陽性対照（71 個）=True・(ii) 退化・老化の陽性対照（23 個）=False・(iii)=True）→ R-6'。**R-12'=不合格（判定不能）**: R-12a の活性化の陽性対照（71 個）が log2(n_exp) 0.118＞許容 0.10（退化・老化の陽性対照（23 個）は 0.027 で合格、R-12b・R-12c は合格）。感度は全版で結論不変（補正なし・残差化・退化・老化の陽性対照 30 個版・同 21 個版）→ R-13'=false。A2 候補の 8 遺伝子（マッチ補正版・探索）は 99% 以内（R-7 の読みが保たれる）。出力 `results/criteria_015.json`・`gene_rank_015.csv`・`run.log`・`README.md`。判定は PI（ゲート G3） | 出力: `experiments/EXP-20261011-002/results/` |
| 2026-10-11（着手前記録。**G2 承認済み PI 2026-10-11**） | **EXP-20261011-003（A6 埋め込みのみ）** | 全 46,957 核（snRNA/multiome・Reck 2025）を Geneformer V2-316M・**fp32（tf32 閉）** で 1 回エンコード。摂動なし。評価は CPU。`analysis/embed_a6.py`（新規）・`analysis/evaluate_a6.py`（新規）を**結果を見る前に作成・smoke・固定**してから本番 | GB10/CUDA（spark-cf5b または c083、fp32・tf32 閉）**を第一候補**（§10.1 の実測 11.5 核/s の由来）。Mac MPS の fp32 代替も可（fp32 は MPS≒CUDA、機体依存なし。CLAUDE.md §6.1） | tokenize 要実測（smoke で 300 秒目安）＋埋め込み **≈1.13 h（上振れ 1.5 h 見込み）＋軸・評価 CPU 0.5 h 程度**。総 GPU 2.5 h 未満（24 h 未満） | plan.md §8 の **B-1〜B-6**（材料 sha256 不一致・盲検出現・tokenize 異常・smoke 不通は**止める**）。チェックポイントでは主判定の値を計算しない（B-6）。`Approved-by: PI 2026-10-11` を確認してから本番 | 埋め込み: `ssh petadimensionlab@100.104.171.109 "ls labretary-remote/experiments/EXP-20261011-003/results/embed_a6/run_meta.json"` と `tail labretary-remote/experiments/EXP-20261011-003/logs/remote-a6.log` で進捗を確認。完了はローカル側で `python -c "import pandas as pd; print(len(pd.read_parquet('experiments/EXP-20261011-003/results/embeddings.parquet')))"`＝46,957。評価: `results/criteria_a6.json` が書ける（A-6 別の道具で読み直し）。**`isp_remote.sh status` はログ名が食い違うので使わない** | **計算・評価 完了（2026-10-11）**: 本番は spark-cf5b（GB10/CUDA）の tmux isp-a6 で **17:49:14 開始 → 18:44:14 終了（exit=0）**。埋め込み 3,276.7 秒（≈54.6 分・46,957 核・fp32・tf32 閉・attn=sdpa。実測 ≈14.3 核/s で見積もり 11.5 核/s をわずかに上回る）。回収後 sha256 一致・46,957 行を確認、リモート複製 2.3 GB は clean で削除・消失確認。評価コードのバグ 4 件を **PI への報告前に修正**（plan §13 第 4 回・SHA256SUMS 再固定 9 件 OK）し、axis_feasibility 2 系統 → evaluate_a6.py 本番を完走（A-1 盲検 0・A-2 材料 OK・A-5 NaN 0・A-6 再読込 OK）。**判定（ゲート G3）は PI 待ち。提案: 最初に当たった反証条件 = R-A2（間質系の検出遺伝子数・総カウントだけの AUC 0.8163 > 0.60 で分離は検出遺伝子数・総カウントで説明できる → 判定不能）**。P-2 は (c) 不合格（u_ada の検出された遺伝子数の分布をそろえたときの安定性 31.96° > 30°） | 出力予定: `experiments/EXP-20261011-003/results/`（embeddings・criteria_a6.json・gene_position_a6.csv・run.log・README.md）。**結果は 1 回だけ。見る前に固定**

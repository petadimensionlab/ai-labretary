# EXP-20260928-002 進行ログ（experiment-planner）

- 対象仮説: HYP-20260928-01（G1 承認済み, PI 2026-09-28）
- 方針: PI の判断（2026-09-28）に従い、案 A を実施。対照を上位集合にして案 B・C を事後の部分集合で評価する。
- G2: 不要（判定条件 `hours > 24 or n_cells > 10000 or g == "all"` のいずれにも当たらない）

## 2026-09-28
- [済] guardrails 50 件（91 行）を全件確認
- [済] HYP-20260928-01、EXP-20260928-001 の plan.md / config.yaml / select_genes_001.py / run_isp.py を確認
- [済] ディレクトリ作成（templates/experiment/ をひな型に）
- [作業中] 経路セット候補（S1〜S3）の作成と、対照の上位集合・部分集合の定義
- [済] `work/tokenized.dataset` を EXP-20260928-001 から `cp -c` で複製（cache-*.arrow は複製していない）。
  3 ファイルとも sha256 が EXP-20260928-001 の事前登録と一致
- [済] `analysis/pathway_candidates_002.tsv`（S1〜S3 の候補 406 個）を作成。摂動の結果は一切見ていない
- [済] `analysis/select_genes_002.py` を実行。**合計 1,136 遺伝子**（経路 120 / 筋線維芽細胞マーカー 24 /
  恒常発現の対照遺伝子 3 個 / 中立対照パネル 7 個 / 参照 2 個 / 対照 980）、延べ 681,600 細胞×遺伝子
  - gene_stats は EXP-20260928-001 と 17,465 行・8 列が一致（最大差 4.4e-16 = CSV の丸め）
  - 候補 406 → S4 除外 22 / 本データに無い 14 → 373 → S5(2) で 185 → S5(3) 185 → S6 で先頭 120
  - 内訳 S1 61 / S2 13 / S3 46。S5(4)（使える対照 20 以上）は差し替え 0 回、最小 275
  - 入力の長さ（検出された遺伝子数にほぼ等しい）の分布をそろえた後の目標セントロイド 588 / 600 / 587 細胞（seed 51 / 52 / 53）
  - 経路セットに盲検・ラベル定義・元の主対象は 0 件（検算済み）
- [注意] 一般機構（プロテアソーム・ユビキチン系）部分セットは 5 個しか入らなかった
  （S6 の辞書順の帰結）。HYP §4.3 の「一般機構 37 個」は成立しないので、F-4 は
  「一般機構 5 個を除いた 115 個でも主予測に合格するか」に置き換える（判定基準ではなく限界の記述）
- [追加指示] PI 2026-09-28: 計算機を DGX 2 台 + Mac に分ける。**seed 単位でのみ分割**する
- [済] `analysis/make_configs_002.py` で設定 10 本を生成（正本 config.yaml ＋ 機体ごと 3 本 ＋ 段階 6 本）
  ＋ 試し実行用 `work/config-smoke.yaml`（事前登録の対象外）＋ 機体一致確認用 `work/config-machinecheck.yaml`
- [済] **`run_isp.py --dry-run` が 12 本すべて通った**。`needs_G2_approval: false`（全件）
  - 正本: forward_passes 681,600 / 実測モデル 21.96 時間。1 seed あたり 7.32 時間
  - `config-seed*` の `WARN seeds が 3 未満` は機体分割の帰結（実験全体では 3 seed）
- [済] `analysis/power_002.py`（最小検出可能効果量を**事前登録の規則そのもの**で計算。guardrails 2026-09-28）
  - K=120: 99 パーセンタイル 0.213、実効の閾値 **0.245**（絶対値条件が効く）、偶然合格率 0.0018
  - **最小検出可能効果量 0.34（効果が均一）／0.44（遺伝子間 SD 2.051 込み）**
  - nonmarker K=76: 0.268 / 0.39 / **0.51**、筋線維芽細胞マーカー K=24: 0.483 / 0.69 / 0.92
  - **復元あり 0.213 対 復元なし 520 版 0.214**（差 0.001 未満）→ 有限母集団の補正はどの版でも効かない
- [済] `analysis/evaluate_002.py`（主判定 980 版／副次 600・520 版／偽経路セット 50 組／
  復元あり・なしの両方を計算／欠測を不合格と数えない）、`merge_results_002.py`、`machine_check_002.py`
- [済] plan.md 執筆（§4 セットの設計 / §6 判定基準 / §7 最小検出可能効果量 / §8 機体の割り当てと段階実行と中断規則 /
  §11 コードの一覧（判定に使う値の出所をたどった） / §13 限界 / §14 答えが出ないこと）
- [済] learnings 5 件を `logs/learnings-2026-09-28.jsonl` に python の json.dumps で追記し、
  全 25 行を読み直して検証（壊れた行 0）
- [作業中] 試し実行（smoke）を `work/config-smoke.yaml` で実行中 → `logs/smoke.log`、出力は `results-smoke/`
- [済] **試し実行（smoke）が最後まで通った**（16:44〜17:20, Mac/MPS/bf16）
  - `smoke: min_cells_per_gene 120 -> 5` がログに出て、**skip 0 件・摂動 41 回**
  - 41 遺伝子（target 3 / positive 24 / housekeeping 3 / negctrl_panel 7 / reference_myof 2 / negative 2）
  - **`emb/{orig,pert}.npy` が 41/41 で実在**（`ls` で確認、形は (50, 1152) float32）
  - `cells.csv` 2,050 行、状態のセントロイド 3 本（各 1,152 次元）、‖Δ‖ 中央値 0.120
  - **本番の `results/` は空のまま**（`output.results_dir: results-smoke` が効いた）
- [済] **事前登録を固定した**: `analysis/freeze_preregistration_002.sh` →
  `work/SHA256SUMS-preregistered.txt`（**38 ファイル、`shasum -c` 全件 OK、不一致 0**）
- [作業中] Mac の補助計算: `analysis/extra_centroids.py` を seed 51 / 52 / 53 で実行中
  → `work/extra_centroids-seed{S}.npz`（**摂動より先に済ませる**。PI 指示 2026-09-28）
- [残] **本番は未実行。** 順序は 補助計算 → 機体の一致確認（H-1〜H-3）→ 本番（seed 単位で 3 機体）
- [済] **Mac の補助計算（追加の状態のセントロイド）が 3 seed とも完了**（17:22〜17:36、13 分/seed）
  - `work/extra_centroids-seed{51,52,53}.npz`（MYOF_neg / FIB_org_unmatched、各 1,152 次元）
  - MYOF_neg 569〜574 細胞 / 33〜34 提供者、**出発点のセントロイドとの細胞の重なり 0**（付け替えが効いている）
  - FIB_org_unmatched 665 細胞 / 64 提供者
- [済] plan.md §8.4 を更新し、変更ログを書いて**再固定**（38 / 38 一致、不一致 0）
- [完了] **EXP-20260928-002 の設計と事前登録は完了。本番は未実行。**
  次にやること: (1) §8.3 の機体の一致確認を 2 台で回す（0.31 時間/機体）→ (2) H-1〜H-3 が合格したら本番

## 2026-09-28（simulation-engineer, 段階 1: 機体の一致確認）
- [済] `shasum -c work/SHA256SUMS-preregistered.txt` 38/38 OK（実行前の再確認）
- [済] `--dry-run`（Mac / config-machinecheck.yaml）OK。`needs_G2_approval: false`
- [済] 両機体 ping 確認（spark-cf5b, spark-c083、hostname で確認）
- [済] push 完了（run_isp.py / capture_env.sh / config-machinecheck.yaml / tokenized.dataset / モデル重み、
  sha256 一致を両機体で確認）
- [済] 追加ファイル `work/state_rows-seed51.json` を両機体へ手動で送付、sha256 一致確認（自動送付対象外。plan.md §8.7）
- [済] 両機体で `--dry-run`（LABRETARY_DEVICE=cuda）OK。`needs_G2_approval: false`（48 遺伝子×200 細胞×seed51、
  0.35 時間/機体の見積もり）
- [開始 17:43:28 ごろ] tmux で両機体同時実行開始。出力先 `results-machinecheck/`（両機体とも同名。回収時に
  `results-machinecheck-cf5b` / `results-machinecheck-c083` に分けて pull する）
- [完了 18:03 ごろ] 両機体とも exit=0、エラーなし。実測 0.33 h（cf5b）/ 0.325 h（c083）（見積もり 0.35 h 以内）
- [済] 回収: `isp_remote.sh pull` は config が `work/` 直下にあるため `exp_dir()` の解決がずれ失敗
  （`output.experiment_dir` と `dirname(CONFIG)` が食い違う。機体一致確認専用の既知の制約として記録）。
  手動 rsync で `results-machinecheck-cf5b/` / `results-machinecheck-c083/` に回収
- [済] 照合: `config.sha256` 両機体一致、48/48 遺伝子、`emb/{orig,pert}.npy` 48/48 実在（両機体）、
  遮断遺伝子 10 個の混入 0（両機体の `cells.csv` で確認）
- [済] `analysis/machine_check_002.py --compare` 実行
  → **H-1〜H-3 すべて合格**。系統成分の絶対値 0.0、最大絶対差 0.0（基準 ≤0.02 / ≤0.05）、
  細胞ごと Spearman ρ = 1.0（基準 ≥0.99, n=9,600）、状態のセントロイドコサイン 最小 1.0（基準 ≥0.9999）。
  **2 機体の埋め込み・cells.csv は sha256 まで完全一致**（バイト同一）。
  両機体 GPU は同一チップ（NVIDIA GB10、driver 595.58.03 / 580.173.02）、bf16 決定的計算のため一致したと考えられる
  （念のため device.json のホスト名が別であること、cells.csv・emb がプル元ごとに別ファイルであることを個別に確認し、
  取り違えでないことを検算済み）
- [判断] **§8.3 の一致確認は合格。本番へ進む（対処 (a)(b) は不要）。**
- [開始 18:08:17] 本番: spark-cf5b が `config-seed51.yaml`（tmux `isp-config-seed51_`, 見積もり 7.32 h）
- [開始 18:08:35] 本番: spark-c083 が `config-seed52.yaml` → 完了後に自動で `config-seed53.yaml`
  （tmux `isp-c083-seed52-53` で連結。見積もり計 14.6 h）
- 事前に両機体で外部ファイル送付・sha256 照合済み（`config-seed{51,52,53}.yaml`、
  `state_rows-seed{51,52,53}.json`）、両機体で `--dry-run` OK（`needs_G2_approval: false`、
  cf5b 227,200 細胞×遺伝子・7.32h、c083 各 seed 同様）
- [監視] A-1〜A-5 の中断規則に従い、最初の seed が終わった時点で G-3・M1〜M4 を確認してから継続する。
  バックグラウンドで進捗を監視中
- [完了 2026-09-29 01:13:12] **c083 の seed 52 が最初に完了**（18:08:35〜01:13:12, 実測 7.08 h, exit=0。
  見積もり 7.32 h の 1.5 倍＝10.98 h 未満で A-2 に抵触せず）。cf5b の seed 51 はまだ実行中
- [済] seed 52 を pull（`results-seed52/`, `cells.csv` 226,845 行, 1,136/1,136 遺伝子、config.sha256 一致、
  遮断遺伝子混入 0、`emb/{orig,pert}.npy` 1,136/1,136 実在）。`logs/remote-config-seed52.log` を回収、
  WARN 2 件のみ（device 上書きと `seeds<3` の想定内警告。`skip` 0 件 → A-4 に抵触せず）
- [済] **A-1 中断規則チェック（seed 52 単独）**: G-3（筋線維芽細胞マーカー24 が筋線維芽細胞の軸で合格）
  set_mean_z=1.448 > 閾値 0.51 → **合格**。M1（校正）n_valid=979, frac_hi=0.092, frac_lo=0.098（基準 ≤0.2）→
  **合格**。M4（相関）uniform ρ=-0.225, placebo ρ=-0.148（基準 ≤0.4 / ≤0.5）→ **合格**。
  M2・M3 は設計上 2 seed 以上の再現を要求するため、この 1 seed だけでは合否を確定できない
  （参考値のみ記録。`work/a1_check_seed52.json`）。**A-1 の判断: 継続してよい**
  （参考: 経路セット単独では G-1 の最終判定基準ではないが、seed 52 だけで見た目安の値は
  set_mean_z=0.139 で 1 seed 分の閾値 0.212 を下回っていた。**これは 3 seed 平均で見る最終判定とは別物**で、
  A-1 の対象にも入っていない。結果を見て設計を変えることはしない）
- 残り: cf5b seed 51 完了待ち、c083 seed 53（自動連結で実行中）完了待ち
- [完了 2026-09-29 01:05:42] **cf5b の seed 51 も完了**（18:08:22〜01:05:42, 実測 6.96 h, exit=0。
  見積もり 7.32 h 以内、A-2 に抵触せず）
- [済] seed 51 を pull・照合（`results-seed51/`, cells.csv 226,845 行, 1,136/1,136 遺伝子、config.sha256 一致、
  遮断遺伝子混入 0、emb 1,136/1,136 実在、WARN 2 件のみ・skip 0 件）
- [失敗→気づき] 監視用の自作スクリプト（`wait_all_done.sh`）が seed51 のログファイル名を
  `remote-config-seed51.log` として ssh 越しに tail していたが、`isp_remote.sh` の `tag()` は
  `basename ... .yaml | tr -c 'A-Za-z0-9_-' '_'` で **末尾の改行を `_` に変換する**ため、
  実際のログ名は `remote-config-seed51_.log`（末尾に `_`）だった。存在しないパスを `tail ... 2>/dev/null`
  で握りつぶしていたため、seed51 が完了しているのに検知できず、**A-2（1.5 倍超過）を誤って発報**した。
  `isp_remote.sh status` 自体は正しい名前を解決するため誤りなし。**教訓: 自作の監視スクリプトで
  リモートのログパスを直接組み立てず、`isp_remote.sh status`（または `tag()` と同じ変換）を使うこと。**
  実際の計算には異常なし（cf5b は見積もり内で正常終了）。cf5b は `clean` 済み（ヒト由来データの複製削除を確認）
- [済] 監視スクリプトを修正（`tag()` 相当のファイル名解決に直す）。c083 の seed 53 完了待ちを再開
- [完了 2026-09-29 08:18:11] **c083 の seed 53 も完了**（01:13:13〜08:18:11, 実測 7.08 h。c083 の 2 seed 計
  14.16 h、見積もり 14.6 h 以内、A-2/A-3 に抵触せず）
- [済] seed 53 を pull・照合（1,136/1,136 遺伝子、config.sha256 一致、遮断遺伝子混入 0、
  emb 1,136/1,136 実在、WARN 2 件のみ・skip 0 件）
- [済] 両機体で `clean` 実行、`removed` を確認（ヒト由来データの複製は両機体とも消失）
- [済] `analysis/merge_results_002.py` で 3 seed を `results/` に統合。自動照合すべて合格
  （`results/merge_report.json`: 3 seed とも 1,136/1,136 遺伝子、役割内訳一致、embedding 欠け 0、
  遮断遺伝子混入 0、`cells.csv`/`summary.csv` の sha256 記録済み）
- [済] `analysis/evaluate_002.py` 実行（9 分 6 秒）→ `results/eval/criteria_002.json`
  - **欠測 1 遺伝子（対照 980 の ADAMTSL2, ENSG00000197859）× 3 seed**。使える対照が 7〜8 個で
    下限 20 未満のため判定不能（規則どおり「不合格」ではなく「欠測」と数えた）。経路セット 120 には
    含まれない（`missing_in_pathway_set: []`）
  - **G-1 主判定（980 版）: 不合格**。経路セット 120 の目標軸遺伝子セットの平均スコア z = 0.157 / 0.139 / 0.227
    （seed 51/52/53）、閾値 0.213 / 0.213 / 0.212（絶対値条件 0.245 も満たさない）。3 seed とも不合格
  - 600 版・520 版でも同じ z 帯・同じ閾値でいずれも不合格（**3 版で判定が一致 = true**）
  - 副次 S-1（nonmarker, K=76）も不合格（z = 0.058 / 0.082 / 0.132、閾値 0.268〜0.271）
  - G-2 交差: 経路セット（目標軸−筋線維芽軸）> 0 が 3 seed とも成立（0.179/0.145/0.288）、
    筋線維芽マーカー（筋線維芽軸−目標軸）> 0 も 3 seed とも成立（1.711/1.923/1.608）→ 合格
  - G-3 アッセイの妥当性: 筋線維芽マーカー 24 が筋線維芽軸で合格（z=1.421/1.448/1.391、閾値 0.51 台）
  - G-4 中立性: 恒常発現・中立パネルとも 3 seed 全部が 5〜95 パーセンタイル内 → 合格
  - G-5 偶然合格率（実測）0.001775（≤0.02）
  - G-6 偽経路セット 50 組・偽 nonmarker セット 50 組とも 0/50 合格 → 合格（経路セットの合格が
    位置・長さ・細胞数だけで説明されない設計であることを確認。もっとも本判定自体は不合格なので
    直接は効かない）
  - 指標の健全性 M1〜M4 すべて合格（`metric_pass: true`）
  - 結論の提案（スクリプトの自動生成）: 「主予測は棄却。ただし『効果がない』ではなく
    『最小検出可能効果量以上の効果は無い』としか読めない（F-1）」
- [ls 確認済み] `results/{cells.csv, summary.csv, machines.json, merge_report.json, eval/criteria_002.json,
  eval/gene_seed_002.csv, eval/set_null_002.csv}`、`results-seed{51,52,53}/`、
  `results-machinecheck-{cf5b,c083}/`、`work/machine_check_002.json`、`work/a1_check_seed52.json`

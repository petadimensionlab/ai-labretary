
## 2026-09-28 EXP-20260928-001 設計（experiment-planner）
- [開始] PI 指示: 陽性対照の確立をやり直す。摂動の実行はしない。
- [済] EXP-20260927-004/validation.md 通読、guardrails.md 69行まで通読
- [済] EXP-004 plan.md / config.yaml 通読、labretary-level.md 通読
- [重要] **run_isp.py は遺伝子を 1 個ずつ回す**（main の `for gene, role in plan: run_one(cfg, f_cells, [gene], ...)`）。
  `perturbation.combos` は Geneformer の anchor_gene 併用（combos=1 は anchor 必須）で、複数遺伝子同時摂動には使えない。
  Geneformer 側は genes_to_perturb にリストを渡すと「まとめて摂動」する（in_silico_perturber.py:218）が、
  run_isp.py がリストを渡していない。**A-1 には run_isp.py への機能追加（遺伝子セット）が必要。**
  さらに select_cells() も単一遺伝子前提（その遺伝子を含む細胞のみ抽出）。セットでは「全員含む細胞」か「和集合」かを決める必要あり。
- [済] experiments/EXP-20260928-001/ を templates/experiment/ から作成
- [作業中] B の調査（SYN3 と RASA1 の共通性質）: analysis/confound_probe_001.py
- [済] B の調査 第1弾（analysis/confound_probe_001.py → work/confound_probe_001.json）
  - 主指標 z は「遺伝子ごとの平均移動方向 d と固定軸 w の内積」でほぼ決まる（r^2 = 0.81〜0.85, Spearman 0.88〜0.91）
  - d の seed 間再現性 cos = 0.985（中央値、最小 0.908）→ ばらつきの 92.5% が遺伝子固有である理由が判明。方向は遺伝子ごとにほぼ決まっている
  - **SYN3 と RASA1 の移動方向の cos = 0.353。全ペア中央値 0.457 より低い。共通方向はない**
  - 上位 8 遺伝子の第1主成分寄与率 43.9%（全 106 では 53.0%）→ 上位は集まっていない
  - z と dnorm / coherence / lp / ll / relpos の相関はすべて |rho| <= 0.14 → **新しい交絡は見つからない**
- [済] 検出力（analysis/power_001.py → work/power_001.json）
  - 単一遺伝子の床: 2.73 標準偏差（変わらない）。偶然 5% の閾値は 1.40 標準偏差
  - **セット単位（K 遺伝子の z の平均）の最小検出可能効果量（minimum detectable effect）: K=10 で 0.54, K=15 で 0.42, K=20 で 0.355, K=24 で 0.32, K=30 で 0.27 標準偏差**
  - EXP-004 の候補 10 個のセット平均 z は 1.17/1.17/1.00（SYN3 を除く 9 個でも 0.94/0.92/0.71）→ K=24 なら余裕で検出可能
  - 軸を替えても床（標準偏差単位）は下がらない（目標の軸 1.40 / 筋線維芽細胞の軸 1.51）→ A-2 は床を下げる方法ではない
  - A-1（同時摂動）の伸びは飽和: k=5 で 1.33 倍、上限 1.48 倍。対照も同じ倍率を受けるので標準偏差単位の利得は A-3 と同じ sqrt(K)
- [決定] 主軸は A-3（検定の単位を遺伝子セットに変える）＋ 軸 × セットの交差設計。A-1・A-2 は採らない（理由は plan.md §5）
- [作業中] 遺伝子セットの選定（analysis/select_genes_001.py）。seed は 41/42/43
- [済] 遺伝子選定（analysis/select_genes_001.py）。合計 156 遺伝子、延べ 93,600 細胞×遺伝子
  - 目標マーカー 24（上位 12 ＋ 無作為 12）/ 筋線維芽細胞マーカー 24 / 恒常発現 3 / 中立パネル 7 / 参照 ACTA2・TAGLN 2 / 対照 96
  - EXP-004 で摂動済みの 106 遺伝子はマーカーセットから除外（独立性）
  - 出発点・目標・対照の状態から MYOF/cycMYOF を除いた（重心の重なりを断つ）。深さマッチ後 595/585/589 細胞
  - 使える対照遺伝子の最小 22（下限 20）
- [済] analysis/extra_centroids.py と select_genes_fib.py を EXP-20260928-001 に複製（事前登録の一覧に入れる）
- [済] config.yaml を生成（遺伝子リストは config に展開。guardrails 2026-09-27）
- [済] `run_isp.py --dry-run` が警告なしで通った（forward_passes 93,600、内蔵見積もり 3.38 時間）
  - 判定に使うのは実測モデルの 6.0 時間（小さいほうを採らない。guardrails 2026-09-27）
- [作業中] 小さな試し実行（smoke）をバックグラウンドで実行中 → logs/smoke.log
- [作業中] plan.md の執筆
- [済] plan.md 執筆（§2 B の調査 / §5 A-1,A-2,A-3 の評価 / §6.3 陽性対照の確立の定義 G-1〜G-6 / §9 検出力 / §10 失敗時の書き方 / §11 コードの一覧と実行順序）
- [済] learnings 4 件を logs/learnings-2026-09-28.jsonl に追記し、全行を JSON として読み直して検証（壊れた行 0）
- [済] evaluate_001.py の自己検査（set_null / sham_sets / compare_same_scale）が通った
- [待ち] smoke 実行中。終わったら SHA256SUMS-preregistered.txt を作って事前登録を固定する
- [発見・重要] smoke が 41 遺伝子すべて skip した。原因 2 件（どちらも run_isp.py の既存の挙動）:
  1. `--smoke` は max_cells_per_gene を 50 に落とすが **min_cells_per_gene は下げない**ので、
     min_cells_per_gene: 120 の設定では全遺伝子が skip され、**摂動が 1 回も実行されない**（smoke が smoke になっていない）
  2. `output.results_dir` を明示していると **`--smoke` でも同じ results/ に書く**
     （`base = cfg["output"]["results_dir"] or ("results-smoke" if smoke ...)`）。
     本番の results/ に 50 細胞の出力が残ると、再開ロジック（済んだ遺伝子は飛ばす）で
     **50 細胞と 200 細胞が混ざった結果**になりうる
  → 対応: results/ を消して作り直し、試し実行専用の `work/config-smoke.yaml`（事前登録の対象外）を作って回す。
     共有スクリプト run_isp.py は**自分では直さず**、Labretary に修正を提案する
- [済] config.yaml と work/genes_001.json の遺伝子集合を役割ごとに照合（6 役割すべて一致、合計 156、重複 0、遮断遺伝子の混入 0）
- [済] extra_centroids.py が読む設定・参照先（state_rows の rows / depth_match_to）が 3 seed すべてで揃っていることを静的に確認
- [修正] FIB_org_unmatched（感度分析用の重心）にも subclass_l2 の絞り込みを入れ、主判定と同じく MYOF/cycMYOF を除いた。dry-run 再実行 OK
- [済] guardrails に 2 件追記（smoke の min_cells_per_gene / results_dir）。learnings は 6 行、壊れた行 0
- [発見・重要] **重心から MYOF を除くだけでは足りなかった。** run_isp.py の start_matrix() は
  状態ラベルだけで摂動する細胞を選ぶので、出発点 5,748 細胞のうち 1,132（19.7%）が筋線維芽細胞のままだった。
  そのままだと筋線維芽細胞マーカーセットの細胞だけ MYOF に偏り、交差（G-2）とアッセイの妥当性（G-3）の比較が成り立たない。
  → 対応: work/tokenized.dataset の MYOF/cycMYOF 細胞の fib_tls_state を `MYOF_<元の状態>` に付け替えた
  （行順と input_ids は不変。state_rows の行番号が正しい状態を指し MYOF 0 件であることを 3 seed で検算）。
  出発点は 4,616 細胞・MYOF 0 件になった。config の extra_centroids.MYOF_neg も新しい状態名に更新。dry-run 再実行 OK
- [作業中] smoke3（付け替え後のデータで）実行中 → logs/smoke3.log
- [済] select_genes_001.py に base_states() を追加（状態ラベルの付け替えの前後どちらでも同じ結果になるようにした）。
  付け替え後のデータで再実行し、**出力 9 ファイルすべてが sha256 で一致**することを確認（再現性の担保）
- [済] **smoke3 が摂動まで到達した（skip 0 件）**。1 遺伝子目 FLVCR2 で
  emb/orig.npy・pert.npy（50 × 1152）、cells.csv（shift_goal / shift_start / shift_alt）、
  state_emb_vectors.npz（FIB_org / FIB_neg_TLShigh / FIB_neg_TLSlow）が正しく出た。移動量 ||Δ|| 中央値 0.124
  → 状態ラベルの付け替え後も run_isp.py の状態解決・細胞選抜・摂動・埋め込み保存が動くことを確認
- [済] **事前登録を固定した**: analysis/freeze_preregistration.sh → work/SHA256SUMS-preregistered.txt（24 ファイル、shasum -c 全件 OK）
  - 含めたもの: plan.md, config.yaml, analysis/*.py 6 件, work/ の中間ファイル 13 件, tokenized.dataset の 3 ファイル
  - 含めないもの: work/config-smoke.yaml（試し実行専用）、results-smoke/（試し実行の出力）
- [残] 本番の実行は PI の承認待ち。実行前に results/ が空であることと、work/ の試し実行の残り
  （cells-seed41-n50-d20 等）が本番の細胞選抜（n200-d20）と別名であることを確認すること
- [済] **smoke3 完了**（02:53:44）。skip 0 件 / WARN・ERROR 0 件。
  summary.csv 10 遺伝子（target / positive / housekeeping / negctrl_panel / reference_myof / negative の 6 役割）、
  cells.csv 500 行（1 遺伝子 50 細胞）。状態重心の抽出も摂動も最後まで動くことを確認した
- [済] 事前登録の照合: `shasum -c work/SHA256SUMS-preregistered.txt` = 24/24 OK、不一致 0 件
- [完了] EXP-20260928-001 の設計と事前登録は完了。**本番は未実行**（PI へ報告し、実行は承認後）

## EXP-20260928-002（HYP-20260928-01 の検定）— experiment-planner / 2026-09-28
- 詳細は `experiments/EXP-20260928-002/logs/progress.md`
- [済] guardrails 50 件を全件確認 → 設計に反映（復元ありの引き方・偽セット 50 組・最小検出可能効果量の較正・
  smoke の min_cells_per_gene・results_dir の分離・外部ファイルの数え上げ・欠測を不合格と数えない）
- [済] **合計 1,136 遺伝子**（経路 120 / 対照 980 / アッセイ対照 36）、延べ 681,600 細胞×遺伝子、**21.9 時間**
- [済] 対照は上位集合 980 で、入れ子の部分集合 600（案 B 相当）・520（案 C 相当）を事前に定義（`work/null_pool_002.json`）
- [済] **G2 は不要**（`needs_G2_approval: false` を全 12 設定の dry-run で確認）
- [済] 最小検出可能効果量を事前登録の規則そのもので計算: K=120 で **0.34（均一）/ 0.44（遺伝子間 SD 2.05 込み）**
- [済] 計算機の割り当て（PI 指示 2026-09-28）: spark-c083 が seed 52・53、spark-cf5b が seed 51、Mac が補助計算。
  **seed 単位でのみ分割**。本番前の機体一致確認（H-1〜H-3）の判定基準を事前登録
- [作業中] 試し実行（smoke）→ 終わり次第 `work/SHA256SUMS-preregistered.txt` で事前登録を固定
- [残] **本番は未実行。** Mac の補助計算（追加の状態重心 3 seed）と機体の一致確認を先に済ませること
- [済] 試し実行（smoke）合格（skip 0 件、摂動 41 回、emb 41/41 実在）、Mac の補助計算（追加の状態重心 3 seed）完了
- [済] **事前登録を固定**（`work/SHA256SUMS-preregistered.txt`、38 ファイル、`shasum -c` 全件 OK）
- [完了] EXP-20260928-002 の設計は完了。**本番は未実行**（次は機体の一致確認 → 本番）

## EXP-20260929-001（HYP-20260929-01, G1 承認 2026-09-29）— experiment-planner
- [済] 背景の読み込み: HYP-20260929-01 / EXP-20260928-002 の plan・validation（末尾の訂正を含む）/ EXP-20260928-001 の枠組み / guardrails 52 件 / labretary-level
- [済] 設計の骨子を確定（結果を見る前）:
  - 状態: 目標 `FIB_orgHigh`（950）/ 出発点 `FIB_orgLow`（18,099）/ 筋線維芽細胞 `MYOF_fib`（第 2 の軸）/ **対照の状態 `FIB_ecmHigh`**（第 3 の軸。目標と同じ作り方で細胞外基質プログラムの上位 5%。**追加の重心として作り、出発点の細胞プールは 18,099 のまま**にして CXCL13 の余裕を削らない）
  - seed: **61 / 62 / 63**（既出と重ならない）
  - 対照遺伝子 300 を主判定、301〜500 は入れ子の感度分析として事前登録
- [作業中] 入力の作成（`analysis/build_input_003.py`）→ トークン化 → 事前チェック P-A（トークン化後の数え直し）/ P-B（深さをそろえた後の AUC）
- [済] `work/input.h5ad`（21,769 細胞 = 線維芽細胞 19,049 + 筋線維芽細胞 2,720）→ トークン化（`work/tokenized.dataset`）。
      目標 950 / 70 提供者、出発点 18,099 / 79 提供者、検出遺伝子数の中央値 1590.5 対 1591.0 を**HYP §2.2 と完全一致で再現**
- [済] 対照の状態のプログラム（細胞外基質）7 個を確定: COL1A1, COL1A2, COL3A1, COL5A1, COL6A3, FBN1, SPARC（900 細胞 / 63 提供者。目標との重なり 0）
- [済] **事前チェック P-A（トークン化後の数え直し）合格 → GO**。`work/precheck_counts_003.json`
      - **LTB だけが不合格**（生カウント 128 → トークン化後 119 細胞。下限 120 を 1 細胞下回る）→ セットから除外、K = 25 → **24**
      - CXCL13 **127**（余裕 +7）/ CCL19 **133** / CCL21 200（上限）→ **3 遺伝子とも合格**。諮り直しには当たらない
      - 目標マーカー 7・筋線維芽細胞マーカー 24・恒常発現 3・中立パネル 7 はすべて合格
- [済] 深さをそろえた状態の行番号（`work/state_rows-seed{61,62,63}.json`）。目標 931 / 出発点 931（1:1 マッチ 100%）、対照の状態 798〜799、筋線維芽細胞 849〜858。4 本の重心で細胞の重なり 0
- [済] 遺伝子の選定（`work/genes_003.json` ほか）。対照候補 9,041 → 対照 300（主判定）/ 500（入れ子）。**摂動する遺伝子 367**
- [作業中] 事前チェック P-B（深さをそろえた後の提供者分割 AUC。**0.65 未満なら実験を行わない**）を Mac で実行中
- [済] **事前チェック P-B 合格**（`work/precheck_auc_extra_003.json`, `results-precheck/precheck/precheck.json`）
      - 深さをそろえた目標 931 対 出発点 931 の **提供者分割 AUC = 0.6803**（5 分割 SD 0.032。**閾値 0.65 以上 → F-6 に当たらない**）
      - 並べ替え検定 **p = 0.005**（< 0.01）。無作為化時の平均 0.529 / SD 0.017
      - **深さ 2 変数だけの AUC = 0.5905**（生カウントでの参照値 0.6025）。トークン長だけでは 0.4627（ほぼ情報なし）
- [済] 最小検出可能効果量（`work/power_003.json`）。事前登録する規則そのもので計算。**主判定 K=24 で 0.91 標準偏差**（効果が均一なら 0.69）。偶然合格率 0.0051（基準 0.02 以下）。絶対値条件 `set_abs_min = 0.46`
- [作業中] 評価コード `analysis/evaluate_003.py`（**3 本の軸すべてでセット検定**／偽セット 200 組＋中心と散らばりの検査）、設定の生成、dry-run、試し実行
- [済] **事前チェック P-C（比べられる対照遺伝子の数）**。対照 300 では IL34 の適格な対照が 12 / 19 / 22（下限 20 割れ）→ **事前に決めた規則どおり対照を 500 に増やした**（入れ子なので 300 版は感度分析に残る）。500 では最小 23 / 37 / 47 で合格。計算量 7.1 → 11.0 時間（上限 24 時間内、G2 不要）
- [済] **dry-run 7 設定すべて通過**（`needs_G2_approval: false`）。正本 11.0 時間 / seed 別 3.67 時間 / 機体一致確認 0.59 時間 / 試し実行 0.10 時間
- [済] **試し実行（smoke）合格**。11 遺伝子 × 50 細胞、**skip 0 件**、`emb/{orig,pert}.npy` **22 / 22 実在**（形 (50, 1152) float32）、`cells.csv` 550 行、移動量の中央値 0.157、本番の `results/` は空のまま
- [済] **Mac の補助計算（追加の状態重心）3 seed 完了**。対照の状態 798〜799 細胞 / 63 提供者、筋線維芽細胞 849〜858 / 66〜69、深さをそろえない目標 931 / 70。**出発点・目標の重心との重なりはすべて 0**
- [済] **事前登録を固定**（`work/SHA256SUMS-preregistered.txt`、**43 ファイル**、`shasum -c` 全件一致・不一致 0）
- [完了] **EXP-20260929-001 の設計は完了。摂動の本番は未実行**（PI 指示どおり）。次は (1) 機体の一致確認（87 遺伝子 × 200 細胞 × 2 台、約 0.59 時間/機体）→ (2) 本番（spark-cf5b: seed 61 / spark-c083: seed 62・63、並行で約 7.3 時間）

## 2026-09-30 着手前チェック（軸が立つかの判定）の新設 — experiment-planner
- [済] `experiments/EXP-20260929-001/validation.md` と guardrails 59 件（特に 2026-09-30 の 7 件）を読了
- [済] 過去 4 実験の重心ファイル（`work/extra_centroids-seed*.npz` + `results/run-*-L-1/state_emb_vectors.npz`）から
      **軸の長さを摂動なしで再計算**。validation.md §2.4 の値を小数第 5 位まで再現
- [発見] **軸の長さの比は 4 実験を分離しない**（0.442〜0.637）。うまくいった EXP-20260928-001 も 0.452〜0.498 で、
      PI 案の閾値 0.70 はこのプロジェクトで一度も到達していない → 閾値の作り直しが必要
- [発見] validation.md §2.4 の「14% 短くなった」は **量の取り違え**（28-001 は向きベクトル長、29-001 は重心間距離）。
      同じ量で比べると −3%〜−9%。**軸長では 4 回の失敗を説明できない**
- [済] 発現データからマーカー群のまとまり・検出率を実測（29-001 の h5ad と 27-002 の h5ad）
- [作業中] `scripts/analysis/axis_feasibility.py` の実装と、過去 4 実験への当てはめ
- [済] `scripts/analysis/axis_feasibility.py` を実装（`py_compile` 通過）。**摂動を 1 回も行わない**。
      合格 / 判定保留 / 不合格 の 3 値。終了コード 0 / 1 / 2（設定・入力のエラーは 3）
- [済] 過去 4 実験に当てはめ（設定と出力は `outputs/axis-feasibility/`）
      - EXP-20260927-004 **不合格**（K = 1）
      - EXP-20260928-001 **判定保留**（K = 24。弱点 4 件）… 「限定つき合格」と一致
      - EXP-20260928-002 **不合格**（K = 0）
      - EXP-20260929-001 **不合格**（K = 7 / 検出率 0.042 / まとまりと検出率が 2 つそろって低い の 3 つが独立に発火）
- [発見] **止める条件は「定義マーカー」ではなく「実際に検定に使う陽性対照セット」に当てる必要がある。**
      28-001 と 28-002 は同じ目標の定義（ケモカイン 3 個）で結果が違った。分離したのは K と陽性対照の検出率
- [済] `templates/experiment/plan.md` §5.5・`templates/experiment/config.yaml` の `axis_feasibility:` 節・
      `docs/用語集.md` §3.5・`CLAUDE.md` §6 に追記。**既存の実験の設定は変更していない**
      （EXP-20260929-001 の事前登録 43 ファイルの `shasum -c` は全件一致のまま）
- [済] learnings 5 件を `logs/learnings-2026-09-30.jsonl` に追記し、全 12 行を読み直して検証（不正 0 件）

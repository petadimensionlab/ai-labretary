#!/usr/bin/env python3
"""2 つの結果を比べる**前に**、実行条件が同じかを突き合わせる門。

## なぜこれがあるか（2026-10-04、PI 指摘）

> 比較検証するのみ、条件違うか確かめず、違って停まったりしてて、おまえら適当に仕事してるやろ？

Mac の機体照合が落ちたとき、Labretary は原因を「MPS と CUDA の数値の違い」と断定して止めた。
**切り分けを一度もしていなかった。** 実際に振ってみると:

| 条件 | 同じ Mac・同じ重み・同じ入力での最大差 |
|---|---|
| bf16 / CPU eager ↔ CPU sdpa | **0.3125**（装置は同じ。注意機構の実装だけが違う） |
| bf16 / CPU ↔ MPS | 0.219 |
| **fp32 / どの組み合わせでも** | **0.0000238** |

**原因は装置ではなく bf16 だった。** しかも当時は `transformers` も 4.46.3 と 4.57.6 で違っていた。
条件を記録していなかったので、後から突き合わせることもできなかった。

→ **数値を比べる前に、まずここを比べる。** 条件が違えば、値が違うのは当たり前であって発見ではない。

## 使い方

    .venv/bin/python scripts/analysis/compare_conditions.py <結果A> <結果B> [<結果C> ...]

## 終了コード（2026-10-04 に設計を直した）

| コード | 意味 | 呼び出し側がすべきこと |
|---|---|---|
| **0** | 判定に効く条件がすべて一致 | そのまま進む |
| **1** | **記録が無く、突き合わせられない** | **止める。** 条件が同じかどうか確かめようがない |
| **2** | 条件が違う。何がどう違うかは出力に書いてある | **止めない。** 数値そのものを比べ、差が出たらこの違いを原因の候補として調べる。差が出なければ、その違いはこの計算に効かなかったという観測として記録する |

**初版は 2 を 1 と同じ扱いにして止めていたが、これは誤りだった**（2026-10-04）。
この手続きが守ろうとしているのは「条件が違うせいで数値がずれていることに気づかないまま進むこと」である。
数値が 1 ビットも違わないのに、条件が違うというだけで止めるのは、守ろうとしている判断から外れている。
実際、spark-cf5b と spark-c083 は `torch` と `transformers` の版が違うが、
埋め込み・摂動に使った細胞・状態の重心・主指標のすべてで差 0.0 だった。
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# 値を動かす条件。**ここが違えば、埋め込みが違っても原因は分かっている**
DECISIVE = [
    ("dtype", "数値型"),
    ("numerics.attn_implementation", "注意機構の実装"),
    ("numerics.cuda_matmul_allow_tf32", "CUDA matmul の tf32"),
    ("numerics.float32_matmul_precision", "fp32 matmul の精度"),
    ("numerics.transformers", "transformers の版"),
    ("numerics.torch", "torch の版"),
]
MISSING = "（記録なし）"
# 記録はするが、それ自体では値を動かさない（fp32 なら機体は一致する）
INFORMATIONAL = [("resolved", "装置"), ("host", "機体"), ("numerics.cuda_device", "GPU")]


def dig(d: dict, path: str):
    cur = d
    for k in path.split("."):
        if not isinstance(cur, dict) or k not in cur:
            return MISSING
        cur = cur[k]
    return cur


def load(p: Path) -> dict:
    f = p / "device.json"
    if not f.exists():
        raise SystemExit(f"{f} がありません。この結果は条件を記録していないので、突き合わせられません")
    return json.loads(f.read_text())


def main(dirs: list[str], quiet: bool = False) -> int:
    paths = [Path(d) for d in dirs]
    recs = [load(p) for p in paths]
    labels = [p.name for p in paths]
    w = max(len(x) for x in labels + ["条件"]) + 2

    bad = []
    for group, title in ((DECISIVE, "判定に効く条件"), (INFORMATIONAL, "参考（fp32 なら値を動かさない）")):
        print(f"\n=== {title} ===")
        for key, name in group:
            vals = [dig(r, key) for r in recs]
            same = len({json.dumps(v, sort_keys=True) for v in vals}) == 1
            # **「記録なし」どうしを一致と見なさない**（2026-10-04）。記録していないものは
            # 突き合わせられていないのであって、同じだと確かめたわけではない。
            # 沈黙で通る検査は、ここで起きた失敗そのものを見逃す。
            missing = any(v == MISSING for v in vals)
            ok = same and not missing
            mark = "  " if ok else ("**記録なし（突き合わせ不能）**" if missing else "**違う**")
            if not ok and group is DECISIVE:
                bad.append((name, dict(zip(labels, vals)),
                            "記録がない" if missing else "値が違う"))
            print(f"{name:<28} {mark}")
            for lab, v in zip(labels, vals):
                print(f"    {lab:<{w}} {v}")

    print()
    missing = [b for b in bad if b[2] == "記録がない"]
    differ = [b for b in bad if b[2] == "値が違う"]
    if missing:
        print("**記録が無く、突き合わせられません。止めてください。**")
        for name, v, _ in missing:
            print(f"  - {name}: " + " / ".join(f"{k}={x}" for k, x in v.items()))
        return 1
    if differ:
        print("**条件が違います。ただしこれは止める理由ではありません。**")
        for name, v, _ in differ:
            print(f"  - {name}: " + " / ".join(f"{k}={x}" for k, x in v.items()))
        print("\n数値そのものを比べてください。差が出たら、上の違いを原因の候補として 1 つずつ確かめること。")
        print("差が出なければ、「この違いはこの計算に効かなかった」という観測として報告に書くこと。")
        print("どちらの場合も「機体が違うから」で片付けないこと。")
        return 2
    print("判定に効く条件はすべて一致しています。値の違いは条件の違いでは説明できません。")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+", help="比べる結果ディレクトリ（device.json を含むもの）")
    a = ap.parse_args()
    sys.exit(main(a.dirs))

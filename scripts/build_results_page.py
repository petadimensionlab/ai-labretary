#!/usr/bin/env python3
"""結果ファイルをそのまま閲覧できる 1 枚の HTML を作る（PI 指示 2026-09-26「結果そのものを閲覧できるように」）。

要約を書き起こさず、既存の .md / .csv / .log の中身をそのまま埋め込む。実験が進んだら再実行すれば最新になる。
  .venv/bin/python scripts/build_results_page.py
  → outputs/results-page.html（Artifact として公開して PI に共有する）
"""
from __future__ import annotations

import html
import json
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "results-page.html"

# (節の見出し, [(表示名, 相対パス), ...])
SECTIONS: list[tuple[str, list[tuple[str, str]]]] = [
    ("通しの報告", [
        ("研究報告（2026-09-26）", "outputs/reports/2026-09-26-腎臓リンパ組織と細胞老化-研究報告.md"),
    ]),
    ("はじめに", [
        ("用語集", "docs/用語集.md"),
        ("仮説トラッカー", "hypotheses/tracker.md"),
    ]),
    ("実験 001 — 修復不全の近位尿細管で老化軸と LTβR を削除", [
        ("計画", "experiments/EXP-20260926-001/plan.md"),
        ("検証", "experiments/EXP-20260926-001/validation.md"),
        ("レポート", "experiments/EXP-20260926-001/report.md"),
    ]),
    ("実験 002 — 一様な削除効果を差し引く指標の開発", [
        ("計画", "experiments/EXP-20260926-002/plan.md"),
        ("結果（本実行）", "experiments/EXP-20260926-002/results/eval/README.md"),
        ("判定基準の厳しさの検討", "experiments/EXP-20260926-002/results/eval/criteria-strictness.md"),
        ("結果（予備実行）", "experiments/EXP-20260926-002/optionB/results-B/eval/README.md"),
        ("判定表", "experiments/EXP-20260926-002/results/eval/criteria.csv"),
    ]),
    ("実験 003 — 追加学習済みモデルでの摂動", [
        ("計画", "experiments/EXP-20260926-003/plan.md"),
        ("結果（追加学習 対 事前学習）", "experiments/EXP-20260926-003/results-ft/eval/README.md"),
        ("判定表（追加学習）", "experiments/EXP-20260926-003/results-ft/eval/criteria.csv"),
        ("判定表（事前学習・同条件）", "experiments/EXP-20260926-003/results-pretrained-matched/eval/criteria.csv"),
    ]),
    ("実験 004 — 比較遺伝子と提供者単位の解析", [
        ("計画", "experiments/EXP-20260926-004/plan.md"),
        ("発現をそろえた比較遺伝子", "experiments/EXP-20260926-004/de_matched/README.md"),
        ("提供者単位の再解析", "experiments/EXP-20260926-004/donor_level/README.md"),
        ("提供者 29 人・34 遺伝子（評価ログ）", "experiments/EXP-20260926-004/expanded/dl_all29_g34.log"),
        ("提供者 29 人・44 遺伝子（評価ログ）", "experiments/EXP-20260926-004/expanded/dl_all29_g44.log"),
        ("提供者 14 人・34 遺伝子（評価ログ）", "experiments/EXP-20260926-004/expanded/dl_high14_g34.log"),
        ("判定表（提供者拡張）", "experiments/EXP-20260926-004/expanded/results-expanded/eval44.log"),
    ]),
    ("実験 005 — 陽性対照の見直し", [
        ("計画", "experiments/EXP-20260926-005/plan.md"),
        ("検証", "experiments/EXP-20260926-005/validation.md"),
        ("結果（削除型）", "experiments/EXP-20260926-005/results-delete/README.md"),
        ("結果（過剰発現型）", "experiments/EXP-20260926-005/results-oe/README.md"),
        ("判定表（削除型）", "experiments/EXP-20260926-005/results-delete/eval/criteria.csv"),
        ("判定表（過剰発現型）", "experiments/EXP-20260926-005/results-oe/eval/criteria.csv"),
    ]),
    ("仮説", [
        ("HYP-20260925-01 老化軸と LTβR", "hypotheses/HYP-20260925-01.md"),
        ("HYP-20260925-02 線維芽細胞", "hypotheses/HYP-20260925-02.md"),
        ("HYP-20260925-03 腫瘍内 B 細胞", "hypotheses/HYP-20260925-03.md"),
        ("HYP-20260925-04 腫瘍内 T 細胞", "hypotheses/HYP-20260925-04.md"),
        ("HYP-20260926-01 補正指標", "hypotheses/HYP-20260926-01.md"),
        ("HYP-20260926-02 追加学習モデル", "hypotheses/HYP-20260926-02.md"),
        ("HYP-20260926-03 比較遺伝子・提供者単位", "hypotheses/HYP-20260926-03.md"),
        ("HYP-20260926-04 陽性対照の見直し", "hypotheses/HYP-20260926-04.md"),
    ]),
    ("調査資料", [
        ("データセット一覧", "outputs/literature/2026-09-25-kidney-TLS-datasets.md"),
        ("状態を作れるかの確認", "outputs/literature/2026-09-25-tls-feasibility/feasibility.md"),
        ("陰性対照の候補", "outputs/literature/2026-09-26-negative-controls.md"),
        ("データの取得元", "data/raw/kidney_tls/SOURCES.md"),
    ]),
]


def collect() -> list[dict]:
    out = []
    for section, items in SECTIONS:
        docs = []
        for label, rel in items:
            p = ROOT / rel
            if not p.is_file():
                continue
            kind = {".md": "md", ".csv": "csv"}.get(p.suffix, "text")
            docs.append({
                "id": rel.replace("/", "-").replace(".", "-"),
                "label": label, "path": rel, "kind": kind,
                "bytes": p.stat().st_size,
                "body": p.read_text(encoding="utf-8", errors="replace"),
            })
        if docs:
            out.append({"section": section, "docs": docs})
    return out


PAGE = """<title>腎臓リンパ組織 研究記録</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans+JP:wght@400;500;600&family=IBM+Plex+Serif:wght@500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root {
  --paper:#F6F5F2; --surface:#FFFFFF; --ink:#16202B; --muted:#5E6C79; --line:#DFDCD5;
  --accent:#0B6E6E; --accent-soft:#E3EFEE; --amber:#B45309; --green:#15803D; --red:#B3261E;
  --sans:"IBM Plex Sans JP",system-ui,sans-serif;
  --serif:"IBM Plex Serif",Georgia,serif;
  --mono:"IBM Plex Mono",ui-monospace,monospace;
}
@media (prefers-color-scheme:dark){ :root:not([data-theme="light"]){
  --paper:#10151A; --surface:#171E25; --ink:#E4EAEF; --muted:#93A3B1; --line:#2A343D;
  --accent:#3FB8AC; --accent-soft:#16302F; --amber:#D98324; --green:#4ADE80; --red:#F87171;
  color-scheme:dark; } }
:root[data-theme="dark"]{
  --paper:#10151A; --surface:#171E25; --ink:#E4EAEF; --muted:#93A3B1; --line:#2A343D;
  --accent:#3FB8AC; --accent-soft:#16302F; --amber:#D98324; --green:#4ADE80; --red:#F87171;
  color-scheme:dark; }
body{background:var(--paper);color:var(--ink);font-family:var(--sans);line-height:1.75;margin:0}
.wrap{display:grid;grid-template-columns:260px minmax(0,1fr);gap:0;min-height:100%}
nav{border-right:1px solid var(--line);background:var(--surface);padding:24px 16px;
  position:sticky;top:env(safe-area-inset-top,0px);align-self:start;max-height:100vh;overflow-y:auto}
nav h1{font-family:var(--serif);font-size:1.05rem;margin:0 0 4px;line-height:1.4}
nav .sub{color:var(--muted);font-size:.74rem;margin:0 0 18px;font-family:var(--mono)}
nav .grp{font-size:.7rem;letter-spacing:.09em;text-transform:uppercase;color:var(--muted);
  margin:18px 0 6px;font-weight:600}
nav a{display:block;padding:5px 10px;border-radius:6px;color:var(--ink);text-decoration:none;
  font-size:.83rem;border-left:2px solid transparent}
nav a:hover{background:var(--accent-soft)}
nav a.on{background:var(--accent-soft);border-left-color:var(--accent);font-weight:600}
main{padding-block:32px;padding-inline:clamp(16px,4vw,56px);max-width:none;min-width:0}
.doc{display:none}.doc.on{display:block}
.meta{font-family:var(--mono);font-size:.72rem;color:var(--muted);border-bottom:1px solid var(--line);
  padding-bottom:10px;margin-bottom:26px;display:flex;gap:14px;flex-wrap:wrap}
.body{max-width:74ch}
.body h1{font-family:var(--serif);font-size:1.6rem;line-height:1.35;text-wrap:balance;margin:0 0 .6em}
.body h2{font-family:var(--serif);font-size:1.22rem;margin:2em 0 .5em;text-wrap:balance;
  padding-bottom:.25em;border-bottom:1px solid var(--line)}
.body h3{font-size:1.02rem;margin:1.6em 0 .4em;color:var(--accent)}
.body h4{font-size:.92rem;margin:1.3em 0 .3em}
.body p,.body li{font-size:.9rem}
.body code{font-family:var(--mono);font-size:.82em;background:var(--accent-soft);padding:.1em .35em;border-radius:4px}
.body pre{background:var(--surface);border:1px solid var(--line);border-radius:8px;padding:14px;
  overflow-x:auto}
.body pre code{background:none;padding:0;font-size:.78rem;line-height:1.6}
.body blockquote{border-left:3px solid var(--accent);margin:1em 0;padding:.1em 0 .1em 14px;color:var(--muted)}
.body a{color:var(--accent)}
.tbl{overflow-x:auto;margin:1em 0;border:1px solid var(--line);border-radius:8px;max-width:100%}
table{border-collapse:collapse;width:100%;font-size:.8rem;font-variant-numeric:tabular-nums;background:var(--surface)}
th,td{text-align:left;padding:7px 11px;border-bottom:1px solid var(--line);white-space:nowrap}
th{background:var(--accent-soft);font-weight:600;font-size:.74rem;position:sticky;top:0}
tr:last-child td{border-bottom:none}
td.num{font-family:var(--mono)}
.raw{white-space:pre;font-family:var(--mono);font-size:.74rem;line-height:1.6;background:var(--surface);
  border:1px solid var(--line);border-radius:8px;padding:14px;overflow-x:auto}
.note{background:var(--accent-soft);border-radius:8px;padding:12px 14px;font-size:.82rem;margin-bottom:24px}
@media (max-width:820px){
  .wrap{grid-template-columns:1fr}
  nav{position:static;max-height:none;border-right:none;border-bottom:1px solid var(--line);padding:16px}
  nav .links{display:flex;gap:6px;overflow-x:auto;padding-bottom:4px}
  nav .grp{margin:14px 0 4px}
  nav a{white-space:nowrap;border-left:none;border-bottom:2px solid transparent}
  nav a.on{border-left:none;border-bottom-color:var(--accent)}
  main{padding-inline:16px}
}
</style>
<div class="wrap">
<nav>
  <h1>腎臓リンパ組織 研究記録</h1>
  <p class="sub" id="stamp"></p>
  <div id="nav"></div>
</nav>
<main>
  <div class="note">各文書は生成時点のファイルの中身をそのまま表示しています。要約ではありません。</div>
  <div id="docs"></div>
</main>
</div>
<script src="https://cdnjs.cloudflare.com/ajax/libs/marked/4.3.0/marked.min.js"></script>
<script id="data" type="application/json">__DATA__</script>
<script>
const DATA = JSON.parse(document.getElementById('data').textContent);
document.getElementById('stamp').textContent = DATA.generated;
const nav = document.getElementById('nav'), docs = document.getElementById('docs');
function csvRows(text){
  const rows=[]; let row=[], cell='', q=false;
  for(let i=0;i<text.length;i++){ const c=text[i];
    if(q){ if(c==='"'){ if(text[i+1]==='"'){cell+='"';i++;} else q=false; } else cell+=c; }
    else if(c==='"') q=true;
    else if(c===','){ row.push(cell); cell=''; }
    else if(c==='\\n'){ row.push(cell); rows.push(row); row=[]; cell=''; }
    else if(c!=='\\r') cell+=c;
  }
  if(cell||row.length){ row.push(cell); rows.push(row); }
  return rows.filter(r=>r.some(x=>x!==''));
}
function renderCsv(text){
  const rows=csvRows(text); if(!rows.length) return '<p>（空）</p>';
  const head=rows[0], body=rows.slice(1);
  const esc=s=>s.replace(/[&<>]/g,m=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[m]));
  const isNum=s=>/^-?\\d*\\.?\\d+(e-?\\d+)?$/i.test(s.trim());
  return '<div class="tbl"><table><thead><tr>'+head.map(h=>'<th>'+esc(h)+'</th>').join('')+
    '</tr></thead><tbody>'+body.map(r=>'<tr>'+head.map((_,i)=>{
      const v=r[i]??''; return '<td'+(isNum(v)?' class="num"':'')+'>'+esc(v)+'</td>';
    }).join('')+'</tr>').join('')+'</tbody></table></div>';
}
let first=true;
for(const grp of DATA.sections){
  const g=document.createElement('div'); g.className='grp'; g.textContent=grp.section; nav.appendChild(g);
  const links=document.createElement('div'); links.className='links'; nav.appendChild(links);
  for(const d of grp.docs){
    const a=document.createElement('a'); a.href='#'+d.id; a.textContent=d.label; a.dataset.id=d.id;
    links.appendChild(a);
    const sec=document.createElement('section'); sec.className='doc'; sec.id=d.id;
    const kb=(d.bytes/1024).toFixed(1);
    let inner='';
    if(d.kind==='md') inner='<div class="body">'+marked.parse(d.body)+'</div>';
    else if(d.kind==='csv') inner=renderCsv(d.body);
    else inner='<div class="raw">'+d.body.replace(/[&<>]/g,m=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[m]))+'</div>';
    sec.innerHTML='<div class="meta"><span>'+d.path+'</span><span>'+kb+' KB</span></div>'+inner;
    for(const t of sec.querySelectorAll('.body table')){
      const w=document.createElement('div'); w.className='tbl';
      t.parentNode.insertBefore(w,t); w.appendChild(t);
    }
    docs.appendChild(sec);
    if(first){ show(d.id); first=false; }
  }
}
function show(id){
  for(const s of docs.children) s.classList.toggle('on', s.id===id);
  for(const a of nav.querySelectorAll('a')) a.classList.toggle('on', a.dataset.id===id);
  window.scrollTo({top:0});
}
nav.addEventListener('click', e=>{
  const a=e.target.closest('a'); if(!a) return;
  e.preventDefault(); show(a.dataset.id); history.replaceState(null,'','#'+a.dataset.id);
});
const h=location.hash.slice(1);
if(h && document.getElementById(h)) show(h);
</script>
"""


def main() -> None:
    sections = collect()
    n = sum(len(s["docs"]) for s in sections)
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    data = {"generated": f"{stamp} 生成 / {n} 文書", "sections": sections}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        PAGE.replace("__DATA__", json.dumps(data, ensure_ascii=False).replace("</", "<\\/")),
        encoding="utf-8")
    print(f"{OUT} ({OUT.stat().st_size/1024:.0f} KB, {n} 文書)")


if __name__ == "__main__":
    main()

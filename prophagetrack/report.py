"""Offline report with escaped data and simple SVG matrix/timelines."""
import html
from collections import Counter


def esc(x):return html.escape(str(x),quote=True)

def table(rows,fields):
    if not rows:return '<p>No entries.</p>'
    return '<div class="scroll"><table><thead><tr>'+''.join(f'<th>{esc(x)}</th>' for x in fields)+'</tr></thead><tbody>'+''.join('<tr>'+''.join(f'<td>{esc(r.get(f,""))}</td>' for f in fields)+'</tr>' for r in rows)+'</tbody></table></div>'

def render(path,samples,regions,events,genes,manifest,recovered):
    groups=sorted({r['exact_group'] for r in regions});sids=[s['sample_id'] for s in samples]
    found={(r['sample_id'],r['exact_group']) for r in recovered}
    display=groups[:60];width=230+22*len(display);height=70+24*len(sids[:80])
    svg=[f'<svg role="img" aria-label="Exact sequence presence matrix" viewBox="0 0 {width} {height}">']
    for j,g in enumerate(display):svg.append(f'<text x="{240+j*22}" y="22" font-size="10" text-anchor="middle">{j+1}</text>')
    for i,s in enumerate(sids[:80]):
        y=40+i*24;svg.append(f'<text x="5" y="{y+15}" font-size="12">{esc(s)}</text>')
        for j,g in enumerate(display):
            color='#087e8b' if (s,g) in found else '#e3e9ef'
            svg.append(f'<rect x="{230+j*22}" y="{y}" width="19" height="19" rx="3" fill="{color}"><title>{esc(s)} / {esc(g)}: {"matched" if (s,g) in found else "not detected; absence unproven"}</title></rect>')
    svg.append('</svg>')
    counts=Counter(e['event'] for e in events)
    body=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ProphageTrack report</title>
<style>body{{font:16px system-ui;color:#192d3b;background:#f3f6f8;margin:0}}main{{max-width:1200px;margin:auto;padding:38px}}h1{{font-size:38px;margin-bottom:8px}}h2{{margin-top:32px}}.note{{background:#fff1d6;padding:18px;border-left:4px solid #ba7b16}}.cards{{display:flex;gap:14px;flex-wrap:wrap}}.card{{background:white;border:1px solid #dbe3e8;padding:20px;border-radius:12px;min-width:160px}}.card b{{font-size:28px;display:block;color:#087e8b}}table{{border-collapse:collapse;width:100%;font-size:13px;background:white}}th,td{{padding:10px;border-bottom:1px solid #dde5eb;text-align:left}}th{{background:#dfeef0}}.scroll{{overflow:auto}}svg{{background:white;max-height:1500px;width:100%}}code{{background:#e3e9ef}}a{{color:#087e8b}}</style>
<main><p>RESEARCH PREVIEW · {esc(manifest['version'])}</p><h1>ProphageTrack</h1><p>Shared prophage sequences, genes, and longitudinal evidence.</p>
<div class="note">{esc(manifest.get('dataset_label','User-supplied data'))}. Predictions and candidate events are not proof of induction, infection, or transmission. Strain assignments are supplied by the user. Empty-site evidence is assembly-based.</div>
<div class="cards"><div class="card"><b>{len(samples)}</b>Assemblies</div><div class="card"><b>{len(regions)}</b>Predicted regions</div><div class="card"><b>{len(groups)}</b>Sequence groups*</div><div class="card"><b>{len(events)}</b>Transitions assessed</div></div>
<p>*Ambiguous-base regions have separate UNCERTAIN identifiers. Exact extracted sequences do not establish complete prophages.</p>
<h2>Exact sequence distribution</h2><p>Teal: full exact sequence recovered. Gray: not detected, not confirmed absence. Columns follow the group list below. Display limited to 60 groups and 80 samples; TSVs contain all data.</p>{''.join(svg)}
{table([dict(column=i+1,exact_group=g) for i,g in enumerate(display)],['column','exact_group'])}
<h2>Longitudinal events</h2><p>{esc(dict(counts))}</p>{table(events[:300],['subject_id','strain_id','from_sample','to_sample','event','reference_region'])}
<h2>Regions and prediction evidence</h2>{table(regions[:300],['region_id','sample_id','contig','start','end','callers','flags','exact_group'])}
<h2>Shared gene families</h2>{table([g for g in genes if g['n_regions']>1][:300],['gene_family','n_proteins','n_regions','regions'])}
<h2>Analysis scope</h2><p>Related-sequence comparison: {esc(manifest['parameters']['related'])}. Gene analysis: {esc(manifest['gene_status'])}. Boundary-uncertain or overlapping regions are excluded from event anchors.</p>
<p>Full results: <a href="regions.tsv">regions</a> · <a href="events.tsv">events</a> · <a href="site_observations.tsv">site evidence</a> · <a href="gene_families.tsv">gene families</a> · <a href="manifest.json">run provenance</a></p></main></html>'''
    path.write_text(body)

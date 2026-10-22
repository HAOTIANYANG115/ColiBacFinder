"""Presentation only: retain original screening rules and all detailed evidence."""
from collections import defaultdict
from pathlib import Path
from html import escape
from .io import read_tsv, write_tsv

STRONG = {'all_model_components_near_reference', 'modeled_components_near_reference',
          'core_near_reference_components_incomplete_or_ambiguous'}
FIELDS = ['sample', 'run_status', 'result', 'main_candidates', 'component_support',
          'immunity_candidates', 'weak_or_partial_loci', 'notes']
LABELS = ['Sample', 'Run status', 'Screening result', 'Main candidates', 'Component support', 'Immunity candidates', 'Weak or partial loci', 'Notes']

def build_rows(summaries, systems):
    by_sample = defaultdict(list)
    for r in systems:
        by_sample[r['sample']].append(r)
    rows = []
    for s in summaries:
        groups, weak = defaultdict(list), set()
        for r in by_sample[s['sample']]:
            anchor = r.get('anchor_locus', '')
            if not anchor:
                continue
            if r['status'] in STRONG:
                groups[(r.get('anchor_contig',''), anchor)].append(r)
            else:
                weak.add(anchor)
        # One physical anchor may have competing system names; never count these as separate genes.
        strong_anchors = {key[1] for key in groups}
        weak -= strong_anchors
        names, support = [], []
        notes = [s.get('warnings',''), s.get('error','')]
        for key, rr in sorted(groups.items()):
            name = '/'.join(sorted({r['system'] for r in rr}))
            names.append(name)
            for r in rr:
                support.append('{}: {}/{} [{}]'.format(r['system'],
                    r.get('nearby_near_components',0),r.get('modeled_components',0),r.get('model_scope','')))
            if len({r['system'] for r in rr}) > 1:
                notes.append('Multiple candidate types at the same locus; subtype unresolved')
        complete = s['status'] == 'COMPLETE'
        near_core = int(s.get('near_core_loci') or 0)
        near_immune = int(s.get('near_nonexport_immunity_loci') or 0)
        if not complete:
            result = 'Analysis failed; no conclusion'
            names, support = [], []
        elif near_core:
            result = 'Near-reference bacteriocin/precursor candidate detected'
            if not names:
                names = [s.get('near_core_candidate_labels','')]
        elif near_immune:
            result = 'Near-reference immunity candidate only'
        elif int(s.get('primary_candidate_loci') or 0):
            result = 'Weak or partial candidates only; review required'
        else:
            result = 'No primary candidate detected at current thresholds'
        if complete and near_core and not groups:
            notes.append('Bacteriocin/precursor candidate has no matching system model; inspect detailed evidence')
        rows.append(dict(sample=s['sample'],run_status=s['status'],result=result,
            main_candidates='; '.join(names),component_support='; '.join(support),
            immunity_candidates=s.get('near_immunity_candidate_labels','') if complete else '',
            weak_or_partial_loci=len(weak) if complete else '',
            notes='; '.join(x for x in notes if x)))
    return rows

def render(run_dir):
    run_dir = Path(run_dir)
    rows = build_rows(read_tsv(run_dir/'01_sample_summary.tsv'),read_tsv(run_dir/'06_system_candidates.tsv'))
    write_tsv(run_dir/'00_results_summary.tsv',FIELDS,rows)
    head = ''.join('<th>'+escape(x)+'</th>' for x in LABELS)
    body = ''.join('<tr>'+''.join('<td>'+escape(str(r.get(k,'')))+'</td>' for k in FIELDS)+'</tr>' for r in rows)
    doc = '<!doctype html><html lang="en"><meta charset="utf-8"><title>ColiBacFinder results</title>'
    doc += '<style>body{font:15px system-ui;margin:32px;color:#172b3a}table{border-collapse:collapse;width:100%}th,td{border:1px solid #ccd5dc;padding:10px;text-align:left;vertical-align:top;overflow-wrap:anywhere}th{background:#edf3f6}tr:nth-child(even){background:#f8fafb}p{max-width:1100px;line-height:1.6}</style>'
    doc += '<h1>ColiBacFinder results overview</h1><p>{} samples; {} completed. Genomic predictions do not establish bacteriocin production or antibacterial activity.</p>'.format(len(rows),sum(r['run_status']=='COMPLETE' for r in rows))
    doc += '<p>Main candidates require near-reference bacteriocin or precursor support. Names separated by a slash are alternatives at one locus, not separate genes. Component support counts refer to the current model, not necessarily all biologically required components. Coding integrity is not assessed. Multiple immunity labels can describe the same locus.</p>'
    doc += '<p>No detection applies only to the current reference database and thresholds. Review weak or partial hits in 03_candidate_loci.tsv, 04_locus_evidence.tsv and 06_system_candidates.tsv. All detailed evidence is retained.</p>'
    doc += '<table><thead><tr>'+head+'</tr></thead><tbody>'+body+'</tbody></table></html>'
    (run_dir/'00_results_report.html').write_text(doc,encoding='utf-8')
    return rows

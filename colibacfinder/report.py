from pathlib import Path
from collections import Counter
from .io import write_tsv, write_fasta, write_json, tokens
from .engine import make_loci, make_pairs, CORE
from .systems import make_systems, MODELS

HIT_FIELDS = ['sample','hsp_id','reference_id','labels','roles','reference_tiers',
    'reference_accessions','source_libraries','contig','start','end','strand','subject_frame',
    'identity_pct','reference_span_pct','alignment_length_aa','reference_start_aa',
    'reference_end_aa','evalue','bitscore','screen_class']
LOCUS_FIELDS = ['sample','locus_id','contig','start','end','strand','subject_frame',
    'locus_category','roles','candidate_labels','near_reference_labels','screen_classes',
    'best_reference_id','best_labels','best_reference_tiers','best_identity_pct',
    'best_reference_span_pct','best_bitscore','best_evalue','alternative_reference_count',
    'second_reference_id','bitscore_gap_to_second','relative_bitscore_gap','multi_label_candidate',
    'mixed_core_immunity','distance_to_sequence_boundary','near_assembly_edge','topology',
    'coding_integrity','functional_status']
PAIR_FIELDS = ['sample','contig','core_locus','immunity_locus','gap_bp','overlap_bp',
    'same_strand','crosses_origin','interpretation']
SYSTEM_FIELDS = ['sample','system','anchor_locus','anchor_contig','anchor_start','anchor_end',
    'status','model_scope','nearby_near_components','modeled_components','near_components',
    'missing_near_components','noncore_near_reference_HSPs','reference_tiers',
    'coding_integrity','functional_status','interpretation']
COMP_FIELDS = ['sample','system','anchor_locus','component','expected_labels','status',
    'location_class','reference_id','labels','reference_tiers','contig','start','end','strand',
    'identity_pct','reference_span_pct','gap_bp','overlap_bp','self_anchor','crosses_origin',
    'nearby_qualifying_reference_count']
SUMMARY_FIELDS = ['sample','status','contigs','total_bp','ambiguous_bp','all_HSPs',
    'primary_candidate_loci','near_reference_primary_loci','near_core_loci',
    'near_nonexport_immunity_loci','auxiliary_candidate_loci','proximity_pairs',
    'near_core_candidate_labels','near_immunity_candidate_labels','system_candidate_names',
    'all_model_components_candidate_names','mixed_role_loci','warnings','error']
TABLES = {
 'hits':('02_all_HSPs.tsv',HIT_FIELDS),
 'loci':('03_candidate_loci.tsv',LOCUS_FIELDS),
 'members':('04_locus_evidence.tsv',HIT_FIELDS+['locus_id','reference_rank','bitscore_gap_from_best']),
 'pairs':('05_proximity_pairs.tsv',PAIR_FIELDS),
 'systems':('06_system_candidates.tsv',SYSTEM_FIELDS),
 'components':('07_system_components.tsv',COMP_FIELDS),
}
REVCOMP = str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')

def extract_context(out, sample, loci, contigs, circular, rules):
    records, coords = [], []
    for r in loci:
        if r['locus_category'] != 'primary':
            continue
        seq = contigs[r['contig']]
        n = len(seq)
        start, end = r['start']-1-rules.flank_bp, r['end']+rules.flank_bp
        if r['contig'] in circular:
            # No duplicate sequence when requested neighborhood exceeds one revolution.
            size = min(n,end-start)
            start %= n
            part = (seq+seq)[start:start+size]
            finish = start+size
            intervals = f'{start+1}..{min(finish,n)}' + (f';1..{finish-n}' if finish>n else '')
        else:
            start, end = max(0,start), min(n,end)
            part = seq[start:end]
            intervals = f'{start+1}..{end}'
        if r['strand']=='-':
            part = part.translate(REVCOMP)[::-1]
        record_id = f'{sample}__{r["locus_id"]}'
        # FASTA ID must not inherit spaces from a filename.
        safe_id = r['locus_id']
        records.append((safe_id, part))
        coords.append(dict(sample=sample,locus_id=r['locus_id'],fasta_id=safe_id,
            contig=r['contig'],source_intervals_1based=intervals,
            reverse_complement=r['strand']=='-',length_bp=len(part),
            topology=r['topology']))
    write_fasta(out/'candidate_neighborhoods.fna', records)
    write_tsv(out/'neighborhood_coordinates.tsv',
        ['sample','locus_id','fasta_id','contig','source_intervals_1based','reverse_complement','length_bp','topology'],coords)

def analyze(sample, hits, contigs, circular, rules, out):
    loci, members = make_loci(hits,contigs,circular,rules)
    pairs = make_pairs(loci,contigs,circular,rules)
    systems, components = make_systems(sample,loci,members,hits,contigs,circular,rules)
    data = dict(hits=hits,loci=loci,members=members,pairs=pairs,systems=systems,components=components)
    for key,(name,fields) in TABLES.items():
        write_tsv(out/name,fields,data[key])
    extract_context(out,sample,loci,contigs,circular,rules)
    primary = [r for r in loci if r['locus_category']=='primary']
    core, immunity = set(),set()
    core_labels, immune_labels = set(),set()
    for h in members:
        if h['screen_class']!='near_reference' or 'export' in tokens(h['roles']):
            continue
        if tokens(h['roles']) & CORE:
            core.add(h['locus_id']); core_labels.update(tokens(h['labels']))
        if 'immunity' in tokens(h['roles']):
            immunity.add(h['locus_id']); immune_labels.update(tokens(h['labels']))
    warnings = []
    ambiguous = sum(sum(c not in 'ACGT' for c in seq) for seq in contigs.values())
    if ambiguous:
        warnings.append('ambiguous_nucleotides_present')
    if len(contigs)>1000:
        warnings.append('high_contig_count_review_input')
    summary = dict(sample=sample,status='COMPLETE',contigs=len(contigs),
        total_bp=sum(map(len,contigs.values())),ambiguous_bp=ambiguous,all_HSPs=len(hits),
        primary_candidate_loci=len(primary),near_reference_primary_loci=sum(bool(r['near_reference_labels']) for r in primary),
        near_core_loci=len(core),near_nonexport_immunity_loci=len(immunity),
        auxiliary_candidate_loci=len(loci)-len(primary),proximity_pairs=len(pairs),
        near_core_candidate_labels='|'.join(sorted(core_labels)),
        near_immunity_candidate_labels='|'.join(sorted(immune_labels)),
        system_candidate_names='|'.join(sorted({r['system'] for r in systems if r['anchor_locus']})),
        all_model_components_candidate_names='|'.join(sorted({r['system'] for r in systems if r['status']=='all_model_components_near_reference'})),
        mixed_role_loci=sum(r['mixed_core_immunity'] for r in primary),warnings=';'.join(warnings),error='')
    return summary,data

def aggregate(out, summaries, results):
    write_tsv(out/'01_sample_summary.tsv',SUMMARY_FIELDS,summaries)
    for key,(name,fields) in TABLES.items():
        write_tsv(out/name,fields,(r for data in results for r in data[key]))
    label_samples = {}
    for data in results:
        for h in data['hits']:
            if h['screen_class']=='near_reference':
                for label in tokens(h['labels']):
                    label_samples.setdefault(label,set()).add(h['sample'])
    n = sum(s['status']=='COMPLETE' for s in summaries)
    write_tsv(out/'08_label_frequency.tsv',['label','samples','completed_denominator','percent'],
        [dict(label=k,samples=len(v),completed_denominator=n,percent=round(100*len(v)/n,3) if n else '')
        for k,v in sorted(label_samples.items(),key=lambda kv:(-len(kv[1]),kv[0]))])
    write_json(out/'system_models.json', MODELS)
    write_json(out/'run_status.json',dict(status='COMPLETE' if n==len(summaries) else 'COMPLETED_WITH_FAILURES',
        samples=len(summaries),completed=n,failed=len(summaries)-n))

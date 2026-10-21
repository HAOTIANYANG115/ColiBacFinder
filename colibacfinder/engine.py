"""Single-HSP classification and fixed-anchor locus grouping.
No HSP chaining, gene calling, protein subtype validation or activity inference.
"""
from dataclasses import dataclass, asdict
from .io import tokens

CORE = {'bacteriocin', 'precursor'}
PRIMARY = CORE | {'immunity'}

@dataclass
class Rules:
    search_evalue: float = 1e-3
    screen_evalue: float = 1e-5
    identity: float = 80.0
    coverage: float = 80.0
    remote_identity: float = 35.0
    remote_coverage: float = 70.0
    remote_length: int = 40
    partial_identity: float = 35.0
    partial_coverage: float = 10.0
    partial_length: int = 60
    overlap_fraction: float = 0.5
    neighbor_bp: int = 10000
    flank_bp: int = 5000
    edge_bp: int = 200
    def validate(self):
        for name in ['identity','coverage','remote_identity','remote_coverage',
                     'partial_identity','partial_coverage']:
            if not 0 < getattr(self, name) <= 100:
                raise ValueError(f'{name} must be in (0,100]')
        if not 0 < self.overlap_fraction <= 1:
            raise ValueError('overlap_fraction must be in (0,1]')
        if not 0 < self.screen_evalue <= self.search_evalue:
            raise ValueError('Require 0 < screen_evalue <= search_evalue')
        if min(self.neighbor_bp, self.flank_bp, self.edge_bp) < 0:
            raise ValueError('Distances must be nonnegative')
    def to_dict(self):
        return asdict(self)

def classify(identity, coverage, length, evalue, rules):
    if evalue <= rules.screen_evalue:
        if identity >= rules.identity and coverage >= rules.coverage:
            return 'near_reference'
        if identity >= rules.remote_identity and coverage >= rules.remote_coverage and length >= rules.remote_length:
            return 'remote'
        if identity >= rules.partial_identity and coverage >= rules.partial_coverage and length >= rules.partial_length:
            return 'partial'
    return 'below_screen_threshold'

def parse_raw(path, db, sample, contigs, id_map, rules):
    hits = []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) != 13:
                raise ValueError('Expected 13 TBLASTN columns')
            rid, sid = p[:2]
            if rid not in db.sequences or sid not in id_map:
                raise ValueError(f'Unknown BLAST sequence ID: {rid}, {sid}')
            contig = id_map[sid]
            identity, aln, qlen, slen = float(p[2]), int(p[3]), int(p[4]), int(p[5])
            qa, qb, sa, sb = map(int, p[6:10])
            ev, score, frame = float(p[10]), float(p[11]), int(p[12])
            if (qlen != len(db.sequences[rid]) or slen != len(contigs[contig])
                or not 1 <= qa <= qb <= qlen or not 1 <= min(sa,sb) <= max(sa,sb) <= slen
                or frame not in {-3,-2,-1,1,2,3}):
                raise ValueError('Invalid BLAST coordinates or lengths')
            coverage = 100.0 * (qb-qa+1) / qlen
            hits.append(dict(
                sample=sample, hsp_id=f'HSP{len(hits)+1:07d}', reference_id=rid,
                labels='|'.join(db.values(rid, 'reporting_label')),
                roles=';'.join(sorted(db.roles(rid))),
                reference_tiers='|'.join(db.values(rid, 'reference_tier')),
                reference_accessions='|'.join(db.values(rid, 'protein_accession')),
                source_libraries='|'.join(db.values(rid, 'source_library')),
                contig=contig, start=min(sa,sb), end=max(sa,sb),
                strand='+' if frame > 0 else '-', subject_frame=frame,
                identity_pct=identity, reference_span_pct=coverage,
                alignment_length_aa=aln, reference_start_aa=qa, reference_end_aa=qb,
                evalue=ev, bitscore=score,
                screen_class=classify(identity, coverage, aln, ev, rules)))
    return hits

def overlap(a, b):
    return max(0, min(a['end'],b['end']) - max(a['start'],b['start']) + 1)

def distance(a, b, length=None):
    """Return minimum interval gap, overlap, whether origin crossing was shorter."""
    ov = overlap(a,b)
    linear = max(0, max(a['start'],b['start']) - min(a['end'],b['end']) - 1)
    if length and not ov:
        first, last = sorted([a,b], key=lambda r:r['start'])
        wrap = length - last['end'] + first['start'] - 1
        if wrap < linear:
            return wrap, 0, True
    return linear, ov, False

def hit_order(h):
    return (-h['bitscore'], h['evalue'], -h['reference_span_pct'],
            -h['identity_pct'], h['reference_id'], h['contig'], h['start'], h['end'])

def group_hits(hits, rules):
    # Separate export-only/dual export-immunity evidence from primary locus counts.
    # Two non-overlapping copies never merge; anchor never expands to bridge loci.
    groups = []
    for h in sorted(hits, key=hit_order):
        if h['screen_class'] == 'below_screen_threshold':
            continue
        primary = bool(tokens(h['roles']) & PRIMARY) and 'export' not in tokens(h['roles'])
        for group in groups:
            a = group[0]
            ap = bool(tokens(a['roles']) & PRIMARY) and 'export' not in tokens(a['roles'])
            if (h['contig'], h['subject_frame'], primary) != (a['contig'], a['subject_frame'], ap):
                continue
            ov = overlap(h, a)
            if ov >= rules.overlap_fraction*(h['end']-h['start']+1) and ov >= rules.overlap_fraction*(a['end']-a['start']+1):
                group.append(h)
                break
        else:
            groups.append([h])
    return sorted(groups, key=lambda g:(g[0]['contig'], min(h['start'] for h in g), g[0]['subject_frame']))

def make_loci(hits, contigs, circular, rules):
    loci, members = [], []
    for i, group in enumerate(group_hits(hits, rules), 1):
        # Unique-reference ranking: a second HSP of the same reference is not a second candidate type.
        by_ref = {}
        for h in sorted(group, key=hit_order):
            by_ref.setdefault(h['reference_id'], h)
        ranked = list(by_ref.values())
        best = ranked[0]
        near = [h for h in ranked if h['screen_class'] == 'near_reference']
        roles = set().union(*(tokens(h['roles']) for h in group))
        primary = bool(roles & PRIMARY) and 'export' not in roles
        start, end = min(h['start'] for h in group), max(h['end'] for h in group)
        rid = f'LOCUS{i:05d}'
        second = ranked[1] if len(ranked)>1 else None
        edge = min(start-1, len(contigs[best['contig']])-end)
        labels = sorted(set().union(*(tokens(h['labels']) for h in ranked)))
        near_labels = sorted(set().union(*(tokens(h['labels']) for h in near))) if near else []
        row = dict(sample=best['sample'], locus_id=rid, contig=best['contig'],
            start=start, end=end, strand=best['strand'], subject_frame=best['subject_frame'],
            locus_category='primary' if primary else 'auxiliary',
            roles=';'.join(sorted(roles)), candidate_labels='|'.join(labels),
            near_reference_labels='|'.join(near_labels),
            screen_classes=';'.join(sorted({h['screen_class'] for h in group})),
            best_reference_id=best['reference_id'], best_labels=best['labels'],
            best_reference_tiers=best['reference_tiers'],
            best_identity_pct=best['identity_pct'], best_reference_span_pct=best['reference_span_pct'],
            best_bitscore=best['bitscore'], best_evalue=best['evalue'],
            alternative_reference_count=len(ranked)-1,
            second_reference_id=second['reference_id'] if second else '',
            bitscore_gap_to_second=best['bitscore']-second['bitscore'] if second else '',
            relative_bitscore_gap=(best['bitscore']-second['bitscore'])/best['bitscore'] if second and best['bitscore'] else '',
            multi_label_candidate=len(labels)>1,
            mixed_core_immunity=bool(roles & CORE) and 'immunity' in roles,
            distance_to_sequence_boundary=edge,
            near_assembly_edge=(edge <= rules.edge_bp and best['contig'] not in circular),
            topology='circular_declared' if best['contig'] in circular else 'linear_or_unknown',
            coding_integrity='not_assessed', functional_status='not_tested')
        loci.append(row)
        ranks = {h['reference_id']:j+1 for j,h in enumerate(ranked)}
        for h in group:
            members.append(dict(h, locus_id=rid, reference_rank=ranks[h['reference_id']],
                                bitscore_gap_from_best=best['bitscore']-h['bitscore']))
    return loci, members

def make_pairs(loci, contigs, circular, rules):
    cores = [r for r in loci if r['locus_category']=='primary' and tokens(r['roles']) & CORE and not r['mixed_core_immunity']]
    immune = [r for r in loci if r['locus_category']=='primary' and 'immunity' in tokens(r['roles']) and not r['mixed_core_immunity']]
    rows = []
    for a in cores:
        for b in immune:
            if a['contig'] != b['contig'] or a['locus_id']==b['locus_id']:
                continue
            gap, ov, wrap = distance(a,b,len(contigs[a['contig']]) if a['contig'] in circular else None)
            if gap <= rules.neighbor_bp:
                rows.append(dict(sample=a['sample'],contig=a['contig'],
                    core_locus=a['locus_id'], immunity_locus=b['locus_id'], gap_bp=gap,
                    overlap_bp=ov, same_strand=a['strand']==b['strand'], crosses_origin=wrap,
                    interpretation='proximity_only_not_cognate_pair_confirmation'))
    return rows

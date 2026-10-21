"""Versioned reporting models, not experimentally calibrated operon definitions."""
from .engine import distance, hit_order
from .io import tokens

# Complete-component wording is restricted to the explicitly enumerated V/J25/B17 models.
MODELS = {
 'Microcin_V': dict(scope='four_component_model', complete=True,
    core=['cvaC'], components={'precursor':['cvaC'],'immunity':['cvi'],
                             'export_A':['cvaA'],'export_B':['cvaB']}),
 'Microcin_J25': dict(scope='four_component_model', complete=True,
    core=['mcjA'], components={'precursor':['mcjA'],'maturation_B':['mcjB'],
                             'maturation_C':['mcjC'],'export_immunity':['mcjD']}),
 'Microcin_B17': dict(scope='seven_component_model', complete=True,
    core=['mcbA'], components={x:[x] for x in ['mcbA','mcbB','mcbC','mcbD','mcbE','mcbF','mcbG']}),
 'Microcin_H47': dict(scope='catalog_components_only', complete=False,
    core=['mchB'], components={'precursor':['mchB'],'immunity':['mchI'],
        'maturation_A':['mchA','mcmL'],'maturation_S1':['mchS1','mcmK'],
        'maturation_C':['mchC'],'maturation_D':['mchD'],'export_E':['mchE'],'export_F':['mchF']}),
 'Microcin_M': dict(scope='precursor_immunity_only', complete=False,
    core=['mcmA'], components={'precursor':['mcmA'],'immunity':['mcmI']}),
 'Microcin_I47': dict(scope='precursor_immunity_only', complete=False,
    core=['mchS2'], components={'precursor':['mchS2'],'immunity':['mchS3']}),
 'Colicin_M': dict(scope='bacteriocin_immunity_only', complete=False,
    core=['Colicin-M'], components={'bacteriocin':['Colicin-M'],
    'immunity':['ColM_immunity','Cmi_like_immunity']}),
 'Colicin_I_family': dict(scope='Ia_Ib_unresolved_family_pair', complete=False,
    core=['Colicin-Ia','Colicin-Ib'], components={'bacteriocin':['Colicin-Ia','Colicin-Ib'],
    'immunity':['ColIa Immunity_type_unverified','Colicin-Ib immunity protein']}),
}

# Additional type-labelled pairs remain candidate associations, never subtype confirmation.
PAIR_LABELS = {
 'A': (['Colicin-A'], ['ColA_immunity','ColA_immunity_N_terminal_extension']),
 'B': (['ColB_type_unverified'], ['ColB Immunity_type_unverified']),
 'D': (['Colicin-D'], ['ColD Immunity_type_unverified']),
 'E1': (['Colicin-E1'], ['Colicin-E1 immunity protein']),
 'E2': (['Colicin-E2'], ['ColE2_immunity','ColE2_like_immunity']),
 'E3': (['Colicin E3'], ['Colicin E3 immunity protein']),
 'E4': (['ColE4_C_terminal_fragment'], ['ColE4_immunity_candidate']),
 'E5': (['ColE5_type_unverified'], ['ColE5_immunity','ColE5_immunity_containing_sequence']),
 'E6': (['Colicin-E6'], ['ColE6_immunity']),
 'E7': (['Colicin-E7'], ['ColE7_immunity']),
 'E8': (['ColE8_type_unverified'], ['ColE8_immunity']),
 'E9': (['Colicin-E9'], ['Colicin-E9 immunity protein']),
 'K': (['Colicin-K'], ['ColK_Col5_shared_immunity']),
 '5': (['Colicin-5'], ['ColK_Col5_shared_immunity']),
 '10': (['Colicin-10'], ['Col10 Immunity_type_unverified']),
 'S4': (['ColS4_type_unverified'], ['ColS4 Immunity_type_unverified']),
 'R': (['ColR_type_unverified'], ['ColR Immunity_type_unverified']),
 'Y': (['ColY_type_unverified'], ['ColY Immunity_type_unverified']),
}
for name, (core, immunity) in PAIR_LABELS.items():
    MODELS['Colicin_'+name] = dict(scope='candidate_label_pair_only', complete=False,
        core=core, components={'bacteriocin':core,'immunity':immunity})
MODELS['Cloacin'] = dict(scope='candidate_label_pair_only', complete=False,
    core=['Cloacin'], components={'bacteriocin':['Cloacin'],
        'immunity':['Cloacin_immunity','Cloacin_like_immunity','Cloacin_family_immunity']})
MODELS['Colicin_Z'] = dict(scope='core_only_no_immunity_model', complete=False,
    core=['ColZ_type_unverified'],components={'bacteriocin':['ColZ_type_unverified']})


def matches(h, labels):
    return bool(tokens(h['labels']) & set(labels))


def make_systems(sample, loci, members, hits, contigs, circular, rules):
    systems, components = [], []
    for system, model in MODELS.items():
        groups = {}
        for h in members:
            if matches(h,model['core']):
                groups.setdefault(h['locus_id'], []).append(h)
        # Every system gets an explicit not-detected row if no screened core anchor exists.
        if not groups:
            n = sum(h['screen_class']=='near_reference' and any(
                matches(h, labels) for labels in model['components'].values()) for h in hits)
            systems.append(dict(sample=sample, system=system, anchor_locus='',
                status='no_screened_core_candidate', model_scope=model['scope'],
                nearby_near_components=0, modeled_components=len(model['components']),
                noncore_near_reference_HSPs=n, coding_integrity='not_assessed',
                functional_status='not_tested'))
            continue
        for locus, core_hits in groups.items():
            # Prefer screened near-reference evidence for component composition.
            core_near = [h for h in core_hits if h['screen_class']=='near_reference']
            anchor = sorted(core_near or core_hits, key=hit_order)[0]
            near_names, missing, selected = [], [], []
            all_model_loci = set()
            for component, labels in model['components'].items():
                choices = [h for h in hits if matches(h, labels) and h['evalue']<=rules.screen_evalue]
                nearby = []
                for h in choices:
                    if h['contig'] != anchor['contig']:
                        continue
                    gap, ov, wrap = distance(h,anchor,len(contigs[h['contig']]) if h['contig'] in circular else None)
                    if gap <= rules.neighbor_bp:
                        nearby.append(h)
                eligible = [h for h in nearby if h['screen_class']=='near_reference']
                candidate = sorted(eligible or nearby or choices, key=hit_order)
                h = candidate[0] if candidate else None
                row = dict(sample=sample,system=system,anchor_locus=locus,component=component,
                    expected_labels='|'.join(labels), status='no_HSP_at_screen_evalue',
                    nearby_qualifying_reference_count=len({x['reference_id'] for x in eligible}))
                if h:
                    same = h['contig']==anchor['contig']
                    gap, ov, wrap = distance(h,anchor,len(contigs[h['contig']]) if h['contig'] in circular else None) if same else ('','',False)
                    is_nearby = same and gap <= rules.neighbor_bp
                    is_self = matches(h,model['core']) and overlap_same(h,anchor)
                    row.update(reference_id=h['reference_id'],labels=h['labels'],
                        reference_tiers=h['reference_tiers'],contig=h['contig'],
                        start=h['start'],end=h['end'],strand=h['strand'],
                        identity_pct=h['identity_pct'],reference_span_pct=h['reference_span_pct'],
                        status=h['screen_class'],location_class='near_anchor' if is_nearby else 'elsewhere',
                        gap_bp='' if is_self else gap,overlap_bp='' if is_self else ov,
                        self_anchor=is_self,crosses_origin=wrap)
                    selected.append(h)
                if eligible:
                    near_names.append(component)
                    for m in members:
                        if m['hsp_id']==eligible[0]['hsp_id']:
                            all_model_loci.add(m['locus_id'])
                else:
                    missing.append(component)
                components.append(row)
            # A single ambiguous locus cannot fill multiple distinct component slots.
            complete = not missing and len(all_model_loci) >= len(model['components'])
            fragment = all('fragment' in tokens(h['reference_tiers']) for h in core_hits)
            if fragment:
                status='fragment_reference_anchor'
            elif complete:
                status='all_model_components_near_reference' if model['complete'] else 'modeled_components_near_reference'
            elif not core_near:
                status='weak_or_partial_core_candidate'
            else:
                status='core_near_reference_components_incomplete_or_ambiguous'
            systems.append(dict(sample=sample,system=system,anchor_locus=locus,
                anchor_contig=anchor['contig'],anchor_start=anchor['start'],anchor_end=anchor['end'],
                status=status,model_scope=model['scope'],nearby_near_components=len(near_names),
                modeled_components=len(model['components']),
                near_components=';'.join(near_names),missing_near_components=';'.join(missing),
                reference_tiers=';'.join(sorted(set().union(*(tokens(h['reference_tiers']) for h in selected)))),
                coding_integrity='not_assessed',functional_status='not_tested',
                interpretation='candidate_association_not_subtype_or_function_confirmation'))
    return systems, components

def overlap_same(a,b):
    return a['contig']==b['contig'] and a['subject_frame']==b['subject_frame'] and max(a['start'],b['start'])<=min(a['end'],b['end'])

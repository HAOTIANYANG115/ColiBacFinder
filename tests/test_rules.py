import tempfile
import unittest
from pathlib import Path
import gzip
import shutil
from colibacfinder.database import Database
from colibacfinder.engine import Rules, classify, group_hits, make_loci, make_pairs, parse_raw
from colibacfinder.systems import make_systems, MODELS
from colibacfinder.io import fasta, discover, tokens
from colibacfinder.report import analyze

DB = Database()
RULES = Rules()

def hit(label, start, end, frame=1, cls='near_reference', score=200, contig='c'):
    rid = next(rid for rid in DB.meta if label in DB.values(rid,'reporting_label'))
    return dict(sample='test',hsp_id=f'{label}_{start}_{score}',reference_id=rid,
        labels='|'.join(DB.values(rid,'reporting_label')),roles=';'.join(sorted(DB.roles(rid))),
        reference_tiers='|'.join(DB.values(rid,'reference_tier')),reference_accessions='',
        source_libraries='|'.join(DB.values(rid,'source_library')),
        contig=contig,start=start,end=end,strand='+' if frame>0 else '-',subject_frame=frame,
        identity_pct=100 if cls=='near_reference' else 40, reference_span_pct=100 if cls=='near_reference' else 30,
        alignment_length_aa=100,reference_start_aa=1,reference_end_aa=100,
        evalue=1e-20,bitscore=score,screen_class=cls)

class Tests(unittest.TestCase):
    def test_database_counts(self):
        self.assertEqual(len(DB.sequences),84)
        self.assertEqual(sum(map(len,DB.meta.values())),92)

    def test_model_labels_all_exist(self):
        labels=set().union(*(set(DB.values(r,'reporting_label')) for r in DB.meta))
        for model in MODELS.values():
            for names in model['components'].values():
                self.assertTrue(set(names)<=labels,names)

    def test_threshold_boundaries(self):
        self.assertEqual(classify(80,80,10,1e-5,RULES),'near_reference')
        self.assertEqual(classify(79.9,80,40,1e-5,RULES),'remote')
        self.assertEqual(classify(99,20,60,1e-5,RULES),'partial')
        self.assertEqual(classify(99,20,59,1e-5,RULES),'below_screen_threshold')
        self.assertEqual(classify(100,100,100,1.1e-5,RULES),'below_screen_threshold')

    def test_shared_reference_count(self):
        h=hit('ColM_immunity',1,423)
        groups=group_hits([h,dict(h,hsp_id='second')],RULES)
        loci,_=make_loci([h,dict(h,hsp_id='second')],{'c':'A'*1000},set(),RULES)
        self.assertEqual(len(groups),1)
        self.assertEqual(loci[0]['alternative_reference_count'],0)

    def test_two_candidates_one_locus(self):
        hh=[hit('Colicin-Ia',1,1800,score=200),hit('Colicin-Ib',1,1800,score=190)]
        loci,members=make_loci(hh,{'c':'A'*2000},set(),RULES)
        self.assertEqual(len(loci),1)
        self.assertTrue(loci[0]['multi_label_candidate'])
        self.assertEqual(loci[0]['bitscore_gap_to_second'],10)
        self.assertEqual(len(members),2)

    def test_distinct_copies_remain_separate(self):
        self.assertEqual(len(group_hits([hit('cvaC',1,309),hit('cvaC',1000,1308)],RULES)),2)

    def test_fixed_anchor_not_chain(self):
        hs=[hit('cvaC',1,100,score=300),hit('cvaC',51,150,score=200),hit('cvaC',101,200,score=100)]
        self.assertEqual(len(group_hits(hs,RULES)),2)

    def test_different_frames_not_merge(self):
        self.assertEqual(len(group_hits([hit('cvaC',1,300,1),hit('cvi',1,300,2)],RULES)),2)

    def test_microcin_overlap_retained(self):
        hs=[hit('cvi',141,374,3),hit('cvaC',355,663,1)]
        contigs={'c':'A'*4000}
        loci,_=make_loci(hs,contigs,set(),RULES)
        pairs=make_pairs(loci,contigs,set(),RULES)
        self.assertEqual(len(pairs),1)
        self.assertEqual(pairs[0]['overlap_bp'],20)
        self.assertEqual(pairs[0]['gap_bp'],0)

    def test_complete_V_model(self):
        hs=[hit('cvi',141,374,3),hit('cvaC',355,663,1),hit('cvaA',1000,2271,-1),hit('cvaB',2300,4393,-2)]
        seq={'c':'A'*6000}
        loci,members=make_loci(hs,seq,set(),RULES)
        systems,components=make_systems('test',loci,members,hs,seq,set(),RULES)
        row=next(r for r in systems if r['system']=='Microcin_V')
        self.assertEqual(row['status'],'all_model_components_near_reference')
        row=next(r for r in components if r['system']=='Microcin_V' and r['component']=='precursor')
        self.assertEqual(row['overlap_bp'],'')
        self.assertTrue(row['self_anchor'])

    def test_immunity_only_not_system(self):
        hs=[hit('cvi',1,234)]
        seq={'c':'A'*2000}; loci,members=make_loci(hs,seq,set(),RULES)
        systems,_=make_systems('test',loci,members,hs,seq,set(),RULES)
        self.assertEqual(next(r for r in systems if r['system']=='Microcin_V')['status'],'no_screened_core_candidate')

    def test_weak_exporter_not_complete(self):
        hs=[hit('cvaC',355,663),hit('cvi',141,374,3),hit('cvaA',1000,2000),hit('cvaB',2500,4000,cls='below_screen_threshold')]
        seq={'c':'A'*6000}; loci,members=make_loci(hs,seq,set(),RULES)
        systems,_=make_systems('test',loci,members,hs,seq,set(),RULES)
        row=next(r for r in systems if r['system']=='Microcin_V')
        self.assertIn('incomplete',row['status'])
        self.assertIn('export_B',row['missing_near_components'])

    def test_nearby_weak_cannot_hide_nearby_near(self):
        hs=[hit('cvaC',355,663),hit('cvi',141,374,3),hit('cvaA',1000,2000),
            hit('cvaB',2500,4000,cls='below_screen_threshold',score=300),hit('cvaB',4500,6000,score=200)]
        seq={'c':'A'*8000}; loci,members=make_loci(hs,seq,set(),RULES)
        systems,_=make_systems('test',loci,members,hs,seq,set(),RULES)
        self.assertEqual(next(r for r in systems if r['system']=='Microcin_V')['status'],'all_model_components_near_reference')

    def test_other_contig_component_not_near(self):
        hs=[hit('cvaC',355,663),hit('cvi',141,374,3),hit('cvaA',1000,2000),hit('cvaB',2500,4000,contig='d')]
        seq={'c':'A'*6000,'d':'A'*6000}; loci,members=make_loci(hs,seq,set(),RULES)
        systems,_=make_systems('test',loci,members,hs,seq,set(),RULES)
        self.assertIn('incomplete',next(r for r in systems if r['system']=='Microcin_V')['status'])

    def test_one_locus_cannot_supply_two_components(self):
        hs=[hit('cvaC',1,309),hit('cvi',500,734),hit('cvaA',1000,2200),hit('cvaB',1000,2200)]
        seq={'c':'A'*4000}; loci,members=make_loci(hs,seq,set(),RULES)
        systems,_=make_systems('test',loci,members,hs,seq,set(),RULES)
        self.assertIn('ambiguous',next(r for r in systems if r['system']=='Microcin_V')['status'])

    def test_fragment_not_complete_system(self):
        hs=[hit('ColE4_C_terminal_fragment',1,531),hit('ColE4_immunity_candidate',550,804)]
        seq={'c':'A'*2000}; loci,members=make_loci(hs,seq,set(),RULES)
        systems,_=make_systems('test',loci,members,hs,seq,set(),RULES)
        self.assertEqual(next(r for r in systems if r['system']=='Colicin_E4')['status'],'fragment_reference_anchor')

    def test_circular_pair(self):
        hs=[hit('Colicin-M',19000,19812),hit('ColM_immunity',20,442)]
        seq={'c':'A'*20000}; loci,_=make_loci(hs,seq,set(),RULES)
        self.assertEqual(make_pairs(loci,seq,set(),RULES),[])
        p=make_pairs(loci,seq,{'c'},RULES)[0]
        self.assertEqual(p['gap_bp'],207)
        self.assertTrue(p['crosses_origin'])

    def test_auxiliary_not_primary(self):
        with tempfile.TemporaryDirectory() as d:
            s,_=analyze('test',[hit('mcjD',1,1740)],{'c':'A'*2000},set(),RULES,Path(d))
            self.assertEqual(s['primary_candidate_loci'],0)
            self.assertEqual(s['near_nonexport_immunity_loci'],0)

    def test_invalid_fasta(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.fa'
            for text in ['abc\n','>x\nACGT\n>x\nACGT\n','>x\n','>x\nACGT*\n']:
                p.write_text(text)
                with self.assertRaises(ValueError):fasta(p,True)

    def test_fasta_gzip(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.fasta.gz'
            with gzip.open(p,'wt') as f:f.write('>a hello\nACGTNN\n')
            self.assertEqual(fasta(p,True),{'a':'ACGTNN'})
            self.assertEqual(discover(Path(d))[0][0],'x')

    def test_duplicate_sample_names(self):
        with tempfile.TemporaryDirectory() as d:
            for name in ['x.fa','x.fasta']:(Path(d)/name).write_text('>a\nACGT\n')
            with self.assertRaises(ValueError):discover(Path(d))

    def test_checksum_failure(self):
        with tempfile.TemporaryDirectory() as d:
            target=Path(d)/'db';shutil.copytree(DB.path,target)
            with (target/'reference_proteins.faa').open('a') as f:f.write('\n')
            with self.assertRaises(ValueError):Database(target)

    def test_bad_rules(self):
        with self.assertRaises(ValueError):Rules(identity=101).validate()
        with self.assertRaises(ValueError):Rules(screen_evalue=1).validate()

    def test_blast_id_mapping(self):
        rid=next(iter(DB.sequences));qlen=len(DB.sequences[rid])
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'hits';p.write_text(f'{rid}\tCBF1\t100\t10\t{qlen}\t100\t1\t10\t1\t30\t1e-20\t100\t1\n')
            rows=parse_raw(p,DB,'test',{'gnl|X|odd':'A'*100},{'CBF1':'gnl|X|odd'},RULES)
            self.assertEqual(rows[0]['contig'],'gnl|X|odd')

    def test_circular_extraction_one_revolution(self):
        with tempfile.TemporaryDirectory() as d:
            analyze('test',[hit('cvaC',1,300)],{'c':'A'*1000},{'c'},RULES,Path(d))
            seq=fasta(Path(d)/'candidate_neighborhoods.fna',True)
            self.assertEqual(len(next(iter(seq.values()))),1000)

if __name__=='__main__':unittest.main()

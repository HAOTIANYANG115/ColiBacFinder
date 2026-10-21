import unittest
from colibacfinder.simple_report import build_rows

class SimpleReportTests(unittest.TestCase):
    def summary(self, **kw):
        return dict(sample='x',status='COMPLETE',near_core_loci=kw.get('core',0),
                    near_nonexport_immunity_loci=kw.get('immune',0),primary_candidate_loci=kw.get('primary',0))
    def test_partial_not_main(self):
        row=build_rows([self.summary(primary=1)], [dict(sample='x',system='E2',anchor_locus='L1',status='weak_or_partial_core_candidate')])[0]
        self.assertEqual(row['main_candidates'],'')
        self.assertEqual(row['weak_or_partial_loci'],1)
    def test_cross_hits_one_anchor(self):
        systems=[dict(sample='x',system=s,anchor_locus='L1',anchor_contig='c',status='modeled_components_near_reference') for s in ['E2','E9']]
        row=build_rows([self.summary(core=1)],systems)[0]
        self.assertEqual(row['main_candidates'],'E2/E9')
        self.assertIn('尚未定型',row['notes'])
    def test_failure_not_negative(self):
        row=build_rows([dict(sample='bad',status='FAILED',error='bad fasta')],[])[0]
        self.assertEqual(row['result'],'分析失败，不能判定')
    def test_immunity_only(self):
        row=build_rows([self.summary(immune=1)],[])[0]
        self.assertEqual(row['result'],'仅检出近参考免疫候选')

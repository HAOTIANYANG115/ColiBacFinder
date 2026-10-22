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
        self.assertIn('subtype unresolved',row['notes'])
    def test_failure_not_negative(self):
        row=build_rows([dict(sample='bad',status='FAILED',error='bad fasta')],[])[0]
        self.assertEqual(row['result'],'Analysis failed; no conclusion')
    def test_immunity_only(self):
        row=build_rows([self.summary(immune=1)],[])[0]
        self.assertEqual(row['result'],'Near-reference immunity candidate only')

    def test_cli_help_formats_defaults(self):
        from colibacfinder.cli import arguments
        import contextlib
        import io
        out = io.StringIO()
        with contextlib.redirect_stdout(out), self.assertRaises(SystemExit) as exc:
            arguments().parse_args(['run', '-h'])
        self.assertEqual(exc.exception.code, 0)
        self.assertIn('default: 80', out.getvalue())
        self.assertIn('Nucleotide FASTA', out.getvalue())

    def test_english_html_escapes_sample_names(self):
        from tempfile import TemporaryDirectory
        from pathlib import Path
        from colibacfinder.io import write_tsv
        from colibacfinder.simple_report import render
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_tsv(root/'01_sample_summary.tsv', ['sample','status'],
                      [dict(sample='<sample>',status='COMPLETE')])
            write_tsv(root/'06_system_candidates.tsv',
                      ['sample','system','anchor_locus','status'], [])
            render(root)
            html = (root/'00_results_report.html').read_text()
            self.assertIn('lang="en"', html)
            self.assertIn('&lt;sample&gt;', html)
            self.assertNotIn('<sample>', html)
            self.assertIn('No primary candidate detected', html)

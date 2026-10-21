"""Actual BLAST integration smoke test. Reference-derived, NOT accuracy validation."""
import os
import sys
import csv
import subprocess
import tempfile
from pathlib import Path
from colibacfinder.database import Database
from colibacfinder.io import fasta, read_tsv, write_tsv

db = Database()
with tempfile.TemporaryDirectory(prefix='cbf_smoke_') as tmp:
    root = Path(tmp)
    inputs = root/'inputs'
    inputs.mkdir()
    cds = fasta(db.path/'available_reference_CDS.fna', nucleotide=True)
    links = read_tsv(db.path/'metadata/CDS_links.tsv')
    parts = []
    for label in ['cvaC','cvi','cvaA','cvaB']:
        rid = next(r for r in db.sequences if label in db.values(r,'reporting_label'))
        link = next(r for r in links if r['protein_sequence_id']==rid)
        parts.append(cds[link['cds_id']])
    (inputs/'V_construct.fasta').write_text('>gnl|test|contig\n'+'NNN'*30+('NNN'*30).join(parts)+'NNN'*30+'\n')
    (inputs/'negative.fasta').write_text('>negative\n'+'ACGT'*1000+'\n')
    (inputs/'invalid.fasta').write_text('>bad\nTHIS_IS_NOT_DNA\n')
    def run(*args, expected=0):
        p = subprocess.run([sys.executable,'-m','colibacfinder',*map(str,args)],capture_output=True,text=True)
        print(p.stdout, end='')
        if p.returncode != expected:
            raise AssertionError((p.returncode,p.stderr))
    out = root/'run'
    run('run','-i',inputs,'-o',out,'-t','2',expected=1)
    rows = {r['sample']:r for r in read_tsv(out/'01_sample_summary.tsv')}
    assert rows['invalid']['status']=='FAILED'
    assert rows['negative']['near_core_loci']=='0'
    assert int(rows['V_construct']['near_core_loci'])>=1
    assert any(r['sample']=='V_construct' and r['system']=='Microcin_V' and r['status']=='all_model_components_near_reference' for r in read_tsv(out/'06_system_candidates.tsv'))
    assert all(r['completed_denominator']=='2' for r in read_tsv(out/'08_label_frequency.tsv'))
    run('report','--from-run',out,'-o',root/'report',expected=1)
    assert (out/'03_candidate_loci.tsv').read_bytes()==(root/'report/03_candidate_loci.tsv').read_bytes()
    truth = root/'truth.tsv'
    write_tsv(truth,['sample','target','truth','evidence'],[
        dict(sample=s,target='system:Microcin_V',truth=n,evidence='Synthetic software test only')
        for s,n in [('V_construct',1),('negative',0),('invalid',0)]])
    run('evaluate','--run',out,'--truth',truth,'-o',root/'eval')
    m=read_tsv(root/'eval/metrics.tsv')[0]
    assert (m['TP'],m['TN'],m['FP'],m['FN'])==('1','1','0','0')
    assert any(r['outcome']=='EXCLUDED_NOT_COMPLETE' for r in read_tsv(root/'eval/decisions.tsv'))
    run('run','-i',inputs,'-o',out,expected=2)
    print('PASS: BLAST, V components, negative, failed sample exclusion, re-report, evaluation, overwrite protection')

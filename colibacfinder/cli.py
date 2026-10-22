import argparse
import csv
import json
import shutil
import subprocess
import sys
import platform
from pathlib import Path
from datetime import datetime, timezone
from . import __version__
from .database import Database
from .engine import Rules, parse_raw
from .io import digest, discover, fasta, read_tsv, write_tsv, write_json, write_fasta
from .report import analyze, aggregate, SUMMARY_FIELDS


def arguments():
    p = argparse.ArgumentParser(prog='colibacfinder', description='Reference-based bacteriocin candidate screening; predictions are not activity tests.')
    p.add_argument('--version', action='version', version=__version__)
    sub = p.add_subparsers(dest='command',required=True)
    brief = sub.add_parser('summarize', help='Create simple TSV and HTML from an existing run; no BLAST rerun')
    brief.add_argument('--run', type=Path, required=True)
    check = sub.add_parser('db-check',help='Verify bundled or external reference database')
    check.add_argument('--db',type=Path)
    for cmd in ['run','report']:
        a = sub.add_parser(cmd,formatter_class=argparse.ArgumentDefaultsHelpFormatter,help='Scan FASTA inputs' if cmd=='run' else 'Reanalyze retained raw hits without rerunning BLAST')
        a.add_argument('-i','--input',type=Path,required=cmd=='run',help='Nucleotide FASTA file or directory; .fa/.fna/.fasta, optionally gzipped')
        a.add_argument('-o','--outdir',type=Path,required=True,help='New output directory (must not exist)')
        a.add_argument('--db',type=Path,help='Reference database directory; omit to use the bundled release')
        a.add_argument('--identity',type=float,default=80,help='Minimum identity (percent) for near-reference hits')
        a.add_argument('--coverage',type=float,default=80,help='Single-HSP reference span percent')
        a.add_argument('--evalue',type=float,default=1e-5,help='Screening E-value')
        a.add_argument('--search-evalue',type=float,default=1e-3,help='Maximum E-value retained by the BLAST search')
        a.add_argument('--remote-identity',type=float,default=35,help='Minimum identity (percent) for remote candidates')
        a.add_argument('--remote-coverage',type=float,default=70,help='Minimum single-HSP reference span (percent) for remote candidates')
        a.add_argument('--partial-identity',type=float,default=35,help='Minimum identity (percent) for partial candidates')
        a.add_argument('--partial-coverage',type=float,default=10,help='Minimum single-HSP reference span (percent) for partial candidates')
        a.add_argument('--neighbor-bp',type=int,default=10000,help='Maximum component-to-anchor gap (bp)')
        a.add_argument('--flank-bp',type=int,default=5000,help='Sequence to extract on each side of a candidate (bp)')
        a.add_argument('--circular-map',type=Path,help='TSV columns sample,contig; explicitly closed molecules only')
        if cmd=='run':
            a.add_argument('-t','--threads',type=int,default=8,help='Threads per sample; samples processed sequentially')
            a.add_argument('--keep-db',action='store_true',help='Keep per-sample BLAST databases and normalized FASTA files')
        else:
            a.add_argument('--from-run',type=Path,required=True,help='ColiBacFinder run or original batch run containing parameters.json')
    e = sub.add_parser('evaluate',help='Score explicit independent truth labels against completed predictions')
    e.add_argument('--run',type=Path,required=True)
    e.add_argument('--truth',type=Path,required=True,help='TSV sample,target,truth,evidence; truth is 0 or 1')
    e.add_argument('-o','--outdir',type=Path,required=True)
    return p


def check_tools():
    tools = {}
    for name in ['makeblastdb','tblastn']:
        executable = shutil.which(name)
        if not executable:
            raise ValueError(f'{name} not found. Install NCBI BLAST+ and activate its environment.')
        r = subprocess.run([executable,'-version'],capture_output=True,text=True,check=True)
        tools[name] = dict(path=executable,version=r.stdout.strip())
    return tools


def execute(cmd, log, env=None):
    with log.open('w') as f:
        subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,env=env)


def run(args):
    db = Database(args.db)
    rules = Rules(search_evalue=args.search_evalue,screen_evalue=args.evalue,
        identity=args.identity,coverage=args.coverage,remote_identity=args.remote_identity,
        remote_coverage=args.remote_coverage,partial_identity=args.partial_identity,
        partial_coverage=args.partial_coverage,neighbor_bp=args.neighbor_bp,flank_bp=args.flank_bp)
    rules.validate()
    legacy = False
    oldparams, oldindex = {}, {}
    if args.command=='report':
        origin = args.from_run.resolve()
        oldparams = json.loads((origin/'parameters.json').read_text())
        if args.input is None:
            args.input = Path(oldparams['input'])
        legacy = oldparams.get('software') != 'ColiBacFinder'
        oldindex = {r['sample']:r for r in read_tsv(origin/'input_manifest.tsv')} if not legacy else {}
        if oldparams.get('database_info') and oldparams['database_info']['manifest_sha256'] != db.fingerprint:
            raise ValueError('Reanalysis database fingerprint differs from original run')
        if legacy:
            previous_db = Path(oldparams['database'])
            if previous_db.exists() and digest(previous_db/'SHA256SUMS.txt') != db.fingerprint:
                raise ValueError('Legacy database fingerprint differs')
        original_search = oldparams.get('rules',oldparams).get('search_evalue',1e-3)
        if rules.screen_evalue > original_search:
            raise ValueError('Reanalysis screen E-value exceeds original search E-value')
        rules.search_evalue = original_search
    files = discover(args.input)
    if args.command=='report':
        oldstates = {r['sample']:r['status'] for r in read_tsv(origin/'01_sample_summary.tsv')}
        if set(s for s,_ in files) != set(oldstates):
            raise ValueError('Reanalysis input sample set differs from original run')
    if args.command=='run' and args.threads < 1:
        raise ValueError('threads must be positive')
    tools = check_tools() if args.command=='run' else {}
    circular_map = {}
    if args.command=='report' and args.circular_map is None:
        circular_map = {k:set(v) for k,v in oldparams.get('circular_contigs',{}).items()}
    if args.circular_map:
        for r in read_tsv(args.circular_map):
            if not r.get('sample') or not r.get('contig'):
                raise ValueError('circular-map requires sample and contig columns')
            circular_map.setdefault(r['sample'],set()).add(r['contig'])
        if set(circular_map)-{s for s,_ in files}:
            raise ValueError('circular-map contains samples absent from input')
    out = args.outdir.resolve()
    if out.exists():
        raise ValueError(f'Output already exists; select a new directory: {out}')
    out.mkdir(parents=True)
    (out/'samples').mkdir()
    # Snapshot the exact database, including original provenance and checksum manifest.
    shutil.copytree(db.path,out/'reference_snapshot')
    params = dict(software='ColiBacFinder',version=__version__,
        started_utc=datetime.now(timezone.utc).isoformat(),input=str(args.input.resolve()),
        database_info=db.info(),rules=rules.to_dict(),python=sys.version,platform=platform.platform(),
        tools=tools,threads=getattr(args,'threads',None),
        circular_contigs={k:sorted(v) for k,v in circular_map.items()},
        mode=args.command,source_run=str(args.from_run.resolve()) if args.command=='report' else None,
        legacy_input_identity_verified=False if legacy else True,
        limitations=['single_HSP_reference_span','no_ORF_integrity_validation','no_activity_inference',
                     'subtype_ambiguity_retained','no_HSP_joining_across_contig_origin'])
    write_json(out/'parameters.json',params)
    summaries, results, inputs = [], [], []
    write_json(out/'run_status.json',dict(status='RUNNING',samples=len(files)))
    print(f'ColiBacFinder {__version__}: {len(files)} samples, {len(db.sequences)} proteins',flush=True)
    print(f'Output: {out}',flush=True)
    for i,(sample,path) in enumerate(files,1):
        work = out/'samples'/f'S{i:05d}'
        work.mkdir()
        print(f'[{i}/{len(files)}] {sample}',flush=True)
        try:
            sha = digest(path)
            input_row = dict(sample=sample,input_path=str(path),sha256=sha,sample_directory=str(work.relative_to(out)))
            inputs.append(input_row)
            seqs = fasta(path,nucleotide=True)
            circular = circular_map.get(sample,set())
            if circular-set(seqs):
                raise ValueError('Declared circular contig IDs absent from FASTA')
            if args.command=='run':
                mapping = {f'CBF{n:07d}':c for n,c in enumerate(seqs,1)}
                write_fasta(work/'normalized_input.fna',((sid,seqs[c]) for sid,c in mapping.items()))
                write_tsv(work/'contig_id_map.tsv',['blast_id','original_id'],
                          [dict(blast_id=k,original_id=v) for k,v in mapping.items()])
                prefix = work/'genome_db'
                cmd1 = [tools['makeblastdb']['path'],'-in',str(work/'normalized_input.fna'),
                        '-dbtype','nucl','-out',str(prefix)]
                raw = work/'tblastn.tsv'
                cmd2 = [tools['tblastn']['path'],'-query',str(out/'reference_snapshot/reference_proteins.faa'),
                    '-db',str(prefix),'-db_gencode','11','-evalue',str(rules.search_evalue),'-seg','yes',
                    '-num_threads',str(args.threads),'-max_target_seqs',str(max(5,len(seqs))),
                    '-outfmt','6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore sframe',
                    '-out',str(raw)]
                write_json(work/'commands.json',[cmd1,cmd2])
                import os
                env = os.environ.copy()
                (work/'tmp').mkdir()
                env['TMPDIR'] = str(work/'tmp')
                execute(cmd1,work/'makeblastdb.log',env)
                execute(cmd2,work/'tblastn.log',env)
            else:
                if oldstates.get(sample)!='COMPLETE':
                    raise ValueError('Source sample did not complete; rerun search')
                if legacy:
                    prevwork = origin/'samples'/sample
                    mapping = {c:c for c in seqs}
                else:
                    old = oldindex[sample]
                    if sha != old['sha256']:
                        raise ValueError('Input checksum differs from original run')
                    prevwork = origin/old['sample_directory']
                    mapping = {r['blast_id']:r['original_id'] for r in read_tsv(prevwork/'contig_id_map.tsv')}
                raw = work/'tblastn.tsv'
                shutil.copyfile(prevwork/'tblastn.tsv',raw)
                write_tsv(work/'contig_id_map.tsv',['blast_id','original_id'],
                          [dict(blast_id=k,original_id=v) for k,v in mapping.items()])
            hits = parse_raw(raw,db,sample,seqs,mapping,rules)
            summary, data = analyze(sample,hits,seqs,circular,rules,work)
            if legacy:
                summary['warnings'] += ';legacy_input_checksum_unavailable'
            results.append(data)
            summaries.append(summary)
            if args.command=='run' and not args.keep_db:
                for p in work.glob('genome_db.*'):
                    p.unlink()
                (work/'normalized_input.fna').unlink()
            print(f"  COMPLETE: {summary['near_core_loci']} near-reference core loci",flush=True)
        except Exception as exc:
            summaries.append(dict(sample=sample,status='FAILED',error=str(exc)))
            write_json(work/'error.json',dict(sample=sample,error=str(exc)))
            print(f'  FAILED: {exc}',file=sys.stderr,flush=True)
        write_tsv(out/'01_sample_summary.tsv',SUMMARY_FIELDS,summaries)
        write_tsv(out/'input_manifest.tsv',['sample','input_path','sha256','sample_directory'],inputs)
    aggregate(out,summaries,results)
    failed = sum(s['status']!='COMPLETE' for s in summaries)
    print(f'Finished: {len(summaries)-failed} complete, {failed} failed. {out}',flush=True)
    return 1 if failed else 0


def evaluate(args):
    run_dir = args.run.resolve()
    states = {r['sample']:r['status'] for r in read_tsv(run_dir/'01_sample_summary.tsv')}
    samples = {r['sample']:r for r in read_tsv(run_dir/'01_sample_summary.tsv')}
    systems = read_tsv(run_dir/'06_system_candidates.tsv')
    from .systems import MODELS
    truth = read_tsv(args.truth)
    if not truth:
        raise ValueError('Truth table is empty')
    seen, rows, counts = set(), [], {}
    for r in truth:
        sample,target,label,evidence = (r.get(k,'').strip() for k in ['sample','target','truth','evidence'])
        if not sample or not target or label not in {'0','1'} or not evidence:
            raise ValueError('Truth needs sample,target,truth(0/1),evidence; negative truth must be supported')
        if (sample,target) in seen:
            raise ValueError('Duplicate sample/target truth row')
        seen.add((sample,target))
        valid_target = target in {'ANY_CORE','ANY_NONEXPORT_IMMUNITY'} or (target.startswith('system:') and target[7:] in MODELS and MODELS[target[7:]]['complete'])
        if not valid_target:
            raise ValueError(f'Unsupported target: {target}; use ANY_CORE, ANY_NONEXPORT_IMMUNITY, system:Microcin_V/J25/B17')
        if states.get(sample)!='COMPLETE':
            rows.append(dict(sample=sample,target=target,truth=label,prediction='',outcome='EXCLUDED_NOT_COMPLETE',evidence=evidence))
            continue
        if target=='ANY_CORE':
            pred = int(samples[sample]['near_core_loci'])>0
        elif target=='ANY_NONEXPORT_IMMUNITY':
            pred = int(samples[sample]['near_nonexport_immunity_loci'])>0
        else:
            pred = any(s['sample']==sample and s['system']==target[7:] and s['status']=='all_model_components_near_reference' for s in systems)
        outcome = ('T' if pred == (label=='1') else 'F') + ('P' if pred else 'N')
        c = counts.setdefault(target,dict(TP=0,FP=0,TN=0,FN=0))
        c[outcome] += 1
        rows.append(dict(sample=sample,target=target,truth=label,prediction=int(pred),outcome=outcome,evidence=evidence))
    out = args.outdir.resolve()
    if out.exists():
        raise ValueError('Evaluation output directory already exists')
    out.mkdir(parents=True)
    metrics=[]
    for target,c in counts.items():
        def ratio(a,b):
            return a/b if b else ''
        tp,fp,tn,fn=(c[k] for k in ['TP','FP','TN','FN'])
        metrics.append(dict(target=target,**c,sensitivity=ratio(tp,tp+fn),specificity=ratio(tn,tn+fp),
            precision=ratio(tp,tp+fp),f1=ratio(2*tp,2*tp+fp+fn)))
    write_tsv(out/'decisions.tsv',['sample','target','truth','prediction','outcome','evidence'],rows)
    write_tsv(out/'metrics.tsv',['target','TP','FP','TN','FN','sensitivity','specificity','precision','f1'],metrics)
    write_json(out/'evaluation.json',dict(run=str(run_dir),truth_sha256=digest(args.truth),
        software_version=__version__,truth_rows=len(truth),excluded=sum(r['outcome'].startswith('EXCLUDED') for r in rows),
        caution='Metrics are only as valid as supplied independent ground truth; independence is not automatically established.'))
    print(f'Evaluation written: {out}')
    return 0


def main(argv=None):
    p = arguments()
    args = p.parse_args(argv)
    try:
        if args.command=='summarize':
            from .simple_report import render
            render(args.run)
            print(f'Simple reports: {args.run}/00_results_summary.tsv and 00_results_report.html')
            return 0
        if args.command=='db-check':
            print(json.dumps(Database(args.db).info(),ensure_ascii=False,indent=2))
            return 0
        if args.command=='evaluate':
            return evaluate(args)
        return run(args)
    except (ValueError,OSError,KeyError,subprocess.SubprocessError) as exc:
        print(f'ERROR: {exc}',file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print('Interrupted; partial output retained, run_status is not COMPLETE.',file=sys.stderr)
        return 130

"""Strict FASTA handling and reproducible tabular output; standard library only."""
import csv
import gzip
import hashlib
import json
from pathlib import Path

SUFFIXES = ('.fasta', '.fna', '.fa')

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()

def fasta(path, nucleotide=False):
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    records, name, parts = {}, None, []
    def add():
        if name is None:
            return
        seq = ''.join(parts).upper()
        if not seq:
            raise ValueError(f'Empty sequence: {name}')
        alphabet = set('ACGTRYSWKMBDHVN') if nucleotide else set('ACDEFGHIKLMNPQRSTVWYBXZJUO')
        bad = set(seq) - alphabet
        if bad:
            raise ValueError(f'Invalid sequence characters in {name}: {sorted(bad)}')
        records[name] = seq
    with opener(path, 'rt', encoding='utf-8-sig') as f:
        for number, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                add()
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f'Empty FASTA header at line {number}')
                name = fields[0]
                if name in records:
                    raise ValueError(f'Duplicate FASTA ID: {name}')
                parts = []
            elif name is None:
                raise ValueError(f'Non-FASTA content at line {number}')
            else:
                parts.append(''.join(line.split()))
        add()
    if not records:
        raise ValueError('No FASTA records')
    return records

def write_fasta(path, records):
    with Path(path).open('w', encoding='utf-8') as f:
        for name, seq in records:
            f.write('>' + name + '\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

def read_tsv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(path, fields, rows):
    with Path(path).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

def tokens(s):
    return {x.strip() for x in s.replace('|', ';').split(';') if x.strip()}

def discover(path):
    path = Path(path).resolve()
    paths = [path] if path.is_file() else sorted(p for p in path.iterdir() if p.is_file())
    result = []
    for p in paths:
        name = p.name[:-3] if p.name.lower().endswith('.gz') else p.name
        if not name.lower().endswith(SUFFIXES):
            if path.is_file():
                raise ValueError('Input must end in .fa/.fna/.fasta, optionally .gz')
            continue
        sample = name.rsplit('.', 1)[0]
        if not sample or any(c in sample for c in '\t\r\n'):
            raise ValueError(f'Invalid sample filename: {p.name}')
        result.append((sample, p))
    if not result:
        raise ValueError('No top-level FASTA input files found')
    if len({s for s, _ in result}) != len(result):
        raise ValueError('Duplicate sample names after removing FASTA extensions')
    return result

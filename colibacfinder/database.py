from pathlib import Path
from collections import defaultdict
import hashlib
from .io import digest, fasta, read_tsv, tokens

RELEASE = 'bacteriocin_unified_v0.1_20261020_194442_905507'
DEFAULT_DB = Path(__file__).parent / 'data' / RELEASE
PINNED_MANIFEST = '55eceee074f5a32343c14e469b2c5d333737c88627eaef3fcfaf8fb89d560cd2'

class Database:
    def __init__(self, path=None):
        self.path = Path(path or DEFAULT_DB).resolve()
        manifest = self.path / 'SHA256SUMS.txt'
        self.fingerprint = digest(manifest)
        if self.path == DEFAULT_DB.resolve() and self.fingerprint != PINNED_MANIFEST:
            raise ValueError('Bundled database differs from pinned release')
        listed = set()
        for line in manifest.read_text().splitlines():
            if not line.strip():
                continue
            expected, name = line.split(None, 1)
            name = name.lstrip('*')
            target = (self.path / name).resolve()
            if self.path not in target.parents:
                raise ValueError('Unsafe checksum manifest path')
            if digest(target) != expected:
                raise ValueError(f'Database checksum mismatch: {name}')
            listed.add(name)
        for required in ['reference_proteins.faa', 'metadata/reference_links.tsv']:
            if required not in listed:
                raise ValueError(f'Database manifest does not cover {required}')
        self.ref = self.path / 'reference_proteins.faa'
        self.sequences = fasta(self.ref)
        self.meta = defaultdict(list)
        for row in read_tsv(self.path / 'metadata/reference_links.tsv'):
            self.meta[row['sequence_id']].append(row)
        if set(self.meta) != set(self.sequences):
            raise ValueError('Protein IDs and metadata IDs differ')
        for rid, seq in self.sequences.items():
            if rid != 'BAC_' + hashlib.sha256(seq.encode()).hexdigest():
                raise ValueError(f'Protein sequence hash ID mismatch: {rid}')

    def values(self, rid, key):
        return sorted({r.get(key, '').strip() for r in self.meta[rid] if r.get(key, '').strip()})

    def roles(self, rid):
        return set().union(*(tokens(r.get('roles', '')) for r in self.meta[rid]))

    def info(self):
        return dict(release=self.path.name, path=str(self.path), proteins=len(self.sequences),
                    annotation_links=sum(map(len, self.meta.values())),
                    manifest_sha256=self.fingerprint,
                    protein_fasta_sha256=digest(self.ref), checksum_status='PASS')

"""Validated, explicit input schemas. Public coordinates are 1-based inclusive."""
import csv
import gzip
import hashlib
import re
from datetime import date
from pathlib import Path
from Bio import SeqIO

DNA = set('ACGTRYSWKMBDHVN')

def read_tsv(path, required=()):
    with open(path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames is None or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f'{path}: missing or duplicate headers')
        missing = set(required) - set(reader.fieldnames)
        if missing:
            raise ValueError(f'{path}: missing columns {sorted(missing)}')
        rows = list(reader)
        for row in rows:
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f'{path}: ragged TSV row')
        return rows

def write_tsv(path, rows, fields):
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def fasta(path, protein=False):
    opener = gzip.open if str(path).endswith('.gz') else open
    records = {}
    with opener(path, 'rt') as handle:
        for rec in SeqIO.parse(handle, 'fasta'):
            seq = str(rec.seq).upper()
            if rec.id in records:
                raise ValueError(f'{path}: duplicate FASTA ID {rec.id}')
            if not seq or (not protein and set(seq) - DNA):
                raise ValueError(f'{path}: empty or invalid sequence {rec.id}')
            if protein and set(seq) - set('ACDEFGHIKLMNPQRSTVWYBXZJUO*'):
                raise ValueError(f'{path}: invalid protein {rec.id}')
            records[rec.id] = seq
    if not records:
        raise ValueError(f'{path}: empty FASTA')
    return records

def write_fasta(path, records):
    with open(path, 'w') as handle:
        for name, seq in records.items():
            handle.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i+80] + '\n')

def load_samples(path):
    path = Path(path).resolve()
    rows = read_tsv(path, ['sample_id', 'assembly'])
    seen = set()
    for r in rows:
        sid = r['sample_id']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sid) or sid in seen:
            raise ValueError(f'Invalid or duplicate sample_id: {sid}')
        seen.add(sid)
        for k in ('subject_id', 'strain_id', 'date', 'lineage_evidence', 'genbank'):
            r.setdefault(k, '')
        if r['date']:
            date.fromisoformat(r['date'])
        for key in ('assembly', 'genbank'):
            if r[key]:
                r[key] = str((path.parent / r[key]).resolve())
                if not Path(r[key]).is_file():
                    raise ValueError(f'Missing {key}: {r[key]}')
    if not rows:
        raise ValueError('Sample sheet is empty')
    return rows

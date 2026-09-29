"""Caller output adapters. Unknown formats fail rather than invent coordinates."""
import csv
from pathlib import Path
from .io import read_tsv
from .regions import CALL_FIELDS


def import_calls(path,caller,sample_id,coordinates='1-inclusive',columns=None):
    rows=[]
    if caller=='phispy':
        for line in Path(path).read_text().splitlines():
            if not line.strip() or line.startswith('#'):continue
            f=line.split('\t')
            if len(f)<4:raise ValueError('PhiSpy coordinates require at least four columns')
            if f[2].lower() in ('start','start_position'):continue
            rows.append(dict(sample_id=sample_id,contig=f[1],start=int(f[2]),end=int(f[3]),caller=caller,call_id=f[0]))
    elif caller=='phageboost':
        for line in Path(path).read_text().splitlines():
            if not line.strip() or line.startswith('#'):continue
            f=line.split('\t')
            if len(f)!=9:raise ValueError('PhageBoost: use phages_*.gff, not genecalls_*.gff3')
            attrs=dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
            rows.append(dict(sample_id=sample_id,contig=f[0],start=int(f[3]),end=int(f[4]),caller=caller,
                call_id=attrs.get('phage_id',attrs.get('ID',f'call{len(rows)+1}'))))
    elif caller in ('vibrant','generic'):
        # Caller versions can vary: explicit column names are supported.
        fields=columns or ('scaffold','start','stop')
        data=read_tsv(path,fields)
        for i,r in enumerate(data,1):
            rows.append(dict(sample_id=sample_id,contig=r[fields[0]],start=int(r[fields[1]]),end=int(r[fields[2]]),caller=caller,call_id=f'call{i}'))
    else:raise ValueError(f'Unknown caller {caller}')
    if coordinates=='0-half-open':
        for r in rows:r['start']+=1
    elif coordinates!='1-inclusive':raise ValueError('Unknown coordinate convention')
    for r in rows:
        if r['start']<1 or r['end']<r['start']:raise ValueError('Invalid coordinates after conversion')
    return rows

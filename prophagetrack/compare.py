"""Exact groups and optional alignment evidence; no transitive identity claims."""
import hashlib
import subprocess
from pathlib import Path
from .io import write_fasta

_COMP = str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')

def reverse_complement(seq):
    return seq.translate(_COMP)[::-1]

def canonical(seq):
    return min(seq, reverse_complement(seq))

def exact_groups(regions):
    for r in regions:
        # An N-to-N match is not a biological exact identity claim.
        if set(r['sequence'])-set('ACGT'):
            r['exact_group'] = 'UNCERTAIN_' + r['region_id']
        else:
            r['exact_group'] = 'PG_' + hashlib.sha256(canonical(r['sequence']).encode()).hexdigest()[:16]
    return regions

def mmseqs_search(records, out, protein=False, threads=1, executable='mmseqs'):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if len(records)<2: return []
    query=out/'sequences.fasta'; hits=out/'alignments.tsv'
    write_fasta(query,records)
    # Single best local alignment per pair: deliberately conservative for mosaics.
    cmd=[executable,'easy-search',str(query),str(query),str(hits),str(out/'tmp'),
         '--search-type','1' if protein else '3','--threads',str(threads),
         '--format-output','query,target,fident,qcov,tcov,alnlen,qstart,qend,tstart,tend',
         '--max-seqs',str(max(300,len(records))),'-e','1e-5']
    if not protein: cmd+=['--strand','2']
    with open(out/'command.txt','w') as h: h.write(repr(cmd)+'\n')
    with open(out/'run.log','w') as log:
        subprocess.run(cmd,check=True,stdout=log,stderr=subprocess.STDOUT)
    result=[]
    for line in hits.read_text().splitlines():
        f=line.split('\t')
        if len(f)!=10: raise ValueError('Unexpected MMseqs2 output schema')
        if f[0]==f[1]: continue
        result.append(dict(query=f[0],target=f[1],identity=float(f[2]),query_coverage=float(f[3]),
            target_coverage=float(f[4]),alignment_length=int(f[5]),qstart=int(f[6]),qend=int(f[7]),tstart=int(f[8]),tend=int(f[9])))
    return result

def related_pairs(regions, hits, identity=0.95, coverage=0.85):
    byid={r['region_id']:r for r in regions}; best={}
    for h in hits:
        a,b=sorted((h['query'],h['target']))
        if a==b or byid[a]['exact_group']==byid[b]['exact_group']:continue
        row=dict(h)
        if h['query']!=a:
            row['query'],row['target']=a,b
            row['query_coverage'],row['target_coverage']=h['target_coverage'],h['query_coverage']
            row['qstart'],row['qend'],row['tstart'],row['tend']=h['tstart'],h['tend'],h['qstart'],h['qend']
        if row['identity']<identity:continue
        row['relationship']='related' if min(row['query_coverage'],row['target_coverage'])>=coverage else 'fragment_match'
        if row['relationship']=='fragment_match' and max(row['query_coverage'],row['target_coverage'])<coverage:continue
        key=(a,b)
        if key not in best or row['alignment_length']>best[key]['alignment_length']:best[key]=row
    return [best[k] for k in sorted(best)]

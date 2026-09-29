"""Interval-level accuracy against independent coordinate truth."""
from collections import defaultdict
from .io import read_tsv


def union(intervals):
    result=[]
    for a,b in sorted(intervals):
        if result and a<=result[-1][1]+1:result[-1]=(result[-1][0],max(b,result[-1][1]))
        else:result.append((a,b))
    return result

def evaluate(truth_path,prediction_path):
    def load(p):
        d=defaultdict(list)
        for r in read_tsv(p,['sample_id','contig','start','end']):
            a,b=int(r['start']),int(r['end'])
            if not 1<=a<=b:raise ValueError('Invalid benchmark coordinates')
            d[r['sample_id'],r['contig']].append((a,b))
        return {k:union(v) for k,v in d.items()}
    t,p=load(truth_path),load(prediction_path)
    truth=sum(b-a+1 for v in t.values() for a,b in v)
    pred=sum(b-a+1 for v in p.values() for a,b in v)
    tp=sum(max(0,min(b,d)-max(a,c)+1) for k in t.keys()|p.keys() for a,b in t.get(k,[]) for c,d in p.get(k,[]))
    return dict(true_positive_bp=tp,false_positive_bp=pred-tp,false_negative_bp=truth-tp,
        nucleotide_precision=tp/pred if pred else None,nucleotide_recall=tp/truth if truth else None,
        note='Truth must be independently curated and exhaustive over evaluated assemblies. This measures nucleotide overlap, not biological event accuracy.')

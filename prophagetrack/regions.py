"""Preserve raw calls; reconcile only unambiguous reciprocal-overlap cliques."""
from collections import defaultdict
from .io import read_tsv

CALL_FIELDS = ['sample_id','contig','start','end','caller','call_id']

def validate_calls(rows, assemblies):
    seen = set()
    for i, r in enumerate(rows):
        if r['sample_id'] not in assemblies or r['contig'] not in assemblies[r['sample_id']]:
            raise ValueError(f'Call references unknown sample/contig: {r}')
        r['start'], r['end'] = int(r['start']), int(r['end'])
        if not 1 <= r['start'] <= r['end'] <= len(assemblies[r['sample_id']][r['contig']]):
            raise ValueError(f'Out-of-range call coordinates: {r}')
        if not r['caller'].strip() or not r['call_id'].strip():
            raise ValueError('caller and call_id cannot be blank')
        key = (r['sample_id'],r['caller'],r['call_id'])
        if key in seen:
            raise ValueError(f'Duplicate caller ID: {key}')
        seen.add(key)
    return rows

def reconcile(rows, assemblies, overlap=0.8):
    grouped = defaultdict(list)
    for r in rows:
        grouped[r['sample_id'], r['contig']].append(r)
    regions = []
    for (sid, ctg), calls in sorted(grouped.items()):
        calls = sorted(calls, key=lambda x:(x['start'],x['end'],x['caller'],x['call_id']))
        adjacency = {i:set() for i in range(len(calls))}
        overlaps = set()
        for i,a in enumerate(calls):
            for j in range(i+1,len(calls)):
                b = calls[j]
                n = max(0, min(a['end'],b['end'])-max(a['start'],b['start'])+1)
                if n:
                    overlaps.update([i,j])
                if n / (a['end']-a['start']+1) >= overlap and n / (b['end']-b['start']+1) >= overlap:
                    adjacency[i].add(j); adjacency[j].add(i)
        pending = set(adjacency)
        while pending:
            seed = min(pending); component = {seed}; stack = [seed]
            while stack:
                for nxt in adjacency[stack.pop()] - component:
                    component.add(nxt); stack.append(nxt)
            pending -= component
            clique = all(component-{i} <= adjacency[i] for i in component)
            unique_callers = len({calls[i]['caller'] for i in component}) == len(component)
            # Do not merge chains or multiple predictions from the same caller.
            partitions = [component] if clique and unique_callers else [{i} for i in sorted(component)]
            for part in partitions:
                members = [calls[i] for i in sorted(part)]
                # Use an actual predicted interval, never an invented union/intersection.
                chosen = min(members,key=lambda a:(sum(abs(a['start']-b['start'])+abs(a['end']-b['end']) for b in members),a['caller'],a['call_id']))
                start,end = chosen['start'],chosen['end']
                bounds = {(x['start'],x['end']) for x in members}
                flags=[]
                if len(bounds)>1: flags.append('boundary_disagreement')
                if any(i in overlaps and (adjacency[i] - part or any(j not in part and max(0,min(calls[i]['end'],calls[j]['end'])-max(calls[i]['start'],calls[j]['start'])+1)>0 for j in range(len(calls)))) for i in part):
                    flags.append('overlap_ambiguity')
                if start==1 or end==len(assemblies[sid][ctg]): flags.append('contig_edge')
                seq=assemblies[sid][ctg][start-1:end]
                if set(seq)-set('ACGT'): flags.append('ambiguous_bases')
                regions.append(dict(sample_id=sid,contig=ctg,start=start,end=end,length=len(seq),
                    callers=';'.join(sorted({x['caller'] for x in members})),n_callers=len({x['caller'] for x in members}),
                    source_calls=';'.join(x['caller']+':'+x['call_id'] for x in members),
                    selected_call=chosen['caller']+':'+chosen['call_id'],flags=';'.join(flags),sequence=seq))
    regions.sort(key=lambda r:(r['sample_id'],r['contig'],r['start'],r['end'],r['selected_call']))
    for i,r in enumerate(regions,1): r['region_id']=f'PT{i:06d}'
    return regions

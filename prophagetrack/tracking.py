"""Conservative exact-flank projection. No strain inference or causal claims."""
from collections import defaultdict, Counter
from .compare import reverse_complement, canonical


def find_all(text, pattern):
    start=0
    while True:
        pos=text.find(pattern,start)
        if pos<0:return
        yield pos
        start=pos+1

def project(left, right, sequence, assembly, max_empty_gap=50):
    """Require each flanking anchor to have exactly one hit over both strands."""
    left_hits=[];right_hits=[]
    for ctg,seq in assembly.items():
        for strand,oriented in [('+',seq),('-',reverse_complement(seq))]:
            left_hits.extend((ctg,strand,p) for p in find_all(oriented,left))
            right_hits.extend((ctg,strand,p) for p in find_all(oriented,right))
    if len(left_hits)!=1 or len(right_hits)!=1:
        return dict(state='unresolved',reason='missing_or_nonunique_flanks',contig='',strand='',start='',end='',gap_length='')
    ctg,strand,l=left_hits[0];rctg,rstrand,r=right_hits[0]
    if (ctg,strand)!=(rctg,rstrand) or r<l+len(left):
        return dict(state='unresolved',reason='flanks_not_collinear',contig='',strand='',start='',end='',gap_length='')
    oriented=assembly[ctg] if strand=='+' else reverse_complement(assembly[ctg])
    begin=l+len(left); gap=oriented[begin:r]
    a,b=(begin+1,r) if strand=='+' else (len(oriented)-r+1,len(oriented)-begin)
    # Empty sites use flanking positions instead of a misleading nonempty interval.
    base=dict(contig=ctg,strand=strand,start=a if gap else '',end=b if gap else '',gap_length=len(gap))
    if set(gap)-set('ACGT'):
        return dict(base,state='unresolved',reason='ambiguous_bases_between_flanks')
    if gap==sequence:
        return dict(base,state='occupied',reason='exact_sequence_at_corresponding_site')
    if len(gap)<=max_empty_gap:
        return dict(base,state='candidate_empty',reason='short_resolved_interval_between_unique_flanks')
    return dict(base,state='different_occupancy',reason='different_sequence_between_corresponding_flanks')

def track(regions,samples,assemblies,flank_length=200,max_empty_gap=50):
    cohorts=defaultdict(list)
    for s in samples:
        if all(s.get(k) for k in ('subject_id','strain_id','date','lineage_evidence')):
            cohorts[s['subject_id'],s['strain_id']].append(s)
    bysample={s['sample_id']:s for s in samples}
    observations=[];events=[];seen=set()
    for r in regions:
        sample=bysample[r['sample_id']];key=(sample['subject_id'],sample['strain_id'])
        cohort=cohorts.get(key,[])
        if len(cohort)<2 or not sample['lineage_evidence']:continue
        seq=assemblies[r['sample_id']][r['contig']]; start=r['start']-1;end=r['end']
        if r['flags'] or r['length']<=max_empty_gap or start<flank_length or len(seq)-end<flank_length:continue
        left=seq[start-flank_length:start];right=seq[end:end+flank_length]
        if set(left+right)-set('ACGT'):continue
        locus=(key,min((left,right),(reverse_complement(right),reverse_complement(left))),r['exact_group'])
        if locus in seen:continue
        seen.add(locus)
        previous=None
        date_counts=Counter(s['date'] for s in cohort)
        for s in sorted(cohort,key=lambda x:(x['date'],x['sample_id'])):
            state=project(left,right,r['sequence'],assemblies[s['sample_id']],max_empty_gap)
            obs=dict(reference_region=r['region_id'],exact_group=r['exact_group'],subject_id=key[0],strain_id=key[1],
                sample_id=s['sample_id'],date=s['date'],lineage_evidence=s['lineage_evidence'],**state)
            observations.append(obs)
            if previous:
                before,after=previous['state'],obs['state']
                pair=(before,after)
                label={('occupied','occupied'):'persistence',('candidate_empty','occupied'):'candidate_gain',
                    ('occupied','candidate_empty'):'candidate_loss',('occupied','different_occupancy'):'candidate_modification_or_replacement',
                    ('different_occupancy','occupied'):'candidate_modification_or_replacement',
                    ('candidate_empty','candidate_empty'):'no_tracked_prophage_at_site'}.get(pair,'unresolved')
                if previous['date']==obs['date']:label='unresolved_same_date'
                elif date_counts[previous['date']]>1 or date_counts[obs['date']]>1:label='unresolved_replicates'
                events.append(dict(reference_region=r['region_id'],exact_group=r['exact_group'],subject_id=key[0],strain_id=key[1],
                    from_sample=previous['sample_id'],to_sample=obs['sample_id'],from_date=previous['date'],to_date=obs['date'],
                    from_state=before,to_state=after,event=label,
                    evidence='unique_exact_flanks;supplied_lineage_assignment;assembly_only'))
            previous=obs
    return observations,events

def recover_exact(regions,assemblies):
    """Search all assemblies, not just caller-positive samples; not proof of integration."""
    representatives={}
    for r in regions:
        if not r['exact_group'].startswith('UNCERTAIN_'):representatives.setdefault(r['exact_group'],r['sequence'])
    rows=[]
    for group,seq in sorted(representatives.items()):
        for sid,assembly in sorted(assemblies.items()):
            for ctg,target in assembly.items():
                for strand,query in [('+',seq),('-',reverse_complement(seq))]:
                    if strand=='-' and query==seq:continue
                    for p in find_all(target,query):
                        rows.append(dict(exact_group=group,sample_id=sid,contig=ctg,start=p+1,end=p+len(seq),strand=strand,
                            evidence='exact_sequence_match;integration_not_established_by_match_alone'))
    return rows

"""Consistent gene prediction or explicit supplied-protein mapping."""
from collections import defaultdict
from .io import fasta,read_tsv,write_fasta,write_tsv
from .compare import mmseqs_search


def get_genes(regions,out,proteins=None,gene_map=None):
    if proteins:
        seqs=fasta(proteins,protein=True)
        rows=read_tsv(gene_map,['protein_id','region_id','start','end','strand'])
        valid={r['region_id']:r for r in regions};seen=set()
        for row in rows:
            pid=row['protein_id'];rid=row['region_id']
            if pid in seen or pid not in seqs or rid not in valid:raise ValueError('Invalid/duplicate supplied protein mapping')
            seen.add(pid)
            row['start'],row['end']=int(row['start']),int(row['end'])
            if not 1<=row['start']<=row['end']<=valid[rid]['length'] or row['strand'] not in ('+','-'):
                raise ValueError('Invalid gene coordinates (relative to extracted region)')
        if seen!=set(seqs):raise ValueError('Every supplied protein needs exactly one mapping')
        for r in rows:r['gene_source']='supplied';r['partial']='unknown'
    else:
        try:import pyrodigal
        except ImportError as exc:raise ValueError('Gene calling needs pyrodigal: install prophagetrack[genes], supply --proteins and --gene-map, or explicitly use --skip-genes') from exc
        finder=pyrodigal.GeneFinder(meta=True)
        seqs={};rows=[]
        for r in regions:
            for i,g in enumerate(finder.find_genes(r['sequence'].encode()),1):
                pid=f"{r['region_id']}_g{i:05d}"
                protein=g.translate().rstrip('*')
                if not protein:continue
                seqs[pid]=protein
                rows.append(dict(protein_id=pid,region_id=r['region_id'],start=g.begin,end=g.end,
                    strand='+' if g.strand==1 else '-',partial=bool(g.partial_begin or g.partial_end),gene_source='pyrodigal_meta'))
    write_fasta(out/'proteins.faa',seqs)
    return seqs,rows

def families(seqs,rows,out,mode='mmseqs',identity=0.5,coverage=0.8,threads=1):
    # Greedy representative groups: every member passes thresholds to representative.
    # Exact proteins always grouped, even when search misses a short sequence.
    neighbors=defaultdict(set)
    if mode=='mmseqs' and len(seqs)>1:
        for h in mmseqs_search(seqs,out/'protein_search',protein=True,threads=threads):
            if h['identity']>=identity and min(h['query_coverage'],h['target_coverage'])>=coverage:
                neighbors[h['query']].add(h['target']);neighbors[h['target']].add(h['query'])
    exact=defaultdict(set)
    for pid,seq in seqs.items():exact[seq].add(pid)
    assigned={};members=[]
    for pid in sorted(seqs,key=lambda p:(-len(seqs[p]),p)):
        if pid in assigned:continue
        family=f'GF{len(members)+1:06d}'
        group=({pid}|neighbors[pid]|exact[seqs[pid]])-set(assigned)
        for g in sorted(group):assigned[g]=family
        members.append((family,pid,sorted(group)))
    for row in rows:row['gene_family']=assigned[row['protein_id']]
    mapping={r['protein_id']:r['region_id'] for r in rows}
    sharing=[]
    for family,representative,group in members:
        prophages=sorted({mapping[p] for p in group})
        sharing.append(dict(gene_family=family,representative=representative,n_proteins=len(group),
            n_regions=len(prophages),regions=';'.join(prophages),method=mode))
    return rows,sharing

def shared_blocks(rows,min_genes=3):
    # Exact family windows in either orientation; gene order, not physical linkage proof.
    byregion=defaultdict(list)
    for r in rows:byregion[r['region_id']].append(r)
    index=defaultdict(set)
    for rid,genes in byregion.items():
        genes.sort(key=lambda r:r['start'])
        for i in range(len(genes)-min_genes+1):
            ids=tuple(g['gene_family'] for g in genes[i:i+min_genes])
            index[min(ids,ids[::-1])].add(rid)
    return [dict(families=';'.join(ids),n_regions=len(rids),regions=';'.join(sorted(rids)))
        for ids,rids in sorted(index.items()) if len(rids)>1]

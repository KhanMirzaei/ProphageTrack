import json
import platform
import shutil
import subprocess
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from . import __version__
from .io import load_samples,fasta,read_tsv,write_tsv,write_fasta,digest
from .regions import CALL_FIELDS,validate_calls,reconcile
from .compare import exact_groups,mmseqs_search,related_pairs
from .tracking import track,recover_exact
from .genes import get_genes,families,shared_blocks
from .report import render

REGION_FIELDS=['region_id','sample_id','contig','start','end','length','callers','n_callers','source_calls','selected_call','flags','exact_group']
OBS_FIELDS=['reference_region','exact_group','subject_id','strain_id','sample_id','date','state','reason','contig','strand','start','end','gap_length','lineage_evidence']
EVENT_FIELDS=['reference_region','exact_group','subject_id','strain_id','from_sample','to_sample','from_date','to_date','from_state','to_state','event','evidence']
PAIR_FIELDS=['query','target','identity','query_coverage','target_coverage','alignment_length','qstart','qend','tstart','tend','relationship']
GENE_FIELDS=['protein_id','region_id','start','end','strand','partial','gene_source','gene_family']
FAMILY_FIELDS=['gene_family','representative','n_proteins','n_regions','regions','method']

def run(args):
    out=Path(args.out).resolve()
    if out.exists() and any(out.iterdir()):raise ValueError(f'Output directory must be empty: {out}')
    if not 0<args.overlap<=1 or not 0<args.identity<=1 or not 0<args.coverage<=1 or not 0<args.protein_identity<=1 or not 0<args.protein_coverage<=1:
        raise ValueError('Identity, coverage and overlap must be in (0, 1]')
    if args.flank_length<20 or args.max_empty_gap<0 or args.threads<1:raise ValueError('Flanks must be >=20 bp; empty gap >=0; threads >=1')
    if bool(args.proteins)!=bool(args.gene_map):raise ValueError('--proteins and --gene-map must be supplied together')
    if args.skip_genes and args.proteins:raise ValueError('Cannot combine --skip-genes and --proteins')
    if (args.related=='mmseqs' or (not args.skip_genes and args.gene_mode=='mmseqs')) and not shutil.which('mmseqs'):
        raise ValueError('MMseqs2 not found. Install it, or explicitly choose --related none --gene-mode exact')
    if not args.skip_genes and not args.proteins:
        try:import pyrodigal
        except ImportError:raise ValueError('Install prophagetrack[genes], supply proteins/mapping, or use --skip-genes')
    samples=load_samples(args.samples)
    assemblies={s['sample_id']:fasta(s['assembly']) for s in samples}
    calls=[]
    for path in args.calls:calls+=read_tsv(path,CALL_FIELDS)
    calls=validate_calls(calls,assemblies)
    regions=exact_groups(reconcile(calls,assemblies,args.overlap))
    out.mkdir(parents=True,exist_ok=True)
    write_tsv(out/'raw_calls.tsv',calls,CALL_FIELDS)
    write_tsv(out/'regions.tsv',regions,REGION_FIELDS)
    write_fasta(out/'prophages.fna',{r['region_id']:r['sequence'] for r in regions})
    with open(out/'prophages.gff3','w') as h:
        from urllib.parse import quote
        h.write('##gff-version 3\n')
        for r in regions:h.write(f"{quote(r['sample_id']+'|'+r['contig'],safe='')}\tProphageTrack\tprophage_region\t{r['start']}\t{r['end']}\t.\t.\t.\tID={r['region_id']};Group={r['exact_group']}\n")
    hits=mmseqs_search({r['region_id']:r['sequence'] for r in regions},out/'nucleotide_search',threads=args.threads) if args.related=='mmseqs' else []
    pairs=related_pairs(regions,hits,args.identity,args.coverage)
    write_tsv(out/'related_pairs.tsv',pairs,PAIR_FIELDS)
    recovered=recover_exact(regions,assemblies)
    write_tsv(out/'exact_recovery.tsv',recovered,['exact_group','sample_id','contig','start','end','strand','evidence'])
    groups=sorted({r['exact_group'] for r in regions})
    called={(r['sample_id'],r['exact_group']) for r in regions}
    found={(r['sample_id'],r['exact_group']) for r in recovered}
    write_tsv(out/'genome_prophage_matrix.tsv',[dict(sample_id=s['sample_id'],**{g:('called_ambiguous' if g.startswith('UNCERTAIN_') and (s['sample_id'],g) in called else 'called_exact' if (s['sample_id'],g) in called else 'recovered_exact' if (s['sample_id'],g) in found else 'not_detected') for g in groups}) for s in samples],['sample_id']+groups)
    observations,events=track(regions,samples,assemblies,args.flank_length,args.max_empty_gap)
    write_tsv(out/'site_observations.tsv',observations,OBS_FIELDS)
    write_tsv(out/'events.tsv',events,EVENT_FIELDS)
    gene_rows=[];family_rows=[]
    if not args.skip_genes and regions:
        proteins,gene_rows=get_genes(regions,out,args.proteins,args.gene_map)
        gene_rows,family_rows=families(proteins,gene_rows,out,args.gene_mode,args.protein_identity,args.protein_coverage,args.threads)
    write_tsv(out/'genes.tsv',gene_rows,GENE_FIELDS)
    write_tsv(out/'gene_families.tsv',family_rows,FAMILY_FIELDS)
    write_tsv(out/'shared_gene_blocks.tsv',shared_blocks(gene_rows),['families','n_regions','regions'])
    fams=sorted({r['gene_family'] for r in gene_rows});counts=Counter((g['region_id'],g['gene_family']) for g in gene_rows)
    write_tsv(out/'prophage_gene_matrix.tsv',[dict(region_id=r['region_id'],**{f:counts[r['region_id'],f] for f in fams}) for r in regions],['region_id']+fams)
    edges=[dict(source='sample:'+r['sample_id'],target=r['region_id'],relationship='predicted_region') for r in regions]
    edges += [dict(source=r['region_id'],target=r['exact_group'],relationship='sequence_group') for r in regions]
    edges += [dict(source=r['region_id'],target=r['gene_family'],relationship='encodes_family') for r in gene_rows]
    edges += [dict(source=p['query'],target=p['target'],relationship=p['relationship']) for p in pairs]
    write_tsv(out/'network_edges.tsv',edges,['source','target','relationship'])
    metadata={s['sample_id']:s for s in samples}
    prevalence=[]
    for g in groups:
        ids={s for s,grp in called|found if grp==g}
        prevalence.append(dict(exact_group=g,n_assemblies=len(ids),total_assemblies=len(samples),
            n_subjects=len({metadata[s]['subject_id'] for s in ids if metadata[s]['subject_id']}),
            n_named_strains=len({(metadata[s]['subject_id'],metadata[s]['strain_id']) for s in ids if metadata[s]['subject_id'] and metadata[s]['strain_id']})))
    write_tsv(out/'group_prevalence.tsv',prevalence,['exact_group','n_assemblies','total_assemblies','n_subjects','n_named_strains'])
    paths=[args.samples]+args.calls+[s['assembly'] for s in samples]+([args.proteins,args.gene_map] if args.proteins else [])
    params={k:v for k,v in vars(args).items() if k not in ('func',)}
    warnings=['Not biologically validated; caller support is not a calibrated probability.',
        'Tracking uses user-supplied lineage evidence and exact unique flanks, not automatic strain inference.',
        'Exact sequence recovery alone does not establish integration, induction, transmission or direction of movement.',
        'No read-based evidence is evaluated in this release. Non-exact recovery in caller-negative assemblies is not implemented.']
    if not observations:warnings.append('No eligible longitudinal loci: inspect metadata, flags and flank availability.')
    versions={'python':platform.python_version()}
    import Bio
    versions['biopython']=Bio.__version__
    try:
        import pyrodigal
        versions['pyrodigal']=pyrodigal.__version__
    except ImportError:pass
    if shutil.which('mmseqs'):
        versions['mmseqs']=subprocess.check_output(['mmseqs','version'],text=True).strip()
    manifest=dict(version=__version__,created_utc=datetime.now(timezone.utc).isoformat(),parameters=params,
        dataset_label=args.label,gene_status='skipped' if args.skip_genes else args.gene_mode,
        software=versions,inputs=[dict(path=str(Path(p).resolve()),sha256=digest(p)) for p in sorted(set(paths))],
        counts=dict(samples=len(samples),regions=len(regions),events=len(events),gene_families=len(family_rows)),warnings=warnings)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    render(out/'report.html',samples,regions,events,family_rows,manifest,recovered)
    (out/'COMPLETE').write_text('Analysis completed successfully. See manifest.json for scope and limitations.\n')
    return manifest

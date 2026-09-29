"""Deterministic, entirely synthetic assemblies and predictions. No biological claims."""
import random
from pathlib import Path
from .io import write_fasta,write_tsv
from .regions import CALL_FIELDS


def make_demo(out):
    out=Path(out)
    if out.exists() and any(out.iterdir()):raise ValueError('Demo destination must be empty')
    out.mkdir(parents=True,exist_ok=True)
    rng=random.Random(1701)
    dna=lambda n:''.join(rng.choice('ACGT') for _ in range(n))
    codons=dict(zip('ACDEFGHIKLMNPQRSTVWY',['GCT','TGT','GAT','GAA','TTT','GGT','CAT','ATT','AAA','CTG','ATG','AAT','CCT','CAA','CGT','TCT','ACT','GTT','TGG','TAT']))
    proteins=['M'+''.join(rng.choice(list(codons)) for _ in range(119)) for _ in range(5)]
    def prophage(indices):return ''.join(''.join(codons[a] for a in proteins[i])+'TAA'+'ACGTACGTACGT' for i in indices)
    p=prophage([0,1,2,3]);q=prophage([0,1,2,4]);left=dna(1200);right=dna(1200)
    specs=[('A0',p,'2025-01-01','person1','strainA'),('A1',p,'2025-02-01','person1','strainA'),
        ('A2','','2025-03-01','person1','strainA'),('A3',p,'2025-04-01','person1','strainA'),('B0',q,'2025-01-01','person2','strainB')]
    samples=[];calls=[];mapping=[];faa={}
    for sid,phage,day,person,strain in specs:
        write_fasta(out/f'{sid}.fna',{'chromosome':left+phage+right})
        samples.append(dict(sample_id=sid,assembly=f'{sid}.fna',subject_id=person,strain_id=strain,date=day,lineage_evidence='synthetic_known_lineage'))
        if sid in ('A0','A1','B0'):
            for caller in ('phispy','phageboost'):
                calls.append(dict(sample_id=sid,contig='chromosome',start=1201,end=1200+len(phage),caller=caller,call_id='synthetic1'))
    for rid,indices in [('PT000001',[0,1,2,3]),('PT000002',[0,1,2,3]),('PT000003',[0,1,2,4])]:
        for j,i in enumerate(indices):
            pid=f'{rid}_g{j+1}';faa[pid]=proteins[i]
            mapping.append(dict(protein_id=pid,region_id=rid,start=j*375+1,end=j*375+363,strand='+'))
    write_tsv(out/'samples.tsv',samples,['sample_id','assembly','subject_id','strain_id','date','lineage_evidence'])
    write_tsv(out/'calls.tsv',calls,CALL_FIELDS)
    write_fasta(out/'proteins.faa',faa)
    write_tsv(out/'gene_map.tsv',mapping,['protein_id','region_id','start','end','strand'])
    (out/'README.txt').write_text('Entirely synthetic sequences, proteins, strain labels and caller predictions. A3 intentionally has no caller prediction. Expected: persistence A0->A1, loss A1->A2, gain A2->A3; related B0 shares three of four genes. This tests software logic, not caller accuracy.\n')

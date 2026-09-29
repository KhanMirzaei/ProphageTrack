"""Explicit execution of locally installed callers, without network setup."""
import json
import re
import shutil
import subprocess
import time
from pathlib import Path
from Bio import SeqIO
from .io import fasta,write_fasta,write_tsv,digest
from .adapters import import_calls
from .regions import CALL_FIELDS,validate_calls


def detect(args):
    out=Path(args.out).resolve()
    if out.exists() and any(out.iterdir()):raise ValueError('Caller output directory must be empty')
    if args.threads<1:raise ValueError('threads must be >=1')
    seqs=fasta(args.assembly)
    if any(not re.fullmatch(r'[A-Za-z0-9_.:-]+', name) for name in seqs):
        raise ValueError('Caller execution requires simple contig IDs (letters, numbers, underscore, period, colon or hyphen)')
    exe=args.executable or {'phispy':'PhiSpy.py','phageboost':'PhageBoost','vibrant':'VIBRANT_run.py'}[args.caller]
    if not shutil.which(exe):raise ValueError(f'Caller executable not found: {exe}')
    if args.caller=='phispy':
        if not args.genbank:raise ValueError('PhiSpy requires --genbank with annotated CDS features; annotate the assembly first')
        records=list(SeqIO.parse(args.genbank,'genbank'))
        if {r.id:str(r.seq).upper() for r in records}!=seqs:raise ValueError('GenBank sequences/IDs must exactly match the assembly')
        if not any(f.type=='CDS' for r in records for f in r.features):raise ValueError('GenBank has no CDS annotations')
    if args.caller=='phageboost' and not args.coordinates:
        raise ValueError('PhageBoost import requires explicit --coordinates after verifying the installed version; see docs/CALLERS.md')
    out.mkdir(parents=True,exist_ok=True)
    # Relative safe filenames protect callers that internally construct shell strings.
    write_fasta(out/'input.fna',seqs)
    raw=out/'raw'
    if args.caller=='phispy':
        shutil.copyfile(args.genbank,out/'input.gbk')
        cmd=[exe,'input.gbk','-o','raw','--threads',str(args.threads)]
    elif args.caller=='phageboost':cmd=[exe,'-f','input.fna','-o','raw','-j',str(args.threads)]
    else:
        cmd=[exe,'-i','input.fna','-folder','raw','-t',str(args.threads),'-no_plot']
        if args.database:
            database=str(Path(args.database).resolve())
            if any(c.isspace() for c in database) or any(c in database for c in "'\";`$|&<>()"):
                raise ValueError('VIBRANT database path must not contain shell-special characters')
            cmd+=['-d',database]
    started=time.monotonic()
    with open(out/'caller.log','w') as log:
        subprocess.run(cmd,cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
    pattern={'phispy':'prophage_coordinates.tsv','phageboost':'phages_*.gff','vibrant':'*integrated_prophage_coordinates*.tsv'}[args.caller]
    paths=sorted(raw.rglob(pattern)) if raw.exists() else []
    if not paths:
        # Missing files must not silently become "no prophages". Require manual
        # review of the caller log, including callers that exit 0 on failure.
        raise ValueError('Caller produced no coordinate file. Inspect caller.log; import an explicit empty normalized table only after confirming a successful no-hit run.')
    if len(paths)!=1:raise ValueError(f'Expected one coordinate file, found {len(paths)}; import explicitly')
    rows=import_calls(paths[0],args.caller,args.sample_id,args.coordinates or '1-inclusive',
        (args.contig_column,args.start_column,args.end_column) if args.contig_column else None)
    validate_calls(rows,{args.sample_id:seqs})
    write_tsv(out/'calls.tsv',rows,CALL_FIELDS)
    (out/'caller_manifest.json').write_text(json.dumps(dict(caller=args.caller,caller_version=args.caller_version or 'not_recorded',
        command=cmd,wall_seconds=time.monotonic()-started,assembly_sha256=digest(args.assembly),
        coordinates=args.coordinates or '1-inclusive',database=args.database or 'caller_default',
        database_version=args.database_version or 'not_recorded'),indent=2)+'\n')
    return len(rows)

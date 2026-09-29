import argparse
import json
import sys
from pathlib import Path
from . import __version__
from .io import write_tsv
from .regions import CALL_FIELDS


def parser():
    p=argparse.ArgumentParser(prog='prophagetrack',description='Shared prophages and conservative longitudinal tracking. Research preview.')
    p.add_argument('--version',action='version',version=__version__)
    sub=p.add_subparsers(dest='command',required=True)
    r=sub.add_parser('run',help='Analyze normalized predictions and assemblies')
    r.add_argument('--samples',required=True);r.add_argument('--calls',required=True,nargs='+');r.add_argument('--out',required=True)
    r.add_argument('--related',choices=['mmseqs','none'],default='mmseqs')
    r.add_argument('--identity',type=float,default=.95);r.add_argument('--coverage',type=float,default=.85)
    r.add_argument('--overlap',type=float,default=.8);r.add_argument('--flank-length',type=int,default=200)
    r.add_argument('--max-empty-gap',type=int,default=50);r.add_argument('--threads',type=int,default=1)
    r.add_argument('--proteins');r.add_argument('--gene-map');r.add_argument('--skip-genes',action='store_true')
    r.add_argument('--gene-mode',choices=['mmseqs','exact'],default='mmseqs')
    r.add_argument('--protein-identity',type=float,default=.5);r.add_argument('--protein-coverage',type=float,default=.8)
    r.add_argument('--label',default='User-supplied data')
    d=sub.add_parser('demo',help='Generate synthetic inputs, not real caller predictions');d.add_argument('--out',required=True)
    a=sub.add_parser('import-calls',help='Normalize an existing caller coordinate file')
    a.add_argument('--input',required=True);a.add_argument('--out',required=True)
    a.add_argument('--sample-id',required=True);a.add_argument('--caller',required=True,choices=['phispy','phageboost','vibrant','generic'])
    a.add_argument('--coordinates',required=True,choices=['1-inclusive','0-half-open'])
    for option in ['contig-column','start-column','end-column']:a.add_argument('--'+option)
    c=sub.add_parser('detect',help='Run one installed caller locally; no downloads')
    c.add_argument('--caller',required=True,choices=['phispy','phageboost','vibrant'])
    c.add_argument('--sample-id',required=True);c.add_argument('--assembly',required=True);c.add_argument('--genbank')
    c.add_argument('--out',required=True);c.add_argument('--threads',type=int,default=1);c.add_argument('--executable')
    c.add_argument('--coordinates',choices=['1-inclusive','0-half-open'])
    c.add_argument('--database');c.add_argument('--caller-version');c.add_argument('--database-version')
    for option in ['contig-column','start-column','end-column']:c.add_argument('--'+option)
    b=sub.add_parser('evaluate',help='Compare predictions with independent coordinate truth')
    b.add_argument('--truth',required=True);b.add_argument('--predictions',required=True);b.add_argument('--out',required=True)
    return p

def main(argv=None):
    p=parser();args=p.parse_args(argv)
    try:
        if args.command in ('detect','import-calls'):
            cols=[args.contig_column,args.start_column,args.end_column]
            if any(cols) and not all(cols):raise ValueError('Provide all three column-name options together')
        if args.command=='run':
            from .pipeline import run
            manifest=run(args);print(json.dumps(manifest['counts']));print(f'Report: {Path(args.out).resolve() / "report.html"}')
        elif args.command=='demo':
            from .demo import make_demo
            make_demo(args.out);print(f'Synthetic inputs: {Path(args.out).resolve()}')
        elif args.command=='import-calls':
            from .adapters import import_calls
            if Path(args.out).exists():raise ValueError('Refusing to overwrite an existing output')
            rows=import_calls(args.input,args.caller,args.sample_id,args.coordinates,tuple(cols) if all(cols) else None)
            write_tsv(args.out,rows,CALL_FIELDS);print(f'Imported {len(rows)} predictions; run validates them against the assembly')
        elif args.command=='detect':
            from .callers import detect
            print(f'Imported {detect(args)} predictions')
        else:
            from .benchmark import evaluate
            if Path(args.out).exists():raise ValueError('Refusing to overwrite an existing output')
            Path(args.out).write_text(json.dumps(evaluate(args.truth,args.predictions),indent=2)+'\n')
    except (ValueError,OSError,RuntimeError) as e:
        p.exit(2,f'prophagetrack: {e}\n')
    except Exception as e:
        # External program errors should be actionable without obscuring failures.
        import subprocess
        if isinstance(e,subprocess.CalledProcessError):p.exit(2,f'External command failed ({e.returncode}): {e.cmd}. Inspect the output log.\n')
        raise

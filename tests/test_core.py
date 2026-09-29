import copy
import json
import tempfile
import unittest
from pathlib import Path
from prophagetrack.io import fasta,write_tsv,write_fasta,load_samples
from prophagetrack.regions import validate_calls,reconcile,CALL_FIELDS
from prophagetrack.compare import reverse_complement,exact_groups,related_pairs
from prophagetrack.tracking import project,track,recover_exact
from prophagetrack.adapters import import_calls
from prophagetrack.genes import families,shared_blocks
from prophagetrack.benchmark import evaluate
from prophagetrack.demo import make_demo
from prophagetrack.pipeline import run
from prophagetrack.cli import parser


def call(start,end,caller='a',cid='1',sid='s'):
    return dict(sample_id=sid,contig='c',start=start,end=end,caller=caller,call_id=cid)

class RegionTests(unittest.TestCase):
    def setUp(self):self.assemblies={'s':{'c':'ACGT'*100}}
    def test_coordinate_bounds(self):
        for a,b in [(0,10),(10,9),(1,401)]:
            with self.assertRaises(ValueError):validate_calls([call(a,b)],self.assemblies)
    def test_duplicate_call(self):
        with self.assertRaises(ValueError):validate_calls([call(10,20),call(40,50)],self.assemblies)
    def test_unknown_contig(self):
        r=call(1,10);r['contig']='missing'
        with self.assertRaises(ValueError):validate_calls([r],self.assemblies)
    def test_agreement(self):
        r=reconcile([call(20,80),call(20,80,'b')],self.assemblies)
        self.assertEqual(len(r),1);self.assertEqual(r[0]['n_callers'],2);self.assertEqual(r[0]['sequence'],self.assemblies['s']['c'][19:80])
    def test_disagreeing_boundaries(self):
        r=reconcile([call(20,100),call(22,99,'b')],self.assemblies)
        self.assertIn('boundary_disagreement',r[0]['flags']);self.assertIn((r[0]['start'],r[0]['end']),[(20,100),(22,99)])
    def test_no_chain_merging(self):
        r=reconcile([call(20,119),call(40,139,'b'),call(60,159,'c')],self.assemblies)
        self.assertEqual(len(r),3);self.assertTrue(all('overlap_ambiguity' in x['flags'] for x in r))
    def test_same_caller_not_merged(self):
        self.assertEqual(len(reconcile([call(20,100),call(22,99,cid='2')],self.assemblies)),2)
    def test_edge_flag(self):self.assertIn('contig_edge',reconcile([call(1,100)],self.assemblies)[0]['flags'])
    def test_exact_reverse(self):
        a='ACGTTGCAAA';r=exact_groups([dict(region_id='a',sequence=a),dict(region_id='b',sequence=reverse_complement(a))])
        self.assertEqual(r[0]['exact_group'],r[1]['exact_group'])
    def test_ambiguous_not_exact(self):
        r=exact_groups([dict(region_id='a',sequence='ACNT'),dict(region_id='b',sequence='ACNT')])
        self.assertNotEqual(r[0]['exact_group'],r[1]['exact_group'])
    def test_fragment_classification(self):
        regions=[dict(region_id='a',exact_group='a'),dict(region_id='b',exact_group='b')]
        h=dict(query='a',target='b',identity=.99,query_coverage=1.,target_coverage=.2,alignment_length=100,qstart=1,qend=100,tstart=1,tend=100)
        self.assertEqual(related_pairs(regions,[h])[0]['relationship'],'fragment_match')
    def test_pair_direction(self):
        regions=[dict(region_id='a',exact_group='a'),dict(region_id='b',exact_group='b')]
        h=dict(query='b',target='a',identity=.99,query_coverage=1.,target_coverage=.2,alignment_length=100,qstart=1,qend=100,tstart=5,tend=104)
        r=related_pairs(regions,[h])[0];self.assertEqual(r['query_coverage'],.2);self.assertEqual(r['qstart'],5)

class SiteTests(unittest.TestCase):
    def setUp(self):self.l='AACCGTAC';self.r='TTCGAGGA';self.p='G'*80
    def test_occupied(self):self.assertEqual(project(self.l,self.r,self.p,{'c':self.l+self.p+self.r})['state'],'occupied')
    def test_empty(self):self.assertEqual(project(self.l,self.r,self.p,{'c':self.l+self.r})['state'],'candidate_empty')
    def test_reverse_site(self):self.assertEqual(project(self.l,self.r,self.p,{'c':reverse_complement(self.l+self.p+self.r)})['state'],'occupied')
    def test_repeated_flank(self):self.assertEqual(project(self.l,self.r,self.p,{'c':self.l+self.p+self.r+self.l})['state'],'unresolved')
    def test_split_contigs(self):self.assertEqual(project(self.l,self.r,self.p,{'c':self.l,'d':self.r})['state'],'unresolved')
    def test_gap_ns(self):self.assertEqual(project(self.l,self.r,self.p,{'c':self.l+'NN'+self.r})['state'],'unresolved')
    def test_other_occupancy(self):self.assertEqual(project(self.l,self.r,self.p,{'c':self.l+'C'*90+self.r})['state'],'different_occupancy')
    def test_missing_flank(self):self.assertEqual(project(self.l,self.r,self.p,{'c':self.l+self.p})['state'],'unresolved')
    def test_exact_recovery(self):
        r=exact_groups([dict(region_id='r',sequence=self.p)])
        rows=recover_exact(r,{'s':{'c':self.l+self.p+self.r},'no_call':{'c':self.p}})
        self.assertEqual({x['sample_id'] for x in rows},{'s','no_call'})

class FileTests(unittest.TestCase):
    def setUp(self):self.temp=tempfile.TemporaryDirectory();self.p=Path(self.temp.name)
    def tearDown(self):self.temp.cleanup()
    def test_duplicate_fasta(self):
        f=self.p/'a.fa';f.write_text('>a\nACGT\n>a\nACGT\n')
        with self.assertRaises(ValueError):fasta(f)
    def test_bad_fasta(self):
        f=self.p/'a.fa';f.write_text('>a\nAC-Q\n')
        with self.assertRaises(ValueError):fasta(f)
    def test_phispy_adapter(self):
        f=self.p/'a.tsv';f.write_text('pp1\tc\t10\t20\textra\n')
        self.assertEqual(import_calls(f,'phispy','s')[0]['end'],20)
    def test_phageboost_adapter_conversion(self):
        f=self.p/'a.gff';f.write_text('##gff-version 3\nc\tPhageBoost\tprophage_prediction\t0\t20\t.9\t.\t.\tphage_id=p1\n')
        self.assertEqual(import_calls(f,'phageboost','s','0-half-open')[0]['start'],1)
    def test_vibrant_explicit_columns(self):
        f=self.p/'v.tsv';f.write_text('host\tbegin\tfinish\nc\t3\t80\n')
        self.assertEqual(import_calls(f,'vibrant','s',columns=('host','begin','finish'))[0]['end'],80)
    def test_benchmark_union(self):
        t=self.p/'t.tsv';p=self.p/'p.tsv'
        write_tsv(t,[call(1,100)],CALL_FIELDS);write_tsv(p,[call(1,50),call(40,110,cid='2')],CALL_FIELDS)
        r=evaluate(t,p);self.assertEqual(r['true_positive_bp'],100);self.assertEqual(r['false_positive_bp'],10)
    def test_gene_exact_families(self):
        rows=[dict(protein_id='a',region_id='r1'),dict(protein_id='b',region_id='r2'),dict(protein_id='c',region_id='r2')]
        genes,fams=families({'a':'MALALA','b':'MALALA','c':'MWWWW'},rows,self.p,mode='exact')
        self.assertEqual(len(fams),2);self.assertEqual(genes[0]['gene_family'],genes[1]['gene_family'])
    def test_shared_blocks_reverse(self):
        rows=[dict(region_id=rid,start=i,gene_family=g) for rid,ids in [('a',['x','y','z']),('b',['z','y','x'])] for i,g in enumerate(ids)]
        self.assertEqual(shared_blocks(rows)[0]['n_regions'],2)

class IntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.p=Path(self.temp.name);self.demo=self.p/'demo';make_demo(self.demo)
    def tearDown(self):self.temp.cleanup()
    def args(self):
        return parser().parse_args(['run','--samples',str(self.demo/'samples.tsv'),'--calls',str(self.demo/'calls.tsv'),
            '--proteins',str(self.demo/'proteins.faa'),'--gene-map',str(self.demo/'gene_map.tsv'),
            '--related','none','--gene-mode','exact','--out',str(self.p/'results'),'--label','SYNTHETIC'])
    def test_end_to_end(self):
        from prophagetrack.io import read_tsv
        m=run(self.args());self.assertEqual(m['counts']['regions'],3)
        events=read_tsv(self.p/'results/events.tsv');self.assertEqual([e['event'] for e in events],['persistence','candidate_loss','candidate_gain'])
        self.assertTrue((self.p/'results/COMPLETE').exists())
        recovery=read_tsv(self.p/'results/exact_recovery.tsv');self.assertIn('A3',{r['sample_id'] for r in recovery})
        self.assertIn('SYNTHETIC',(self.p/'results/report.html').read_text())
    def test_no_overwrite(self):
        run(self.args())
        with self.assertRaises(ValueError):run(self.args())
    def test_missing_lineage_no_events(self):
        from prophagetrack.io import read_tsv
        rows=read_tsv(self.demo/'samples.tsv')
        for r in rows:r['lineage_evidence']=''
        write_tsv(self.demo/'samples.tsv',rows,list(rows[0]))
        run(self.args());self.assertEqual(read_tsv(self.p/'results/events.tsv'),[])
    def test_same_date_not_event(self):
        from prophagetrack.io import read_tsv
        rows=read_tsv(self.demo/'samples.tsv');rows[1]['date']=rows[0]['date']
        write_tsv(self.demo/'samples.tsv',rows,list(rows[0]))
        run(self.args());self.assertEqual(read_tsv(self.p/'results/events.tsv')[0]['event'],'unresolved_same_date')
    def test_html_escaped(self):
        a=self.args();a.label='<script>alert(1)</script>';run(a)
        html=(self.p/'results/report.html').read_text();self.assertNotIn('<script>',html);self.assertIn('&lt;script&gt;',html)
    def test_out_of_range_parameter(self):
        a=self.args();a.coverage=1.1
        with self.assertRaises(ValueError):run(a)
    def test_empty_predictions(self):
        from prophagetrack.io import read_tsv
        write_tsv(self.demo/'calls.tsv',[],CALL_FIELDS)
        a=self.args();a.skip_genes=True;a.proteins=None;a.gene_map=None
        m=run(a);self.assertEqual(m['counts']['regions'],0)

if __name__=='__main__':unittest.main()

class CallerRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.p=Path(self.temp.name)
        write_fasta(self.p/'assembly.fna',{'c':'ACGT'*100})
    def tearDown(self):self.temp.cleanup()
    def executable(self,body):
        import sys
        path=self.p/'mock_caller'
        path.write_text('#!'+sys.executable+'\n'+body);path.chmod(0o755)
        return str(path)
    def args(self,exe):
        return parser().parse_args(['detect','--caller','phageboost','--assembly',str(self.p/'assembly.fna'),
            '--sample-id','S1','--out',str(self.p/'out'),'--executable',exe,'--coordinates','1-inclusive','--caller-version','mock-fixture'])
    def test_mock_runner_manifest_and_calls(self):
        from prophagetrack.callers import detect
        exe=self.executable("from pathlib import Path\nPath('raw').mkdir()\nPath('raw/phages_input.gff').write_text('c\\tPhageBoost\\tprophage_prediction\\t10\\t100\\t0.9\\t.\\t.\\tphage_id=p1\\n')\n")
        self.assertEqual(detect(self.args(exe)),1)
        self.assertEqual(json.loads((self.p/'out/caller_manifest.json').read_text())['caller_version'],'mock-fixture')
    def test_exit_zero_without_output_fails(self):
        from prophagetrack.callers import detect
        with self.assertRaises(ValueError):detect(self.args(self.executable('print("configuration failure")\n')))
    def test_missing_coordinate_convention_fails(self):
        from prophagetrack.callers import detect
        a=self.args(self.executable('pass\n'));a.coordinates=None
        with self.assertRaises(ValueError):detect(a)

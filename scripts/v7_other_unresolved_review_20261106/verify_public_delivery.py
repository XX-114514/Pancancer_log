"""Validate the real, aggregate-only review delivery and explicit SVG row/panel labels."""
import argparse,json,hashlib,re,ast,xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd,numpy as np
p=argparse.ArgumentParser();p.add_argument('result_root',type=Path);a=p.parse_args();r=a.result_root
models={'Erythroid':(9098,18),'Platelet':(2318,11),'HSPC':(14072,16),'ILC':(45290,14)}
manifest=json.loads((r/'artifact_manifest.json').read_text())
for z in manifest['files']:
 f=r/z['path'];assert f.stat().st_size==z['bytes'] and hashlib.sha256(f.read_bytes()).hexdigest()==z['sha256']
allrows=pd.read_csv(r/'all_59_cluster_review_with_evidence.tsv',sep='\t');assert len(allrows)==59 and not allrows.duplicated(['lineage','cluster']).any()
assert int(allrows.n_cells.sum())==70778 and int(allrows.n_identity_review_required.sum())==43704
markers=0;svg_files=0
for name,(n,k) in models.items():
 d=r/name;t=pd.read_csv(d/'markers/measured_marker_expression.tsv',sep='\t');lab=pd.read_csv(d/(name+'_cluster_decisions.tsv'),sep='\t');spec=json.loads((d/(name+'_plot_panels.json')).read_text())
 assert len(t)==260*k and len(lab)==k and not t.duplicated(['cluster','gene']).any()
 assert t.groupby('cluster').n_cells.first().sum()==n and t.groupby('cluster').n_cells.nunique().eq(1).all()
 assert t.n_expressing.between(0,t.n_measured).all() and t.n_measured.between(0,t.n_cells).all()
 assert np.allclose(t.coverage_fraction,t.n_measured/t.n_cells)
 assert set(t.gene)==set(pd.read_csv(r/'Erythroid/markers/measured_marker_expression.tsv',sep='\t').gene)
 assert not lab.final_fine_identity_accepted.any()
 assert json.loads((d/'markers/independent_postfit_audit.json').read_text())['status']=='PASS'
 assert json.loads((d/'model/terminal.json').read_text())['exit_code']==0
 res=pd.read_csv(d/'model/cluster_resolution_counts.tsv',sep='\t');assert res.groupby('graph').n_cells.sum().eq(n).all()
 plots=json.loads((d/'plots/plot_manifest.json').read_text());assert plots['n_cells']==n and plots['n_clusters']==k
 assert len(list((d/'plots').glob('*.png')))==4 and len(list((d/'plots').glob('*.svg')))==4
 for z in plots['outputs']:assert hashlib.sha256((d/'plots'/z['name']).read_bytes()).hexdigest()==z['sha256']
 for f in (d/'plots').glob('*.svg'):
  text=' '.join(ET.parse(f).getroot().itertext());flat=re.sub(r'\s+',' ',text)
  for row in lab.itertuples():assert row.display_annotation in flat and f'n={int(t[t.cluster.eq(row.cluster)].n_cells.iloc[0]):,}' in flat
  assert 'source coverage <80%' in flat and 'not accepted fine subtypes' in flat
  svg_files+=1
 markers+=len(t)
for f in r.rglob('*.tsv'):
 h=f.read_text().splitlines()[0].split('\t');assert not set(h)&{'cell_id','patient_id','patientID','sampleID','source_obs_name','barcode'}
for f in r.rglob('*'):
 if f.is_file() and f.suffix not in {'.png'}:
  assert not re.search(r'/data\d/|/Users/[^/]+/|BE[G]IN (?:OPENSSH|RSA) PRIVATE KEY|gh[pousr]_[A-Za-z0-9]{20,}',f.read_text())
for f in Path(__file__).parent.glob('*.py'):ast.parse(f.read_text(),filename=str(f))
print(json.dumps({'status':'PASS','n_clusters':len(allrows),'n_cells':int(allrows.n_cells.sum()),'n_marker_records':markers,'n_SVG_row_and_caption_checks':svg_files,'n_PNG_figures':16,'all_manifest_hashes_match':True,'no_cell_patient_sample_row_keys_in_public_TSV':True,'new_Python_AST_parse':'PASS'}))

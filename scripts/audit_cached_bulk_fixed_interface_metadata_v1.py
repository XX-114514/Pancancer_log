"""Independent complete feature-axis reads and clinical aggregate count replay."""
from pathlib import Path
import os,json,gzip,csv,hashlib,math
from datetime import datetime,timezone
from collections import Counter
import numpy as np,pandas as pd

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def nval(v):
 try:return float(v)
 except (TypeError,ValueError):return float('nan')
def tabcounts(path,tcga=False):
 with path.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 key=list(rows[0])[0] if rows else None
 if rows:assert len({r[key] for r in rows})==len(rows)
 counts={}
 for ep in (['OS','PFI'] if tcga else ['OS','PFS']):
  ee=ep if tcga else ep+'_event';tt=ep+'.time' if tcga else ep
  pairs=[(nval(r[tt]),nval(r[ee])) for r in rows]
  valid=[e for t,e in pairs if math.isfinite(t) and t>0 and e in [0,1]]
  counts[ep]=dict(known_positive_time_units=len(valid),events=int(sum(valid)))
 if tcga:return dict(cancer=path.name.split('_')[1],mapped_primary_tumor_units=len(rows),endpoint_counts=counts)
 y=[nval(r['response']) for r in rows]
 return dict(old_mapped_units=len(rows),response_known=sum(v in [0,1] for v in y),response_positive=y.count(1.),response_negative=y.count(0.),endpoint_counts=counts)
def main():
 root=Path(os.environ['PANCANCER_BULK_ROOT']);old=Path(os.environ['PANCANCER_OLD_BULK_RUN']);out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
 terminal=json.loads((out/'private_terminal.json').read_text())
 for n,v in terminal['files'].items():assert (out/n).stat().st_size==v['bytes'] and sha(out/n)==v['SHA256']
 axes=json.loads((out/'private_all7_physical_bulk_feature_axis_qualification.json').read_text())
 c=json.loads((out/'private_all6_ICI_and33_TCGA_old_clinical_unit_counts.json').read_text())
 paths=[root/'ici_bulk/imvigor210/data_mrna_seq_tpm.txt',root/'ici_bulk/riaz2017/GSE91061_raw_counts.csv.gz',root/'ici_bulk/hugo2016/GSE78220_FPKM.xlsx',root/'ici_bulk/jung2019_nsclc/GSE135222_exp.tsv.gz',root/'ici_bulk/rose2021_urothelial/GSE176307_salmon_tpm.tsv.gz',old/'SU2C_MARK_panel_ready_v1/fixed_panel_TPM_private.tsv.gz',root/'tcga/tcga_RSEM_gene_tpm.gz']
 for row,p in zip(axes,paths):
  if p.suffix=='.xlsx':a=pd.read_excel(p,usecols=[0],dtype=str).iloc[:,0].astype(str).tolist()
  else:
   chunks=pd.read_csv(p,sep=',' if 'raw_counts.csv' in p.name else '\t',usecols=[0],dtype=str,keep_default_na=False,chunksize=4096)
   a=[str(v) for ch in chunks for v in ch.iloc[:,0]]
  assert len(a)==row['raw_axis_rows']
  assert hashlib.sha256('\n'.join(a).encode()).hexdigest()==row['physical_axis_SHA256']
  assert p.stat().st_size==row['raw_file_bytes'] and sha(p)==row['raw_file_SHA256']
  del a
  if row['cohort']=='TCGA':
   cc=[tabcounts(p,True) for p in sorted((old/'analysis_v2').glob('TCGA_*_patient_scores_private.tsv'))]
   assert cc==next(v for v in c if v['cohort']=='TCGA')['source_cancer_counts']
  else:
   name=row['cohort'];path=old/('SU2C_analysis_v1' if name=='SU2C_MARK' else 'analysis_v2')/(name+'_patient_scores_private.tsv')
   got=tabcounts(path);want=next(v for v in c if v['cohort']==name)
   assert all(want[k]==v for k,v in got.items())
 result=dict(status='ACCEPTED_ALL7_FULL_FEATURE_AXES_AND_ALL6_PLUS33_SOURCE_CLINICAL_COUNTS',completed_UTC=datetime.now(timezone.utc).isoformat(),audit_source_SHA256=sha(__file__),producer_terminal_SHA256=sha(out/'private_terminal.json'),
  raw_feature_axis_alternate_parser_complete=True,all_raw_files_size_SHA256_verified=True,all_clinical_source_counts_csv_vs_pandas_equal=True,
  qualification_not_biological_expression_or_patient_truth=True,no_expression_or_patient_rows_or_model_arrays_persisted=True,no_new_outcomes_or_effects=True)
 dest=out/'private_full_metadata_independent_audit.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result))
if __name__=='__main__':main()

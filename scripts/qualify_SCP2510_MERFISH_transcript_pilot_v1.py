"""Full public transcript pilot: conservative FOV/segmentation fragments, no P."""
import argparse,hashlib,json
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter,defaultdict
import numpy as np,pandas as pd
from scipy.sparse import coo_matrix,save_npz

EXCLUDED={'','0','0.0','-1','-1.0','nan','NaN','NA','None'}

def main(a):
 receipt=json.loads((a.source/'download_receipt_public.json').read_text());raw=a.source/'tx_merfish_Pt43.csv'
 md=hashlib.md5();ss=hashlib.sha256()
 with raw.open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):md.update(b);ss.update(b)
 assert md.hexdigest()==receipt['md5'] and ss.hexdigest()==receipt['sha256'] and raw.stat().st_size==receipt['bytes']
 assert not a.output.exists();a.output.mkdir(mode=0o700)
 panels=json.loads(a.plan.read_text())['panels'];gene_count=Counter();fovs=set();baysor_fov=defaultdict(set);cellpose_fov=defaultdict(set);sentinel=Counter()
 cells={};genes={};counts=Counter();xy_sum=defaultdict(lambda:np.zeros(3,float));totals=Counter();n=assigned=bad_coordinates=0
 low=np.full(3,np.inf);high=-low
 for chunk in pd.read_csv(raw,chunksize=100000,dtype={'gene':str,'fov':str,'cellpose_segmentation':str,'baysor_segmentation':str},keep_default_na=False):
  assert list(chunk.columns)==['global_x','global_y','global_z','fov','gene','cellpose_segmentation','baysor_segmentation']
  xyz=chunk[['global_x','global_y','global_z']].apply(pd.to_numeric,errors='coerce').to_numpy(float);finite=np.isfinite(xyz).all(1)
  bad_coordinates+=int((~finite).sum());low=np.minimum(low,np.nanmin(xyz,axis=0));high=np.maximum(high,np.nanmax(xyz,axis=0))
  assert finite.all() and chunk.gene.ne('').all();n+=len(chunk);gene_count.update(chunk.gene);fovs.update(chunk.fov)
  for segmentation,registry in [('baysor_segmentation',baysor_fov),('cellpose_segmentation',cellpose_fov)]:
   for cell,fov in chunk[[segmentation,'fov']].drop_duplicates().itertuples(index=False,name=None):
    if cell not in EXCLUDED:registry[cell].add(fov)
  invalid=chunk.baysor_segmentation.isin(EXCLUDED);sentinel.update(chunk.loc[invalid,'baysor_segmentation']);selected=chunk.loc[~invalid].copy();assigned+=len(selected)
  # FOV is always part of the key: avoid merging reused labels. This may split
  # real cross-FOV cells; fragments must not be promoted to independent cells.
  for (fov,label),frame in selected.groupby(['fov','baysor_segmentation'],sort=False):
   key=(fov,label)
   if key not in cells:cells[key]=len(cells)
   i=cells[key];xy_sum[i]+=frame[['global_x','global_y','global_z']].to_numpy(float).sum(0);totals[i]+=len(frame)
   for gene,value in frame.gene.value_counts(sort=False).items():
    if gene not in genes:genes[gene]=len(genes)
    counts[i,genes[gene]]+=int(value)
  if n//1000000!=(n-len(chunk))//1000000:print(json.dumps(dict(phase='full_transcript_scan',rows=n)),flush=True)
 row=np.fromiter((k[0] for k in counts),np.int32);col=np.fromiter((k[1] for k in counts),np.int32);val=np.fromiter(counts.values(),np.int64)
 x=coo_matrix((val,(row,col)),shape=(len(cells),len(genes))).tocsr();assert int(x.sum())==assigned
 g=sorted(genes,key=genes.get);gene_index={v:i for i,v in enumerate(g)};tot=np.array([totals[i] for i in range(len(cells))],np.int64);centers=np.array([xy_sum[i]/totals[i] for i in range(len(cells))])
 assert np.array_equal(np.asarray(x.sum(1)).ravel(),tot) and np.isfinite(centers).all()
 save_npz(a.output/'private_all_observed_gene_fragment_counts.npz',x)
 np.savez_compressed(a.output/'private_fragment_axes_and_transcript_centroids.npz',fov=np.array([v[0] for v in cells],str),segmentation_label=np.array([v[1] for v in cells],str),gene=np.array(g,str),transcript_centroid_xyz=centers,total_transcripts=tot)
 pd.DataFrame(sorted(gene_count.items()),columns=['gene','observed_transcripts']).to_csv(a.output/'private_all_observed_gene_counts.tsv',sep='\t',index=False)
 coverage={k:dict(required=len(gg),observed_in_raw=len(set(gg)&set(gene_count)),not_observed_measurement_status_uncertain=sorted(set(gg)-set(gene_count))) for k,gg in panels.items()}
 target=set(panels['pan_T']+panels['regulatory_associated']+panels['cytotoxic']);missing=sorted(target-set(gene_count));complete=not missing
 fixed_counts={}
 def detected(gg):return np.asarray((x[:,[gene_index[z] for z in gg]]>0).sum(1)).ravel()
 if complete and target<=set(g):
  T=detected(panels['pan_T'])>=2;R=detected(panels['regulatory_associated'])>=2;C=detected(panels['cytotoxic'])>=2
  fixed_counts=dict(T_capture_fragments=int(T.sum()),R_exclusive_fragments=int((T&R&~C).sum()),C_exclusive_fragments=int((T&C&~R).sum()),double_positive_fragments=int((T&R&C).sum()))
 report=dict(status='FULL_SELECTED_PUBLIC_TRANSCRIPT_FILE_QUALIFIED_AS_TECHNICAL_FRAGMENT_INPUT_ONLY',finished_utc=datetime.now(timezone.utc).isoformat(),source_record=receipt['source_record'],source_file='tx_merfish_Pt43.csv',source_bytes=raw.stat().st_size,source_md5=md.hexdigest(),source_sha256=ss.hexdigest(),
  all_transcript_rows=n,all_observed_gene_labels=len(gene_count),n_FOV_codes=len(fovs),nonfinite_coordinate_rows=bad_coordinates,excluded_or_unverified_assignment_code_rows=sum(sentinel.values()),retained_transcript_rows=assigned,
  conservative_FOV_Baysor_fragments=len(cells),observed_assigned_genes=len(g),Baysor_labels_observed_in_multiple_FOVs=sum(len(v)>1 for v in baysor_fov.values()),Cellpose_labels_observed_in_multiple_FOVs=sum(len(v)>1 for v in cellpose_fov.values()),
  complete_fixed_CM02_14_genes_observed=complete,CM02_genes_not_observed_not_proven_unmeasured=missing,HLAII_group_genes_not_observed=sorted({'HLA-DRA','CD74'}-set(gene_count)),coverage=coverage,
  fragment_capture_counts_private=fixed_counts,spatial_or_clinical_effects_not_tested=True,independent_patient_validation=False,
  limits=['Only smallest-file technical pilot, not a four-patient analysis or representative effect selection','Zero/unassigned-code meaning and cross-FOV segmentation identity require author verification','Conservative FOV+label keys may split genuine cells; no fragment is counted as independent patient or validated cell identity','Coordinate means are transcript centroids, not nuclear/boundary cell centroids or contact','Unobserved genes may be unassayed or below capture; no zero-filled missing genes','Published pre-treatment cohort design does not by itself provide the four MERFISH cases matched clinical outcomes'])
 (a.output/'private_complete_pilot_qualification.json').write_text(json.dumps(report,indent=2)+'\n')
 public={k:v for k,v in report.items() if k not in ['fragment_capture_counts_private']};public['code_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 (a.output/'technical_qualification_public.json').write_text(json.dumps(public,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['status','all_transcript_rows','complete_fixed_CM02_14_genes_observed','CM02_genes_not_observed_not_proven_unmeasured','HLAII_group_genes_not_observed']}),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser()
 for name in ['source','plan','output']:p.add_argument('--'+name,type=Path,required=True)
 main(p.parse_args())

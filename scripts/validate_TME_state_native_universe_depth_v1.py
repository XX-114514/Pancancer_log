"""Development-only official UCell rank-universe and independent-gene probes.

Balanced probes are measurement diagnostics, never functional validation.
No held-out score arrays or outcomes are opened. No score is imputed as zero.
"""
import argparse,datetime,hashlib,json,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr

def load(p):return json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def bound(b):
 p=Path(b['path']);assert sha(p)==b['SHA256'];return p
def binding(p):return dict(path=str(p),SHA256=sha(p))
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def corr(x,y):
 if len(x)<3 or np.ptp(x)==0 or np.ptp(y)==0:return None
 v=float(spearmanr(x,y).statistic);return v if np.isfinite(v) else None

def main(R):
 cfg=load(R/'config/spec_v1.json');Q=Path(cfg['candidates']);assert load(Q/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 qp=Q/'outputs/accepted_development_state_measurement_candidates_v1.json';assert sha(qp)==load(Q/'provenance/terminal_v1.json')['acceptance_SHA256'];cand=load(qp)
 readouts=load(bound(cand['readout_acceptance']));app=load(bound(cand['applicability']));roles=load(bound(cand['source_roles']))
 audit=load(bound(readouts['score_numerical_audit']));idx=load(bound(audit['producer_index']));P=Path(cfg['score_producer']);pc=load(P/'config/spec_v1.json')
 J=Path(pc['joined']);H=Path(pc['model']);sys.path[:0]=[str(J/'scripts'),str(pc['vendor_parent'])]
 from native_archive_collection_v1 import NativeArchiveCollection
 from native_annotation_review_v1 import native_measured_genes
 from native_panel_aggregate_v1 import gather_csr_rows
 from audit_TME_standard_Harmony_v2 import gather_nonnegative_float_delta
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 import anndata,pyucell
 assert pyucell.__version__=='0.8.0';storage_gate(R,cfg)
 with np.load(H/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as f:axes={k:f[k] for k in f.files}
 totals=np.load(H/'outputs/private_analysis_RNA_totals_v1.npy',mmap_mode='r');ctx=np.load(bound(app['context_codes']))
 lineage={x['source']:x for x in load(H/'outputs/source_input_lineage_v1.json')};common=load(H/'outputs/scope_and_common_genes_v1.json')['genes'];commonset=set(common)
 native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest'])
 catalogue=load(bound(pc['catalog']));names=idx['program_ids'];genes={x['program_id']:sorted(set(x['genes'])) for x in catalogue['programs']};columns={v:i for i,v in enumerate(names)}
 role={r['source']:r['role'] for r in roles['source_roles']};allow=np.load(bound(roles['development_mask']));ap=np.load(bound(app['context_program_applicability']))
 wanted={(r['source'],int(r['context_code']),r['program_id']) for r in pd.read_csv(bound(cand['diagnostics']),sep='\t').to_dict('records') if r['measurement_candidate']}
 rng=np.random.default_rng(cfg['seed']);sources=[];diagnostics=[];probe_count=0
 for record in idx['sources']:
  source=record['source']
  if role[source]!='development':continue
  with np.load(bound(record['cell_axis']),allow_pickle=False) as f:rows=f['fit_row'];units=f['qualification_unit_code'];detected=f['detected_native_genes']
  assert allow[rows].all();contexts=ctx[rows];strata=[]
  for key in np.unique(units.astype(np.int64)*app['n_contexts']+contexts):
   pos=np.flatnonzero(units.astype(np.int64)*app['n_contexts']+contexts==key)
   if len(pos)<cfg['minimum_group_cells']:continue
   strata.append(rng.permutation(pos)[:cfg['maximum_cells_per_unit_context']])
  chosen=np.array([a[k] for k in range(cfg['maximum_cells_per_unit_context']) for a in strata if len(a)>k][:cfg['maximum_cells_per_source']],int)
  if not len(chosen):continue
  chosen=np.sort(chosen);rr=rows[chosen];cc=contexts[chosen];uu=units[chosen]
  gn,ph=native_measured_genes(native,source);gn=list(gn);gset=set(gn);pos={g:i for i,g in enumerate(gn)};assert len(pos)==len(gn) and commonset<=gset
  z=native.store(source);x=gather_csr_rows(z[lineage[source]['counts_key']],axes['native_row'][rr],max_span=8192).astype(np.float64).tocsr()
  if 'ambient_removed_counts' in z['layers']:x=(x-gather_nonnegative_float_delta(z['layers/ambient_removed_counts'],axes['native_row'][rr])).tocsr()
  x.eliminate_zeros();assert np.isfinite(x.data).all() and x.data.min(initial=0)>=-1e-8;x.data[x.data<0]=0;x.eliminate_zeros()
  np.testing.assert_allclose(np.asarray(x.sum(1)).ravel(),totals[rr],rtol=1e-10,atol=1e-7);x=x[:,ph].astype(np.float32).tocsr()
  full={};shared={};halfa={};halfb={}
  for j,n in enumerate(names):
   gg=genes[n];measured=sorted(set(gg)&gset)
   if len(measured)>=pc['minimum_genes'] and len(measured)/len(gg)>=pc['minimum_fraction']:full[n]=measured
   ggshared=sorted(set(gg)&commonset)
   if len(ggshared)>=pc['minimum_genes'] and len(ggshared)/len(gg)>=pc['minimum_fraction']:shared[n]=ggshared
   # Split the complete external gene list before intersecting measured genes.
   ordered=sorted(gg,key=lambda g:hashlib.sha256((n+'|'+g+'|fixed_gene_halves_v1').encode()).digest())
   a=ordered[::2];b=ordered[1::2];aa=sorted(set(a)&gset);bb=sorted(set(b)&gset)
   if min(len(aa),len(bb))>=pc['minimum_genes'] and len(aa)/len(a)>=pc['minimum_fraction'] and len(bb)/len(b)>=pc['minimum_fraction']:halfa[n]=aa;halfb[n]=bb
  scores={}
  settings=[('native_1500',x,gn,full,1500),('native_1000',x,gn,full,1000),('native_2000',x,gn,full,2000),('native_shared_signature_1500',x,gn,shared,1500),('common_shared_signature_1500',x[:,[pos[g] for g in common]].tocsr(),common,shared,1500),('disjoint_gene_A_1500',x,gn,halfa,1500),('disjoint_gene_B_1500',x,gn,halfb,1500)]
  for setting,matrix,gene_axis,sig,rank in settings:
   y=np.full((len(chosen),len(names)),np.nan,np.float32)
   if sig:
    a=anndata.AnnData(matrix,var=pd.DataFrame(index=gene_axis));pyucell.compute_ucell_scores(a,sig,max_rank=rank,missing_genes='skip',ties_method='average',chunk_size=cfg['score_chunk_size'],n_jobs=cfg['threads'],device='cpu')
    js=[columns[n] for n in sig];y[:,js]=a.obs[[n+'_UCell' for n in sig]].to_numpy(np.float32);assert np.isfinite(y[:,js]).all()
   scores[setting]=y
  baseline=np.load(bound(record['scores']),mmap_mode='r')[chosen];np.testing.assert_allclose(scores['native_1500'],baseline,atol=3e-6,rtol=1e-5,equal_nan=True)
  for c in np.unique(cc):
   use=np.flatnonzero(cc==c)
   for j,n in enumerate(names):
    if (source,int(c),n) not in wanted:continue
    yy=scores['native_1500'][use,j]
    assert np.isfinite(yy).all() and ap[c,j]
    comparisons=[('rank_1000',yy,scores['native_1000'][use,j]),('rank_2000',yy,scores['native_2000'][use,j]),('rank_universe_shared_signature',scores['native_shared_signature_1500'][use,j],scores['common_shared_signature_1500'][use,j]),('native_signature_rescue',yy,scores['native_shared_signature_1500'][use,j]),('disjoint_gene_halves',scores['disjoint_gene_A_1500'][use,j],scores['disjoint_gene_B_1500'][use,j])]
    for kind,a,b in comparisons:
     good=np.isfinite(a)&np.isfinite(b);rho=corr(a[good],b[good]) if good.any() else None;means=[]
     for u in np.unique(uu[use][good]):
      take=good&(uu[use]==u)
      if take.sum()>=cfg['minimum_probe_cells_per_unit']:means.append([a[take].mean(),b[take].mean()])
     means=np.asarray(means)
     diagnostics.append(dict(source=source,context_code=int(c),program_id=n,comparison=kind,n_probe_cells=int(good.sum()),n_probe_units=len(means),cell_Spearman=rho,unit_mean_Spearman=corr(means[:,0],means[:,1]) if len(means)>=cfg['minimum_probe_units'] else None,mean_absolute_difference=float(np.abs(a[good]-b[good]).mean()) if good.any() else None,RNA_Spearman=corr(yy,np.log1p(totals[rr[use]])),detected_gene_Spearman=corr(yy,detected[chosen[use]]),functional_validation_claimed=False))
  out=R/'outputs'/('private_'+source+'_native_state_sensitivity_probe_v1.npz');np.savez_compressed(out,fit_row=rr,qualification_unit_code=uu,context_code=cc,**scores)
  sources.append(dict(source=source,n_cells=len(chosen),probe=binding(out)));probe_count+=len(chosen);native._fields.clear();storage_gate(R,cfg)
  save(R/'checkpoints'/('source_'+str(len(sources))+'_complete_v1.json'),dict(n_sources=len(sources),n_probe_cells=probe_count))
 table=R/'outputs/private_development_state_native_sensitivity_v1.tsv';pd.DataFrame(diagnostics).to_csv(table,sep='\t',index=False)
 op=R/'outputs/accepted_development_state_native_sensitivity_v1.json';save(op,dict(schema='TME_development_state_rank_universe_disjoint_gene_diagnostics_v1',candidate_acceptance=binding(qp),source_roles=cand['source_roles'],applicability=cand['applicability'],n_sources=len(sources),n_probe_cells=probe_count,n_settings=7,n_diagnostic_records=len(diagnostics),sources=sources,diagnostics=binding(table),official_pyucell_version=pyucell.__version__,original_native_baseline_reproduced=True,source_roles_development_only=True,n_lockbox_arrays_opened=0,no_molecule_thinning_applied_to_fractional_ambient_corrected_counts=True,gene_disjoint_coherence_not_biological_function_validation=True,missing_coverage_or_unidentifiable_correlation_is_NaN=True,clinical_outcomes_used=False))
 save(R/'provenance/terminal_v1.json',dict(status='DEVELOPMENT_NATIVE_RANK_UNIVERSE_AND_DISJOINT_GENE_SENSITIVITY_COMPLETE',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(op)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise

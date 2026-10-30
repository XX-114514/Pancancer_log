"""Complete four frozen feature-coverage-gated sources with official inference.

The strict project-threshold view is retained. This branch uses the official
matched-feature inference API on the predeclared positive-RNA sources, with
unchanged model numbers, genes, normalization, scope and block size. It is
not proof of biological identity accuracy or a new source-qualification gate.
"""
import argparse,copy,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def main(R):
    cfg=load(R/'config/spec_v1.json');A=Path(cfg['automatic']);D=Path(cfg['automatic_auditor']);old=load(A/'config/spec_v1.json');H=Path(old['model']);J=Path(old['joined']);sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(A/'scripts'));sys.path.insert(0,str(R/'scripts'))
    from native_archive_collection_v1 import NativeArchiveCollection
    from native_annotation_review_v1 import native_measured_genes,predict_reference,model_fingerprint
    from native_panel_aggregate_v1 import gather_csr_rows
    from audit_TME_standard_Harmony_v2 import gather_nonnegative_float_delta
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    from CellTypist_official_unicode_axis_v1 import official_annotate_unicode_axis,independent_exact_membership_check
    import anndata,celltypist,logging
    from celltypist.models import Model
    logging.getLogger('celltypist').setLevel(logging.WARNING);storage_gate(R,cfg)
    for root in [A,D]:assert load(root/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    ap=D/'outputs/accepted_TME_CellTypist_cluster_crosscheck_v1.json';assert sha(ap)==load(D/'provenance/terminal_v1.json')['acceptance_SHA256'];accepted=load(ap);auto=load(A/'outputs/CellTypist_cluster_crosscheck_v1.json');assert sha(A/'outputs/CellTypist_cluster_crosscheck_v1.json')==accepted['automatic_crosscheck']['SHA256']
    eq=load(old['metadata_equivalence_acceptance']['path']);assert sha(R/'scripts/CellTypist_official_unicode_axis_v1.py')==eq['helper_SHA256'];assert celltypist.__version__=='1.7.1'
    with np.load(H/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as f:axes={k:f[k] for k in f.files}
    N=len(axes['source_code']);assert N==auto['n_fit'];lineage=load(H/'outputs/source_input_lineage_v1.json');totals=np.load(H/'outputs/private_analysis_RNA_totals_v1.npy',mmap_mode='r');native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest']);models={};parts={};base={};before={};probes=[];source_rows=[]
    for key in cfg['core_models']:
        m=old['models'][key];assert sha(m['path'])==m['sha256'];models[key]=Model.load(m['path']);assert models[key].cell_types.tolist()==auto['models'][key]['labels'];base[key]={};parts[key]={}
        for field in ['best_code','n_positive','max_probability']:
            path=A/'outputs'/('private_'+key+'_'+field+'_v1.npy');before[str(path)]=sha(path);v=np.load(path,mmap_mode='r');assert v.shape==(N,);base[key][field]=v;dest=np.lib.format.open_memmap(R/'outputs'/('private_'+key+'_'+field+'_v1.npy'),mode='w+',dtype=v.dtype,shape=v.shape);dest[:]=v;parts[key][field]=dest
    target=np.zeros(N,bool)
    for source in cfg['predeclared_sources']:
        code=next(i for i,r in enumerate(lineage) if r['source']==source);r=lineage[code];wh=np.flatnonzero(axes['source_code']==code);prior=A/'outputs'/('private_source_coverage_'+str(code).zfill(2)+'_v1.json');previous=load(prior);assert sha(prior)==cfg['baseline_source_coverage_SHA256'][source];assert previous['source']==source and previous['n_fit']==len(wh) and previous['n_zero_RNA']==0;assert (totals[wh]>0).all()
        for key in models:assert previous['reference_predictions'][key]==0 and previous['measured_feature_fractions'][key]<old['minimum_feature_fraction'];assert (base[key]['best_code'][wh]==-1).all()
        target[wh]=True;store=native.store(source);genes,columns=native_measured_genes(native,source);by=dict(zip(columns.tolist(),genes));all_names=[by.get(j,'__UNRESOLVED_NATIVE_COLUMN_'+str(j)) for j in range(int(store[r['counts_key']].attrs['shape'][1]))];tested=set();fingerprints={key:model_fingerprint(model) for key,model in models.items()}
        for key,model in models.items():assert independent_exact_membership_check(all_names,model.features).any()
        for lo in range(0,len(wh),old['read_block_size']):
            ix=wh[lo:lo+old['read_block_size']];x=gather_csr_rows(store[r['counts_key']],axes['native_row'][ix],max_span=8192).astype(np.float64).tocsr();np.testing.assert_array_equal(np.asarray(x.sum(axis=1)).ravel(),axes['native_total_RNA'][ix])
            if 'ambient_removed_counts' in store['layers']:x=(x-gather_nonnegative_float_delta(store['layers/ambient_removed_counts'],axes['native_row'][ix])).tocsr()
            x.eliminate_zeros();assert x.data.min(initial=0)>=-1e-8;x.data[x.data<0]=0;x.eliminate_zeros();total=np.asarray(x.sum(axis=1)).ravel();np.testing.assert_allclose(total,totals[ix],rtol=1e-10,atol=1e-7);assert (total>0).all();x=x.multiply((10000/total)[:,None]).tocsr();x.data=np.log1p(x.data);measured=x[:,columns].tocsr()
            for key,model in models.items():
                ad=anndata.AnnData(x.copy(),var=pd.DataFrame(index=all_names));ad.obs_names=ix.astype(str);pred=official_annotate_unicode_axis(ad,model,mode='best match');labels=pred.predicted_labels.predicted_labels.astype(str).to_numpy();classes=auto['models'][key]['labels'];codes=np.asarray([classes.index(v) for v in labels],np.int16);prob=pred.probability_matrix.to_numpy()
                if key not in tested:
                    q=predict_reference(measured[:16],genes,model,0.0);np.testing.assert_array_equal(codes[:16],q['best_code']);np.testing.assert_allclose(prob[:16],q['probability_scores'],rtol=1e-6,atol=1e-7)
                    ordinary=celltypist.annotate(ad.copy(),model=copy.deepcopy(model),mode='best match',majority_voting=False);np.testing.assert_array_equal(ordinary.probability_matrix.to_numpy(),prob);np.testing.assert_array_equal(ordinary.predicted_labels.predicted_labels.astype(str),labels);tested.add(key);probes.append(dict(source=source,model=key,n_numpy_probe_cells=min(16,len(ix)),n_same_batch_ordinary_probe_cells=len(ix),ordinary_official_probabilities_bitwise_equal=True,independent_numpy_replay_equal=True))
                parts[key]['best_code'][ix]=codes;parts[key]['n_positive'][ix]=(prob>.5).sum(axis=1);parts[key]['max_probability'][ix]=prob.max(axis=1)
            if lo==0 or lo//old['read_block_size']%50==0:save(R/'checkpoints'/(str(time.time_ns())+'.json'),dict(source=source,processed_source_cells=min(lo+old['read_block_size'],len(wh)),n_source_cells=len(wh)))
        for key,model in models.items():assert model_fingerprint(model)==fingerprints[key]
        for fields in parts.values():
            for a in fields.values():a.flush()
        source_rows.append(dict(source=source,n_completed=len(wh),baseline_source_coverage_SHA256=sha(prior),baseline_matching_feature_fractions=previous['measured_feature_fractions'],old_project_fraction_gate=old['minimum_feature_fraction'],official_inference_mode='best_match_matched_features_no_project_percentage_cutoff',no_zero_RNA_cells=True,biological_accuracy_not_claimed=True));native._fields.clear();storage_gate(R,cfg)
    summaries=[]
    for key in models:
        fields={}
        for name,a in parts[key].items():
            a.flush();np.testing.assert_array_equal(a[~target],base[key][name][~target]);path=R/'outputs'/('private_'+key+'_'+name+'_v1.npy');fields[name]=dict(path=str(path),SHA256=sha(path))
        code=np.load(fields['best_code']['path'],mmap_mode='r');p=np.load(fields['max_probability']['path'],mmap_mode='r');pos=np.load(fields['n_positive']['path'],mmap_mode='r');assert (code[target]>=0).all() and (code<len(auto['models'][key]['labels'])).all();valid=code>=0;assert np.isfinite(p[valid]).all() and ((p[valid]>=0)&(p[valid]<=1)).all() and ((pos[valid]>=0)&(pos[valid]<=len(auto['models'][key]['labels']))).all();summaries.append(dict(model=key,fields=fields,n_predicted=int(valid.sum()),n_added=int(target.sum()),labels=auto['models'][key]['labels'],model_SHA256=auto['models'][key]['model_SHA256']))
    assert all(sha(path)==digest for path,digest in before.items());np.save(R/'outputs/private_completed_source_cell_mask_v1.npy',target)
    result=dict(schema='accepted_TME_core_reference_coverage_completion_v1',created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_fit=N,n_added_cells=int(target.sum()),sources=source_rows,models=summaries,baseline_execution_acceptance=dict(path=str(ap),SHA256=sha(ap)),baseline_automatic=str(A),all_existing_predictions_and_model_numbers_unchanged=True,original_strict_coverage_view_retained=True,input_axis_SHA256=sha(H/'outputs/private_cell_axes_v1.npz'),independent_probes=probes,default_official_inference_not_calibrated_biological_accuracy=True,no_new_source_QC_or_patient_certification=True,final_identity_or_CM_gate_not_granted=True)
    save(R/'outputs/accepted_core_coverage_completion_v1.json',result);storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='PREDECLARED_SOURCE_CORE_REFERENCE_COVERAGE_COMPLETED_WITH_OFFICIAL_INFERENCE',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(R/'outputs/accepted_core_coverage_completion_v1.json')))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise

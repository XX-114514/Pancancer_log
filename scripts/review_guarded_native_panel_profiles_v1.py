"""Quantitative multigene review of frozen guarded subsets, without relabeling.

Reuses audited complete-native-gene sufficient statistics. Detection co-support
is described separately from RNA intensity. Group means do not establish
same-cell coexpression, contamination or a doublet diagnosis. No threshold is
learned from clinical outcomes and no inferential test is performed.
"""
import argparse,collections,datetime,gc,json,resource,sys,time
from pathlib import Path
import numpy as np
from annotation_representation_pilot_v1 import save
from audit_full_source_local_graphs_v1 import index_at
from extract_full_source_local_marker_evidence_v1 import storage_gate
from full_context_atlas_v1 import sha256_file
from native_panel_aggregate_v1 import LINEAGE_PANELS,STATE_PANELS
from native_annotation_review_v1 import TYPE_PANELS
from native_lineage_routing_v1 import FAMILIES
from summarize_full_native_annotation_evidence_v1 import read_source


def panel_profiles(n,cp,log,detected,genes,panels):
    n=np.asarray(n);genes=list(genes)
    if len(set(genes))!=len(genes) or any(v.shape!=(len(n),len(genes)) for v in [cp,log,detected]) or (n<=0).any():raise ValueError('Native statistic axes invalid')
    if any(not np.isfinite(v).all() or (v<0).any() for v in [cp,log,detected]) or np.any(detected>n[:,None]):raise ValueError('Native statistic domain invalid')
    lookup={g:j for j,g in enumerate(genes)};k=len(panels)
    out={name:np.full((len(n),k),np.nan) for name in ['mean_CP10K','mean_log1p_CP10K','mean_gene_detection_fraction']};coverage=[]
    for j,(name,panel) in enumerate(panels.items()):
        measured=[g for g in panel if g in lookup];coverage.append({'panel':name,'requested_genes':panel,'measured_genes':measured,'n_measured':len(measured),'n_requested':len(panel)})
        if not measured:continue
        q=[lookup[g] for g in measured]
        for key,values in [('mean_CP10K',cp),('mean_log1p_CP10K',log),('mean_gene_detection_fraction',detected)]:out[key][:,j]=values[:,q].sum(axis=1)/(n*len(q))
    return out,coverage


def tests():
    n=np.array([2,4]);cp=np.array([[4.,0],[0,8]]);log=cp/2;det=np.array([[2.,0],[0,2.]])
    p,c=panel_profiles(n,cp,log,det,['A','B'],{'one':['A'],'unmeasured':['X'],'partial':['A','X']})
    np.testing.assert_array_equal(p['mean_CP10K'][:,0],[2,0]);np.testing.assert_array_equal(p['mean_CP10K'][:,2],[2,0]);assert np.isnan(p['mean_CP10K'][:,1]).all() and c[2]['n_measured']==1
    # Weighted descriptive means use cells, not an unweighted mean of groups.
    assert np.average(p['mean_CP10K'][:,0],weights=n)!=np.mean(p['mean_CP10K'][:,0])
    for bad in [np.array([[3.,0],[0,2]]),np.array([[np.nan,0],[0,2]])]:
        try:panel_profiles(n,cp,log,bad,['A','B'],{'x':['A']})
        except ValueError:pass
        else:raise AssertionError('Impossible/missing sufficient statistics accepted')


def main(root):
    start=time.monotonic();cfg=json.loads((root/'config/guarded_panel_review_spec_v1.json').read_text())
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Panel review input changed')
    tests();save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':3,'no_skips':True})
    markers,cat,evidence=map(Path,[cfg['markers'],cfg['catalog'],cfg['native_evidence']])
    mi=index_at(markers,'full_source_local_marker_evidence_index_v1.json');ci=index_at(cat,'full_guarded_identity_catalog_index_v1.json');ei=index_at(evidence,'full_archive_reference_index_v1.json')
    if set(mi)!=set(ci) or set(ci)!=set(ei):raise ValueError('All-source panel review scope differs')
    panels={**{'lineage__'+k:v for k,v in LINEAGE_PANELS.items()},**{'type__'+k:v for k,v in TYPE_PANELS.items()}}
    index=[];overall=collections.Counter();type_counts=collections.Counter();type_groups=collections.Counter();weighted_profiles={}
    for ds in sorted(mi):
        storage_gate(root,cfg);p=markers/'outputs'/mi[ds]['metadata']
        if sha256_file(p)!=mi[ds]['metadata_sha256']:raise ValueError('Audited marker metadata changed')
        m=json.loads(p.read_text());sp=markers/'outputs'/m['group_native_statistics']
        if sha256_file(sp)!=m['group_native_statistics_sha256']:raise ValueError('Audited native sufficient statistics changed')
        with np.load(sp,allow_pickle=False) as f:stat={k:f[k] for k in f.files}
        genes=m['native_measured_gene_names'];profiles,coverage=panel_profiles(stat['n'],stat['cp_sum'],stat['log_sum'],stat['detected'],genes,panels)
        cm,ca=read_source(cat,ci[ds]);_,group=read_source(markers,m['source_graph_record']);em,ea=read_source(evidence,ei[ds]);fit=group['fit_pool_candidate'];code=group['source_local_group_code'];ng=len(stat['n'])
        for side in [ca,ea]:np.testing.assert_array_equal(side['native_row'],group['native_row'])
        np.testing.assert_array_equal(np.bincount(code[fit],minlength=ng),stat['n'])
        bits=ca['native_family_bitmask'];t=(bits&(1<<FAMILIES.index('T_NK')))!=0;my=(bits&(1<<FAMILIES.index('Myeloid')))!=0
        tmy=np.bincount(code[fit&t&my],minlength=ng).astype(np.int64)
        if len(m['guarded_subsets'])!=ng:raise ValueError('Guarded descriptor count differs')
        types=[r['guarded_celltype_candidate'] for r in m['guarded_subsets']]
        for j,r in enumerate(m['guarded_subsets']):
            if r['guarded_subset_code']!=j or r['n_subset_fit_rows']!=stat['n'][j]:raise ValueError('Guarded exact membership/count differs')
        state=np.full((ng,len(em['state_order'])),np.nan);state_n=np.zeros_like(state,dtype=np.int64)
        for j,name in enumerate(em['state_order']):
            v=ea['state_native_logCP10K'][:,j];ok=fit&np.isfinite(v);state_n[:,j]=np.bincount(code[ok],minlength=ng);total=np.bincount(code[ok],weights=v[ok],minlength=ng);np.divide(total,state_n[:,j],out=state[:,j],where=state_n[:,j]>0)
        output=root/'outputs'/(ds+'_guarded_native_multigene_profiles_v1.npz')
        with output.open('xb') as h:np.savez_compressed(h,n_cells=stat['n'],native_T_NK_Myeloid_binary_cosupport_count=tmy,**profiles,native_RNA_state_mean=state,native_RNA_state_measured_cells=state_n)
        descriptions=[]
        for j,t in enumerate(types):
            n=int(stat['n'][j]);type_counts[t]+=n;type_groups[t]+=1
            descriptions.append({'source_subset_code':j,'candidate_celltype':t,'n_cells':n,'native_T_NK_Myeloid_binary_cosupport_count':int(tmy[j]),'own_marker':m['groups'][j]['own_group_markers'][0]['gene'] if m['groups'][j]['own_group_markers'] else None})
            key=(t,tuple(tuple(c['measured_genes']) for c in coverage));acc=weighted_profiles.setdefault(key,{'n_cells':0,'cp_sum':np.zeros(len(panels)),'log_sum':np.zeros(len(panels)),'detection_sum':np.zeros(len(panels)),'n_groups':0,'sources':set()})
            acc['n_cells']+=n;acc['n_groups']+=1;acc['sources'].add(ds)
            for name,label in [('cp_sum','mean_CP10K'),('log_sum','mean_log1p_CP10K'),('detection_sum','mean_gene_detection_fraction')]:acc[name]+=np.nan_to_num(profiles[label][j])*n
        meta=output.with_suffix('.json');save(meta,{'source':ds,'n_groups':ng,'native_ids_sha256':cm['native_ids_sha256'],'guarded_subset_descriptors':descriptions,'panels':coverage,'RNA_state_order':em['state_order'],'state_measurement_coverage':em['state_measurement_coverage'],'arrays':output.name,'arrays_sha256':sha256_file(output),'input_marker_metadata_sha256':mi[ds]['metadata_sha256'],'input_statistics_sha256':m['group_native_statistics_sha256'],'group_means_do_not_establish_same_cell_coexpression':True,'binary_cosupport_not_doublet_or_ambient_diagnosis':True,'no_automatic_relabeling':True})
        index.append({'source':ds,'n_groups':ng,'n_fit_rows':int(stat['n'].sum()),'metadata':meta.name,'metadata_sha256':sha256_file(meta),'arrays':output.name,'arrays_sha256':sha256_file(output)});overall.update(n_groups=ng,n_fit_rows=int(stat['n'].sum()),n_T_NK_Myeloid_binary_cosupport=int(tmy.sum()))
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'n_complete_sources':len(index)})
        del stat,profiles,ca,ea,group,state,state_n;gc.collect()
    aggregate=[]
    for (t,coverage),v in sorted(weighted_profiles.items()):
        means={label:[float(x/v['n_cells']) if measured else None for x,measured in zip(v[key],coverage)] for key,label in [('cp_sum','cell_weighted_mean_CP10K'),('log_sum','cell_weighted_mean_log1p_CP10K'),('detection_sum','cell_weighted_mean_detection_fraction')]}
        aggregate.append({'candidate_celltype':t,'panel_measured_genes':[list(q) for q in coverage],'n_cells':v['n_cells'],'n_groups':v['n_groups'],'n_sources':len(v['sources']),**means,'missing_panels_remain_null':True})
    ip=root/'outputs/guarded_native_panel_review_index_v1.json';save(ip,index)
    summary=root/'outputs/guarded_native_panel_review_summary_v1.json';save(summary,{'n_sources':len(index),**dict(overall),'candidate_type_fit_cell_counts':dict(type_counts),'candidate_type_group_counts':dict(type_groups),'panel_order':list(panels),'coverage_stratified_candidate_profiles':aggregate,'clinical_labels_not_used':True,'no_identity_or_state_labels_changed':True,'not_new_independent_identity_validation':True})
    save(root/'provenance/terminal_v1.json',{'status':'ALL_GUARDED_NATIVE_MULTIGENE_INTENSITY_DETECTION_STATE_REVIEW_PROFILE_COMPLETE_NO_RELABELING_AUDIT_REQUIRED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(index),**dict(overall),'index_sha256':sha256_file(ip),'summary_sha256':sha256_file(summary),'native_counts_reread':False,'identity_labels_changed':False,'clinical_labels_used':False,'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as e:save(root/'provenance/failure_v1.json',{'error_class':type(e).__name__,'error':str(e)});raise

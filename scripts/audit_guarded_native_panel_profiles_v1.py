"""Independently replay all frozen group/panel statistics and RNA-state means."""
import argparse,collections,datetime,gc,json,resource,time
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from annotation_representation_pilot_v1 import save
from audit_full_source_local_graphs_v1 import index_at
from extract_full_source_local_marker_evidence_v1 import storage_gate
from full_context_atlas_v1 import sha256_file
from joined_native_atlas_delivery_v1 import completed
from native_panel_aggregate_v1 import LINEAGE_PANELS
from native_annotation_review_v1 import TYPE_PANELS
from native_lineage_routing_v1 import FAMILIES
from summarize_full_native_annotation_evidence_v1 import read_source

def main(root):
    start=time.monotonic();cfg=json.loads((root/'config/guarded_panel_audit_spec_v1.json').read_text());producer=Path(cfg['producer']);t=completed(producer)
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Native panel audit bound input changed')
    pc=json.loads((producer/'config/guarded_panel_review_spec_v1.json').read_text());markers,evidence=map(Path,[pc['markers'],pc['native_evidence']]);mi=index_at(markers,'full_source_local_marker_evidence_index_v1.json');ei=index_at(evidence,'full_archive_reference_index_v1.json');pi=index_at(producer,'guarded_native_panel_review_index_v1.json')
    if set(mi)!=set(ei) or set(pi)!=set(mi):raise ValueError('Full panel audit scope differs')
    panels={**{'lineage__'+k:v for k,v in LINEAGE_PANELS.items()},**{'type__'+k:v for k,v in TYPE_PANELS.items()}}
    total=collections.Counter();types=collections.Counter();type_groups=collections.Counter();aggregate={};records=[]
    for ds in sorted(pi):
        storage_gate(root,cfg);mr=mi[ds];p=markers/'outputs'/mr['metadata']
        if sha256_file(p)!=mr['metadata_sha256']:raise ValueError('Marker metadata changed')
        m=json.loads(p.read_text());p=markers/'outputs'/m['group_native_statistics']
        if sha256_file(p)!=m['group_native_statistics_sha256']:raise ValueError('Native statistic input changed')
        with np.load(p,allow_pickle=False) as f:stat={k:f[k] for k in f.files}
        p=producer/'outputs'/pi[ds]['metadata']
        if sha256_file(p)!=pi[ds]['metadata_sha256']:raise ValueError('Profile metadata changed')
        pm=json.loads(p.read_text());p=producer/'outputs'/pi[ds]['arrays']
        if sha256_file(p)!=pi[ds]['arrays_sha256']:raise ValueError('Profile arrays changed')
        with np.load(p,allow_pickle=False) as f:out={k:f[k] for k in f.files}
        genes=m['native_measured_gene_names'];lookup={g:i for i,g in enumerate(genes)};n=stat['n'];ng=len(n);np.testing.assert_array_equal(out['n_cells'],n)
        measured=[]
        for j,(name,panel) in enumerate(panels.items()):
            got=[g for g in panel if g in lookup];measured.append(tuple(got));expected={'panel':name,'requested_genes':panel,'measured_genes':got,'n_measured':len(got),'n_requested':len(panel)}
            if pm['panels'][j]!=expected:raise ValueError('Panel native measurement coverage changed')
            for k,label in [('cp_sum','mean_CP10K'),('log_sum','mean_log1p_CP10K'),('detected','mean_gene_detection_fraction')]:
                value=np.full(ng,np.nan)
                if got:value=sum(stat[k][:,lookup[g]] for g in got)/(n*len(got))
                np.testing.assert_allclose(out[label][:,j],value,atol=1e-10,rtol=1e-12,equal_nan=True)
        _,group=read_source(markers,m['source_graph_record']);em,ea=read_source(evidence,ei[ds]);fit=group['fit_pool_candidate'];code=group['source_local_group_code'];row=np.flatnonzero(fit);ng=len(n)
        incidence=sp.csr_matrix((np.ones(len(row)),(code[row],row)),shape=(ng,len(fit)))
        bit=ea['native_family_bitmask'];both=((bit&(1<<FAMILIES.index('T_NK')))!=0)&((bit&(1<<FAMILIES.index('Myeloid')))!=0)
        cosupport=np.asarray(incidence@both.astype(float)).ravel().astype(np.int64);np.testing.assert_array_equal(cosupport,out['native_T_NK_Myeloid_binary_cosupport_count'])
        if pm['RNA_state_order']!=em['state_order']:raise ValueError('State axis changed')
        for j in range(len(em['state_order'])):
            values=ea['state_native_logCP10K'][:,j];ok=np.isfinite(values);number=np.asarray(incidence@ok.astype(float)).ravel();sums=np.asarray(incidence@np.where(ok,values,0).astype(float)).ravel();expected=np.full(ng,np.nan);np.divide(sums,number,out=expected,where=number>0)
            np.testing.assert_array_equal(out['native_RNA_state_measured_cells'][:,j],number);np.testing.assert_allclose(out['native_RNA_state_mean'][:,j],expected,rtol=1e-9,atol=1e-10,equal_nan=True)
        for j,desc in enumerate(m['guarded_subsets']):
            kind=desc['guarded_celltype_candidate'];number=int(n[j]);types[kind]+=number;type_groups[kind]+=1
            if pm['guarded_subset_descriptors'][j]['candidate_celltype']!=kind or pm['guarded_subset_descriptors'][j]['n_cells']!=number:raise ValueError('Profile candidate membership changed')
            key=(kind,tuple(measured));v=aggregate.setdefault(key,{'n_cells':0,'n_groups':0,'sources':set(),'cp':np.zeros(len(panels)),'log':np.zeros(len(panels)),'detection':np.zeros(len(panels))});v['n_cells']+=number;v['n_groups']+=1;v['sources'].add(ds)
            for k,label in [('cp_sum','cp'),('log_sum','log'),('detected','detection')]:
                for col,mg in enumerate(measured):
                    if mg:v[label][col]+=sum(float(stat[k][j,lookup[g]]) for g in mg)/len(mg)
        records.append({'source':ds,'n_groups':ng,'n_fit_rows':int(n.sum()),'n_T_NK_Myeloid_binary_cosupport':int(cosupport.sum())});total.update({k:v for k,v in records[-1].items() if k!='source'});save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'n_complete_sources':len(records)})
        del stat,group,ea,out,incidence;gc.collect()
    p=producer/'outputs/guarded_native_panel_review_summary_v1.json'
    if sha256_file(p)!=t['summary_sha256']:raise ValueError('Panel review summary changed')
    summary=json.loads(p.read_text())
    if summary['candidate_type_fit_cell_counts']!=dict(types) or summary['candidate_type_group_counts']!=dict(type_groups):raise ValueError('Candidate-type denominators differ')
    keys=set()
    for a in summary['coverage_stratified_candidate_profiles']:
        key=(a['candidate_celltype'],tuple(tuple(g) for g in a['panel_measured_genes']));v=aggregate[key];keys.add(key)
        if a['n_cells']!=v['n_cells'] or a['n_groups']!=v['n_groups'] or a['n_sources']!=len(v['sources']):raise ValueError('Coverage-stratified denominator differs')
        for k,label in [('cp','cell_weighted_mean_CP10K'),('log','cell_weighted_mean_log1p_CP10K'),('detection','cell_weighted_mean_detection_fraction')]:
            for i,genes in enumerate(key[1]):
                if genes:np.testing.assert_allclose(a[label][i],v[k][i]/v['n_cells'],rtol=1e-10,atol=1e-10)
                elif a[label][i] is not None:raise ValueError('Unmeasured panel filled with biological zero')
    if keys!=set(aggregate) or total['n_groups']!=t['n_groups'] or total['n_fit_rows']!=t['n_fit_rows']:raise ValueError('Complete profile/aggregate scope differs')
    out=root/'outputs/guarded_native_panel_review_audit_v1.json';save(out,{'sources':records,'totals':dict(total),'all_panels_all_groups_all_state_denominators_and_coverage_stratified_means_replayed':True,'biological_identity_not_validated_by_reusing_annotation_markers':True})
    save(root/'provenance/terminal_v1.json',{'status':'ALL_GUARDED_GROUP_MULTIGENE_PANEL_STATE_AND_AGGREGATE_PROFILES_AUDITED_NO_RELABELING','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(records),**dict(total),'output_sha256':sha256_file(out),'producer_terminal_sha256':sha256_file(producer/'provenance/terminal_v1.json'),'identity_labels_changed':False,'new_independent_biological_validation':False,'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as e:save(root/'provenance/failure_v1.json',{'error_class':type(e).__name__,'error':str(e)});raise

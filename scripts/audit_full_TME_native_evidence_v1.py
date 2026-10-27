"""Independent complete native replay, consuming completed source checkpoints.

No new production normalizer, marker selector or CSR gather helper is imported.
The archive reader supplies only the frozen storage/identity/gene-lookup contract.
Acceptance requires successful real production termination and the full scope.
"""
import argparse,collections,datetime,gc,hashlib,json,re,resource,sys,time
from pathlib import Path
import numpy as np
import scipy.sparse as sp

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def live(root):
    p=root/'provenance/child_launch_v1.json'
    if not p.exists():return False
    d=load(p);p=Path('/proc')/str(d['PID'])/'stat'
    try:s=p.read_text().rsplit(')',1)[1].split();return s[19]==str(d['startticks']) and s[0]!='Z'
    except FileNotFoundError:return False
def gather(csr,rows):
    unique,inverse=np.unique(rows,return_inverse=True);parts=[];low=0
    while low<len(unique):
        high=int(np.searchsorted(unique,unique[low]+2048));first=int(unique[low]);last=int(unique[high-1])+1
        pointers=np.asarray(csr['indptr'][first:last+1],np.int64);start=int(pointers[0]);end=int(pointers[-1]);data=np.asarray(csr['data'][start:end],np.float64);indices=np.asarray(csr['indices'][start:end],np.int64)
        if not np.isfinite(data).all() or (data<0).any() or not np.equal(data,np.floor(data)).all():raise ValueError('Invalid original native values')
        block=sp.csr_matrix((data,indices,pointers-start),shape=(last-first,int(csr.attrs['shape'][1])));block.sum_duplicates();block.eliminate_zeros()
        parts.append(block[unique[low:high]-first]);low=high
    return sp.vstack(parts,format='csr')[inverse]
def family(major):
    return 'T_NK' if major in ['T','NK','T_NK'] else 'B_Plasma' if major in ['B','Plasma','B_Plasma'] else major

def audit_source(root,P,reader,source_item,source_groups,cfg):
    from anndata._io.specs import read_elem
    ds=source_item['source'];m=load(P/'outputs'/(ds+'_TME_native_evidence_v1.json'));assert m['source']==ds
    workbase=Path(cfg['working_index']).parent;workmeta=load(workbase/source_item['metadata']);assert sha(workbase/source_item['metadata'])==source_item['metadata_SHA256'];assert sha(workbase/source_item['arrays'])==source_item['arrays_SHA256']
    with np.load(workbase/source_item['arrays'],allow_pickle=False) as f:q={k:f[k] for k in ['source_native_position','native_row','annotation_group_code','native_total_RNA','fit_pool_candidate_not_final_QC']}
    np.testing.assert_array_equal(reader.rows(ds)[q['source_native_position']],q['native_row']);n=len(q['native_row']);ng=len(source_groups)
    z=reader.store(ds);var=read_elem(z['var']);lookup=reader.feature_lookup(ds);genes=[];columns=[]
    for j,index in enumerate(var.index):
        name=str(var.iloc[j]['gene_symbol']) if 'gene_symbol' in var.columns else str(index)
        if lookup.get(name,set())=={j}:genes.append(name);columns.append(j)
    assert genes==m['native_measured_gene_names'] and columns==m['native_columns'];assert len(set(genes))==len(genes)
    expected={'n':np.zeros(ng,np.int64),**{k:np.zeros((ng,len(genes)),np.float64) for k in ['cp_sum','log_sum','detected']}}
    contextp=Path(m['native_context']['path']);statp=Path(m['group_statistics']['path']);assert sha(contextp)==m['native_context']['SHA256'] and sha(statp)==m['group_statistics']['SHA256']
    with np.load(contextp,allow_pickle=False) as f:c={k:f[k] for k in f.files}
    for name in ['source_native_position','native_row','annotation_group_code','native_total_RNA']:np.testing.assert_array_equal(c[name],q[name])
    by={g:col for g,col in zip(genes,columns)};physical=[by.get(g) for g in m['context_genes']];mask=np.array([v is not None for v in physical]);np.testing.assert_array_equal(mask,c['measured']);assert mask.tolist()==m['context_gene_measured']
    assert np.isnan(c['native_gene_values'][:,~mask]).all()
    key=reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key'];csr=z[key];canonical_nnz=0
    for low in range(0,n,1024):
        high=min(n,low+1024);raw=gather(csr,q['native_row'][low:high]);canonical_nnz+=raw.nnz;totals=np.asarray(raw.sum(axis=1,dtype=np.float64)).ravel();np.testing.assert_array_equal(totals,q['native_total_RNA'][low:high])
        if mask.any():np.testing.assert_array_equal(raw[:,[v for v in physical if v is not None]].toarray(),c['native_gene_values'][low:high,mask])
        codes=q['annotation_group_code'][low:high];fit=codes>=0;assert (totals[fit]>0).all()
        if not fit.any():continue
        measured=raw[fit][:,columns].tocsr();which=codes[fit];expected['n']+=np.bincount(which,minlength=ng)
        coo=measured.tocoo();values=coo.data*(10000/totals[fit])[coo.row]
        norm=sp.csr_matrix((values,(coo.row,coo.col)),shape=measured.shape);logged=norm.copy();logged.data=np.log1p(logged.data)
        link=sp.csr_matrix((np.ones(len(which)),(which,np.arange(len(which)))),shape=(ng,len(which)))
        for name,x in [('cp_sum',norm),('log_sum',logged),('detected',measured>0)]:
            v=(link@x).tocoo();np.add.at(expected[name],(v.row,v.col),v.data)
    assert canonical_nnz==m['canonical_native_nnz']
    with np.load(statp,allow_pickle=False) as f:saved={k:f[k] for k in f.files}
    np.testing.assert_array_equal(saved['native_columns'],columns)
    for name in ['n','detected']:np.testing.assert_array_equal(expected[name],saved[name])
    for name in ['cp_sum','log_sum']:np.testing.assert_allclose(expected[name],saved[name],rtol=2e-10,atol=1e-7)
    families=[]
    for g in source_groups:
        names=list(g['major_distribution']);families.append(family(names[0]) if len(names)==1 else 'mixed_major_review')
    families=np.array(families,object);assert len(m['groups'])==ng
    for j,(g,record) in enumerate(zip(source_groups,m['groups'])):
        assert record['source_annotation_group_code']==j and record['stable_group_sha256']==g['stable_group_sha256'];assert record['n_exact_group_cells']==int(saved['n'][j]) and record['current_analysis_family']==families[j]
        compare=families==families[j];compare[j]=False;count=int(saved['n'][compare].sum());scope='other_current_TME_source_groups_in_same_analysis_family'
        if not count:compare=np.arange(ng)!=j;count=int(saved['n'][compare].sum());scope='other_current_TME_source_groups_no_within_family_comparator'
        assert count==record['n_comparator_cells'] and scope==record['comparator_scope']
        cp=saved['cp_sum'][j]/saved['n'][j];det=saved['detected'][j]/saved['n'][j]
        eligible=np.array([not re.match(r'^(MT-|RPL|RPS|ENSG|IG[HKL][VDJC])',gene) for gene in genes])&(det>=cfg['minimum_detection'])&(cp>0)
        other=saved['cp_sum'][compare].sum(axis=0)/count if count else None;otherdet=saved['detected'][compare].sum(axis=0)/count if count else None
        effect=np.log2((cp+.1)/(other+.1)) if count else None
        if count:
            good=np.flatnonzero(eligible&(det>otherdet)&(effect>0));score=effect*np.maximum(det-otherdet,0);order=sorted(good,key=lambda k:(-score[k],genes[k]))[:cfg['marker_limit']];kind='own_subset_expression_specificity_not_inferential_DEG'
        else:order=[];kind='own_subset_expression_without_comparator_not_inferential_DEG'
        if not order:
            order=sorted(np.flatnonzero(eligible),key=lambda k:(-cp[k]*det[k],genes[k]))[:cfg['marker_limit']]
            if count:kind='own_subset_expression_without_positive_specificity_not_inferential_DEG'
        candidates=record['own_native_expression_marker_candidates'];assert [v['gene'] for v in candidates]==[genes[k] for k in order]
        for k,v in zip(order,candidates):
            assert v['marker_evidence_type']==kind and v['native_mean_CP10K']==float(cp[k]) and v['native_detection_fraction']==float(det[k]);assert v['comparator_detection_fraction']==(float(otherdet[k]) if count else None);assert v['descriptive_log2_mean_ratio']==(float(effect[k]) if count else None)
    out={'source':ds,'n_scope_rows_checked':n,'n_exact_groups_checked':ng,'n_context_values_checked':int(c['native_gene_values'].size),'n_measured_native_genes':len(genes),'metadata_SHA256':sha(P/'outputs'/(ds+'_TME_native_evidence_v1.json')),'all_native_context_and_full_gene_group_statistics_and_marker_ranking_checked':True}
    save(root/'outputs'/(ds+'_independent_native_evidence_audit_v1.json'),out);reader._fields.clear();gc.collect();return out

def main(root):
    start=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['producer']);pc=load(P/'config/spec_v1.json');sys.path.insert(0,str(Path(pc['joined'])/'scripts'))
    from native_archive_collection_v1 import NativeArchiveCollection
    accepted=load(pc['working_acceptance']);manifest=load(accepted['producer_manifest']['path']);index=load(manifest['index']['path']);cfg.update(working_index=manifest['index']['path'],minimum_detection=pc['minimum_detection'],marker_limit=pc['marker_limit'])
    allgroups=load(manifest['group_catalog']['path'])['groups'];expected={r['source']:r for r in index};groups={s:sorted([g for g in allgroups if g['source']==s],key=lambda g:g['source_annotation_group_code']) for s in expected}
    reader=NativeArchiveCollection(load(Path(pc['joined'])/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest']);done={}
    while len(done)<len(expected):
        ready=[s for s in expected if (P/'checkpoints'/('source_'+s+'.json')).exists() and s not in done]
        for source in ready:
            done[source]=audit_source(root,P,reader,expected[source],groups[source],cfg)
            save(root/'checkpoints'/('source_'+source+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources_checked':len(done),'n_scope_rows_checked':sum(v['n_scope_rows_checked'] for v in done.values())})
        if len(done)<len(expected):
            terminal=P/'provenance/supervisor_terminal_v1.json'
            if terminal.exists():
                assert load(terminal)['actual_child_exit_code']==0,'Original producer actually failed; partial audits preserved'
                assert ready,'Producer ended without all expected completed sources'
            elif not live(P):
                time.sleep(2)
                assert terminal.exists() and load(terminal)['actual_child_exit_code']==0,'Producer handle disappeared without successful termination'
            if time.monotonic()-start>cfg['maximum_producer_wait_seconds']:raise TimeoutError('Bounded native audit time exceeded')
            time.sleep(15)
    while not (P/'provenance/supervisor_terminal_v1.json').exists():
        if not live(P):
            time.sleep(2);assert (P/'provenance/supervisor_terminal_v1.json').exists()
        time.sleep(2)
    assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    terminal=load(P/'provenance/terminal_v1.json');pi=P/'outputs/TME_full_native_evidence_index_v1.json';assert sha(pi)==terminal['index_SHA256']
    for entry in load(pi):assert done[entry['source']]['metadata_SHA256']==entry['metadata_SHA256']
    totals={'n_sources':len(done),'n_scope_rows_checked':sum(v['n_scope_rows_checked'] for v in done.values()),'n_exact_groups_checked':sum(v['n_exact_groups_checked'] for v in done.values()),'n_context_values_checked':sum(v['n_context_values_checked'] for v in done.values())}
    assert totals['n_scope_rows_checked']==terminal['totals']['n_current_scope_records'] and totals['n_exact_groups_checked']==terminal['totals']['n_exact_groups']
    out={'status':'COMPLETE_CURRENT_TME_NATIVE_CONTEXT_AND_ALL_GENE_GROUP_STATISTICS_AND_MARKER_RANKING_INDEPENDENTLY_ACCEPTED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':totals,'producer_terminal_SHA256':sha(P/'provenance/terminal_v1.json'),'production_normalizer_marker_selector_or_CSR_gather_imported':False,'numerical_acceptance_not_biological_identity_validation':True,'no_labels_QC_or_CM_changed':True,'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024}
    op=root/'outputs/full_TME_native_evidence_audit_v1.json';save(op,out);save(root/'provenance/terminal_v1.json',dict(out,output_SHA256=sha(op)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except BaseException:
        import traceback
        save(root/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise

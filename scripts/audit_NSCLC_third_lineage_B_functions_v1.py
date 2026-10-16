"""Independent sparse pooling and pinned official ULM for every B-lineage readout."""
import ast,hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.sparse import csr_matrix


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=r/'NSCLC_third_lineage_B_functions_20261016_v1';job=r/'HNSCC_receiver_followup_job_control_20261016_v2/measure_NSCLC_third_lineage_B_functions_v1'
    receipt=json.loads((job/'launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    t=json.loads((out/'private_terminal.json').read_text());plan=json.loads((out/'plan_before_any_B_function_measurement.json').read_text());assert plan['code_sha256']==receipt['source_sha256']
    for n,h in t['files'].items():assert sha(out/n)==h
    native=r/'NSCLC_full_native_gene_units_20261016_v1';pre=r/'NSCLC_disjoint_hallmark_preflight_20261016_v1'
    paths=dict(question=r/'NSCLC_third_lineage_background_plan_20261016_v1/PRIVATE_NEXT_QUESTION_AND_SCOPE.json',units=native/'private_all_native_units_and_quality.tsv',matrix=native/'mean_logCP10K.npy',QC=r/'NSCLC_author_RNA_composition_20261016_v1/private_all_native_unit_RNA_composition.tsv',masks=pre/'private_physical_column_masks.json',references=pre/'private_reference_axis.tsv',coverage=pre/'private_all1200_functional_scope.tsv',official=r/'fullgene_receiver_target_code_v1/collectri_reference_v1/official_decoupler192_ULM_source.txt')
    for k,p in paths.items():assert sha(p)==plan['input_sha256'][k]
    u=pd.read_csv(paths['units'],sep='\t').sort_values('unit_index');b=u[u.major_cell_type.eq('B cell')].sort_values('unit_index');axis=pd.read_csv(out/'private_all_B_case_axis_and_quality.tsv',sep='\t');mapping=pd.read_csv(out/'private_all_B_native_unit_pool_mapping.tsv',sep='\t');q=pd.read_csv(paths['QC'],sep='\t').set_index('unit_index')
    cases=sorted(b.sampleID.unique());assert axis.sampleID.tolist()==cases and axis.pool_index.tolist()==list(range(len(cases))) and len(axis)==t['author_B_cases']
    totals=u.groupby('sampleID').n_cells.sum();n=b.groupby('sampleID').n_cells.sum().reindex(cases);assert np.array_equal(axis.n_cells,n)
    index=pd.Index(cases).get_indexer(b.sampleID);weights=b.n_cells.to_numpy(float)/n.to_numpy()[index];P=csr_matrix((weights,(index,np.arange(len(b)))),shape=(len(cases),len(b)));assert np.allclose(np.asarray(P.sum(axis=1)).ravel(),1,atol=1e-14)
    expected_mapping=pd.DataFrame(dict(pool_index=index,unit_index=b.unit_index.to_numpy(int),cell_weight=weights)).sort_values(['pool_index','unit_index']).reset_index(drop=True)
    pd.testing.assert_frame_equal(mapping.sort_values(['pool_index','unit_index']).reset_index(drop=True),expected_mapping,check_exact=False,rtol=1e-12,atol=1e-12)
    assert np.allclose(axis.fraction,n.to_numpy()/totals.reindex(cases).to_numpy()) and np.allclose(axis.logN,np.log1p(n.to_numpy()))
    assert np.array_equal(axis.author_subtypes,b.groupby('sampleID').sub_cell_type.nunique().reindex(cases))
    for field in ['mean_logRNA','mean_loggenes','mean_author_pct_mito','mean_author_pct_ribo']:
        values=b[field].to_numpy(float) if field in b else q.reindex(b.unit_index)[field].to_numpy(float)
        expected=np.asarray(P@values).ravel();assert np.allclose(axis[field],expected,rtol=1e-12,atol=1e-12)
    m=np.load(paths['matrix'],mmap_mode='r');pooled=np.empty((len(cases),m.shape[1]),float);units=b.unit_index.to_numpy(int)
    for start in range(0,m.shape[1],512):pooled[:,start:start+512]=P@np.asarray(m[np.ix_(units,np.arange(start,min(start+512,m.shape[1])))])
    ns={'np':np};funcs=[n for n in ast.parse(paths['official'].read_text()).body if isinstance(n,ast.FunctionDef) and n.name in {'mat_cov','mat_cor','t_val'}];assert len(funcs)==3
    exec(compile(ast.Module(body=funcs,type_ignores=[]),str(paths['official']),'exec'),ns)
    masks=json.loads(paths['masks'].read_text());refs=pd.read_csv(paths['references'],sep='\t').sort_values('reference_column').reference.tolist();mask_axis=pd.read_csv(out/'private_all16_mask_axis.tsv',sep='\t');assert mask_axis.mask_id.tolist()==sorted(masks)
    npz=np.load(out/'private_all_B_function_scores_and_QA.npz',allow_pickle=False);arrays={k:npz[k] for k in ['ULM_score','target_mean','nontarget_mean','expressed_target_genes']};assert all(a.shape==(16,len(cases),50) for a in arrays.values());maximum=0.;checked=0
    for row in mask_axis.itertuples():
        mask=masks[row.mask_id];excluded=set(mask['forbidden_native_columns']);bg=np.array([i for i in range(m.shape[1]) if i not in excluded],int);position={int(c):i for i,c in enumerate(bg)};w=np.zeros((len(bg),50),float);ix=[]
        for j,ref in enumerate(refs):
            targets=np.array([position[int(c)] for c in mask['targets'][ref]],int);w[targets,j]=1;ix.append(targets)
        e=pooled[:,bg].T
        with np.errstate(divide='ignore',invalid='ignore'):expected=ns['t_val'](ns['mat_cor'](w,e),len(bg)-2)
        observed=arrays['ULM_score'][row.mask_index];assert np.allclose(observed,expected,atol=1e-8,rtol=1e-8,equal_nan=True);maximum=max(maximum,float(np.nanmax(np.abs(observed-expected))))
        total=e.sum(axis=0)
        for j,targets in enumerate(ix):
            x=e[targets];assert np.allclose(arrays['target_mean'][row.mask_index,:,j],x.mean(axis=0),atol=1e-12,rtol=1e-10)
            assert np.allclose(arrays['nontarget_mean'][row.mask_index,:,j],(total-x.sum(axis=0))/(len(bg)-len(targets)),atol=1e-12,rtol=1e-10)
            assert np.array_equal(arrays['expressed_target_genes'][row.mask_index,:,j],np.count_nonzero(x>0,axis=0))
        checked+=expected.size
        (out/'private_independent_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),official_function_scores_checked=checked))+'\n')
    assert checked==t['functional_score_values'] and set(mapping.unit_index)==set(b.unit_index)
    proof=dict(status='PASS_ALL_AUTHOR_B_NATIVE_UNIT_POOLS_QUALITY_KEYS_AND_ALL_OFFICIAL_MASKED_FUNCTION_SCORES_QA',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),pooled_cases=len(cases),official_scores_checked=checked,maximum_score_difference=maximum,raw_replay_scope='All new pool/score outputs independently reconstructed from previously audited native-unit cache, not all original cells independently reread.',new_P_BH=False)
    p=out/'private_independent_B_function_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()

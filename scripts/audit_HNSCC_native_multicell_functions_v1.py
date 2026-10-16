"""Replay every selected raw cell and official function score independently."""
import ast,hashlib,json,os,time
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.sparse import csr_matrix
from HNSCC_native_multicell_function_inputs_v1 import inputs


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=r/'HNSCC_native_multicell_functions_20261016_v1';job=r/'HNSCC_receiver_followup_job_control_20261016_v2/measure_HNSCC_native_multicell_functions_v1';rc=json.loads((job/'launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(rc['pid'])/'stat';assert p.exists();a=p.read_text().rsplit(')',1)[1].split();assert a[0]!='Z' and a[19]==str(rc['proc_start_ticks']) and time.monotonic()<deadline;time.sleep(5)
    terminal=json.loads((out/'private_terminal.json').read_text());plan=json.loads((out/'plan_before_any_multicell_function_measurement.json').read_text());assert rc['source_sha256']==terminal['code_sha256']==plan['code_sha256']
    for name,digest in terminal['files'].items():assert sha(out/name)==digest
    data=inputs(r);assert data['input_sha256']==plan['input_sha256'] and data['coverage_input_sha256']==plan['coverage_input_sha256']
    assert sha(Path(__file__).with_name('HNSCC_native_multicell_function_inputs_v1.py'))==plan['input_helper_sha256']
    axis=pd.read_csv(out/'private_all228_pool_library_axis_and_quality.tsv',sep='\t',dtype={'library_id':str,'patient_block':str},na_values=['NA']);wanted=data['axis'];pd.testing.assert_frame_equal(axis[wanted.columns],wanted,check_dtype=False,rtol=1e-12,atol=1e-12)
    mapping=data['mapping'];counts=data['zarr']['layers']['counts'];N,G=tuple(counts.attrs['shape']);U=len(axis);means=np.zeros((U,G));quality=np.zeros((U,7));symbols=data['symbols'];mito=np.array([str(g).startswith('MT-') for g in symbols]);ribo=np.array([str(g).startswith(('RPS','RPL')) for g in symbols]);visited=0
    # Reverse, unequal blocks, no producer sparse group projection.
    for start in reversed(range(0,N,3072)):
        stop=min(N,start+3072);keep=mapping[start:stop]>=0
        if not keep.any():continue
        pointers=np.asarray(counts['indptr'][start:stop+1],np.int64);matrix=csr_matrix((np.asarray(counts['data'][int(pointers[0]):int(pointers[-1])]),np.asarray(counts['indices'][int(pointers[0]):int(pointers[-1])]),pointers-pointers[0]),shape=(stop-start,G))[keep].copy();ids=mapping[start:stop][keep]
        matrix.sum_duplicates();matrix.eliminate_zeros();assert np.isfinite(matrix.data).all() and (matrix.data>0).all() and np.equal(matrix.data,np.floor(matrix.data)).all() and data['native_mask'][matrix.indices].all()
        total=np.asarray(matrix.sum(axis=1)).ravel().astype(float);detected=np.diff(matrix.indptr);logged=matrix.astype(float);logged.data=np.log1p(logged.data*np.repeat(10000/total,detected));m=100*np.asarray(matrix[:,mito].sum(axis=1)).ravel()/total;b=100*np.asarray(matrix[:,ribo].sum(axis=1)).ravel()/total
        for index in np.unique(ids):
            selected=ids==index;means[index]+=np.asarray(logged[selected].sum(axis=0)).ravel();quality[index]+=np.array([selected.sum(),total[selected].sum(),np.log1p(total[selected]).sum(),np.log1p(detected[selected]).sum(),m[selected].sum(),b[selected].sum(),detected[selected].sum()])
        visited+=int(keep.sum())
        (out/'private_independent_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),raw_selected_cells_replayed=visited))+'\n')
    assert visited==terminal['selected_cells'] and np.array_equal(quality[:,0],axis.n_cells);valid=np.flatnonzero(axis.n_cells.to_numpy()>0);means[valid]/=quality[valid,:1]
    for j,key in enumerate(['mean_full_RNA','mean_logRNA','mean_loggenes','mean_symbol_pct_mito','mean_symbol_pct_ribo','mean_detected_genes'],start=1):
        expected=np.divide(quality[:,j],quality[:,0],out=np.full(U,np.nan),where=quality[:,0]>0);assert np.allclose(axis[key],expected,rtol=2e-11,atol=1e-11,equal_nan=True)
    masks=pd.read_csv(out/'private_all16_mask_axis.tsv',sep='\t');pd.testing.assert_frame_equal(masks,data['masks'].sort_values('mask_id').reset_index(drop=True),check_dtype=False);refs=pd.read_csv(out/'private_all50_reference_axis.tsv',sep='\t');assert refs.reference.tolist()==sorted(data['references'])
    saved=json.loads((out/'private_physical_function_masks.json').read_text());values=np.load(out/'private_all_multicell_function_scores_and_QA.npz',allow_pickle=False);assert all(v.shape==(16,228,50) for v in values.values())
    namespace={'np':np};tree=ast.parse(data['official_path'].read_text());funcs=[f for f in tree.body if isinstance(f,ast.FunctionDef) and f.name in {'mat_cov','mat_cor','t_val'}];assert len(funcs)==3;exec(compile(ast.Module(body=funcs,type_ignores=[]),str(data['official_path']),'exec'),namespace)
    physical=data['physical'];params=data['parameters'];maximum=0.;checked=0
    for i,mask in enumerate(masks.itertuples()):
        learned=set(params.loc[(params.context==mask.context)&(params.program==mask.program),'gene']);excluded=set()
        for gene in learned:excluded.update(physical['lookup'].get(gene,set()))
        background=np.array([c for c in sorted(physical['native']) if c not in excluded],int);assert background.tolist()==saved[str(mask.mask_id)]['background'];position={c:j for j,c in enumerate(background)};weights=np.zeros((len(background),50));targets=[]
        for j,reference in enumerate(refs.reference):
            owners=defaultdict(list)
            for gene in data['references'][reference]-learned:
                hits=physical['lookup'].get(gene,set())
                if len(hits)==1:
                    column,=hits
                    if column not in excluded:owners[column].append(gene)
            columns=sorted(c for c,g in owners.items() if len(g)==1);assert columns==saved[str(mask.mask_id)]['targets'][reference];ix=np.array([position[c] for c in columns],int);weights[ix,j]=1;targets.append(ix)
        expression=means[np.ix_(valid,background)].T
        with np.errstate(divide='ignore',invalid='ignore'):expected=namespace['t_val'](namespace['mat_cor'](weights,expression),len(background)-2)
        observed=values['ULM_score'][i,valid,:];assert np.allclose(observed,expected,rtol=1e-8,atol=1e-8,equal_nan=True);maximum=max(maximum,float(np.nanmax(abs(observed-expected))));sums=expression.sum(axis=0)
        for j,ix in enumerate(targets):
            foreground=expression[ix];assert np.allclose(values['target_mean'][i,valid,j],foreground.mean(axis=0),rtol=1e-10,atol=1e-12)
            assert np.allclose(values['nontarget_mean'][i,valid,j],(sums-foreground.sum(axis=0))/(len(background)-len(ix)),rtol=1e-10,atol=1e-12)
            assert np.array_equal(values['expressed_target_genes'][i,valid,j],np.count_nonzero(foreground>0,axis=0))
        absent=np.flatnonzero(axis.n_cells.to_numpy()==0)
        for key in ['ULM_score','target_mean','nontarget_mean']:assert np.isnan(values[key][i,absent,:]).all()
        assert (values['expressed_target_genes'][i,absent,:]==-1).all();checked+=expected.size
    proof=dict(status='PASS_ALL_SELECTED_RAW_CELLS_ALL228_POOL_LIBRARY_QC_ALL16_MASKS50_OFFICIAL_FUNCTIONS_AND_QA',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),selected_raw_cells_independently_replayed=visited,nonempty_pool_library_slots=len(valid),official_function_scores_checked=checked,maximum_score_difference=maximum,raw_replay='Every selected original atlas count row and all native features; reverse unequal blocks and per-pool sparse sums, independent official ULM and physical reverse mapping. No full expression copy.',new_P_BH=False)
    p=out/'private_independent_multicell_function_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()

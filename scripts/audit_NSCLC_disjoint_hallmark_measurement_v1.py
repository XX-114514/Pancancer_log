"""Every native functional score versus pinned dense official ULM and RNA QA."""
import ast,hashlib,json,os,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);pre=r/'NSCLC_disjoint_hallmark_preflight_20261016_v1';out=r/'NSCLC_disjoint_hallmark_measurement_20261016_v1'
    receipt=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2/measure_NSCLC_disjoint_hallmark_v1/launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    terminal=json.loads((out/'private_terminal.json').read_text())
    for n,h in terminal['files'].items():assert sha(out/n)==h
    plan=json.loads((pre/'plan_before_functional_RNA_measurements.json').read_text());mp=json.loads((out/'plan_before_real_functional_measurements.json').read_text())
    assert sha(pre/'plan_before_functional_RNA_measurements.json')==mp['original_plan_sha256']
    assert sha(pre/'private_independent_scope_audit.json')==mp['scope_audit_sha256']
    p=r/'fullgene_receiver_target_code_v1/collectri_reference_v1/official_decoupler192_ULM_source.txt';assert sha(p)==plan['input_sha256']['official_ULM']
    ns={'np':np};funcs=[n for n in ast.parse(p.read_text()).body if isinstance(n,ast.FunctionDef) and n.name in {'mat_cov','mat_cor','t_val'}];assert len(funcs)==3
    exec(compile(ast.Module(body=funcs,type_ignores=[]),str(p),'exec'),ns)
    axis=pd.read_csv(pre/'private_score_row_axis.tsv.gz',sep='\t');scope=pd.read_csv(pre/'private_all1200_functional_scope.tsv',sep='\t');refs=pd.read_csv(pre/'private_reference_axis.tsv',sep='\t').reference.tolist();masks=json.loads((pre/'private_physical_column_masks.json').read_text())
    p=r/'NSCLC_full_native_gene_units_20261016_v1/mean_logCP10K.npy';assert sha(p)==plan['input_sha256']['matrix'];m=np.load(p,mmap_mode='r')
    path=out/'private_all_native_functional_scores_and_QA.npz';z=np.load(path,allow_pickle=False);values={k:z[k] for k in ['ULM_score','target_mean','nontarget_mean','expressed_target_genes']}
    assert all(x.shape==(len(axis),50) for x in values.values());maximum=0.;rows_checked=0;score_values=0
    for (context,program),f in axis.groupby(['context','program'],sort=True):
        mask=masks[f'C{context}_{program}'];excluded=set(mask['forbidden_native_columns']);bg=np.array([i for i in range(m.shape[1]) if i not in excluded],int)
        position={int(g):i for i,g in enumerate(bg)};s=scope[scope.context.eq(context)&scope.program.eq(program)].drop_duplicates('reference').set_index('reference')
        names=[name for name in refs if s.loc[name,'status']=='MEASURABLE_INPUT_NO_RNA_SCORE'];columns=[refs.index(name) for name in names];missing=sorted(set(range(50))-set(columns));dest=f.score_row.to_numpy(int)
        for key in ['ULM_score','target_mean','nontarget_mean']:assert np.isnan(values[key][np.ix_(dest,missing)]).all()
        assert (values['expressed_target_genes'][np.ix_(dest,missing)]==-1).all()
        w=np.zeros((len(bg),len(names)),float);target_positions=[]
        for j,name in enumerate(names):
            pos=np.array([position[int(g)] for g in mask['targets'][name]],int);target_positions.append(pos);w[pos,j]=1
        if not names:rows_checked+=len(f);continue
        for start in range(0,len(f),128):
            part=f.iloc[start:start+128];ids=part.score_row.to_numpy(int);e=np.asarray(m[np.ix_(part.unit_index.to_numpy(int),bg)]).T
            with np.errstate(divide='ignore',invalid='ignore'):expected=ns['t_val'](ns['mat_cor'](w,e),len(bg)-2)
            actual=values['ULM_score'][np.ix_(ids,columns)]
            assert np.array_equal(np.isfinite(actual),np.isfinite(expected)) and np.allclose(actual,expected,atol=1e-8,rtol=1e-8,equal_nan=True)
            if np.isfinite(expected).any():maximum=max(maximum,float(np.nanmax(np.abs(actual-expected))))
            total=e.sum(axis=0)
            for j,pos in enumerate(target_positions):
                x=e[pos];mean=x.mean(axis=0);other=(total-x.sum(axis=0))/(len(bg)-len(pos));detected=np.count_nonzero(x>0,axis=0)
                assert np.allclose(values['target_mean'][ids,columns[j]],mean,atol=1e-12,rtol=1e-10)
                assert np.allclose(values['nontarget_mean'][ids,columns[j]],other,atol=1e-12,rtol=1e-10)
                assert np.array_equal(values['expressed_target_genes'][ids,columns[j]],detected)
                if start==0:
                    keep=np.ones(len(bg),bool);keep[pos]=False
                    assert np.isclose(e[keep,0].mean(),other[0],atol=1e-12)
            rows_checked+=len(part);score_values+=expected.size
            (out/'private_independent_audit_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),native_rows_checked=rows_checked,official_score_values_checked=score_values))+'\n')
    summary=pd.read_csv(out/'private_all1200_measurement_summary.tsv',sep='\t');assert len(summary)==1200
    for row in summary.itertuples():
        cov=scope[scope.readout_id.eq(row.readout_id)].iloc[0];ids=axis.loc[axis.context.eq(cov.context)&axis.type.eq(cov.type)&axis.program.eq(cov.program),'score_row'].to_numpy(int);j=int(cov.reference_column)
        x=values['ULM_score'][ids,j];finite=x[np.isfinite(x)]
        assert row.native_units==len(ids) and row.finite_scores==len(finite) and row.coverage_status==cov.status
        assert row.positive_foreground_units==int((values['target_mean'][ids,j]>0).sum())
        assert np.isclose(row.mean_ULM,float(finite.mean()) if len(finite) else np.nan,equal_nan=True)
        assert np.isclose(row.SD_ULM,float(finite.std()) if len(finite) else np.nan,equal_nan=True)
        assert np.isclose(row.median_expressed_target_genes,float(np.median(values['expressed_target_genes'][ids,j])) if len(finite) else np.nan,equal_nan=True)
    assert rows_checked==len(axis)==terminal['source_score_rows']
    proof=dict(status='PASS_ALL1200_SCOPE_ALL_NATIVE_ROWS_ALL_OFFICIAL_ULM_SCORES_MEANS_DETECTION_AND_MISSINGNESS_QA',finished_UTC=datetime.now(timezone.utc).isoformat(),
               audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),native_rows_checked=rows_checked,official_score_values_checked=score_values,
               maximum_official_score_difference=maximum,all_raw_cells_independently_replayed=False,new_association_P=False)
    p=out/'private_independent_measurement_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()

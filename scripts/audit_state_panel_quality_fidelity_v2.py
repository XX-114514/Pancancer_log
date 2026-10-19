"""Independent source aggregation/triage and actual scalar weighted-LS checks."""
import hashlib,json,os,time,shutil
from datetime import datetime,timezone
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem
from scipy.stats import rankdata,spearmanr


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def actual_case(row, V, manifest, score_dir):
    info=manifest['private_native_files'][int(row.native_file_ordinal)]
    with h5py.File(V/'state_spatial_proxy_native_cells_20261016_v1'/info['relative_file'],'r') as h:
        obs=pd.DataFrame({k:read_elem(h['obs'][k]) for k in ['node_id','stratum','effective_resampling_block','full_RNA_total','full_detected_genes']})
    groups=pd.factorize(pd.MultiIndex.from_frame(obs[['node_id','stratum']]),sort=True)[0]
    ix=np.flatnonzero(groups==row.identity_stratum_ordinal)
    blocks=obs.effective_resampling_block.iloc[ix].astype(str).to_numpy()
    _,inv,cnt=np.unique(blocks,return_inverse=True,return_counts=True)
    assert len(ix)==row.cells and len(cnt)==row.analysis_blocks
    p=score_dir/'private_cell_scores'/(Path(info['relative_file']).stem+'_'+row.axis_id[:10]+'.npz')
    with np.load(p,allow_pickle=False) as z:
        j=z['programme_ids'].tolist().index(row.programme_id)
        a=z['proxy'][ix,j].astype(float);b=z['outside_reference'][ix,j].astype(float)
    assert np.isfinite(a).all() and np.isfinite(b).all()
    w=1/(len(cnt)*cnt[inv]);sw=np.sqrt(w)
    q=rankdata(obs[['full_RNA_total','full_detected_genes']].iloc[ix].to_numpy(float),axis=0)
    q-=np.average(q,axis=0,weights=w);q/=np.sqrt(np.average(q*q,axis=0,weights=w))
    X=np.column_stack([np.ones(len(ix)),q]);ar=rankdata(a);br=rankdata(b)
    ar-=X@np.linalg.lstsq(sw[:,None]*X,sw*ar,rcond=None)[0]
    br-=X@np.linalg.lstsq(sw[:,None]*X,sw*br,rcond=None)[0]
    cell=float(np.average(ar*br,weights=w)/np.sqrt(np.average(ar*ar,weights=w)*np.average(br*br,weights=w)))
    means=np.array([[ar[inv==k].mean(),br[inv==k].mean()] for k in range(len(cnt))])
    varying=(means.var(0)>np.array([np.average(ar*ar,weights=w),np.average(br*br,weights=w)])*1e-10).all()
    block=float(spearmanr(means[:,0],means[:,1]).statistic) if varying else np.nan
    np.testing.assert_allclose(cell,row.adjusted_cell_rho,atol=2e-8,rtol=1e-7)
    np.testing.assert_allclose(block,row.adjusted_block_rho,atol=2e-8,rtol=1e-7,equal_nan=True)
    return dict(case='positive' if cell>0 else 'negative',cell_error=abs(cell-row.adjusted_cell_rho),
                block_error=abs(block-row.adjusted_block_rho) if np.isfinite(block) else None)


def main():
    parent=Path('/dev/shm/pancancer-state-proxy-quality-20261020-v2/private_output')
    old=Path('/dev/shm/pancancer-state-proxy-quality-20261020-v1/private_output')
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT'])
    assert out.parent.parent==Path('/dev/shm') and not out.exists() and shutil.disk_usage('/dev/shm').free>512*1024**2
    for _ in range(360):
        if (parent/'private_terminal.json').is_file() and (old/'private_terminal.json').is_file():break
        time.sleep(5)
    else:raise TimeoutError('Same producer dependency, no restart performed')
    terminal=json.loads((parent/'private_terminal.json').read_text())
    assert terminal['producer_sha256']=='47a8774f645b0a759c31dce3626c1dfdff0fa9d1cf6e74d5c632dcac431692b2'
    assert terminal['plan_sha256']=='0839f00e7f3b3e81d72d648f7d9042fb805914e92bea5370c6a0855f37e97310'
    for folder in [old,parent]:
        for name,d in json.loads((folder/'private_terminal.json').read_text())['files'].items():assert sha(folder/name)==d['sha256']
    out.mkdir(mode=0o700)
    d=pd.read_csv(parent/'private_all_identity_quality_adjusted_fidelity.tsv.gz',sep='\t')
    earlier=pd.read_csv(old/'private_all_identity_quality_adjusted_fidelity.tsv.gz',sep='\t')
    assert len(d)==len(earlier)==725004
    for col in d:
        if col in ['status','adjusted_block_rho']:continue
        pd.testing.assert_series_equal(d[col],earlier[col],check_names=True)
    changed=earlier.adjusted_block_rho.notna()&d.adjusted_block_rho.isna()
    assert not (earlier.adjusted_block_rho.isna()&d.adjusted_block_rho.notna()).any()
    common=d.adjusted_block_rho.notna()
    np.testing.assert_array_equal(d.loc[common,'adjusted_block_rho'],earlier.loc[common,'adjusted_block_rho'])
    status_changed=d.status.ne(earlier.status)
    assert (changed==status_changed).all()
    assert earlier.loc[changed,'status'].eq('ESTIMABLE_TECHNICAL_DESCRIPTION').all()
    assert d.loc[changed,'status'].eq('CONSTANT_OR_QUALITY_EXPLAINED_SCORE').all()
    metrics=['raw_block_weighted_rho','adjusted_cell_rho','adjusted_block_rho','proxy_RNA_rho','proxy_detected_rho']
    source=pd.read_csv(parent/'private_all_source_major_programme_axis_fidelity.tsv.gz',sep='\t')
    keys=['gse_id','major','programme_id','axis_id'];index=source.set_index(keys);assert index.index.is_unique and len(index)==127636
    checked=0
    for key,ix in d.groupby(keys,sort=True).indices.items():
        g=d.iloc[ix];use=g[g.status.eq('ESTIMABLE_TECHNICAL_DESCRIPTION')];r=index.loc[key]
        assert r.identity_strata==len(g) and r.estimable_identity_strata==len(use)
        assert r.partial_signature_genes==g.partial_signature_genes.iloc[0]
        for col in metrics:
            value=float(np.median(use[col].to_numpy(float))) if len(use) else np.nan
            np.testing.assert_allclose(value,r[col],atol=2e-12,rtol=0,equal_nan=True)
        checked+=1
    plan=json.loads((parent/'plan_before_quality_adjusted_scores.json').read_text())
    dictionary=json.loads((V/plan['inputs']['dictionary']['directory']/plan['inputs']['dictionary']['manifest']).read_text())['programmes']
    bypid={r['programme_id']:r for r in dictionary}
    assert all(type(r['algorithm_stable_candidate']) is bool for r in dictionary)
    triage=pd.read_csv(parent/'private_complete_programme_major_axis_triage.tsv',sep='\t')
    seen=set();passed=0
    for row in triage.itertuples(index=False):
        key=(row.programme_id,row.major,row.axis_id);assert key not in seen;seen.add(key)
        rec=bypid[row.programme_id];assert row.major in rec['allowed_majors']
        g=source[source.programme_id.eq(row.programme_id)&source.major.eq(row.major)&source.axis_id.eq(row.axis_id)]
        s=g[g.estimable_identity_strata.gt(0)&g.partial_signature_genes.ge(5)]
        assert len(g)==row.available_sources and len(s)==row.qualifiable_sources and s.gse_id.is_unique
        a=s.adjusted_cell_rho.to_numpy();b=s.adjusted_block_rho.to_numpy();n=len(s)
        am=np.median(a) if n else np.nan;bm=np.median(b) if n else np.nan
        frac=np.mean((a>=.2)&(b>=.2)) if n else np.nan
        loo=min(min(np.median(np.delete(a,k)),np.median(np.delete(b,k))) for k in range(n)) if n>1 else np.nan
        for observed,expected in [(am,row.source_median_adjusted_cell_rho),(bm,row.source_median_adjusted_block_rho),(frac,row.source_fraction_both_ge_point2),(loo,row.leave_one_source_min_of_two_medians)]:
            np.testing.assert_allclose(observed,expected,atol=2e-12,rtol=0,equal_nan=True)
        result=bool(rec['algorithm_stable_candidate'] and n>=3 and am>=.3 and bm>=.3 and frac>=.75 and loo>=.3)
        assert result==row.exploratory_proxy_triage_pass;passed+=result
    expected={(r['programme_id'],major,axis) for r in dictionary for major in r['allowed_majors'] for axis in source.axis_id.unique()}
    assert seen==expected
    native=json.loads((V/plan['inputs']['native']['directory']/plan['inputs']['native']['manifest']).read_text())
    cases=[]
    for positive in [True,False]:
        c=d[d.status.eq('ESTIMABLE_TECHNICAL_DESCRIPTION')&((d.adjusted_cell_rho>.05) if positive else (d.adjusted_cell_rho<-.05))]
        if len(c):cases.append(actual_case(c.iloc[0],V,native,V/plan['inputs']['scores']['directory']))
    # Explicitly inspect a real excluded coverage case without interpreting it as negative.
    excluded=d[d.status.eq('FEWER_THAN50CELLS_OR10BLOCKS')]
    assert len(excluded)>0 and ((excluded.cells<50)|(excluded.analysis_blocks<10)).all()
    (out/'private_actual_scalar_cases.json').write_text(json.dumps(cases,indent=2)+'\n')
    report=dict(status='ACCEPTED_ALL_QUALITY_FIDELITY_SOURCE_AGGREGATES_TRIAGE_AND_ACTUAL_SCALAR_CHECKS',
                completed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),
                producer_terminal_sha256=sha(parent/'private_terminal.json'),old_unaccepted_terminal_sha256=sha(old/'private_terminal.json'),
                identity_positions_preserved=len(d),source_groups_independently_recomputed=checked,triage_positions_checked=len(seen),
                original_cell_metrics_and_all_supported_block_values_exact=True,real_scalar_cases_checked=len(cases),real_insufficient_coverage_positions_checked=len(excluded),
                full_RNA_or_score_arrays_copied_or_persisted=False,formal_state_or_interaction_validation=False,
                scientific_outputs_public=False,files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(v['bytes'] for v in report['files'].values())<1024**2
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['status','identity_positions_preserved','source_groups_independently_recomputed','triage_positions_checked','real_scalar_cases_checked']}))


if __name__=='__main__':main()

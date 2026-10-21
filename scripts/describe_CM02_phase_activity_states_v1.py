"""Outcome-free phase-specific frozen CM02/observed-state correspondence.

Read existing private model exports; learned arrays and new individual arrays
stay in RAM. This is descriptive source evidence, not a new clinical family.
"""
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import subprocess
import urllib.request

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.optimize import nnls, lsq_linear
from scipy.stats import rankdata

import measure_CM02_native_marker_author_correspondence_v3 as native
import observed_T_marker_evidence_v1 as evidence
import frozen_CM_network_contexts_v1 as network

SCOPE = 'T_whole_TME_primary_tumor'
STATES = ['R_TNFRSF9', 'R_LAG3', 'R_cycling',
          'E_dysfunction', 'E_cycling', 'E_memory']
PHASES = ['pre', 'post']
CONTEXTS = ['full_current_source', 'CM02_network_member_union']
B = 2000


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def manual_rank(x):
    values = sorted(set(x))
    ranks = {}
    position = 0
    for value in values:
        count = int(np.count_nonzero(x == value))
        ranks[value] = position + (count + 1) / 2
        position += count
    return np.array([ranks[value] for value in x])


def rho(x, y, independent=False):
    if len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        return None
    ranking = manual_rank if independent else rankdata
    a, b = ranking(x), ranking(y)
    a -= a.mean(); b -= b.mean()
    return float(np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b)))


def describe(x, y, slot):
    result = dict(reviewed_cases=len(x), unique_activity_values=len(np.unique(x)),
                  unique_state_fraction_values=len(np.unique(y)))
    if len(x) < 20:
        return dict(result, status='INSUFFICIENT_20_REVIEWED_CASES')
    effect = rho(x, y)
    if effect is None:
        return dict(result, status='CONSTANT_ACTIVITY_OR_STATE')
    assert abs(effect - rho(x, y, True)) < 1e-12
    seed = int.from_bytes(hashlib.sha256(slot.encode()).digest()[:8], 'little')
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(B):
        indices = rng.integers(len(x), size=len(x))
        first = rho(x[indices], y[indices])
        second = rho(x[indices], y[indices], True)
        assert (first is None) == (second is None)
        if first is not None:
            assert abs(first - second) < 1e-12
            boots.append(first)
    loo = [rho(np.delete(x, j), np.delete(y, j)) for j in range(len(x))] if len(x) - 1 >= 20 else []
    finite = [v for v in loo if v is not None]
    return dict(result, status='DESCRIPTIVE_SOURCE_PHASE_CORRESPONDENCE',
                spearman_rho=effect, bootstrap_valid=len(boots),
                bootstrap_CI=np.quantile(boots, [.025, .975]).tolist() if len(boots) >= 1800 else None,
                leave_one_case_min_max=[min(finite), max(finite)] if finite else None,
                leave_one_case_same_sign=sum(v * effect > 0 for v in finite),
                leave_one_case_valid=len(finite),
                bootstrap_manual_tie_rank_replay=True, no_P_or_FDR=True)


R_CHECK = r'''
args<-commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(NMF))
scope<-args[[1]];old<-args[[2]];helper<-args[[3]]
expressions<-parse(helper)
hit<-vapply(expressions,function(x)is.call(x)&&is.symbol(x[[1]])&&
 as.character(x[[1]])%in%c('<-','=')&&is.symbol(x[[2]])&&
 as.character(x[[2]])=='freq_normalize',logical(1))
stopifnot(sum(hit)==1L);eval(expressions[[which(hit)]])
f<-readRDS(file.path(scope,'discovery_v3','NMF_discovery.rds'))
W<-as.matrix(basis(f));H<-as.matrix(coef(f));S<-smoothing(fit(f));A<-W%*%S
table<-function(root,name)as.matrix(read.delim(file.path(root,name),row.names=1,check.names=FALSE))
fq<-read.delim(file.path(scope,'inputs','frequency_private.tsv'),row.names=1,check.names=FALSE)
norm<-freq_normalize(fq,normalize='minmax')
meta<-read.delim(file.path(scope,'inputs','sample_metadata_private.tsv'),check.names=FALSE)
meta<-meta[match(colnames(norm),meta$sampleID),,drop=FALSE]
stopifnot(!anyNA(meta$sampleID))
balanced<-sapply(unique(meta$independent_unit),function(b)rowMeans(norm[,meta$independent_unit==b,drop=FALSE]))
rownames(balanced)<-rownames(norm)
stopifnot(identical(dimnames(balanced),dimnames(fitted(f))))
stopifnot(identical(dimnames(W),dimnames(table(old,'W_private.tsv'))),
 identical(dimnames(H),dimnames(table(old,'H_private.tsv'))))
errors<-c(W=max(abs(W-table(old,'W_private.tsv'))),H=max(abs(H-table(old,'H_private.tsv'))),
 A=max(abs(A-table(old,'A_private.tsv'))),F=max(abs(fitted(f)-A%*%H)),
 norm=max(abs(norm-table(old,'norm_private.tsv'))))
stopifnot(all(errors<1e-8))
cat(jsonlite::toJSON(list(model_class=class(fit(f)),NMF_version=as.character(packageVersion('NMF')),
 rank=ncol(W),nodes=nrow(W),original_unit_order_verified=TRUE,
 original_discovery_RSS=sum((balanced-fitted(f))^2),max_replay_errors=as.list(errors),
 no_arrays_or_new_fits_written=TRUE),auto_unbox=TRUE))
'''


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    atlas = root.parents[1] / 'runs/20261013_immunotherapy_atlas_v1'
    scope = atlas / 'CM' / SCOPE
    old = root / 'HNSCC_new_CM_frozen_model_export_20261016_v2' / SCOPE
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    out.mkdir(mode=0o700)
    def save(name, value):
        (out / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    files = {name: old / name for name in ['A_private.tsv', 'H_private.tsv', 'F_private.tsv', 'norm_private.tsv', 'scaling_private.tsv']}
    files.update(frequency=scope/'inputs/frequency_private.tsv', metadata=scope/'inputs/sample_metadata_private.tsv',
                 original_model=scope/'discovery_v3/NMF_discovery.rds',
                 prior_baseline=root/'HNSCC_new_CM_baseline_preparation_20261016_v1/private_patient_baseline_component_scores.tsv')
    previous=json.loads((old.parent/'private_terminal.json').read_text())
    accepted=next(r for r in previous['scopes'] if r['scope']==SCOPE)
    assert sha(files['A_private.tsv'])==accepted['effective_basis_sha256']=='8fa2e933d2fd109da7c0d2e4a1a88bd1c31edaf3c93fb2ab09ed012792295da0'
    assert sha(files['norm_private.tsv'])==accepted['normalization_sha256']=='916efbace913881d2d3f9b683e0dda17356baaabf2c24edbb8a56f62aef740af'
    fixture=np.array([0.,0.,1.,1.,1.,2.])
    np.testing.assert_array_equal(manual_rank(fixture),rankdata(fixture))
    assert rho(fixture,np.ones(len(fixture))) is None
    assert abs(rho(fixture,-fixture)+1)<1e-12
    for y in [np.array([0.,0.]),np.array([0.,2.]),np.array([1.,3.])]:
        np.testing.assert_allclose(nnls(np.eye(2),y)[0],lsq_linear(np.eye(2),y,bounds=(0,np.inf),method='bvls').x,atol=1e-12)
    launch = json.loads((atlas/'manifests/CM_T_whole_TME_primary_tumor_recovery_v2_started_private.json').read_text())
    command = launch['command']
    assert command[-1] == 'discovery_v3' and Path(command[2]).resolve() == scope.resolve()
    assert sha(command[1]) == '4cf9af42279924f071ec1cffdbf49f1467cd4693aa85732bad6de9ad08485feb'
    helper = Path(command[3])/'R/CoVarNet.R'
    save('plan_before_phase_projection_and_state_correspondence.json', dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),
        scope=SCOPE,component='CM02',source='GSE286827',phases=PHASES,contexts=CONTEXTS,
        parents=native.PARENTS,states=STATES,descriptive_slots=48,
        input_SHA256={k:sha(p) for k,p in files.items()},actual_discovery_code_SHA256=sha(command[1]),
        original_helper_SHA256=sha(helper),
        activity='Original full-discovery minmax; reviewed patient-phase material-mean normalized frequencies; fixed original A=W S NNLS, full rank/condition<=1e8, no refit or new threshold.',
        state_fraction='All current source observations or original seven-member network union; state/own marker parent >=20 observations, >=20 reviewed cases per descriptive rho.',
        context='Activity always original full15-node model, not recomputed from the seven network members.',
        uncertainty='2000 whole-case bootstrap reranking; >=1800 valid for interval, leave-one-case requires >=20 remaining. Independently replay every bootstrap rank with manual ties.',
        sensitivity='State/observed current T count and parent/observed T count, same eligible cases. Cross-lineage marker parent observations may exceed original T annotation; never assert true T proportions.',
        interpretation='Same RNA and original composition discovery; correspondence is circular descriptive characterization, not independent state validation or clinical benefit.',
        old_baseline_family='Already executed/audited20261016, retained without new clinical tests.',
        no_outcomes_or_sequences_read=True,no_new_individual_or_learned_arrays_written=True))
    checked = subprocess.run([command[0], '--vanilla', '-e', R_CHECK, str(scope), str(old), str(helper)],
                             capture_output=True, text=True, timeout=120)
    assert checked.returncode == 0, 'R original model parity failed'
    contract = json.loads(checked.stdout)
    save('private_original_model_parity.json', contract)
    A = pd.read_csv(files['A_private.tsv'], sep='\t', index_col=0)
    H = pd.read_csv(files['H_private.tsv'], sep='\t', index_col=0)
    F = pd.read_csv(files['F_private.tsv'], sep='\t', index_col=0)
    norm = pd.read_csv(files['norm_private.tsv'], sep='\t', index_col=0)
    scale = pd.read_csv(files['scaling_private.tsv'], sep='\t', index_col=0)
    fq = pd.read_csv(files['frequency'], sep='\t', index_col=0)
    meta = pd.read_csv(files['metadata'], sep='\t', dtype=str, keep_default_na=False).set_index('sampleID')
    assert all(A.index.equals(x.index) for x in [fq,norm,scale])
    assert list(A.columns) == list(H.index) and list(norm.columns) == list(fq.columns[1:])
    assert set(norm.columns) == set(meta.index)
    np.testing.assert_allclose((fq.iloc[:,1:].to_numpy()-scale.minimum.to_numpy()[:,None])/scale['range'].to_numpy()[:,None], norm, atol=1e-12, rtol=1e-12)
    a=A.to_numpy(); condition=float(np.linalg.cond(a))
    assert np.linalg.matrix_rank(a)==a.shape[1] and condition<=1e8
    for j in range(H.shape[1]):
        np.testing.assert_allclose(nnls(a,F.iloc[:,j].to_numpy())[0],H.iloc[:,j],atol=1e-6,rtol=1e-7)
    fine = Path(os.environ['PANCANCER_FINE_LOADER'])
    ft = Path(os.environ['PANCANCER_FINE_SCOPE_OUTPUT'])/'private_terminal.json'
    loader = native.module_at('pinned_phase_reference',fine,json.loads(ft.read_text())['source_SHA256'])
    urllib.request.install_opener(urllib.request.build_opener(urllib.request.ProxyHandler({})))
    tab,frame,unique,shared,matrices,table=loader.load_inputs()
    frame['patient']=frame.library_id.map(tab.set_index('library_id').patient_registry)
    assert frame.patient.notna().all()
    source_meta=meta[meta.dataset.eq('GSE286827')].copy()
    source_meta['library']=source_meta.index.str.removeprefix('GSE286827::')
    assert set(source_meta.library)==set(tab.library_id)
    registry=tab.set_index('library_id')
    source_meta['patient']=source_meta.library.map(registry.patient_registry)
    for library,g in frame.groupby('library_id'):
        assert set(g.timepoint)=={source_meta.set_index('library').loc[library,'timepoint']}
    assert source_meta.groupby('independent_unit').patient.nunique().max()==1
    assert source_meta.groupby('patient').independent_unit.nunique().max()==1
    legacy=pd.read_csv(files['prior_baseline'],sep='\t')
    legacy=legacy[legacy.scope.eq(SCOPE)].set_index(['patient_id','component'])
    activities={};projection=[];j=list(A.columns).index('CM02')
    for (patient,phase),g in source_meta.groupby(['patient','timepoint'],sort=True):
        y=norm.loc[:,g.index].mean(axis=1).to_numpy()
        h,res=nnls(a,y,maxiter=10000)
        other=lsq_linear(a,y,bounds=(0,np.inf),method='bvls',tol=1e-12,max_iter=10000)
        assert other.success
        np.testing.assert_allclose(h,other.x,atol=1e-7,rtol=1e-7)
        grad=a.T@(a@h-y);tol=1e-7*(1+np.linalg.norm(a,2)*np.linalg.norm(y))
        assert grad.min()>=-tol and np.max(np.abs(h*grad))<=tol*(1+np.linalg.norm(h))
        if phase=='pre':
            np.testing.assert_allclose(h,[legacy.loc[(patient,c),'activity'] for c in A.columns],atol=1e-7,rtol=1e-7)
        activities[(patient,phase)]=float(h[j])
        projection.append(dict(phase=phase,activity=float(h[j]),residual=float(res/max(np.linalg.norm(y),1e-12))))
    contexts,membership=network.resolve(atlas,native.CONTEXTS,native.ANCHORS)
    z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r')
    selected=np.flatnonzero(np.asarray(read_elem(z['obs/dataset'])).astype(str)=='GSE286827')
    assert len(selected)==len(frame) and np.array_equal(frame.position,np.arange(len(frame)))
    var=read_elem(z['var']);sources=list(map(str,read_elem(z['uns/source_order'])))
    measured=np.asarray(read_elem(z['varm/measured_in_source']),bool)[:,sources.index('GSE286827')]
    human=measured&var.feature_identity_status.astype(str).eq('stable_Ensembl_ID').to_numpy()
    lookup=defaultdict(set)
    for col in np.flatnonzero(human):
        r=var.iloc[col]
        for value in {str(var.index[col]),str(r.gene_symbol)}|set(str(r.source_gene_aliases).split('|')):
            lookup[value].add(int(col))
    assert all(len(lookup[g])==1 for g in evidence.GENES)
    columns=[next(iter(lookup[g])) for g in evidence.GENES]
    counts,totals,detected=native.native_counts(z,selected,columns,human,measured,lambda v:save('private_progress.json',v))
    results=[]
    for parent in native.PARENTS:
        flags=evidence.vector_evidence(counts,np.ones(len(counts),bool),parent)
        np.testing.assert_array_equal(flags,evidence.scalar_evidence(counts,np.ones(len(counts),bool),parent))
        for context in CONTEXTS:
            mask=np.ones(len(frame),bool) if contexts[context] is None else frame.node_id.isin(contexts[context]).to_numpy()
            work=frame.loc[mask,['patient','timepoint','major']].copy()
            for k,name in enumerate(evidence.METRICS): work[name]=flags[mask,k].astype(int)
            work['current_T_count']=work.major.eq('T').astype(int)
            sums=work.groupby(['patient','timepoint'])[evidence.METRICS+['current_T_count']].sum()
            # Separate scalar dictionary count audit for all patient/phase/metric cells.
            independent=defaultdict(lambda:np.zeros(len(evidence.METRICS)+1,np.int64))
            for row,values in zip(work.itertuples(index=False,name=None),flags[mask]):
                independent[(row[0],row[1])]+=np.r_[values.astype(np.int64),int(row[2]=='T')]
            assert set(independent)==set(sums.index)
            for key in sums.index: np.testing.assert_array_equal(sums.loc[key].to_numpy(),independent[key])
            for phase in PHASES:
                current=sums[sums.index.get_level_values('timepoint')==phase]
                for state in STATES:
                    parent_name='regulatory_only' if state.startswith('R_') else 'CD8_only'
                    use=current[(current[parent_name]>=20)&(current.current_T_count>=20)]
                    x=np.array([activities[key] for key in use.index])
                    y=(use[state]/use[parent_name]).to_numpy()
                    slot='|'.join([parent,context,phase,state])
                    result=describe(x,y,slot)
                    result.update(source='GSE286827',phase=phase,parent=parent,context=context,state=state,
                                  projected_phase_cases=int(len(current)),
                                  parent_or_T_support_excluded=int(len(current)-len(use)),
                                  current_T_state_abundance_rho=rho(x,(use[state]/use.current_T_count).to_numpy()) if len(x)>=20 else None,
                                  current_T_parent_abundance_rho=rho(x,(use[parent_name]/use.current_T_count).to_numpy()) if len(x)>=20 else None)
                    results.append(result)
    assert len(results)==48
    save('private_all48_phase_activity_state_correspondences.json',results)
    profiles=[]
    for phase in native.PHASES:
        rows=[r for r in projection if r['phase']==phase]
        profiles.append(dict(source='GSE286827',phase=phase,reviewed_cases=len(rows),
            positive_CM02_activity_cases=sum(r['activity']>0 for r in rows),
            activity_min_q25_median_q75_max=np.quantile([r['activity'] for r in rows],[0,.25,.5,.75,1]).tolist(),
            relative_projection_residual_min_q25_median_q75_max=np.quantile([r['residual'] for r in rows],[0,.25,.5,.75,1]).tolist()))
    save('private_source_phase_projection_profiles.json',profiles)
    save('private_terminal.json',dict(status='COMPLETE_OUTCOME_FREE_PHASE_ACTIVITY_STATE_CORRESPONDENCE',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),scope=SCOPE,
        actual_original_producer_and_original_unit_order_verified=True,effective_condition=condition,
        legacy20261016_baseline_scores_all_exact_replayed=True,all_phase_NNLS_BVLS_KKT_passed=True,
        original_fitted_H_all_replayed=True,all_source_marker_CSR_and_scalar_truth_table_replayed=True,
        all_case_phase_state_counts_independently_replayed=True,all_descriptive_bootstrap_manual_rank_replayed=True,
        descriptive_slots=len(results),no_clinical_outcomes_read=True,no_clinical_or_spatial_tests_reopened=True,
        no_original_models_changed=True,no_new_individual_or_learned_arrays_written=True,
        files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status='COMPLETE_OUTCOME_FREE_PHASE_ACTIVITY_STATE_CORRESPONDENCE',terminal_SHA256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

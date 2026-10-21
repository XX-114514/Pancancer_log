"""Complete fixed observed-state pair descriptions, with RNA capture controls.

All individual arrays stay in server RAM. Existing patients, source RNA and
network membership are read-only. No new clinical endpoint is selected.
"""
from collections import defaultdict
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
import hashlib
import json
import os
import traceback
import urllib.request

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem

import observed_T_marker_evidence_v1 as evidence
import frozen_CM_network_contexts_v1 as network
import measure_CM02_native_marker_author_correspondence_v3 as native
import observed_T_pair_statistics_v1 as statistics

STATES = ['R_TNFRSF9','R_LAG3','R_cycling','E_dysfunction','E_cycling','E_memory']
CONTEXTS = ['full_current_source','CM02_network_member_union']
DEPTHS = ['all_native_filtered','low_native_capture','mid_native_capture','high_native_capture']
PHASES = ['pre','post']
PHASE_COLUMNS = ['patient','timepoint']


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def count_scope(frame, flags, mask, bg_umi, bg_genes):
    positions = np.flatnonzero(mask)
    w = frame.iloc[positions][PHASE_COLUMNS].copy()
    w[evidence.METRICS] = flags[positions].astype(np.int64)
    sums = w.groupby(PHASE_COLUMNS)[evidence.METRICS].sum()
    independent = defaultdict(lambda: np.zeros(len(evidence.METRICS),np.int64))
    for key,values in zip(w[PHASE_COLUMNS].itertuples(index=False,name=None),flags[positions]):
        independent[key] += values.astype(np.int64)
    assert set(independent) == set(sums.index)
    for key in sums.index:
        np.testing.assert_array_equal(sums.loc[key].to_numpy(),independent[key])
    qualities = {}
    for parent in ['regulatory_only','CD8_only']:
        take = positions[flags[positions,evidence.METRICS.index(parent)]]
        q = frame.iloc[take][PHASE_COLUMNS].copy()
        q['log1p_background_UMI'] = np.log1p(bg_umi[take])
        q['log1p_background_genes'] = np.log1p(bg_genes[take])
        medians = q.groupby(PHASE_COLUMNS)[['log1p_background_UMI','log1p_background_genes']].median()
        groups = defaultdict(list)
        for key,i in zip(q[PHASE_COLUMNS].itertuples(index=False,name=None),take):
            groups[key].append(int(i))
        assert set(groups) == set(medians.index)
        for key,ids in groups.items():
            raw = np.column_stack([np.log1p(bg_umi[ids]),np.log1p(bg_genes[ids])])
            np.testing.assert_allclose(np.median(raw,axis=0),medians.loc[key].to_numpy(),atol=1e-12,rtol=1e-12)
            qualities[(key,parent)] = medians.loc[key].to_numpy()
    return sums,qualities


def materialize(sums, qualities, metadata, phase):
    current = sums[sums.index.get_level_values('timepoint')==phase]
    keys=[];excluded=defaultdict(int)
    for key,row in current.iterrows():
        if min(row.regulatory_only,row.CD8_only,row.T_support)<20:
            excluded['parent_or_marker_T_below20']+=1
            continue
        arm,hpv=metadata.loc[key,['neo_type','HPV/p16']]
        if arm not in ['D','D+T'] or hpv not in ['positive','negative']:
            excluded['drug_or_HPV_not_qualified']+=1
            continue
        assert (key,'regulatory_only') in qualities and (key,'CD8_only') in qualities
        keys.append(key)
    keys.sort()
    n=len(keys)
    counts=current.loc[keys] if n else current.iloc[:0]
    nums=counts[STATES].to_numpy(float)
    nR=counts.regulatory_only.to_numpy(float);nE=counts.CD8_only.to_numpy(float)
    T=counts.T_support.to_numpy(float)
    g=metadata.loc[keys] if n else metadata.iloc[:0]
    tx=g.neo_type.eq('D+T').to_numpy(int);hpv=g['HPV/p16'].eq('positive').to_numpy(int)
    quality=np.array([np.r_[qualities[(k,'regulatory_only')],qualities[(k,'CD8_only')]] for k in keys]).reshape(n,4)
    # Input quality order R_UMI,R_genes,E_UMI,E_genes -> declared control order.
    c=np.column_stack([np.ones(n),np.log1p(nR),np.log1p(nE),tx,hpv,quality[:,0],quality[:,2],quality[:,1],quality[:,3]])
    y=nums/np.column_stack([np.repeat(nR[:,None],3,axis=1),np.repeat(nE[:,None],3,axis=1)])
    sensitivity=dict(state_per_marker_T=nums/T[:,None],
        parent_per_marker_T=np.column_stack([np.repeat((nR/T)[:,None],3,axis=1),np.repeat((nE/T)[:,None],3,axis=1)]))
    assert np.isfinite(c).all() and np.isfinite(y).all() and np.all((y>=0)&(y<=1))
    report=dict(source_phase_cases=len(current),eligible_reviewed_cases=n,exclusions=dict(excluded),
        eligible_by_arm=[dict(arm=a,reviewed_cases=int(g.neo_type.eq(a).sum())) for a in ['D','D+T']],
        background_capture_median_profiles={name:np.quantile(quality[:,j],[0,.25,.5,.75,1]).tolist() if n>=5 else None
            for j,name in enumerate(['R_logUMI','R_logGenes','E_logUMI','E_logGenes'])},
        control_columns=9,independent_patient_not_library_or_cell_unit=True)
    return y,c,2*tx+hpv,sensitivity,report


def main():
    out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent==Path('/dev/shm') and not out.exists()
    out.mkdir(mode=0o700)
    def save(name,value):
        (out/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    fine=Path(os.environ['PANCANCER_FINE_LOADER'])
    scope=Path(os.environ['PANCANCER_FINE_SCOPE_OUTPUT'])
    terminal=json.loads((scope/'private_terminal.json').read_text())
    loader=native.module_at('pinned_T_pair_patient_loader',fine,terminal['source_SHA256'])
    assert sha(evidence.__file__)=='e8950f189eb399a4aa647e9e36c96dd4f2d0733404d47aa5d5f72b9d02922ae6'
    save('plan_before_observed_state_pair_effects.json',dict(frozen_UTC=datetime.now(timezone.utc).isoformat(),
        producer_SHA256=sha(__file__),statistics_SHA256=sha(statistics.__file__),recipe_SHA256=sha(evidence.__file__),
        source='GSE286827',component='2026-10-14 immune atlas CM02 v1 T_whole_TME_primary_tumor',
        phases=PHASES,contexts=CONTEXTS,parents=native.PARENTS,depth=DEPTHS,states=STATES,all_positions=288,
        main='pre/full_current_source/dual_regulatory_evidence/all_native_filtered; all9 state pairs',
        secondary='post same cohort, all predetermined contexts/parents/depth; never promote a sensitivity',
        unit='Reviewed patient-phase, all libraries pooled. Both parents and marker-T>=20 observations, >=20 patients.',
        controls=['intercept','rank_log1p_R_count','rank_log1p_E_count','DplusT','HPVpositive','rank_R_background_logUMI','rank_E_background_logUMI','rank_R_background_logGenes','rank_E_background_logGenes'],
        capture='Patient/context/depth/parent median cellular log1p human RNA counts/detection, excluding all22 recipe genes; endogenous proxies, not pure technical truth.',
        design='All9 columns, full rank, condition<1e6, >=10 residual df including tested state. No control dropping.',
        uncertainty='2000 within drug x HPV whole-patient bootstrap, every draw reranked, >=1800 valid for pointwiseCI; leave-one still>=20 patients.',
        numeric='Every point/draw/LOO/sensitivity SciPy rank+lstsq+QR versus independent Pandas average rank+SVD.',
        source_scope='Descriptive characterization within original discovery data; neither independent validation nor direct communication.',
        external='GSE288199 clinical bridge deferred, no state-pair effects in unreviewed coded units.',
        no_new_P_FDR_endpoints_or_clinical_models=True,no_new_individual_arrays_saved=True))
    save('private_synthetic_statistics_proof.json',statistics.synthetic())
    urllib.request.install_opener(urllib.request.build_opener(urllib.request.ProxyHandler({})))
    tab,frame,unique,shared,matrices,table=loader.load_inputs()
    assert sha(table)==terminal['original_author_join_SHA256']
    frame['patient']=frame.library_id.map(tab.set_index('library_id').patient_registry)
    assert frame.patient.notna().all()
    assert (frame.timepoint==frame.library_id.map(tab.set_index('library_id').timepoint_registry)).all()
    metadata=tab[['patient_registry','timepoint_registry','neo_type','HPV/p16']].drop_duplicates()
    assert not metadata.duplicated(['patient_registry','timepoint_registry']).any()
    metadata=metadata.rename(columns=dict(patient_registry='patient',timepoint_registry='timepoint')).set_index(PHASE_COLUMNS)
    assert set(zip(frame.patient,frame.timepoint))<=set(metadata.index)
    atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1'
    z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r')
    selected=np.flatnonzero(np.asarray(read_elem(z['obs/dataset'])).astype(str)=='GSE286827')
    assert len(selected)==len(frame) and np.array_equal(frame.position,np.arange(len(frame)))
    var=read_elem(z['var']);sources=list(map(str,read_elem(z['uns/source_order'])))
    measured=np.asarray(read_elem(z['varm/measured_in_source']),bool)[:,sources.index('GSE286827')]
    human=measured&var.feature_identity_status.astype(str).eq('stable_Ensembl_ID').to_numpy()
    lookup=defaultdict(set)
    for col in np.flatnonzero(human):
        row=var.iloc[col]
        for alias in {str(var.index[col]),str(row.gene_symbol)}|set(str(row.source_gene_aliases).split('|')):
            lookup[alias].add(int(col))
    assert all(len(lookup[g])==1 for g in evidence.GENES)
    columns=[next(iter(lookup[g])) for g in evidence.GENES]
    assert len(set(columns))==len(columns)
    contexts,membership=network.resolve(atlas,native.CONTEXTS,native.ANCHORS)
    assert set(CONTEXTS)<=set(contexts)
    save('private_input_proofs.json',dict(author_join_SHA256=sha(table),fine_terminal_SHA256=sha(scope/'private_terminal.json'),
        authority_release_SHA256=sha(atlas/'manifests/immune_authority_release_public.json'),network=membership,
        no_new_response_TCR_or_clinical_endpoint_values_selected=True))
    counts,totals,detected=native.native_counts(z,selected,columns,human,measured,lambda v:save('private_progress.json',v))
    bg_umi=totals-counts.sum(axis=1);bg_genes=detected-(counts>0).sum(axis=1)
    assert (bg_umi>=0).all() and (bg_genes>=0).all()
    depth=np.full(len(frame),-1,np.int8)
    for lib in frame.library_id.unique():
        ix=np.flatnonzero(frame.library_id.eq(lib).to_numpy())
        cuts=np.unique(np.quantile(totals[ix],[1/3,2/3]))
        depth[ix]=np.searchsorted(cuts,totals[ix],side='left')
    results=[];profiles=[];completed=0
    for parent in native.PARENTS:
        flags=evidence.vector_evidence(counts,np.ones(len(counts),bool),parent)
        np.testing.assert_array_equal(flags,evidence.scalar_evidence(counts,np.ones(len(counts),bool),parent))
        assert (flags[:,1:5].sum(axis=1)==1).all()
        for context,di in product(CONTEXTS,range(4)):
            mask=np.ones(len(frame),bool) if contexts[context] is None else frame.node_id.isin(contexts[context]).to_numpy()
            if di:mask&=depth==(di-1)
            sums,qualities=count_scope(frame,flags,mask,bg_umi,bg_genes)
            for phase in PHASES:
                y,c,strata,sensitivity,coverage=materialize(sums,qualities,metadata,phase)
                label='|'.join([parent,context,DEPTHS[di],phase])
                values=statistics.profile(y,c,strata,sensitivity,label)
                for row in values:
                    row.update(source='GSE286827',phase=phase,parent=parent,context=context,capture_scope=DEPTHS[di],
                        regulatory_state=STATES[row.pop('i')],effector_state=STATES[3+row.pop('j')])
                    results.append(row)
                profiles.append(dict(source='GSE286827',phase=phase,parent=parent,context=context,capture_scope=DEPTHS[di],**coverage))
                completed+=1
                save('private_progress.json',dict(stage='complete_scope_patient_resampling_and_independent_rank_SVD_replay',completed=completed,total=32,observed_UTC=datetime.now(timezone.utc).isoformat()))
    assert len(results)==288 and len(profiles)==32
    save('private_all288_observed_T_pair_capture_descriptions.json',results)
    save('private_all32_source_scope_and_quality_profiles.json',profiles)
    save('private_numeric_proofs.json',dict(all_native_raw_CSR_and_marker_scalar_counts=True,
        all_patient_phase_counts_dictionary_checked=True,all_parent_capture_medians_Pandas_and_Numpy_checked=True,
        all_actual_points_bootstrap_LOO_and_sensitivity_independent_ranks_QR_SVD_checked=True,
        no_effect_below20_reviewed_cases_or_without_full_controls=True,all288_opportunities_retained=True,
        no_P_FDR_clinical_or_spatial_test=True,no_individual_arrays_saved=True))
    save('private_terminal.json',dict(status='COMPLETE_ALL288_OBSERVED_T_PAIR_CAPTURE_DESCRIPTIONS',
        completed_UTC=datetime.now(timezone.utc).isoformat(),producer_SHA256=sha(__file__),statistics_SHA256=sha(statistics.__file__),
        all288_complete_numeric_and_scope_checks=True,no_clinical_endpoints_or_P_FDR=True,
        no_cell_patient_expression_control_draw_or_learned_arrays_saved=True,
        original_CM_models_authorities_and_closed_families_unchanged=True,
        files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status='COMPLETE_ALL288_OBSERVED_T_PAIR_CAPTURE_DESCRIPTIONS',terminal_SHA256=sha(out/'private_terminal.json'))))


if __name__=='__main__':
    try:
        main()
    except Exception as e:
        out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
        if out.is_dir():
            proof=dict(status='FAILED_NOT_SCIENTIFICALLY_ACCEPTED',observed_UTC=datetime.now(timezone.utc).isoformat(),
                error_class=type(e).__name__,frames=[dict(file=Path(f.filename).name,line=f.lineno,function=f.name) for f in traceback.extract_tb(e.__traceback__)],
                no_individual_values_in_error_record=True)
            (out/'private_failed_terminal.json').write_text(json.dumps(proof,indent=2)+'\n')
            print(json.dumps(dict(status=proof['status'],error_class=proof['error_class'])))
        raise SystemExit(1)

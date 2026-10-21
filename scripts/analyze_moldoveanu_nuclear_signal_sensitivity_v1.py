"""Fixed nuclear-proxy sensitivity, distinct from the closed constant-control run.

Histone H3 and the equal mean of standardized Ir191/Ir193 rank signals are
declared nuclear covariates, not pure technical negative controls. Never infer
artifact or absence of real cycling from attenuation. No fitted coefficients
or individual rows leave RAM. Twelve positions, one common set per relation.
"""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import importlib.util
import io
import json
import numpy as np
import pandas as pd

PATH=Path(__file__).with_name('analyze_moldoveanu_matched_representation_null_v1.py')
assert hashlib.sha256(PATH.read_bytes()).hexdigest()=='9698a399640f7be594b487c28647ed95b3b2a06dcdec1b82dd720c530c5c8302'
SPEC=importlib.util.spec_from_file_location('fixed_matched_nuclear_parent',PATH)
core=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(core)
base=core.base
METRICS=['raw_Ki67_rank_covariance','nuclear_residual_Ki67_covariance','DNA_proxy_covariance']
SEED=2026102203
B=2000


def nuclear_scores(part,rank):
    if len(part)<20:return None,'parent_below_minimum'
    values=part[['Ki67','Histone.H3','X191Ir','X193Ir']].to_numpy(dtype=float)
    if not np.isfinite(values).all():return None,'nonfinite_required_channel'
    standardized,varied=rank(values)
    if not varied.all():return None,'constant_required_channel_no_automatic_dropping'
    dna=standardized[:,2:4].mean(axis=1)
    dna_rms=float(np.sqrt(np.mean(dna**2)))
    if dna_rms<1e-10:return None,'DNA_proxy_without_variation'
    dna=dna/dna_rms
    design=np.column_stack([np.ones(len(part)),standardized[:,1],dna])
    if np.linalg.matrix_rank(design)!=3:return None,'nuclear_design_rank_deficient'
    if np.linalg.cond(design)>1e8:return None,'nuclear_design_ill_conditioned'
    coefficients=np.linalg.lstsq(design,standardized[:,0],rcond=None)[0]
    residual=standardized[:,0]-design@coefficients
    for column in design.T:
        assert abs(sum(float(x)*float(y) for x,y in zip(column,residual)))<1e-7*len(part)
    scale=float(np.sqrt(np.mean(residual**2)))
    if scale<=1e-10:return None,'no_residual_Ki67_variation'
    return np.column_stack([standardized[:,0],residual/scale,dna]),None


def prepare(frame,rank):
    prepared=[];counts=[Counter() for _ in range(3)]
    for _,part in frame.groupby('sample.id',sort=True):
        part=part.reset_index(drop=True)
        positions=[np.flatnonzero(part['Cluster'].isin(labels)) for labels in base.PARENTS]
        scores=[];reasons=[]
        for role,position in enumerate(positions):
            value,reason=nuclear_scores(part.iloc[position],rank)
            counts[role][reason or 'all_declared_nuclear_scores_identifiable']+=1
            scores.append(value);reasons.append(reason)
        prepared.append(dict(coords=part[['coord.x','coord.y']].to_numpy(dtype=float),
            tumor=part['Cluster'].eq('melano').to_numpy(),positions=positions,scores=scores,reasons=reasons))
    profile=dict(status='COMPLETE_DECLARED_NUCLEAR_DESIGN_QUALIFICATION_BEFORE_EFFECTS',
        roles=[dict(role=name,reason_counts=dict(c)) for name,c in zip(['regulatory','conventional_CD4','CD8'],counts)],
        raw_channel_names=['Ki67','Histone.H3','X191Ir','X193Ir'],
        fixed_model='rank(Ki67) ~ intercept + rank(Histone.H3) + unit_RMS(mean(rank(Ir191),rank(Ir193)))',
        no_learned_coefficients_or_individual_arrays_exported=True,no_state_effects_yet=True)
    return prepared,profile


def covariance_scores(left,right,a,b):
    return np.array([base.covariance(left[:,j],right[:,j],a,b) for j in range(3)])


def statistics(prepared,knn,grid):
    records=[[],[]];coverage=[Counter(),Counter()]
    rngs=[np.random.default_rng(SEED),np.random.default_rng(SEED+100000)]
    for entry in prepared:
        coords=entry['coords'];positions=entry['positions'];values=entry['scores'];reasons=entry['reasons']
        if values[2] is None or all(v is None for v in values[:2]):
            for parent in range(2):coverage[parent][reasons[parent] or reasons[2] or 'unidentifiable_other_parent']+=1
            continue
        edges=knn(coords);near=base.tumor_neighbors(edges,entry['tumor'])
        tiles=grid(coords,(coords.min(0),coords.max(0)))
        pairs=[base.pair_indices(edges,positions[i],positions[2],len(coords)) for i in range(2)]
        valid=[]
        for parent in range(2):
            reason=reasons[parent] or reasons[2]
            if reason is None and len(pairs[parent][0])<20:reason='cross_edges_below_minimum'
            valid.append(reason is None)
            if reason is not None:coverage[parent][reason]+=1
        if not any(valid):continue
        observed=[covariance_scores(values[i],values[2],*pairs[i]) if valid[i] else None for i in range(2)]
        expected=[np.zeros((2,3)) for _ in range(2)];draws=[np.zeros((2,B,3)) for _ in range(2)]
        for mode in range(2):
            strata=tiles if mode==0 else tiles*2+near.astype(int)
            groups=[strata[position] for position in positions]
            means=[np.column_stack([base.means_by_group(v[:,j],g) for j in range(3)]) if v is not None else None
                   for v,g in zip(values,groups)]
            for parent in range(2):
                if valid[parent]:expected[parent][mode]=covariance_scores(means[parent],means[2],*pairs[parent])
            for j in range(B):
                permutations=[core.permute_indices(g,rngs[mode]) if values[i] is not None else None for i,g in enumerate(groups)]
                for parent in range(2):
                    if valid[parent]:draws[parent][mode,j]=covariance_scores(values[parent][permutations[parent]],
                        values[2][permutations[2]],*pairs[parent])
        for parent in range(2):
            if valid[parent]:records[parent].append(dict(observed=observed[parent],expected=expected[parent],draws=draws[parent]))
    # Reuse the unchanged scalar-replayed generic three-column summary. Only
    # its labels are configured; all gates, P/BH rules and twelve slots persist.
    old=core.METRICS
    try:
        core.METRICS=METRICS
        result=core.summarize(records,coverage)
    finally:core.METRICS=old
    return dict(family=result,source_materials_considered=len(prepared),all_three_scores_jointly_permuted=True,
        each_relationship_has_identical_materials_graph_and_edges_across_six_positions=True)


def analyze(bundle,parameters):
    assert set(parameters)=={'scope','permutations','seed','checkpoint_directory'}
    assert parameters['scope']=='fixed_nuclear_proxy_sensitivity_no_clinical'
    assert parameters['seed']==SEED and parameters['permutations']==B
    folder=Path(__file__).parent
    qualifier=base.load(folder/'qualify_moldoveanu_ICI_cells_in_RAM_v1.py',
        'a498faeb15053036e72003585ba82d5c9b44b2f8f9d3c16d392c3fcf8e693125','nuclear_qualified_cells')
    source=qualifier.analyze(bundle,dict(coordinate_tolerance=1e-6,minimum_parent_cells=20))
    archive=bundle['archive'];frame=archive._own_ICI_qualified_cells
    raw=qualifier.decoded(qualifier.table_bytes(archive,'ICI_meanIntensity.tsv.gz'))
    fields=['Cell_ID','Sample_ID','Cluster.v2','prob.FoxP3','prob.CD45RO','Histone.H3','X191Ir','X193Ir']
    extra=pd.read_csv(io.BytesIO(raw),sep='\t',usecols=fields,dtype={'Cell_ID':str,'Sample_ID':str})
    extra=extra.drop(columns=[x for x in ['Histone.H3','X191Ir','X193Ir'] if x in frame])
    checked=frame.merge(extra,left_on=['sample.id','obj.id'],right_on=['Sample_ID','Cell_ID'],
        how='left',validate='one_to_one',suffixes=('','_raw'),indicator='_nuclear_join')
    assert len(checked)==len(frame) and checked['_nuclear_join'].eq('both').all()
    assert hashlib.sha256(qualifier.table_bytes(archive,'cell_type_clustering.R')).hexdigest()=='2d38c485cfdc5e5292d92d5f4ce4f11a08c95767a2eaa382aa92a077c840487f'
    audit=base.audit_final_author_phenotypes(checked);assert audit['author_probability_recipes_all_passed']
    helper=base.load(folder/'analyze_melanoma_IMC_continuous_state_neighbors_v2.py',
        '51b8a4323d806ae011ecb7d9419d56548eba63ca947a6756d3ab1dc5f9baf7bf','nuclear_rank_grid')
    neighbors=base.load(folder/'melanoma_IMC_numpy_neighbors_v1.py',
        'acbde69dae387ee79d6d8397899681bedf00553dd3da08878abd16b392858471','nuclear_neighbors')
    prepared,profile=prepare(checked,helper.rank_z);profile['completed_UTC']=datetime.now(timezone.utc).isoformat()
    out=Path(parameters['checkpoint_directory']);assert out.is_dir() and not out.is_symlink()
    with (out/'private_nuclear_design_qualification_v1.json').open('x') as handle:json.dump(profile,handle,indent=2,allow_nan=False)
    result=statistics(prepared,lambda c:base.exact_edges(c,neighbors),helper.tile_groups)
    return dict(status='COMPLETE_DECLARED_NUCLEAR_SIGNAL_SENSITIVITY_NOT_CAUSAL_OR_CLINICAL_VALIDATION',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_qualification=source,final_type_audit=audit,
        preeffect_nuclear_design_qualification=profile,results=result,family_size=12,permutations=B,seed=SEED,
        all_required_nuclear_covariates_fixed_no_automatic_dropping=True,
        constant_X190BCKG_closed_family_not_recomputed_or_silently_changed=True,
        raw_and_residual_and_DNA_same_eligible_materials=True,
        nuclear_signals_may_contain_real_biology_not_pure_negative_controls=True,
        residual_attenuation_not_proof_of_artifact_or_absent_proliferation=True,
        effect_differences_descriptive_no_formal_superiority_test=True,
        no_patient_CI_clinical_model_or_original_CM02_projection=True,
        arrays_and_fitted_coefficients_RAM_only=True,producer_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())

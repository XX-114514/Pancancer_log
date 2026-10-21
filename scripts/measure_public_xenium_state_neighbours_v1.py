"""Frozen observed-state spatial question; no clinical outcomes or label learning.

Raw source bytes, cell arrays, and author-case mapping live in process RAM only.
Source-only summaries may be reviewed privately. Slot-level checks stay server.
"""
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
import hashlib
import json
import os

import numpy as np
from scipy.spatial import cKDTree
import xenium_observed_state_context_v2 as context
import xenium_state_spatial_kernel_v1 as kernel

MODES = list(product([20,50], ['dual_regulatory_evidence','costate_regulatory_sensitivity']))
RADII = [25,50,100]
NULLS = ['parent_only','local_matched']
B = 2000


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def matched_groups(data, masks, parent):
    xy = data['coordinates'][['x_centroid','y_centroid']].to_numpy(float)
    choose = masks[parent]
    positions = xy[choose]
    # All observed T positions (including unresolved parents) are fixed context.
    density = cKDTree(xy[masks['T_support']]).query_ball_point(
        positions, 50, return_length=True) - 1
    assert (density >= 0).all()
    ix = {g:i for i,g in enumerate(context.REQUIRED)}
    epi = masks['qc'] & ~masks['T_support'] & (data['counts'][:,ix['EPCAM']]>=1)
    if epi.any():
        distance = cKDTree(xy[epi]).query(positions,k=1)[0]
        close = distance<=50
    else:
        close = np.zeros(len(positions),bool)
    keys = np.column_stack([np.floor(positions/500).astype(np.int64),
        kernel.quantile_bins(density),
        kernel.quantile_bins(data['human_gene_molecules'][choose]),close.astype(int)])
    groups = kernel.groups_from_keys(keys)
    return groups


def eligibility(masks):
    nr = int(masks['regulatory'].sum()); ne = int(masks['CD8'].sum())
    rs = masks['R_states'].sum(axis=0); es = masks['E_states'].sum(axis=0)
    allowed = np.ones((3,3),bool)
    statuses = np.full((3,3),'OBSERVED_STATE_COUNT_SCOPE_READY',object)
    for i,j in product(range(3),range(3)):
        s = 'OBSERVED_STATE_COUNT_SCOPE_READY'
        if nr<20 or ne<20: s='INSUFFICIENT_PARENT_OBSERVATIONS'
        elif rs[i]<5 or es[j]<5: s='INSUFFICIENT_OBSERVED_CHILD_SUPPORT'
        elif rs[i]==nr or es[j]==ne: s='NO_WITHIN_PARENT_STATE_VARIATION'
        statuses[i,j]=s; allowed[i,j]=s=='OBSERVED_STATE_COUNT_SCOPE_READY'
    return allowed,statuses


def one_section(ctx, ordinal, support):
    metadata = ctx['pairs'][ordinal]['source_metadata']
    if metadata['diagnosis']!='HNSCC':
        assert support['excluded']
        return dict(ordinal=ordinal,excluded=True,metadata=metadata)
    d = context.load_pair(ctx,ordinal)
    xy = d['coordinates'][['x_centroid','y_centroid']].to_numpy(float)
    all_modes=[]; technical=[]
    for mi,(minimum,parent) in enumerate(MODES):
        masks = context.observed_masks(d,minimum,parent)
        sm = support['modes'][mi]
        measured={k:int(masks[k].sum()) for k in ['qc','T_support','regulatory','CD8','ambiguous','B_cycling','epithelial_cycling']}
        assert measured==sm['counts']
        assert masks['R_states'].sum(axis=0).tolist()==sm['R_observed_state_counts']
        assert masks['E_states'].sum(axis=0).tolist()==sm['E_observed_state_counts']
        count_ready,statuses = eligibility(masks)
        assert statuses.ravel().tolist()==[v['status'] for v in sm['opportunities']]
        nr=int(masks['regulatory'].sum()); ne=int(masks['CD8'].sum())
        mode=dict(results=[],count_statuses=statuses)
        if count_ready.any():
            left=masks['R_states'][masks['regulatory']].astype(np.int64)
            right=masks['E_states'][masks['CD8']].astype(np.int64)
            graphs=[kernel.graph_for_radius(xy[masks['regulatory']],xy[masks['CD8']],r) for r in RADII]
            for ni,null in enumerate(NULLS):
                lg=[np.arange(nr)] if null=='parent_only' else matched_groups(d,masks,'regulatory')
                rg=[np.arange(ne)] if null=='parent_only' else matched_groups(d,masks,'CD8')
                result=kernel.measure_graphs(left,right,lg,rg,graphs,B,[20261021,ordinal,mi,ni])
                ready=count_ready[None,...]&result['positive_expected']&result['empirical_null_varies']
                denominator=result['expected']+.5
                result['effect']=(result['observed']-result['expected'])/denominator
                result['null_effect']=(result['draws']-result['expected'][None,...])/denominator[None,...]
                result['ready']=ready
                # Draws and arrays are RAM-only, never a cell/case table on disk.
                mode['results'].append(result)
                technical.append(dict(mode_index=mi,null=null,groups_left=len(lg),groups_right=len(rg),
                    singleton_left=sum(len(g)==1 for g in lg),singleton_right=sum(len(g)==1 for g in rg),
                    edges_by_radius=[g[0].nnz for g in graphs],
                    independent_observed_and_expected_and_first8_draws=True,
                    joint_vectors_and_source_measurement_counts_checked=True))
        all_modes.append(mode)
    return dict(ordinal=ordinal,excluded=False,metadata=metadata,modes=all_modes,technical=technical)


def source_summary(sections, mode, null, radius, i,j,scope):
    subset=[s for s in sections if scope=='all_HNSCC' or s['metadata']['recurrence']==scope]
    counts=Counter(); cases={}
    for s in subset:
        m=s['modes'][mode]
        status=m['count_statuses'][i,j]
        if status!='OBSERVED_STATE_COUNT_SCOPE_READY':counts[status]+=1;continue
        r=m['results'][null]
        if not r['positive_expected'][radius,i,j]:counts['NO_EXPECTED_CROSS_PARENT_STATE_EDGES']+=1;continue
        if not r['empirical_null_varies'][radius,i,j]:counts['NO_NULL_VARIATION_IN_FIXED_MONTE_CARLO']+=1;continue
        counts['ESTIMABLE_SECTION']+=1
        cases.setdefault(s['metadata']['case'],[]).append((
            r['effect'][radius,i,j],r['null_effect'][:,radius,i,j],
            r['log2_enrichment'][radius,i,j]))
    result=dict(source='GSE300147',scope=scope,minimum_human_gene_molecules=MODES[mode][0],
        regulatory_parent=MODES[mode][1],null=NULLS[null],radius_um=RADII[radius],
        R_observed_state=context.R_NAMES[i],E_observed_state=context.E_NAMES[j],
        source_cases=len({s['metadata']['case'] for s in subset}),
        supporting_cases=len(cases),section_status_counts=dict(counts))
    if len(cases)<5:
        result['status']='INSUFFICIENT_SOURCE_CASES_FOR_FROZEN_SPATIAL_SUMMARY'
        return result
    ordered=[cases[k] for k in sorted(cases)]
    effects=np.array([np.mean([s[0] for s in c]) for c in ordered])
    nulls=np.column_stack([np.mean([s[1] for s in c],axis=0) for c in ordered])
    logs=np.array([np.mean([s[2] for s in c]) for c in ordered])
    point=float(effects.mean()); pooled_null=nulls.mean(axis=1)
    # Under this fixed abundance/stratum null, the analytic excess expectation is zero.
    p=float((1+np.sum(np.abs(pooled_null)>=abs(point)))/(B+1))
    seed=[20261021,mode,null,radius,i,j, ['all_HNSCC','primary','recurrent'].index(scope)]
    boot=np.random.default_rng(seed).integers(0,len(effects),size=(B,len(effects)))
    ci=np.quantile(effects[boot].mean(axis=1),[.025,.975])
    log_ci=np.quantile(logs[boot].mean(axis=1),[.025,.975])
    loo=np.array([np.delete(effects,k).mean() for k in range(len(effects))])
    result.update(status='ESTIMATED_SOURCE_STRUCTURAL_OBSERVED_RNA_PATTERN',
        conditional_excess=point,CI95_case_bootstrap=ci.tolist(),
        log2_enrichment_equal_case_mean=float(logs.mean()),log2_CI95_case_bootstrap=log_ci.tolist(),
        permutation_P=p,positive_case_effects=int(np.sum(effects>0)),negative_case_effects=int(np.sum(effects<0)),
        zero_case_effects=int(np.sum(effects==0)),leave_one_case_range=[float(loo.min()),float(loo.max())],
        leave_one_all_same_sign_as_full=bool(np.all(np.sign(loo)==np.sign(point))),
        mean_of_independent_conditional_null=float(pooled_null.mean()),
        SE_of_conditional_null_mean=float(pooled_null.std(ddof=1)/np.sqrt(B)),
        no_clinical_or_functional_interpretation=True)
    return result


def fixed_family_bh(rows):
    """Nine prespecified slots; missing tests stay missing, family size stays nine."""
    main=[r for r in rows if r['scope']=='all_HNSCC' and r['radius_um']==50
          and r['minimum_human_gene_molecules']==20 and r['regulatory_parent']==MODES[0][1]
          and r['null']=='local_matched']
    assert len(main)==9
    valid=sorted([r for r in main if 'permutation_P' in r],key=lambda r:r['permutation_P'])
    running=1.
    for k in range(len(valid)-1,-1,-1):
        running=min(running,9*valid[k]['permutation_P']/(k+1))
        valid[k]['BH_nine_primary_opportunities']=min(1.,running)
    for r in main:r['primary_family_size']=9


def main():
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']); assert out.parent.parent==Path('/dev/shm') and not out.exists();out.mkdir(mode=0o700)
    save=lambda n,v:(out/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    question=Path(os.environ['PANCANCER_SPATIAL_QUESTION'])
    assert sha(question)==os.environ['PANCANCER_SPATIAL_QUESTION_SHA256']
    parameters=json.loads(question.read_text())
    assert parameters['radii_um']==RADII and parameters['primary_radius_um']==50
    assert parameters['nulls']==NULLS and parameters['primary_null']=='local_matched'
    assert parameters['primary_QC']==20 and parameters['primary_parent']==MODES[0][1]
    assert parameters['source_summary_min_cases']==5 and parameters['bootstraps']==B
    assert parameters['seed']==20261021
    support_root=Path(os.environ['PANCANCER_OBSERVED_SUPPORT']);accepted=context.verified(support_root)
    assert accepted['all18_source_slots_and36_source_opportunities_preserved']
    assert accepted['context_source_SHA256']==sha(context.__file__)
    plan=dict(started_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),
        context_SHA256=sha(context.__file__),kernel_SHA256=sha(kernel.__file__),
        question_SHA256=sha(question),observed_support_terminal_SHA256=sha(support_root/'private_terminal.json'),
        parameters=parameters,no_clinical_outcomes=True)
    save('plan_before_observed_state_spatial_effects.json',plan)
    save('private_synthetic_spatial_fixture.json',kernel.synthetic_fixture())
    supports=json.loads((support_root/'private_all18_source_slot_marker_support.json').read_text())
    ctx=context.load_context();rows=[]
    with ThreadPoolExecutor(max_workers=2) as pool:
        work={pool.submit(one_section,ctx,k,supports[k]):k for k in range(18)}
        for f in as_completed(work):
            rows.append(f.result())
            save('private_checkpoint.json',dict(completed_source_slots=len(rows)))
    sections=sorted([s for s in rows if not s['excluded']],key=lambda s:s['ordinal'])
    assert len(rows)==18 and len(sections)==17 and len({s['metadata']['case'] for s in sections})==10
    summaries=[source_summary(sections,mi,ni,ri,i,j,scope)
        for mi,ni,ri,i,j,scope in product(range(4),range(2),range(3),range(3),range(3),['all_HNSCC','primary','recurrent'])]
    assert len(summaries)==648
    fixed_family_bh(summaries)
    save('private_all648_source_spatial_opportunities.json',summaries)
    save('private_all17_slot_technical_checks.json',[dict(ordinal=s['ordinal'],technical=s['technical']) for s in sections])
    save('private_terminal.json',dict(status='COMPLETE_ALL648_SOURCE_OBSERVED_STATE_SPATIAL_OPPORTUNITIES',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),
        context_SHA256=sha(context.__file__),kernel_SHA256=sha(kernel.__file__),
        all18_source_slots_preserved=True,all_observed_counts_match_prior_support=True,
        all_graph_observed_expectations_and_first8_permutations_independent_edge_checked=True,
        no_payload_cells_or_case_rows_saved=True,no_clinical_validation=True,
        files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status='COMPLETE_ALL648_SOURCE_OBSERVED_STATE_SPATIAL_OPPORTUNITIES',terminal_SHA256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

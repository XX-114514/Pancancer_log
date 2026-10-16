"""Independent cell-code totals and case-set intersections; semantic roles stay candidates."""
import itertools
import json
from collections import defaultdict
from datetime import datetime,timezone

import numpy as np
import pandas as pd

from NSCLC_external_lineage_inputs_v2 import load,cell_arms,STRATA,matched_quota,sha,VIEWS,CHANNELS,PARENTS


def main():
    data=load();root=data['root'];src=root/'NSCLC_native_MFA_external_lineage_qualification_20261016_v2'
    terminal=json.loads((src/'private_terminal.json').read_text())
    for name,item in terminal['files'].items():assert sha(src/name)==item['sha256']
    proof=src/'private_independent_lineage_case_audit.json';assert not proof.exists();matched_quota(root)
    with np.load(src/'private_all1690_three_channel_role_counts.npz',allow_pickle=False) as z:observed=z['counts'];parent_observed=z['parent_counts']
    n=len(data['lib']);expected=np.zeros((n,3,4),np.int64);parents=np.zeros((n,4),np.int64);total=0
    for arm in cell_arms(data):
        idx=arm['library_index'];total+=len(idx)
        for ci,roles in enumerate(arm['roles']):
            valid=roles>=0;codes=idx[valid]*4+roles[valid];unique,counts=np.unique(codes,return_counts=True)
            expected[unique//4,ci,unique%4]=counts
        for j,view in enumerate(VIEWS):
            use=np.isin(arm['major'],list(PARENTS[view]))[arm['node'].codes]
            unique,counts=np.unique(idx[use],return_counts=True);parents[unique,j]=counts
        record=next(x for x in terminal['arm_input_proofs'] if x['arm']==arm['arm'])
        assert record['author_values_sha256']==arm['author_values_sha256'] and record['cells']==arm['cell_count']
    assert np.array_equal(observed,expected) and np.array_equal(parents,parent_observed)
    lib=pd.read_csv(src/'private_all1690_library_metadata_and_label_diagnostics.tsv',sep='\t',keep_default_na=False)
    pd.testing.assert_frame_equal(lib[data['lib'].columns],data['lib'],check_dtype=False)
    key={(r.arm,r.dataset,r.library_id):i for i,r in enumerate(lib.itertuples())}
    long=pd.read_csv(src/'private_all20280_library_role_capture.tsv.gz',sep='\t')
    assert len(long)==20280
    for r in long.itertuples():
        i=key[(r.arm,r.dataset,r.library_id)];ci=CHANNELS.index(r.channel);j=VIEWS.index(r.view)
        assert r.n_cells==expected[i,ci,j] and r.parent_cells==parents[i,j]
        assert bool(r.captured)==bool(expected[i,ci,j]>=20) and bool(r.parent_captured)==bool(parents[i,j]>=20)
        assert bool(r.parent_present_role_unresolved)==bool(parents[i,j]>=20 and expected[i,ci,j]<20)
    scopes={v:[j] for j,v in enumerate(VIEWS)}
    scopes.update({'pair::'+VIEWS[a]+'::'+VIEWS[b]:[a,b] for a,b in itertools.combinations(range(4),2)})
    scopes['all4']=list(range(4))
    source_groups={k:set(g.index) for k,g in lib.groupby(['arm','dataset'],sort=True)}
    strata_groups={k:set(g.index) for k,g in lib.groupby(['arm','dataset']+STRATA,sort=True,dropna=False)}
    def check(r,universe):
        ci=CHANNELS.index(r.channel);ends=scopes[r.capture_scope_id]
        eligible={i for i in universe if all(expected[i,ci,j]>=20 for j in ends)}
        reviewed={i for i in eligible if lib.iloc[i].reviewed_patient};people={lib.iloc[i].patient_block for i in reviewed}
        assert r.captured_libraries==len(eligible) and r.reviewed_libraries==len(reviewed) and r.reviewed_patient_blocks==len(people)
        assert r.unreviewed_library_fallback==len(eligible)-len(reviewed) and r.repeat_reviewed_library_rows==len(reviewed)-len(people)
        reviewed_endpoint=0
        for i in eligible:
            text=str(lib.iloc[i].clinical_endpoint_review).casefold()
            reviewed_endpoint+=any(s in text for s in ['reviewed','author_definition','official']) and not any(s in text for s in ['unreviewed','not_reviewed','not reviewed','pending','unknown','not_new_endpoint'])
        assert r.clinical_endpoint_reviewed_libraries==reviewed_endpoint
        phasegroups=defaultdict(set)
        for i in reviewed:
            m=lib.iloc[i];kk=tuple(m[k] for k in ['patient_block','cancer_label','tissue','material_class','treatment_arm','capture_scope','ex_vivo'])
            phasegroups[kk].add(m.phase_class)
        assert r.comparable_pre_post_patient_groups==sum({'pre','post'}<=v for v in phasegroups.values())
        assert r.comparable_pre_on_treatment_patient_groups==sum({'pre','on_treatment'}<=v for v in phasegroups.values())
        return eligible
    source=pd.read_csv(src/'private_all2013_source_role_pair_capture.tsv',sep='\t',keep_default_na=False)
    strata=pd.read_csv(src/'private_all_source_stratum_role_pair_capture.tsv.gz',sep='\t',keep_default_na=False)
    assert len(source)==2013 and not source.duplicated(['arm','dataset','channel','capture_scope_id']).any()
    for r in source.itertuples():
        check(r,source_groups[(r.arm,r.dataset)]);assert r.all_source_libraries==len(source_groups[(r.arm,r.dataset)])
        assert not r.exact_native_identity_accepted and not r.independent_validation_accepted
    for r in strata.itertuples():
        universe=strata_groups[(r.arm,r.dataset,*[getattr(r,k) for k in STRATA])]
        check(r,universe);assert r.all_stratum_libraries==len(universe)
    genes=pd.read_csv(root/'NSCLC_native_MFA_external_gene_qualification_20261016_v1/private_all2440_native_factor_source_gene_coverage.tsv',sep='\t').set_index(['arm','dataset','view','factor'])
    opportunities=pd.read_csv(src/'private_all20130_factor_source_capture_and_gene_opportunities.tsv.gz',sep='\t')
    ss=source.set_index(['arm','dataset','channel','capture_scope_id']);assert len(opportunities)==20130
    for r in opportunities.itertuples():
        s=ss.loc[(r.arm,r.dataset,r.channel,r.capture_scope_id)];ends=scopes[r.capture_scope_id]
        g=all(bool(genes.loc[(r.arm,r.dataset,VIEWS[j],r.factor),'physical_coverage_pass']) for j in ends)
        assert r.all_endpoints_gene_gate==g and r.captured_libraries==s.captured_libraries and r.reviewed_patient_blocks==s.reviewed_patient_blocks
        assert bool(r.gene_and_capture_library_opportunity)==bool(g and s.captured_libraries>0)
    assert source.captured_libraries.eq(0).any() and long.parent_present_role_unresolved.any()
    report=dict(status='PASS_ALL61_SOURCES1690_LIBRARIES_THREE_ROLE_CHANNELS_ALL_CAPTURE_PATIENT_PHASE_AND20130_FACTOR_OPPORTUNITIES',finished_UTC=datetime.now(timezone.utc).isoformat(),
        audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(src/'private_terminal.json'),fresh_metadata_cells_checked=total,
        library_role_rows_checked=len(long),source_rows_checked=len(source),stratum_rows_checked=len(strata),factor_opportunities_checked=len(opportunities),
        semantic_rules_are_candidates_not_independent_biological_validation=True,new_expression_or_P=False,identity_or_patient_independence_promoted=False)
    proof.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],audit_sha256=sha(proof))))


if __name__=='__main__':main()

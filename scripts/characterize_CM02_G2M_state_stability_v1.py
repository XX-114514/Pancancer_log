"""Finite descriptive RNA follow-up, closed primary family remains unchanged.

Exact original references and native counts are re-read with all individual
arrays RAM-only. Writes only aggregate fine-label coverage, LOO ranges, and
fixed odd/even score sensitivities. No new P, CI, clinical fit or identities.
"""
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import resource
import traceback
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
import analyze_CM02_RNA_multigene_in_RAM_v1 as baseline
import CM02_RNA_state_stability_helpers_v1 as helpers


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_all_scores(z,selected,human,measured,programmes,checkpoint):
    raw=z['layers/counts'];shape=tuple(raw.attrs['shape']);ptr=raw['indptr'][:]
    assert len(ptr)==shape[0]+1 and ptr[0]==0 and np.all(np.diff(ptr)>=0)
    assert ptr[-1]==raw['data'].shape[0]==raw['indices'].shape[0] and np.all(np.diff(selected)>0)
    scores=np.zeros((len(selected),len(programmes)));total=np.zeros(len(selected),np.int64);detected=total.copy();completed=0;replayed=0
    for start in range(0,shape[0],5000):
        stop=min(start+5000,shape[0]);lo,hi=np.searchsorted(selected,[start,stop])
        if lo==hi:continue
        a,b=int(ptr[start]),int(ptr[stop]);indices=raw['indices'][a:b];values=raw['data'][a:b];localptr=ptr[start:stop+1]-a
        rows=selected[lo:hi]-start
        ss,tt,dd=baseline.scoring.score_selected_CSR_rows(localptr,indices,values,rows,human,measured,programmes)
        for position in sorted({0,len(rows)-1}):
            row=int(rows[position]);begin,end=map(int,localptr[row:row+2]);dense=np.zeros(shape[1],np.int64)
            np.add.at(dense,np.asarray(indices[begin:end],int),np.asarray(values[begin:end],np.int64))
            ranks=np.full(shape[1],np.nan);ranks[human]=pd.Series(-dense[human]).rank(method='average').to_numpy()
            expected=[]
            for columns in programmes:
                n=len(columns);minimum=n*(n+1)/2;maximum=min(1500,int(human.sum()))
                expected.append(1-(np.minimum(ranks[columns],maximum).sum()-minimum)/(n*maximum-minimum))
            np.testing.assert_allclose(ss[position],expected,atol=1e-12,rtol=1e-12)
            assert dense[human].sum()==tt[position] and np.count_nonzero(dense[human])==dd[position];replayed+=1
        scores[lo:hi]=ss;total[lo:hi]=tt;detected[lo:hi]=dd;completed+=hi-lo
        checkpoint(dict(stage='native_four_programme_scores',completed_source_cells=int(completed),independent_dense_score_rows=replayed))
        assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024<=8*1024**3
    assert completed==len(selected)
    return scores,total,detected,dict(source_cells_scored=int(completed),independent_dense_boundary_rows=replayed,
        all_selected_rows_canonicalized_before_rank=True,no_rank_or_individual_array_written=True)


def main():
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert out.parent.parent==Path('/dev/shm') and not out.exists();out.mkdir(mode=0o700)
    def save(name,value):(out/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    try:
        root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1'
        previous=Path(os.environ['PANCANCER_BASELINE_OUTPUT']);terminal=previous/'private_terminal.json'
        assert sha(terminal)==os.environ['PANCANCER_BASELINE_TERMINAL_SHA256']
        accepted=json.loads(terminal.read_text());assert accepted['status']=='COMPLETE_FIXED32_NATIVE_RNA_PROGRAMME_ASSOCIATIONS'
        for name,expected in accepted['files'].items():assert sha(previous/name)==expected['SHA256']
        previous_plan=json.loads((previous/'plan_before_native_RNA_effects.json').read_text())
        for name,expected in previous_plan['script_SHA256'].items():assert sha(Path(__file__).parent/name)==expected
        original=json.loads((previous/'private_all32_RNA_programme_associations.json').read_text())['results']
        signature=Path(os.environ['PANCANCER_PROGRAMME_DEFINITION']);assert sha(signature)==os.environ['PANCANCER_PROGRAMME_SHA256']
        definition=json.loads(signature.read_text());programmes=helpers.split_programmes(definition['programs'])
        save('plan_before_state_stability.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source='GSE286827',
            scope='2026-10-14_CM02_v1_T_whole_TME_primary_tumor',previous_terminal_SHA256=sha(terminal),
            code_SHA256={p.name:sha(p) for p in Path(__file__).parent.glob('*.py')},external_programme_SHA256=sha(signature),
            new_programmes={name:genes for name,genes in programmes.items() if name.startswith('g2m.') and name!='g2m.genes'},
            source_coverage_minimum=.8,minimum_parent=20,minimum_cases=20,minimum_residual_df=10,maximum_condition=1e8,
            all32_previous_slots_LOO=True,LOO_original_categorical_levels_fixed=True,no_case_removed_from_primary=True,
            split_sensitivity='all8 original G2M/G2M scopes x4 odd/even endpoint combinations; descriptive rho only',
            fine_profiles='all13 fixed author fine labels xpre/post x2 original contexts; >=20 cells per case, distributions only>=5 cases',
            post2_coverage_only=True,no_new_P_CI_or_clinical_family=True,
            state_threshold_or_new_classifier_not_fitted=True,no_individual_or_learned_arrays_written=True))
        fine=Path(os.environ['PANCANCER_FINE_LOADER']);ft=Path(os.environ['PANCANCER_FINE_SCOPE_OUTPUT'])/'private_terminal.json'
        loader=baseline.load(fine,json.loads(ft.read_text())['source_SHA256'],'accepted_fine_reference')
        eligibility=baseline.load(os.environ['PANCANCER_ELIGIBILITY_SOURCE'],'a05c1d06b3c873b8abc0e78da20e420b20a0768e241487a43b17439fbfa2d3bd','accepted_RNA_eligibility')
        network=baseline.load(os.environ['PANCANCER_NETWORK_RESOLVER'],os.environ['PANCANCER_NETWORK_RESOLVER_SHA256'],'accepted_CM02_network')
        tab,cur,unique,shared,matrices,table=loader.load_inputs()
        contexts,membership=network.resolve(atlas,['all','regulatory_anchor','cytotoxic_anchor'],['39697d077b6cb8bf','c3f761a9da61f578']);nodes=contexts['CM02_network_member_union']
        input_source=Path(os.environ['PANCANCER_RNA_INPUT_OUTPUT'])
        assert sha(input_source/'private_terminal.json')==os.environ['PANCANCER_RNA_INPUT_TERMINAL_SHA256']
        original_input=json.loads((input_source/'private_terminal.json').read_text())
        for name,expected in original_input['files'].items():assert sha(input_source/name)==expected['SHA256']
        assert eligibility.reference_coverage(tab,unique,nodes)==json.loads((input_source/'private_author_parent_phase_coverage.json').read_text())
        assert membership==json.loads((input_source/'private_original_network_context_qualification.json').read_text())
        z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r')
        selected=np.flatnonzero(np.asarray(read_elem(z['obs/dataset'])).astype(str)=='GSE286827')
        assert len(selected)==len(cur) and np.array_equal(cur.position,np.arange(len(cur)))
        var=read_elem(z['var']);source_order=list(map(str,read_elem(z['uns/source_order'])))
        measured=np.asarray(read_elem(z['varm/measured_in_source']),bool)[:,source_order.index('GSE286827')]
        human=measured&var.feature_identity_status.astype(str).eq('stable_Ensembl_ID').to_numpy()
        coverage=eligibility.feature_coverage(var,human,programmes);assert all(r['source_coverage_accepted'] for r in coverage)
        columns=baseline.programme_columns(var,human,definition['programs'])
        for name in ['g2m.odd_positions','g2m.even_positions']:
            columns.append(baseline.programme_columns(var,human,{'s.genes':programmes['s.genes'],'g2m.genes':programmes[name]})[1])
        assert set(columns[2]).isdisjoint(columns[3]) and set(columns[2])|set(columns[3])==set(columns[1])
        unique=unique.loc[unique.major.eq('T')].sort_values('position').reset_index(drop=True)
        scores,total,detected,score_receipt=read_all_scores(z,selected[unique.position.to_numpy(int)],human,measured,columns,
            lambda v:save('private_progress.json',dict(observed_UTC=datetime.now(timezone.utc).isoformat(),**v)))
        previous_audit=json.loads((previous/'private_RNA_input_and_score_audit.json').read_text())
        assert score_receipt==previous_audit['score_audit']
        units,case_coverage=baseline.inputs.aggregate(tab,cur,unique,nodes,scores[:,:2],total,detected)
        assert case_coverage==previous_audit['case_coverage']
        loo=[]
        for slot,definition_slot in enumerate(baseline.statistics.fixed_slots()):
            x,y,c,cat,cases=baseline.inputs.pair_arrays(units,definition_slot)
            result=helpers.describe_fit(x,y,c,cat,baseline.statistics.categorical_spec(cat));old=original[slot]
            assert old['n_cases']==len(cases) and old['status']=='ESTIMABLE' and result['status']=='DESCRIPTIVE_ESTIMABLE'
            np.testing.assert_allclose([result['rho'],result['unadjusted_rho'],result['condition']],
                [old['rho'],old['unadjusted_rho'],old['design_condition']],atol=1e-11,rtol=1e-11)
            loo.append(dict(slot=slot+1,**definition_slot,original_rho=old['rho'],**helpers.leave_one_case(x,y,c,cat)))
        save('private_all32_LOO_descriptions.json',dict(source='GSE286827',descriptions=loo,all32_previous_observed_scalars_replayed=True,
            no_new_P_CI_or_modified_primary_results=True))
        fine_rows=helpers.fine_profiles(tab,cur,unique,nodes,scores)
        save('private_all52_author_fine_profiles.json',dict(source='GSE286827',profiles=fine_rows,
            post2_only_reviewed_case_coverage=int(tab.loc[tab.timepoint_registry.eq('post2')].patient_registry.nunique()),
            post2_expression_not_described=True,source_coverage=coverage))
        half_units,half_coverage=baseline.inputs.aggregate(tab,cur,unique,nodes,scores[:,2:],total,detected)
        assert half_coverage==case_coverage
        halves=[]
        for definition_slot in baseline.statistics.fixed_slots():
            if definition_slot['left_programme']!=definition_slot['right_programme'] or definition_slot['left_programme']!='g2m.genes':continue
            for left,labelleft in [('s.genes','g2m.odd_positions'),('g2m.genes','g2m.even_positions')]:
                for right,labelright in [('s.genes','g2m.odd_positions'),('g2m.genes','g2m.even_positions')]:
                    slot=dict(definition_slot,left_programme=left,right_programme=right)
                    x,y,c,cat,cases=baseline.inputs.pair_arrays(half_units,slot)
                    result=helpers.describe_fit(x,y,c,cat,baseline.statistics.categorical_spec(cat))
                    halves.append(dict(**dict(definition_slot,left_programme=labelleft,right_programme=labelright),n_cases=len(cases),**result))
        assert len(halves)==32
        save('private_all32_split_programme_descriptions.json',dict(source='GSE286827',descriptions=halves,
            no_new_P_CI_or_independent_validation=True,posthoc_finite_measurement_sensitivity=True))
        save('private_score_and_scope_audit.json',dict(source='GSE286827',score_audit=score_receipt,
            previous_observed32_scalars_and_case_coverage_replayed=True,source_programme_coverage=coverage,
            original_network_and_reference_keys_verified=True,no_new_individual_or_learned_arrays_written=True))
        assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024<=8*1024**3
        save('private_terminal.json',dict(status='COMPLETE_DESCRIPTIVE_FINE_LOO_AND_SPLIT_MEASUREMENT',completed_UTC=datetime.now(timezone.utc).isoformat(),
            source_SHA256=sha(__file__),peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
            no_new_P_CI_or_clinical_tests=True,no_individual_or_learned_arrays_written=True,
            files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()}))
        print(json.dumps(dict(status='COMPLETE_DESCRIPTIVE_FINE_LOO_AND_SPLIT_MEASUREMENT',terminal_SHA256=sha(out/'private_terminal.json'))))
    except Exception as error:
        save('private_failed_terminal.json',dict(status='DESCRIPTIVE_STATE_STABILITY_NOT_ACCEPTED',observed_UTC=datetime.now(timezone.utc).isoformat(),error_class=type(error).__name__,
            frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)],no_full_error_or_individual_values_written=True))
        raise SystemExit(1)


if __name__=='__main__':main()

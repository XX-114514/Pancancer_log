"""Outcome-free source eligibility for external RNA programmes.

Existing reviewed patient/library and exact author-reference joins are reused.
Only source-level coverage is saved; no cell rows, rank arrays, expression or
learned model arrays are read out or persisted. This step does not score cells.
"""
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import os
import resource
import traceback
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem

SOURCE='GSE286827'
CONTEXTS=['current_T_author_reference','CM02_network_member_union_author_reference']
PARENTS=['author_regulatory_reference','author_conventional_CD4_reference','author_CD8_reference']
REGULATORY=['4-1BB+ Treg','HLA-DR+ Treg','LAG3+ Treg','Prolif Treg']
CONVENTIONAL=['CXCL13+ CD4T','IFN-stim CD4T','Tn/Tmem']
CD8=['CXCL13+Tex','Tdys','GZMK+T','Tmem','Tprolif','Tstr']


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path,expected,name):
    assert sha(path)==expected
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def feature_coverage(var,human,programs):
    lookup=defaultdict(set)
    for column in np.flatnonzero(human):
        row=var.iloc[column]
        names={str(var.index[column]),str(row.gene_symbol)}|set(str(row.source_gene_aliases).split('|'))
        for name in names:lookup[name].add(int(column))
    results=[]
    for name,genes in programs.items():
        assert len(genes)==len(set(genes)) and not {'MKI67','TOP2A','UBE2C'}.intersection(genes)
        mapped={gene:next(iter(lookup[gene])) for gene in genes if len(lookup[gene])==1}
        per_column=Counter(mapped.values())
        duplicate={gene for gene,column in mapped.items() if per_column[column]>1}
        accepted={gene:column for gene,column in mapped.items() if gene not in duplicate}
        missing=[gene for gene in genes if not lookup[gene]]
        ambiguous=[gene for gene in genes if len(lookup[gene])>1]
        results.append(dict(program=name,requested_external_genes=len(genes),uniquely_mapped_measured_genes=len(accepted),
            missing_gene_names=missing,ambiguous_gene_names=ambiguous,same_feature_multiple_gene_names=sorted(duplicate),
            source_coverage_fraction=len(accepted)/len(genes),minimum_fraction=.8,
            source_coverage_accepted=len(accepted)>=math.ceil(.8*len(genes)),
            missing_not_imputed_to_zero=True))
    return results


def reference_coverage(tab,unique,network_nodes):
    registry=tab.set_index('library_id')
    frame=unique.copy()
    assert frame.position.is_unique and frame.library_id.isin(registry.index).all()
    frame['reviewed_case']=frame.library_id.map(registry.patient_registry)
    expected=frame.library_id.map(registry.timepoint_registry)
    assert expected.eq(frame.timepoint).all() and frame.reviewed_case.notna().all()
    def parent(label):
        if label in REGULATORY:return PARENTS[0]
        if label in CONVENTIONAL:return PARENTS[1]
        assert label in CD8;return PARENTS[2]
    frame['reference_parent']=frame.label.map(parent)
    out=[];pairs=[]
    for context in CONTEXTS:
        selected=frame.loc[frame.major.eq('T')]
        if context==CONTEXTS[1]:selected=selected.loc[selected.node_id.isin(network_nodes)]
        counts=selected.groupby(['reviewed_case','timepoint','reference_parent']).size()
        for phase in ['pre','post','post2']:
            phase_tab=tab.loc[tab.timepoint_registry.eq(phase)]
            cases=sorted(phase_tab.patient_registry.unique())
            matrix=np.zeros((len(cases),3),dtype=int)
            for i,case in enumerate(cases):
                for j,role in enumerate(PARENTS):matrix[i,j]=counts.get((case,phase,role),0)
            for j,role in enumerate(PARENTS):
                out.append(dict(source=SOURCE,phase=phase,context=context,parent=role,
                    reviewed_cases=len(cases),source_cell_observations=int(matrix[:,j].sum()),
                    cases_with_reference_cells=int(np.sum(matrix[:,j]>0)),cases_parent_atleast20=int(np.sum(matrix[:,j]>=20)),
                    reference_absence_not_asserted_biological_absence=True))
            for left in range(2):
                eligible=(matrix[:,left]>=20)&(matrix[:,2]>=20)
                arms=[];hpvs=[];complete=0
                for case in np.asarray(cases)[eligible]:
                    rows=phase_tab.loc[phase_tab.patient_registry.eq(case)]
                    arm=set(rows.neo_type);hpv=set(rows['HPV/p16'])
                    if len(arm)==1 and arm<=set(['D','D+T']) and len(hpv)==1 and hpv<=set(['positive','negative']):
                        complete+=1;arms.extend(arm);hpvs.extend(hpv)
                pairs.append(dict(source=SOURCE,phase=phase,context=context,left_parent=PARENTS[left],right_parent=PARENTS[2],
                    both_parents_atleast20_reviewed_cases=int(eligible.sum()),
                    cases_with_concordant_known_arm_and_HPV=complete,
                    arm_category_counts=dict(Counter(arms)),HPV_category_counts=dict(Counter(hpvs)),
                    future_minimum20_case_coverage_only=complete>=20,
                    expression_variation_and_covariate_rank_not_yet_qualified=True))
    other=frame.loc[~frame.major.eq('T')].groupby(['reference_parent','major']).size()
    phase_summary=[dict(phase=phase,reviewed_libraries=len(g),reviewed_cases=g.patient_registry.nunique(),
        duplicate_patient_phase_libraries=int(g.patient_registry.duplicated(keep=False).sum()))
        for phase,g in tab.groupby('timepoint_registry',sort=True)]
    return dict(parent_coverage=out,pair_coverage=pairs,reviewed_phase_scope=phase_summary,
        author_reference_current_nonT_exclusions=[dict(parent=a,current_major=b,source_cells=int(n)) for (a,b),n in other.items()],
        primary_scope_current_major_T_and_unique_author_reference=True,
        fine_reference_subset_not_whole_parent_truth=True,no_new_annotation_or_patient_identity_assignment=True)


def main():
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert out.parent.parent==Path('/dev/shm') and not out.exists();out.mkdir(mode=0o700)
    def save(name,value):(out/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    try:
        root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1'
        signature=Path(os.environ['PANCANCER_PROGRAMME_DEFINITION']);assert sha(signature)==os.environ['PANCANCER_PROGRAMME_SHA256']
        definition=json.loads(signature.read_text());assert definition['excluded_prior_marker_genes']==['MKI67','TOP2A','UBE2C']
        fine=Path(os.environ['PANCANCER_FINE_LOADER']);terminal=Path(os.environ['PANCANCER_FINE_SCOPE_OUTPUT'])/'private_terminal.json'
        accepted=json.loads(terminal.read_text());loader=load(fine,accepted['source_SHA256'],'accepted_author_reference_inputs')
        network_path=Path(os.environ['PANCANCER_NETWORK_RESOLVER']);network=load(network_path,os.environ['PANCANCER_NETWORK_RESOLVER_SHA256'],'accepted_original_network')
        save('plan_before_RNA_programme_qualification.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),
            source_SHA256=sha(__file__),source=SOURCE,scope='2026-10-14_CM02_v1_T_whole_TME_primary_tumor',
            public_programme_definition_SHA256=sha(signature),prior_fine_scope_terminal_SHA256=sha(terminal),
            fine_loader_SHA256=sha(fine),regulatory_labels=REGULATORY,conventional_labels=CONVENTIONAL,CD8_labels=CD8,
            current_major_T_required=True,minimum_parent=20,minimum_future_case_coverage=20,
            minimum_measured_unique_signature_fraction=.8,missing_not_filled_negative=True,
            no_expression_or_outcome_effects_computed=True,no_individual_or_learned_arrays_written=True))
        tab,frame,unique,shared,matrices,table=loader.load_inputs()
        contexts,membership=network.resolve(atlas,['all','regulatory_anchor','cytotoxic_anchor'],['39697d077b6cb8bf','c3f761a9da61f578'])
        nodes=contexts['CM02_network_member_union']
        z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r')
        var=read_elem(z['var']);source_order=list(map(str,read_elem(z['uns/source_order'])))
        measured=np.asarray(read_elem(z['varm/measured_in_source']),bool)[:,source_order.index(SOURCE)]
        human=measured & var.feature_identity_status.astype(str).eq('stable_Ensembl_ID').to_numpy()
        gene=feature_coverage(var,human,definition['programs'])
        coverage=reference_coverage(tab,unique,nodes)
        save('private_source_gene_programme_qualification.json',dict(source=SOURCE,
            measured_stable_human_features=int(human.sum()),programmes=gene))
        save('private_author_parent_phase_coverage.json',coverage)
        save('private_original_network_context_qualification.json',membership)
        assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024<=8*1024**3
        save('private_terminal.json',dict(status='COMPLETE_OUTCOME_FREE_RNA_PROGRAMME_INPUT_QUALIFICATION',
            completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),
            all_programme_coverage_gates_accepted=all(x['source_coverage_accepted'] for x in gene),
            exact_prior_reference_key_and_crosscount_checks_replayed=True,
            unique_author_reference_overlap_excluded=True,
            source_count_layer_not_read_no_expression_or_P_or_clinical_effects=True,
            no_new_individual_or_learned_arrays_written=True,
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
            files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()}))
        print(json.dumps(dict(status='COMPLETE_OUTCOME_FREE_RNA_PROGRAMME_INPUT_QUALIFICATION',terminal_SHA256=sha(out/'private_terminal.json'))))
    except Exception as error:
        save('private_failed_terminal.json',dict(status='RNA_PROGRAMME_INPUT_QUALIFICATION_NOT_ACCEPTED',
            observed_UTC=datetime.now(timezone.utc).isoformat(),error_class=type(error).__name__,
            frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)],
            individual_values_and_full_errors_not_returned=True))
        raise SystemExit(1)


if __name__=='__main__':main()

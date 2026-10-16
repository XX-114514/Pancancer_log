"""Join capture and physical coverage with all original programme/source roles."""
import hashlib,itertools,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd

STRATA=['cancer_label','tissue','sample_class','timepoint','treatment_arm','capture_scope','ex_vivo','official_stage']
ATLAS_MAJOR={'B cell':'B','T cell':'T','Macrophage':'Myeloid','Epithelial cell':'Epithelial','Fibroblast':'Fibroblast','Endothelial cell':'Endothelial','Mural cell':'Pericyte_SMC'}


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);capture=r/'complete_external_lineage_opportunities_20261016_v2';genes=r/'complete_external_functional_gene_coverage_20261016_v1'
    proofs=[capture/'private_independent_opportunity_audit.json',genes/'private_independent_gene_coverage_audit.json']
    assert all(json.loads(p.read_text())['status'].startswith('PASS_ALL') for p in proofs)
    paths={'libraries':capture/'private_all_library_latest_metadata.tsv','counts':capture/'private_all_original_node_library_counts.tsv.gz','pool_counts':capture/'private_all_library_pool_counts.tsv.gz','pool_definitions':capture/'private_all_pool_node_definitions.tsv','programme_coverage':genes/'private_all3416_programme_native_coverage.tsv','functional_coverage':genes/'private_all48800_functional_native_coverage.tsv','original_units':r/'all_CM_state_native_input_preflight_20261014_v1/private_all_frozen_sample_node_units.tsv.gz','register':r/'state_discovery_evidence_register_20261014_v1.json'}
    for directory in [capture,genes]:
        terminal=json.loads((directory/'private_terminal.json').read_text())
        for name,digest in terminal['files'].items():assert sha(directory/name)==digest
    lib=pd.read_csv(paths['libraries'],sep='\t',keep_default_na=False,dtype={'dataset':str,'library_id':str,'patient_block':str});units=pd.read_csv(paths['counts'],sep='\t',dtype={'dataset':str,'library_id':str,'node_id':str});pc=pd.read_csv(paths['pool_counts'],sep='\t',dtype={'dataset':str,'library_id':str})
    original=pd.read_csv(paths['original_units'],sep='\t',dtype={'node_id':str});labels=original[['node_id','subCluster','major']].drop_duplicates();assert labels.subCluster.is_unique and labels.node_id.is_unique
    mapping=labels.set_index('subCluster');coverage=pd.read_csv(paths['programme_coverage'],sep='\t',keep_default_na=False);functions=pd.read_csv(paths['functional_coverage'],sep='\t');slots=coverage[['context','program','type','original_endpoint']].drop_duplicates();assert len(slots)==56
    parent=Path(json.loads(paths['register'].read_text())['model_output_root'])/'state_after_DIALOGUE_inputs_20261016_v3';roles=[];role_hashes={}
    for context,group in slots.groupby('context'):
        p=parent/f'context_{context:04d}'/'private_input_inventory.json';x=json.loads(p.read_text());assert int(x['context'])==context
        ends=[group.loc[group.type.eq(t),'original_endpoint'].unique().item() for t in ['typeA','typeB']];assert ends==x['endpoints']
        roles.append(dict(context=context,training_dataset=x['gse_id'],training_stratum=x['stratum'],training_mode=x['mode'],endpoint_A=ends[0],endpoint_B=ends[1],major_A=mapping.loc[ends[0],'major'],major_B=mapping.loc[ends[1],'major']))
        role_hashes[str(context)]=sha(p)
    bycontext={x['context']:x for x in roles};assert len(roles)==14
    count_index={(x.arm,x.dataset,x.library_id,x.node_id):int(x.n_cells) for x in units.itertuples()};pool_index={(x.arm,x.dataset,x.library_id,x.pool):int(x.n_cells) for x in pc.itertuples()}
    cf=coverage.set_index(['arm','dataset','context','program','type']);functional=functions.groupby(['arm','dataset','context','program']).technical_measurement_gate.sum().to_dict()
    out=r/'complete_external_programme_readiness_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),input_sha256={k:sha(p) for k,p in paths.items()},audit_sha256={p.name:sha(p) for p in proofs},training_inventory_sha256=role_hashes,scope='Every original14 context x2 programme x61 source-arm; exact frozen172 endpoints in original arm, all inherited major extensions in both arms and prespecified multigene candidate pairs in atlas. All zero/unsupported/unknown strata retained.',limits='Broad-major and candidate transfer are not exact old172 or NSCLC identity. Source roles refer to programme training, not independent original-CM validation. No new effects, endpoint acceptance, P or downloads.')
    (out/'plan_before_source_question_selection.json').write_text(json.dumps(plan,indent=2)+'\n')
    rows=[];strata_rows=[];Hrows=[]
    for (arm,source),libraries in lib.groupby(['arm','dataset'],sort=True):
        for context,program in slots[['context','program']].drop_duplicates().itertuples(index=False,name=None):
            role=bycontext[context];a=role['endpoint_A'];b=role['endpoint_B'];ma=role['major_A'];mb=role['major_B'];specs=[]
            if arm=='frozen_original_annotation':specs.append(('exact_frozen172_nodes',str(mapping.loc[a,'node_id']),str(mapping.loc[b,'node_id'])))
            broad=(ma,mb) if arm=='frozen_original_annotation' else (ATLAS_MAJOR.get(ma),ATLAS_MAJOR.get(mb))
            assert all(broad);specs.append(('broad_inherited_major_extension','major::'+broad[0],'major::'+broad[1]))
            candidate=lambda endpoint:'candidate::CD8' if endpoint.startswith('CD8_T_cell_') else 'candidate::Treg' if endpoint.startswith('Regulatory_T_cell_') else 'candidate::Macrophage' if endpoint.startswith('Macrophage_') else None
            if arm=='archived_immune_atlas' and candidate(a) and candidate(b):specs.append(('multigene_candidate_extension',candidate(a),candidate(b)))
            ca=cf.loc[(arm,source,context,program,'typeA')];cb=cf.loc[(arm,source,context,program,'typeB')]
            for scope,pa,pb in specs:
                counts=count_index if scope=='exact_frozen172_nodes' else pool_index
                members=set();metadata={}
                for m in libraries.itertuples():
                    key=(arm,source,m.library_id);metadata[m.library_id]=m
                    if counts.get((*key,pa),0)>=20 and counts.get((*key,pb),0)>=20:members.add(m.library_id)
                static=dict(arm=arm,dataset=source,context=context,program=program,scope=scope,endpoint_A=a,endpoint_B=b,pool_A=pa,pool_B=pb,training_dataset=role['training_dataset'],same_programme_training_source=source==role['training_dataset'],original_CM_discovery_arm=arm=='frozen_original_annotation',exact_old172_identity=scope=='exact_frozen172_nodes',exact_NSCLC_identity_accepted=False,untouched_independent_cohort_accepted=False,new_effect_or_clinical_model_accepted=False,programme_A_95=bool(ca.programme95_gate),programme_B_95=bool(cb.programme95_gate),both_programmes_95=bool(ca.programme95_gate and cb.programme95_gate),native_function_slots_pass=int(functional.get((arm,source,context,program),0)),original_NSCLC_function_mask_available=(arm,source,context,program) in functional)
                def summary(selected):
                    people={metadata[k].patient_block for k in selected if metadata[k].reviewed_patient};return dict(libraries=len(selected),reviewed_patient_blocks=len(people),unreviewed_library_fallback=sum(not metadata[k].reviewed_patient for k in selected))
                rows.append(dict(**static,all_source_libraries=len(libraries),**summary(members)))
                for values,g in libraries.groupby(STRATA,dropna=False,sort=True):
                    selected=members&set(g.library_id);strata_rows.append(dict(**static,**dict(zip(STRATA,values)),all_stratum_libraries=len(g),**summary(selected)))
                if source=='GSE286827' and scope=='multigene_candidate_extension':
                    for require_B in [False,True]:
                        selected=members if not require_B else {k for k in members if pool_index.get((arm,source,k,'extension::B_plus_Plasma'),0)>=20};phases={}
                        for key in selected:
                            m=metadata[key]
                            if m.reviewed_patient and m.timepoint in ['pre','post']:phases.setdefault(m.patient_block,set()).add(m.timepoint)
                        Hrows.append(dict(**static,B_Plasma_required=require_B,reviewed_pre_post_patient_pairs=sum({'pre','post'}<=x for x in phases.values())))
    table=pd.DataFrame(rows);ss=pd.DataFrame(strata_rows);hp=pd.DataFrame(Hrows)
    assert len(table)==52*28*2+9*28+9*8 and len(hp)==16
    assert not table.duplicated(['arm','dataset','context','program','scope']).any() and (table.libraries<=table.all_source_libraries).all()
    assert (ss.libraries<=ss.all_stratum_libraries).all() and not table.new_effect_or_clinical_model_accepted.any()
    for n,t in [('private_all_source_programme_readiness.tsv',table),('private_all_source_stratum_programme_readiness.tsv.gz',ss),('private_all_HNSCC_programme_pair_capture.tsv',hp),('private_all14_programme_discovery_source_roles.tsv',pd.DataFrame(roles))]:t.to_csv(out/n,sep='\t',index=False)
    terminal=dict(status='COMPLETE_ALL_ORIGINAL_CONTEXT_PROGRAMME_SOURCE_ROLES_EXACT_AND_EXTENSION_READINESS_NO_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),source_rows=len(table),stratum_rows=len(ss),HNSCC_pair_capture_rows=len(hp),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_expression_P_or_clinical_acceptance=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()

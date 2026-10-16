"""All sources, three explicit label channels and all native view/pair captures."""
import itertools
import json
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_external_lineage_inputs_v1 import load,cell_arms,STRATA,matched_quota,sha,VIEWS,CHANNELS,PARENTS


def scopes():
    return [(v,[j]) for j,v in enumerate(VIEWS)]+[('pair::'+VIEWS[a]+'::'+VIEWS[b],[a,b]) for a,b in itertools.combinations(range(4),2)]+[('all4',list(range(4)))]


def summary(lib,take):
    g=lib.loc[take];reviewed=g[g.reviewed_patient]
    pairgroup=['patient_block','cancer_label','tissue','material_class','treatment_arm','capture_scope','ex_vivo']
    phases=reviewed.groupby(pairgroup,dropna=False).phase_class.agg(set)
    endpoint=g.clinical_endpoint_review.astype(str)
    endpoint_reviewed=endpoint.str.contains('reviewed|author_definition|official',case=False)&~endpoint.str.contains('unreviewed|not_reviewed|not reviewed|pending|unknown|not_new_endpoint',case=False)
    return dict(captured_libraries=len(g),reviewed_libraries=len(reviewed),reviewed_patient_blocks=int(reviewed.patient_block.nunique()),
        unreviewed_library_fallback=int((~g.reviewed_patient).sum()),repeat_reviewed_library_rows=int(len(reviewed)-reviewed.patient_block.nunique()),
        comparable_pre_post_patient_groups=int(phases.map(lambda s:{'pre','post'}<=s).sum()),
        comparable_pre_on_treatment_patient_groups=int(phases.map(lambda s:{'pre','on_treatment'}<=s).sum()),
        clinical_endpoint_reviewed_libraries=int(endpoint_reviewed.sum()))


def main():
    data=load();root=data['root'];lib=data['lib'];out=root/'NSCLC_native_MFA_external_lineage_qualification_20261016_v1'
    assert not out.exists();quota=matched_quota(root);out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        input_helper_sha256=sha(Path(__file__).with_name('NSCLC_external_lineage_inputs_v1.py')),
        role_rules_sha256=sha(Path(__file__).with_name('NSCLC_external_lineage_role_rules_v1.py')),input_sha256=data['hashes'],quota=quota,
        scope='All61 existing source-arm/1690 library entries, all4 native view roles, all6 pairs and all4 capture. Three label channels remain separate; all10 factors retained in gene-qualified opportunities.',
        channels=CHANNELS,frozen='Original172 current label prefixes CD8_T_cell/Regulatory_T_cell, Macrophage major and whole B-cell major; atlas270 supported multigene CD8_Effector/Treg/Macro_* candidates and whole B+Plasma. Compatible subsets, never exact native full-view or formal state bridge.',
        source_author='Use retained ev_author_label/source_author_label, conservative explicit lexical roles; mixed CD4/CD8, mono/mac, macro/DC or multiple roles unresolved. Source-author candidate, not a new accepted cell annotation.',
        concordance='Same unambiguous author roles additionally compatible with current broad parent sets; all conflicts separately counted, no automatic correction.',
        parent_diagnostic='Broad parent capture only; overlapping unresolved lymphoid parent sets are not disjoint native denominators or fractions.',
        capture='Each endpoint>=20 cells within the SAME physical library; zero subgroup counts remain zero. Distinguish parent-not-captured and parent-present-role-unresolved.',
        patients='Reviewed unique patient blocks, library rows, duplicate reviewed library rows and unknown sample substitutes separately. Pre/post and pre/on-treatment groups require same patient,cancer,tissue,material,treatment,capture and exvivo. Do not sum source-arm or stratum counts as independent cohorts/patients.',
        material='Preserve original strata and separately classify exvivo/product/blood/marrow/LN/normal/recorded solid tumor/anatomical-only/unknown; no pooled effect.',
        role_or_clinical_identity_promoted=False,new_expression_scores_or_P=False,external_independence_or_model_accepted=False)
    (out/'plan_before_complete_native_lineage_case_capture.json').write_text(json.dumps(plan,indent=2)+'\n')
    n=len(lib);values=np.zeros((n,3,4),np.int64);parent=np.zeros((n,4),np.int64)
    ambiguous=np.zeros(n,np.int64);conflict=np.zeros(n,np.int64);author_present=np.zeros(n,np.int64);vocab=[];arm_proofs=[]
    for arm in cell_arms(data):
        idx=arm['library_index'];nodes=arm['node'];authors=arm['author']
        for ci,roles in enumerate(arm['roles']):
            for vi in range(4):values[:,ci,vi]+=np.bincount(idx[roles==vi],minlength=n)
        for vi,view in enumerate(VIEWS):
            allowed=np.isin(arm['major'],list(PARENTS[view]))[nodes.codes]
            parent[:,vi]+=np.bincount(idx[allowed],minlength=n)
        ambiguous+=np.bincount(idx[arm['roles'][1]==-2],minlength=n)
        conflict+=np.bincount(idx[(arm['roles'][1]>=0)&(arm['roles'][2]<0)],minlength=n)
        nonempty=np.array([str(s).strip().casefold() not in {'','nan','none','unknown','not_reviewed','not reviewed','unavailable','unassigned'} for s in authors.categories])
        present=np.zeros(len(idx),bool);valid=authors.codes>=0;present[valid]=nonempty[authors.codes[valid]]
        author_present+=np.bincount(idx[present],minlength=n)
        # Vocab stays private and contains metadata only, never model-learned genes.
        encoded=arm['source'].codes.astype(np.int64)*(len(authors.categories)+1)+(authors.codes+1)
        unique,first,counts=np.unique(encoded,return_index=True,return_counts=True)
        for i,count in zip(first,counts):
            ac=int(authors.codes[i]);role=int(arm['author_mapping'][ac]) if ac>=0 else -1
            vocab.append(dict(arm=arm['arm'],dataset=str(arm['source'][i]),author_label=str(authors[i]) if ac>=0 else '',parsed_role=VIEWS[role] if role>=0 else 'AMBIGUOUS' if role==-2 else 'UNRESOLVED',cells=int(count),author_field=arm['author_field']))
        arm_proofs.append(dict(arm=arm['arm'],cells=arm['cell_count'],author_field=arm['author_field'],author_values_sha256=arm['author_values_sha256']))
    assert (values.sum(axis=2)<=lib.annotated_scope_cells.to_numpy()[:,None]).all()
    lib['author_label_present_cells']=author_present;lib['author_role_ambiguous_cells']=ambiguous;lib['author_role_parent_conflict_cells']=conflict
    lib.to_csv(out/'private_all1690_library_metadata_and_label_diagnostics.tsv',sep='\t',index=False)
    pd.DataFrame(vocab).to_csv(out/'private_all_source_author_label_role_vocabulary.tsv.gz',sep='\t',index=False,compression='gzip')
    definitions=[]
    for arm,catalog in data['catalogs'].items():
        for row in catalog.itertuples():definitions.append(dict(arm=arm,node_id=row.node_id,major=row.major,label=row.label,frozen_role=VIEWS[row.frozen_role] if row.frozen_role>=0 else 'UNRESOLVED'))
    pd.DataFrame(definitions).to_csv(out/'private_all442_frozen_node_role_definitions.tsv',sep='\t',index=False)
    np.savez_compressed(out/'private_all1690_three_channel_role_counts.npz',counts=values,parent_counts=parent,channels=np.array(CHANNELS),views=np.array(VIEWS))
    counts=[]
    for i,row in enumerate(lib.itertuples()):
        for ci,ch in enumerate(CHANNELS):
            for vi,view in enumerate(VIEWS):counts.append(dict(arm=row.arm,dataset=row.dataset,library_id=row.library_id,channel=ch,view=view,n_cells=int(values[i,ci,vi]),parent_cells=int(parent[i,vi]),captured=bool(values[i,ci,vi]>=20),parent_captured=bool(parent[i,vi]>=20),parent_present_role_unresolved=bool(parent[i,vi]>=20 and values[i,ci,vi]<20)))
    pd.DataFrame(counts).to_csv(out/'private_all20280_library_role_capture.tsv.gz',sep='\t',index=False,compression='gzip')
    source_rows=[];stratum_rows=[]
    for (arm,source),g in lib.groupby(['arm','dataset'],sort=True):
        ids=g.index.to_numpy()
        for ci,ch in enumerate(CHANNELS):
            for name,ends in scopes():
                ok=(values[:,ci,ends]>=20).all(axis=1)
                record=dict(arm=arm,dataset=source,channel=ch,capture_scope_id=name,all_source_libraries=len(g),**summary(lib,lib.index.isin(ids)&ok),exact_native_identity_accepted=False,independent_validation_accepted=False)
                source_rows.append(record)
                for levels,sg in g.groupby(STRATA,sort=True,dropna=False):
                    stratum_rows.append(dict(arm=arm,dataset=source,channel=ch,capture_scope_id=name,**dict(zip(STRATA,levels)),all_stratum_libraries=len(sg),**summary(lib,lib.index.isin(sg.index)&ok)))
    sources=pd.DataFrame(source_rows);assert len(sources)==2013
    sources.to_csv(out/'private_all2013_source_role_pair_capture.tsv',sep='\t',index=False)
    pd.DataFrame(stratum_rows).to_csv(out/'private_all_source_stratum_role_pair_capture.tsv.gz',sep='\t',index=False,compression='gzip')
    genes=pd.read_csv(root/'NSCLC_native_MFA_external_gene_qualification_20261016_v1/private_all2440_native_factor_source_gene_coverage.tsv',sep='\t')
    gene=genes.set_index(['arm','dataset','view','factor']).physical_coverage_pass
    opportunities=[]
    for row in sources.itertuples():
        ends=dict(scopes())[row.capture_scope_id]
        for k in range(10):
            factor='NSCLC_MFA_20261016_F%02d'%(k+1);coverage=all(bool(gene.loc[(row.arm,row.dataset,VIEWS[j],factor)]) for j in ends)
            opportunities.append(dict(arm=row.arm,dataset=row.dataset,channel=row.channel,capture_scope_id=row.capture_scope_id,factor=factor,all_endpoints_gene_gate=coverage,
                captured_libraries=row.captured_libraries,reviewed_patient_blocks=row.reviewed_patient_blocks,unreviewed_library_fallback=row.unreviewed_library_fallback,
                gene_and_capture_library_opportunity=bool(coverage and row.captured_libraries>0),exact_identity_or_external_validation_accepted=False))
    assert len(opportunities)==20130
    pd.DataFrame(opportunities).to_csv(out/'private_all20130_factor_source_capture_and_gene_opportunities.tsv.gz',sep='\t',index=False,compression='gzip')
    terminal=dict(status='COMPLETE_ALL61_SOURCE1690_LIBRARIES_THREE_ROLE_CHANNELS_ALL40_AND60_OPPORTUNITIES_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),
        source_sha256=sha(__file__),arm_input_proofs=arm_proofs,files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
        new_expression_or_P=False,identity_or_patient_independence_promoted=False)
    assert sum(v['bytes'] for v in terminal['files'].values())<268435456
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

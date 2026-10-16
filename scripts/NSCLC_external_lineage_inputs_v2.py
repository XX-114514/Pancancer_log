"""Read original author text with frozen cell rows; no annotation writes."""
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem

from NSCLC_native_MFA_clinical_inputs_v1 import matched_quota,sha
from NSCLC_external_lineage_role_rules_v1 import VIEWS,CHANNELS,PARENTS,author_role,phase,material,rule_fixture

STRATA=['cancer_label','tissue','sample_class','timepoint','treatment_arm','capture_scope','ex_vivo','official_stage','material_class','phase_class']


def load():
    rule_fixture();root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);capture=root/'complete_external_lineage_opportunities_20261016_v2'
    proof=capture/'private_independent_opportunity_audit.json';j=json.loads(proof.read_text())
    assert j['status'].startswith('PASS_ALL61_SOURCE_ARM_ALL1690_')
    terminal=capture/'private_terminal.json';assert j['producer_terminal_sha256']==sha(terminal)
    t=json.loads(terminal.read_text())
    names=['private_all_library_latest_metadata.tsv','private_all_original_node_library_counts.tsv.gz']
    for name in names:assert sha(capture/name)==t['files'][name]
    lib=pd.read_csv(capture/names[0],sep='\t',keep_default_na=False,dtype={'dataset':str,'library_id':str,'patient_block':str})
    lib=lib.sort_values(['arm','dataset','library_id']).reset_index(drop=True)
    counts=pd.read_csv(capture/names[1],sep='\t',keep_default_na=False,dtype={'dataset':str,'library_id':str,'node_id':str})
    assert len(lib)==1690 and pd.api.types.is_bool_dtype(lib.reviewed_patient)
    assert not lib.duplicated(['arm','dataset','library_id']).any()
    lib['phase_class']=lib.timepoint.map(phase);lib['material_class']=[material(r) for r in lib.itertuples()]
    old=root/'all_CM_state_native_input_preflight_20261014_v1';inventory=json.loads((old/'inventory.json').read_text())
    annotations=root.parents[2]/'Project_v3/annotation_rebuild'
    dictionary=annotations/'20261010_expression_annotation_downstream_v1/downstream_review_v2/number_to_annotation_dictionary.tsv'
    contract_path=annotations/'20261010_uniform_expression_annotation_v1/input_contract.json'
    assert sha(dictionary)==inventory['input_sha256']['current_catalogue'] and sha(contract_path)==inventory['input_sha256']['input_contract']
    assert sha(old/'private_canonical_row_selection.npz')==inventory['private_row_selection_sha256']
    oldcat=pd.read_csv(dictionary,sep='\t',keep_default_na=False)[['node_id','subCluster','majorCluster']].rename(columns={'subCluster':'label','majorCluster':'major'})
    oldcat['frozen_role']=-1
    oldcat.loc[oldcat.label.str.startswith('CD8_T_cell_'),'frozen_role']=0
    oldcat.loc[oldcat.label.str.startswith('Regulatory_T_cell_'),'frozen_role']=1
    oldcat.loc[oldcat.major.eq('Macrophage'),'frozen_role']=2
    oldcat.loc[oldcat.major.eq('B cell'),'frozen_role']=3
    atlas=root.parents[2]/'Project_v5/runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'
    z=zarr.open_group(str(atlas),mode='r');text=read_elem(z['uns']['subcluster_catalog_json']);text=text.item() if isinstance(text,np.ndarray) else text
    newcat=pd.DataFrame(json.loads(text));candidate=newcat.subtype_evidence_status.eq('multigene_subtype_candidate')
    newcat['frozen_role']=-1
    newcat.loc[newcat.major.eq('T')&newcat.subtype_candidate.eq('CD8_Effector')&candidate,'frozen_role']=0
    newcat.loc[newcat.major.eq('T')&newcat.subtype_candidate.eq('Treg')&candidate,'frozen_role']=1
    newcat.loc[newcat.major.eq('Myeloid')&newcat.subtype_candidate.str.startswith('Macro_')&candidate,'frozen_role']=2
    newcat.loc[newcat.major.isin(['B','Plasma']),'frozen_role']=3
    newcat=newcat[['node_id','subcluster','major','frozen_role']].rename(columns={'subcluster':'label'})
    assert len(oldcat)==172 and len(newcat)==270
    gene=root/'NSCLC_native_MFA_external_gene_qualification_20261016_v1'
    assert sha(gene/'private_independent_gene_qualification_audit.json')=='096d5166f1016fac39ca4feea462a3399eb7b120c3adfa1e8166d20f74c09a22'
    question=gene/'PRIVATE_NEXT_COMPLETE_LINEAGE_CASE_AND_CALIBRATION_SCOPE.json'
    assert sha(question)=='15423fd5ccaff0eefdc80ee8c39a876a9c6db59de2b4ae16b9d124ce0c1bb879'
    catalogs={'frozen_original_annotation':oldcat,'archived_immune_atlas':newcat}
    return dict(root=root,lib=lib,counts=counts,catalogs=catalogs,old=old,old_zarr=json.loads(contract_path.read_text())['core_zarr'],atlas=atlas,
        hashes=dict(capture_proof=sha(proof),capture_terminal=sha(terminal),libraries=sha(capture/names[0]),counts=sha(capture/names[1]),
            old_dictionary=sha(dictionary),old_row_selection=sha(old/'private_canonical_row_selection.npz'),atlas_catalog=hashlib.sha256(text.encode()).hexdigest(),next_question=sha(question)))


def cell_arms(data):
    """Yield encoded cells; original metadata stay exclusively on the server."""
    lookup={(r.arm,r.dataset,r.library_id):i for i,r in enumerate(data['lib'].itertuples())}
    for arm in ['frozen_original_annotation','archived_immune_atlas']:
        original=arm=='frozen_original_annotation'
        z=zarr.open_group(str(data['old_zarr'] if original else data['atlas']),mode='r')
        rows=None
        if original:
            with np.load(data['old']/'private_canonical_row_selection.npz',allow_pickle=False) as f:rows=f['global_rows']
        fields=['gse_id','sample_uid','expression_v1_node_id','ev_author_label'] if original else ['dataset','library_id','node_id','source_author_label']
        arrays=[]
        for field in fields:
            values=read_elem(z['obs'][field]);values=values[rows] if rows is not None else values
            arrays.append(pd.Categorical(values).remove_unused_categories())
        source,library,node,author=arrays;n=len(node)
        assert all(len(x)==n for x in arrays) and (source.codes>=0).all() and (library.codes>=0).all() and (node.codes>=0).all()
        packed=source.codes.astype(np.int64)*len(library.categories)+library.codes
        _,first,inverse=np.unique(packed,return_index=True,return_inverse=True)
        unique_lib=np.array([lookup[(arm,str(source[i]),str(library[i]))] for i in first])
        lib_index=unique_lib[inverse]
        catalog=data['catalogs'][arm].set_index('node_id');names=np.array(node.categories).astype(str)
        assert set(names)<=set(catalog.index)
        major=catalog.loc[names,'major'].to_numpy(str);frozen=catalog.loc[names,'frozen_role'].to_numpy(int)[node.codes]
        assert len(author.categories)<=20000,'Author field resembles observation identifiers, do not infer roles'
        mapping=np.array([author_role(s) for s in author.categories],dtype=np.int8)
        parsed=np.full(n,-1,np.int8);valid=author.codes>=0;parsed[valid]=mapping[author.codes[valid]]
        concordant=parsed.copy()
        for j,view in enumerate(VIEWS):
            allowed=np.isin(major,list(PARENTS[view]))[node.codes]
            concordant[(parsed==j)&~allowed]=-1
        fresh=np.bincount(lib_index*len(names)+node.codes,minlength=len(data['lib'])*len(names)).reshape(len(data['lib']),len(names))
        inherited=data['counts'][data['counts'].arm.eq(arm)];ni={s:j for j,s in enumerate(names)}
        assert np.count_nonzero(fresh)==len(inherited)
        for row in inherited.itertuples():assert fresh[lookup[(arm,row.dataset,row.library_id)],ni[row.node_id]]==row.n_cells
        assert fresh.sum()==n==int(data['lib'].loc[data['lib'].arm.eq(arm),'annotated_scope_cells'].sum())
        digest=hashlib.sha256(json.dumps([str(s) for s in author.categories],ensure_ascii=False,separators=(',',':')).encode()+author.codes.astype('<i4').tobytes()).hexdigest()
        yield dict(arm=arm,source=source,library_index=lib_index,node=node,major=major,author=author,author_mapping=mapping,
            roles=[frozen,parsed,concordant],author_field=fields[-1],author_values_sha256=digest,cell_count=n)

"""Freeze exploratory annotation memberships before any human-data CM fit.

All arms share the accepted compartment denominator and numerical graph.
Reference conflicts remain explicit denominator-retained fallback nodes.
SingleR proposals are group-level, not new per-cell reference predictions.
"""
import argparse, datetime, hashlib, json, sys
from pathlib import Path
import numpy as np
import pandas as pd

def load(p): return json.loads(Path(p).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    with Path(p).open('x') as f: json.dump(v,f,indent=2,allow_nan=False)
def binding(p): return dict(path=str(p),SHA256=sha(p))
def bound(b):
    assert sha(b['path'])==b['SHA256']; return Path(b['path'])

HPCA_MAIN = {
    'B_cell':'B_Plasma','T_cells':'T_NK','NK_cell':'T_NK',
    'DC':'Myeloid','Macrophage':'Myeloid','Monocyte':'Myeloid',
    'Myelocyte':'Myeloid','Neutrophils':'Myeloid','Pro-Myelocyte':'Myeloid',
    'Epithelial_cells':'Epithelial','Keratinocytes':'Epithelial',
    'Hepatocytes':'Epithelial','Endothelial_cells':'Endothelial',
    'Fibroblasts':'Mesenchymal_stromal','MSC':'Mesenchymal_stromal',
    'Smooth_muscle_cells':'Mesenchymal_stromal',
    'Chondrocytes':'Mesenchymal_stromal','Osteoblasts':'Mesenchymal_stromal',
    'Astrocyte':'Neural_related','Neuroepithelial_cell':'Neural_related',
    'Neurons':'Neural_related','Erythroblast':'Erythroid','Platelets':'Platelet',
    'BM':None,'BM & Prog.':None,'CMP':None,'Embryonic_stem_cells':None,
    'GMP':None,'Gametocytes':None,'HSC_-G-CSF':None,'HSC_CD34+':None,
    'MEP':None,'Pre-B_cell_CD34-':None,'Pro-B_cell_CD34+':None,
    'Tissue_stem_cells':None,'iPS_cells':None}
MONACO_MAIN = {'B cells':'B_Plasma','Basophils':'Myeloid',
    'CD4+ T cells':'T_NK','CD8+ T cells':'T_NK','T cells':'T_NK',
    'NK cells':'T_NK','Dendritic cells':'Myeloid','Monocytes':'Myeloid',
    'Neutrophils':'Myeloid','Progenitors':None}

def main(R):
    cfg=load(R/'config/spec_v1.json')
    sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(R,cfg)
    scope_root=Path(cfg['scope']); sem_root=Path(cfg['semantics'])
    sr_root=Path(cfg['SingleR']); audit_root=Path(cfg['SingleR_auditor'])
    B=Path(cfg['naming']); H=Path(cfg['model']); closure=Path(cfg['unit_closure'])
    for root in [scope_root,sem_root,sr_root,audit_root,B,closure]:
        assert load(root/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    scope_p=scope_root/'outputs/accepted_shared_compartment_scope_v1.json'
    assert sha(scope_p)==load(scope_root/'provenance/terminal_v1.json')['acceptance_SHA256']
    scope=load(scope_p); assert scope['accepted_purpose']==cfg['accepted_purpose']
    sem_p=bound(scope['reference_semantics_acceptance']); assert sem_p.parent.parent==sem_root
    sem=load(sem_p); assert sem['core_coverage_completion_acceptance'], 'Primary core coverage must be complete'
    assert scope['input_axis_SHA256']==sha(H/'outputs/private_cell_axes_v1.npz')
    p=np.load(bound(scope['codes'])); N=len(p); parents=scope['compartments']
    groups=np.load(B/'outputs/private_full_fit_cluster_code_v1.npy')
    catalog_p=B/'outputs/marker_first_cluster_catalog_v1.json'; cat=load(catalog_p)
    assert sha(catalog_p)==scope['marker_catalog_SHA256'] and len(groups)==N
    cluster_rows={r['cluster']:r for r in cat['clusters']}
    axis=np.load(H/'outputs/private_cell_axes_v1.npz'); source=axis['source_code']
    sources=load(H/'outputs/scope_and_common_genes_v1.json')['sources']
    unit_record=load(closure/'outputs/accepted_fit_qualification_unit_bridge_v3.json')
    unit=np.load(bound(unit_record['unit_codes'])); assert unit.shape==(N,)
    nu=int(unit.max())+1
    ontology=load(bound(sem['exact_ontology']))
    policy=load(bound(cfg['compartment_policy']))
    family_of={a:k for k,v in policy['compartments'].items() for a in v['identity_parents']}
    immune={parents.index(x) for x in ['T_NK','B_Plasma','Myeloid','Mast']}
    semantics_models={x['model']:x for x in sem['models']}
    proposals=[]
    # Fine labels are reference phenotypes in one view, identity coarsenings
    # in another. Neither raw state words nor tissue words become state assays.
    low=semantics_models['immune_low']; raw=np.load(bound(low['raw_prediction_codes']))
    assert (raw>=0).all(), 'Completed core must cover all fitted cells'
    rows=ontology['models']['immune_low']['rows']
    assert len(rows)>int(raw.max())
    low_family=np.array([parents.index(family_of[r['coarse_reference_parent']])
        if r['coarse_reference_parent'] in family_of else -1 for r in rows],np.int16)
    for scheme,field in [('CT_reference_phenotype_cluster','label'),
                         ('CT_identity_cluster','canonical_identity')]:
        vocabulary=sorted({r[field] for r in rows}); by={x:i for i,x in enumerate(vocabulary)}
        label_map=np.array([by[r[field]] for r in rows],np.int32)
        fine=np.full(N,-1,np.int32)
        in_immune=np.isin(p,list(immune)); compatible=in_immune&(p==low_family[raw])
        fine[compatible]=label_map[raw[compatible]]
        proposals.append((scheme,fine,vocabulary,'per_cell_official_CellTypist_best_match'))
    proposals.insert(0,('shared_compartment_cluster',np.full(N,-1,np.int32),[],
                         'common_marker_reference_native_scope_numeric_cluster_control'))
    # Use the actual training-reference main/fine relationship. In particular,
    # no string fragment in a fine name is used to infer its biological parent.
    proposal_p=sr_root/'outputs/private_all_matched_SingleR_proposals_v1.tsv'
    df=pd.read_csv(proposal_p,sep='\t',keep_default_na=False)
    assert not df.duplicated(['source','cluster','reference','label_level']).any()
    for reference,mainmap in [('HPCA',HPCA_MAIN),('Monaco_immune',MONACO_MAIN)]:
        vocab_p=sr_root/'outputs'/('private_reference_'+reference+'_vocabulary_v1.tsv')
        voc=pd.read_csv(vocab_p,sep='\t',keep_default_na=False)
        assert set(voc.label_main)<=set(mainmap),'Unreviewed reference main label'
        fineparents={}
        for label,x in voc.groupby('label_fine'):
            fam={mainmap[v] for v in x.label_main}
            fineparents[label]=next(iter(fam)) if len(fam)==1 else None
        vocabulary=sorted(fineparents); label_index={x:i for i,x in enumerate(vocabulary)}
        lookup=np.full((len(sources),int(groups.max())+1),-1,np.int32)
        parent_lookup=np.full_like(lookup,-1)
        sub=df[(df.reference==reference)&(df.label_level=='fine')]
        assert len(sub)>0
        for row in sub.itertuples():
            label=row.pruned_label
            if label in ['', 'NA', 'NaN']: continue
            assert label in fineparents
            family=fineparents[label]
            if family is None: continue
            sc=sources.index(row.source); gr=int(row.cluster)
            lookup[sc,gr]=label_index[label];parent_lookup[sc,gr]=parents.index(family)
        candidate=lookup[source,groups]; proposed_parent=parent_lookup[source,groups]
        fine=np.where((p>=0)&np.isin(p,list(immune))&(p==proposed_parent),candidate,-1).astype(np.int32)
        proposals.append(('SingleR_'+reference+'_group_cluster',fine,vocabulary,
                          'source_by_cluster_group_profile_official_SingleR_pruned_fine'))
    output=[]
    for scheme,fine,vocabulary,origin in proposals:
        # Nonimmune classes use the identical numeric-cluster partition in all
        # primary arms. Organ-specific references have separate assay domains.
        numeric=(p>=0)&(~np.isin(p,list(immune)))
        if scheme=='shared_compartment_cluster': numeric=p>=0
        fine[numeric]=-2
        good=p>=0
        keys=np.stack([p[good],groups[good],fine[good]],axis=1)
        unique,inverse=np.unique(keys,axis=0,return_inverse=True)
        node=np.full(N,-1,np.int32);node[good]=inverse
        nodes=[]
        for k,(parent,cluster,fid) in enumerate(unique):
            parent=int(parent);cluster=int(cluster);fid=int(fid)
            description='numeric_cluster' if fid==-2 else ('reference_unresolved_or_conflict' if fid==-1 else vocabulary[fid])
            # Display is deliberately decoupled from membership. Own cluster
            # markers are not borrowed across a rescued compartment boundary.
            row=cluster_rows[cluster]; marker_family=family_of.get(row['major_marker_first'])
            marker=row['marker'] if marker_family==parents[parent] else None
            nodes.append(dict(code=k,name='N_'+str(k),parent_code=parent,
                numeric_cluster=cluster,reference_label=description,
                discovery_eligible=fid!=-1,own_node_DEG_marker_not_yet_assigned=True,
                parent_cluster_marker=marker,parent_marker_borrowed=marker is not None,
                reference_or_numeric_definition=origin,
                reference_fallback_retained_in_full_parent_denominator=fid==-1))
        # All primary core models actually ran in every source. SingleR ran in
        # every observed source/group; an unpruned call is an observed fallback,
        # not a missing assay. The opportunities are explicit, never inferred
        # from which labels happened to be present in a work unit.
        opportunity=np.ones((len(nodes),nu),bool)
        npth=R/'outputs'/('private_'+scheme+'_node_codes_v1.npy')
        ppth=R/'outputs'/('private_'+scheme+'_parent_codes_v1.npy')
        opth=R/'outputs'/('private_'+scheme+'_node_measurement_v1.npy')
        np.save(npth,node);np.save(ppth,p);np.save(opth,opportunity)
        np.testing.assert_array_equal(np.load(npth),node)
        mapping=np.asarray([n['parent_code'] for n in nodes]); np.testing.assert_array_equal(mapping[node[good]],p[good])
        # Independent scalar label/cluster replay samples every realized node.
        first=np.full(len(nodes),N,np.int64);np.minimum.at(first,inverse,np.flatnonzero(good))
        for k,i in enumerate(first):
            assert tuple(unique[k])==(int(p[i]),int(groups[i]),int(fine[i]))
        a=dict(schema='accepted_TME_identity_scheme_membership_v1',
            accepted_purpose=cfg['accepted_purpose'],scheme_id=scheme,n_fit=N,
            final_membership_codes_not_display_names=True,
            unresolved_fine_identity_has_explicit_fallback_node=True,
            parents=[dict(code=i,name=n) for i,n in enumerate(parents)],nodes=nodes,
            parent_codes=binding(ppth),node_codes=binding(npth),node_measurement_mask=binding(opth),
            input_axis_SHA256=scope['input_axis_SHA256'],scope_acceptance=binding(scope_p),
            reference_semantics_acceptance=binding(sem_p),SingleR_proposals=binding(proposal_p),
            numeric_cluster_catalog=binding(catalog_p),n_assessed_cells=int(good.sum()),
            n_nodes=len(nodes),n_scalar_membership_replays=len(nodes),
            nonimmune_numeric_cluster_memberships_shared=True,
            SingleR_group_assignments_not_per_cell_inference=scheme.startswith('SingleR'),
            state_hints_not_multigene_state_scores=True,
            no_clinical_outcome_or_CM_result_used=True)
        ap=R/'outputs'/('private_'+scheme+'_accepted_membership_v1.json');save(ap,a)
        output.append(dict(scheme_id=scheme,acceptance=binding(ap),n_nodes=len(nodes),
            n_discovery_eligible_nodes=sum(x['discovery_eligible'] for x in nodes),
            n_fallback_cells=int(np.sum(good&(fine==-1)))))
    # The default discovery view uses complete work units. It is exploratory;
    # source/patient-held-out validation must have a separately frozen mask.
    train=R/'outputs/private_all_unit_discovery_cell_mask_v1.npy';np.save(train,np.ones(N,bool))
    storage_gate(R,cfg)
    mp=R/'outputs/accepted_annotation_membership_collection_v1.json'
    save(mp,dict(schema='accepted_TME_annotation_membership_collection_v1',
        accepted_purpose=cfg['accepted_purpose'],n_fit=N,schemes=output,
        training_cell_mask=dict(**binding(train),frozen_before_module_fitting=True),
        no_heldout_validation_claim=True,input_axis_SHA256=scope['input_axis_SHA256']))
    save(R/'provenance/terminal_v1.json',dict(status='PRE_CM_ANNOTATION_MEMBERSHIP_VIEWS_FROZEN',
        completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        acceptance_SHA256=sha(mp),n_schemes=len(output)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try: main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise

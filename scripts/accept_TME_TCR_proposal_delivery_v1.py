"""Accept the full proposal audit, exercise the reader and export aggregates."""
import argparse,collections,csv,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
def table(p,rows):
    with Path(p).open('x',newline='',encoding='utf-8-sig') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def main(root):
    begin=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['producer']);A=Path(cfg['auditor']);J=Path(cfg['joined'])
    for s in [P,A]:assert load(s/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    pc=load(P/'config/spec_v1.json');pt=P/'provenance/terminal_v1.json';at=A/'provenance/terminal_v1.json';assert load(at)['producer_terminal_SHA256']==sha(pt)
    ip=P/'outputs/TCR_boundary_proposal_index_v1.json';assert sha(ip)==load(pt)['index_SHA256'];sp=P/'outputs/TCR_boundary_proposal_summary_v1.json';assert sha(sp)==load(pt)['summary_SHA256']
    sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(root/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(root,pc)
    from TME_TCR_boundary_policy_v1 import REASONS
    from TME_TCR_proposal_view_v1 import TMETCRProposalView
    bindings={k:{'path':str(p),'SHA256':sha(p)} for k,p in {'working_acceptance':Path(pc['working_acceptance']),'proposal_index':ip,'proposal_terminal':pt,'independent_audit_terminal':at}.items()}
    contract={'schema':'accepted_TME_TCR_boundary_proposal_v1','accepted_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'bindings':bindings,'evidence_version':'TME_additional_TCR_RNA_coarse_proposal_20261027_v1','review_dictionary':REASONS,
        'acceptance_scope':'all_scope_numeric_RNA_proposal_and_reader_not_biological_accuracy_or_current_consensus_release',
        'changed_membership_full_gene_marker_gate_pending':True,'current_consensus_and_states_preserved':True,'final_module_input_ready':False,'new_CM_started':False}
    ap=root/'outputs/accepted_TME_TCR_boundary_proposal_v1.json';save(ap,contract);view=TMETCRProposalView(ap)
    E=Path(pc['native_evidence']);ei={r['source']:r for r in load(E/'outputs/TME_full_native_evidence_index_v1.json')}
    lmroot=Path(load(J/'outputs/joined_native_atlas_delivery_v1.json')['local_review']);li={r['source']:r for r in load(lmroot/'outputs/full_local_identity_review_index_v1.json')}
    catalog=view.base.groups;bygroup={(g['source'],g['source_annotation_group_code']):g for g in catalog}
    cases={};coverage=[];splits=[];gaprows=[];gapcounts=collections.Counter();totals=collections.Counter()
    for source in view.sources:
        s=view.source(source);a=s['arrays'];m=s['metadata'];b=view.base.source(source);q=b['arrays'];groups=m['affected_previous_exact_groups']
        for role,mask in [('proposed_T',a['proposed_T']),('competing_core',a['review_code']==3),('unmeasured_competitor',a['review_code']==4),('local_group_discordant',a['review_code']==7)]:
            if role not in cases and mask.any():cases[role]=(source,np.flatnonzero(mask)[:2])
        major=np.asarray(b['metadata']['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code']]
        if 'baseline_NK_control' not in cases and (major=='NK').any():cases['baseline_NK_control']=(source,np.flatnonzero(major=='NK')[:2])
        coverage.append({'source':source,'n_scope_records':m['n_scope_records'],'n_baseline_unresolved':m['n_baseline_unresolved'],'n_additional_T_coarse_proposed':m['n_proposed_T'],
            'n_TCR_core':m['n_TCR_core_records'],'n_local_group_discordant':m['reason_counts_in_baseline_unresolved'].get(REASONS[7],0),
            'n_affected_previous_groups':len(groups),'proposals_applied_to_current_consensus':False})
        for r in groups:
            old=bygroup[(source,r['source_annotation_group_code'])]
            splits.append({'source':source,'previous_group_number':old['global_subcluster_number'],'previous_group_SHA256':old['stable_group_sha256'],
                'n_previous_group':r['n_previous_exact_group'],'n_additional_T_proposed':r['n_proposed_T'],'n_unresolved_remainder':r['n_remainder'],
                'new_exact_membership_and_marker_release_pending':True})
        er=ei[source];ep=E/'outputs'/er['metadata'];assert sha(ep)==er['metadata_SHA256'];em=load(ep);markers={g['source_annotation_group_code']:[r['gene'] for r in g['own_native_expression_marker_candidates'][:8]] for g in em['groups']}
        lr=li[source];lp=lmroot/'outputs'/lr['output'];assert sha(lp)==lr['output_sha256'];local={r['source_local_group_code']:r['local_major_proposal'] for r in load(lp)['groups']}
        held=a['review_code']==7
        for code in np.unique(q['annotation_group_code'][held]):
            sel=q['annotation_group_code']==code;old=bygroup[(source,int(code))];nheld=int((sel&held).sum());parent=old['parent_source_local_group_code'];label=local.get(parent,'not_available')
            gapcounts[label]+=nheld
            gaprows.append({'source':source,'previous_group_number':old['global_subcluster_number'],'parent_local_major':label,'n_scope_group':int(sel.sum()),
                'n_TCR_core':int((sel&a['TCR_core']).sum()),'n_TCR_and_competing_core':int((sel&a['TCR_core']&(a['competing_core_mask']!=0)).sum()),
                'n_local_context_discordant':nheld,'current_group_marker_candidates':'|'.join(markers[int(code)]),
                'review_purpose':'source_local_group_reassessment_not_automatic_T_assignment'})
        totals.update(n_sources=1,n_scope_records=m['n_scope_records'],n_baseline_unresolved=m['n_baseline_unresolved'],n_proposed_T=m['n_proposed_T'],n_affected_previous_groups=len(groups),n_split_previous_groups=sum(r['n_remainder']>0 for r in groups))
    assert dict(totals)==load(pt)['totals'];assert set(cases)=={'proposed_T','competing_core','unmeasured_competitor','local_group_discordant','baseline_NK_control'}
    checks=[]
    for role,(source,pos) in cases.items():
        result=view.metadata(source,pos);base=view.base.metadata(source,pos)
        for key in ['native_row','source_native_position','consensus_major','consensus_celltype','label_version','subcluster','global_subcluster_number','stable_group_sha256']:
            np.testing.assert_array_equal(result['columns'][key],base['columns'][key])
        if role=='proposed_T':assert all(v=='T' for v in result['columns']['additional_TCR_coarse_proposal']) and all(v=='unresolved' for v in result['columns']['consensus_major'])
        else:assert all(v is None for v in result['columns']['additional_TCR_coarse_proposal'])
        actual=view.read_counts(source,pos,['CD3D','TRAC','__UNMEASURED_SENTINEL__']);expected=view.base.read_counts(source,pos,['CD3D','TRAC','__UNMEASURED_SENTINEL__'])
        assert (actual['counts']!=expected['counts']).nnz==0
        for k in ['native_row','native_total_counts','measured']:np.testing.assert_array_equal(actual[k],expected[k])
        assert not actual['measured'][-1];checks.append({'role':role,'n_rows':len(pos),'baseline_labels_native_counts_and_missingness_preserved':True})
    try:view.metadata(view.sources[0],[-1])
    except IndexError:pass
    else:raise AssertionError('Negative position accepted')
    try:view.source('__out_of_scope__')
    except KeyError:pass
    else:raise AssertionError('Out-of-scope source accepted')
    table(root/'outputs/TCR_proposal_source_coverage_v1.tsv',coverage);table(root/'outputs/TCR_proposed_group_split_preparation_v1.tsv',splits);table(root/'outputs/TCR_local_group_context_review_priority_v1.tsv',gaprows)
    summary=load(sp);summary.update(local_group_discordance_by_parent_local_major=dict(gapcounts),n_local_context_discordant_exact_groups=len(gaprows),bounded_reader_cases=checks,all_existing_working_labels_and_native_counts_preserved=True,
        accepted_proposal_contract_SHA256=sha(ap),full_group_tables_saved_server_only=True)
    save(root/'outputs/TCR_proposal_delivery_summary_v1.json',summary)
    save(root/'provenance/terminal_v1.json',{'status':'FULL_TCR_RNA_PROPOSAL_LAYER_READER_AND_SOURCE_GROUP_SUMMARIES_DELIVERED_NOT_APPLIED',
        'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':dict(totals),'contract_SHA256':sha(ap),
        'files':{p.name:sha(p) for p in (root/'outputs').glob('*') if p.is_file()},'changed_membership_marker_release_pending':True,'final_module_input_ready':False,'wall_seconds':time.monotonic()-begin})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
    try:main(r)
    except BaseException:
        import traceback
        save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise

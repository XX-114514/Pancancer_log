"""Replay all preflight counts/LR means from original native cell caches.

No producer import. Reconstruct diagnostic design from per-cell RNA and gene
counts, and rank via pivoted QR; fullgene targets use previously audited arrays.
"""
from pathlib import Path
from datetime import datetime, timezone
import os,json,hashlib,gc
import numpy as np
import pandas as pd
from scipy.linalg import qr


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for v in iter(lambda:f.read(8388608),b''):h.update(v)
    return h.hexdigest()


def truth(v):return v is True or v==True or v=='True'


def main():
    start=datetime.now(timezone.utc).isoformat();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';O=P/'fullgene_receiver_response_preflight_v1';F=P/'fullgene_receiver_target_readiness_v1';N=P/'supported_LR_native_cell_panel_v1'
    assert not (O/'independent_audit.json').exists();inv=json.loads((O/'inventory.json').read_text());plan=json.loads((R/'fullgene_receiver_target_code_v1/fullgene_receiver_response_preflight_plan_v1.json').read_text())
    assert sha(__file__)==plan['code_sha256'][Path(__file__).name]
    inputs=dict(full_inventory=F/'inventory.json',full_audit=F/'independent_audit.json',native_inventory=N/'inventory.json',provenance_inventory=P/'patient_block_provenance_overlay_v1/inventory.json',provenance_audit=P/'patient_block_provenance_overlay_v1/independent_audit.json',membership=P/'membership_review_v2/source_node_analysis_membership.tsv',fixed_candidates=N/'fixed336_reference_supported_candidates.tsv')
    for name,path in inputs.items():assert sha(path)==plan['input_sha256'][name]

    for n,v in inv['outputs'].items():assert sha(O/n)==v
    assert inv['plan_sha256']==sha(R/'fullgene_receiver_target_code_v1/fullgene_receiver_response_preflight_plan_v1.json')
    scope=pd.read_csv(P/'patient_block_provenance_overlay_v1/private_all1290_patient_provenance_overlay.tsv',sep='\t',keep_default_na=False).set_index('sampleID')
    provenance=json.loads((P/'patient_block_provenance_overlay_v1/inventory.json').read_text());assert sha(P/'patient_block_provenance_overlay_v1/private_all1290_patient_provenance_overlay.tsv')==provenance['outputs']['private_all1290_patient_provenance_overlay.tsv']
    fixed=pd.read_csv(N/'fixed336_reference_supported_candidates.tsv',sep='\t',keep_default_na=False);fixed['model_key']=fixed[['context_index','program','sender_family','receiver_family','ligand','receptor']].astype(str).agg('|'.join,axis=1)
    cols=['gse_id','stratum','sender_family','receiver_family','ligand','receptor'];mapped=pd.read_csv(O/'all336_to_unique_interfaces.tsv',sep='\t',keep_default_na=False);d=pd.read_csv(O/'all_interface_original_node_preflight.tsv',sep='\t',keep_default_na=False)
    assert len(mapped)==336 and mapped.model_key.is_unique
    pd.testing.assert_frame_equal(mapped[['model_key']+cols].sort_values('model_key').reset_index(drop=True),fixed[['model_key']+cols].sort_values('model_key').reset_index(drop=True),check_dtype=False)
    interfaces=mapped[['interface_id']+cols].drop_duplicates();assert interfaces.interface_id.is_unique
    metadata=interfaces.set_index('interface_id')
    for r in d.itertuples():
        for name in cols:assert getattr(r,name)==metadata.loc[r.interface_id,name]
    member=pd.read_csv(P/'membership_review_v2/source_node_analysis_membership.tsv',sep='\t');catalogue=member.loc[member.analysis_family_main_membership,['subCluster','family']].drop_duplicates()
    expected={(r.interface_id,n) for r in interfaces.itertuples() for n in catalogue[catalogue.family.eq(r.receiver_family)].subCluster}
    assert set(zip(d.interface_id,d.receiver_original_node))==expected and not d[['interface_id','receiver_original_node']].duplicated().any()
    priors=pd.read_csv(F/'all336_prior_coverage.tsv',sep='\t',keep_default_na=False);checked=0;designs=0;source_counts={}
    for rec in json.loads((N/'inventory.json').read_text())['source_records']:
        sid=rec['gse_id'];assert sha(rec['selected_cache'])==rec['cache_sha256'];z=dict(np.load(rec['selected_cache'],allow_pickle=False));b=dict(np.load(F/'private_source_unit_arrays'/(sid+'.npz'),allow_pickle=False));samples=z['samples'].astype(str);ns=len(samples)
        family=z['families'][z['family_code']];rawgi={g:i for i,g in enumerate(z['genes'])};gi={g:i for i,g in enumerate(b['genes'])};nodes=list(b['nodes']);Xgenes=b['mean_log1p_CP10K'].reshape(len(nodes),ns,len(gi));det=b['positive_cell_counts'].reshape(len(nodes),ns,len(gi))
        s=scope.loc[samples];counts={};moments={};joint={}
        for row in d[d.gse_id.eq(sid)].itertuples():
            node=row.receiver_original_node;captured=node in nodes;assert truth(row.receiver_node_captured)==captured
            if not captured:assert row.status=='RECEIVER_ORIGINAL_NODE_NOT_CAPTURED';checked+=1;continue
            family_mask=family==row.sender_family;node_mask=z['subCluster']==node
            sender_n=np.bincount(z['sample_code'][family_mask],minlength=ns);receiver_n=np.bincount(z['sample_code'][node_mask],minlength=ns)
            context=np.array([r.stratum==row.stratum and bool(r.analysis_eligible) for r in s.itertuples()]);ok=context & (sender_n>=20) & (receiver_n>=20);picked=s.iloc[np.flatnonzero(ok)]
            nb=picked.effective_resampling_block.nunique();npb=picked.loc[picked.unit_kind.eq('reviewed_patient'),'effective_resampling_block'].nunique()
            assert int(row.context_libraries)==int(context.sum()) and int(row.n_libraries_endpoint_cells_ge20)==int(ok.sum()) and int(row.n_analysis_blocks)==nb and int(row.n_reviewed_patient_blocks)==npb
            assert truth(row.library_gate20)==(ok.sum()>=20) and truth(row.block_support_ge10)==(nb>=10) and truth(row.block_support_ge20)==(nb>=20)
            L=row.ligand.split('_');Q=row.receptor.split('_');measured=all(g in gi for g in L+Q);assert truth(row.all_LR_constituents_in_native_axis)==measured
            raw_meas=all(g in rawgi for g in L+Q);assert truth(row.joint_RNA_panel_measurable)==raw_meas
            if raw_meas:
                lc=np.bincount(z['sample_code'][family_mask & np.all(z['raw_counts'][:,[rawgi[g] for g in L]]>0,axis=1)],minlength=ns)
                rc=np.bincount(z['sample_code'][node_mask & np.all(z['raw_counts'][:,[rawgi[g] for g in Q]]>0,axis=1)],minlength=ns)
                assert int(float(row.n_libraries_joint_LR_RNA_positive_cells_ge20))==int((ok & (lc>=20) & (rc>=20)).sum())
            pr=priors[(priors.gse_id.eq(sid)) & priors.sender_family.eq(row.sender_family) & priors.receiver_family.eq(row.receiver_family) & priors.ligand.eq(row.ligand) & priors.receptor.eq(row.receptor)].iloc[0]
            names=pr.top200_after_exclusions.split(';') if pr.top200_after_exclusions else [];indices=[gi[g] for g in names if g in gi];assert int(row.top_targets_fixed)==len(names) and int(row.top_targets_in_native_axis)==len(indices)
            k=nodes.index(node);var=0
            if indices and ok.any():
                expression=Xgenes[k,ok][:,indices];positive=det[k,ok][:,indices].sum(axis=0)/receiver_n[ok].sum();var=int((np.var(expression,axis=0)>1e-12).sum()) if len(expression)>1 else 0
                assert int(float(row.targets_detected1pct))==int((positive>=.01).sum()) and int(float(row.targets_detected10pct))==int((positive>=.1).sum()) and int(float(row.variable_target_genes))==var
            lvary=rvary=fullrank=False
            if measured and ok.sum()>=2:
                assert raw_meas,'LR is native-measured but missing accepted cell cache: requires separate raw replay'
                features=[]
                for i in np.flatnonzero(ok):
                    sp=np.flatnonzero(family_mask & (z['sample_code']==i));rp=np.flatnonzero(node_mask & (z['sample_code']==i))
                    ls=np.mean(z['raw_counts'][sp][:,[rawgi[g] for g in L]]*(10000/z['full_RNA_total'][sp])[:,None],axis=0)
                    rs=np.mean(z['raw_counts'][rp][:,[rawgi[g] for g in Q]]*(10000/z['full_RNA_total'][rp])[:,None],axis=0)
                    features.append([np.log1p(ls).mean(),np.log1p(rs).mean(),np.log1p(len(sp)),np.log1p(len(rp)),np.log1p(z['full_RNA_total'][sp].mean()),np.log1p(z['full_RNA_total'][rp].mean()),np.log1p(z['native_detected_genes'][sp].mean()),np.log1p(z['native_detected_genes'][rp].mean())])
                a=np.asarray(features);sd=a.std(axis=0,ddof=1);v=np.divide(a-a.mean(axis=0),sd,out=np.zeros_like(a),where=sd>1e-12);X=np.column_stack([np.ones(len(a)),v]);q,tri,piv=qr(X,mode='economic',pivoting=True);rank=int((np.abs(np.diag(tri))>1e-10*np.linalg.norm(X,2)).sum());assert rank==int(float(row.diagnostic_design_rank))
                fullrank=rank==9;lvary=sd[0]>1e-12;rvary=sd[1]>1e-12
                assert truth(row.diagnostic_design_full_rank)==fullrank and truth(row.ligand_varies)==lvary and truth(row.receptor_varies)==rvary
                if fullrank:assert np.isclose(float(row.diagnostic_design_condition),np.linalg.cond(X),rtol=1e-6,atol=1e-6)
                designs+=1
            status=('NO_EXACT_LIGAND_TARGET_PRIOR' if not truth(pr.exact_ligand_target_column) else 'LR_GENE_NOT_MEASURED' if not measured else 'INSUFFICIENT_PAIRED_LIBRARIES_OR_BLOCKS' if not (ok.sum()>=20 and nb>=10) else 'L_OR_R_CONSTANT' if not lvary or not rvary else 'NO_VARIABLE_NATIVE_TARGETS' if not var else 'DIAGNOSTIC_DESIGN_NOT_FULL_RANK' if not fullrank else 'COUNT_VARIATION_DIAGNOSTIC_DESIGN_SCREEN_PASSED_NOT_INFERENCE')
            assert row.status==status;checked+=1
        source_counts[sid]=int(d.gse_id.eq(sid).sum());del z,b,Xgenes,det;gc.collect()
    assert checked==len(d);result=dict(status='PASS_ALL336_MAPPING_ALL_ORIGINAL_RECEIVER_NODE_COUNT_TARGET_AND_NATIVE_CELL_DESIGN_DIAGNOSTICS',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(O/'inventory.json'),frozen_record_mapping=336,unique_interfaces=len(interfaces),checked_receiver_node_records=checked,cell_reconstructed_designs=designs,source_record_counts=source_counts,new_models_fitted=False,final_inference_accepted=False,scientific_release=False,axis_validation_complete=False)
    (O/'independent_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()

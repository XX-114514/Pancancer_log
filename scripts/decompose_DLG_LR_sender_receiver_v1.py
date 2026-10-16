"""All14 selected associations: descriptive sender/receiver RNA decomposition.

This post-selection diagnostic supplies no new p values or confirmatory FDR.
Conditional RNA associations do not establish causal direction or mediation.
"""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata

KEYS=['context','source','target','MCP','ligand','receptor']
BASE=['log_source_n','log_target_n','source_within_major','target_within_major','log_source_mean_RNA','log_target_mean_RNA','log_source_mean_detected','log_target_mean_detected']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rho(x,y,D=None):
    a=rankdata(x).astype(float);b=rankdata(y).astype(float);a-=a.mean();b-=b.mean()
    if D is not None:
        Q=np.linalg.qr(D,mode='reduced')[0];a-=Q@(Q.T@a);b-=Q@(Q.T@b)
    d=np.linalg.norm(a)*np.linalg.norm(b)
    return float(a@b/d) if d>1e-12 else np.nan


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);V=Path(json.loads((R/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    O=R/'DLG_LR_sender_receiver_decomposition_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
        scope='Every original14 BH-selected record, two fixed conditional RNA diagnostics each; retain reference-rejected records as flagged controls.',
        X='Source ligand or target receptor minimum across literal subunit mean(per-cell log1pCP10K), every native endpoint cell. Not surface protein or soluble concentration.',
        Y='Original frozen receiver MCP weighted components, remove all own ligand/receptor genes exactly as original.',
        nuisance='All original8 controls, plus each literal opposite-side subunit mean(per-cell log1pCP10K); no automatic control removal.',
        estimation='Rank X/Y, residualize against standardized nuisance; n>=20, full design rank, residualdf>=10, condition<1e8.',
        sensitivity='Complete leave-one-original-block-out, original interface conditional on receptor RNA, internal holdout with nuisance fitted only on discovery.',
        inference='Post-selection descriptive diagnostics only, no p/BH or causal/mediation claims, no new or exchanged endpoints.',
        inputs='Require independent all-native RNA panel audit and all4012 point/sensitivity/scope/BH replay before effects.')
    (O/'plan_before_conditional_RNA_diagnostics.json').write_text(json.dumps(plan,indent=2)+'\n')
    A=R/'DLG_LR_native_receptor_panel_audit_20261016_v1';receipt=R/'HNSCC_receiver_followup_job_control_20261016_v2/audit_DLG_LR_native_receptor_panels_v1/launch_receipt.json'
    launch=json.loads(receipt.read_text());deadline=time.monotonic()+7200
    while not (A/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists()
        st=p.read_text().rsplit(')',1)[1].split();assert st[19]==str(launch['proc_start_ticks']) and st[0]!='Z'
        assert time.monotonic()<deadline;time.sleep(15)
    assert json.loads((A/'private_terminal.json').read_text())['status'].startswith('PASS_ALL_NATIVE_PANEL_CELLS_')
    assert json.loads((R/'DLG_LR_point_binding_review_20261016_v1/private_terminal.json').read_text())['status'].startswith('PASS_ALL4012')
    original=pd.read_csv(V/'state_after_DIALOGUE_LR_response_20261016_v1/all_estimable_exploratory_LR_receiver_associations.tsv',sep='\t')
    original=original[original.BH_q.lt(.05)].sort_values(KEYS);assert len(original)==14
    N=R/'DLG_LR_native_receptor_panels_20261016_v3';expr=pd.read_csv(N/'private_all_native_unit_panel_expression.tsv',sep='\t').set_index(['unit_index','gene'])
    comp=pd.read_csv(V/'state_after_DIALOGUE_gene_components_20261016_v1/private_all_block_gene_components.tsv.gz',sep='\t')
    alignment=pd.read_csv(V/'state_after_DIALOGUE_all_LR_preflight_20261016_v1/all28_cache_alignment_checks.tsv',sep='\t').set_index(['context','source','target'])
    Q=V/'all_source_native_UCell_unit_sensitivity_20261014_v1';qi=json.loads((Q/'inventory.json').read_text());quality={}
    for info in qi['private_outputs']:
        p=Q/info['relative_file'];assert sha(p)==info['sha256'];z=np.load(p,allow_pickle=False)
        quality.update({int(k):(float(r),float(g)) for k,r,g in zip(z['unit_index'],z['mean_full_RNA'],z['mean_full_detected_genes'])})
    ref=pd.read_csv(R/'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv',sep='\t').set_index(['ligand','receptor'])
    result=[];loo=[];stored=O/'private_fixed_designs';stored.mkdir(mode=0o700)
    for number,r in enumerate(original.itertuples(),1):
        key={k:getattr(r,k) for k in KEYS};root=V/'state_after_DIALOGUE_inputs_20261016_v3'/f'context_{r.context:04d}'
        ci=json.loads((root/'private_input_inventory.json').read_text());part=pd.read_csv(root/'private_block_library_partition.tsv',sep='\t')
        reverse=r.source==ci['endpoints'][1];typ='typeA' if reverse else 'typeB'
        su=part['unit_index_B' if reverse else 'unit_index_A'].to_numpy(int);tu=part['unit_index_A' if reverse else 'unit_index_B'].to_numpy(int)
        ar=alignment.loc[(r.context,r.source,r.target)];z=np.load(R/'immune_same_name_family_source_review_v1/CM_edge_LR_extraction_v2/private_group_arrays'/(ar.old_group+'.npz'),allow_pickle=False)
        ix=pd.Index(z['samples'].astype(str)).get_indexer(part.sample_uid.astype(str));assert (ix>=0).all()
        j=list(zip(z['ligands'],z['receptors'])).index((r.ligand,r.receptor));interface=np.log1p(z['scores'][ix,j,0]);mask=np.isfinite(interface)
        assert int(mask.sum())==int(r.blocks) and list(z['patient_blocks'][ix])==list(part.effective_resampling_block)
        assert part.effective_resampling_block.is_unique
        qs=np.array([quality[k] for k in su]);qt=np.array([quality[k] for k in tu])
        base=np.column_stack([np.log1p(z['source_n'][ix]),np.log1p(z['target_n'][ix]),z['source_within_major'][ix],z['target_within_major'][ix],np.log1p(qs[:,0]),np.log1p(qt[:,0]),np.log1p(qs[:,1]),np.log1p(qt[:,1])])[mask]
        c=comp[comp.context.eq(r.context)&comp.type.eq(typ)&comp.program.eq(r.MCP)]
        B=c.pivot(index='sample',columns='gene',values='contribution').reindex(part.effective_resampling_block)
        weights=c.groupby('gene').weight.first().reindex(B.columns).to_numpy();remove=set(r.ligand.split('_')+r.receptor.split('_'))
        keep=~B.columns.isin(remove);Y=(B.loc[:,keep].to_numpy()@weights[keep])[mask]
        ls=r.ligand.split('_');rs=r.receptor.split('_')
        L=np.array([[expr.loc[(k,g),'mean_logCP10K'] for g in ls] for k in su])[mask]
        Rec=np.array([[expr.loc[(k,g),'mean_logCP10K'] for g in rs] for k in tu])[mask]
        assert np.isfinite(L).all() and np.isfinite(Rec).all()
        for us,expected in [(su,z['source_n'][ix]),(tu,z['target_n'][ix])]:
            assert np.array_equal(np.array([expr.loc[(k,ls[0]),'n_cells'] for k in us]),expected)
        baseD=np.column_stack([np.ones(len(Y)),(base-base.mean(0))/base.std(0)])
        assert np.isclose(rho(interface[mask],Y,baseD),r.adjusted_masked_rho,atol=1e-8)
        for mode,X,opposite,names in [('sender_L_given_receiver_R',L.min(1),Rec,['target_'+g for g in rs]),('receiver_R_given_sender_L',Rec.min(1),L,['source_'+g for g in ls])]:
            nuisance=np.column_stack([base,opposite]);sd=nuisance.std(0);n=len(Y)
            D=np.column_stack([np.ones(n),(nuisance-nuisance.mean(0))/np.where(sd>1e-12,sd,1)])
            rank=int(np.linalg.matrix_rank(D));condition=float(np.linalg.cond(D));df=n-rank-1
            row=dict(**key,mode=mode,blocks=n,original_BH_q=r.BH_q,original_interface_adjusted_rho=r.adjusted_masked_rho,
                reference_status=ref.loc[(r.ligand,r.receptor),'reference_status'],design_columns=D.shape[1],design_rank=rank,residual_df=df,design_condition=condition,
                extra_controls=';'.join(names),base_side_rho=rho(X,Y,baseD),conditional_rho=np.nan,interface_conditional_on_opposite_rho=np.nan,
                internal_holdout_rho=np.nan,LOO_estimable=0,LOO_same_sign=0,LOO_minimum=np.nan,LOO_maximum=np.nan,
                status='UNESTIMABLE',new_p_values=False,post_selection_diagnostic=True)
            identifier=f'R{number:02d}_'+mode
            tr=part.partition.eq('discovery').to_numpy()[mask];te=part.partition.eq('holdout').to_numpy()[mask]
            np.savez_compressed(stored/(identifier+'.npz'),X=X,Y=Y,nuisance=nuisance,interface=interface[mask],blocks=part.effective_resampling_block.to_numpy(str)[mask],train=tr,test=te)
            if n<20:row['status']='INSUFFICIENT_BLOCKS'
            elif (sd<=1e-12).any() or rank!=D.shape[1]:row['status']='CONSTANT_OR_COLLINEAR_FULL_CONTROLS'
            elif df<10 or condition>=1e8:row['status']='LOW_RESIDUAL_DF_OR_ILL_CONDITIONED'
            elif not np.isfinite(rho(X,Y,D)):row['status']='NO_RESIDUAL_X_OR_Y_VARIATION'
            else:
                row['status']='ESTIMABLE_POST_SELECTION_RNA_DESCRIPTION';row['conditional_rho']=rho(X,Y,D)
                row['interface_conditional_on_opposite_rho']=rho(interface[mask],Y,D)
                if tr.sum()>=20 and te.sum()>=10 and np.linalg.matrix_rank(D[tr])==D.shape[1]:
                    xh=X[te]-D[te]@np.linalg.lstsq(D[tr],X[tr],rcond=None)[0];yh=Y[te]-D[te]@np.linalg.lstsq(D[tr],Y[tr],rcond=None)[0]
                    row['internal_holdout_rho']=rho(xh,yh)
                values=[]
                for k in range(n):
                    use=np.arange(n)!=k;v=rho(X[use],Y[use],D[use]) if np.linalg.matrix_rank(D[use])==D.shape[1] else np.nan
                    loo.append(dict(record=identifier,removed_block_index=k,rho=v))
                    if np.isfinite(v):values.append(v)
                row.update(LOO_estimable=len(values),LOO_same_sign=int(np.sum(np.sign(values)==np.sign(row['conditional_rho']))),LOO_minimum=min(values) if values else np.nan,LOO_maximum=max(values) if values else np.nan)
            row['record']=identifier;result.append(row)
    assert len(result)==28
    result=pd.DataFrame(result);result.to_csv(O/'all28_private_sender_receiver_RNA_diagnostics.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(loo).to_csv(O/'private_all_estimable_leave_one_block_out.tsv',sep='\t',index=False,na_rep='NA')
    report=dict(status='ALL14_SELECTED_LR_TWO_CONDITIONAL_RNA_DIAGNOSTICS_COMPLETE_NOT_CONFIRMATORY_TESTS',finished_UTC=datetime.now(timezone.utc).isoformat(),records=len(result),
        status_counts=result.status.value_counts().to_dict(),native_panel_audit_sha256=sha(A/'private_terminal.json'),
        result_sha256=sha(O/'all28_private_sender_receiver_RNA_diagnostics.tsv'),new_p_values=False,independent_validation=False,causal_direction_established=False)
    (O/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    print('MIF_DIAGNOSTICS',result[result.ligand.eq('MIF')].to_json(orient='records'))


if __name__=='__main__':main()

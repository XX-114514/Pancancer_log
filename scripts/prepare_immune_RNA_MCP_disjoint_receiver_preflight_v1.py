"""All frozen candidate endpoints: disjoint readouts and design coverage, no effects."""
import gzip
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import rankdata

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from prepare_NSCLC_fixed_programme_continuous_priors_v1 import read_rds
from decompose_immune_RNA_MCP_original_nodes_v1 import node_baseline


def design_gate(x,y,C):
    n=len(x)
    if n<20:return 'FEWER_THAN20_COMPLETE_UNITS',np.nan,np.nan
    if not np.isfinite(np.c_[x,y,C]).all():return 'NONFINITE_DESIGN',np.nan,np.nan
    if min(np.ptp(x),np.ptp(y))==0:return 'CONSTANT_PREDICTOR_OR_READOUT',np.nan,np.nan
    r=rankdata(C,axis=0);s=r.std(axis=0,ddof=1)
    if np.any(s==0):return 'CONSTANT_REQUIRED_COVARIATE',np.nan,np.nan
    D=np.c_[np.ones(n),(r-r.mean(axis=0))/s]
    rank=int(np.linalg.matrix_rank(D));condition=float(np.linalg.cond(D));df=n-D.shape[1]-1
    if rank!=D.shape[1]:return 'COLLINEAR_REQUIRED_DESIGN',condition,df
    if condition>=1e8:return 'ILL_CONDITIONED_REQUIRED_DESIGN',condition,df
    if df<10:return 'FEWER_THAN10_FULL_MODEL_RESIDUAL_DF',condition,df
    Q=np.linalg.qr(D,mode='reduced')[0]
    for v in [x,y]:
        v=rankdata(v);residual=v-Q@(Q.T@v)
        if np.linalg.norm(residual)<1e-12:return 'SCORE_EXPLAINED_BY_REQUIRED_CONTROLS',condition,df
    return 'ESTIMABLE',condition,df


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,128*1024**2)
    qroot=root/'immune_RNA_MCP_all_ligand_qualification_20261018_v2'
    accepted=root/'immune_RNA_MCP_all_ligand_qualification_audit_20261018_v1/private_terminal.json'
    assert sha(accepted)=='ff2928361a945ec06bfd3e2827d149405351b54954d8b4406b90bc91eae5b659'
    proof=json.loads(accepted.read_text());assert sha(qroot/'private_terminal.json')==proof['producer_terminal_sha256']
    for name,item in json.loads((qroot/'private_terminal.json').read_text())['files'].items():assert sha(qroot/name)==item['sha256']
    out=root/'immune_RNA_MCP_disjoint_receiver_preflight_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        candidate_acceptance_sha256=sha(accepted),scope='All120 fixed selection slots; deduplicate programme/receiver/ligand/receptor/sender across signed-set selections, no replacement. All3 sources and66 case/sample units.',
        readout='Original weights/scaling/discovery nuisance transform, excluding ligand and all flat/reference receptors/components/physical aliases. Full=remaining+removed every cell and case. All3 source-centred Pearson>=.95 and RMSE/original source SD<=.25.',
        composition='Discovery38 only equal-case original-node baseline refit to each disjoint receiver; sender original programme baseline. Unseen nodes missing.',
        native='Limiting-component min(mean logCP10K) for ligand/receptor; source-measured true ligand zero retained. Receptor every-component detection>=1pct main,10pct coverage sensitivity. No assembled-complex claim.',
        controls='Both endpoint fraction/logN/cellQ/geneQ=8, two fixed node compositions, receptor limiting-component RNA=11 required. >=20 complete units, full rank/condition<1e8 and >=10 full residual df. No deletion.',
        no_new_communication_association_effect_P_or_FDR=True,identity_or_clinical_label_changes=False,
        privacy='Candidate/weights/case rows remain private. Only technical counts or authorized source aggregates may be reviewed.')
    (out/'plan_before_any_disjoint_readout_or_model_coverage.json').write_text(json.dumps(plan,indent=2)+'\n')
    shortlist=pd.read_csv(qroot/'private_frozen_top10_per_programme_view_direction.tsv',sep='\t');assert len(shortlist)==120
    edge=pd.read_csv(qroot/'private_all_reference_RNA_qualifications.tsv',sep='\t',keep_default_na=False)
    units=pd.read_csv(qroot/'private_all198_original_case_view_units.tsv',sep='\t')
    assays=np.load(qroot/'private_all198_native_ligand_receptor_assays.npz',allow_pickle=False)
    gene_index={g:i for i,g in enumerate(assays['genes'].tolist())}
    unit_index=pd.MultiIndex.from_frame(units[['fit_sample','view']]);assert unit_index.is_unique
    qmap=pd.read_csv(qroot/'private_all_prior_native_gene_columns.tsv',sep='\t').set_index('gene')
    sources=sorted(units.dataset.unique())
    flat=read_rds(root/'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1/verified_prior_inputs_v2/lr_network_human_21122021.rds')
    flat['from']=flat['from'].astype(str);flat['to']=flat['to'].astype(str)
    reference=pd.read_csv(root/'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv',sep='\t',keep_default_na=False)
    projection=root/'immune_RNA_MCP_frozen_projection_20261018_v1'
    assert sha(projection/'private_terminal.json')=='06c5637d9acf562a9266436da9f814af403fceed330a6bfd951ff45dc1d920a3'
    for name,digest in json.loads((projection/'private_terminal.json').read_text())['files'].items():assert sha(projection/name)==digest
    weights=pd.read_csv(projection/'private_R_projection/private_frozen_gene_weights.tsv',sep='\t')
    weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')]
    fullcell=pd.read_csv(projection/'private_R_projection/private_all_cell_score_reconstruction.tsv.gz',sep='\t')
    fullcell=fullcell[fullcell.fit_id.eq('K2_seed20261018_source_adjusted')].set_index(['cell','program']).score
    inputs=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1'
    assert sha(inputs/'private_terminal.json')=='2474355a003d96232b752ad2c77ba92a506f4c68829c7e0adf434045cf03d656'
    for name,item in json.loads((inputs/'private_terminal.json').read_text())['files'].items():assert sha(inputs/name)==item['sha256']
    selected=pd.read_csv(inputs/'private_frozen_original_cell_selection.tsv',sep='\t')
    old=root/'existing_immune_fine_native_measurement_scope_20261016_v1'
    assert sha(old/'private_terminal.json')==json.loads((inputs/'plan_before_selected_native_expression_read.json').read_text())['raw_original_cell_scope_sha256']
    for name,digest in json.loads((old/'private_terminal.json').read_text())['files'].items():assert sha(old/name)==digest
    oldunits=pd.read_csv(old/'private_all_frozen_T_myeloid_node_units.tsv',sep='\t').set_index('unit_index')
    oldmap=np.load(old/'private_original_cell_to_fine_node_unit.npy',mmap_mode='r')
    selected['node_id']=oldunits.node_id.loc[oldmap[selected.global_row.to_numpy(int)]].to_numpy()
    nodes=selected.set_index('cell').node_id
    native=pd.read_csv(root/'immune_RNA_MCP_discovery_scope_20261018_v1/private_marker_excluded_native_gene_map.tsv',sep='\t').set_index('gene')
    metadata=pd.read_csv(root/'immune_RNA_MCP_source_coordination_20261018_v1/private_all_case_view_controls.tsv',sep='\t')
    coordination=root/'immune_RNA_MCP_source_coordination_20261018_v1'
    assert sha(coordination/'private_terminal.json')=='72fdf3928ac53318e2c5b4db162d4e9c3c432b76f5105191e463a02ee73d7b23'
    control_name='private_all_case_view_controls.tsv'
    assert sha(coordination/control_name)==json.loads((coordination/'private_terminal.json').read_text())['files'][control_name]
    decomposition=root/'immune_RNA_MCP_original_node_decomposition_20261018_v1'
    assert sha(decomposition/'private_terminal.json')=='785f3a4cf7eba7c56e791abfddfd8040de7b9cc93dc2506a1905e00c9f745f5d'
    baseline_name='private_all396_case_view_decompositions.tsv'
    assert sha(decomposition/baseline_name)==json.loads((decomposition/'private_terminal.json').read_text())['files'][baseline_name]
    sender_baseline=pd.read_csv(root/'immune_RNA_MCP_original_node_decomposition_20261018_v1/private_all396_case_view_decompositions.tsv',sep='\t').set_index(['sample','view','program']).composition_prediction
    proto=[];fidelity=[];parameters=[];case_scores=[];cell_checks=[];baseline_proofs=[]
    for view,sv in shortlist.groupby('view',sort=True):
        data={}
        for part in ['discovery','TNBC_holdout']:
            prefix=inputs/'private_inputs'/(view+'_'+part)
            md=pd.read_csv(str(prefix)+'_metadata.tsv',sep='\t')
            genes=Path(str(prefix)+'_genes.tsv').read_text().splitlines()
            with gzip.open(str(prefix)+'_tpm.mtx.gz','rb') as handle:matrix=mmread(handle).tocsr()
            data[part]=(md,matrix,genes)
        td,tm,genes=data['discovery'];hd,hm,hgenes=data['TNBC_holdout'];assert genes==hgenes
        D=np.c_[np.ones(len(td)),td[['cellQ','geneQ','sourceCC','fraction','logN']].to_numpy(float)]
        Dh=np.c_[np.ones(len(hd)),hd[['cellQ','geneQ','fraction','logN']].to_numpy(float)];ci=[0,1,2,4,5]
        for program,sp in sv.groupby('program',sort=True):
            w=weights[weights.view.eq(view)&weights.program.eq(program)].copy()
            positions=pd.Index(genes).get_indexer(w.gene);assert (positions>=0).all()
            v=w.weight.to_numpy()/w.discovery_sd.to_numpy();mu=w.discovery_mean.to_numpy()
            X=tm[positions];Xh=hm[positions]
            total_raw=np.asarray(v@X).ravel()-v@mu;total_coef=np.linalg.lstsq(D,total_raw,rcond=None)[0]
            total_scores={'discovery':total_raw-D@total_coef,'TNBC_holdout':np.asarray(v@Xh).ravel()-v@mu-Dh@total_coef[ci]}
            for part,(md,_,_) in data.items():
                expected=fullcell.loc[pd.MultiIndex.from_arrays([md.cell,np.repeat(program,len(md))])].to_numpy()
                assert np.max(np.abs(total_scores[part]-expected))<1e-8
            physical=native.loc[w.gene,sources].to_numpy(int)
            for ligand,sl in sp.groupby('ligand',sort=True):
                pid=len(proto);directions='|'.join(sorted(sl.direction.unique()))
                receptor_names=set(flat.loc[flat['from'].eq(ligand),'to'])|set(reference.loc[reference.ligand.eq(ligand),'receptor'])
                symbols=set(str(ligand).split('_'))
                for name in receptor_names:symbols.update(str(name).split('_'))
                excluded=np.zeros(len(w),bool)
                for si,source in enumerate(sources):
                    columns=qmap.loc[sorted(symbols),source].to_numpy(int);columns=columns[columns>=0]
                    excluded|=np.isin(physical[:,si],columns)
                remain=v.copy();remain[excluded]=0
                raw=np.asarray(remain@X).ravel()-remain@mu;coef=np.linalg.lstsq(D,raw,rcond=None)[0]
                values={'discovery':raw-D@coef,'TNBC_holdout':np.asarray(remain@Xh).ravel()-remain@mu-Dh@coef[ci]}
                allcells=[]
                for part,(md,matrix,_) in data.items():
                    removed=v-remain;R=np.asarray(removed@(X if part=='discovery' else Xh)).ravel()-removed@mu
                    Rc=np.linalg.lstsq(D,np.asarray(removed@X).ravel()-removed@mu,rcond=None)[0]
                    R-=D@Rc if part=='discovery' else Dh@Rc[ci]
                    error=float(np.max(np.abs(values[part]+R-total_scores[part])));assert error<1e-8
                    cell_checks.append(dict(prototype_id=pid,partition=part,cells=len(md),max_complement_error=error))
                    f=md[['sample','dataset','cell']].copy();f['score']=values[part];f['original_score']=total_scores[part];f['node_id']=nodes.loc[md.cell].to_numpy();allcells.append(f)
                cf=pd.concat(allcells,ignore_index=True)
                means=cf.groupby(['sample','dataset','node_id']).agg(mean_score=('score','mean'),selected_cells=('score','size')).reset_index()
                beta,proof=node_baseline(means);baseline_proofs.append(dict(prototype_id=pid,**proof))
                if beta is not None:parameters.append(beta.assign(prototype_id=pid));beta=beta.set_index('node_id').coefficient
                rows=[]
                for (sample,source),g in cf.groupby(['sample','dataset'],sort=True):
                    contribution=g.node_id.map(beta) if beta is not None else pd.Series(np.nan,index=g.index)
                    rows.append(dict(prototype_id=pid,sample=sample,dataset=source,receiver_view=view,program=program,
                        receiver_score=float(g.score.mean()),original_receiver_score=float(g.original_score.mean()),
                        receiver_composition=float(contribution.mean()) if contribution.notna().all() else np.nan,
                        unseen_node_cells=int(contribution.isna().sum())))
                cs=pd.DataFrame(rows);assert len(cs)==66
                case_scores.append(cs);passes=[]
                for source,g in cs.groupby('dataset',sort=True):
                    a=g.receiver_score.to_numpy(copy=True);b=g.original_receiver_score.to_numpy(copy=True);a-=a.mean();b-=b.mean()
                    correlation=float(a@b/(np.linalg.norm(a)*np.linalg.norm(b))) if min(np.std(a),np.std(b))>1e-12 else np.nan
                    rmse=float(np.sqrt(np.mean((a-b)**2))/np.std(b)) if np.std(b)>1e-12 else np.nan
                    ok=np.isfinite(correlation) and correlation>=.95 and rmse<=.25;passes.append(ok)
                    fidelity.append(dict(prototype_id=pid,dataset=source,units=len(g),Pearson=correlation,centred_RMSE_over_original_SD=rmse,faithful=bool(ok)))
                proto.append(dict(prototype_id=pid,receiver_view=view,program=program,ligand=ligand,selection_directions=directions,
                    selection_slots=len(sl),excluded_active_weight_genes=int(excluded.sum()),all3_source_fidelity_pass=bool(all(passes))))
    prototypes=pd.DataFrame(proto);case_score=pd.concat(case_scores,ignore_index=True);assert len(case_score)==66*len(prototypes)
    slots=[]
    for r in prototypes.itertuples(index=False):
        edges=edge[edge.receiver_view.eq(r.receiver_view)&edge.ligand.eq(r.ligand)&edge.candidate_qualified]
        for e in edges.itertuples(index=False):
            for sender in str(e.eligible_senders).split('|'):
                assert sender and sender!=r.receiver_view
                slots.append(dict(prototype_id=r.prototype_id,program=r.program,receiver_view=r.receiver_view,ligand=r.ligand,receptor=e.receptor,
                    sender_view=sender,reference_status=e.reference_status,selection_directions=r.selection_directions,all3_source_fidelity_pass=r.all3_source_fidelity_pass))
    registry=pd.DataFrame(slots).drop_duplicates(['program','receiver_view','ligand','receptor','sender_view']).sort_values(['program','receiver_view','ligand','receptor','sender_view']).reset_index(drop=True)
    registry['slot']=np.arange(len(registry));assert len(registry)>0
    all_inputs=[];coverage=[];controls=[c+'_'+side for side in ['S','R'] for c in ['fraction','logN','cellQ','geneQ']]+['sender_composition','receiver_composition','receptor_RNA']
    for r in registry.itertuples(index=False):
        ms=metadata[metadata.view.eq(r.sender_view)].drop(columns='view');mr=metadata[metadata.view.eq(r.receiver_view)].drop(columns='view')
        d=ms.merge(mr,on=['sample','dataset','partition','reviewed_patient'],suffixes=('_S','_R'),validate='one_to_one')
        cs=case_score[case_score.prototype_id.eq(r.prototype_id)]
        d=d.merge(cs[['sample','receiver_score','receiver_composition','unseen_node_cells']],on='sample',validate='one_to_one')
        d['sender_composition']=sender_baseline.loc[pd.MultiIndex.from_arrays([d['sample'],np.repeat(r.sender_view,len(d)),np.repeat(r.program,len(d))])].to_numpy()
        si=unit_index.get_indexer(pd.MultiIndex.from_arrays([d['sample'],np.repeat(r.sender_view,len(d))]));ri=unit_index.get_indexer(pd.MultiIndex.from_arrays([d['sample'],np.repeat(r.receiver_view,len(d))]));assert min(si.min(),ri.min())>=0
        lc=[gene_index[g] for g in str(r.ligand).split('_')];rc=[gene_index[g] for g in str(r.receptor).split('_')]
        d['ligand_RNA']=assays['mean_logCP10K'][np.ix_(si,lc)].min(axis=1);d['receptor_RNA']=assays['mean_logCP10K'][np.ix_(ri,rc)].min(axis=1)
        rh=assays['positive_cells'][np.ix_(ri,rc)];rn=assays['cells'][ri]
        d['receptor_support1pct']=(rh*100>=rn[:,None]).all(axis=1);d['receptor_support10pct']=(rh*10>=rn[:,None]).all(axis=1)
        d['slot']=r.slot;d['prototype_id']=r.prototype_id
        d['complete_design']=np.isfinite(d[['ligand_RNA','receiver_score']+controls].to_numpy(float)).all(axis=1)
        d['primary_eligible']=d.receptor_support1pct&d.complete_design
        all_inputs.append(d[['slot','prototype_id','sample','dataset','partition','reviewed_patient','ligand_RNA','receiver_score']+controls+['receptor_support1pct','receptor_support10pct','complete_design','primary_eligible','unseen_node_cells']])
        for source,g in d.groupby('dataset',sort=True):
            eligible=g[g.primary_eligible]
            reason,condition,df=design_gate(eligible.ligand_RNA.to_numpy(),eligible.receiver_score.to_numpy(),eligible[controls].to_numpy(float))
            status=reason if r.all3_source_fidelity_pass else 'DISJOINT_RECEIVER_NOT_FAITHFUL'
            coverage.append(dict(slot=r.slot,prototype_id=r.prototype_id,dataset=source,source_units=len(g),complete_units=int(g.complete_design.sum()),
                receptor_supported_units=int(g.receptor_support1pct.sum()),receptor_supported10pct_units=int(g.receptor_support10pct.sum()),
                primary_units=len(eligible),zero_ligand_units=int(eligible.ligand_RNA.eq(0).sum()),reviewed_patient_units=int(eligible.reviewed_patient.sum()),
                status=status,design_status_before_fidelity=reason,design_condition=condition,residual_df=df,no_communication_effect_calculated=True))
    prototypes.to_csv(out/'private_all_disjoint_receiver_prototypes.tsv',sep='\t',index=False)
    registry.to_csv(out/'private_all_unique_directional_endpoint_slots.tsv',sep='\t',index=False)
    case_score.to_csv(out/'private_all_case_disjoint_receiver_scores.tsv.gz',sep='\t',index=False)
    pd.DataFrame(fidelity).to_csv(out/'private_all_source_readout_fidelity.tsv',sep='\t',index=False)
    pd.concat(parameters,ignore_index=True).to_csv(out/'private_disjoint_discovery_node_coefficients.tsv',sep='\t',index=False)
    pd.DataFrame(baseline_proofs).to_csv(out/'private_all_disjoint_node_baseline_designs.tsv',sep='\t',index=False)
    pd.DataFrame(cell_checks).to_csv(out/'private_all_cell_complement_checks.tsv',sep='\t',index=False)
    pd.concat(all_inputs,ignore_index=True).to_csv(out/'private_all_slot_case_model_inputs.tsv.gz',sep='\t',index=False)
    cov=pd.DataFrame(coverage);assert len(cov)==3*len(registry)
    cov.to_csv(out/'private_all_source_model_coverage.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()};assert sum(p['bytes'] for p in files.values())<128*1024**2
    (out/'private_terminal.json').write_text(json.dumps(dict(status='COMPLETE_ALL_DISJOINT_READOUTS_AND_ENDPOINT_COVERAGE_INDEPENDENT_ACCEPTANCE_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,
        frozen_selection_slots=120,receiver_prototypes=len(prototypes),unique_directional_endpoint_slots=len(registry),source_coverage_positions=len(cov),
        faithful_prototypes=int(prototypes.all3_source_fidelity_pass.sum()),estimable_TNBC_slots=int((cov.dataset.eq('GSE305078')&cov.status.eq('ESTIMABLE')).sum()),
        no_new_communication_effect_P_or_FDR=True,no_identity_or_clinical_label_change=True),indent=2)+'\n')
    print('ALL_DISJOINT_READOUTS_AND_ENDPOINT_COVERAGE_COMPLETE_ACCEPTANCE_PENDING')


if __name__=='__main__':main()

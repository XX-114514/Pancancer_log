"""Independent QR reconstruction, endpoint inventory and every design gate."""
import gzip
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from prepare_NSCLC_fixed_programme_continuous_priors_v1 import read_rds


def same(a,b):
    assert (pd.isna(a) and pd.isna(b)) or (np.isfinite(a) and np.isfinite(b) and abs(a-b)<1e-8)


def gate(x,y,C):
    n=len(x)
    if n<20:return 'FEWER_THAN20_COMPLETE_UNITS',np.nan,np.nan
    if not np.isfinite(np.column_stack([x,y,C])).all():return 'NONFINITE_DESIGN',np.nan,np.nan
    if min(np.ptp(x),np.ptp(y))==0:return 'CONSTANT_PREDICTOR_OR_READOUT',np.nan,np.nan
    ranks=pd.DataFrame(C).rank().to_numpy();sd=ranks.std(axis=0,ddof=1)
    if np.any(sd==0):return 'CONSTANT_REQUIRED_COVARIATE',np.nan,np.nan
    D=np.column_stack([np.ones(n),(ranks-ranks.mean(axis=0))/sd])
    values=np.linalg.svd(D,compute_uv=False);rank=np.sum(values>values[0]*max(D.shape)*np.finfo(float).eps)
    condition=float(values[0]/values[-1]);df=n-D.shape[1]-1
    if rank!=D.shape[1]:return 'COLLINEAR_REQUIRED_DESIGN',condition,df
    if condition>=1e8:return 'ILL_CONDITIONED_REQUIRED_DESIGN',condition,df
    if df<10:return 'FEWER_THAN10_FULL_MODEL_RESIDUAL_DF',condition,df
    for v in [x,y]:
        a=pd.Series(v).rank().to_numpy();residual=a-D@np.linalg.lstsq(D,a,rcond=None)[0]
        if np.linalg.norm(residual)<1e-12:return 'SCORE_EXPLAINED_BY_REQUIRED_CONTROLS',condition,df
    return 'ESTIMABLE',condition,df


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,8*1024**2)
    source=root/'immune_RNA_MCP_disjoint_receiver_preflight_20261018_v1'
    terminal=json.loads((source/'private_terminal.json').read_text())
    assert terminal['source_sha256']=='c49b5ff607053017892c5dfb18e7170b9723d117c77b219345bf162d02fa6ba9'
    for name,item in terminal['files'].items():assert sha(source/name)==item['sha256']
    out=root/'immune_RNA_MCP_disjoint_receiver_preflight_audit_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_all_readout_and_design_replay.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        producer_terminal_sha256=sha(source/'private_terminal.json'),
        scope='All frozen selection/prototype/explicit receptor-sender positions, all disjoint case scores from complete candidate-cell matrices, all discovery-node coefficients by independent full QR, all source fidelities, every case model input and all source design gates. No new ligand association effect.'),indent=2)+'\n')
    qroot=root/'immune_RNA_MCP_all_ligand_qualification_20261018_v2'
    assert sha(root/'immune_RNA_MCP_all_ligand_qualification_audit_20261018_v1/private_terminal.json')=='ff2928361a945ec06bfd3e2827d149405351b54954d8b4406b90bc91eae5b659'
    shortlist=pd.read_csv(qroot/'private_frozen_top10_per_programme_view_direction.tsv',sep='\t')
    edges=pd.read_csv(qroot/'private_all_reference_RNA_qualifications.tsv',sep='\t',keep_default_na=False)
    gene_columns=pd.read_csv(qroot/'private_all_prior_native_gene_columns.tsv',sep='\t').set_index('gene')
    native=pd.read_csv(root/'immune_RNA_MCP_discovery_scope_20261018_v1/private_marker_excluded_native_gene_map.tsv',sep='\t').set_index('gene')
    prototype=pd.read_csv(source/'private_all_disjoint_receiver_prototypes.tsv',sep='\t')
    registry=pd.read_csv(source/'private_all_unique_directional_endpoint_slots.tsv',sep='\t')
    keys=['receiver_view','program','ligand']
    expected=set(map(tuple,shortlist.rename(columns={'view':'receiver_view'})[keys].to_numpy()))
    assert expected==set(map(tuple,prototype[keys].to_numpy())) and len(prototype)==len(expected)
    actual=set(map(tuple,registry[['program','receiver_view','ligand','receptor','sender_view']].to_numpy()));expected=set()
    for p in prototype.itertuples(index=False):
        for e in edges[edges.receiver_view.eq(p.receiver_view)&edges.ligand.eq(p.ligand)&edges.candidate_qualified].itertuples(index=False):
            for sender in e.eligible_senders.split('|'):expected.add((p.program,p.receiver_view,p.ligand,e.receptor,sender))
        selections=shortlist[shortlist.view.eq(p.receiver_view)&shortlist.program.eq(p.program)&shortlist.ligand.eq(p.ligand)]
        assert p.selection_slots==len(selections) and p.selection_directions=='|'.join(sorted(selections.direction.unique()))
    assert actual==expected and len(actual)==len(registry)
    assert registry.slot.tolist()==list(range(len(registry)))
    flatpath=root/'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1/verified_prior_inputs_v2/lr_network_human_21122021.rds'
    assert sha(flatpath)=='47c971d2fbba4ecd0ba7485d1846a74054432a0d1a97ea3e6a79ae27d0da8094'
    flat=read_rds(flatpath);flat['from']=flat['from'].astype(str);flat['to']=flat['to'].astype(str)
    refpath=root/'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv'
    assert sha(refpath)=='df574533d866438e646423f29bd555cb283173dbceadfe9f8264e1622a86fd06'
    reference=pd.read_csv(refpath,sep='\t',keep_default_na=False)
    projection=root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection'
    weights=pd.read_csv(projection/'private_frozen_gene_weights.tsv',sep='\t');weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')]
    oldcases=pd.read_csv(projection/'private_all_case_mean_scores.tsv',sep='\t');oldcases=oldcases[oldcases.fit_id.eq('K2_seed20261018_source_adjusted')].set_index(['sample','view','program']).score
    inputs=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1'
    selected=pd.read_csv(inputs/'private_frozen_original_cell_selection.tsv',sep='\t')
    old=root/'existing_immune_fine_native_measurement_scope_20261016_v1'
    units=pd.read_csv(old/'private_all_frozen_T_myeloid_node_units.tsv',sep='\t').set_index('unit_index')
    mapping=np.load(old/'private_original_cell_to_fine_node_unit.npy',mmap_mode='r')
    selected['node_id']=units.node_id.loc[mapping[selected.global_row.to_numpy(int)]].to_numpy();cell_node=selected.set_index('cell').node_id
    saved_cases=pd.read_csv(source/'private_all_case_disjoint_receiver_scores.tsv.gz',sep='\t').set_index(['prototype_id','sample'])
    saved_beta=pd.read_csv(source/'private_disjoint_discovery_node_coefficients.tsv',sep='\t').set_index(['prototype_id','node_id'])
    saved_fidelity=pd.read_csv(source/'private_all_source_readout_fidelity.tsv',sep='\t').set_index(['prototype_id','dataset'])
    sources=sorted(gene_columns.columns);reconstructed=[];max_case_error=0.;node_parameter_checks=0
    for view,group in prototype.groupby('receiver_view',sort=True):
        matrices={}
        for part in ['discovery','TNBC_holdout']:
            prefix=inputs/'private_inputs'/(view+'_'+part);md=pd.read_csv(str(prefix)+'_metadata.tsv',sep='\t')
            genes=Path(str(prefix)+'_genes.tsv').read_text().splitlines()
            with gzip.open(str(prefix)+'_tpm.mtx.gz','rb') as f:matrix=mmread(f).tocsc()
            matrices[part]=(md,matrix,genes)
        td,tm,names=matrices['discovery'];hd,hm,hnames=matrices['TNBC_holdout'];assert names==hnames
        D=np.column_stack([np.ones(len(td)),td[['cellQ','geneQ','sourceCC','fraction','logN']].to_numpy()]);Q,R=np.linalg.qr(D,mode='reduced')
        Dh=np.column_stack([np.ones(len(hd)),hd[['cellQ','geneQ','fraction','logN']].to_numpy()])
        for p in group.itertuples(index=False):
            w=weights[weights.view.eq(view)&weights.program.eq(p.program)].copy();where=pd.Index(names).get_indexer(w.gene)
            all_receptors=set(flat.loc[flat['from'].eq(p.ligand),'to'])|set(reference.loc[reference.ligand.eq(p.ligand),'receptor'])
            symbols=set(g for s in {p.ligand,*all_receptors} for g in str(s).split('_'))
            keep=np.ones(len(w),bool)
            for s in sources:
                excluded=set(int(v) for v in gene_columns.loc[sorted(symbols),s] if v>=0)
                keep&=~np.isin(native.loc[w.gene,s].to_numpy(int),list(excluded))
            assert p.excluded_active_weight_genes==int((~keep).sum())
            b=w.weight.to_numpy()/w.discovery_sd.to_numpy();b[~keep]=0
            raw=tm[where].T.dot(b)-b@w.discovery_mean.to_numpy()
            coef=np.linalg.solve(R,Q.T@raw);values={'discovery':raw-Q@(Q.T@raw),
                'TNBC_holdout':hm[where].T.dot(b)-b@w.discovery_mean.to_numpy()-Dh@coef[[0,1,2,4,5]]}
            cell=[]
            for part,(md,_,_) in matrices.items():
                f=md[['sample','dataset','cell']].copy();f['score']=values[part];f['node_id']=cell_node.loc[md.cell].to_numpy();cell.append(f)
            cell=pd.concat(cell,ignore_index=True);m=cell.groupby(['sample','dataset','node_id']).agg(mean_score=('score','mean'),n=('score','size')).reset_index()
            train=m[~m.dataset.eq('GSE305078')].copy();cases=sorted(train['sample'].unique());node_ids=sorted(train.node_id.unique())
            reference_node=train.groupby('node_id').n.sum().sort_values(ascending=False,kind='stable').index[0];other=[n for n in node_ids if n!=reference_node]
            weight=train.n.to_numpy()/train.groupby('sample').n.transform('sum').to_numpy()
            full=np.column_stack([pd.get_dummies(pd.Categorical(train['sample'],categories=cases),dtype=float),
                pd.get_dummies(pd.Categorical(train.node_id,categories=other),dtype=float)])*np.sqrt(weight[:,None])
            assert np.linalg.matrix_rank(full)==full.shape[1];q,r=np.linalg.qr(full,mode='reduced')
            effect=np.linalg.solve(r,q.T@(train.mean_score.to_numpy()*np.sqrt(weight)))
            beta=pd.Series(0.,index=node_ids);beta.loc[other]=effect[len(cases):]
            beta-=float(np.sum(weight*train.node_id.map(beta))/38)
            saved=saved_beta.loc[p.prototype_id];assert set(saved.index)==set(beta.index)
            assert np.allclose(beta.loc[saved.index],saved.coefficient,atol=1e-8,rtol=0);node_parameter_checks+=len(beta)
            current=[]
            for (sample,dataset),g in cell.groupby(['sample','dataset'],sort=True):
                contribution=g.node_id.map(beta);score=float(g.score.mean());comp=float(contribution.mean()) if contribution.notna().all() else np.nan
                saved=saved_cases.loc[(p.prototype_id,sample)];same(score,saved.receiver_score);same(comp,saved.receiver_composition)
                assert saved.unseen_node_cells==contribution.isna().sum()
                original=float(oldcases.loc[(sample,view,p.program)]);same(original,saved.original_receiver_score)
                max_case_error=max(max_case_error,abs(score-saved.receiver_score))
                current.append(dict(prototype_id=p.prototype_id,sample=sample,dataset=dataset,receiver_score=score,receiver_composition=comp,original=original))
            current=pd.DataFrame(current);passes=[]
            for dataset,g in current.groupby('dataset'):
                a=g.receiver_score.to_numpy(copy=True);b=g.original.to_numpy(copy=True);a-=a.mean();b-=b.mean()
                corr=float(np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b))) if min(np.std(a),np.std(b))>1e-12 else np.nan
                error=float(np.sqrt(np.average((a-b)**2))/np.std(b)) if np.std(b)>1e-12 else np.nan
                saved=saved_fidelity.loc[(p.prototype_id,dataset)];same(corr,saved.Pearson);same(error,saved.centred_RMSE_over_original_SD)
                ok=np.isfinite(corr) and corr>=.95 and error<=.25;assert ok==saved.faithful;passes.append(ok)
            assert p.all3_source_fidelity_pass==all(passes);reconstructed.append(current)
    newcases=pd.concat(reconstructed,ignore_index=True).set_index(['prototype_id','sample'])
    data=pd.read_csv(source/'private_all_slot_case_model_inputs.tsv.gz',sep='\t')
    coverage=pd.read_csv(source/'private_all_source_model_coverage.tsv',sep='\t').set_index(['slot','dataset'])
    meta=pd.read_csv(root/'immune_RNA_MCP_source_coordination_20261018_v1/private_all_case_view_controls.tsv',sep='\t').set_index(['sample','view'])
    oldcomp=pd.read_csv(root/'immune_RNA_MCP_original_node_decomposition_20261018_v1/private_all396_case_view_decompositions.tsv',sep='\t').set_index(['sample','view','program']).composition_prediction
    nativeunits=pd.read_csv(qroot/'private_all198_original_case_view_units.tsv',sep='\t')
    native_index=pd.MultiIndex.from_frame(nativeunits[['fit_sample','view']]);native_assay=np.load(qroot/'private_all198_native_ligand_receptor_assays.npz',allow_pickle=False)
    gi={g:i for i,g in enumerate(native_assay['genes'].tolist())}
    controls=[c+'_'+side for side in ['S','R'] for c in ['fraction','logN','cellQ','geneQ']]+['sender_composition','receiver_composition','receptor_RNA']
    source_checks=0
    for slot in registry.itertuples(index=False):
        d=data[data.slot.eq(slot.slot)].copy();assert len(d)==66 and d['sample'].is_unique
        for r in d.itertuples(index=False):
            sc=newcases.loc[(slot.prototype_id,r.sample)];same(r.receiver_score,sc.receiver_score);same(r.receiver_composition,sc.receiver_composition)
            for view,side in [(slot.sender_view,'S'),(slot.receiver_view,'R')]:
                md=meta.loc[(r.sample,view)];assert md.dataset==r.dataset and md.reviewed_patient==r.reviewed_patient
                for c in ['fraction','logN','cellQ','geneQ']:same(getattr(r,c+'_'+side),md[c])
            same(r.sender_composition,oldcomp.loc[(r.sample,slot.sender_view,slot.program)])
            si=native_index.get_loc((r.sample,slot.sender_view));ri=native_index.get_loc((r.sample,slot.receiver_view))
            lc=[gi[g] for g in slot.ligand.split('_')];rc=[gi[g] for g in slot.receptor.split('_')]
            same(r.ligand_RNA,float(np.min(native_assay['mean_logCP10K'][si,lc])));same(r.receptor_RNA,float(np.min(native_assay['mean_logCP10K'][ri,rc])))
            hits=native_assay['positive_cells'][ri,rc];n=native_assay['cells'][ri]
            assert r.receptor_support1pct==bool(np.all(hits*100>=n)) and r.receptor_support10pct==bool(np.all(hits*10>=n))
        finite=np.isfinite(d[['ligand_RNA','receiver_score']+controls].to_numpy(float)).all(axis=1)
        assert np.array_equal(finite,d.complete_design) and np.array_equal(finite&d.receptor_support1pct,d.primary_eligible)
        for dataset,g in d.groupby('dataset'):
            saved=coverage.loc[(slot.slot,dataset)];eligible=g[g.primary_eligible]
            reason,condition,df=gate(eligible.ligand_RNA.to_numpy(),eligible.receiver_score.to_numpy(),eligible[controls].to_numpy(float))
            assert saved.design_status_before_fidelity==reason and saved.status==(reason if slot.all3_source_fidelity_pass else 'DISJOINT_RECEIVER_NOT_FAITHFUL')
            same(condition,saved.design_condition);same(df,saved.residual_df)
            assert saved.primary_units==len(eligible) and saved.source_units==len(g) and saved.reviewed_patient_units==eligible.reviewed_patient.sum()
            assert saved.complete_units==g.complete_design.sum() and saved.receptor_supported_units==g.receptor_support1pct.sum()
            assert saved.receptor_supported10pct_units==g.receptor_support10pct.sum() and saved.zero_ligand_units==eligible.ligand_RNA.eq(0).sum()
            source_checks+=1
    assert source_checks==3*len(registry) and len(data)==66*len(registry)
    (out/'private_terminal.json').write_text(json.dumps(dict(status='ACCEPTED_ALL_DISJOINT_READOUTS_NODE_BASELINES_ENDPOINTS_AND_COVERAGE',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),
        receiver_prototypes=len(prototype),unique_endpoint_slots=len(registry),source_designs_replayed=source_checks,
        case_model_rows_replayed=len(data),node_parameter_checks=node_parameter_checks,max_case_score_error=max_case_error,
        no_new_ligand_association_effect_P_or_FDR=True),indent=2)+'\n')
    print('ALL_DISJOINT_READOUTS_AND_ENDPOINT_DESIGNS_INDEPENDENTLY_ACCEPTED')


if __name__=='__main__':main()

"""All frozen own-view vectors vs unsigned prior, no external effect selection."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadr
from scipy.stats import rankdata

from NSCLC_native_MFA_interpretation_inputs_v1 import load, sha, VIEWS
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

PRIORS={'lr_network_human_21122021.rds':'47c971d2fbba4ecd0ba7485d1846a74054432a0d1a97ea3e6a79ae27d0da8094',
        'ligand_target_matrix_nsga2r_final.rds':'699fce17ff65e2511277359696306bb130aa833bd7f30dd9290ad9efd9dc9c5d'}
ACCEPTED={'simplex_reference_supported','all_constituents_supported_complex_assembly_unverified'}


def read_rds(path):
    result=pyreadr.read_r(str(path));assert len(result)==1
    return next(iter(result.values()))


def calibrate(prediction,original,train):
    scale=original[train].std(axis=0);values={}
    for name,take in [('discovery',train),('original_holdout',~train)]:
        a=prediction[take]-prediction[take].mean(axis=0);b=original[take]-original[take].mean(axis=0)
        with np.errstate(divide='ignore',invalid='ignore'):
            values[name+'_rho']=np.sum(a*b,axis=0)/(np.linalg.norm(a,axis=0)*np.linalg.norm(b,axis=0))
            values[name+'_RMSE_over_discovery_SD']=np.sqrt(np.mean((prediction[take]-original[take])**2,axis=0))/scale
    gate=np.ones(10,bool)
    for key,value in values.items():gate &= np.isfinite(value)&(value>=.95 if key.endswith('_rho') else value<=.25)
    return values,gate


def main():
    data=load();root=data['root'];budget=resource(root,64*1024**2)
    prior=root/'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1/verified_prior_inputs_v2'
    for name,digest in PRIORS.items():assert sha(prior/name)==digest
    refpath=root/'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv'
    assert sha(refpath)=='df574533d866438e646423f29bd555cb283173dbceadfe9f8264e1622a86fd06'
    positive_path=data['native']/'positive_cell_counts.npy'
    native_terminal=json.loads((data['native']/'private_terminal.json').read_text())
    assert sha(positive_path)==native_terminal['outputs']['positive_cell_counts.npy']['sha256']
    out=root/'NSCLC_fixed_programme_continuous_prior_inputs_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),budget=budget,
              scope='All10 frozen factors x4 own views x1226 prior ligands; all80 signed top20 slots retained. No external coordination/clinical or function-based factor selection.',
              background='Original selected features explicitly rechecked >=100 raw counts in >=ceil(.25*145) discovery cases, exact prior target symbols. No extra cell-fraction filter that would remove rare-state genes. Contributions outside prior background omitted, never missing-RNA zero fill.',
              coefficient='Fixed posterior-ridge A=(I+W.T Tau W)^-1 W.T Tau; continuous prior comparison uses A/original view RMS. Full original prepared-X projection checked against saved own-view scores.',
              fidelity='Base background AND each ligand/receptor-component-excluded background must preserve original score: both discovery and original63 rho>=.95, RMSE<=.25 original discovery SD. At least200 retained background genes. No gate rescue or threshold relaxation.',
              exclusion='Literal ligand and all components; all flat-prior and complete-reference receptors and their components omitted from both continuous vectors before alignment.',
              reference='At least one accepted literal whole-ligand reference edge, all receptor components in native axis and each pooled-discovery positive>=1pct, with10pct sensitivity retained. Candidate-only marginal RNA, not same-cell/assembled protein or signalling proof.',
              alignment='Pearson continuous coefficient vs unsigned prior potential, Spearman sensitivity; not official NicheNet activity or P/FDR.',
              selection='Each fixed orientation: directional Pearson>0, fidelity and receptor marginal compatibility, top20 decreasing directional Pearson/lexical ligand. No rank21 replacement after later RNA/case qualification.',
              original_holdout_not_untouched_validation=True,clinical_values_or_external_effects_loaded=False,prior_sha256=PRIORS,
              reference_sha256=sha(refpath),native_positive_cache_sha256=sha(positive_path),model_sha256=sha(data['model_path']))
    (out/'plan_before_any_continuous_prior_alignment.json').write_text(json.dumps(plan,indent=2)+'\n')
    genes=data['genes'];units=data['units'];train=data['train'];cases=data['cases'].sampleID.to_numpy(str)
    hits=np.load(positive_path,mmap_mode='r');assert hits.shape==(len(units),31831)
    selections=[];totals=[]
    for view in VIEWS:
        select=units.sampleID.astype(str).isin(set(cases[train])) & units.sub_cell_type.isin(data['labels'][view])
        selections.append(units.loc[select,'unit_index'].to_numpy(int));totals.append(int(units.loc[select,'n_cells'].sum()))
        assert totals[-1]>0 and set(units.loc[select,'sampleID'].astype(str))==set(cases[train])
    detected=np.zeros((4,31831),np.int64)
    for start in range(0,31831,512):
        columns=np.arange(start,min(start+512,31831))
        for vi,rows in enumerate(selections):detected[vi,columns]=np.asarray(hits[np.ix_(rows,columns)]).sum(axis=0,dtype=np.int64)
    assert (detected<=np.asarray(totals)[:,None]).all()
    np.savez_compressed(out/'private_discovery_only_gene_detection.npz',positive=detected,cells=np.array(totals),views=np.array(VIEWS))
    target=read_rds(prior/'ligand_target_matrix_nsga2r_final.rds');target.index=target.index.astype(str);target.columns=target.columns.astype(str)
    assert target.index.is_unique and target.columns.is_unique and len(target.columns)==1226
    assert np.isfinite(target.to_numpy()).all() and (target.to_numpy()>=0).all()
    flat=read_rds(prior/'lr_network_human_21122021.rds');flat['from']=flat['from'].astype(str);flat['to']=flat['to'].astype(str)
    reference=pd.read_csv(refpath,sep='\t',keep_default_na=False);assert len(reference)==4624 and not reference.duplicated(['ligand','receptor']).any()
    gene_index={str(g):i for i,g in enumerate(genes)};ligands=sorted(target.columns.tolist())
    exclusions={};edges={ligand:reference[reference.ligand.eq(ligand)] for ligand in ligands}
    for ligand in ligands:
        receptors=set(flat.loc[flat['from'].eq(ligand),'to'])|set(edges[ligand].receptor)
        strings={ligand,*receptors};parts=set(strings)
        for s in strings:parts.update(str(s).split('_'))
        exclusions[ligand]=parts
    ranks=[];base_rows=[];edge_rows=[];shortlist=[];factors=['NSCLC_MFA_20261016_F%02d'%k for k in range(1,11)]
    for vi,view in enumerate(VIEWS):
        p=data['prepared'][vi];w=data['model']['weights_'+str(vi)];tau=data['model']['precision_'+str(vi)]
        a=np.linalg.solve(np.eye(10)+w.T@(tau[:,None]*w),w.T*tau)
        original=data['model']['own_view_scores_'+str(vi)]
        assert np.allclose(p['X']@a.T,original,atol=1e-10,rtol=1e-10)
        selected=p['selected_native_columns'];names=genes[selected]
        raw=pd.read_csv(root/'NSCLC_native_MFA_inputs_20261016_v1'/('private_'+view+'_pseudobulk.tsv.gz'),sep='\t',index_col=0)
        assert raw.columns.tolist()==cases.tolist() and raw.index.astype(str).tolist()==genes.tolist()
        assert ((raw.to_numpy()[selected][:,train]>=100).sum(axis=1)>=int(np.ceil(.25*train.sum()))).all()
        expressed=detected[vi]*100>=totals[vi];expressed10=detected[vi]*10>=totals[vi]
        bg=np.array([j for j,g in enumerate(names) if g in target.index],int)
        bg_names=names[bg];base=p['X'][:,bg]@a[:,bg].T
        base_values,base_gate=calibrate(base,original,train);base_gate &= len(bg)>=200
        for k,factor in enumerate(factors):
            base_rows.append(dict(view=view,factor=factor,selected_features=len(names),discovery_expressed_prior_features=len(bg),
                                  background_preserves_original_score=bool(base_gate[k]),
                                  **{field:float(value[k]) for field,value in base_values.items()}))
        compatible={}
        for ligand in ligands:
            valid=[]
            for edge in edges[ligand].itertuples():
                components=str(edge.receptor).split('_');measured=all(g in gene_index for g in components)
                marginal=measured and all(expressed[gene_index[g]] for g in components)
                marginal10=measured and all(expressed10[gene_index[g]] for g in components)
                accepted=edge.reference_status in ACCEPTED and marginal
                if accepted:valid.append(str(edge.receptor))
                edge_rows.append(dict(view=view,ligand=ligand,receptor=edge.receptor,reference_status=edge.reference_status,
                                      all_components_native=measured,all_component_discovery_marginals1pct=marginal,
                                      all_component_discovery_marginals10pct=marginal10,
                                      reference_and_RNA_marginal_compatible=accepted,same_cell_complex_verified=False))
            compatible[ligand]=sorted(set(valid))
        abs_all=np.abs(a).sum(axis=1);pos_all=np.maximum(a,0).sum(axis=1);neg_all=np.maximum(-a,0).sum(axis=1)
        local=[]
        prior_values=target.loc[bg_names,ligands].to_numpy(float)
        for j,ligand in enumerate(ligands):
            keep=np.array([str(g) not in exclusions[ligand] for g in bg_names]);remain=bg[keep];drop=bg[~keep]
            prediction=base-p['X'][:,drop]@a[:,drop].T if len(drop) else base.copy()
            cv,gate=calibrate(prediction,original,train);gate &= base_gate & (len(remain)>=200)
            coefficient=a[:,remain]/float(p['training_view_RMS']);prior_vector=prior_values[keep,j]
            pearson=np.full(10,np.nan);spearman=np.full(10,np.nan)
            if len(remain)>=200:
                y=prior_vector-prior_vector.mean();x=coefficient-coefficient.mean(axis=1,keepdims=True)
                den=np.linalg.norm(x,axis=1)*np.linalg.norm(y);good=gate&(den>1e-15)
                pearson[good]=(x[good]@y)/den[good]
                yr=rankdata(prior_vector);yr-=yr.mean();xr=rankdata(coefficient,axis=1);xr-=xr.mean(axis=1,keepdims=True)
                dr=np.linalg.norm(xr,axis=1)*np.linalg.norm(yr);valid=good&(dr>1e-15)
                spearman[valid]=(xr[valid]@yr)/dr[valid]
            for k,factor in enumerate(factors):
                status='CONTINUOUS_PRIOR_ALIGNMENT_MEASURABLE'
                if not base_gate[k]:status='BASE_PRIOR_BACKGROUND_DOES_NOT_PRESERVE_PROGRAMME'
                elif not gate[k]:status='LR_EXCLUDED_BACKGROUND_DOES_NOT_PRESERVE_PROGRAMME'
                elif not np.isfinite(pearson[k]):status='CONSTANT_PRIOR_OR_PROGRAMME_VECTOR'
                local.append(dict(view=view,factor=factor,ligand=ligand,background_genes=len(bg),LR_excluded_background_genes=len(remain),
                                  absolute_coefficient_mass_retained=float(np.abs(a[k,remain]).sum()/abs_all[k]) if abs_all[k]>0 else np.nan,
                                  positive_coefficient_mass_retained=float(np.maximum(a[k,remain],0).sum()/pos_all[k]) if pos_all[k]>0 else np.nan,
                                  negative_coefficient_mass_retained=float(np.maximum(-a[k,remain],0).sum()/neg_all[k]) if neg_all[k]>0 else np.nan,
                                  base_background_gate=bool(base_gate[k]),LR_excluded_background_gate=bool(gate[k]),
                                  prior_Pearson=float(pearson[k]),prior_Spearman=float(spearman[k]),
                                  receiver_reference_RNA_compatible=bool(compatible[ligand]),compatible_receptors='|'.join(compatible[ligand]),
                                  status=status,**{field:float(value[k]) for field,value in cv.items()}))
        table=pd.DataFrame(local);assert len(table)==12260
        for factor in factors:
            for direction in [1,-1]:
                group=table[table.factor.eq(factor)].copy();group['directional_prior_alignment']=direction*group.prior_Pearson
                choose=group[group.status.eq('CONTINUOUS_PRIOR_ALIGNMENT_MEASURABLE')&group.receiver_reference_RNA_compatible&group.directional_prior_alignment.gt(0)]
                choose=choose.sort_values(['directional_prior_alignment','ligand'],ascending=[False,True]).head(20).to_dict('records')
                for pos in range(20):
                    row=dict(view=view,factor=factor,direction=direction,prospective_position=pos+1)
                    if pos<len(choose):row.update(choose[pos]);row['slot_status']='FROZEN_PRIOR_CANDIDATE_REQUIRES_SENDER_RECEIVER_CASE_QUALIFICATION'
                    else:row.update(ligand='',slot_status='NO_ELIGIBLE_PRIOR_CANDIDATE_UNDER_FIXED_RULE')
                    shortlist.append(row)
        ranks.extend(local)
        (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),completed_views=vi+1,total_views=4))+'\n')
    pd.DataFrame(base_rows).to_csv(out/'private_all40_base_background_fidelity.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(ranks).to_csv(out/'private_all49040_continuous_prior_alignments.tsv.gz',sep='\t',index=False,compression='gzip',na_rep='NA')
    pd.DataFrame(edge_rows).to_csv(out/'private_all_discovery_receptor_reference_qualifications.tsv.gz',sep='\t',index=False,compression='gzip')
    slots=pd.DataFrame(shortlist);assert len(slots)==1600
    slots.to_csv(out/'private_all1600_prospective_ligand_slots.tsv.gz',sep='\t',index=False,compression='gzip',na_rep='NA')
    terminal=dict(status='COMPLETE_ALL40_PROGRAMME_VECTORS49040_PRIOR_ALIGNMENTS_AND1600_PROSPECTIVE_SLOTS_AUDIT_PENDING',
                  finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
                  files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
                  external_effects_or_clinical_values_loaded=False,official_NicheNet_activity_or_direct_communication_validation=False)
    assert sum(v['bytes'] for v in terminal['files'].values())<budget['output_cap_bytes']
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

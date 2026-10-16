"""All frozen-vector prior inputs: sparse counts, QR and independent metrics."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.sparse import csr_matrix
from scipy.stats import pearsonr,spearmanr

from NSCLC_native_MFA_interpretation_inputs_v1 import load,sha,VIEWS
from prepare_NSCLC_fixed_programme_continuous_priors_v1 import read_rds,PRIORS,ACCEPTED
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def fidelity(predicted,original,train):
    values={};gate=np.ones(10,bool);sd=np.std(original[train],axis=0)
    for tag,mask in [('discovery',train),('original_holdout',~train)]:
        correlations=np.array([pearsonr(predicted[mask,k],original[mask,k]).statistic if np.ptp(predicted[mask,k])>0 and np.ptp(original[mask,k])>0 else np.nan for k in range(10)])
        errors=np.sqrt(np.mean(np.square(predicted[mask]-original[mask]),axis=0))/sd
        values[tag+'_rho']=correlations;values[tag+'_RMSE_over_discovery_SD']=errors
        gate &= np.isfinite(correlations)&np.isfinite(errors)&(correlations>=.95)&(errors<=.25)
    return values,gate


def main():
    data=load();root=data['root'];resource(root,1048576)
    src=root/'NSCLC_fixed_programme_continuous_prior_inputs_20261016_v1'
    terminal=json.loads((src/'private_terminal.json').read_text())
    for name,item in terminal['files'].items():assert sha(src/name)==item['sha256']
    plan=json.loads((src/'plan_before_any_continuous_prior_alignment.json').read_text())
    assert terminal['external_effects_or_clinical_values_loaded'] is False
    positive_path=data['native']/'positive_cell_counts.npy';assert sha(positive_path)==plan['native_positive_cache_sha256']
    native=data['units'];cases=data['cases'].sampleID.astype(str).to_numpy();train=data['train']
    rows=[];columns=[]
    for vi,view in enumerate(VIEWS):
        take=native.sampleID.astype(str).isin(set(cases[train]))&native.sub_cell_type.isin(data['labels'][view])
        ids=native.loc[take,'unit_index'].to_numpy(int);rows.extend([vi]*len(ids));columns.extend(ids.tolist())
    mixing=csr_matrix((np.ones(len(rows),np.int64),(np.array(rows),np.array(columns))),shape=(4,len(native)))
    totals=np.asarray(mixing@native.n_cells.to_numpy(np.int64)).ravel()
    hits=np.load(positive_path,mmap_mode='r');pooled=np.zeros((4,31831),np.int64)
    for start in range(0,31831,512):pooled[:,start:start+512]=np.asarray(mixing@hits[:,start:start+512])
    with np.load(src/'private_discovery_only_gene_detection.npz',allow_pickle=False) as z:
        assert np.array_equal(z['positive'],pooled) and np.array_equal(z['cells'],totals)
    prior=root/'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1/verified_prior_inputs_v2'
    for name,digest in PRIORS.items():assert sha(prior/name)==digest
    target=read_rds(prior/'ligand_target_matrix_nsga2r_final.rds');target.index=target.index.astype(str);target.columns=target.columns.astype(str)
    flat=read_rds(prior/'lr_network_human_21122021.rds');flat['from']=flat['from'].astype(str);flat['to']=flat['to'].astype(str)
    refpath=root/'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv';assert sha(refpath)==plan['reference_sha256']
    ref=pd.read_csv(refpath,sep='\t',keep_default_na=False)
    base=pd.read_csv(src/'private_all40_base_background_fidelity.tsv',sep='\t').set_index(['view','factor'])
    ranking=pd.read_csv(src/'private_all49040_continuous_prior_alignments.tsv.gz',sep='\t',keep_default_na=False)
    numeric=['absolute_coefficient_mass_retained','positive_coefficient_mass_retained','negative_coefficient_mass_retained','prior_Pearson','prior_Spearman','discovery_rho','discovery_RMSE_over_discovery_SD','original_holdout_rho','original_holdout_RMSE_over_discovery_SD']
    for field in numeric:ranking[field]=pd.to_numeric(ranking[field],errors='coerce')
    assert len(ranking)==49040 and not ranking.duplicated(['view','factor','ligand']).any()
    reference=pd.read_csv(src/'private_all_discovery_receptor_reference_qualifications.tsv.gz',sep='\t',keep_default_na=False)
    assert not reference.duplicated(['view','ligand','receptor']).any()
    slots=pd.read_csv(src/'private_all1600_prospective_ligand_slots.tsv.gz',sep='\t',keep_default_na=False)
    assert len(slots)==1600 and not slots.duplicated(['view','factor','direction','prospective_position']).any()
    index=ranking.set_index(['view','factor','ligand']);ref_index=reference.set_index(['view','ligand','receptor'])
    gene={str(g):i for i,g in enumerate(data['genes'])};ligands=sorted(target.columns.astype(str));factors=['NSCLC_MFA_20261016_F%02d'%k for k in range(1,11)]
    expected_edges=set();maximum=0.;checked=0
    checkpoint=root/'HNSCC_receiver_followup_job_control_20261016_v2'/Path(__file__).stem/'private_audit_checkpoint.json'
    for vi,view in enumerate(VIEWS):
        p=data['prepared'][vi];w=data['model']['weights_'+str(vi)];sw=np.sqrt(data['model']['precision_'+str(vi)])
        q,r=np.linalg.qr(np.vstack([sw[:,None]*w,np.eye(10)]),mode='reduced')
        aq=solve_triangular(r,q[:len(w)].T*sw)
        # Canonical float coefficients define exact ties; QR validates them.
        canonical=np.linalg.solve(np.eye(10)+w.T@((sw*sw)[:,None]*w),w.T*(sw*sw))
        assert np.allclose(aq,canonical,atol=1e-10,rtol=1e-10)
        original=data['model']['own_view_scores_'+str(vi)];assert np.allclose(p['X']@aq.T,original,atol=1e-10,rtol=1e-10)
        selected=p['selected_native_columns'];names=data['genes'][selected]
        raw=pd.read_csv(root/'NSCLC_native_MFA_inputs_20261016_v1'/('private_'+view+'_pseudobulk.tsv.gz'),sep='\t',index_col=0)
        assert raw.columns.tolist()==cases.tolist()
        assert ((raw.to_numpy()[selected][:,train]>=100).sum(axis=1)>=int(np.ceil(.25*train.sum()))).all()
        bg=np.flatnonzero(np.array([str(name) in target.index for name in names]));bg_names=names[bg]
        base_prediction=p['X'][:,bg]@aq[:,bg].T
        bs,bg_gate=fidelity(base_prediction,original,train);bg_gate &= len(bg)>=200
        for k,factor in enumerate(factors):
            observed=base.loc[(view,factor)]
            assert observed.background_preserves_original_score==bg_gate[k] and observed.discovery_expressed_prior_features==len(bg)
            for field,value in bs.items():assert np.isclose(observed[field],value[k],atol=1e-9,rtol=1e-9,equal_nan=True)
        marginal=pooled[vi]*100>=totals[vi];marginal10=pooled[vi]*10>=totals[vi]
        for ligand in ligands:
            relevant=ref[ref.ligand.eq(ligand)];compatible=[]
            excluded={ligand,*flat.loc[flat['from'].eq(ligand),'to'].tolist(),*relevant.receptor.tolist()}
            removed=set(excluded)
            for name in excluded:removed.update(str(name).split('_'))
            for edge in relevant.itertuples():
                key=(view,ligand,edge.receptor);expected_edges.add(key);observed=ref_index.loc[key]
                parts=str(edge.receptor).split('_');measured=all(part in gene for part in parts)
                m1=measured and all(marginal[gene[part]] for part in parts);m10=measured and all(marginal10[gene[part]] for part in parts)
                accepted=edge.reference_status in ACCEPTED and m1
                assert observed.all_components_native==measured and observed.all_component_discovery_marginals1pct==m1
                assert observed.all_component_discovery_marginals10pct==m10 and observed.reference_and_RNA_marginal_compatible==accepted
                assert observed.same_cell_complex_verified is False or observed.same_cell_complex_verified==False
                if accepted:compatible.append(str(edge.receptor))
            keep=np.array([str(g) not in removed for g in bg_names]);remain=bg[keep]
            prediction=p['X'][:,remain]@aq[:,remain].T
            cs,gate=fidelity(prediction,original,train);gate &= bg_gate & (len(remain)>=200)
            values=target.loc[names[remain],ligand].to_numpy(float)
            for k,factor in enumerate(factors):
                observed=index.loc[(view,factor,ligand)];status='CONTINUOUS_PRIOR_ALIGNMENT_MEASURABLE'
                x=aq[k,remain]/float(p['training_view_RMS']);canon=canonical[k,remain]/float(p['training_view_RMS'])
                den=np.linalg.norm(x-x.mean())*np.linalg.norm(values-values.mean()) if len(remain) else 0.
                corr=float(pearsonr(x,values).statistic) if gate[k] and den>1e-15 else np.nan
                sensitivity=float(spearmanr(canon,values).statistic) if np.isfinite(corr) else np.nan
                if not bg_gate[k]:status='BASE_PRIOR_BACKGROUND_DOES_NOT_PRESERVE_PROGRAMME'
                elif not gate[k]:status='LR_EXCLUDED_BACKGROUND_DOES_NOT_PRESERVE_PROGRAMME'
                elif not np.isfinite(corr):status='CONSTANT_PRIOR_OR_PROGRAMME_VECTOR'
                assert observed.status==status and observed.base_background_gate==bg_gate[k] and observed.LR_excluded_background_gate==gate[k]
                assert observed.LR_excluded_background_genes==len(remain) and observed.background_genes==len(bg)
                assert observed.receiver_reference_RNA_compatible==bool(compatible)
                assert observed.compatible_receptors=='|'.join(sorted(set(compatible)))
                assert np.isclose(observed.prior_Pearson,corr,equal_nan=True,atol=1e-9,rtol=1e-9)
                assert np.isclose(observed.prior_Spearman,sensitivity,equal_nan=True,atol=1e-9,rtol=1e-9)
                for field,value in cs.items():assert np.isclose(observed[field],value[k],equal_nan=True,atol=1e-9,rtol=1e-9)
                for field,all_weights,part in [('absolute_coefficient_mass_retained',np.abs(aq[k]),np.abs(aq[k,remain])),('positive_coefficient_mass_retained',np.maximum(aq[k],0),np.maximum(aq[k,remain],0)),('negative_coefficient_mass_retained',np.maximum(-aq[k],0),np.maximum(-aq[k,remain],0))]:
                    value=float(part.sum()/all_weights.sum()) if all_weights.sum()>0 else np.nan
                    assert np.isclose(observed[field],value,equal_nan=True,atol=1e-9,rtol=1e-9)
                if np.isfinite(corr):maximum=max(maximum,abs(corr-observed.prior_Pearson))
                checked+=1
        checkpoint.write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),completed_views=vi+1,total_views=4))+'\n')
    assert set(ref_index.index)==expected_edges
    # Rank the checked canonical floating measurements, preserving exact tie policy.
    for (view,factor,direction),group in slots.groupby(['view','factor','direction']):
        candidates=ranking[ranking.view.eq(view)&ranking.factor.eq(factor)].copy()
        candidates['directional_prior_alignment']=int(direction)*candidates.prior_Pearson
        candidates=candidates[candidates.status.eq('CONTINUOUS_PRIOR_ALIGNMENT_MEASURABLE')&candidates.receiver_reference_RNA_compatible&candidates.directional_prior_alignment.gt(0)]
        candidates=candidates.sort_values(['directional_prior_alignment','ligand'],ascending=[False,True]).head(20)
        observed=group.sort_values('prospective_position');assert observed.prospective_position.tolist()==list(range(1,21))
        assert observed.ligand.tolist()==candidates.ligand.tolist()+['']*(20-len(candidates))
        expected=['FROZEN_PRIOR_CANDIDATE_REQUIRES_SENDER_RECEIVER_CASE_QUALIFICATION']*len(candidates)+['NO_ELIGIBLE_PRIOR_CANDIDATE_UNDER_FIXED_RULE']*(20-len(candidates))
        assert observed.slot_status.tolist()==expected
    proof=src/'private_independent_continuous_prior_audit.json';assert not proof.exists()
    report=dict(status='PASS_ALL40_QR_VECTORS_DISCOVERY_DETECTION49040_ALIGNMENTS_FIDELITY_REFERENCE_EDGES_AND1600_SLOTS',
                finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(src/'private_terminal.json'),
                all_alignment_records_checked=checked,maximum_continuous_Pearson_error=maximum,
                canonical_ties='QR verifies original continuous coefficients; discrete coefficient-rank and shortlist ties use the unchanged canonical float coefficients/checked stored Pearson values.',
                clinical_or_external_effects_loaded=False,direct_communication_or_clinical_validation=False)
    proof.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],audit_sha256=sha(proof))))


if __name__=='__main__':main()

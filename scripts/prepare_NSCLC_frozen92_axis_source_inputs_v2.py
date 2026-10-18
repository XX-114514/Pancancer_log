"""Freeze all92 hypotheses and prepare native RNA, controls, blocks, no effects."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from preflight_NSCLC_native_MFA_clinical_coverage_v1 import COMMON
from NSCLC_frozen92_axis_contract_v1 import CONTRACT
from NSCLC_frozen_axis_block_math_v2 import fixture
from NSCLC_frozen92_axis_io_v2 import prepare,controls


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,40*1024**2);out=root/'NSCLC_frozen92_axis_source_inputs_20261019_v2';assert not out.exists();out.mkdir(mode=0o700)
    original=root/'NSCLC_communication_tensor_family_20261019_v1';blindroot=root/'NSCLC_communication_tensor_inputs_20261019_v1';native=root/'NSCLC_full_native_gene_units_20261016_v1';cat=root/'NSCLC_tensor_joint_candidate_expression_qualification_20261019_v1/private_CPG04_two_stable_ranks_all_common_joint_candidate_names_and_coverage_v1.tsv.gz'
    assert sha(cat)=='a63e70e488239ad337c1699c2afdcfb6cbef301646d7c4b42e39a5dd3ee759b6';assert sha(root/'NSCLC_complete_native_LR_statistics_audit_20261019_v1/private_terminal.json')=='64e1604175e43b3816ab3dc626ec60300b03a35a5d6b09a2d3e1244b2d4d8c2e';assert sha(root/'NSCLC_tensor_joint_candidate_expression_qualification_20261019_v1/private_terminal.json')=='82344b0f9d8204796b3bea185038889eb68c9b45322dd0e170a388b85a1ccd49'
    save=lambda n,d:(out/n).write_text(json.dumps(d,indent=2)+'\n')
    modules=['NSCLC_frozen92_axis_contract_v1.py','NSCLC_frozen_axis_block_math_v2.py','NSCLC_frozen92_axis_io_v2.py'];plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),module_sha256={n:sha(Path(__file__).with_name(n)) for n in modules},total_resource_guard=guard,CONTRACT=CONTRACT,catalogue_sha256=sha(cat),new_source_effects_read=False)
    save('plan_before_all92_source_axis_readouts_and_effects.json',plan);save('private_synthetic_block_math_gate.json',fixture())
    original_t=json.loads((original/'private_terminal.json').read_text());axespath=original/'private_tensor_axes.npz';vocabpath=original/'private_selected_LR_vocabulary.json'
    for p in [axespath,vocabpath]:assert sha(p)==original_t['files'][p.name]['sha256']
    axes=dict(np.load(axespath));cases=axes['discovery_ids'].tolist();states=axes['states'].tolist();N=len(cases);S=len(states);assert (N,S)==(129,38);vocab=json.loads(vocabpath.read_text());L=len(vocab);assert L==876
    candidates=pd.read_csv(cat,sep='\t',keep_default_na=False);assert len(candidates)==92 and not candidates.duplicated(['canonical_resource_LR_id','sender_original_author_type','receiver_original_author_type']).any();candidates.insert(0,'axis_id',['AX%03d'%(j+1) for j in range(92)]);vi={d['canonical_resource_LR_id']:i for i,d in enumerate(vocab)};si={s:i for i,s in enumerate(states)};li=np.array([vi[v] for v in candidates.canonical_resource_LR_id]);sender=np.array([si[v] for v in candidates.sender_original_author_type]);receiver=np.array([si[v] for v in candidates.receiver_original_author_type]);candidates.to_csv(out/'private_all92_frozen_axis_registry.tsv',sep='\t',index=False)
    blind=pd.read_csv(blindroot/'private_all243_blind_record_partitions.tsv',sep='\t',keep_default_na=False).set_index('sampleID');assert blind.loc[cases].partition.eq('DISCOVERY').all() and not set(cases)&set(axes['heldout_ids']);assert blind.loc[cases].profile_block.nunique()==128
    units=pd.read_csv(native/'private_all_native_units_and_quality.tsv',sep='\t');assert units.unit_index.tolist()==list(range(len(units)));lookup={(r.sampleID,r.sub_cell_type):r.unit_index for r in units.itertuples()};cols=axes['native_gene_columns'];rev={int(c):j for j,c in enumerate(cols)};M=np.zeros((N,S,len(cols)));det=np.zeros_like(M,dtype=bool);cells=np.zeros((N,S),int);quality=np.full((N,S,2),np.nan);unit_id=np.full((N,S),-1,int);major={s:units.loc[units.sub_cell_type.eq(s),'major_cell_type'].unique().tolist() for s in states};assert all(len(v)==1 for v in major.values())
    mean=np.load(native/'mean_logCP10K.npy',mmap_mode='r');positive=np.load(native/'positive_cell_counts.npy',mmap_mode='r')
    for i,case in enumerate(cases):
        for j,state in enumerate(states):
            unit=lookup.get((case,state))
            if unit is not None:
                unit_id[i,j]=unit;n=int(units.iloc[unit].n_cells);cells[i,j]=n;quality[i,j]=units.iloc[unit][['mean_logRNA','mean_loggenes']].to_numpy(float);M[i,j]=mean[unit,cols];det[i,j]=positive[unit,cols]*10>=n
    del mean,positive;available=cells>=20;np.testing.assert_array_equal(available,axes['availability'][:129]);lig=np.zeros((N,L,S),np.float32);rec=np.zeros_like(lig);left=np.zeros((N,L,S),bool);right=np.zeros_like(left)
    for j,d in enumerate(vocab):
        lc=[rev[c] for c in d['native_ligand_columns']];rc=[rev[c] for c in d['native_receptor_columns']];lig[:,j]=M[:,:,lc].min(2);rec[:,j]=M[:,:,rc].min(2);left[:,j]=det[:,:,lc].all(2);right[:,j]=det[:,:,rc].all(2)
    del M,det;raw=np.log1p(lig[:,li,sender].astype(float)*rec[:,li,receiver].astype(float)).astype(np.float32).astype(float);measurable=available[:,sender]&available[:,receiver];supported=left[:,li,sender]&right[:,li,receiver]&measurable;raw[~measurable]=np.nan;standardized=np.full_like(raw,np.nan);scales=[]
    for j in range(92):
        values=raw[measurable[:,j],j];mu=float(values.mean()) if len(values) else np.nan;sd=float(values.std()) if len(values) else np.nan;scales.append([mu,sd])
        if np.isfinite(sd) and sd>1e-12:standardized[measurable[:,j],j]=(values-mu)/sd
    burden=np.zeros(N);rms=float(axes['global_RMS'])
    for i in range(N):
        score=np.log1p(lig[i].astype(float)[:,:,None]*rec[i].astype(float)[:,None,:]).astype(np.float32).astype(float)/rms;a=available[i];burden[i]=np.log1p(score[:,a][:,:,a].mean())
    del lig,rec,left,right
    all_totals=units.groupby('sampleID').n_cells.sum();totals=all_totals.reindex(cases).to_numpy()
    global_quality=np.column_stack([np.log1p(totals)]+[(units.assign(weighted=units.n_cells*units[f]).groupby('sampleID').weighted.sum()/all_totals).reindex(cases).to_numpy() for f in ['mean_logRNA','mean_loggenes']])
    denom=units.groupby(['sampleID','major_cell_type']).n_cells.sum();fractions=np.zeros_like(cells,float)
    for i,case in enumerate(cases):
        for j,state in enumerate(states):
            n=denom.get((case,major[state][0]),0);fractions[i,j]=cells[i,j]/n if n else np.nan
    levels=sorted(blind.loc[blind.partition.isin(['DISCOVERY','HELDOUT']),'agent_identity'].unique());agent=np.array([[float(v==level) for level in levels[1:]] for v in blind.loc[cases,'agent_identity']]);agent_names=['agent:'+v for v in levels[1:]]
    # First clinical outcome read: hypotheses, scales and covariate definitions fixed.
    save('private_before_first_new_outcome_read.json',dict(observed_UTC=datetime.now(timezone.utc).isoformat(),plan_sha256=sha(out/'plan_before_all92_source_axis_readouts_and_effects.json'),axis_scores_and_blind_scales_fixed=True,new_clinical_effects_read=False))
    clinicalpath=root/'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1/private_all243_current_RNA_clinical_scope.tsv';clinical=pd.read_csv(clinicalpath,sep='\t',keep_default_na=False).set_index('sampleID').loc[cases]
    base=np.c_[clinical[COMMON].apply(pd.to_numeric,errors='coerce').to_numpy(float),agent];response=pd.to_numeric(clinical.response,errors='coerce').to_numpy(float)
    arrays=dict(author_record_ids=np.asarray(cases),profile_blocks=blind.loc[cases,'profile_block'].to_numpy(str),cancer_type=blind.loc[cases,'cancer_type'].to_numpy(str),response=response,clinical_controls=base,clinical_control_names=np.asarray(COMMON+agent_names),global_quality=global_quality,unit_cells=cells,unit_quality=quality,within_major_fraction=fractions,author_states=np.asarray(states),original_unit_indices=unit_id,LR_index=li,sender_index=sender,receiver_index=receiver,axis_ids=candidates.axis_id.to_numpy(str),raw_axis_score=raw,standardized_axis_score=standardized,endpoint_measurable=measurable,all_complex_components_detection_supported=supported,blind_score_mean_SD=np.asarray(scales),log1p_mean_observed_normalized_LRscore=burden,captured_author_types=available.sum(1),all_agent_levels=np.asarray(levels))
    np.savez_compressed(out/'private_all129_all92_native_axis_and_design_inputs.npz',**arrays);coverage=[];excluded=[]
    for j,axis in enumerate(candidates.axis_id):
        model,g,status,valid=prepare(arrays,j);x=response[valid];coverage.append(dict(axis_id=axis,source='GSE243013',partition='ORIGINAL_DISCOVERY_SEPARATE_LABEL_SET',author_records=int(valid.sum()),nonMPR=int((x==0).sum()),pCR_or_MPR=int((x==1).sum()),profile_blocks=len(set(arrays['profile_blocks'][valid])),reviewed_independent_patients=0,endpoint_measurable=int(measurable[:,j].sum()),strict_detection_supported=int(supported[:,j].sum()),measured_true_zero_scores=int(np.sum(measurable[:,j]&(raw[:,j]==0))),status=status,nuisance_rank=model['nuisance_rank'] if model else np.nan,residual_df=model['df'] if model else np.nan,max_block_leverage=g['max_block_leverage'] if g else np.nan))
        C,names=controls(arrays,j)
        for i,case in enumerate(cases):
            reason=[]
            for label,t in [('SENDER',sender[j]),('RECEIVER',receiver[j])]:
                if unit_id[i,t]<0:reason.append(label+'_AUTHOR_UNIT_ABSENT')
                elif cells[i,t]<20:reason.append(label+'_FEWER_THAN20_CELLS')
            if not np.isfinite(response[i]):reason.append('MISSING_ORIGINAL_RESPONSE')
            if not np.isfinite(C[i]).all():reason.append('MISSING_FIXED_CONTROLS')
            if not np.isfinite(standardized[i,j]) and measurable[i,j]:reason.append('CONSTANT_OR_INVALID_BLIND_SCORE_SCALE')
            excluded.append(dict(axis_id=axis,author_record_id=case,in_main_design=bool(valid[i]),reasons='|'.join(reason) if reason else 'INCLUDED'))
    pd.DataFrame(coverage).to_csv(out/'private_all92_source_design_coverage.tsv',sep='\t',index=False);pd.DataFrame(excluded).to_csv(out/'private_all11868_record_design_exclusions.tsv.gz',sep='\t',index=False,compression='gzip')
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<2*1024**2
    t=dict(status='COMPLETE_ALL92_FROZEN_SOURCE_INPUTS_AND_DESIGN_COVERAGE_NO_TRUE_EFFECTS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,candidate_slots=92,author_records=129,source_profile_blocks=128,clinical_source_sha256=sha(clinicalpath),native_unit_source_sha256=sha(native/'private_all_native_units_and_quality.tsv'),true_clinical_effects_read=False,scientific_outputs_public=False)
    save('private_terminal.json',t);print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

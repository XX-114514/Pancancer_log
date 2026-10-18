"""Every frozen axis native score, control, exclusion, block and QR geometry."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from preflight_NSCLC_native_MFA_clinical_coverage_v1 import COMMON
from NSCLC_frozen92_axis_io_v2 import prepare as production_prepare,controls
from NSCLC_frozen_axis_block_audit_math_v1 import prepare_axis,synthetic_stream_gate


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,2*1024**2);out=root/'NSCLC_frozen92_axis_source_input_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700);src=root/'NSCLC_frozen92_axis_source_inputs_20261019_v2';t=json.loads((src/'private_terminal.json').read_text());assert t['source_sha256']=='698e06c0b2f8e224d1c7ede17a41c0054ee811ca033675f429f31af6ddf6162d' and not t['true_clinical_effects_read']
    for name,d in t['files'].items():assert sha(src/name)==d['sha256']
    save=lambda name,data:(out/name).write_text(json.dumps(data,indent=2)+'\n')
    save('private_independent_synthetic_stream_gate.json',synthetic_stream_gate());a=dict(np.load(src/'private_all129_all92_native_axis_and_design_inputs.npz'));cases=a['author_record_ids'].tolist();states=a['author_states'].tolist();original=root/'NSCLC_communication_tensor_family_20261019_v1';axes=dict(np.load(original/'private_tensor_axes.npz'));assert cases==axes['discovery_ids'].tolist() and states==axes['states'].tolist() and not set(cases)&set(axes['heldout_ids']);vocab=json.loads((original/'private_selected_LR_vocabulary.json').read_text());registry=pd.read_csv(src/'private_all92_frozen_axis_registry.tsv',sep='\t',keep_default_na=False);cat=root/'NSCLC_tensor_joint_candidate_expression_qualification_20261019_v1/private_CPG04_two_stable_ranks_all_common_joint_candidate_names_and_coverage_v1.tsv.gz';assert sha(cat)=='a63e70e488239ad337c1699c2afdcfb6cbef301646d7c4b42e39a5dd3ee759b6';old=pd.read_csv(cat,sep='\t',keep_default_na=False);pd.testing.assert_frame_equal(registry.drop(columns='axis_id'),old);assert registry.axis_id.tolist()==a['axis_ids'].tolist()==['AX%03d'%(j+1) for j in range(92)]
    blind=pd.read_csv(root/'NSCLC_communication_tensor_inputs_20261019_v1/private_all243_blind_record_partitions.tsv',sep='\t',keep_default_na=False).set_index('sampleID');np.testing.assert_array_equal(a['profile_blocks'],blind.loc[cases,'profile_block'].to_numpy(str));np.testing.assert_array_equal(a['cancer_type'],blind.loc[cases,'cancer_type'].to_numpy(str));assert len(set(a['profile_blocks']))==128
    raw=root/'NSCLC_full_native_gene_units_20261016_v1';u=pd.read_csv(raw/'private_all_native_units_and_quality.tsv',sep='\t');lookup=u.set_index(['sampleID','sub_cell_type']);mean=np.load(raw/'mean_logCP10K.npy',mmap_mode='r');positive=np.load(raw/'positive_cell_counts.npy',mmap_mode='r');N=len(cases);S=len(states);M=np.zeros((N,S,len(axes['native_gene_columns'])));ids=np.full((N,S),-1,int);cells=np.zeros((N,S),int);quality=np.full((N,S,2),np.nan);fractions=np.full((N,S),np.nan);global_quality=[]
    for i,case in enumerate(cases):
        source=u[u.sampleID.eq(case)];total=int(source.n_cells.sum());global_quality.append([np.log1p(total),np.average(source.mean_logRNA,weights=source.n_cells),np.average(source.mean_loggenes,weights=source.n_cells)])
        for j,state in enumerate(states):
            major=u.loc[u.sub_cell_type.eq(state),'major_cell_type'].unique();assert len(major)==1;denominator=int(source.loc[source.major_cell_type.eq(major[0]),'n_cells'].sum());fractions[i,j]=0. if denominator else np.nan
            if (case,state) in lookup.index:
                row=lookup.loc[(case,state)];unit=int(row.unit_index);ids[i,j]=unit;cells[i,j]=int(row.n_cells);quality[i,j]=[row.mean_logRNA,row.mean_loggenes];fractions[i,j]=cells[i,j]/denominator;M[i,j]=mean[unit,axes['native_gene_columns']]
    for key,value in [('original_unit_indices',ids),('unit_cells',cells),('unit_quality',quality),('within_major_fraction',fractions),('global_quality',np.asarray(global_quality))]:np.testing.assert_allclose(a[key],value,atol=1e-10,rtol=1e-10,equal_nan=True)
    available=cells>=20;np.testing.assert_array_equal(available,axes['availability'][:129]);np.testing.assert_array_equal(a['captured_author_types'],available.sum(1));vindex={d['canonical_resource_LR_id']:i for i,d in enumerate(vocab)};sindex={s:i for i,s in enumerate(states)};score=np.full((N,92),np.nan);supported=np.zeros((N,92),bool);endpoint=np.zeros((N,92),bool)
    for j,row in enumerate(registry.itertuples()):
        li=vindex[row.canonical_resource_LR_id];s=sindex[row.sender_original_author_type];r=sindex[row.receiver_original_author_type];d=vocab[li];assert (a['LR_index'][j],a['sender_index'][j],a['receiver_index'][j])==(li,s,r);assert row.ligand_complex=='_'.join(d['ligand_components']) and row.receptor_complex=='_'.join(d['receptor_components'])
        for i in range(N):
            if cells[i,s]<20 or cells[i,r]<20:continue
            endpoint[i,j]=True;ls=np.float32(min(float(mean[ids[i,s],c]) for c in d['native_ligand_columns']));rs=np.float32(min(float(mean[ids[i,r],c]) for c in d['native_receptor_columns']));score[i,j]=float(np.float32(np.log1p(float(ls)*float(rs))));supported[i,j]=all(int(positive[ids[i,s],c])*10>=cells[i,s] for c in d['native_ligand_columns']) and all(int(positive[ids[i,r],c])*10>=cells[i,r] for c in d['native_receptor_columns'])
    np.testing.assert_allclose(score,a['raw_axis_score'],atol=0,rtol=0,equal_nan=True);np.testing.assert_array_equal(endpoint,a['endpoint_measurable']);np.testing.assert_array_equal(supported,a['all_complex_components_detection_supported']);standard=np.full_like(score,np.nan);scales=[]
    for j in range(92):
        v=score[endpoint[:,j],j];mu=v.mean() if len(v) else np.nan;sd=v.std() if len(v) else np.nan;scales.append([mu,sd])
        if np.isfinite(sd) and sd>1e-12:standard[endpoint[:,j],j]=(v-mu)/sd
    np.testing.assert_allclose(a['standardized_axis_score'],standard,atol=1e-12,rtol=1e-12,equal_nan=True);np.testing.assert_allclose(a['blind_score_mean_SD'],scales,atol=1e-12,rtol=1e-12,equal_nan=True)
    reverse={int(c):j for j,c in enumerate(axes['native_gene_columns'])};rms=float(axes['global_RMS']);burden=[]
    for i in range(N):
        mask=available[i];total=0.
        for d in vocab:
            l=M[i][:,[reverse[c] for c in d['native_ligand_columns']]].min(1).astype(np.float32)[mask];r=M[i][:,[reverse[c] for c in d['native_receptor_columns']]].min(1).astype(np.float32)[mask];total+=float(np.log1p(np.outer(l.astype(float),r.astype(float))).astype(np.float32).astype(float).sum())/rms
        burden.append(np.log1p(total/(876*int(mask.sum())**2)))
    np.testing.assert_allclose(a['log1p_mean_observed_normalized_LRscore'],burden,atol=1e-12,rtol=1e-12);del mean,positive,M
    clinical=pd.read_csv(root/'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1/private_all243_current_RNA_clinical_scope.tsv',sep='\t',keep_default_na=False).set_index('sampleID').loc[cases];np.testing.assert_allclose(a['response'],pd.to_numeric(clinical.response,errors='coerce'),equal_nan=True);levels=sorted(blind.loc[blind.partition.isin(['DISCOVERY','HELDOUT']),'agent_identity'].unique());agent=np.asarray([[float(value==level) for level in levels[1:]] for value in blind.loc[cases,'agent_identity']]);expected=np.column_stack([clinical[COMMON].apply(pd.to_numeric,errors='coerce').to_numpy(float),agent]);np.testing.assert_allclose(a['clinical_controls'],expected,equal_nan=True);assert a['clinical_control_names'].tolist()==COMMON+['agent:'+v for v in levels[1:]] and a['all_agent_levels'].tolist()==levels
    coverage=pd.read_csv(src/'private_all92_source_design_coverage.tsv',sep='\t');ex=pd.read_csv(src/'private_all11868_record_design_exclusions.tsv.gz',sep='\t',keep_default_na=False);assert len(coverage)==92 and len(ex)==11868;maximum=0.
    for j,axis in enumerate(a['axis_ids']):
        p,g,status,valid=production_prepare(a,j);m,h,other,v=prepare_axis(a,j);assert status==other==coverage.iloc[j].status;np.testing.assert_array_equal(valid,v);row=coverage.iloc[j];assert row.author_records==int(v.sum()) and row.nonMPR==int((a['response'][v]==0).sum()) and row.pCR_or_MPR==int((a['response'][v]==1).sum()) and row.profile_blocks==len(set(a['profile_blocks'][v])) and row.strict_detection_supported==int(supported[:,j].sum()) and row.measured_true_zero_scores==int(np.sum(endpoint[:,j]&(score[:,j]==0)))
        if h is not None:
            np.testing.assert_allclose(g['row'],h['row'],atol=1e-9,rtol=1e-9);np.testing.assert_allclose(g['score'],-h['delta_rows'],atol=1e-8,rtol=1e-8);maximum=max(maximum,float(np.max(np.abs(g['score']+h['delta_rows']))));assert row.residual_df==m['df'] and row.nuisance_rank==m['nuisance_rank']
        z=ex[ex.axis_id.eq(axis)];assert z.author_record_id.tolist()==cases;np.testing.assert_array_equal(z.in_main_design,v);C,names=controls(a,j)
        for i,reason in enumerate(z.reasons):
            expected=[]
            for label,k in [('SENDER',a['sender_index'][j]),('RECEIVER',a['receiver_index'][j])]:
                if ids[i,k]<0:expected.append(label+'_AUTHOR_UNIT_ABSENT')
                elif cells[i,k]<20:expected.append(label+'_FEWER_THAN20_CELLS')
            if not np.isfinite(a['response'][i]):expected.append('MISSING_ORIGINAL_RESPONSE')
            if not np.isfinite(C[i]).all():expected.append('MISSING_FIXED_CONTROLS')
            if not np.isfinite(standard[i,j]) and endpoint[i,j]:expected.append('CONSTANT_OR_INVALID_BLIND_SCORE_SCALE')
            assert reason==('|'.join(expected) if expected else 'INCLUDED')
    receipt=dict(status='ACCEPTED_ALL92_NATIVE_SCORES_FULL_BACKGROUND_CONTROLS_BLOCKS_EXCLUSIONS_AND_DELETED_QR_GEOMETRY_NO_TRUE_EFFECTS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(src/'private_terminal.json'),resource_guard=guard,candidate_slots=92,all_record_axis_positions=11868,maximum_independent_deleted_QR_geometry_difference=maximum,true_clinical_effects_read=False,scientific_outputs_public=False)
    save('private_terminal.json',receipt);print(json.dumps(dict(status=receipt['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

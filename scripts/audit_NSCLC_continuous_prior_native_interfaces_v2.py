"""All fixed interface inputs; independent sums, intersections and QR scores."""
import gzip
import json
import time
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.stats import pearsonr

from NSCLC_native_MFA_interpretation_inputs_v1 import load,sha,VIEWS
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from prepare_NSCLC_fixed_programme_continuous_priors_v1 import read_rds


def main():
    data=load();root=data['root'];resource(root,1048576)
    out=root/'NSCLC_continuous_prior_native_interfaces_20261016_v2'
    job=root/'HNSCC_receiver_followup_job_control_20261016_v2/prepare_NSCLC_continuous_prior_native_interfaces_v2'
    receipt=json.loads((job/'launch_receipt.json').read_text());deadline=time.monotonic()+21600
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists(),'Input producer missing without terminal'
        fields=p.read_text().rsplit(')',1)[1].split();assert fields[0]!='Z' and fields[19]==str(receipt['proc_start_ticks'])
        assert time.monotonic()<deadline,'Observation expired, not a termination claim';time.sleep(5)
    terminal=json.loads((out/'private_terminal.json').read_text());plan=json.loads((out/'plan_before_native_interface_measurement.json').read_text())
    assert terminal['source_sha256']==receipt['source_sha256']
    for name,item in terminal['files'].items():assert sha(out/name)==item['sha256']
    native=data['native'];units=data['units'];gene=pd.Index(data['genes']);train=data['train'];nc=len(data['cases'])
    z=np.load(out/'private_all_native_interface_measurements.npz',allow_pickle=False)
    assert np.array_equal(z['cases'],data['cases'].sampleID.to_numpy(str)) and np.array_equal(z['train'],train)
    assert z['views'].tolist()==VIEWS
    case_index={s:i for i,s in enumerate(z['cases'])};expected_unit=np.full(len(units),-1,np.int64)
    for vi,view in enumerate(VIEWS):
        for u in units[units.sampleID.isin(case_index)&units.sub_cell_type.isin(data['labels'][view])].itertuples():expected_unit[u.unit_index]=vi*nc+case_index[str(u.sampleID)]
    assert np.array_equal(expected_unit,z['native_unit_to_case_pool'])
    chosen=np.flatnonzero(expected_unit>=0);group=expected_unit[chosen];totals=np.bincount(group,weights=units.n_cells.to_numpy()[chosen],minlength=4*nc).astype(np.int64)
    assert np.array_equal(totals.reshape(4,nc),z['cells']) and (totals>=20).all()
    hits=np.load(native/'positive_cell_counts.npy',mmap_mode='r');means=np.load(native/'mean_logCP10K.npy',mmap_mode='r')
    assert sha(native/'positive_cell_counts.npy')==plan['native_positive_sha256']
    axis=pd.read_csv(out/'private_ligand_axis.tsv',sep='\t');assert axis.gene.tolist()==z['ligands'].tolist() and axis.native_column.tolist()==gene.get_indexer(axis.gene).tolist()
    positives=np.zeros((4*nc,len(axis)),np.int64);sums=np.zeros_like(positives,dtype=float)
    measured=axis.native_column.to_numpy()>=0
    assert np.array_equal(measured,z['ligand_native_measured'])
    observed_positive=np.zeros((4*nc,int(measured.sum())),np.int64)
    observed_sums=np.zeros_like(observed_positive,dtype=float)
    np.add.at(observed_positive,group,np.asarray(hits[np.ix_(chosen,axis.native_column.to_numpy()[measured])]))
    np.add.at(observed_sums,group,np.asarray(means[np.ix_(chosen,axis.native_column.to_numpy()[measured])])*units.n_cells.to_numpy()[chosen,None])
    assert np.array_equal(observed_positive.reshape(4,nc,-1),z['ligand_positive'][:,:,measured])
    assert np.allclose((observed_sums/totals[:,None]).reshape(4,nc,-1),z['ligand_mean'][:,:,measured],atol=1e-10,rtol=1e-10)
    assert np.isnan(z['ligand_positive'][:,:,~measured]).all() and np.isnan(z['ligand_mean'][:,:,~measured]).all()
    panel=pd.read_csv(out/'private_composite_panel_axis.tsv',sep='\t');assert panel.native_column.tolist()==gene.get_indexer(panel.gene).tolist()
    cell_units=np.load(native/'private_cell_unit_index.npy',mmap_mode='r');assert sha(native/'private_cell_unit_index.npy')==plan['native_cell_units_sha256']
    one=np.unpackbits(np.load(out/'private_positive_panel_packed.npy'),axis=1,count=len(panel),bitorder='little').astype(bool)
    two=np.unpackbits(np.load(out/'private_positive2_panel_packed.npy'),axis=1,count=len(panel),bitorder='little').astype(bool)
    assert one.shape==two.shape==(len(cell_units),len(panel)) and (two<=one).all()
    for j,col in enumerate(panel.native_column):
        counts=np.bincount(cell_units[one[:,j]],minlength=len(units)).astype(np.int64)
        assert np.array_equal(counts,hits[:,col])
    fixed=np.load(out/'private_fixed_prefix128_audit_cells.npy');assert np.array_equal(fixed,np.sort(np.random.default_rng(202610162120).choice(min(4096,len(cell_units)),128,replace=False)))
    lookup={int(c):j for j,c in enumerate(fixed)};cols={int(c):j for j,c in enumerate(panel.native_column)};raw_expected=np.zeros((len(fixed),len(panel)),np.int64)
    download=json.loads((root/'count_matrix_download_v1/download_state.json').read_text());assert terminal['source_compressed_sha256']==download['sha256']==plan['matrix_expected_sha256']
    with gzip.open(download['path'],'rt') as stream:
        assert stream.readline().startswith('%%MatrixMarket matrix coordinate');line=stream.readline()
        while line.startswith('%'):line=stream.readline()
        assert tuple(map(int,line.split()))==(len(cell_units),len(gene),terminal['source_coordinates'])
        reader=pd.read_csv(stream,sep=r'\s+',header=None,names=['row','col','count'],dtype=np.int64,chunksize=1000000,engine='c')
        for f in reader:
            selected=f[f.row.isin(fixed+1)&f.col.isin(panel.native_column.to_numpy()+1)]
            for r in selected.itertuples(index=False):raw_expected[lookup[int(r.row)-1],cols[int(r.col)-1]]=int(r.count)
            if int(f.row.iloc[-1])>int(fixed.max())+1:break
        reader.close()
    assert np.array_equal(raw_expected>0,one[fixed]) and np.array_equal(raw_expected>=2,two[fixed])
    cp=expected_unit[cell_units];inside=cp>=0;pi={g:j for j,g in enumerate(panel.gene)}
    for j,receptor in enumerate(z['receptors']):
        parts=str(receptor).split('_')
        if len(parts)==1:
            expected=np.zeros(4*nc,np.int64);np.add.at(expected,group,hits[chosen,gene.get_loc(parts[0])])
            assert np.array_equal(expected.reshape(4,nc),z['receptor_positive'][:,:,j]) and np.isnan(z['receptor_positive2'][:,:,j]).all()
        else:
            positions=[pi[g] for g in parts]
            for arr,key in [(one,'receptor_positive'),(two,'receptor_positive2')]:
                selected=inside.copy()
                for col in positions:selected &= arr[:,col]
                expected=pd.Series(cp[selected]).value_counts().reindex(range(4*nc),fill_value=0).to_numpy(np.int64).reshape(4,nc)
                assert np.array_equal(expected,z[key][:,:,j])
    prior=root/'NSCLC_fixed_programme_continuous_prior_inputs_20261016_v1'
    assert sha(prior/'private_independent_continuous_prior_audit.json')==plan['prior_audit_sha256']
    original=pd.read_csv(prior/'private_all1600_prospective_ligand_slots.tsv.gz',sep='\t',keep_default_na=False);original['slot_index']=np.arange(1600)
    slots=pd.read_csv(out/'private_all1600_original_slots.tsv.gz',sep='\t',keep_default_na=False)
    pd.testing.assert_frame_equal(original,slots)
    flat=read_rds(root/'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1/verified_prior_inputs_v2/lr_network_human_21122021.rds')
    reference=pd.read_csv(root/'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv',sep='\t',keep_default_na=False)
    fids=pd.read_csv(out/'private_all_occupied_response_fidelity.tsv.gz',sep='\t').set_index('slot_index');maximum=0.
    for vi,view in enumerate(VIEWS):
        p=data['prepared'][vi];w=data['model']['weights_'+str(vi)];sw=np.sqrt(data['model']['precision_'+str(vi)])
        q,r=np.linalg.qr(np.vstack([sw[:,None]*w,np.eye(10)]),mode='reduced');a=solve_triangular(r,q[:len(w)].T*sw)
        names=data['genes'][p['selected_native_columns']];original_score=data['model']['own_view_scores_'+str(vi)]
        for ligand,g in slots[slots.view.eq(view)&slots.ligand.ne('')].groupby('ligand'):
            strings={ligand,*flat.loc[flat['from'].astype(str).eq(ligand),'to'].astype(str),*reference.loc[reference.ligand.eq(ligand),'receptor']};omit=set(strings)
            for text in strings:omit.update(text.split('_'))
            keep=np.array([str(s) not in omit for s in names]);predicted=p['X'][:,keep]@a[:,keep].T
            for slot in g.itertuples():
                k=int(slot.factor.rsplit('F',1)[1])-1;observed=fids.loc[slot.slot_index];ok=True
                assert observed.retained_selected_features==keep.sum()
                for tag,take in [('discovery',train),('original_holdout',~train)]:
                    rho=float(pearsonr(predicted[take,k],original_score[take,k]).statistic)
                    rmse=float(np.sqrt(np.mean((predicted[take,k]-original_score[take,k])**2))/np.std(original_score[train,k]))
                    assert np.isclose(observed[tag+'_rho'],rho,atol=1e-9,rtol=1e-9) and np.isclose(observed[tag+'_RMSE_over_discovery_SD'],rmse,atol=1e-9,rtol=1e-9)
                    ok &= np.isfinite(rho) and np.isfinite(rmse) and rho>=.95 and rmse<=.25
                assert observed.response_preserves_original_score==ok and z['response_ok'][slot.slot_index]==ok
                assert np.allclose(predicted[:,k],z['response'][slot.slot_index],atol=1e-9,rtol=1e-9)
                maximum=max(maximum,float(np.max(abs(predicted[:,k]-z['response'][slot.slot_index]))))
    empty=slots.ligand.eq('').to_numpy();assert np.isnan(z['response'][empty]).all() and not z['response_ok'][empty].any()
    interfaces=pd.read_csv(out/'private_all4800_native_interface_readiness.tsv.gz',sep='\t',keep_default_na=False)
    edges=pd.read_csv(out/'private_all_fixed_reference_edge_joint_RNA.tsv.gz',sep='\t',keep_default_na=False)
    li={g:j for j,g in enumerate(z['ligands'])};ri={r:j for j,r in enumerate(z['receptors'])};vi={v:j for j,v in enumerate(VIEWS)}
    totals=z['cells'][:,train].sum(axis=1);lp=z['ligand_positive'][:,train,:].sum(axis=1);rp=z['receptor_positive'][:,train,:].sum(axis=1)
    expected_edges=[];expected_interfaces=[]
    for slot in slots.itertuples():
        recv=vi[slot.view];valid=[]
        if slot.ligand:
            for receptor in slot.compatible_receptors.split('|'):
                count=int(rp[recv,ri[receptor]]);gate=count*100>=totals[recv]
                expected_edges.append((slot.slot_index,receptor,count,int(totals[recv]),gate,count*10>=totals[recv]))
                if gate:valid.append(receptor)
        for sender in VIEWS:
            if sender==slot.view:continue
            reasons=[]
            if not slot.ligand:reasons.append('EMPTY_FROZEN_PRIOR_SLOT')
            else:
                if not measured[li[slot.ligand]]:reasons.append('LIGAND_NOT_PRESENT_IN_NATIVE_RNA_AXIS')
                elif lp[vi[sender],li[slot.ligand]]*100<totals[vi[sender]]:reasons.append('DISCOVERY_SENDER_RNA_BELOW1PCT')
                if not valid:reasons.append('DISCOVERY_SAME_CELL_RECEPTOR_RNA_BELOW1PCT')
                if not z['response_ok'][slot.slot_index]:reasons.append('LR_EXCLUDED_RESPONSE_FIDELITY_FAILED')
            expected_interfaces.append((slot.slot_index,sender,slot.view,slot.factor,int(slot.direction),int(slot.prospective_position),slot.ligand,li.get(slot.ligand,-1),'|'.join(valid),';'.join(reasons) if reasons else 'RNA_INPUT_READY_EFFECT_MODEL_NOT_RUN'))
    assert len(interfaces)==4800 and interfaces.interface_index.tolist()==list(range(4800))
    cols=['slot_index','sender_view','receiver_view','factor','direction','prospective_position','ligand','ligand_axis_index','receiver_RNA_compatible_receptors','status']
    assert list(interfaces[cols].itertuples(index=False,name=None))==expected_interfaces
    cols=['slot_index','receptor','pooled_discovery_joint_positive','pooled_discovery_cells','joint1pct','joint10pct']
    assert list(edges[cols].itertuples(index=False,name=None))==expected_edges
    assert not interfaces.per_case_expression_filter_applied.any() and interfaces.discovery_cases.eq(145).all() and interfaces.original_holdout_cases.eq(63).all()
    report=dict(status='PASS_ALL4800_NATIVE_INTERFACE_INPUTS_SPARSE_SUMS_JOINT_RNA_QR_RESPONSES_AND_FIXED_PREFIX128_RAW_CELLS',finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),all_interface_slots_checked=4800,maximum_QR_response_error=maximum,raw_replay_cells=128,raw_replay_is_prefix_bounded=True,all_cells_independently_raw_replayed=False,no_new_association_or_clinical_results=True,direct_communication_validated=False)
    proof=out/'private_independent_native_interface_audit.json';assert not proof.exists();proof.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],audit_sha256=sha(proof))))


if __name__=='__main__':main()

"""All fixed prior positions: native RNA and disjoint receiver input only.

No association, clinical value, source effect or selected gene is printed.
"""
import gzip
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from NSCLC_native_MFA_interpretation_inputs_v1 import load, sha, VIEWS
from prepare_NSCLC_fixed_programme_continuous_priors_v1 import read_rds, PRIORS, calibrate
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from prepare_NSCLC_full_native_gene_units_v1 import TrackedReader, chunks


def main():
    data=load();root=data['root'];budget=resource(root,128*1024**2)
    prior=root/'NSCLC_fixed_programme_continuous_prior_inputs_20261016_v1'
    audit=prior/'private_independent_continuous_prior_audit.json'
    assert sha(audit)=='d87696f9800cf51ef735f6799cfccb0e43ca4c1d032ee9dc4e641274b2cf001a'
    accepted=json.loads(audit.read_text());terminal=json.loads((prior/'private_terminal.json').read_text())
    assert accepted['producer_terminal_sha256']==sha(prior/'private_terminal.json')
    for name,item in terminal['files'].items():assert sha(prior/name)==item['sha256']
    slot_path=prior/'private_all1600_prospective_ligand_slots.tsv.gz'
    slots=pd.read_csv(slot_path,sep='\t',keep_default_na=False);assert len(slots)==1600
    slots['slot_index']=np.arange(len(slots));filled=slots[slots.ligand.ne('')]
    ligands=sorted(set(filled.ligand));receptors=sorted({r for v in filled.compatible_receptors for r in v.split('|') if r})
    multicomponent=[r for r in receptors if '_' in r];panel=sorted({g for r in multicomponent for g in r.split('_')})
    out=root/'NSCLC_continuous_prior_native_interfaces_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    native=data['native'];units=data['units'];genes=pd.Index(data['genes']);cases=data['cases'].sampleID.to_numpy(str);train=data['train']
    positive=native/'positive_cell_counts.npy';native_terminal=json.loads((native/'private_terminal.json').read_text())
    assert sha(positive)==native_terminal['outputs']['positive_cell_counts.npy']['sha256']
    cell_units_path=native/'private_cell_unit_index.npy';cell_units=np.load(cell_units_path,mmap_mode='r');ncells=len(cell_units)
    assert np.array_equal(np.bincount(cell_units,minlength=len(units)),units.n_cells)
    download=json.loads((root/'count_matrix_download_v1/download_state.json').read_text());matrix=Path(download['path'])
    assert matrix.stat().st_size==download['bytes']
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
              prior_audit_sha256=sha(audit),prior_slot_sha256=sha(slot_path),native_terminal_sha256=sha(native/'private_terminal.json'),
              native_positive_sha256=sha(positive),native_cell_units_sha256=sha(cell_units_path),matrix_expected_sha256=download['sha256'],
              scope='All1600 fixed prior positions x3 distinct sender pools; all208 original cases, unchanged145/63 partitions. No rank21 replacement.',
              sender='Cell-weighted original mean per-cell log1p(CP10K) plus positive cell counts; native full RNA denominator. Real zeros retained.',
              receptor='All components in same cell count>=1 main and >=2 composite sensitivity; singleton positives from accepted native cache. Every compatible reference edge retained.',
              response='Original frozen own-view score after deleting ligand and all literal flat/reference receptor components. Same discovery and63 fidelity rho>=.95/RMSE<=.25SD gates.',
              RNA_gate='Discovery-only pooled sender positivity>=1pct; at least one fixed compatible receiver edge jointly positive>=1pct in discovery. Report10pct; no per-case expression-based exclusion.',
              no_correlations_P_or_clinical_values=True,clinical_or_direct_communication_validation=False)
    (out/'plan_before_native_interface_measurement.json').write_text(json.dumps(plan,indent=2)+'\n')
    pd.DataFrame(dict(gene=ligands,native_column=genes.get_indexer(ligands))).to_csv(out/'private_ligand_axis.tsv',sep='\t',index=False)
    panel_columns=genes.get_indexer(panel);assert (panel_columns>=0).all() and (genes.get_indexer(ligands)>=0).all()
    available=int(next(s for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith('MemAvailable:')).split()[1])*1024
    assert available>ncells*len(panel)*3+2*1024**3
    pd.DataFrame(dict(gene=panel,native_column=panel_columns)).to_csv(out/'private_composite_panel_axis.tsv',sep='\t',index=False)
    pd.DataFrame(dict(receptor=receptors)).to_csv(out/'private_receptor_axis.tsv',sep='\t',index=False)
    slots.to_csv(out/'private_all1600_original_slots.tsv.gz',sep='\t',index=False,compression='gzip')
    # Integer sparse membership independently auditable without cell weights.
    case_index={s:i for i,s in enumerate(cases)};row=[];col=[]
    native_to_pool=np.full(len(units),-1,np.int64)
    for vi,view in enumerate(VIEWS):
        take=units.sampleID.isin(case_index)&units.sub_cell_type.isin(data['labels'][view])
        for u in units[take].itertuples():
            ri=vi*len(cases)+case_index[str(u.sampleID)];row.append(ri);col.append(int(u.unit_index))
            assert native_to_pool[u.unit_index]==-1;native_to_pool[u.unit_index]=ri
    membership=csr_matrix((np.ones(len(row),np.int64),(row,col)),shape=(4*len(cases),len(units)))
    totals=np.asarray(membership@units.n_cells.to_numpy(np.int64)).ravel();assert (totals>=20).all()
    weights=membership.multiply(1/totals[:,None]).multiply(units.n_cells.to_numpy()[None,:]).tocsr()
    hits=np.load(positive,mmap_mode='r');means=np.load(native/'mean_logCP10K.npy',mmap_mode='r')
    lc=genes.get_indexer(ligands)
    ligand_positive=np.asarray(membership@hits[:,lc]).reshape(4,len(cases),len(ligands))
    ligand_mean=np.asarray(weights@means[:,lc]).reshape(4,len(cases),len(ligands))
    assert np.isfinite(ligand_mean).all() and (ligand_positive<=totals.reshape(4,len(cases),1)).all()
    # Full source SHA and all-cell panel read, no new network download.
    one=np.zeros((ncells,len(panel)),bool);two=np.zeros_like(one);lookup=np.full(len(genes),-1,int);lookup[panel_columns]=np.arange(len(panel))
    seen=0;previous=-1
    with matrix.open('rb') as raw:
        tracked=TrackedReader(raw)
        with gzip.GzipFile(fileobj=tracked,mode='rb') as stream:
            assert stream.readline().startswith(b'%%MatrixMarket matrix coordinate');line=stream.readline()
            while line.startswith(b'%'):line=stream.readline()
            nr,ng,nnz=map(int,line.split());assert nr==ncells and ng==len(genes)
            for block in chunks(stream):
                assert np.isfinite(block).all() and np.array_equal(block,np.floor(block)) and (block[:,2]>0).all()
                rr=block[:,0].astype(np.int64)-1;cc=block[:,1].astype(np.int64)-1;vv=block[:,2]
                key=rr*ng+cc;assert key[0]>previous and (np.diff(key)>0).all();previous=int(key[-1]);seen+=len(block)
                local=lookup[cc];take=local>=0;one[rr[take],local[take]]=True;two[rr[take],local[take]]=vv[take]>=2
                if seen//50000000!=(seen-len(block))//50000000:
                    (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),stage='native_joint_panel',coordinates_read=seen,total_coordinates=nnz))+'\n')
        assert seen==nnz and tracked.bytes==download['bytes'] and tracked.digest.hexdigest()==download['sha256']
    np.save(out/'private_positive_panel_packed.npy',np.packbits(one,axis=1,bitorder='little'))
    np.save(out/'private_positive2_panel_packed.npy',np.packbits(two,axis=1,bitorder='little'))
    fixed_cells=np.sort(np.random.default_rng(202610162120).choice(min(4096,ncells),128,replace=False))
    np.save(out/'private_fixed_prefix128_audit_cells.npy',fixed_cells)
    panel_index={g:i for i,g in enumerate(panel)};cell_pool=native_to_pool[cell_units];in_case=cell_pool>=0
    receptor_positive=np.zeros((4*len(cases),len(receptors)),np.int64);receptor_positive2=np.full(receptor_positive.shape,np.nan)
    for j,receptor in enumerate(receptors):
        parts=receptor.split('_')
        if len(parts)==1:receptor_positive[:,j]=np.asarray(membership@hits[:,genes.get_loc(parts[0])]).ravel()
        else:
            columns=[panel_index[p] for p in parts]
            for arr,dest in [(one,receptor_positive),(two,receptor_positive2)]:
                joint=arr[:,columns].all(axis=1)&in_case
                dest[:,j]=np.bincount(cell_pool[joint],minlength=4*len(cases))
    receptor_positive=receptor_positive.reshape(4,len(cases),len(receptors));receptor_positive2=receptor_positive2.reshape(receptor_positive.shape)
    assert (receptor_positive<=totals.reshape(4,len(cases),1)).all()
    resource_dir=root/'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1/verified_prior_inputs_v2'
    assert sha(resource_dir/'lr_network_human_21122021.rds')==PRIORS['lr_network_human_21122021.rds']
    flat=read_rds(resource_dir/'lr_network_human_21122021.rds');flat['from']=flat['from'].astype(str);flat['to']=flat['to'].astype(str)
    refpath=root/'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv'
    assert sha(refpath)=='df574533d866438e646423f29bd555cb283173dbceadfe9f8264e1622a86fd06'
    reference=pd.read_csv(refpath,sep='\t',keep_default_na=False)
    response=np.full((1600,len(cases)),np.nan);response_ok=np.zeros(1600,bool);fidelity_rows=[]
    for vi,view in enumerate(VIEWS):
        p=data['prepared'][vi];w=data['model']['weights_'+str(vi)];tau=data['model']['precision_'+str(vi)]
        a=np.linalg.solve(np.eye(10)+w.T@(tau[:,None]*w),w.T*tau);names=data['genes'][p['selected_native_columns']]
        original=data['model']['own_view_scores_'+str(vi)]
        for ligand,g in filled[filled.view.eq(view)].groupby('ligand'):
            strings={ligand,*flat.loc[flat['from'].eq(ligand),'to'],*reference.loc[reference.ligand.eq(ligand),'receptor']}
            omitted=set(strings)
            for text in strings:omitted.update(str(text).split('_'))
            keep=np.array([str(gene) not in omitted for gene in names]);predicted=p['X'][:,keep]@a[:,keep].T
            values,gate=calibrate(predicted,original,train)
            for slot in g.itertuples():
                k=int(slot.factor.rsplit('F',1)[1])-1;response[slot.slot_index]=predicted[:,k];response_ok[slot.slot_index]=gate[k]
                fidelity_rows.append(dict(slot_index=slot.slot_index,view=view,factor=slot.factor,ligand=ligand,retained_selected_features=int(keep.sum()),response_preserves_original_score=bool(gate[k]),**{field:float(v[k]) for field,v in values.items()}))
    li={g:i for i,g in enumerate(ligands)};ri={r:i for i,r in enumerate(receptors)};view_index={v:i for i,v in enumerate(VIEWS)}
    total=totals.reshape(4,len(cases));discovery_total=total[:,train].sum(axis=1)
    lp=ligand_positive[:,train,:].sum(axis=1);rp=receptor_positive[:,train,:].sum(axis=1)
    interfaces=[];edge_rows=[]
    for slot in slots.itertuples():
        receiver=view_index[slot.view];compatible=[]
        if slot.ligand:
            for receptor in slot.compatible_receptors.split('|'):
                j=ri[receptor];positive_count=int(rp[receiver,j]);ok=positive_count*100>=discovery_total[receiver]
                edge_rows.append(dict(slot_index=slot.slot_index,receiver_view=slot.view,receptor=receptor,pooled_discovery_joint_positive=positive_count,pooled_discovery_cells=int(discovery_total[receiver]),joint1pct=bool(ok),joint10pct=bool(positive_count*10>=discovery_total[receiver]),same_cell_RNA_not_protein_assembly=True))
                if ok:compatible.append(receptor)
        for sender,view in enumerate(VIEWS):
            if sender==receiver:continue
            reasons=[]
            if not slot.ligand:reasons.append('EMPTY_FROZEN_PRIOR_SLOT')
            else:
                if lp[sender,li[slot.ligand]]*100<discovery_total[sender]:reasons.append('DISCOVERY_SENDER_RNA_BELOW1PCT')
                if not compatible:reasons.append('DISCOVERY_SAME_CELL_RECEPTOR_RNA_BELOW1PCT')
                if not response_ok[slot.slot_index]:reasons.append('LR_EXCLUDED_RESPONSE_FIDELITY_FAILED')
            interfaces.append(dict(interface_index=len(interfaces),slot_index=slot.slot_index,sender_view=view,receiver_view=slot.view,factor=slot.factor,direction=slot.direction,prospective_position=slot.prospective_position,ligand=slot.ligand,ligand_axis_index=li.get(slot.ligand,-1),receiver_RNA_compatible_receptors='|'.join(compatible),status=';'.join(reasons) if reasons else 'RNA_INPUT_READY_EFFECT_MODEL_NOT_RUN',discovery_cases=int(train.sum()),original_holdout_cases=int((~train).sum()),per_case_expression_filter_applied=False))
    assert len(interfaces)==4800
    pd.DataFrame(interfaces).to_csv(out/'private_all4800_native_interface_readiness.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(edge_rows).to_csv(out/'private_all_fixed_reference_edge_joint_RNA.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(fidelity_rows).to_csv(out/'private_all_occupied_response_fidelity.tsv.gz',sep='\t',index=False,compression='gzip')
    np.savez_compressed(out/'private_all_native_interface_measurements.npz',ligand_mean=ligand_mean,ligand_positive=ligand_positive,receptor_positive=receptor_positive,receptor_positive2=receptor_positive2,response=response,response_ok=response_ok,original_scores=np.stack([data['model']['own_view_scores_'+str(v)] for v in range(4)]),cells=total,cases=cases,train=train,views=np.array(VIEWS),ligands=np.array(ligands),receptors=np.array(receptors),native_unit_to_case_pool=native_to_pool)
    result=dict(status='COMPLETE_ALL4800_NATIVE_INTERFACE_INPUTS_NO_ASSOCIATIONS_INDEPENDENT_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),prior_audit_sha256=sha(audit),source_compressed_sha256=tracked.digest.hexdigest(),source_coordinates=seen,
                files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},clinical_or_direct_communication_validation=False,new_association_or_P_calculated=False)
    assert sum(v['bytes'] for v in result['files'].values())<budget['output_cap_bytes']
    (out/'private_terminal.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],finished_UTC=result['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

"""Single passed RNA arm, conditional on original author T labels, no inference.

This never overrides the failed lineage bridge or admits dual-arm CM02 spatial
tests. Original author labels plus measured pan-T support define description.
"""
import argparse
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from fit_CM02_panel_bridge_LOSO_v1 import choose, fit, prepare, save_bundle, weights
from prepare_CM02_panel_bridge_reference_scope_v1 import budget

AUTHOR_T=['CD4 T','CD8 T','CXCL13+ CD4 T','CXCL13+ CD8 T','TCF7+ CD4 T',
          'TCF7+ CD8 T','TCF7+PD1+ CD8 T','Treg']


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()


def main(a):
    budget(a.root,2*1024**3)
    validation=a.root/'model_validation_v1'
    audit=a.root/'model_independent_audit_v2/independent_model_audit_technical_public.json'
    assert json.loads(audit.read_text())['producer_terminal_sha256']==sha(validation/'model_validation_technical_public.json')
    gates=pd.read_csv(validation/'private_all_frozen_engineering_gates.tsv',sep='\t')
    assert gates[gates.target.eq('cytotoxic')].passed.all()
    assert not gates[gates.target.eq('regulatory')].passed.all()
    assert not gates[gates.target.eq('T')].passed.all()
    accepted=json.loads((a.source/'download_technical_public.json').read_text())
    for row in accepted['files']:
        p=a.source/row['relative_path'];assert p.stat().st_size==row['bytes'] and sha(p)==row['sha256']
    out=a.root/'cytotoxic_author_T_description_v1';assert not out.exists();out.mkdir(mode=0o700)
    contract=dict(created_utc=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(Path(__file__)),
                  original_author_fine_T_whitelist=AUTHOR_T,pan_T_support='>=2 of CD3D/CD3E/TRAC with count>=1',
                  target='only cytotoxic-associated RNA, passed native and stress reference gates',
                  failed_lineage_bridge_not_deployed=True,failed_regulatory_arm_not_projected=True,
                  author_identity_is_conditional_not_independent_identity_validation=True,
                  final_ridge_selection='same joint-target inner3-case-fold alpha procedure on all reference T; export only passed cytotoxic output',
                  reference_domain='case-weighted99th percentile of standardized permitted-feature squared norms in fixed matched-depth reference T',
                  out_of_range_scores='retained raw; excluded from display, no clipping into apparent qualification',
                  no_new_spatial_statistic_P_or_clinical_model=True,no_TMA_or_cell_contact_claim=True)
    (out/'plan_before_target_predictions.json').write_text(json.dumps(contract,indent=2)+'\n')
    c=json.loads((a.root/'reference_scope_v1/private_frozen_reference_contract.json').read_text())
    ref=a.root/'reference_native_inputs_v2';meta=pd.read_csv(ref/'private_selected_reference_metadata.tsv.gz',sep='\t',keep_default_na=False)
    raw=np.load(ref/'private_panel_counts.npy',mmap_mode='r');target=np.load(ref/'private_full_native_target_arrays.npz')['original_full_signature_UCell1500']
    cols=[c['shared_panel_genes'].index(g) for g in c['programme_feature_genes']]
    x,positive=prepare(raw,cols);ids=np.flatnonzero(meta.lineage.eq('T').to_numpy() & positive & np.isfinite(target).all(1))
    groups=meta.patient_block.astype(str).to_numpy()[ids]
    alpha,inner=choose(x[ids],target[ids],groups,'programme');assert alpha is not None
    bundle=fit(x[ids],target[ids],groups,alpha,'programme');assert bundle[3]
    save_bundle(out/'private_final_programme_model_only_cytotoxic_used.npz',bundle,cols,ids,np.array([],int))
    pd.DataFrame(inner).to_csv(out/'private_final_reference_only_alpha_selection.tsv',sep='\t',index=False)
    model,scaler,keep,_=bundle
    stress=np.load(validation/'private_fixed_capture_stress_counts.npy',mmap_mode='r')
    sx,_=prepare(stress[ids],cols);norm=np.mean(scaler.transform(sx[:,keep])**2,axis=1)
    order=np.argsort(norm);w=weights(groups);cw=np.cumsum(w[order]);limit=float(norm[order[np.searchsorted(cw,.99*cw[-1])]])
    source=a.source/'single_cell_portal'
    names=pd.read_csv(source/'counts/barcodes.tsv.gz',header=None,dtype=str)[0]
    genes=pd.read_csv(source/'counts/features.tsv.gz',header=None,dtype=str)[0].tolist()
    cells=pd.read_csv(source/'metadata.txt.gz',sep='\t',skiprows=[1],dtype=str,keep_default_na=False).set_index('NAME').loc[names].reset_index()
    author=cells.type_lvl1.eq('Lymphocyte') & cells.type_lvl2.isin(AUTHOR_T)
    selected=np.flatnonzero(author);inverse=np.full(len(cells),-1,np.int32);inverse[selected]=np.arange(len(selected))
    data=np.lib.format.open_memmap(out/'private_all_author_T_original479_counts.npy',mode='w+',dtype=np.int32,shape=(len(selected),len(genes)))
    with gzip.open(source/'counts/counts.mtx.gz','rt') as stream:
        stream.readline()
        while True:
            line=stream.readline()
            if not line.startswith('%'):break
        for part in pd.read_csv(stream,sep=r'\s+',header=None,names=['gene','cell','value'],chunksize=1000000):
            ri=inverse[part.cell.to_numpy(int)-1];take=ri>=0
            np.add.at(data,(ri[take],part.gene.to_numpy(int)[take]-1),part.value.to_numpy(np.int32)[take])
    data.flush();table=cells.iloc[selected].copy().reset_index(drop=True)
    assert np.array_equal(np.asarray(data.sum(1)).ravel(),pd.to_numeric(table.ntx))
    assert np.array_equal(np.asarray((data>0).sum(1)).ravel(),pd.to_numeric(table.ngene))
    positions=[genes.index(g) for g in c['programme_feature_genes']]
    tx,nonzero=prepare(data,positions);z=scaler.transform(tx[:,keep]);pred=model.predict(z)[:,1]
    domain=np.mean(z*z,axis=1);pan=np.asarray((data[:,[genes.index(g) for g in ['CD3D','CD3E','TRAC']]]>0).sum(1)).ravel()
    state=np.ones(len(selected),np.uint8);state[pan>=2]=2
    state[(pan>=2)&nonzero]=3
    state[(pan>=2)&nonzero&(domain<=limit)]=4
    valid=(pan>=2)&nonzero&(domain<=limit)&(pred>=0)&(pred<=1);state[valid]=5
    table['original_matrix_column']=selected;table['pan_T_genes_detected']=pan;table['permitted_features_positive']=nonzero
    table['reference_domain_squared_norm']=domain;table['reference_domain_limit']=limit
    table['raw_cytotoxic_programme_prediction']=pred
    table['cytotoxic_RNA_description']=np.where(valid,pred,np.nan);table['description_status_code']=state
    table['identity_source']='original_author_T_only_not_failed_lineage_model'
    table.to_csv(out/'private_author_T_only_cytotoxic_RNA_description.tsv.gz',sep='\t',index=False,compression='gzip')
    allstates=np.zeros(len(cells),np.uint8);allstates[selected]=state
    np.save(out/'private_all_original_cell_description_state.npy',allstates)
    summary=table.groupby('donor_id',sort=True).agg(author_T_cells=('NAME','size'),display_eligible=('cytotoxic_RNA_description','count'))
    summary.to_csv(out/'private_all4_author_T_conditional_description_coverage.tsv',sep='\t')
    note=dict(status='CONDITIONAL_SINGLE_CYTOTOXIC_ARM_DESCRIPTION_NO_SPATIAL_ADMISSION',
              finished_utc=datetime.now(timezone.utc).isoformat(),all_original_cells=len(cells),
              author_T_scope_cells=len(selected),case_count=4,
              coordinate_units='published_coordinate_units_not_new_physical_contact_measurement',
              author_T_label_is_not_independently_validated_by_this_projection=True,
              failed_lineage_and_regulatory_gates_kept=True,no_new_spatial_P=True,
              producer_sha256=sha(Path(__file__)),frozen_contract_sha256=sha(out/'plan_before_target_predictions.json'),
              input_download_receipt_sha256=sha(a.source/'download_technical_public.json'),
              files={p.name:sha(p) for p in out.iterdir() if p.is_file()})
    (out/'description_technical_public.json').write_text(json.dumps(note,indent=2)+'\n');print(json.dumps({k:v for k,v in note.items() if k!='files'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--source',type=Path,required=True);main(p.parse_args())

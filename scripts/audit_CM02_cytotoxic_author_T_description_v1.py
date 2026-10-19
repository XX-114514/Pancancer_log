"""Independent coefficient replay and complete conditional projection masks."""
import argparse
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd

T_LABELS={'CD4 T','CD8 T','CXCL13+ CD4 T','CXCL13+ CD8 T','TCF7+ CD4 T',
          'TCF7+ CD8 T','TCF7+PD1+ CD8 T','Treg'}

def main(a):
    out=a.root/'cytotoxic_author_T_description_v1'
    producer=out/'description_technical_public.json'
    proof=json.loads(producer.read_text())
    for name,digest in proof['files'].items():
        h=hashlib.sha256()
        with (out/name).open('rb') as f:
            for b in iter(lambda:f.read(4194304),b''):h.update(b)
        assert h.hexdigest()==digest
    source=a.source/'single_cell_portal'
    names=pd.read_csv(source/'counts/barcodes.tsv.gz',header=None,dtype=str)[0]
    genes=pd.read_csv(source/'counts/features.tsv.gz',header=None,dtype=str)[0].tolist()
    metadata=pd.read_csv(source/'metadata.txt.gz',sep='\t',skiprows=[1],dtype=str,keep_default_na=False).set_index('NAME').loc[names].reset_index()
    mask=metadata.type_lvl1.eq('Lymphocyte') & metadata.type_lvl2.isin(T_LABELS)
    rows=np.flatnonzero(mask);inverse=np.full(len(metadata),-1,np.int32);inverse[rows]=np.arange(len(rows))
    original=np.load(out/'private_all_author_T_original479_counts.npy',mmap_mode='r')
    # All raw COO entries belonging to these cells are checked against the saved
    # matrix; the accepted source has independently proven no duplicate pairs.
    found=0
    with gzip.open(source/'counts/counts.mtx.gz','rt') as f:
        line=f.readline()
        while True:
            line=f.readline()
            if not line.startswith('%'):break
        for block in pd.read_csv(f,sep=r'\s+',header=None,names=['g','c','v'],chunksize=1000000):
            r=inverse[block.c.to_numpy(int)-1];take=r>=0
            assert np.array_equal(original[r[take],block.g.to_numpy(int)[take]-1],block.v.to_numpy(int)[take])
            found+=int(take.sum())
    assert found==int(np.count_nonzero(original))
    table=pd.read_csv(out/'private_author_T_only_cytotoxic_RNA_description.tsv.gz',sep='\t',float_precision='round_trip')
    assert np.array_equal(table.NAME.to_numpy(),metadata.NAME.iloc[rows].to_numpy())
    assert np.array_equal(table.original_matrix_column.to_numpy(),rows)
    assert np.array_equal(original.sum(1),pd.to_numeric(metadata.ntx.iloc[rows]))
    assert np.array_equal((original>0).sum(1),pd.to_numeric(metadata.ngene.iloc[rows]))
    contract=json.loads((a.root/'reference_scope_v1/private_frozen_reference_contract.json').read_text())
    model=np.load(out/'private_final_programme_model_only_cytotoxic_used.npz')
    allowed=[contract['shared_panel_genes'][int(i)] for i in model['denominator_positions']]
    selected=[contract['shared_panel_genes'][int(i)] for i in model['selected_feature_positions']]
    assert not set(allowed)&set(contract['core14']+contract['LR_excluded_genes'])
    x=np.asarray(original[:,[genes.index(g) for g in selected]],float)
    denom=np.asarray(original[:,[genes.index(g) for g in allowed]],float).sum(1)
    good=denom>0;x[good]=np.log1p(10000*x[good]/denom[good,None]);z=(x-model['mean'])/model['scale']
    pred=z@model['coef'][1]+model['intercept'][1];norm=(z*z).mean(1)
    np.testing.assert_allclose(pred,table.raw_cytotoxic_programme_prediction,atol=2e-12,rtol=0)
    np.testing.assert_allclose(norm,table.reference_domain_squared_norm,atol=2e-10,rtol=0)
    refmeta=pd.read_csv(a.root/'reference_native_inputs_v2/private_selected_reference_metadata.tsv.gz',sep='\t',keep_default_na=False)
    tr=model['training_rows'];assert refmeta.lineage.iloc[tr].eq('T').all()
    stress=np.load(a.root/'model_validation_v1/private_fixed_capture_stress_counts.npy',mmap_mode='r')
    xx=np.asarray(stress[tr][:,model['selected_feature_positions']],float)
    total=stress[tr][:,model['denominator_positions']].sum(1);nonzero=total>0
    xx[nonzero]=np.log1p(10000*xx[nonzero]/total[nonzero,None]);zz=(xx-model['mean'])/model['scale'];rn=(zz*zz).mean(1)
    groups=refmeta.patient_block.iloc[tr];w=1/groups.map(groups.value_counts()).to_numpy(float)
    ix=np.argsort(rn);limit=rn[ix[np.searchsorted(np.cumsum(w[ix]),.99*w.sum())]]
    np.testing.assert_allclose(table.reference_domain_limit,limit,rtol=0,atol=2e-10)
    pan=(original[:,[genes.index(g) for g in ['CD3D','CD3E','TRAC']]]>0).sum(1)
    state=np.ones(len(rows),np.uint8);state[pan>=2]=2;state[(pan>=2)&good]=3
    state[(pan>=2)&good&(norm<=limit)]=4
    valid=(pan>=2)&good&(norm<=limit)&(pred>=0)&(pred<=1);state[valid]=5
    assert np.array_equal(table.description_status_code,state)
    assert np.array_equal(table.cytotoxic_RNA_description.notna(),valid)
    allstates=np.zeros(len(metadata),np.uint8);allstates[rows]=state
    assert np.array_equal(allstates,np.load(out/'private_all_original_cell_description_state.npy'))
    summary=table.groupby('donor_id').agg(author_T_cells=('NAME','size'),display_eligible=('cytotoxic_RNA_description','count'))
    pd.testing.assert_frame_equal(summary,pd.read_csv(out/'private_all4_author_T_conditional_description_coverage.tsv',sep='\t',index_col=0))
    assert len(summary)==4
    receipt=dict(status='ALL_CONDITIONAL_AUTHOR_T_COUNTS_COEFFICIENT_PREDICTIONS_DOMAIN_AND_NA_MASKS_ACCEPTED',
                 finished_utc=datetime.now(timezone.utc).isoformat(),producer_receipt_sha256=hashlib.sha256(producer.read_bytes()).hexdigest(),
                 auditor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 no_new_spatial_inference=True,no_independent_author_identity_validation_claim=True)
    dest=a.root/'cytotoxic_description_independent_audit_v1';assert not dest.exists();dest.mkdir(mode=0o700)
    (dest/'independent_description_audit_technical_public.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--source',type=Path,required=True);main(p.parse_args())

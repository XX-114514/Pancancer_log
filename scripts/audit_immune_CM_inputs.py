#!/usr/bin/env python3
"""Independently recount saved CM frequencies from all annotated cells."""
import argparse,hashlib,json
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem
import zarr
from prepare_immune_CM import material_classes


def main(args):
    obs=read_elem(zarr.open_group(str(args.scratch/'immune_refined.zarr'),mode='r')['obs'])
    yes=lambda v:v.astype(str).str.lower().isin(['true','1'])
    eligible=yes(obs.CM_composition_eligible)&~yes(obs.annotation_exclude_from_primary)&~yes(obs.ex_vivo_stimulation)&~obs.major.astype(str).isin(['Unresolved','T_NK_boundary'])
    car=obs.analysis_role.astype(str).eq('CAR_product');capture=obs.capture_scope.astype(str)
    families={'whole_TME':~car&capture.str.contains('whole_tumor|serial_tumor'),
              'immune_enriched':~car&capture.str.contains('CD45_enriched|lymphocyte_enriched|T_enriched'),
              'mixed_enriched_nuclei':~car&capture.eq('mixed_tumor_nuclei_NFAT_enrichment')}
    masks={k:eligible&v for k,v in families.items()}
    masks['CAR_product']=eligible&car;masks['T_lineage_biopsy']=eligible&~car&obs.major.astype(str).eq('T')
    material=material_classes(obs)
    for family,family_mask in families.items():
        for label in material.unique():
            mask=eligible&family_mask&material.eq(label)
            masks[family+'_'+label]=mask
            masks['T_'+family+'_'+label]=mask&obs.major.astype(str).eq('T')
    scopes=[]
    for name in ['CM_preparation_public.json','CM_material_preparation_public.json']:
        scopes+=json.loads((args.run/'manifests'/name).read_text())['scopes']
    rows=[]
    for scope in scopes:
        if not scope['ready']:continue
        name=scope['scope'];folder=args.run/'CM'/name/'inputs'
        fq=pd.read_csv(folder/'frequency_private.tsv',sep='\t',index_col=0)
        frame=obs.loc[masks[name],['sampleID','subcluster','major','dataset','patient_block']].astype(str)
        sizes=frame.groupby(['subcluster','sampleID'],observed=True).size()
        denominators=frame.groupby(['major','sampleID'],observed=True).size()
        actual=fq.drop(columns=['majorCluster']);expected=np.empty(actual.shape)
        for i,node in enumerate(actual.index):
            major=str(fq.loc[node,'majorCluster'])
            local=frame.loc[frame.subcluster.eq(node),'major'].unique()
            if len(local)!=1 or local[0]!=major:raise ValueError('Saved CM major/node identity differs')
            for j,unit in enumerate(actual.columns):
                denominator=int(denominators.get((major,unit),0))
                if denominator<20:raise ValueError('Unmeasured major entered saved CM frequency')
                expected[i,j]=int(sizes.get((str(node),str(unit)),0))/denominator
        error=float(np.max(np.abs(expected-actual.to_numpy())))
        if error>1e-12:raise ValueError('Saved CM frequencies disagree with actual cells')
        metadata=pd.read_csv(folder/'sample_metadata_private.tsv',sep='\t',keep_default_na=False)
        if not set(actual.columns)==set(metadata.sampleID):raise ValueError('Saved CM metadata units differ')
        for row in metadata.itertuples():
            group=frame.loc[frame.sampleID.eq(row.sampleID)]
            if len(group)!=row.n_cells:raise ValueError('CM input cell count differs')
            if group.dataset.nunique()!=1 or str(group.dataset.iloc[0])!=row.dataset:raise ValueError('CM source differs')
            blocks=group.patient_block.unique()
            if len(blocks)!=1 or blocks[0]!=row.patient_block:raise ValueError('CM patient-block mapping differs')
        rows.append(dict(scope=name,n_nodes=actual.shape[0],n_units=actual.shape[1],max_frequency_error=error,
            input_frequency_sha256=hashlib.sha256((folder/'frequency_private.tsv').read_bytes()).hexdigest(),
            input_metadata_sha256=hashlib.sha256((folder/'sample_metadata_private.tsv').read_bytes()).hexdigest()))
    report=dict(observed_utc=datetime.now(timezone.utc).isoformat(),status='PASS_SAVED_CM_INPUTS_RECOUNTED_FROM_ALL_CELL_METADATA',
        n_atlas_cells=len(obs),all_ready_scopes=len(rows),source_and_patient_blocks_exact=True,scopes=rows,
        not_biological_or_clinical_validation=True)
    with (args.run/'manifests/CM_input_recount_v1_public.json').open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--scratch',type=Path,required=True)
    main(p.parse_args())

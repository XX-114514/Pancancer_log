"""Qualify remaining author namespace candidates; never adopt a guessed alias."""
import json
import os
import re
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from restore_GSE176021_exact_author_annotation_sidecar_v2 import original_frames,author_frame
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2);out=root/'GSE176021_unmatched_author_namespace_qualification_20261019_v2';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_unmatched_namespace_candidates.json').write_text(json.dumps({'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'scope':'Every remaining unmatched original CD3/CD8 author record. Previous exact-patient search has no candidate; search ALL raw namespaces for exact suffix-barcode and cached-native-QC candidates, explicitly retaining patient disagreement. Candidate qualification only, no alias/label or patient mapping adoption; independent original metadata/full native fingerprint proof required before any mapping.','resource_guard':guard},indent=2)+'\n')
    original,raw=original_frames(root);keys=pd.Index(raw.imid+'_'+raw.barcode);side=np.load(root/'GSE176021_exact_original_author_annotation_sidecar_20261019_v2/private_all_atlas_cells_exact_author_bridge.npz');a=np.load(root/'GSE176021_author_lineage_native_RNA_20261019_v1/private_selected_cell_quality_scores.npz');d=pd.DataFrame({'gr':a['global_row'],'RNA':a['fullRNA'],'det':a['detected']});assert d.groupby('gr').RNA.nunique().eq(1).all() and d.groupby('gr').det.nunique().eq(1).all();d=d.drop_duplicates('gr').set_index('gr');rawrow=pd.Series(side['original_raw_row'],index=side['global_atlas_row']);d['rawrow']=rawrow.reindex(d.index).to_numpy();cache=d.set_index('rawrow');candidates=[];aggregate=[]
    norm=lambda s:re.sub('[^a-z0-9]','',s.lower())
    for branch in ['CD3','CD8']:
        author=author_frame(original/('GSE176021_'+branch+'_annotations.rds.gz'));unmatched=author[keys.get_indexer(author.barcode)<0].copy();splits=unmatched.barcode.str.rsplit('_',n=1,expand=True);unmatched['namespace']=splits[0];unmatched['suffix_barcode']=splits[1]
        for namespace,g in unmatched.groupby('namespace',sort=True):
            assert g.imid.nunique()==1;patient=g.imid.iloc[0];rr=raw;full=[]
            for raw_namespace,h in rr.groupby('imid',sort=True):
                assert h.barcode.is_unique;position=pd.Index(h.barcode).get_indexer(g.suffix_barcode);good=position>=0;n=int(good.sum());raw_positions=raw.index.get_indexer(h.index.to_numpy())[position[good]];author_positions=np.flatnonzero(good);cached=np.isin(raw_positions,cache.index);qc_count=int(cached.sum());RNA_equal=0;det_equal=0
                if qc_count:
                    found=cache.loc[raw_positions[cached]];target=g.iloc[author_positions[cached]];RNA_equal=int(np.equal(found.RNA.to_numpy(),pd.to_numeric(target.nCount_RNA).to_numpy()).sum());det_equal=int(np.equal(found.det.to_numpy(),pd.to_numeric(target.nFeature_RNA).to_numpy()).sum())
                row={'branch':branch,'private_author_namespace':namespace,'private_raw_namespace':raw_namespace,'author_rows':len(g),'matching_raw_barcodes':n,'all_suffix_barcodes_match':n==len(g),'normalized_namespace_equal':norm(namespace)==norm(raw_namespace),'raw_patient_agrees_with_author':bool(h.patient_id.eq(patient).all()),'cached_native_QC_cells':qc_count,'RNA_count_exact':RNA_equal,'detected_count_exact':det_equal};candidates.append(row)
                if n==len(g):full.append({k:v for k,v in row.items() if not k.startswith('private_')})
            aggregate.append({'branch':branch,'unmatched_namespace_rows':len(g),'candidate_raw_namespaces_with_full_barcode_match':len(full),'fully_matched_candidate_metrics_no_ids':full})
    pd.DataFrame(candidates).to_csv(out/'private_all_unmatched_namespace_candidates.tsv',sep='\t',index=False);p=out/'private_terminal.json';p.write_text(json.dumps({'status':'COMPLETE_ALL_UNMATCHED_NAMESPACE_CANDIDATES_NO_ALIAS_OR_LABEL_ADOPTED','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'candidate_table_sha256':sha(out/'private_all_unmatched_namespace_candidates.tsv'),'aggregate_no_identifiers':aggregate,'new_mapping_or_effect':False},indent=2)+'\n');print(json.dumps({'status':'COMPLETE_ALL_UNMATCHED_NAMESPACE_CANDIDATES_NO_ALIAS_OR_LABEL_ADOPTED','terminal_sha256':sha(p),'aggregate_no_identifiers':aggregate}))


if __name__=='__main__':main()

"""Calibrate author QC representation on exact matches before judging orphan candidates."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from restore_GSE176021_exact_author_annotation_sidecar_v2 import original_frames,author_frame
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def summary(d):
    result={'cells':len(d)}
    for native,author,key in [('native_RNA','author_RNA','RNA'),('native_detected','author_detected','detected')]:
        x=d[native].to_numpy(float);y=d[author].to_numpy(float);assert np.isfinite(x).all() and np.isfinite(y).all() and (x>0).all() and (y>0).all();ratio=x/y
        result.update({key+'_exact_cells':int(np.equal(x,y).sum()),key+'_ratio_quantiles':np.quantile(ratio,[0,.25,.5,.75,1]).tolist(),key+'_difference_quantiles':np.quantile(x-y,[0,.25,.5,.75,1]).tolist(),key+'_Pearson':float(np.corrcoef(x,y)[0,1]) if len(d)>=3 and min(np.ptp(x),np.ptp(y))>0 else None,key+'_Spearman':float(pd.Series(x).corr(pd.Series(y),method='spearman')) if len(d)>=3 and min(np.ptp(x),np.ptp(y))>0 else None})
    return result


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2);assert sha(root/'GSE176021_author_lineage_native_RNA_audit_20261019_v1/private_terminal.json')=='6950e1c15fda41b93962cfd763348898802d8875f1eef6719ec395b25392730d';assert sha(root/'GSE176021_unmatched_author_namespace_qualification_20261019_v2/private_terminal.json')=='4d319bb64f7037b14c7bbd32be73596970f5b75446ca88f59b6d191ee007e3d8'
    out=root/'GSE176021_author_native_QC_calibration_20261019_v1';assert not out.exists();out.mkdir(mode=0o700);(out/'plan_before_QC_representation_calibration.json').write_text(json.dumps({'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'resource_guard':guard,'scope':'Every cached unique selected GSE176021 cell with exact existing author metadata; compare author nCount/nFeature to native fullRNA/detected. All unmatched namespace full-barcode candidates separately, no mapping adoption or programme association. Exact matches serve as representation controls; low equality alone does not prove false candidate.','individual_QC_arrays':'Server private only; chat/source summaries omit raw positions, barcodes and patient IDs.'},indent=2)+'\n')
    original,raw=original_frames(root);keys=pd.Index(raw.imid+'_'+raw.barcode);s=np.load(root/'GSE176021_exact_original_author_annotation_sidecar_20261019_v2/private_all_atlas_cells_exact_author_bridge.npz');q=np.load(root/'GSE176021_author_lineage_native_RNA_20261019_v1/private_selected_cell_quality_scores.npz');qc=pd.DataFrame({'global_row':q['global_row'],'native_RNA':q['fullRNA'],'native_detected':q['detected']});assert qc.groupby('global_row').native_RNA.nunique().eq(1).all() and qc.groupby('global_row').native_detected.nunique().eq(1).all();qc=qc.drop_duplicates('global_row').set_index('global_row');rr=pd.Series(s['original_raw_row'],index=s['global_atlas_row']);qc['raw_row']=rr.reindex(qc.index).to_numpy();qc=qc.set_index('raw_row');candidate=pd.read_csv(root/'GSE176021_unmatched_author_namespace_qualification_20261019_v2/private_all_unmatched_namespace_candidates.tsv',sep='\t');pairs=[];summaries=[];by_namespace=[]
    for branch in ['CD3','CD8']:
        d=author_frame(original/('GSE176021_'+branch+'_annotations.rds.gz'));loc=keys.get_indexer(d.barcode);valid=loc>=0;cached=valid&np.isin(loc,qc.index);positions=np.flatnonzero(cached);known=d.iloc[positions];found=qc.loc[loc[cached]];controls=pd.DataFrame({'branch':branch,'role':'exact_key_patient_controls','raw_row':loc[cached],'author_row':positions,'native_RNA':found.native_RNA.to_numpy(),'author_RNA':pd.to_numeric(known.nCount_RNA).to_numpy(),'native_detected':found.native_detected.to_numpy(),'author_detected':pd.to_numeric(known.nFeature_RNA).to_numpy(),'private_raw_namespace':raw.imid.iloc[loc[cached]].to_numpy()});pairs.append(controls);summaries.append({'branch':branch,'role':'exact_key_patient_controls',**summary(controls)})
        for _,g in controls.groupby('private_raw_namespace'):
            if len(g)>=20:by_namespace.append({'branch':branch,'role':'control_namespace_summary',**summary(g)})
        for c in candidate[candidate.branch.eq(branch)&candidate.all_suffix_barcodes_match].itertuples(index=False):
            unmatched=d[~valid];take=unmatched.barcode.str.rsplit('_',n=1).str[0].eq(c.private_author_namespace);unmatched=unmatched[take];bare=unmatched.barcode.str.rsplit('_',n=1).str[-1];h=raw[raw.imid.eq(c.private_raw_namespace)];ix=pd.Index(h.barcode).get_indexer(bare);assert (ix>=0).all();rawpos=raw.index.get_indexer(h.index.to_numpy())[ix];good=np.isin(rawpos,qc.index);target=unmatched[good];found=qc.loc[rawpos[good]];dt=pd.DataFrame({'branch':branch,'role':'unadopted_namespace_candidate','raw_row':rawpos[good],'author_row':target.index.to_numpy(),'native_RNA':found.native_RNA.to_numpy(),'author_RNA':pd.to_numeric(target.nCount_RNA).to_numpy(),'native_detected':found.native_detected.to_numpy(),'author_detected':pd.to_numeric(target.nFeature_RNA).to_numpy(),'private_raw_namespace':c.private_raw_namespace});pairs.append(dt);summaries.append({'branch':branch,'role':'unadopted_namespace_candidate','patient_field_agreement':False,**summary(dt)})
    allpairs=pd.concat(pairs,ignore_index=True);allpairs.to_csv(out/'private_all_cached_author_native_QC_pairs.tsv.gz',sep='\t',index=False);pd.DataFrame(by_namespace).to_json(out/'private_all_control_namespace_summary_no_identifiers.json',orient='records',indent=2);p=out/'private_terminal.json';d={'status':'COMPLETE_EXACT_MATCH_CONTROL_QC_CALIBRATION_AND_UNADOPTED_CANDIDATE_DIAGNOSTIC','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'all_QC_pair_rows':len(allpairs),'private_QC_pairs_sha256':sha(out/'private_all_cached_author_native_QC_pairs.tsv.gz'),'aggregate_no_identifiers':summaries,'no_new_alias_or_effect':True};p.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'status':d['status'],'terminal_sha256':sha(p),'aggregate_no_identifiers':summaries}))


if __name__=='__main__':main()

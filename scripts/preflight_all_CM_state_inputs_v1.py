"""All frozen172 /1290 sample metadata/row alignment, no expression copies."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd,zarr
from anndata._io.specs import read_elem

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()
def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);P=Path(os.environ['PANCANCER_CM_STATE_INPUT_PLAN']);plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256'];U=AR/'20261010_uniform_expression_annotation_v1';B=AR/'20261010_expression_annotation_downstream_v1';S=R/'immune_same_name_family_source_review_v1'
    files={'input_contract':U/'input_contract.json','original_catalogue':U/'all_node_number_catalog.tsv','current_catalogue':B/'downstream_review_v2/number_to_annotation_dictionary.tsv','counts':B/'inputs/sample_major_subcluster_counts.tsv','scope':S/'patient_block_provenance_overlay_v1/private_all1290_patient_provenance_overlay.tsv'}
    fingerprints={k:sha(v) for k,v in files.items()};scope=pd.read_csv(files['scope'],sep='\t',keep_default_na=False).set_index('sampleID');assert scope.index.is_unique and len(scope)==1290;current=pd.read_csv(files['current_catalogue'],sep='\t',keep_default_na=False).set_index('node_id');assert len(current)==172 and current.subCluster.is_unique;original=pd.read_csv(files['original_catalogue'],sep='\t',keep_default_na=False);counts=pd.read_csv(files['counts'],sep='\t',keep_default_na=False);counts=counts[counts.sampleID.isin(scope.index)];expected=counts.set_index(['sampleID','subCluster']);assert expected.index.is_unique
    contract=json.loads(files['input_contract'].read_text());assert contract['raw_layer']=='X';g=zarr.open_group(contract['core_zarr'],mode='r');sc=pd.Categorical(read_elem(g['obs']['sample_uid']));nc=pd.Categorical(read_elem(g['obs']['expression_v1_node_id']));assert len(sc)==len(nc)==contract['raw_shape'][0];si={str(v):i for i,v in enumerate(sc.categories)};ni={str(v):i for i,v in enumerate(nc.categories)};key_to_unit={};rows=[];sources=[];seen=set()
    for root in sorted((U/'by_gse').iterdir()):
        if not (root/'manifest.json').exists():continue
        sid=root.name;manifest=json.loads((root/'manifest.json').read_text());assert manifest['status']=='COMPLETE_EXACT_MEMBER_EXPRESSION_EVIDENCE'
        for n in ['units.tsv','genes.tsv']:assert sha(root/n)==manifest['output_sha256'][n]
        units=pd.read_csv(root/'units.tsv',sep='\t',keep_default_na=False);genes=pd.read_csv(root/'genes.tsv',sep='\t',keep_default_na=False);assert genes.gene.is_unique and genes.core_column.is_unique;sel=units[units.candidate_eligible & units.sample_uid.isin(scope.index)].copy();sel['node_id']=original.node_id.to_numpy()[sel.node_index];assert sel.node_id.isin(current.index).all();sel['subCluster']=current.loc[sel.node_id,'subCluster'].to_numpy();sel['major']=current.loc[sel.node_id,'majorCluster'].to_numpy()
        for row in sel.itertuples():
            pair=(row.sample_uid,row.subCluster);assert pair in expected.index and pair not in seen;seen.add(pair);assert int(row.n_cells)==int(expected.loc[pair,'n_cells']);assert row.major==expected.loc[pair,'majorCluster'];key=si[row.sample_uid]*len(nc.categories)+ni[row.node_id];assert key not in key_to_unit;key_to_unit[key]=len(rows);meta=scope.loc[row.sample_uid];rows.append(dict(unit_index=len(rows),unit_key=key,gse_id=sid,sample_uid=row.sample_uid,node_id=row.node_id,subCluster=row.subCluster,major=row.major,n_cells=int(row.n_cells),native_unit_row=int(row.unit_row),effective_resampling_block=meta.effective_resampling_block,unit_kind=meta.unit_kind,stratum=meta.stratum,analysis_eligible=meta.analysis_eligible))
        if len(sel):sources.append(dict(gse_id=sid,selected_cells=int(sel.n_cells.sum()),selected_libraries=sel.sample_uid.nunique(),selected_nodes=sel.node_id.nunique(),native_genes=len(genes),source_manifest_sha256=sha(root/'manifest.json'),native_gene_axis_sha256=sha(root/'genes.tsv')))
    assert seen==set(expected.index);E=pd.DataFrame(rows);allkeys=sc.codes.astype(np.int64)*len(nc.categories)+nc.codes;valid=(sc.codes>=0)&(nc.codes>=0);allowed=np.zeros((len(sc.categories),len(nc.categories)),bool)
    for k in key_to_unit:allowed.ravel()[k]=True
    mask=np.zeros(len(sc),bool);mask[valid]=allowed.ravel()[allkeys[valid]];positions=np.flatnonzero(mask);keys=allkeys[positions];uniq,num=np.unique(keys,return_counts=True);assert set(uniq)==set(key_to_unit)
    for k,n in zip(uniq,num):assert E.iloc[key_to_unit[int(k)]].n_cells==n
    unit=np.array([key_to_unit[int(k)] for k in keys],np.int32);assert len(positions)==E.n_cells.sum()==counts.n_cells.sum();assert E.sample_uid.nunique()==1290
    O=R/'all_CM_state_native_input_preflight_20261014_v1';assert not O.exists();O.mkdir();E.to_csv(O/'private_all_frozen_sample_node_units.tsv.gz',sep='\t',index=False,compression='gzip');np.savez_compressed(O/'private_canonical_row_selection.npz',global_rows=positions,unit_index=unit)
    pd.DataFrame(sources).to_csv(O/'all_source_native_scope.tsv',sep='\t',index=False);major=E.groupby(['gse_id','major'],observed=True,as_index=False).agg(cells=('n_cells','sum'),libraries=('sample_uid','nunique'),nodes=('node_id','nunique'));major.to_csv(O/'all_source_major_native_scope.tsv',sep='\t',index=False)
    out=dict(status='COMPLETE_ALL172_ALL1290_NATIVE_ROW_PREFLIGHT_PENDING_AUDIT',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),input_sha256=fingerprints,canonical_shape=contract['raw_shape'],sources=len(sources),libraries=1290,frozen_nodes=len(current),selected_cells=len(positions),selected_sample_node_units=len(E),all_counts_keys_and_rows_exact=True,no_expression_matrix_copied=True,original_labels_changed=False,formal_state_labels=False,source_conditions_not_independent_patients=True,private_row_selection_sha256=sha(O/'private_canonical_row_selection.npz'),private_unit_table_sha256=sha(O/'private_all_frozen_sample_node_units.tsv.gz'),public_output_sha256={p.name:sha(p) for p in O.glob('*.tsv')});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()

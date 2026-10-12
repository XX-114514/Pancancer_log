"""All frozen CM-edge/consensus LR observations and explicit missingness.

No LR rank selection, hypothesis tests, identity reassignment, or zero filling.
Private library arrays stay on the data host; public tables are context counts.
"""
from pathlib import Path
from datetime import datetime, timezone
from collections import defaultdict
import csv, gzip, hashlib, json, os
import numpy as np
import pandas as pd
from scipy.optimize import nnls
import liana as li

METRICS = ['expr_prod', 'lr_means', 'specificity_rank', 'magnitude_rank']

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda: f.read(8388608), b''):
            h.update(block)
    return h.hexdigest()

def main():
    started = datetime.now(timezone.utc).isoformat()
    R = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    AR = Path(os.environ['PANCANCER_ANNOTATION_ROOT'])
    P = R/'immune_same_name_family_source_review_v1'
    B = AR/'20261010_expression_annotation_downstream_v1'
    C = P/'CM_edge_coverage_v1'
    planpath = Path(os.environ['PANCANCER_CM_LR_PLAN'])
    plan = json.loads(planpath.read_text())
    O = P/'CM_edge_LR_extraction_v1'
    assert not O.exists() and sha(__file__) == plan['code_sha256'][Path(__file__).name]
    assert li.__version__ == plan['liana_version'] == '1.7.1'
    for key, name in [('coverage_inventory_sha256','inventory.json'), ('coverage_audit_sha256','independent_audit.json'), ('coverage_table_sha256','all60_pair_context_coverage.tsv'), ('mapping_sha256','all240_directed_CM_edge_mapping.tsv')]:
        assert sha(C/name) == plan[key], key
    old = json.loads((C/'inventory.json').read_text())
    U = P/'patient_block_provenance_overlay_v1'
    paths = dict(counts=B/'inputs/sample_major_subcluster_counts.tsv', sample_metadata=B/'inputs/sample_metadata.tsv', W=B/'covarnet/tables/nmf_W_subcluster_by_cm.tsv', H=B/'covarnet/tables/nmf_H_cm_by_sample.tsv', F=B/'covarnet/tables/frequency_normalized.tsv', catalogue=B/'downstream_review_v2/number_to_annotation_dictionary.tsv', liana_status=B/'liana/sample_status.tsv', liana_audit=B/'liana/independent_audit.json', provenance_scope=U/'private_all1290_patient_provenance_overlay.tsv')
    for k,p in paths.items():
        assert sha(p) == old['input_sha256'][k]['sha256'], k
    resource = li.resource.select_resource('consensus')[['ligand','receptor']].sort_values(['ligand','receptor']).reset_index(drop=True)
    resource_text = resource.to_csv(sep='\t', index=False)
    assert hashlib.sha256(resource_text.encode()).hexdigest() == plan['resource_sha256']
    assert len(resource) == 4624 and not resource.duplicated().any()
    LR = list(resource.itertuples(index=False, name=None))
    lrindex = {key:i for i,key in enumerate(LR)}
    scope = pd.read_csv(paths['provenance_scope'], sep='\t', keep_default_na=False).set_index('sampleID')
    status = pd.read_csv(paths['liana_status'], sep='\t', keep_default_na=False).set_index('sampleID')
    catalogue = pd.read_csv(paths['catalogue'], sep='\t', keep_default_na=False).set_index('subCluster')
    counts = pd.read_csv(paths['counts'], sep='\t')
    counts = counts[counts.sampleID.isin(scope.index)]
    nodes = counts.pivot_table(index='sampleID', columns='subCluster', values='n_cells', aggfunc='sum', fill_value=0).reindex(scope.index, fill_value=0)
    majors = counts.groupby(['sampleID','majorCluster']).n_cells.sum().unstack(fill_value=0).reindex(scope.index, fill_value=0)
    W = pd.read_csv(paths['W'], sep='\t', index_col=0)
    H = pd.read_csv(paths['H'], sep='\t', index_col=0)
    F = pd.read_csv(paths['F'], sep='\t', index_col=0)
    assert W.shape == (172,9) and H.shape == (9,1290) and F.shape == (172,1290)
    assert list(W.index) == list(F.index) and list(W.columns) == list(H.index) and list(H.columns) == list(F.columns)
    assert set(H.columns) == set(scope.index) and scope.index.is_unique
    scope = scope.loc[H.columns]
    coverage = pd.read_csv(C/'all60_pair_context_coverage.tsv', sep='\t')
    q = coverage[(coverage.scope_mode=='reviewed_context') & coverage.coverage_ge20libraries_ge10blocks].sort_values(['stratum','pair_id']).reset_index(drop=True)
    assert len(q)==34 and q.stratum.nunique()==10 and q.gse_id.nunique()==8
    mapping = pd.read_csv(C/'all240_directed_CM_edge_mapping.tsv', sep='\t')
    native = {}
    for sid, expected in plan['native_gene_inputs'].items():
        path = AR/'20261010_uniform_expression_annotation_v1/by_gse'/sid/'genes.tsv'
        assert sha(path)==expected['sha256']
        g = pd.read_csv(path, sep='\t', keep_default_na=False).gene
        assert len(g)==expected['genes'] and g.is_unique
        native[sid]=set(g)
    O.mkdir(); (O/'private_group_arrays').mkdir()
    groups=[]; queries=defaultdict(list)
    for ix,row in q.iterrows():
        mask = scope.stratum.eq(row.stratum) & scope.analysis_eligible & nodes.loc[scope.index,row.source].ge(20) & nodes.loc[scope.index,row.target].ge(20) & status.loc[scope.index].status.eq('success')
        samples = scope.index[mask].to_numpy(dtype=str)
        assert len(samples)==row.n_endpoint_eligible_LIANA_success
        meta = scope.loc[samples]
        assert meta.effective_resampling_block.nunique()==row.n_analysis_blocks
        pair=(row.source,row.target)
        keep=~W.index.isin(pair); basis=W.loc[keep].to_numpy(float)
        rank=int(np.linalg.matrix_rank(basis)); cond=float(np.linalg.cond(basis))
        activity=np.full((len(samples),9),np.nan)
        activity_status='ENDPOINT_EXCLUDED_BASIS_NOT_IDENTIFIABLE'
        if rank==9 and cond<1e8:
            try:
                activity=np.stack([nnls(basis,F.loc[keep,s].to_numpy(float),maxiter=5000)[0] for s in samples])
                activity_status='NNLS_POINT_ACTIVITY_AVAILABLE_NOT_CM_DRIVER_EVIDENCE'
            except RuntimeError:
                activity_status='NNLS_MAXITER_REACHED_NO_ACTIVITY'
                activity[:]=np.nan
        for direction,(s,t) in enumerate([pair,pair[::-1]]):
            gid='GRP'+str(2*ix+direction+1).zfill(3)
            measured=np.array([set(l.split('_')+r.split('_')).issubset(native[row.gse_id]) for l,r in LR])
            a=dict(group_id=gid,pair_id=row.pair_id,stratum=row.stratum,gse_id=row.gse_id,cancerType=row.cancerType,sample_class=row.sample_class,source=s,target=t,samples=samples,meta=meta,measured=measured,values=np.full((len(samples),len(LR),4),np.nan),activity=activity,activity_status=activity_status,basis_rank=rank,basis_condition=cond,CMs=';'.join(sorted(mapping[(mapping.pair_id==row.pair_id)&(mapping.source==s)&(mapping.target==t)].cm.unique())))
            assert a['CMs']
            groups.append(a)
            for k,sample in enumerate(samples): queries[sample].append((len(groups)-1,k))
    assert len(groups)==68
    audit_records={r['sampleID']:r for r in json.loads(paths['liana_audit'].read_text())['sample_reports']}
    libraryproof=[]; scanned=0; extracted=0
    private = O/'private_all_observed_LR_rows.tsv.gz'
    with gzip.open(private,'wt') as target:
        writer=csv.DictWriter(target,fieldnames=['group_id','sampleID','gse_id','source','target','ligand_complex','receptor_complex']+METRICS,delimiter='\t'); writer.writeheader()
        for sample in sorted(queries):
            st=status.loc[sample]; expected=audit_records[sample]
            assert st.status==expected['status']=='success' and st.input_signature==expected['input_signature'] and st.output_sha256==expected['output_sha256']
            path=Path(st.output_path); assert path.parent==B/'liana/per_sample' and path.name==hashlib.sha256(sample.encode()).hexdigest()[:20]+'.tsv.gz'
            assert sha(path)==st.output_sha256
            assert int(st.n_native_measured_genes)==len(native[st.gse_id])
            df=pd.read_csv(path,sep='\t',keep_default_na=False); assert len(df)==int(st.n_interactions)
            scanned+=len(df)
            assert df.sampleID.eq(sample).all() and df.gse_id.eq(st.gse_id).all()
            for gi,k in queries[sample]:
                group=groups[gi]; sub=df[df.source.eq(group['source'])&df.target.eq(group['target'])]
                assert not sub.duplicated(['ligand_complex','receptor_complex']).any()
                for r in sub.itertuples(index=False):
                    key=(r.ligand_complex,r.receptor_complex); assert key in lrindex,key
                    j=lrindex[key]; assert group['measured'][j]
                    v=np.array([getattr(r,m) for m in METRICS],float)
                    assert np.isfinite(v).all() and (v>=0).all() and (v[2:]<=1).all()
                    assert np.isnan(group['values'][k,j]).all()
                    group['values'][k,j]=v
                    writer.writerow(dict(group_id=group['group_id'],sampleID=sample,gse_id=st.gse_id,source=group['source'],target=group['target'],ligand_complex=key[0],receptor_complex=key[1],**dict(zip(METRICS,v))))
                    extracted+=1
            libraryproof.append(dict(sampleID=sample,gse_id=st.gse_id,input_signature=st.input_signature,output_sha256=st.output_sha256,rows=len(df)))
            if len(libraryproof)%50==0: print('scanned_libraries',len(libraryproof),'rows',scanned,'extracted',extracted,flush=True)
    (O/'private_library_input_proofs.json').write_text(json.dumps(libraryproof,indent=2)+'\n')
    casepath=O/'all314432_consensus_case_coverage.tsv.gz'; summary=[]; cases=0
    with gzip.open(casepath,'wt') as file:
        writer=None
        for group in groups:
            samples=group['samples']; meta=group['meta']; values=group['values']; n=len(samples); observed=np.isfinite(values[:,:,0]); assert np.array_equal(np.isfinite(values),np.repeat(observed[:,:,None],4,axis=2))
            sm=catalogue.loc[group['source'],'majorCluster']; tm=catalogue.loc[group['target'],'majorCluster']
            np.savez_compressed(O/'private_group_arrays'/(group['group_id']+'.npz'),samples=samples,patient_blocks=meta.effective_resampling_block.to_numpy(dtype=str),reviewed_patient=meta.unit_kind.eq('reviewed_patient').to_numpy(),ligands=resource.ligand.to_numpy(dtype=str),receptors=resource.receptor.to_numpy(dtype=str),native_measured=group['measured'],scores=values,source_n=nodes.loc[samples,group['source']].to_numpy(),target_n=nodes.loc[samples,group['target']].to_numpy(),source_within_major=(nodes.loc[samples,group['source']]/majors.loc[samples,sm]).to_numpy(),target_within_major=(nodes.loc[samples,group['target']]/majors.loc[samples,tm]).to_numpy(),endpoint_excluded_activity=group['activity'],original_H=H.loc[:,samples].to_numpy().T,CMs=W.columns.to_numpy(dtype=str))
            gate=0
            for j,(ligand,receptor) in enumerate(LR):
                obs=observed[:,j]; good=meta.iloc[np.flatnonzero(obs)]; no=int(obs.sum()); nb=good.effective_resampling_block.nunique(); reviewed=good[good.unit_kind.eq('reviewed_patient')]; measured=bool(group['measured'][j]); missing=sorted(set(ligand.split('_')+receptor.split('_'))-native[group['gse_id']]); data_gate=measured and no>=20 and nb>=10 and no/n>=.8; gate+=int(data_gate)
                row={k:group[k] for k in ['group_id','pair_id','stratum','gse_id','cancerType','sample_class','source','target','CMs']}
                row.update(case_id=group['group_id']+'_LR'+str(j+1).zfill(4),ligand=ligand,receptor=receptor,n_endpoint_libraries=n,n_endpoint_blocks=meta.effective_resampling_block.nunique(),native_subunits_measured=measured,missing_native_subunits=';'.join(missing),n_observed_libraries=no,n_observed_blocks=nb,n_observed_reviewed_libraries=len(reviewed),n_observed_reviewed_source_patient_blocks=reviewed.effective_resampling_block.nunique(),n_observed_sample_fallback_libraries=no-len(reviewed),n_native_missing_libraries=n if not measured else 0,n_measured_no_LIANA_output=n-no if measured else 0,observed_coverage=no/n,data_coverage_gate_only=data_gate,expr_prod_variable=bool(np.unique(values[obs,j,0]).size>1),activity_status=group['activity_status'],basis_rank=group['basis_rank'],basis_condition=group['basis_condition'],source_identity_status=catalogue.loc[group['source'],'identity_status'],target_identity_status=catalogue.loc[group['target'],'identity_status'],binding_status='symbol_warning_ACTR2_not_ACVR2A' if (ligand,receptor)==('INHBA','ACTR2') else 'disputed_descriptive_C1QA_CD93' if (ligand,receptor)==('C1QA','CD93') else 'resource_entry_not_mechanistically_reviewed')
                for ixmetric,m in enumerate(METRICS): row['mean_observed_'+m]=float(values[obs,j,ixmetric].mean()) if no else None
                if writer is None: writer=csv.DictWriter(file,fieldnames=list(row),delimiter='\t');writer.writeheader()
                writer.writerow(row);cases+=1
            summary.append({k:group[k] for k in ['group_id','pair_id','stratum','gse_id','cancerType','source','target','CMs','activity_status','basis_rank','basis_condition']}|dict(n_endpoint_libraries=n,n_endpoint_blocks=meta.effective_resampling_block.nunique(),consensus_cases=4624,observed_library_LR_records=int(observed.sum()),native_measured_cases=int(group['measured'].sum()),data_coverage_gate_cases=gate))
    assert cases==314432 and extracted==sum(r['observed_library_LR_records'] for r in summary)
    pd.DataFrame(summary).to_csv(O/'all68_directed_group_summary.tsv',sep='\t',index=False)
    inv=dict(schema='original_CM_complete_consensus_LR_extraction_v1',status='COMPLETE_EXTRACTION_PENDING_INDEPENDENT_AUDIT',started_UTC=started,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planpath),coverage_inventory_sha256=plan['coverage_inventory_sha256'],resource_sha256=plan['resource_sha256'],native_gene_inputs=plan['native_gene_inputs'],selected_unique_libraries=len(queries),unique_GSE=8,reviewed_contexts=10,pair_contexts=34,directed_group_contexts=68,consensus_resource_rows=4624,case_coverage_rows=cases,full_LIANA_rows_scanned=scanned,observed_library_LR_records=extracted,data_coverage_gate_cases=sum(r['data_coverage_gate_cases'] for r in summary),scores_selected_by_rank_or_p=False,unemitted_interactions_are_zero=False,new_associations=False,independent_biological_validation=False,outputs={str(p.relative_to(O)):sha(p) for p in O.rglob('*') if p.is_file()})
    (O/'inventory.json').write_text(json.dumps(inv,indent=2)+'\n')
    print(json.dumps({k:v for k,v in inv.items() if k!='outputs'}),flush=True)
if __name__=='__main__':main()

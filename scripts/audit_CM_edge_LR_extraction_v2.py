"""Independent raw-CSV score replay and all consensus missingness audit.

Uses row counters, not producer pivot/filter logic. Checks every selected raw
score and each public case, plus NNLS KKT conditions for remaining-CM activity.
"""
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter,defaultdict
import csv,gzip,json,hashlib,os
import numpy as np,pandas as pd
from scipy.linalg import qr

METRICS=['expr_prod','lr_means','specificity_rank','magnitude_rank']
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()
def read(p):return list(csv.DictReader(Path(p).open(),delimiter='\t'))
def main():
    start=datetime.now(timezone.utc).isoformat();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);P=R/'immune_same_name_family_source_review_v1';B=AR/'20261010_expression_annotation_downstream_v1';C=P/'CM_edge_coverage_v1';O=P/'CM_edge_LR_extraction_v2';planpath=Path(os.environ['PANCANCER_CM_LR_PLAN']);plan=json.loads(planpath.read_text());assert sha(__file__)==plan['code_sha256'][Path(__file__).name] and not (O/'independent_audit.json').exists();inv=json.loads((O/'inventory.json').read_text());assert inv['plan_sha256']==sha(planpath)
    for n,h in inv['outputs'].items():assert sha(O/n)==h
    coverage_inv=json.loads((C/'inventory.json').read_text());assert sha(C/'inventory.json')==plan['coverage_inventory_sha256'] and sha(C/'independent_audit.json')==plan['coverage_audit_sha256']
    U=P/'patient_block_provenance_overlay_v1'
    paths=dict(counts=B/'inputs/sample_major_subcluster_counts.tsv',sample_metadata=B/'inputs/sample_metadata.tsv',W=B/'covarnet/tables/nmf_W_subcluster_by_cm.tsv',H=B/'covarnet/tables/nmf_H_cm_by_sample.tsv',F=B/'covarnet/tables/frequency_normalized.tsv',catalogue=B/'downstream_review_v2/number_to_annotation_dictionary.tsv',liana_status=B/'liana/sample_status.tsv',liana_audit=B/'liana/independent_audit.json',provenance_scope=U/'private_all1290_patient_provenance_overlay.tsv')
    for k,p in paths.items():assert sha(p)==coverage_inv['input_sha256'][k]['sha256'],k
    resourcepath=planpath.parent/'frozen_LIANA171_consensus4624_v1.tsv';assert sha(resourcepath)==plan['resource_sha256'];LR=[(r['ligand'],r['receptor']) for r in read(resourcepath)];assert len(LR)==4624 and len(set(LR))==4624;lrindex={k:i for i,k in enumerate(LR)}
    scope={r['sampleID']:r for r in read(paths['provenance_scope'])};status={r['sampleID']:r for r in read(paths['liana_status'])};catalogue={r['subCluster']:r for r in read(paths['catalogue'])};count=Counter();major=Counter()
    for r in read(paths['counts']):
        if r['sampleID'] in scope:count[r['sampleID'],r['subCluster']]+=int(r['n_cells']);major[r['sampleID'],r['majorCluster']]+=int(r['n_cells'])
    cov=read(C/'all60_pair_context_coverage.tsv');q=sorted([r for r in cov if r['scope_mode']=='reviewed_context' and r['coverage_ge20libraries_ge10blocks']=='True'],key=lambda r:(r['stratum'],r['pair_id']));assert len(q)==34
    native={}
    for sid,expected in plan['native_gene_inputs'].items():
        p=AR/'20261010_uniform_expression_annotation_v1/by_gse'/sid/'genes.tsv';assert sha(p)==expected['sha256'];genes=[r['gene'] for r in read(p)];assert len(genes)==len(set(genes))==expected['genes'];native[sid]=set(genes)
    W=pd.read_csv(paths['W'],sep='\t',index_col=0);F=pd.read_csv(paths['F'],sep='\t',index_col=0);H=pd.read_csv(paths['H'],sep='\t',index_col=0);summ={r['group_id']:r for r in read(O/'all68_directed_group_summary.tsv')};maps=read(C/'all240_directed_CM_edge_mapping.tsv');groups=[];queries=defaultdict(list);kkt=0
    for ix,c in enumerate(q):
        eligible={s for s,r in scope.items() if r['stratum']==c['stratum'] and r['analysis_eligible']=='True' and count[s,c['source']]>=20 and count[s,c['target']]>=20 and status[s]['status']=='success'}
        assert len(eligible)==int(c['n_endpoint_eligible_LIANA_success'])
        for direction,(source,target) in enumerate([(c['source'],c['target']),(c['target'],c['source'])]):
            gid='GRP'+str(2*ix+direction+1).zfill(3);saved=dict(np.load(O/'private_group_arrays'/(gid+'.npz'),allow_pickle=False));samples=list(saved['samples']);assert len(samples)==len(set(samples)) and set(samples)==eligible;blocks=[scope[s]['effective_resampling_block'] for s in samples];reviewed=np.array([scope[s]['unit_kind']=='reviewed_patient' for s in samples]);assert list(saved['patient_blocks'])==blocks and np.array_equal(saved['reviewed_patient'],reviewed);assert list(zip(saved['ligands'],saved['receptors']))==LR
            for label,node in [('source',source),('target',target)]:
                n=np.array([count[s,node] for s in samples]);den=np.array([major[s,catalogue[node]['majorCluster']] for s in samples]);assert (den>=n).all() and np.array_equal(saved[label+'_n'],n) and np.allclose(saved[label+'_within_major'],n/den,rtol=1e-12,atol=1e-12)
            measured=np.array([all(g in native[c['gse_id']] for g in l.split('_')+r.split('_')) for l,r in LR]);assert np.array_equal(saved['native_measured'],measured)
            mask=~W.index.isin([source,target]);A=W.loc[mask].to_numpy(float);_,rr,_=qr(A,mode='economic',pivoting=True);rank=int(np.sum(np.abs(np.diag(rr))>np.finfo(float).eps*max(A.shape)*abs(rr[0,0])));condition=np.linalg.cond(A);assert int(summ[gid]['basis_rank'])==rank and np.isclose(float(summ[gid]['basis_condition']),condition)
            activity=saved['endpoint_excluded_activity'];assert activity.shape==(len(samples),9) and np.array_equal(saved['original_H'],H.loc[:,samples].to_numpy().T) and list(saved['CMs'])==list(W.columns)
            if summ[gid]['activity_status']=='NNLS_POINT_ACTIVITY_AVAILABLE_NOT_CM_DRIVER_EVIDENCE':
                assert rank==9 and condition<1e8 and np.isfinite(activity).all() and (activity>=0).all()
                for i,sample in enumerate(samples):
                    residual=A@activity[i]-F.loc[mask,sample].to_numpy(float);gradient=A.T@residual;scale=max(1.,np.linalg.norm(A.T@F.loc[mask,sample].to_numpy(float),np.inf));assert (gradient>=-1e-8*scale).all();assert np.max(np.abs(activity[i]*gradient))<1e-7*scale;kkt+=1
            else:assert np.isnan(activity).all()
            CMs=';'.join(sorted({r['cm'] for r in maps if r['pair_id']==c['pair_id'] and r['source']==source and r['target']==target}));assert CMs==summ[gid]['CMs']
            group=dict(gid=gid,coverage=c,source=source,target=target,samples=samples,blocks=blocks,reviewed=reviewed,measured=measured,saved=saved,expected=np.full_like(saved['scores'],np.nan),CMs=CMs);groups.append(group)
            for i,sample in enumerate(samples):queries[sample].append((len(groups)-1,i))
    parent=json.loads(paths['liana_audit'].read_text());assert sha(B/'liana/manifest.json')==parent['source_manifest_sha256'] and sha(B/'run_contract.json')==parent['input_contract_sha256']
    proofs=json.loads((O/'private_library_input_proofs.json').read_text());assert {r['sampleID'] for r in proofs}==set(queries) and len(proofs)==len(queries);proofs={r['sampleID']:r for r in proofs};accepted={r['sampleID']:r for r in json.loads(paths['liana_audit'].read_text())['sample_reports']};rawrows=0;observed=0
    for sample in sorted(queries):
        statusfile=B/'liana/status'/(hashlib.sha256(sample.encode()).hexdigest()[:20]+'.json');sidecar=json.loads(statusfile.read_text());assert sha(statusfile)==proofs[sample]['status_sidecar_sha256'];assert sidecar['sampleID']==sample and sidecar['gse_id']==status[sample]['gse_id'] and sidecar['input_signature']==status[sample]['input_signature'] and sidecar['output_sha256']==status[sample]['output_sha256'] and sidecar['status']=='success'
        s=status[sample];assert s['output_sha256']==accepted[sample]['output_sha256']==proofs[sample]['output_sha256'] and s['input_signature']==proofs[sample]['input_signature'] and accepted[sample]['signature_independently_reconstructed'] is True;path=Path(s['output_path']);assert sha(path)==s['output_sha256'];local={(groups[k]['source'],groups[k]['target']):(k,i) for k,i in queries[sample]};assert len(local)==len(queries[sample]);seen=set();nr=0
        with gzip.open(path,'rt') as f:
            for row in csv.DictReader(f,delimiter='\t'):
                nr+=1;assert row['sampleID']==sample and row['gse_id']==s['gse_id'];pair=(row['source'],row['target'])
                if pair not in local:continue
                key=(row['ligand_complex'],row['receptor_complex']);assert key in lrindex;j=lrindex[key];k,i=local[pair];assert groups[k]['measured'][j] and (k,j) not in seen;seen.add((k,j));groups[k]['expected'][i,j]=[float(row[m]) for m in METRICS];observed+=1
        assert nr==int(s['n_interactions'])==proofs[sample]['rows'];rawrows+=nr
    assert rawrows==inv['full_LIANA_rows_scanned'] and observed==inv['observed_library_LR_records']
    scores_checked=0
    for g in groups:
        assert np.allclose(g['expected'],g['saved']['scores'],rtol=1e-10,atol=1e-12,equal_nan=True);scores_checked+=g['expected'].size
    private_seen=set()
    with gzip.open(O/'private_all_observed_LR_rows.tsv.gz','rt') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            k=int(row['group_id'][3:])-1;g=groups[k];i=g['samples'].index(row['sampleID']);j=lrindex[row['ligand_complex'],row['receptor_complex']];key=(k,i,j);assert key not in private_seen;private_seen.add(key);assert row['gse_id']==g['coverage']['gse_id'] and row['source']==g['source'] and row['target']==g['target'];assert np.allclose([float(row[m]) for m in METRICS],g['expected'][i,j],rtol=1e-12,atol=1e-12)
    assert len(private_seen)==observed
    cases=0;gates=Counter();obs_per_group=Counter()
    with gzip.open(O/'all314432_consensus_case_coverage.tsv.gz','rt') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            k=cases//4624;j=cases%4624;g=groups[k];ligand,receptor=LR[j];n=len(g['samples']);ob=np.isfinite(g['expected'][:,j,0]);indices=np.flatnonzero(ob);no=len(indices);blocks={g['blocks'][i] for i in indices};review=[i for i in indices if g['reviewed'][i]];measured=bool(g['measured'][j]);gate=measured and no>=20 and len(blocks)>=10 and no/n>=.8
            assert row['group_id']==g['gid'] and row['case_id']==g['gid']+'_LR'+str(j+1).zfill(4) and row['ligand']==ligand and row['receptor']==receptor
            for field in ['pair_id','stratum','gse_id','cancerType','sample_class']:assert row[field]==g['coverage'][field]
            assert row['source']==g['source'] and row['target']==g['target'] and row['CMs']==g['CMs'];assert row['missing_native_subunits']==';'.join(sorted(set(ligand.split('_')+receptor.split('_'))-native[g['coverage']['gse_id']]))
            expected=dict(n_endpoint_libraries=n,n_endpoint_blocks=len(set(g['blocks'])),n_observed_libraries=no,n_observed_blocks=len(blocks),n_observed_reviewed_libraries=len(review),n_observed_reviewed_source_patient_blocks=len({g['blocks'][i] for i in review}),n_observed_sample_fallback_libraries=no-len(review),n_native_missing_libraries=n if not measured else 0,n_measured_no_LIANA_output=n-no if measured else 0)
            for field,value in expected.items():assert int(row[field])==value
            for field,value in dict(native_subunits_measured=measured,data_coverage_gate_only=gate,expr_prod_variable=np.unique(g['expected'][ob,j,0]).size>1).items():assert row[field]==str(bool(value))
            assert np.isclose(float(row['observed_coverage']),no/n)
            for ix,m in enumerate(METRICS):
                if no:assert np.isclose(float(row['mean_observed_'+m]),g['expected'][ob,j,ix].mean(),rtol=1e-10,atol=1e-12)
                else:assert row['mean_observed_'+m]==''
            for label,node in [('source',g['source']),('target',g['target'])]:assert row[label+'_identity_status']==catalogue[node]['identity_status']
            binding='symbol_warning_ACTR2_not_ACVR2A' if (ligand,receptor)==('INHBA','ACTR2') else 'disputed_descriptive_C1QA_CD93' if (ligand,receptor)==('C1QA','CD93') else 'resource_entry_not_mechanistically_reviewed';assert row['binding_status']==binding
            gates[g['gid']]+=int(gate);obs_per_group[g['gid']]+=no;cases+=1
    assert cases==314432 and sum(gates.values())==inv['data_coverage_gate_cases']
    for gid,s in summ.items():assert int(s['data_coverage_gate_cases'])==gates[gid] and int(s['observed_library_LR_records'])==obs_per_group[gid]
    audit=dict(status='PASS_ALL314432_CASES_FULL_SELECTED_RAW_LIANA_SCORES_NATIVE_MISSINGNESS_BLOCKS_AND_NNLS_KKT',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(O/'inventory.json'),unique_libraries=len(queries),raw_LIANA_rows_replayed=rawrows,observed_library_LR_records_replayed=observed,private_score_values_including_missing_checked=scores_checked,public_cases_replayed=cases,group_activity_KKT_checks=kkt,primary_association_tests=False,independent_biological_validation=False)
    (O/'independent_audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit),flush=True)
if __name__=='__main__':main()

"""All measured proxies against all fixed functions on disjoint reference RNA."""
import hashlib,json,os,shutil
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from anndata._io.specs import read_elem
from state_proxy_disjoint_function_math_v1 import fit
from qualify_native_spatial_state_probe_contract_v1 import resolved,physical_keys


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def read_axis(path):
    with h5py.File(path,'r') as h:
        k=h['var'].attrs['_index'];k=k.decode() if isinstance(k,bytes) else k
        return np.asarray(read_elem(h['var'][k]),str)


def normalize(g):return str(g).strip().upper()


def make_contract(root,V):
    p=Path('/dev/shm/pancancer-native-state-probe-match-20261020-v1/private_output/private_frozen_matched_measurement_and_target_contract.json');assert sha(p)=='0d6f8e2bac99780b7a7b211c33468ae243ece7920167c3ce85fcf465d994345e'
    nativep=V/'state_spatial_proxy_native_cells_20261016_v1/inventory.json';assert sha(nativep)=='91be1ae48a338aa6bfcd2d89107fd62bf3a00566a90ec5aa953fe0a3f37478f0'
    parent=json.loads(p.read_text());used={m['proxy_group_id'] for t in parent['targets'] for m in t['proxy_mappings']};groups=[g for g in parent['groups'] if g['proxy_group_id'] in used];assert len(groups)==31
    dictionaryp=V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json';assert sha(dictionaryp)=='6a3c4311d89a9a02d55a19f2d249c960c508b5973368c06c2e5c43bb605edec9'
    dictionary={r['programme_id']:r for r in json.loads(dictionaryp.read_text())['programmes']}
    panelp=V/'state_spatial_proxy_calibration_inputs_20261016_v1/fixed_all14_panel_gene_axes.json';panels={}
    for panel in json.loads(panelp.read_text()):panels.setdefault(panel['gene_axis_sha256'],panel)
    ni=json.loads(nativep.read_text());lookup={(r['gse_id'],r['major']):r for r in ni['private_native_files']}
    pairs=sorted({(s,g['major']) for g in groups for s in g['reference_GSEs']});files=[];unique_axes=[]
    for i,key in enumerate(pairs):
        info=lookup[key];path=V/'state_spatial_proxy_native_cells_20261016_v1'/info['relative_file'];assert sha(path)==info['sha256']
        genes=read_axis(path);upper=[normalize(g) for g in genes];counts=Counter(upper);unique_axes.append({g for g in upper if g and counts[g]==1})
        files.append(dict(native_ordinal=i,gse_id=key[0],major=key[1],path=str(path),sha256=info['sha256'],cells=info['cells'],relative_file=info['relative_file']))
    common=sorted(set.intersection(*unique_axes));assert len(files)==12 and sum(int(x['cells']) for x in files)==4349
    gene_path=root/'GSE243013_genes.csv.gz';assert sha(gene_path)=='ac325c9ce6bfaf2270eb931e7d947bf31d1fb151b14cb166ed9f6c8f4b34a104'
    native_symbols=pd.read_csv(gene_path,keep_default_na=False).geneSymbol.map(normalize).tolist();lrpath=root/'NSCLC_communication_tensor_family_20261019_v1/private_tensor_axes.npz'
    with np.load(lrpath,allow_pickle=False) as z:LR={native_symbols[int(j)] for j in z['native_gene_columns']}
    assert len(LR)==1840
    referencep=root/'NSCLC92_receiver_function_reference_scope_20261019_v1/private_all50_LR_disjoint_reference_definitions.json';assert sha(referencep)=='89ad355c5c6f02f0d0e68b2b0f35cd7e1cc217097512c99e2e306e5f9c92f4fa'
    refs=json.loads(referencep.read_text());assert len(refs)==50
    definitions=[];parents=[]
    for g in groups:
        blocked=set(LR);panel=panels[g['axis_id']];pk=physical_keys(panel)
        for gene,j in panel['symbol_to_feature_index'].items():
            if pk[j] in g['signature_probe_keys']:blocked.add(normalize(gene))
        for pid in g['parent_programme_ids']:
            rec=dictionary[pid];blocked.update(map(normalize,rec['primary_genes']))
            parents.append(dict(proxy_group_id=g['proxy_group_id'],major=g['major'],parent_programme_id=pid,primary_gene_count=len(rec['primary_genes']),
              algorithm_stable_candidate=rec['algorithm_stable_candidate'],state_naming_status=rec['state_naming_status'],best_reference_descriptor=rec['best_reference_descriptor'],best_reference_jaccard=rec['best_reference_jaccard'],
              foreign_RNA_sensitive=rec['foreign_RNA_sensitive'],basic_expression_review_flag=rec['basic_expression_review_flag'],formal_state_labels=rec['formal_state_labels']))
        bg=[s for s in common if s not in blocked];assert len(bg)>1500;at={s:i for i,s in enumerate(bg)};functions=[]
        for ref in refs:
            full=set(map(normalize,ref['full_reference_genes']));remaining=full-blocked;measured=sorted(remaining&set(bg));coverage=len(measured)/len(remaining) if remaining else 0.;retention=len(remaining)/len(full)
            eligible=len(measured)>=20 and coverage>=.8
            functions.append(dict(reference=ref['reference'],indices=[at[s] for s in measured],original_genes=len(full),remaining_genes=len(remaining),mapped_genes=len(measured),remaining_coverage=coverage,retained_fraction=retention,measurement_eligible=eligible,strong_reference_representation=bool(eligible and retention>=.8)))
        definitions.append(dict(proxy=g,background_symbols=bg,blocked_symbols=sorted(blocked),functions=functions))
    return dict(created_UTC=datetime.now(timezone.utc).isoformat(),all_proxy_groups=31,all_references=50,source_proxy_positions=203,source_function_positions=10150,
                common_unique_native_symbols=common,native_inputs=files,definitions=definitions,parent_descriptors=parents,LR_universe_symbols=sorted(LR),
                source_SHA=dict(physical_contract=sha(p),dictionary=sha(dictionaryp),native_inventory=sha(nativep),panel_definition=sha(panelp),reference_definition=sha(referencep),LR_axes=sha(lrpath)),
                scoring='Fixed common native background minus all1840 LR/native symbols, all fullparent primary genes and partialsignature physicalaliases. Same background/function genes in every reference source per proxy; no constants removed or zero-imputed.')


def ucells(x,bg,signatures,independent=False):
    out=np.full((x.shape[0],len(signatures)),np.nan);cap=min(1500,len(bg));positions=0
    for start in range(0,x.shape[0],128):
        end=min(start+128,x.shape[0]);a=x[start:end][:,bg].toarray().astype(float)
        assert np.isfinite(a).all() and (a>=0).all()
        if not independent:
            ranks=np.minimum(cap,rankdata(-a,axis=1,method='average'));ranks[a<=0]=cap
        else:
            ranks=np.full(a.shape,float(cap))
            for i,row in enumerate(a):
                values,inv,counts=np.unique(row,return_inverse=True,return_counts=True);greater=len(row)-np.cumsum(counts)
                rr=np.minimum(cap,1+greater+(counts-1)/2);rr[values<=0]=cap;ranks[i]=rr[inv]
        positions+=a.size
        for j,sig in enumerate(signatures):
            if sig is None:continue
            k=len(sig);assert 0<k<cap;b=k*(k+1)/2
            out[start:end,j]=1-(ranks[:,sig].sum(1)-b)/(k*cap-b)
    return out,positions


def load_source_and_scores(info,definition,V,panels,score_files,independent=False):
    assert sha(info['path'])==info['sha256'];a=ad.read_h5ad(info['path']);assert a.n_obs==info['cells']
    obs=a.obs;assert np.isfinite(obs[['full_RNA_total','full_detected_genes']].to_numpy(float)).all()
    np.testing.assert_allclose(np.asarray(a.X.sum(1)).ravel(),obs.full_RNA_total.to_numpy(),atol=0,rtol=0)
    np.testing.assert_array_equal(np.asarray((a.X>0).sum(1)).ravel(),obs.full_detected_genes.to_numpy())
    genes=np.asarray(a.var_names,str);upper=[normalize(s) for s in genes];cnt=Counter(upper);lookup={s:j for j,s in enumerate(upper) if cnt[s]==1};literal={s:j for j,s in enumerate(genes)}
    g=definition['proxy'];panel=panels[g['axis_id']];pk=physical_keys(panel);names=resolved(panel,set(genes));keygene={pk[panel['symbol_to_feature_index'][s]]:s for s in names}
    assert set(g['background_probe_keys'])==set(keygene)
    cols=[literal[keygene[k]] for k in g['background_probe_keys']];sig=[g['background_probe_keys'].index(k) for k in g['signature_probe_keys']]
    p,partial_positions=ucells(a.X,cols,[sig],independent);p=p[:,0]
    rel='private_cell_scores/'+Path(info['relative_file']).stem+'_'+g['axis_id'][:10]+'.npz';entry=score_files[rel];cache=V/'state_spatial_partial_panel_calibration_20261016_v2'/rel;assert sha(cache)==entry['sha256']
    with np.load(cache,allow_pickle=False) as z:
        np.testing.assert_array_equal(z['global_rows'],obs.core_position);np.testing.assert_array_equal(z['unit_index'],obs.unit_index)
        for pid in g['parent_programme_ids']:
            j=z['programme_ids'].tolist().index(pid);np.testing.assert_allclose(p,z['proxy'][:,j].astype(float),atol=1e-6,rtol=0)
    background=[lookup[s] for s in definition['background_symbols']];fs=definition['functions'];signatures=[f['indices'] if f['measurement_eligible'] else None for f in fs]
    function,positions=ucells(a.X,background,signatures,independent)
    xx=a.X[:,background];maxima=np.asarray(xx.max(0).toarray()).ravel();minima=np.asarray(xx.min(0).toarray()).ravel();variable=maxima>minima
    variable_counts=[int(variable[f['indices']].sum()) for f in fs]
    quality=np.column_stack([obs[['full_RNA_total','full_detected_genes']].to_numpy(float),np.asarray(a.X[:,cols].sum(1)).ravel(),np.asarray((a.X[:,cols]>0).sum(1)).ravel()])
    return obs,p,function,quality,variable_counts,dict(function_rank_positions=positions,partial_rank_positions=partial_positions,proxy_cache_all_parent_values_checked=len(p)*len(g['parent_programme_ids']),
                 proxy_score_sha256=hashlib.sha256(p.astype('<f8').tobytes()).hexdigest(),function_scores_sha256=hashlib.sha256(function.astype('<f8').tobytes()).hexdigest())


def describe(obs,p,function,quality,variable_counts,info,definition,fitter=fit):
    groups=pd.factorize(pd.MultiIndex.from_frame(obs[['node_id','stratum']]),sort=True)[0];rows=[];g=definition['proxy'];fs=definition['functions']
    for st in sorted(set(groups)):
        ix=np.flatnonzero(groups==st);blocks=obs.effective_resampling_block.iloc[ix].astype(str).to_numpy();nb=len(set(blocks));ready=len(ix)>=50 and nb>=10
        structural=np.array([f['measurement_eligible'] and variable_counts[j]>=20 and variable_counts[j]>=.8*f['mapped_genes'] for j,f in enumerate(fs)])
        result=None;status='FEWER_THAN50_CELLS_OR10_BLOCKS'
        if ready and structural.any():result,status=fitter(p[ix],function[ix][:,structural],quality[ix],blocks)
        where={j:k for k,j in enumerate(np.flatnonzero(structural))}
        for j,f in enumerate(fs):
            reason=status if structural[j] else 'INSUFFICIENT_FIXED_OR_VARIABLE_FUNCTION_GENES'
            rr=dict(gse_id=info['gse_id'],major=info['major'],proxy_group_id=g['proxy_group_id'],reference=f['reference'],native_ordinal=info['native_ordinal'],identity_stratum_ordinal=int(st),cells=len(ix),analysis_blocks=nb,
                    mapped_function_genes=f['mapped_genes'],variable_function_genes=variable_counts[j],retained_fraction=f['retained_fraction'],remaining_coverage=f['remaining_coverage'],strong_reference_representation=f['strong_reference_representation'],
                    status=reason,raw_rho=np.nan,adjusted_cell_rho=np.nan,adjusted_block_rho=np.nan)
            if result is not None and structural[j]:
                k=where[j];rr.update(raw_rho=float(result['raw'][k]),adjusted_cell_rho=float(result['cell'][k]),adjusted_block_rho=float(result['block'][k]))
                if not np.isfinite(rr['adjusted_cell_rho']) or not np.isfinite(rr['adjusted_block_rho']):rr['status']='CONSTANT_OR_QUALITY_EXPLAINED_SCORE'
            rows.append(rr)
    return rows


def summarize(table,contract):
    metrics=['raw_rho','adjusted_cell_rho','adjusted_block_rho'];source=[]
    for key,b in table.groupby(['gse_id','major','proxy_group_id','reference'],sort=True):
        z=b[b.status.eq('ESTIMABLE_DESCRIPTION')];r=dict(zip(['gse_id','major','proxy_group_id','reference'],key));r.update(identity_strata=len(b),estimable_identity_strata=len(z),strong_reference_representation=bool(b.strong_reference_representation.iloc[0]))
        for col in metrics:r[col]=float(z[col].median()) if len(z) else np.nan
        source.append(r)
    source=pd.DataFrame(source);assert len(source)==10150;combined=[]
    for key,b in source.groupby(['major','proxy_group_id','reference'],sort=True):
        z=b[b.estimable_identity_strata.gt(0)];a=z.adjusted_cell_rho.to_numpy();v=z.adjusted_block_rho.to_numpy();n=len(z)
        lo=min(min(np.median(np.delete(a,i)),np.median(np.delete(v,i))) for i in range(n)) if n>1 else np.nan
        r=dict(zip(['major','proxy_group_id','reference'],key));r.update(predefined_sources=len(b),estimable_sources=n,median_adjusted_cell_rho=float(np.median(a)) if n else np.nan,median_adjusted_block_rho=float(np.median(v)) if n else np.nan,
             fraction_sources_both_ge_point2=float(np.mean((a>=.2)&(v>=.2))) if n else np.nan,leave_one_min_two_medians=float(lo),strong_reference_representation=bool(b.strong_reference_representation.iloc[0]))
        r['descriptive_concordance_flag']=bool(n>=3 and r['median_adjusted_cell_rho']>=.3 and r['median_adjusted_block_rho']>=.3 and r['fraction_sources_both_ge_point2']>=.75 and lo>=.3)
        combined.append(r)
    combined=pd.DataFrame(combined);assert len(combined)==1550
    return source,combined


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);V=Path(json.loads((root/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root']);out=Path(os.environ['PANCANCER_RAM_OUTPUT']);control=out.parent;plan=json.loads((control/'plan.json').read_text())
    assert not out.exists() and control.parent==Path('/dev/shm') and shutil.disk_usage(control).free>512*1024**2 and sha(__file__)==plan['producer_sha256']
    for n,digest in plan['helper_sha256'].items():assert sha(control/n)==digest
    out.mkdir(mode=0o700);contract=make_contract(root,V);(out/'private_frozen_gene_disjoint_function_contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    panels={}
    for p in json.loads((V/'state_spatial_proxy_calibration_inputs_20261016_v1/fixed_all14_panel_gene_axes.json').read_text()):panels.setdefault(p['gene_axis_sha256'],p)
    smp=V/'state_spatial_partial_panel_calibration_20261016_v2/inventory.json';assert sha(smp)=='f4802973cf3ca061b35c0659389572ed513ae012d8e5957bf0bfeae684d515fd';sf={p['relative_file']:p for p in json.loads(smp.read_text())['private_outputs']}
    allrows=[];proofs=[];processed=0
    for info in contract['native_inputs']:
        defs=[d for d in contract['definitions'] if info['major']==d['proxy']['major'] and info['gse_id'] in d['proxy']['reference_GSEs']]
        for definition in defs:
            obs,p,f,q,var,proof=load_source_and_scores(info,definition,V,panels,sf)
            allrows.extend(describe(obs,p,f,q,var,info,definition));proofs.append(dict(native_ordinal=info['native_ordinal'],proxy_group_id=definition['proxy']['proxy_group_id'],**proof));processed+=1
            (out/'private_progress.json').write_text(json.dumps(dict(stage='ACTUAL_DISJOINT_REFERENCE_RNA_FUNCTIONS',source_proxy_positions_completed=processed,total_positions=203,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    assert processed==203;table=pd.DataFrame(allrows);source,combined=summarize(table,contract)
    table.to_csv(out/'private_all_identity_disjoint_function_descriptions.tsv.gz',sep='\t',index=False,compression='gzip');source.to_csv(out/'private_all10150_source_function_descriptions.tsv.gz',sep='\t',index=False,compression='gzip');combined.to_csv(out/'private_all1550_cross_source_function_descriptions.tsv',sep='\t',index=False)
    (out/'private_all203_score_proofs.json').write_text(json.dumps(proofs,indent=2)+'\n');pd.DataFrame(contract['parent_descriptors']).to_csv(out/'private_all_parent_dictionary_descriptors.tsv',sep='\t',index=False)
    terminal=dict(status='COMPLETE_ALL31_PROXIES50_DISJOINT_FUNCTIONS203_SOURCE_MEASUREMENTS_PENDING_AUDIT',completed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(control/'plan.json'),native_reference_cells=4349,
        source_proxy_positions=203,source_function_positions=10150,proxy_function_positions=1550,identity_positions=len(table),common_unique_native_symbols=len(contract['common_unique_native_symbols']),
        function_rank_positions=sum(p['function_rank_positions'] for p in proofs),partial_rank_positions=sum(p['partial_rank_positions'] for p in proofs),all_cached_parent_proxy_values_checked=sum(p['proxy_cache_all_parent_values_checked'] for p in proofs),
        no_clinical_or_spatial_effect_selection=True,no_new_P_FDR_or_formal_state_or_interaction=True,no_RNA_or_score_arrays_persisted=True,scientific_outputs_public=False,
        files={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in out.iterdir() if p.is_file()})
    assert sum(d['bytes'] for d in terminal['files'].values())<32*1024**2;(out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps({k:terminal[k] for k in ['status','source_proxy_positions','source_function_positions','identity_positions']}))


if __name__=='__main__':main()

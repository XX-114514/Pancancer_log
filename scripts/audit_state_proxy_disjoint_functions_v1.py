"""Complete independent contract, frequency ranks, weighted LS and summaries."""
import hashlib,json,os,sys,shutil
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem

CONTROL=Path(os.environ['PANCANCER_PARENT_CONTROL']);sys.path.insert(0,str(CONTROL))
from measure_state_proxy_disjoint_functions_v1 import load_source_and_scores,describe
from state_proxy_disjoint_function_math_v1 import independent_fit


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def check_contract(c,root,V):
    parentpath=Path('/dev/shm/pancancer-native-state-probe-match-20261020-v1/private_output/private_frozen_matched_measurement_and_target_contract.json')
    assert sha(parentpath)==c['source_SHA']['physical_contract']=='0d6f8e2bac99780b7a7b211c33468ae243ece7920167c3ce85fcf465d994345e'
    p=json.loads(parentpath.read_text());used={m['proxy_group_id'] for t in p['targets'] for m in t['proxy_mappings']};groups={g['proxy_group_id']:g for g in p['groups'] if g['proxy_group_id'] in used}
    assert len(groups)==len(c['definitions'])==31
    axes=[]
    for info in c['native_inputs']:
        assert sha(info['path'])==info['sha256']
        with h5py.File(info['path'],'r') as h:
            k=h['var'].attrs['_index'];k=k.decode() if isinstance(k,bytes) else k;genes=[str(g).strip().upper() for g in read_elem(h['var'][k])]
            cnt=Counter(genes);axes.append({g for g in genes if g and cnt[g]==1});assert int(h['X'].attrs['shape'][0])==info['cells']
    common=sorted(set.intersection(*axes));assert common==c['common_unique_native_symbols']
    dictionaryp=V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json';assert sha(dictionaryp)==c['source_SHA']['dictionary']
    dictionary={r['programme_id']:r for r in json.loads(dictionaryp.read_text())['programmes']}
    pp=V/'state_spatial_proxy_calibration_inputs_20261016_v1/fixed_all14_panel_gene_axes.json';assert sha(pp)==c['source_SHA']['panel_definition'];panels={}
    for x in json.loads(pp.read_text()):panels.setdefault(x['gene_axis_sha256'],x)
    gpath=root/'GSE243013_genes.csv.gz';assert sha(gpath)=='ac325c9ce6bfaf2270eb931e7d947bf31d1fb151b14cb166ed9f6c8f4b34a104'
    names=pd.read_csv(gpath,keep_default_na=False).geneSymbol.str.strip().str.upper().to_numpy();lp=root/'NSCLC_communication_tensor_family_20261019_v1/private_tensor_axes.npz';assert sha(lp)==c['source_SHA']['LR_axes']
    with np.load(lp,allow_pickle=False) as z:LR=set(names[z['native_gene_columns'].astype(int)])
    assert sorted(LR)==c['LR_universe_symbols'] and len(LR)==1840
    rp=root/'NSCLC92_receiver_function_reference_scope_20261019_v1/private_all50_LR_disjoint_reference_definitions.json';assert sha(rp)==c['source_SHA']['reference_definition']=='89ad355c5c6f02f0d0e68b2b0f35cd7e1cc217097512c99e2e306e5f9c92f4fa'
    references={r['reference']:set(str(s).strip().upper() for s in r['full_reference_genes']) for r in json.loads(rp.read_text())};assert len(references)==50
    source_positions=0
    for d in c['definitions']:
        g=d['proxy'];assert g==groups[g['proxy_group_id']];source_positions+=len(g['reference_GSEs']);blocked=set(LR)
        for pid in g['parent_programme_ids']:blocked.update(str(s).strip().upper() for s in dictionary[pid]['primary_genes'])
        p=panels[g['axis_id']];aliases={}
        for s,i in p['symbol_to_feature_index'].items():aliases.setdefault(i,[]).append(s)
        fingerprint={i:hashlib.sha256(json.dumps(sorted(ss),sort_keys=True,separators=(',',':')).encode()).hexdigest() for i,ss in aliases.items()}
        for i,ss in aliases.items():
            if fingerprint[i] in g['signature_probe_keys']:blocked.update(str(s).strip().upper() for s in ss)
        assert sorted(blocked)==d['blocked_symbols'];bg=sorted(set(common)-blocked);assert bg==d['background_symbols']
        assert len(d['functions'])==50 and {f['reference'] for f in d['functions']}==set(references)
        for f in d['functions']:
            full=references[f['reference']];remaining=full-blocked;measured=sorted(remaining&set(bg));assert [bg[j] for j in f['indices']]==measured
            assert len(set(f['indices']))==len(f['indices'])==f['mapped_genes'] and not set(measured)&blocked
            cov=len(measured)/len(remaining) if remaining else 0.;ret=len(remaining)/len(full)
            np.testing.assert_allclose([cov,ret],[f['remaining_coverage'],f['retained_fraction']],atol=0,rtol=0)
            eligible=len(measured)>=20 and cov>=.8
            assert f['measurement_eligible']==eligible and f['strong_reference_representation']==bool(eligible and ret>=.8)
    assert source_positions==203
    expected={(s,g['major']) for g in groups.values() for s in g['reference_GSEs']}
    assert expected=={(r['gse_id'],r['major']) for r in c['native_inputs']}
    return panels


def compare_frames(actual,expected):
    assert len(actual)==len(expected) and list(actual.columns)==list(expected.columns);maxerr=0.
    for col in actual:
        if pd.api.types.is_numeric_dtype(actual[col]) and not pd.api.types.is_bool_dtype(actual[col]):
            a=actual[col].to_numpy(float);b=expected[col].to_numpy(float);np.testing.assert_allclose(a,b,atol=2e-8,rtol=1e-7,equal_nan=True)
            finite=np.isfinite(a)&np.isfinite(b)
            if finite.any():maxerr=max(maxerr,float(np.max(np.abs(a[finite]-b[finite]))))
        else:np.testing.assert_array_equal(actual[col].to_numpy(),expected[col].to_numpy())
    return maxerr


def independent_summaries(table):
    keys=['gse_id','major','proxy_group_id','reference'];valid=table.status.eq('ESTIMABLE_DESCRIPTION')
    base=table.groupby(keys,sort=True).agg(identity_strata=('status','size'),strong_reference_representation=('strong_reference_representation','first'))
    good=table[valid].groupby(keys,sort=True).agg(estimable_identity_strata=('status','size'),raw_rho=('raw_rho','median'),adjusted_cell_rho=('adjusted_cell_rho','median'),adjusted_block_rho=('adjusted_block_rho','median'))
    base=base.join(good);base['estimable_identity_strata']=base.estimable_identity_strata.fillna(0).astype(int)
    source=base.reset_index()[keys+['identity_strata','estimable_identity_strata','strong_reference_representation','raw_rho','adjusted_cell_rho','adjusted_block_rho']]
    result=[]
    for key,b in source.groupby(['major','proxy_group_id','reference'],sort=True):
        z=b[b.estimable_identity_strata.gt(0)];n=len(z);cell=z.adjusted_cell_rho.to_numpy();block=z.adjusted_block_rho.to_numpy()
        ca=float(pd.Series(cell).median()) if n else np.nan;ba=float(pd.Series(block).median()) if n else np.nan
        fraction=float(sum(a>=.2 and v>=.2 for a,v in zip(cell,block))/n) if n else np.nan
        values=[]
        for i in range(n):
            if n>1:values.extend([pd.Series([v for j,v in enumerate(cell) if j!=i]).median(),pd.Series([v for j,v in enumerate(block) if j!=i]).median()])
        lo=float(min(values)) if values else np.nan
        r=dict(zip(['major','proxy_group_id','reference'],key));r.update(predefined_sources=len(b),estimable_sources=n,median_adjusted_cell_rho=ca,median_adjusted_block_rho=ba,fraction_sources_both_ge_point2=fraction,leave_one_min_two_medians=lo,strong_reference_representation=bool(b.strong_reference_representation.iloc[0]),descriptive_concordance_flag=bool(n>=3 and ca>=.3 and ba>=.3 and fraction>=.75 and lo>=.3));result.append(r)
    return source,pd.DataFrame(result)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=Path(os.environ['PANCANCER_RAM_OUTPUT']);parent=CONTROL/'private_output';plan=json.loads((CONTROL/'plan.json').read_text())
    assert not out.exists() and out.parent.parent==Path('/dev/shm') and shutil.disk_usage(out.parent).free>512*1024**2
    term=json.loads((parent/'private_terminal.json').read_text());assert term['status']=='COMPLETE_ALL31_PROXIES50_DISJOINT_FUNCTIONS203_SOURCE_MEASUREMENTS_PENDING_AUDIT'
    assert sha(CONTROL/'measure_state_proxy_disjoint_functions_v1.py')==plan['producer_sha256']==term['producer_sha256']
    for name,digest in plan['helper_sha256'].items():assert sha(CONTROL/name)==digest
    for name,d in term['files'].items():assert sha(parent/name)==d['sha256']
    V=Path(json.loads((root/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root']);c=json.loads((parent/'private_frozen_gene_disjoint_function_contract.json').read_text());panels=check_contract(c,root,V)
    manifest=V/'state_spatial_partial_panel_calibration_20261016_v2/inventory.json';assert sha(manifest)=='f4802973cf3ca061b35c0659389572ed513ae012d8e5957bf0bfeae684d515fd';sf={r['relative_file']:r for r in json.loads(manifest.read_text())['private_outputs']}
    oldproof={(p['native_ordinal'],p['proxy_group_id']):p for p in json.loads((parent/'private_all203_score_proofs.json').read_text())};assert len(oldproof)==203
    out.mkdir(mode=0o700);rows=[];processed=0;rank_positions=0;partial_positions=0
    for info in c['native_inputs']:
        for d in c['definitions']:
            g=d['proxy']
            if g['major']!=info['major'] or info['gse_id'] not in g['reference_GSEs']:continue
            obs,p,f,q,var,proof=load_source_and_scores(info,d,V,panels,sf,independent=True)
            for k,v in proof.items():assert v==oldproof[(info['native_ordinal'],g['proxy_group_id'])][k]
            rank_positions+=proof['function_rank_positions'];partial_positions+=proof['partial_rank_positions'];rows.extend(describe(obs,p,f,q,var,info,d,fitter=independent_fit));processed+=1
            (out/'private_progress.json').write_text(json.dumps(dict(stage='ALL_FREQUENCY_RANKS_AND_INDEPENDENT_WEIGHTED_LS',source_proxy_positions=processed,total_positions=203,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    assert processed==203 and rank_positions==term['function_rank_positions'] and partial_positions==term['partial_rank_positions']
    table=pd.DataFrame(rows);original=pd.read_csv(parent/'private_all_identity_disjoint_function_descriptions.tsv.gz',sep='\t');err=compare_frames(table,original)
    source,combined=independent_summaries(table);err=max(err,compare_frames(source,pd.read_csv(parent/'private_all10150_source_function_descriptions.tsv.gz',sep='\t')),compare_frames(combined,pd.read_csv(parent/'private_all1550_cross_source_function_descriptions.tsv',sep='\t')))
    assert len(source)==10150 and len(combined)==1550
    receipt=dict(status='ACCEPTED_ALL203_NATIVE_SCORE_RECONSTRUCTIONS_ALL_IDENTITY_LS_AND10150_SOURCE1550_AGGREGATES',completed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_terminal_sha256=sha(parent/'private_terminal.json'),
                 all_source_proxy_scores_checked=processed,function_rank_positions_reconstructed=rank_positions,partial_rank_positions_reconstructed=partial_positions,identity_positions_recomputed=len(table),source_function_positions_recomputed=len(source),proxy_function_positions_recomputed=len(combined),
                 maximum_numeric_error=err,all_score_vector_sha256_exact=True,common_background_and_all_reference_gene_exclusions_independently_rebuilt=True,
                 reused_production_input_loader_and_stratum_enumeration=True,independent_frequency_ranking_and_weighted_least_squares=True,independent_raw_archive_replay=False,
                 no_arrays_persisted=True,scientific_review_pending=True,formal_state_or_interaction_or_clinical_validation=False,scientific_outputs_public=False)
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:receipt[k] for k in ['status','all_source_proxy_scores_checked','identity_positions_recomputed','maximum_numeric_error']}))


if __name__=='__main__':main()

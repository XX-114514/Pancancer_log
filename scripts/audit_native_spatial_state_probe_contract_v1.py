"""Independently rebuild probe profiles, selection and all target capture counts."""
import hashlib,json,os,shutil
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def code(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def main():
    parent=Path('/dev/shm/pancancer-native-state-probe-match-20261020-v1/private_output')
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    assert not out.exists() and out.parent.parent==Path('/dev/shm') and shutil.disk_usage('/dev/shm').free>512*1024**2
    assert sha(parent/'private_terminal.json')=='7f275021fc2137b11b589b7fe560a8c8385c51059ea94a3dad5a4b557b2114e0'
    term=json.loads((parent/'private_terminal.json').read_text())
    for name,d in term['files'].items():assert sha(parent/name)==d['sha256']
    V=Path(json.loads((root/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    panels=json.loads((V/'state_spatial_proxy_calibration_inputs_20261016_v1/fixed_all14_panel_gene_axes.json').read_text())
    byaxis={};bysample={x['sample_id']:x for x in panels}
    for p in panels:byaxis.setdefault(p['gene_axis_sha256'],p)
    recs={x['programme_id']:x for x in json.loads((V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json').read_text())['programmes']}
    ni=json.loads((V/'state_spatial_proxy_native_cells_20261016_v1/inventory.json').read_text())
    native={(x['gse_id'],x['major']):x for x in ni['private_native_files']}
    profiles=pd.read_csv(parent/'private_all_reference_source_measurement_profiles.tsv',sep='\t')
    expected=pd.read_csv(parent/'private_all_previously_qualified_T_Mac_matched_profiles.tsv',sep='\t')
    sources=pd.read_csv('/dev/shm/pancancer-state-proxy-quality-20261020-v2/private_output/private_all_source_major_programme_axis_fidelity.tsv.gz',sep='\t')
    q=pd.read_csv('/dev/shm/pancancer-state-proxy-quality-20261020-v2/private_output/private_complete_programme_major_axis_triage.tsv',sep='\t')
    original=q[q.exploratory_proxy_triage_pass&q.major.isin(['T cell','Macrophage'])]
    assert set(zip(expected.programme_id,expected.major,expected.axis_id))==set(zip(original.programme_id,original.major,original.axis_id))
    gene_cache={};backgrounds={};counts=defaultdict(list);definitions={}
    for row in profiles.itertuples(index=False):
        nk=(row.gse_id,row.major)
        if nk not in gene_cache:
            inf=native[nk];p=V/'state_spatial_proxy_native_cells_20261016_v1'/inf['relative_file'];assert sha(p)==inf['sha256']
            with h5py.File(p,'r') as h:
                k=h['var'].attrs['_index'];k=k.decode() if isinstance(k,bytes) else k;gene_cache[nk]=set(np.asarray(read_elem(h['var'][k]),str))
        bk=(*nk,row.axis_id)
        if bk not in backgrounds:
            panel=byaxis[row.axis_id];tmp=defaultdict(set)
            for name,j in panel['symbol_to_feature_index'].items():tmp[j].add(name)
            chosen=[];mapping={}
            for j,names in sorted(tmp.items()):
                hits=names&gene_cache[nk];preferred=panel['preferred_symbol_by_feature'][j]
                selected=preferred if preferred in hits else next(iter(hits)) if len(hits)==1 else None
                if selected is not None:chosen.append(selected);mapping[selected]=code(sorted(names))
            usable={x for x in chosen if chosen.count(x)==1}
            backgrounds[bk]=(usable,mapping,sorted(mapping[x] for x in usable))
        usable,mapping,bg=backgrounds[bk]
        sig=sorted(mapping[x] for x in set(recs[row.programme_id]['primary_genes'])&usable)
        definition=dict(background_probe_keys=bg,signature_probe_keys=sig,maxrank=min(1500,len(bg)))
        assert code(definition)==row.measurement_fingerprint and len(bg)==row.background_probes and len(sig)==row.signature_probes
        key=(row.programme_id,row.major,row.axis_id);counts[key].append((row.gse_id,row.measurement_fingerprint));definitions[(key,row.measurement_fingerprint)]=definition
    passed=set()
    for row in expected.itertuples(index=False):
        key=(row.programme_id,row.major,row.axis_id);g=sources[sources.programme_id.eq(row.programme_id)&sources.major.eq(row.major)&sources.axis_id.eq(row.axis_id)&sources.estimable_identity_strata.gt(0)&sources.partial_signature_genes.ge(5)]
        recorded=counts[key];assert len(recorded)==len(set(x[0] for x in recorded))==len(g)
        assert set(x[0] for x in recorded)==set(g.gse_id)
        group=defaultdict(set)
        for src,fingerprint in recorded:group[fingerprint].add(src)
        winner=sorted(group,key=lambda k:(-len(group[k]),k))[0]
        assert winner==row.chosen_profile and len(group)==row.distinct_measurement_profiles
        selected=g[g.gse_id.isin(group[winner])];a=selected.adjusted_cell_rho.to_numpy();b=selected.adjusted_block_rho.to_numpy();n=len(a)
        am=float(np.median(a));bm=float(np.median(b));fraction=float(np.mean((a>=.2)&(b>=.2)))
        loo=min(min(np.median(np.delete(a,k)),np.median(np.delete(b,k))) for k in range(n)) if n>1 else np.nan
        result=bool(recs[row.programme_id]['algorithm_stable_candidate'] and n>=3 and am>=.3 and bm>=.3 and fraction>=.75 and loo>=.3)
        assert result==row.matched_profile_fidelity_pass and n==row.profile_sources
        np.testing.assert_allclose([am,bm,fraction,loo],[row.median_cell_rho,row.median_block_rho,row.source_fraction_both_ge_point2,row.leave_one_min_median],atol=2e-12,rtol=0,equal_nan=True)
        if result:passed.add(key)
    contract=json.loads((parent/'private_frozen_matched_measurement_and_target_contract.json').read_text())
    mapped=set();measurements=set()
    for g in contract['groups']:
        key=(g['axis_id'],g['major'],g['measurement_fingerprint']);assert key not in measurements;measurements.add(key)
        d={k:g[k] for k in ['background_probe_keys','signature_probe_keys','maxrank']};assert code(d)==g['measurement_fingerprint']
        for pid in g['parent_programme_ids']:
            pk=(pid,g['major'],g['axis_id']);assert pk not in mapped;mapped.add(pk);assert definitions[(pk,g['measurement_fingerprint'])]==d
    assert mapped==passed
    fields=pd.read_csv(parent/'private_all_author_FOV_capture_counts.tsv.gz',sep='\t',dtype={'sample_id':str,'FOV':str})
    field_checked=0;mapping_checked=0
    for target in contract['targets']:
        assert sha(target['analysis_path'])==target['input_sha256'] and sha(target['cache_path'])==target['cache_sha256']
        with np.load(target['cache_path'],allow_pickle=False) as z:
            obs=pd.DataFrame({'label':z['author_labels'].astype(str),'fov':z['fov'].astype(str),'qc':z['qc'].astype(bool)})
        assert len(obs)==target['native_cells']
        table=pd.crosstab(obs.loc[obs.qc,'fov'],obs.loc[obs.qc,'label'])
        for row in fields[fields.sample_id.eq(target['sample_id'])].itertuples(index=False):
            cnt=table.loc[row.FOV] if row.FOV in table.index else pd.Series(dtype=int)
            nt=sum(int(cnt.get(k,0)) for k in contract['T_author_labels']);nm=sum(int(cnt.get(k,0)) for k in contract['strict_macrophage_author_labels'])
            assert nt==row.author_T_QC_cells and nm==row.strict_author_Mac_QC_cells and (nt>=30 and nm>=30)==row.both_atleast30;field_checked+=1
        panel=bysample[target['sample_id']];tmp=defaultdict(list)
        for name,j in panel['symbol_to_feature_index'].items():tmp[j].append(name)
        keybyidx={j:code(sorted(names)) for j,names in tmp.items()}
        for m in target['proxy_mappings']:
            g=next(x for x in contract['groups'] if x['proxy_group_id']==m['proxy_group_id'])
            assert [keybyidx[j] for j in m['background_indices']]==g['background_probe_keys']
            assert [keybyidx[j] for j in m['signature_indices']]==g['signature_probe_keys']
            assert target['axis_id']==g['axis_id'] and m['maxrank']==g['maxrank'];mapping_checked+=1
    assert field_checked==len(fields) and len(contract['targets'])==4
    out.mkdir(mode=0o700)
    report=dict(status='ACCEPTED_ALL_PROBE_PROFILES_WINNERS_FIXED_GATES_DEDUP_TARGET_MAPPINGS_AND_AUTHOR_FIELD_COUNTS',completed_UTC=datetime.now(timezone.utc).isoformat(),
                auditor_sha256=sha(__file__),producer_terminal_sha256=sha(parent/'private_terminal.json'),
                programme_positions_checked=len(expected),source_profile_positions_checked=len(profiles),native_gene_axes_independently_read=len(gene_cache),
                author_target_images_checked=4,target_mapping_positions_checked=mapping_checked,technical_field_counts_checked=field_checked,
                new_RNA_scores_P_or_spatial_effects=False,scientific_outputs_public=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()

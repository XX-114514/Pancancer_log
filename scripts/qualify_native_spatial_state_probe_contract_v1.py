"""Freeze identical physical-probe measurements before new spatial effects."""
import hashlib,json,os,shutil
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem

T_AUTHOR=['CD4+ T cells','CD8+ T cells','T cell proliferating','T cells','Treg']
MAC_AUTHOR=['Alveolar macrophages','Alveolar mph proliferating','Interstitial mph prevascular','Macrophages']


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def resolved(panel,native):
    aliases=defaultdict(list)
    for gene,j in panel['symbol_to_feature_index'].items():
        if gene in native:aliases[j].append(gene)
    chosen=[]
    for j,names in sorted(aliases.items()):
        pref=panel['preferred_symbol_by_feature'][j]
        if pref in names:chosen.append(pref)
        elif len(names)==1:chosen.append(names[0])
    repeat={g for g in chosen if chosen.count(g)>1}
    return sorted(g for g in chosen if g not in repeat)


def physical_keys(panel):
    aliases=defaultdict(list)
    for gene,j in panel['symbol_to_feature_index'].items():aliases[j].append(gene)
    out={j:digest(sorted(names)) for j,names in aliases.items()}
    assert len(set(out.values()))==len(out)
    return out


def fidelity(g,stable):
    a=g.adjusted_cell_rho.to_numpy(float);b=g.adjusted_block_rho.to_numpy(float);n=len(g)
    am=float(np.median(a)) if n else np.nan;bm=float(np.median(b)) if n else np.nan
    fraction=float(np.mean((a>=.2)&(b>=.2))) if n else np.nan
    loo=min(min(np.median(np.delete(a,k)),np.median(np.delete(b,k))) for k in range(n)) if n>1 else np.nan
    okay=bool(stable and n>=3 and am>=.3 and bm>=.3 and fraction>=.75 and loo>=.3)
    return dict(profile_sources=n,median_cell_rho=am,median_block_rho=bm,source_fraction_both_ge_point2=fraction,
                leave_one_min_median=float(loo),matched_profile_fidelity_pass=okay)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
    planp=Path(os.environ['PANCANCER_STATE_NATIVE_PROBE_PLAN']);plan=json.loads(planp.read_text())
    assert out.parent.parent==Path('/dev/shm') and not out.exists() and shutil.disk_usage('/dev/shm').free>512*1024**2
    assert sha(__file__)==plan['producer_sha256']
    V=Path(json.loads((root/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    qroot=Path('/dev/shm/pancancer-state-proxy-quality-20261020-v2/private_output')
    audit=Path('/dev/shm/pancancer-state-proxy-quality-audit-20261020-v2/private_output/private_terminal.json')
    assert sha(audit)=='d0221b7719999d7678d79fc33bdf8f8da1e3f0602f4954ef941ac956b9b13759'
    qt=json.loads((qroot/'private_terminal.json').read_text())
    for name in ['private_complete_programme_major_axis_triage.tsv','private_all_source_major_programme_axis_fidelity.tsv.gz']:
        assert sha(qroot/name)==qt['files'][name]['sha256']
    q=pd.read_csv(qroot/'private_complete_programme_major_axis_triage.tsv',sep='\t')
    sources=pd.read_csv(qroot/'private_all_source_major_programme_axis_fidelity.tsv.gz',sep='\t')
    q=q[q.exploratory_proxy_triage_pass&q.major.isin(['T cell','Macrophage'])].copy()
    dp=V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json'
    assert sha(dp)=='6a3c4311d89a9a02d55a19f2d249c960c508b5973368c06c2e5c43bb605edec9'
    recs={r['programme_id']:r for r in json.loads(dp.read_text())['programmes']}
    pp=V/'state_spatial_proxy_calibration_inputs_20261016_v1/fixed_all14_panel_gene_axes.json'
    panels=json.loads(pp.read_text());assert len(panels)==14;byaxis={}
    for p in panels:byaxis.setdefault(p['gene_axis_sha256'],p)
    ni=V/'state_spatial_proxy_native_cells_20261016_v1/inventory.json'
    assert sha(ni)=='91be1ae48a338aa6bfcd2d89107fd62bf3a00566a90ec5aa953fe0a3f37478f0'
    native={(r['gse_id'],r['major']):r for r in json.loads(ni.read_text())['private_native_files']}
    assert len(native)==709
    out.mkdir(mode=0o700);save=lambda name,d:(out/name).write_text(json.dumps(d,indent=2)+'\n')
    save('plan_before_physical_measurement_matching.json',plan)
    profiles={};checked_files={};allrows=[];definitions=[];source_profiles=[]
    for row in q.sort_values(['axis_id','major','programme_id']).itertuples(index=False):
        p=byaxis[row.axis_id];pk=physical_keys(p)
        eligible=sources[sources.programme_id.eq(row.programme_id)&sources.major.eq(row.major)&sources.axis_id.eq(row.axis_id)&sources.estimable_identity_strata.gt(0)&sources.partial_signature_genes.ge(5)]
        assert eligible.gse_id.is_unique
        groups=defaultdict(list);groupdef={}
        for src in eligible.itertuples(index=False):
            key=(src.gse_id,src.major,row.axis_id)
            if key not in profiles:
                info=native[(src.gse_id,src.major)];path=V/'state_spatial_proxy_native_cells_20261016_v1'/info['relative_file']
                if str(path) not in checked_files:
                    assert sha(path)==info['sha256'];checked_files[str(path)]=info['sha256']
                with h5py.File(path,'r') as h:
                    index=h['var'].attrs['_index'];index=index.decode() if isinstance(index,bytes) else index
                    genes=np.asarray(read_elem(h['var'][index]),str)
                assert len(set(genes))==len(genes)
                names=resolved(p,set(genes));background=sorted(pk[p['symbol_to_feature_index'][g]] for g in names)
                assert len(background)==len(set(background))
                profiles[key]=(names,background)
            names,background=profiles[key]
            sig=sorted(pk[p['symbol_to_feature_index'][g]] for g in recs[row.programme_id]['primary_genes'] if g in names)
            assert len(sig)==len(set(sig))==src.partial_signature_genes
            definition=dict(background_probe_keys=background,signature_probe_keys=sig,maxrank=min(1500,len(background)))
            code=digest(definition);groups[code].append(src.gse_id);groupdef[code]=definition
            source_profiles.append(dict(gse_id=src.gse_id,major=src.major,programme_id=row.programme_id,axis_id=row.axis_id,
                                        measurement_fingerprint=code,background_probes=len(background),signature_probes=len(sig)))
        codes=sorted(groups,key=lambda k:(-len(groups[k]),k))
        if codes:
            chosen=codes[0];use=eligible[eligible.gse_id.isin(groups[chosen])]
            metrics=fidelity(use,bool(recs[row.programme_id]['algorithm_stable_candidate']))
            definition=dict(axis_id=row.axis_id,major=row.major,programme_id=row.programme_id,
                            measurement_fingerprint=chosen,**groupdef[chosen],reference_GSEs=sorted(groups[chosen]))
            definitions.append(definition)
        else:chosen='';metrics=fidelity(eligible.iloc[:0],False)
        allrows.append(dict(programme_id=row.programme_id,major=row.major,axis_id=row.axis_id,previous_eligible_sources=len(eligible),
                            distinct_measurement_profiles=len(groups),chosen_profile=chosen,**metrics))
    triage=pd.DataFrame(allrows);assert len(triage)==len(q)
    triage.to_csv(out/'private_all_previously_qualified_T_Mac_matched_profiles.tsv',sep='\t',index=False)
    pd.DataFrame(source_profiles).to_csv(out/'private_all_reference_source_measurement_profiles.tsv',sep='\t',index=False)
    passing=set(zip(triage.loc[triage.matched_profile_fidelity_pass,'programme_id'],triage.loc[triage.matched_profile_fidelity_pass,'axis_id']))
    passed=[x for x in definitions if (x['programme_id'],x['axis_id']) in passing]
    dedup={}
    for x in passed:
        key=(x['axis_id'],x['major'],x['measurement_fingerprint'])
        if key not in dedup:dedup[key]=dict(x,parent_programme_ids=[])
        dedup[key]['parent_programme_ids'].append(x['programme_id'])
    groups=[]
    for j,(key,x) in enumerate(sorted(dedup.items())):
        x.pop('programme_id');x['proxy_group_id']='P%03d'%j;x['parent_programme_ids']=sorted(x['parent_programme_ids']);groups.append(x)
    manifest=root/'immune_same_name_family_source_review_v1/spatial_native_input_repair_v2/analysis_input_manifest.tsv'
    assert sha(manifest)=='6fda649c8b2bf4abe010ba8d53d50e495bae29ce8d30648675d9be9bf91dcf83'
    mf=pd.read_csv(manifest,sep='\t',keep_default_na=False).set_index('sample_id')
    auth=root/'all14_author_spatial_lineage_qualification_20261019_v1';assert sha(auth/'private_terminal.json')=='e68ea3caf60fc05a09abd0072ed9155aa8e09b4d27a19b72fe82f010a5ec1c9a'
    at=json.loads((auth/'private_terminal.json').read_text());name='private_complete14_author_lineage_and_assay_scope.tsv';assert sha(auth/name)==at['files'][name]['sha256']
    authors=pd.read_csv(auth/name,sep='\t');assert len(authors)==14
    old=root/'immune_same_name_family_source_review_v1/spatial_contact_native_preflight_v3';oi=json.loads((old/'inventory.json').read_text())
    assert sha(old/'inventory.json')==json.loads((old/'independent_audit.json').read_text())['producer_inventory_sha256']
    panels={p['sample_id']:p for p in panels};image_rows=[];field_rows=[];target_contract=[]
    for src in authors.itertuples(index=False):
        p=panels[src.sample_id];axis=p['gene_axis_sha256']
        record=dict(dataset_id=src.dataset_id,cancer=src.cancer,sample_id=src.sample_id,axis_id=axis,
                    author_and_geometry_available=bool(src.author_labels_available and src.selected_geometry_has_source_proof),
                    native_cells=int(src.native_cells),T_QC_cells=np.nan,strict_macrophage_QC_cells=np.nan,
                    FOVs_T_Mac_ge30=np.nan,matched_T_proxy_groups=0,matched_Mac_proxy_groups=0,
                    status='AUTHOR_OR_GEOMETRY_UNAVAILABLE_NOT_CELL_ABSENCE')
        if record['author_and_geometry_available']:
            assert src.source_input_sha256==p['input_sha256']
            path=Path(mf.loc[src.sample_id,'analysis_path']);assert sha(path)==src.source_input_sha256
            rel='private_native_cell_caches/'+src.sample_id+'.npz';cache=old/rel;assert sha(cache)==oi['outputs'][rel]
            with np.load(cache,allow_pickle=False) as z:
                labels=z['author_labels'].astype(str);qc=z['qc'].astype(bool);fov=z['fov'].astype(str)
            assert len(labels)==len(qc)==len(fov)==record['native_cells']
            T=np.isin(labels,T_AUTHOR)&qc;M=np.isin(labels,MAC_AUTHOR)&qc;assert not (T&M).any()
            for field in sorted(set(fov)):
                use=fov==field;nt=int((T&use).sum());nm=int((M&use).sum())
                field_rows.append(dict(dataset_id=src.dataset_id,cancer=src.cancer,sample_id=src.sample_id,FOV=field,
                                       author_T_QC_cells=nt,strict_author_Mac_QC_cells=nm,both_atleast30=nt>=30 and nm>=30))
            target_key_to_index={k:j for j,k in physical_keys(p).items()};mapping=[]
            for g in groups:
                if g['axis_id']!=axis:continue
                assert set(g['background_probe_keys'])<=set(target_key_to_index)
                mapping.append(dict(proxy_group_id=g['proxy_group_id'],major=g['major'],
                                    background_indices=[target_key_to_index[k] for k in g['background_probe_keys']],
                                    signature_indices=[target_key_to_index[k] for k in g['signature_probe_keys']],
                                    maxrank=g['maxrank']))
            nt=sum(g['major']=='T cell' for g in mapping);nm=sum(g['major']=='Macrophage' for g in mapping)
            record.update(T_QC_cells=int(T.sum()),strict_macrophage_QC_cells=int(M.sum()),
                          FOVs_T_Mac_ge30=sum(x['both_atleast30'] for x in field_rows if x['sample_id']==src.sample_id),
                          matched_T_proxy_groups=nt,matched_Mac_proxy_groups=nm,status='NATIVE_METADATA_MATCHED_MEASUREMENT_SCOPE_ONLY')
            target_contract.append(dict(sample_id=src.sample_id,dataset_id=src.dataset_id,cancer=src.cancer,axis_id=axis,
                                        analysis_path=str(path),input_sha256=sha(path),cache_path=str(cache),cache_sha256=sha(cache),
                                        native_cells=record['native_cells'],proxy_mappings=mapping))
        image_rows.append(record)
    assert len(image_rows)==14 and len(target_contract)==4
    pd.DataFrame(image_rows).to_csv(out/'private_complete14_native_target_scope.tsv',sep='\t',index=False)
    fields=pd.DataFrame(field_rows);fields.to_csv(out/'private_all_author_FOV_capture_counts.tsv.gz',sep='\t',index=False,compression='gzip')
    save('private_frozen_matched_measurement_and_target_contract.json',dict(groups=groups,targets=target_contract,T_author_labels=T_AUTHOR,strict_macrophage_author_labels=MAC_AUTHOR,
            parent_dictionary_sha256=sha(dp),panel_definition_sha256=sha(pp),source_quality_terminal_sha256=sha(qroot/'private_terminal.json'),
            source_quality_audit_sha256=sha(audit),selection='Largest eligible-GSE group of identical physical background/signature probes; fingerprint lexical tie. Original fixed fidelity thresholds unchanged.',
            next_target_score_readout='Same measured background/maxrank and positive signature UCell; no phenotype, spatial effect or formal state accepted here.'))
    images=pd.DataFrame(image_rows)
    images.groupby(['dataset_id','cancer'],sort=True).agg(source_inputs=('sample_id','size'),known_author_geometry_inputs=('author_and_geometry_available','sum'),
              T_QC_cells=('T_QC_cells',lambda z:z.sum(min_count=1)),strict_Mac_QC_cells=('strict_macrophage_QC_cells',lambda z:z.sum(min_count=1)),
              technical_FOVs_T_Mac_ge30=('FOVs_T_Mac_ge30',lambda z:z.sum(min_count=1))).reset_index().to_csv(out/'private_complete_source_target_capture_summary.tsv',sep='\t',index=False)
    triage.groupby(['major','axis_id'],sort=True).agg(previous_programme_positions=('programme_id','size'),
                matched_measurement_pass_positions=('matched_profile_fidelity_pass','sum')).reset_index().to_csv(out/'private_major_axis_matched_measurement_counts.tsv',sep='\t',index=False)
    terminal=dict(status='COMPLETE_ALL_PREQUALIFIED_PROFILES_AND14_NATIVE_AUTHOR_TARGETS_NO_NEW_EXPRESSION_OR_SPATIAL_EFFECT',
                  completed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planp),
                  parent_quality_audit_sha256=sha(audit),all_target_inputs=14,known_author_geometry_inputs=4,
                  profile_programme_positions=len(triage),current_native_reference_files_SHA_checked=len(checked_files),
                  no_new_expression_scores_P_or_clinical_outcomes=True,no_formal_identity_or_full_signature_gate_changed=True,
                  scientific_outputs_public=False,files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(x['bytes'] for x in terminal['files'].values())<8*1024**2
    save('private_terminal.json',terminal);print(json.dumps({k:terminal[k] for k in ['status','all_target_inputs','known_author_geometry_inputs','profile_programme_positions','current_native_reference_files_SHA_checked']}))


if __name__=='__main__':main()

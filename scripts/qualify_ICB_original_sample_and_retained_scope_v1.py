"""Current metadata-only donor/phase versus original sample-unit qualification."""
import hashlib
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem

FIELDS=['Study_name','PMID_donor_id','pre_post','Sample','outcome',
        'Cancer_type_update','Treatment.or.Mode.of.Action','tissue_type',
        'author_cell_type_update','PMID']
T_LABELS=['CD4_Tcells','CD8_Tcells','CD8act_Tcells','CD8eff_Tcell','CD8ex_Tcells',
          'CD8mem_Tcells','CD8naive_Tcell','Naive_Tcell','Prolif_Tcell','T_cell','Tfh','Th17','Treg']
M_LABELS=['Macrophage','TAM','Monocyte','Myeloid_cell','DCs','pDCs','Microglia']
MISSING={'','nan','none','na','n/a','unknown','not_reviewed','__source_value_missing__'}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def valid(s):
    return ~s.str.strip().str.casefold().isin(MISSING)


def summarize(f):
    f=f.copy()
    f['t']=f.author_cell_type_update.isin(T_LABELS).astype(int)
    f['m']=f.author_cell_type_update.isin(M_LABELS).astype(int)
    f['one']=1
    pre=f[f.pre_post.eq('Pre') & valid(f.PMID_donor_id)]
    groups=['PMID_donor_id','pre_post']
    g=pre.groupby(groups,sort=True).agg(cells=('one','sum'),T=('t','sum'),M=('m','sum'),samples=('Sample',lambda s:s[valid(s)].nunique()))
    sample=pre[valid(pre.Sample)].groupby(groups+['Sample'],sort=True)[['t','m']].sum()
    good=sample[(sample.t>=20)&(sample.m>=20)].reset_index()
    same=set(zip(good.PMID_donor_id,good.pre_post))
    both=(g.T>=20)&(g.M>=20)
    if len(pre):
        coded=pd.DataFrame({'d':pd.Categorical(pre.PMID_donor_id).codes,'p':pd.Categorical(pre.pre_post).codes})
        assert len(np.unique(coded.to_numpy(),axis=0))==len(g)
    conflicts={}
    for field in ['outcome','Cancer_type_update','Treatment.or.Mode.of.Action','tissue_type']:
        counts=pre.groupby(groups)[field].apply(lambda s:s[valid(s)].nunique())
        conflicts[field+'_within_donor_phase_conflicts']=int((counts>1).sum())
    literal=pre.groupby(groups).outcome.apply(lambda s:s[valid(s)].nunique())
    ambiguous=pre[valid(pre.Sample)].groupby('Sample').PMID_donor_id.nunique()
    # Distinct author tokens are not automatically reviewed independent patients.
    return dict(available_rows=len(f),baseline_rows=len(pre),baseline_author_identifier_tokens=len(g),
                baseline_original_sample_tokens=int(pre.loc[valid(pre.Sample),'Sample'].nunique()),
                baseline_donor_phase_with_multiple_original_samples=int((g.samples>1).sum()),
                baseline_donor_phase_T_and_broad_myeloid_ge20=int(both.sum()),
                baseline_donor_phase_with_one_original_sample_both_lineages_ge20=len(same),
                baseline_pair_eligibility_only_after_merging_samples=sum(bool(ok) and key not in same for key,ok in both.items()),
                original_sample_token_assigned_multiple_donors=int((ambiguous>1).sum()),
                baseline_rows_missing_original_sample_token=int((~valid(pre.Sample)).sum()),
                baseline_author_units_with_one_literal_outcome=int((literal==1).sum()),
                baseline_author_units_missing_literal_outcome=int((literal==0).sum()),
                **conflicts,clinical_endpoint_or_patient_identity_accepted=False,
                natural_population_composition_eligible=False)


def old_capture(f):
    records=[]
    for key,g in f.groupby(['PMID_donor_id','pre_post'],sort=True):
        if key[1]!='Pre':continue
        counts=g.author_cell_type_update.value_counts()
        records.append((sum(int(counts.get(v,0)) for v in ['Macrophage','TAM']),int(counts.get('Myeloid_cell',0)),sum(int(n) for v,n in counts.items() if v.startswith('CD8')),int(counts.get('T_cell',0))))
    a=np.asarray(records,dtype=int).reshape(-1,4)
    return dict(baseline_patients=len(a),baseline_atleast20_author_macrophage=int((a[:,0]>=20).sum()),baseline_atleast20_unsplit_myeloid=int((a[:,1]>=20).sum()),baseline_atleast20_both_macrophage_CD8=int(((a[:,0]>=20)&(a[:,2]>=20)).sum()),baseline_atleast20_both_unsplit_myeloid_T=int(((a[:,1]>=20)&(a[:,3]>=20)).sum()))


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent==Path('/dev/shm') and not out.exists() and shutil.disk_usage('/dev/shm').free>512*1024**2
    contract=root.parents[2]/'Project_v5/runs/20261013_immunotherapy_atlas_v1/config/source_contract_v2_private.json'
    assert sha(contract)=='8d7f54a385389d0ecbc825155df727ed065cae2ed1c3f3ed584c1a7896647e71'
    records=[r for r in json.loads(contract.read_text())['records'] if r['dataset'].startswith('ICB_')]
    assert len(records)==8 and len({r['path'] for r in records})==1
    old=root/'ICB_clinical_programme_scope_20261014_v1/aggregate_terminal.json'
    old_scope={r['source']:r for r in json.loads(old.read_text())['sources']}
    out.mkdir(mode=0o700);save=lambda name,value:(out/name).write_text(json.dumps(value,indent=2)+'\n')
    save('plan_before_new_original_sample_unit_counts.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),source_contract_sha256=sha(contract),old_scope_sha256=sha(old),all_source_positions=8,populations=['contract_retained','all_available_same_study_in_cache'],metadata_fields=FIELDS,T_author_labels=T_LABELS,broad_myeloid_author_labels=M_LABELS,missing_tokens=sorted(MISSING),prior_closed_Bassez_family_unchanged=True,Combined_outcome_not_read=True,scope='Metadata key aliases/endpoint consistency and retained-subset effects only. Same-cache complete study rows remain compiler-filtered, not natural populations. No RNA, P, new clinical models, patient identity promotion or source-quarantine removal.',output_cap_bytes=1024**2))
    path=Path(records[0]['path']);before=(path.stat().st_size,path.stat().st_mtime_ns)
    with h5py.File(path,'r') as h:
        data=pd.DataFrame({k:pd.Series(read_elem(h['obs'][k])).astype(str) for k in FIELDS})
    rows=[];checks=[];benchmarks=[]
    for record in records:
        selected=np.load(record['row_positions_ref'],allow_pickle=False)
        assert len(selected)==record['n_cells'] and len(np.unique(selected))==len(selected)
        study=record['dataset'][4:]
        full=np.flatnonzero(data.Study_name.eq(study).to_numpy())
        assert len(full)>0 and np.isin(selected,full).all()
        retained=data.iloc[selected].reset_index(drop=True)
        original=old_capture(retained)
        assert all(original[k]==old_scope[record['dataset']][k] for k in original)
        benchmarks.append(dict(source=record['dataset'],old5_capture_numbers_exact=True))
        for population,positions in [('contract_retained',selected),('all_available_same_study_in_cache',full)]:
            f=data.iloc[positions].reset_index(drop=True)
            row=dict(source=record['dataset'],role=record['role'],population=population,**summarize(f))
            row['source_overlap_quarantine_preserved']=record['role']=='overlap_quarantine'
            row['closed_Bassez_family_not_reopened']=record['dataset']=='ICB_Bassez'
            rows.append(row)
            checks.append(dict(source=record['dataset'],population=population,selected_metadata_rows=len(f),row_selection_sha256=hashlib.sha256(positions.astype('<i8').tobytes()).hexdigest(),all_selected_metadata_fields_sha256=hashlib.sha256(pd.util.hash_pandas_object(f,index=False).to_numpy().tobytes()).hexdigest(),public_primary_PMID_literals=sorted(f.PMID.unique().tolist()),current_count_gate_vs_groupby_checked=True))
        save('private_progress.json',dict(stage='SOURCE_UNIT_ALIAS_AND_RETAINED_SCOPE',completed_sources=len(benchmarks),total_sources=8,observed_UTC=datetime.now(timezone.utc).isoformat()))
    assert before==(path.stat().st_size,path.stat().st_mtime_ns)
    pd.DataFrame(rows).to_csv(out/'private_complete16_source_population_unit_and_capture_scope.tsv',sep='\t',index=False)
    save('private_current_source_metadata_and_selector_proofs.json',dict(source_contract_sha256=sha(contract),current_source_stat_unchanged=True,source_metadata_fields_only=True,checks=checks,all8_original5_scope_numbers_reproduced=benchmarks))
    t=dict(status='COMPLETE_ALL8_CURRENT_METADATA_ORIGINAL_SAMPLE_ALIASES_AND_BOTH_RETAINED_PARENT_SCOPES',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),source_positions=8,source_population_positions=16,old_capture_numbers_reproduced=40,no_new_RNA_or_effects=True,no_patient_identity_or_clinical_endpoint_promotion=True,quarantines_and_closed_families_preserved=True,scientific_outputs_public=False,files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(v['bytes'] for v in t['files'].values())<1024**2
    save('private_terminal.json',t);print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':
    main()

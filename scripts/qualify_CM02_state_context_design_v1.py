"""Outcome-free HNSCC state-count and exact-label opportunities; no RNA fit."""
from datetime import datetime, timezone
from pathlib import Path
from itertools import product
import hashlib
import json
import os
import shutil
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem

REG = ['4-1BB+ Treg','HLA-DR+ Treg','LAG3+ Treg','Prolif Treg']
OTHER = ['Tn/Tmem','IFN-stim CD4T','GZMK+T','Tmem','Teff','Tstr','Tprolif']
MISSING = {'','nan','na','none','unknown','n/a','__source_value_missing__'}


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def qualify(f):
    counts = f[REG+OTHER].to_numpy(float)
    assert np.isfinite(counts).all() and (counts >= 0).all()
    np.testing.assert_array_equal(counts, np.floor(counts))
    parent = np.column_stack([counts[:,:4].sum(axis=1),counts[:,4:].sum(axis=1)])
    tx = f.neo_type.isin(['D','D+T']).to_numpy()
    hpv = f['HPV/p16'].isin(['positive','negative']).to_numpy()
    unique = ~f.patient_registry.duplicated(keep=False).to_numpy()
    eligible = (parent >= 20).all(axis=1) & tx & hpv & unique
    records = []
    for i,j in product(range(4),range(7)):
        a = counts[eligible,i]/parent[eligible,0]
        b = counts[eligible,4+j]/parent[eligible,1]
        status = ('COUNT_SCOPE_READY_NOT_STATE_OR_CLINICAL_VALIDATION' if eligible.sum()>=20 and np.ptp(a)>0 and np.ptp(b)>0 else 'INSUFFICIENT_COVERAGE_OR_VARIATION')
        records.append(dict(reg_author_state=REG[i],other_author_T_state=OTHER[j],eligible_author_patients=int(eligible.sum()),reg_positive_patients=int((counts[eligible,i]>0).sum()),other_positive_patients=int((counts[eligible,4+j]>0).sum()),reg_unique_fractions=len(np.unique(a)),other_unique_fractions=len(np.unique(b)),zero_child_counts_retained=True,status=status,other_T_not_automatically_CD8_or_cytotoxic=True,exact_CM02_per_cell_membership_not_established=True))
    summary = dict(author_materials=len(f),unique_author_patients=int(f.patient_registry.nunique()),duplicate_phase_patients=int(f.loc[~unique,'patient_registry'].nunique()),both_disjoint_parents_atleast20=int((parent>=20).all(axis=1).sum()),missing_or_unsupported_treatment=int((~tx).sum()),missing_or_unsupported_HPV=int((~hpv).sum()),eligible_for_count_structure=int(eligible.sum()),by_treatment=[dict(arm=arm,author_patients=int((eligible&f.neo_type.eq(arm).to_numpy()).sum())) for arm in ['D','D+T']],clinical_outcome_values_not_read=True)
    return records,summary


def fixture():
    f = pd.DataFrame({k:[0,10,20] for k in REG+OTHER})
    f['patient_registry']=['P1','P2','P3'];f['neo_type']=['D','D','D+T'];f['HPV/p16']='positive'
    rows,s = qualify(f)
    assert len(rows)==28 and s['both_disjoint_parents_atleast20']==2
    assert all(r['status']=='INSUFFICIENT_COVERAGE_OR_VARIATION' for r in rows)
    f.loc[1,REG[0]]=0
    rows,s = qualify(f);assert s['eligible_for_count_structure']==2 and rows[0]['reg_positive_patients']==1
    return dict(zero_child_retained_when_parent_present=True,absent_parent_not_negative=True,constant_and_small_scope_not_promoted=True)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent==Path('/dev/shm') and not out.exists()
    assert shutil.disk_usage('/dev/shm').free>512*1024**2
    old=root/'HNSCC_author_state_count_layer_preflight_20261016_v1'
    table=root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv'
    plan=json.loads((old/'plan_before_any_author_state_clinical_effect.json').read_text())
    assert sha(table)==plan['input_sha256']['author_join']
    atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1'
    contract=atlas/'config/source_contract_v2_private.json'
    out.mkdir(mode=0o700);save=lambda name,value:(out/name).write_text(json.dumps(value,indent=2)+'\n')
    save('plan_before_new_state_opportunities.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),author_join_SHA256=sha(table),prior_plan_SHA256=sha(old/'plan_before_any_author_state_clinical_effect.json'),source_contract_SHA256=sha(contract),module='2026-10-14 immune atlas v1 T_whole_TME_primary_tumor CM02',source='GSE286827',primary_phase='pre',secondary_phase='post_not_independent_replication',all_count_opportunity_positions=56,reg_author_labels=REG,other_author_T_labels=OTHER,clinical_outcome_columns_not_read=True,no_RNA_or_new_association_tests=True,other_T_labels_not_promoted_to_CD8=True))
    save('private_synthetic_fixture.json',fixture())
    cols=['library_id','patient_registry','timepoint_registry','neo_type','HPV/p16']+REG+OTHER
    f=pd.read_csv(table,sep='\t',usecols=cols,keep_default_na=False,dtype={'library_id':str,'patient_registry':str})
    assert f.library_id.is_unique
    old_counts=pd.read_csv(old/'private_all57_original_author_counts.tsv',sep='\t',keep_default_na=False,dtype={'library_id':str,'patient_registry':str})
    assert set(f.library_id)==set(old_counts.library_id)
    for field in REG+OTHER:
        np.testing.assert_array_equal(pd.to_numeric(f.set_index('library_id')[field].sort_index()),pd.to_numeric(old_counts.set_index('library_id')[field].sort_index()))
    rows=[];coverage=[]
    for phase in ['pre','post']:
        g=f[f.timepoint_registry.eq(phase)].reset_index(drop=True)
        pairs,summary=qualify(g);rows.extend([dict(source='GSE286827',phase=phase,**r) for r in pairs]);coverage.append(dict(source='GSE286827',phase=phase,**summary))
    assert len(rows)==56
    z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r')
    source=np.asarray(read_elem(z['obs/dataset'])).astype(str)
    author=np.asarray(read_elem(z['obs/source_author_label'])).astype(str)[source=='GSE286827']
    nonmissing=sum(x.strip().casefold() not in MISSING for x in author)
    mapping=dict(source='GSE286827',current_author_label_cells=len(author),nonmissing_author_label_cells=int(nonmissing),complete_original_author_fine_state_mapping_accepted=False,sample_state_counts_cannot_reverse_assign_cell_labels=True,original_inputs_read_only=True)
    save('private_all56_state_pair_opportunities.json',rows);save('private_all2_phase_scope.json',coverage);save('private_exact_state_mapping_scope.json',mapping)
    files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()}
    assert sum(x['bytes'] for x in files.values())<1024**2
    terminal=dict(status='COMPLETE_CM02_CONTEXT_AND_ALL_AUTHOR_STATE_COUNT_OPPORTUNITIES_NO_NEW_EFFECTS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),all56_positions_preserved=True,prior_author_counts_reproduced=True,clinical_outcome_values_expression_weights_or_individual_rows_exported=False,files=files)
    save('private_terminal.json',terminal);print(json.dumps(dict(status=terminal['status'],terminal_SHA256=sha(out/'private_terminal.json'))))


if __name__=='__main__': main()

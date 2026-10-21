"""All frozen observed-marker support opportunities; no spatial effects or P.

Vector masks are checked against an independent scalar truth-table Counter on
every measured cell observation. Case maps and all individual arrays stay RAM.
"""
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
import hashlib
import json
import math
import os

import numpy as np
import pandas as pd
import xenium_observed_state_context_v1 as context_module

MODES = list(product([20,50], ['dual_regulatory_evidence','costate_regulatory_sensitivity']))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def scalar_counts(data, minimum, parent):
    ix={g:i for i,g in enumerate(context_module.REQUIRED)}
    counts=Counter(); rp=Counter(); ep=Counter()
    area=data['coordinates'].cell_area.to_numpy(float)
    for row,h,a in zip(data['counts'],data['human_gene_molecules'],area):
        if h<minimum or not math.isfinite(a) or a<=0:continue
        counts['qc']+=1
        t=any(row[ix[g]]>0 for g in ['CD3D','CD3E','TRAC'])
        counts['T_support']+=int(t)
        reg=t and row[ix['FOXP3']]>0 and (row[ix['IL2RA']]>0 if parent=='dual_regulatory_evidence' else (row[ix['IL2RA']]>0 or row[ix['CTLA4']]>0))
        cd8=t and row[ix['CD8A']]>0
        if reg and cd8:
            counts['ambiguous']+=1;reg=cd8=False
        cycle=row[ix['MKI67']]>0 and (row[ix['TOP2A']]>0 or row[ix['UBE2C']]>0)
        if reg:
            counts['regulatory']+=1
            flags=[row[ix['TNFRSF9']]>0,row[ix['LAG3']]>0,cycle]
            rp[sum((1<<j) for j,b in enumerate(flags) if b)]+=1
            for j,b in enumerate(flags):counts['R'+str(j)]+=int(b)
        if cd8:
            counts['CD8']+=1
            flags=[row[ix['PDCD1']]>0 and (row[ix['TOX']]>0 or row[ix['HAVCR2']]>0),cycle,row[ix['TCF7']]>0 or row[ix['IL7R']]>0]
            ep[sum((1<<j) for j,b in enumerate(flags) if b)]+=1
            for j,b in enumerate(flags):counts['E'+str(j)]+=int(b)
        counts['B_cycling']+=int(not t and row[ix['MS4A1']]>0 and row[ix['CD79A']]>0 and cycle)
        counts['epithelial_cycling']+=int(not t and row[ix['EPCAM']]>0 and cycle)
    return counts,[rp[k] for k in range(8)],[ep[k] for k in range(8)]


def measure_one(context, ordinal):
    metadata=context['pairs'][ordinal]['source_metadata']
    if metadata['diagnosis']!='HNSCC':return dict(ordinal=ordinal,metadata=metadata,excluded=True,modes=[])
    d=context_module.load_pair(context,ordinal);results=[]
    for minimum,parent in MODES:
        masks=context_module.observed_masks(d,minimum,parent)
        counts={k:int(masks[k].sum()) for k in ['qc','T_support','regulatory','CD8','ambiguous','B_cycling','epithelial_cycling']}
        rs=masks['R_states'].sum(axis=0).astype(int).tolist();es=masks['E_states'].sum(axis=0).astype(int).tolist()
        rp=np.bincount(masks['R_states'][masks['regulatory']].astype(np.int64)@np.array([1,2,4]),minlength=8).tolist()
        ep=np.bincount(masks['E_states'][masks['CD8']].astype(np.int64)@np.array([1,2,4]),minlength=8).tolist()
        independent,irp,iep=scalar_counts(d,minimum,parent)
        assert all(independent[k]==v for k,v in counts.items())
        assert rs==[independent['R'+str(k)] for k in range(3)] and es==[independent['E'+str(k)] for k in range(3)]
        assert rp==irp and ep==iep and sum(rp)==counts['regulatory'] and sum(ep)==counts['CD8']
        assert not (masks['regulatory']&masks['CD8']).any()
        statuses=[]
        for i,j in product(range(3),range(3)):
            nr,ne=counts['regulatory'],counts['CD8']
            status='OBSERVED_STATE_COUNT_SCOPE_READY'
            if nr<20 or ne<20:status='INSUFFICIENT_PARENT_OBSERVATIONS'
            elif rs[i]<5 or es[j]<5:status='INSUFFICIENT_OBSERVED_CHILD_SUPPORT'
            elif rs[i]==nr or es[j]==ne:status='NO_WITHIN_PARENT_STATE_VARIATION'
            statuses.append(dict(R_state=context_module.R_NAMES[i],E_state=context_module.E_NAMES[j],status=status))
        results.append(dict(minimum_human_gene_molecules=minimum,regulatory_parent_definition=parent,counts=counts,R_observed_state_counts=rs,E_observed_state_counts=es,R_joint_state_patterns=rp,E_joint_state_patterns=ep,opportunities=statuses,all_cell_observations_independent_scalar_checked=True))
    return dict(ordinal=ordinal,metadata=metadata,excluded=False,matrix_cells=d['matrix_cells'],canonical_human_features=d['canonical_human_features'],raw_CSC_marker_columns_and_human_totals_independent_bincount_checked=d['raw_CSC_marker_columns_and_human_totals_independent_bincount_checked'],modes=results)


def synthetic_fixture():
    x=np.zeros((7,len(context_module.REQUIRED)),np.int64);ix={g:i for i,g in enumerate(context_module.REQUIRED)}
    rows=[['CD3E','FOXP3','IL2RA','TNFRSF9','MKI67','TOP2A'],['CD3D','CD8A','PDCD1','TOX','MKI67','UBE2C'],['TRAC','FOXP3','IL2RA','CD8A'],['CD3E','FOXP3','CTLA4','LAG3'],['CD3E','FOXP3','IL2RA'],['MS4A1','CD79A','MKI67','TOP2A'],['EPCAM','MKI67','UBE2C']]
    for i,genes in enumerate(rows):
        for g in genes:x[i,ix[g]]=1
    d=dict(counts=x,human_gene_molecules=np.array([20,20,20,20,19,50,50]),coordinates=pd.DataFrame({'cell_area':[10.]*7}))
    primary=context_module.observed_masks(d,20,'dual_regulatory_evidence');sensitivity=context_module.observed_masks(d,20,'costate_regulatory_sensitivity');high=context_module.observed_masks(d,50,'dual_regulatory_evidence')
    assert primary['regulatory'].sum()==primary['CD8'].sum()==primary['ambiguous'].sum()==1
    assert primary['R_states'].sum(axis=0).tolist()==[1,0,1] and primary['E_states'].sum(axis=0).tolist()==[1,1,0]
    assert sensitivity['regulatory'].sum()==2 and sensitivity['R_states'].sum(axis=0).tolist()==[1,1,1]
    assert high['regulatory'].sum()==high['CD8'].sum()==high['ambiguous'].sum()==0
    assert high['B_cycling'].sum()==high['epithelial_cycling'].sum()==1
    return dict(overlapping_lineage_unresolved=True,zero_children_retained=True,parent_definition_sensitivity_distinct=True,QC_threshold_exact=True,background_control_counts_exact=True)


def main():
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert out.parent.parent==Path('/dev/shm') and not out.exists();out.mkdir(mode=0o700)
    save=lambda n,v:(out/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    save('plan_before_observed_T_state_support.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),context_source_SHA256=sha(context_module.__file__),prior_private_question_SHA256=os.environ['PANCANCER_PRIOR_QUESTION_SHA256'],all4_QC_parent_modes=[dict(minimum=minimum,parent=parent) for minimum,parent in MODES],all9_pairs_each_mode=True,primary_mode=MODES[0],minimum_parent_observations=20,minimum_observed_child=5,within_parent_state_variation_required=True,no_new_spatial_effects_P_or_clinical_tests=True,definitions='All known-marker recipes from the prior question; molecular counts not formal identity or function.',all_cell_observations_checked_by_independent_scalar_Counter=True,all_payloads_cell_and_case_arrays_RAM_only=True,source_anatomical_site_recurrence_and_ICI_timing_distinct=True))
    save('private_synthetic_fixture.json',synthetic_fixture())
    context=context_module.load_context();rows=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        work={pool.submit(measure_one,context,i):i for i in range(18)}
        for future in as_completed(work):
            rows.append(future.result())
            temp=out/'private_checkpoint.tmp';temp.write_text(json.dumps(dict(completed_source_slots=len(rows)),indent=2)+'\n');temp.replace(out/'private_checkpoint.json')
    rows.sort(key=lambda r:r['ordinal']);assert len(rows)==18 and sum(not r['excluded'] for r in rows)==17
    sections=[r for r in rows if not r['excluded']];cases={r['metadata']['case'] for r in sections};assert len(cases)==10
    technical=[{k:v for k,v in r.items() if k!='metadata'} for r in rows]
    save('private_all18_source_slot_marker_support.json',technical)
    opportunities=[];source_counts=[]
    for mi,(minimum,parent) in enumerate(MODES):
        pooled=Counter();r_counts=np.zeros(3,np.int64);e_counts=np.zeros(3,np.int64);rp=np.zeros(8,np.int64);ep=np.zeros(8,np.int64)
        for r in sections:
            m=r['modes'][mi];pooled.update(m['counts']);r_counts+=m['R_observed_state_counts'];e_counts+=m['E_observed_state_counts'];rp+=m['R_joint_state_patterns'];ep+=m['E_joint_state_patterns']
        source_counts.append(dict(minimum_human_gene_molecules=minimum,regulatory_parent_definition=parent,cell_observation_slots_across_sections_not_patient_prevalence=dict(pooled),R_observed_counts=r_counts.tolist(),E_observed_counts=e_counts.tolist(),R_patterns=rp.tolist(),E_patterns=ep.tolist()))
        for pi in range(9):
            eligible=[r for r in sections if r['modes'][mi]['opportunities'][pi]['status']=='OBSERVED_STATE_COUNT_SCOPE_READY']
            statuses=Counter(r['modes'][mi]['opportunities'][pi]['status'] for r in sections)
            item=sections[0]['modes'][mi]['opportunities'][pi]
            opportunities.append(dict(source='GSE300147',minimum_human_gene_molecules=minimum,regulatory_parent_definition=parent,R_observed_state=item['R_state'],E_observed_state=item['E_state'],source_cases=len(cases),HNSCC_sections=len(sections),supporting_sections=len(eligible),supporting_source_cases=len({r['metadata']['case'] for r in eligible}),supporting_primary_source_cases=len({r['metadata']['case'] for r in eligible if r['metadata']['recurrence']=='primary'}),supporting_recurrent_source_cases=len({r['metadata']['case'] for r in eligible if r['metadata']['recurrence']=='recurrent'}),section_status_counts=dict(statuses),no_spatial_effect_or_identity_validation=True))
    sites={s:len({r['metadata']['case'] for r in sections if r['metadata']['anatomical_site']==s}) for s in sorted({r['metadata']['anatomical_site'] for r in sections})}
    save('private_all36_source_state_support_opportunities.json',opportunities);save('private_all4_source_marker_observation_counts.json',source_counts)
    save('private_source_anatomical_scope.json',dict(source='GSE300147',HNSCC_source_cases=10,HNSCC_sections=17,excluded_non_HNSCC_source_slots=1,case_counts_by_anatomical_site=sites,cases_may_overlap_site_categories=True,source_recurrence_not_anatomical_class_or_ICI_timing=True,no_raw_source_case_or_sample_rows_saved=True))
    terminal=dict(status='COMPLETE_ALL36_SOURCE_OBSERVED_STATE_SUPPORT_OPPORTUNITIES_WITH_SCALAR_CHECKS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),context_source_SHA256=sha(context_module.__file__),joint_audit_SHA256=context['joint_audit_SHA256'],unit_audit_SHA256=context['unit_audit_SHA256'],all18_source_slots_and36_source_opportunities_preserved=True,all17_measured_sections_all_cells_four_modes_scalar_checked=True,no_cell_case_arrays_or_payloads_saved=True,no_new_spatial_or_clinical_effects=True,files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()})
    save('private_terminal.json',terminal);print(json.dumps(dict(status=terminal['status'],terminal_SHA256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

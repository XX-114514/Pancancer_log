"""Freeze all native ligand-programme tests before any association or P.

Reads only the eight covariates from the existing clinical source, not its
outcome or clinical eligibility field. Raw measurements remain private.
"""
import json
import os
import time
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import load,sha,VIEWS
from NSCLC_native_MFA_clinical_inputs_v1 import COMMON
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_external_programme_coordination_math_v2 import prepare


def design_gate(p):
    s=np.linalg.svd(p['design'],compute_uv=False)
    machine=int((s>s[0]*np.finfo(float).eps*max(p['design'].shape)).sum())
    if p['design_rank']!=machine or s[p['design_rank']-1]/s[0]<1e-8:
        raise ValueError('UNSTABLE_FULL_NUISANCE_SPAN')


def read_controls(data):
    root=data['root'];src=root/'NSCLC_native_MFA_clinical_coverage_20261016_v1'
    terminal=src/'private_terminal.json'
    assert sha(terminal)=='76bf844a948caab15a05b2eb5a67a3a0e4e1ac04cf92e193ff981566e543f224'
    records=json.loads(terminal.read_text());path=src/'private_all208_clinical_eligibility.tsv'
    assert sha(path)==records['files'][path.name]['sha256']
    c=pd.read_csv(path,sep='\t',usecols=['sampleID']+COMMON).set_index('sampleID')
    assert c.index.is_unique and set(c.index)==set(data['cases'].sampleID)
    c=c.loc[data['cases'].sampleID,COMMON].apply(pd.to_numeric,errors='coerce')
    return c.to_numpy(float),sha(path)


def make_case_partitions(data,clinical):
    result={};train=data['train'];blocks=data['cases'].holdout_block.to_numpy(str)
    for sv,sender in enumerate(VIEWS):
        for rv,receiver in enumerate(VIEWS):
            if sv==rv:continue
            quality=np.column_stack([clinical,data['prepared'][sv]['quality'],data['prepared'][rv]['quality']])
            complete=np.isfinite(quality).all(axis=1)
            for part,mask in [('internal_holdout',~train),('discovery',train),('all_frozen_cases',np.ones(len(train),bool))]:
                take=np.flatnonzero(mask&complete);counts=pd.Series(blocks[take]).value_counts()
                take=take[np.array([counts[b]>=2 for b in blocks[take]],bool)]
                result[(sender,receiver,part)]=(take,quality[take],blocks[take])
    return result


def qualify(data,z,slots,interfaces,clinical):
    partitions=make_case_partitions(data,clinical);views={v:i for i,v in enumerate(VIEWS)}
    receptor={r:i for i,r in enumerate(z['receptors'])};by_slot=slots.set_index('slot_index');rows=[];edges=[]
    for interface in interfaces.itertuples():
        slot=by_slot.loc[interface.slot_index];sv=views[interface.sender_view];rv=views[interface.receiver_view]
        record='NATIVE_PRIOR_INTERFACE|%04d'%interface.interface_index
        item=dict(record=record,interface_index=interface.interface_index,slot_index=interface.slot_index,
                  sender_view=interface.sender_view,receiver_view=interface.receiver_view,factor=interface.factor,
                  direction=int(interface.direction),prospective_position=int(interface.prospective_position),
                  ligand=interface.ligand,ligand_axis_index=int(interface.ligand_axis_index),status='ELIGIBLE_FOR_FIXED_PRIMARY_TEST',
                  active_cases=0,discovery_cases=0,nuisance_rank=0,residual_df=0,qualifying_reference_receptors='',test_index=-1)
        reasons=[]
        hp,hq,hb=partitions[(interface.sender_view,interface.receiver_view,'internal_holdout')]
        dp,dq,db=partitions[(interface.sender_view,interface.receiver_view,'discovery')]
        item.update(active_cases=len(hp),discovery_cases=len(dp))
        if interface.status!='RNA_INPUT_READY_EFFECT_MODEL_NOT_RUN':reasons.append('INPUT:'+interface.status)
        compatible=[]
        if not reasons:
            j=interface.ligand_axis_index
            if min(int((z['ligand_positive'][sv,hp,j]>0).sum()),int((z['ligand_positive'][sv,dp,j]>0).sum()))<10:
                reasons.append('SENDER_POSITIVE_IN_FEWER_THAN10_CASES_IN_EITHER_PARTITION')
            for rec in interface.receiver_RNA_compatible_receptors.split('|'):
                k=receptor[rec];hn=int((z['receptor_positive'][rv,hp,k]>0).sum());dn=int((z['receptor_positive'][rv,dp,k]>0).sum())
                ok=min(hn,dn)>=10
                edges.append(dict(interface_index=interface.interface_index,receptor=rec,holdout_joint_positive_cases=hn,discovery_joint_positive_cases=dn,both_partition10case_gate=ok))
                if ok:compatible.append(rec)
            if not compatible:reasons.append('NO_FIXED_RECEPTOR_WITH10_CASES_IN_BOTH_PARTITIONS')
            item['qualifying_reference_receptors']='|'.join(compatible)
            for name,indices,quality,blocks in [('primary',hp,hq,hb),('discovery',dp,dq,db)]:
                try:
                    p=prepare(z['ligand_mean'][sv,indices,j],int(interface.direction)*z['response'][interface.slot_index,indices],quality,blocks,min_n=50,min_df=10)
                    design_gate(p)
                    if name=='primary':item.update(nuisance_rank=p['design_rank'],residual_df=p['residual_df'])
                except ValueError as error:reasons.append(name.upper()+':'+str(error))
        if reasons:item['status']=';'.join(reasons)
        rows.append(item)
    counter=0
    for item in rows:
        if item['status']=='ELIGIBLE_FOR_FIXED_PRIMARY_TEST':item['test_index']=counter;counter+=1
    return pd.DataFrame(rows),pd.DataFrame(edges),partitions


def main():
    data=load();root=data['root'];budget=resource(root,32*1024**2)
    out=root/'NSCLC_continuous_prior_native_association_inputs_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    src=root/'NSCLC_continuous_prior_native_interfaces_20261016_v2';proof=src/'private_independent_native_interface_audit.json'
    receipt=json.loads((root/'HNSCC_receiver_followup_job_control_20261016_v2/audit_NSCLC_continuous_prior_native_interfaces_v2/launch_receipt.json').read_text())
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
              scope='All4800 fixed ligand-receiver positions; native63 internal-holdout primary,145 and208 effect-only descriptions. Not untouched external or reviewed-independent-patient validation.',
              controls=COMMON,quality='All original sender and receiver quality variables plus intercept and complete author histology-agent block dummies; rank all numeric columns, retain entire nuisance span.',
              primary='Partial rank correlation of sender mean-logCP10K ligand and prior orientation times LR-excluded frozen receiver programme. No clinical outcome loaded.',
              gates='Input RNA/fidelity readiness;>=10 actual sender-positive cases AND>=1 fixed receptor jointly positive in>=10 cases in BOTH original partitions. Primary/discovery >=50 informative cases,df>=10,full numerical span stable and both residuals variable. No per-case expression filtering.',
              nuisance_rank_relative_cutoff=1e-10,minimum_retained_relative_singular=1e-8,
              permutations=199999,bootstrap=2000,min_valid_bootstrap=1900,
              inferential_family='All pre-effect eligible primary records, all4800 planned slots retained. BH,BY,simultaneous99pct MonteCarlo-upper-BH. Two-sided plus-one P; no after-result precision increases.',
              permutations_method='Within author category blocks Freedman-Lane residual permutation; approximate conditional exchangeability, not universal finite-sample guarantee.',
              bootstrap_method='Within original category blocks case bootstrap, rerank and refit all original controls each draw; same fixed streams and explicit invalids. Case identity not asserted to be independent patient.',
              seed_namespace='native_continuous_prior_association_20261016_v1',
              sensitivities='Discovery/all frozen effect only; original unmasked receiver response; add sender own same factor; add both full-native-immune cell fractions; add complete both-pool author subtype proportions. No extra P or rescue.',
              direct_binding_causality_or_clinical_validation=False,associations_or_P_computed=False)
    (out/'plan_before_any_native_interface_association.json').write_text(json.dumps(plan,indent=2)+'\n')
    deadline=time.monotonic()+21600
    while not proof.exists():
        process=Path('/proc')/str(receipt['pid'])/'stat';assert process.exists(),'Bound input audit missing without proof'
        fields=process.read_text().rsplit(')',1)[1].split();assert fields[0]!='Z' and fields[19]==str(receipt['proc_start_ticks'])
        assert time.monotonic()<deadline,'Observation expired; no termination inferred';time.sleep(5)
    accepted=json.loads(proof.read_text());assert accepted['status'].startswith('PASS_ALL4800_NATIVE_INTERFACE_INPUTS_')
    assert accepted['producer_terminal_sha256']==sha(src/'private_terminal.json')
    terminal=json.loads((src/'private_terminal.json').read_text())
    for name,item in terminal['files'].items():assert sha(src/name)==item['sha256']
    clinical,clinical_sha=read_controls(data)
    z=np.load(src/'private_all_native_interface_measurements.npz',allow_pickle=False)
    assert np.array_equal(z['cases'],data['cases'].sampleID.to_numpy(str)) and np.array_equal(z['train'],data['train'])
    slots=pd.read_csv(src/'private_all1600_original_slots.tsv.gz',sep='\t',keep_default_na=False)
    interfaces=pd.read_csv(src/'private_all4800_native_interface_readiness.tsv.gz',sep='\t',keep_default_na=False)
    scope,edges,parts=qualify(data,z,slots,interfaces,clinical)
    assert len(scope)==4800 and scope.record.is_unique
    scope.to_csv(out/'private_all4800_pre_effect_eligibility.tsv.gz',sep='\t',index=False,compression='gzip')
    edges.to_csv(out/'private_all_fixed_receptor_case_coverage.tsv.gz',sep='\t',index=False,compression='gzip')
    np.savez_compressed(out/'private_all208_covariates_only.npz',clinical=clinical,clinical_names=np.array(COMMON),cases=z['cases'],train=z['train'])
    partition_rows=[]
    for (sender,receiver,partition),(indices,quality,blocks) in parts.items():
        partition_rows.append(dict(sender_view=sender,receiver_view=receiver,partition=partition,informative_cases=len(indices),category_blocks=len(set(blocks))))
    pd.DataFrame(partition_rows).to_csv(out/'private_all36_covariate_partition_coverage.tsv',sep='\t',index=False)
    frozen=dict(input_audit_sha256=sha(proof),clinical_covariate_source_sha256=clinical_sha,eligible_primary_records=int(scope.test_index.ge(0).sum()),all_planned_positions=4800,
                no_response_clinical_eligibility_or_pathology_columns_loaded=True,associations_or_P_calculated=False)
    (out/'private_frozen_family_before_effects.json').write_text(json.dumps(frozen,indent=2)+'\n')
    terminal=dict(status='COMPLETE_ALL4800_NATIVE_INTERFACE_ASSOCIATION_ELIGIBILITY_FROZEN_NO_EFFECTS_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
                  files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},associations_or_P_calculated=False)
    assert sum(v['bytes'] for v in terminal['files'].values())<budget['output_cap_bytes']
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

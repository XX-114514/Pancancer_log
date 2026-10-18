"""Independent complete scalar ULM/OLS and same-cell aggregate replay; no P."""
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from NSCLC_native_MFA_interpretation_inputs_v1 import load, pool_native, sha, VIEWS
from test_NSCLC_native_interface_association_v1 import inputs, context
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from describe_NSCLC_native_interface_mechanistic_followup_v1 import read_scope, exclusions, OUT, SCOPE_SHA


def independent_ols(x, y, q, blocks):
    q = rankdata(q, axis=0)
    sd = q.std(axis=0)
    q = (q - q.mean(axis=0)) / np.where(sd > 0, sd, 1)
    design = np.column_stack([np.ones(len(x)), *[(blocks == b).astype(float) for b in sorted(set(blocks))], q])
    s = np.linalg.svd(design, compute_uv=False)
    rank = int((s > s[0]*1e-10).sum())
    if len(x) < 50 or len(x)-rank < 10:
        return np.nan
    assert rank == int((s > s[0]*np.finfo(float).eps*max(design.shape)).sum())
    assert s[rank-1]/s[0] >= 1e-8
    response = np.column_stack([rankdata(x), rankdata(y)])
    residual = response - design @ np.linalg.lstsq(design, response, rcond=1e-10)[0]
    norms = np.linalg.norm(residual, axis=0)
    if norms.min() < 1e-10: return np.nan
    return float(residual[:, 0] @ residual[:, 1] / np.prod(norms))


def equal(expected, actual, tolerance=1e-9):
    if np.isnan(expected): assert pd.isna(actual)
    else: assert np.isfinite(actual) and abs(expected-float(actual)) < tolerance, (expected, actual)


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); resource(root, 1024**2)
    job = root / 'HNSCC_receiver_followup_job_control_20261016_v2/describe_NSCLC_native_interface_mechanistic_followup_v1'
    receipt = json.loads((job/'launch_receipt.json').read_text())
    out = root / OUT; terminal_path = out/'private_terminal.json'
    deadline = time.monotonic()+3600
    while not terminal_path.exists():
        proc = Path('/proc')/str(receipt['pid'])/'stat'
        assert proc.exists(), 'Bound producer absent without terminal; no restart inferred'
        fields = proc.read_text().rsplit(')', 1)[1].split()
        assert fields[0] != 'Z' and fields[19] == str(receipt['proc_start_ticks'])
        assert time.monotonic() < deadline, 'Observation deadline, not scientific failure'
        time.sleep(3)
    terminal = json.loads(terminal_path.read_text())
    assert terminal['code_sha256'] == receipt['source_sha256'] and terminal['scope_sha256'] == SCOPE_SHA
    for name, item in terminal['files'].items(): assert sha(out/name) == item['sha256']
    data = load(); scope, candidates = read_scope(data)
    _, _, _, z, clinical, _ = inputs(root); state = context(data, z, clinical)
    tables = {key: pd.read_csv(out/name, sep='\t', keep_default_na=False, na_values=['']) for key, name in [
        ('functional', 'private_all_functional_descriptions.tsv'),
        ('control', 'private_all_original_programme_control_descriptions.tsv'),
        ('receptor', 'private_all_canonical_receptor_support.tsv'),
        ('sender', 'private_all_sender_subtype_contributions.tsv')]}
    assert not tables['functional'].duplicated(['record','partition','function','readout','control']).any()
    assert not tables['control'].duplicated(['record','partition','control']).any()
    assert len(tables['functional']) == scope['reference_slots']*2*3*3
    assert len(tables['control']) == len(candidates)*2*5
    packed_root = root/'NSCLC_continuous_prior_native_interfaces_20261016_v2'
    panel = pd.read_csv(packed_root/'private_composite_panel_axis.tsv', sep='\t').gene.tolist()
    bits = np.unpackbits(np.load(packed_root/'private_positive_panel_packed.npy'), axis=1, bitorder='little')[:, :len(panel)].astype(bool)
    cell_units = np.load(data['native']/'private_cell_unit_index.npy', mmap_mode='r')
    cell_pool = z['native_unit_to_case_pool'][cell_units]
    hits = np.load(data['native']/'positive_cell_counts.npy', mmap_mode='r')
    genes = pd.Index(data['genes']); pooled, _ = pool_native(data); omit = exclusions(root, scope)
    max_error = 0.; n_function = n_control = n_receptor = n_sender = 0
    for item in candidates.itertuples():
        sv, rv = VIEWS.index(item.sender_view), VIEWS.index(item.receiver_view)
        li = int(item.ligand_axis_index); k = int(item.factor.rsplit('F',1)[1])-1
        p = data['prepared'][sv]; w = data['model']['weights_'+str(sv)]; tau = data['model']['precision_'+str(sv)]
        # Independent augmented system, rather than producer's normal equations.
        augmented = np.vstack([np.eye(10), np.sqrt(tau)[:,None]*w])
        rhs = np.vstack([np.zeros((10,2000)), np.diag(np.sqrt(tau))])
        a = np.linalg.lstsq(augmented, rhs, rcond=None)[0]
        retained = data['genes'][p['selected_native_columns']] != item.ligand
        sender = p['X'][:, retained] @ a[k, retained]
        remain = np.array([j for j in data['remaining'] if str(data['genes'][j]) not in omit[item.ligand]])
        for part in ['internal_holdout','discovery']:
            ix, q, blocks = state['parts'][(item.sender_view,item.receiver_view,part)]
            qs = dict(base=q, recipient_own_ligand=np.column_stack([q,z['ligand_mean'][rv,ix,li]]),
                sender_original_programme=np.column_stack([q,z['original_scores'][sv,ix,k]]),
                sender_ligand_removed_programme=np.column_stack([q,sender[ix]]),
                native_immune_fractions=np.column_stack([q,state['fractions'][[sv,rv]][:,ix].T]))
            x = z['ligand_mean'][sv,ix,li]
            for label, quality in qs.items():
                row = tables['control'].query('record == @item.record and partition == @part and control == @label')
                assert len(row)==1; row=row.iloc[0]
                value=independent_ols(x,item.direction*z['response'][item.slot_index,ix],quality,blocks)
                equal(value,row.rho)
                assert row.n_cases==len(ix) and row.sender_ligand_features_removed==int((~retained).sum())
                if np.isfinite(value):max_error=max(max_error,abs(value-row.rho));assert row.status=='ESTIMABLE'
                n_control+=1
            for complex_ in scope['ligand_specs'][item.ligand]['complexes']:
                rec='_'.join(complex_)
                selected=np.isin(z['native_unit_to_case_pool'],rv*208+ix)
                marginal={g:int(hits[selected,genes.get_loc(g)].sum()) for g in complex_}
                row=tables['receptor'].query('record == @item.record and partition == @part and receptor == @rec')
                assert len(row)==1;row=row.iloc[0]
                assert json.loads(row.marginal_counts_json)==marginal and row.joint_upper_bound==min(marginal.values())
                assert row.receiver_cells==int(z['cells'][rv,ix].sum())
                if set(complex_).issubset(panel):
                    joint=np.all(bits[:,[panel.index(g) for g in complex_]],axis=1)
                    counts=np.array([np.count_nonzero(joint & (cell_pool==rv*208+i)) for i in ix])
                    assert row.joint_positive_cells==int(counts.sum()) and row.positive_cases==int((counts>0).sum())
                    assert row.joint_status=='ACTUAL_SAME_CELL_RNA'
                else:assert pd.isna(row.joint_positive_cells) and row.joint_status=='COMPONENT_NOT_IN_ACCEPTED_PACKED_PANEL'
                n_receptor+=1
            expected_labels=set(data['units'].loc[np.isin(z['native_unit_to_case_pool'],sv*208+ix),'sub_cell_type'])
            source=tables['sender'].query('record == @item.record and partition == @part')
            assert set(source.author_subtype)==expected_labels and not source.author_subtype.duplicated().any()
            assert int(source.ligand_positive_cells.sum())==int(z['ligand_positive'][sv,ix,li].sum())
            assert int(source.cells.sum())==int(z['cells'][sv,ix].sum())
            for row in source.itertuples():
                take=np.isin(z['native_unit_to_case_pool'],sv*208+ix)&data['units'].sub_cell_type.eq(row.author_subtype).to_numpy()
                assert row.cells==int(data['units'].loc[take,'n_cells'].sum())
                assert row.ligand_positive_cells==int(hits[take,genes.get_loc(item.ligand)].sum())
            n_sender+=len(source)
        for function in scope['ligand_specs'][item.ligand]['functions']:
            target=np.isin(data['genes'][remain],list(data['refs'][function])); nt=int(target.sum());nb=int((~target).sum())
            xall=pooled[rv][:,remain]; qualified=nt>=20 and nb>=500
            if qualified:
                tm=xall[:,target].mean(axis=1);bm=xall[:,~target].mean(axis=1)
                center=xall-xall.mean(axis=1)[:,None]
                corr=center[:,target].sum(axis=1)/np.sqrt((nt-nt*nt/len(target))*np.square(center).sum(axis=1))
                ulm=corr*np.sqrt((len(target)-2)/((1-corr+1e-16)*(1+corr+1e-16)))
                values=dict(ULM_score=ulm,target_mean=tm,nontarget_mean=bm)
                detected=(xall[:,target]>0).sum(axis=1)
            for part in ['internal_holdout','discovery']:
                ix,q,blocks=state['parts'][(item.sender_view,item.receiver_view,part)]
                qs=dict(base=q,recipient_own_ligand=np.column_stack([q,z['ligand_mean'][rv,ix,li]]),
                    native_immune_fractions=np.column_stack([q,state['fractions'][[sv,rv]][:,ix].T]))
                for readout in ['ULM_score','target_mean','nontarget_mean']:
                    for label,quality in qs.items():
                        row=tables['functional'].query('record == @item.record and partition == @part and function == @function and readout == @readout and control == @label')
                        assert len(row)==1;row=row.iloc[0]
                        assert bool(row.qualified)==qualified and row.remaining_targets==nt and row.remaining_background==nb and row.n_cases==len(ix)
                        value=independent_ols(z['ligand_mean'][sv,ix,li],values[readout][ix],quality,blocks) if qualified else np.nan
                        equal(value,row.rho)
                        if qualified:assert bool(row.low_foreground_support)==bool(np.median(detected[ix])<10 or np.median(tm[ix])<=0)
                        if np.isfinite(value):max_error=max(max_error,abs(value-row.rho));assert row.status=='ESTIMABLE'
                        n_function+=1
    assert n_control==len(tables['control']) and n_function==len(tables['functional']) and n_receptor==len(tables['receptor']) and n_sender==len(tables['sender'])
    result=dict(status='PASS_COMPLETE_FIXED_FUNCTIONAL_ULM_OLS_CONTROLS_RECEPTOR_AND_SENDER_AGGREGATE_REPLAY',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        producer_terminal_sha256=sha(terminal_path),scope_sha256=SCOPE_SHA,
        rows_checked=dict(functional=n_function,control=n_control,receptor=n_receptor,sender=n_sender),
        maximum_effect_error=max_error,
        limits='Shared previously accepted native pooling and exclusion I/O; independent ULM scalar formula, augmented sender projection, full rank OLS effects and unpacked Boolean receptor counts. Not a new raw-cell RNA audit or biological validation.')
    path=out/'private_independent_mechanistic_description_audit.json';assert not path.exists()
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=result['status'],completed_UTC=result['completed_UTC'],proof_sha256=sha(path))))


if __name__=='__main__':main()

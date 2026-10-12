"""Full frozen-candidate coverage and expression/count/depth design diagnostics.

No association model or p-value. The nine-column diagnostic design omits
major fractions, sender subtype mixture and shared programmes; it is not
acceptance of the eventual inferential model or receptor-specific mechanism.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import gc
import numpy as np
import pandas as pd


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''):
            h.update(b)
    return h.hexdigest()


def main():
    start = datetime.now(timezone.utc).isoformat()
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    parent = root / 'immune_same_name_family_source_review_v1'
    full = parent / 'fullgene_receiver_target_readiness_v1'
    native = parent / 'supported_LR_native_cell_panel_v1'
    provenance = parent / 'patient_block_provenance_overlay_v1'
    out = parent / 'fullgene_receiver_response_preflight_v1'
    assert not out.exists()
    plan_path = root / 'fullgene_receiver_target_code_v1/fullgene_receiver_response_preflight_plan_v1.json'
    plan = json.loads(plan_path.read_text())
    assert sha(__file__) == plan['code_sha256'][Path(__file__).name]
    paths = dict(full_inventory=full / 'inventory.json', full_audit=full / 'independent_audit.json',
                 native_inventory=native / 'inventory.json', provenance_inventory=provenance / 'inventory.json',
                 provenance_audit=provenance / 'independent_audit.json',
                 membership=parent / 'membership_review_v2/source_node_analysis_membership.tsv', fixed_candidates=native / 'fixed336_reference_supported_candidates.tsv')
    for name, p in paths.items():
        assert sha(p) == plan['input_sha256'][name], name
    inventory = json.loads(paths['full_inventory'].read_text())
    for n, h in inventory['outputs'].items():
        assert sha(full / n) == h
    prv = json.loads(paths['provenance_inventory'].read_text())
    scope_path = provenance / 'private_all1290_patient_provenance_overlay.tsv'
    assert sha(scope_path) == prv['outputs'][scope_path.name]
    scope = pd.read_csv(scope_path, sep='\t', keep_default_na=False)
    assert len(scope) == 1290 and scope.sampleID.is_unique and scope.analysis_eligible.dtype==bool
    fixed = pd.read_csv(native / 'fixed336_reference_supported_candidates.tsv', sep='\t', keep_default_na=False)
    fixed['model_key'] = fixed[['context_index', 'program', 'sender_family', 'receiver_family', 'ligand', 'receptor']].astype(str).agg('|'.join, axis=1)
    priors = pd.read_csv(full / 'all336_prior_coverage.tsv', sep='\t', keep_default_na=False)
    assert len(fixed) == 336 and set(fixed.model_key) == set(priors.model_key)
    member = pd.read_csv(paths['membership'], sep='\t')
    catalogue = member.loc[member.analysis_family_main_membership, ['subCluster', 'family']].drop_duplicates()
    assert len(catalogue) == 13 and catalogue.subCluster.is_unique
    keycols = ['gse_id', 'stratum', 'sender_family', 'receiver_family', 'ligand', 'receptor']
    interfaces = fixed[keycols].drop_duplicates().sort_values(keycols).reset_index(drop=True)
    interfaces['interface_id'] = ['IF%04d' % (i + 1) for i in range(len(interfaces))]
    mapping = fixed[['model_key'] + keycols].merge(interfaces, on=keycols, validate='many_to_one')
    rows = []; design_names = ['L', 'R', 'log_sender_n', 'log_receiver_n', 'log_sender_RNA', 'log_receiver_RNA', 'log_sender_nGenes', 'log_receiver_nGenes']
    records = json.loads(paths['native_inventory'].read_text())['source_records']
    for rec in records:
        sid = rec['gse_id']
        assert sha(rec['selected_cache']) == rec['cache_sha256']
        z = dict(np.load(rec['selected_cache'], allow_pickle=False))
        b = dict(np.load(full / 'private_source_unit_arrays' / (sid + '.npz'), allow_pickle=False))
        samples = b['samples'].astype(str); genes = list(b['genes']); gi = {g: i for i, g in enumerate(genes)}
        nodes = list(b['nodes']); groups = len(nodes); ns = len(samples)
        assert np.array_equal(samples, z['samples'])
        cell_family = z['families'][z['family_code']]; rawgi = {g: i for i, g in enumerate(z['genes'])}
        frame = scope.set_index('sampleID').loc[samples].reset_index()
        n = b['n_cells'].reshape(groups, ns)
        cp = b['mean_CP10K'].reshape(groups, ns, len(genes))
        lg = b['mean_log1p_CP10K'].reshape(groups, ns, len(genes))
        pos = b['positive_cell_counts'].reshape(groups, ns, len(genes))
        rna = b['full_RNA_total_sum'].reshape(groups, ns)
        ng = b['full_detected_gene_sum'].reshape(groups, ns)
        joint_cache = {}
        def joint(symbols, family, node=None):
            key = (tuple(symbols), family, node)
            if key not in joint_cache:
                if not all(g in rawgi for g in symbols):
                    joint_cache[key] = None
                else:
                    mask = cell_family == family
                    if node is not None: mask &= z['subCluster'] == node
                    good = np.all(z['raw_counts'][:, [rawgi[g] for g in symbols]] > 0, axis=1)
                    joint_cache[key] = np.bincount(z['sample_code'][mask & good], minlength=ns)
            return joint_cache[key]
        for interface in interfaces[interfaces.gse_id.eq(sid)].itertuples():
            pr = priors[(priors.gse_id.eq(sid)) & priors.sender_family.eq(interface.sender_family) & priors.receiver_family.eq(interface.receiver_family) & priors.ligand.eq(interface.ligand) & priors.receptor.eq(interface.receptor)].iloc[0]
            target_names = pr.top200_after_exclusions.split(';') if pr.top200_after_exclusions else []
            target_idx = [gi[g] for g in target_names if g in gi]
            Lgenes = interface.ligand.split('_'); Rgenes = interface.receptor.split('_')
            measured = all(g in gi for g in Lgenes + Rgenes)
            sender_indices = np.flatnonzero(b['node_family'] == interface.sender_family)
            sn = n[sender_indices].sum(axis=0)
            srna = rna[sender_indices].sum(axis=0); sng = ng[sender_indices].sum(axis=0)
            sender_cp = np.divide(np.nansum(cp[sender_indices] * n[sender_indices, :, None], axis=0), sn[:, None], out=np.full((ns,len(genes)),np.nan), where=sn[:, None]>0)
            context = frame.stratum.eq(interface.stratum).to_numpy() & frame.analysis_eligible.to_numpy(dtype=bool)
            for node in sorted(catalogue[catalogue.family.eq(interface.receiver_family)].subCluster):
                row = dict(interface_id=interface.interface_id, gse_id=sid, stratum=interface.stratum,
                           sender_family=interface.sender_family, receiver_family=interface.receiver_family,
                           ligand=interface.ligand, receptor=interface.receptor, receiver_original_node=node,
                           receiver_node_captured=node in nodes, ligand_target_prior_present=bool(pr.exact_ligand_target_column),
                           LR_status=pr.LR_status, all_LR_constituents_in_native_axis=measured,
                           complex_RNA_does_not_prove_assembled_LR=len(Lgenes)>1 or len(Rgenes)>1,
                           top_targets_fixed=len(target_names), top_targets_in_native_axis=len(target_idx),
                           context_libraries=int(context.sum()), diagnostic_design_columns=9,
                           final_inferential_design_not_evaluated=True)
                if node not in nodes:
                    row.update(status='RECEIVER_ORIGINAL_NODE_NOT_CAPTURED', n_libraries_endpoint_cells_ge20=0, n_analysis_blocks=0, n_reviewed_patient_blocks=0)
                    rows.append(row);continue
                k=nodes.index(node);rn=n[k];ok=context & (sn>=20) & (rn>=20)
                chosen=frame.loc[ok]
                blocks=chosen.effective_resampling_block.nunique()
                row.update(n_libraries_endpoint_cells_ge20=int(ok.sum()),n_analysis_blocks=int(blocks),
                           n_reviewed_patient_blocks=int(chosen.loc[chosen.unit_kind.eq('reviewed_patient'),'effective_resampling_block'].nunique()),
                           library_gate20=bool(ok.sum()>=20),block_support_ge10=bool(blocks>=10),block_support_ge20=bool(blocks>=20),
                           receiver_cells_in_eligible_units=int(rn[ok].sum()))
                jl=joint(Lgenes,interface.sender_family);jr=joint(Rgenes,interface.receiver_family,node)
                row['n_libraries_joint_LR_RNA_positive_cells_ge20']=int((ok & (jl>=20) & (jr>=20)).sum()) if jl is not None and jr is not None else None
                row['joint_RNA_panel_measurable']=jl is not None and jr is not None
                if target_idx and ok.any():
                    detected=pos[k,ok][:,target_idx].sum(axis=0)/rn[ok].sum()
                    expr=lg[k,ok][:,target_idx]
                    variable=np.var(expr,axis=0)>1e-12 if len(expr)>1 else np.zeros(len(target_idx),bool)
                    row.update(targets_detected1pct=int((detected>=.01).sum()),targets_detected10pct=int((detected>=.1).sum()),variable_target_genes=int(variable.sum()))
                else:row.update(targets_detected1pct=None,targets_detected10pct=None,variable_target_genes=None)
                if measured and ok.sum()>=2:
                    lv=np.mean(np.log1p(sender_cp[ok][:,[gi[g] for g in Lgenes]]),axis=1)
                    rv=np.mean(np.log1p(cp[k,ok][:,[gi[g] for g in Rgenes]]),axis=1)
                    a=np.column_stack([lv,rv,np.log1p(sn[ok]),np.log1p(rn[ok]),np.log1p(srna[ok]/sn[ok]),np.log1p(rna[k,ok]/rn[ok]),np.log1p(sng[ok]/sn[ok]),np.log1p(ng[k,ok]/rn[ok])])
                    sd=a.std(axis=0,ddof=1);v=np.divide(a-a.mean(axis=0),sd,out=np.zeros_like(a),where=sd>1e-12)
                    X=np.column_stack([np.ones(len(a)),v]);s=np.linalg.svd(X,compute_uv=False);rank=int((s>1e-10*s[0]).sum())
                    row.update(ligand_varies=bool(sd[0]>1e-12),receptor_varies=bool(sd[1]>1e-12),diagnostic_design_rank=rank,
                               diagnostic_design_full_rank=rank==9,constant_diagnostic_covariates=';'.join(np.asarray(design_names)[sd<=1e-12]),
                               diagnostic_design_condition=float(s[0]/s[-1]) if rank==9 else None)
                else:row.update(ligand_varies=False,receptor_varies=False,diagnostic_design_rank=None,diagnostic_design_full_rank=False,constant_diagnostic_covariates='NOT_EVALUABLE',diagnostic_design_condition=None)
                row['status']=('NO_EXACT_LIGAND_TARGET_PRIOR' if not row['ligand_target_prior_present'] else 'LR_GENE_NOT_MEASURED' if not measured else 'INSUFFICIENT_PAIRED_LIBRARIES_OR_BLOCKS' if not (ok.sum()>=20 and blocks>=10) else 'L_OR_R_CONSTANT' if not row['ligand_varies'] or not row['receptor_varies'] else 'NO_VARIABLE_NATIVE_TARGETS' if not row['variable_target_genes'] else 'DIAGNOSTIC_DESIGN_NOT_FULL_RANK' if not row['diagnostic_design_full_rank'] else 'COUNT_VARIATION_DIAGNOSTIC_DESIGN_SCREEN_PASSED_NOT_INFERENCE')
                rows.append(row)
        del z,b,cp,lg,pos,sender_cp;gc.collect()
    out.mkdir();result=pd.DataFrame(rows)
    for name in ['n_libraries_joint_LR_RNA_positive_cells_ge20','targets_detected1pct','targets_detected10pct','variable_target_genes','diagnostic_design_rank','receiver_cells_in_eligible_units']:
        result[name]=pd.array(result[name],dtype='Int64')
    result.to_csv(out/'all_interface_original_node_preflight.tsv',sep='\t',index=False,na_rep='NA')
    mapping.to_csv(out/'all336_to_unique_interfaces.tsv',sep='\t',index=False)
    summary=dict(frozen_records=336,unique_interfaces=len(interfaces),original_receiver_node_records=len(result),
                 source_context_node_status_counts=result.status.value_counts().to_dict(),
                 scope='Current3 audited sources,all336 mapping,all13 catalogue nodes separately; missing original nodes retained',
                 inference_not_accepted=True)
    inventory=dict(status='COMPLETE_PREFLIGHT_PENDING_INDEPENDENT_AUDIT',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),
                   producer_sha256=sha(__file__),plan_sha256=sha(plan_path),input_sha256=plan['input_sha256'],summary=summary,
                   outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_models_fitted=False,
                   formal_states=False,canonical_Zarr_changed=False,scientific_release=False,axis_validation_complete=False)
    (out/'inventory.json').write_text(json.dumps(inventory,indent=2)+'\n');print(json.dumps(inventory),flush=True)


if __name__=='__main__':main()

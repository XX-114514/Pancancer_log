"""Audit existing post-treatment NSCLC programme measurement without outcomes."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    source=root/'native_unit_expression_v1';family=root/'GSE243013_macrophage_CD8_family_cache_v1'
    accepted=json.loads((family/'independent_audit.json').read_text())
    assert accepted['status'].startswith('PASS_ALL200_SELECTION_KEYS_ALL21_STATE_MEMBERSHIPS400_NATIVE_WEIGHTED_FAMILIES')
    assert sha(family/'native_family_cache.npz')==accepted['cache_sha256']
    native_audit=json.loads((source/'independent_expression_audit.json').read_text())
    assert native_audit['status'].startswith('PASS_')
    V=Path(json.loads((root/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    projected=V/'state_after_DIALOGUE_frozen_projection_20261016_v1'
    out=root/'NSCLC_existing_frozen_programme_measurement_20261016_v1';assert not out.exists();out.mkdir()
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(Path(__file__)),
        question='Can frozen expression programmes be measured from existing native cell-summary extraction, before any new clinical effects?',
        source='Existing GSE243013 native all-cell/author-state extraction and frozen selected200 anti-PD1 cases; post-treatment immune-only, previously examined clinically.',
        old_programme_scope='Both original CRC MCP1 endpoints; no unsupported old MCP2 score manufactured.',
        new_scope='All56 fixed context/type/programme slots, including unsupported slots; no selection by internal-heldout or clinical P.',
        measurement='Unique exact target-cache gene identifiers and provided=True only; missing cache genes are unextracted, not proven native absence.',
        transfer_gate='For subsequent partial observed programme use, >=95% frozen absolute weight retained; complete measurement separately reported.',
        required_summary='Mean per-cell log1pCP10K supports a linear frozen gene projection; log1p(meanCP10K) does not. Cell quality means and native identities still require review.',
        family_limits='Coarse native lineage/author states need an explicit mapping; original172 fine identities not automatically accepted.',
        clinical_limits='Native case identifiers not all formally linked to independent patients; post-treatment association, not baseline prediction. No clinical test or new data download.')
    (out/'plan_before_measurement_coverage.json').write_text(json.dumps(plan,indent=2)+'\n')
    with np.load(source/'native_unit_expression.npz',allow_pickle=False) as z:
        keys=z.files;genes=z['genes'].astype(str);provided=z['provided'].astype(bool)
    assert len(genes)==len(provided) and len(set(genes))==len(genes)
    available=set(genes[provided]);targeted=set(genes)
    native_log_keys=[key for key in keys if 'log' in key.lower() and 'mean' in key.lower()]
    inventory=json.loads((source/'inventory.json').read_text())
    source_plan=json.loads((source/'plan.json').read_text())
    rows=[]
    def coverage(gene_weights,fields):
        assert gene_weights.gene.is_unique
        w=np.abs(gene_weights.weight.to_numpy(float));assert np.isfinite(w).all() and w.sum()>0
        present=gene_weights.gene.astype(str).isin(available).to_numpy()
        in_target=gene_weights.gene.astype(str).isin(targeted).to_numpy()
        frac=float(w[present].sum()/w.sum())
        return dict(**fields,programme_genes=len(w),provided_gene_count=int(present.sum()),
            targeted_but_not_provided=int((in_target&~present).sum()),not_in_existing_target_union=int((~in_target).sum()),
            absolute_weight_in_existing_cache=frac,complete_cache_measurement=bool(present.all()),
            fixed95pct_weight_gate=bool(frac>=.95),native_gene_absence_inferred=False)
    for endpoint in ['Tnode9','Macnode18']:
        path=root/'CRC_fine_pair_DIALOGUE_20261014_v1/private_transfer_fixed_weights_v1'/(endpoint+'_MCP1_genes.tsv')
        w=pd.read_csv(path,sep='\t')
        rows.append(coverage(w,dict(model='original_CRC_fixed_pair_MCP1',context='original_pair',endpoint=endpoint,program='MCP1',supported=True)))
    scope=pd.read_csv(projected/'public_slot_scope.tsv',sep='\t')
    weights=pd.read_csv(projected/'public_fixed_gene_weights.tsv',sep='\t')
    assert len(scope)==56
    for slot in scope.itertuples():
        selected=weights[weights.context.eq(slot.context)&weights.type.eq(slot.type)&weights.program.eq(slot.program)]
        fields=dict(model='2026-10-16_CONTEXT_DIALOGUE',context=str(slot.context),endpoint=slot.type,program=slot.program,
            supported=slot.status=='FROZEN_RECONSTRUCTION_PASS')
        if len(selected) and fields['supported']:
            rows.append(coverage(selected[['gene','weight']],fields))
        else:
            rows.append(dict(**fields,programme_genes=0,provided_gene_count=0,targeted_but_not_provided=0,
                not_in_existing_target_union=0,absolute_weight_in_existing_cache=np.nan,complete_cache_measurement=False,
                fixed95pct_weight_gate=False,native_gene_absence_inferred=False))
    table=pd.DataFrame(rows);assert len(table)==58
    table.to_csv(out/'all58_private_fixed_programme_cache_coverage.tsv',sep='\t',index=False,na_rep='NA')
    metadata=pd.read_csv(root/'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz',nrows=0)
    required=['sampleID','sub_cell_type','total_counts','n_genes_by_counts','cancer_type','anti-PD1_therapy',
        'chemotherapy','pathological_response','pathological_response_rate','gender','age','smoking_history']
    report=dict(status='ALL58_FROZEN_PROGRAMME_EXISTING_NSCLC_CACHE_MEASUREMENT_PREFLIGHT_NO_EFFECTS',
        finished_UTC=datetime.now(timezone.utc).isoformat(),native_summary_keys=keys,mean_log_summary_keys=native_log_keys,
        target_union_genes=len(genes),provided_target_genes=int(provided.sum()),frozen_programme_slots=len(table),
        fixed95pct_cache_eligible=int(table.fixed95pct_weight_gate.sum()),complete_cache_eligible=int(table.complete_cache_measurement.sum()),
        original_two_endpoints=table[table.model.eq('original_CRC_fixed_pair_MCP1')].to_dict('records'),
        missing_clinical_column_names=sorted(set(required)-set(metadata.columns)),
        native_extraction_status=inventory.get('status'),normalization=source_plan.get('normalization'),
        source_family_audit_sha256=sha(family/'independent_audit.json'),native_extraction_audit_sha256=sha(source/'independent_expression_audit.json'),
        programme_weights_sha256=sha(projected/'public_fixed_gene_weights.tsv'),
        coverage_sha256=sha(out/'all58_private_fixed_programme_cache_coverage.tsv'),
        independent_patient_linkage_review_pending=True,source_not_part_of_current_immune_atlas=True,
        clinical_effects_executed=False,new_expression_download=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()

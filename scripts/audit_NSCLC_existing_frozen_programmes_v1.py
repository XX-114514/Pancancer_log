"""Independently sum frozen absolute weights over existing target-cache genes."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=R/'NSCLC_existing_frozen_programme_measurement_20261016_v1'
    target=O/'private_independent_audit.json';assert not target.exists()
    report=json.loads((O/'private_terminal.json').read_text());table_path=O/'all58_private_fixed_programme_cache_coverage.tsv'
    assert hashlib.sha256(table_path.read_bytes()).hexdigest()==report['coverage_sha256']
    table=pd.read_csv(table_path,sep='\t',dtype={'context':str})
    with np.load(R/'native_unit_expression_v1/native_unit_expression.npz',allow_pickle=False) as z:
        all_names=list(map(str,z['genes']));flags=list(map(bool,z['provided']));keys=z.files
    provided={gene for gene,ok in zip(all_names,flags) if ok};targeted=set(all_names)
    V=Path(json.loads((R/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    P=V/'state_after_DIALOGUE_frozen_projection_20261016_v1'
    weights_path=P/'public_fixed_gene_weights.tsv';assert hashlib.sha256(weights_path.read_bytes()).hexdigest()==report['programme_weights_sha256']
    weights=pd.read_csv(weights_path,sep='\t');scope=pd.read_csv(P/'public_slot_scope.tsv',sep='\t')
    assert len(table)==58 and not table.duplicated(['model','context','endpoint','program']).any()
    for row in table.itertuples():
        if row.model=='original_CRC_fixed_pair_MCP1':
            w=pd.read_csv(R/'CRC_fine_pair_DIALOGUE_20261014_v1/private_transfer_fixed_weights_v1'/(row.endpoint+'_MCP1_genes.tsv'),sep='\t')
            supported=True
        else:
            match=scope[scope.context.eq(int(row.context))&scope.type.eq(row.endpoint)&scope.program.eq(row.program)]
            assert len(match)==1;supported=match.iloc[0].status=='FROZEN_RECONSTRUCTION_PASS'
            w=weights[weights.context.eq(int(row.context))&weights.type.eq(row.endpoint)&weights.program.eq(row.program)]
        assert supported==bool(row.supported)
        if not supported:
            assert row.programme_genes==0 and not row.fixed95pct_weight_gate and pd.isna(row.absolute_weight_in_existing_cache)
            continue
        pairs=list(zip(w.gene.astype(str),w.weight.astype(float)))
        total=sum(abs(value) for _,value in pairs)
        fraction=sum(abs(value) for gene,value in pairs if gene in provided)/total
        assert abs(fraction-row.absolute_weight_in_existing_cache)<1e-12
        assert row.programme_genes==len(pairs)
        assert row.provided_gene_count==sum(gene in provided for gene,_ in pairs)
        assert row.targeted_but_not_provided==sum(gene in targeted and gene not in provided for gene,_ in pairs)
        assert row.not_in_existing_target_union==sum(gene not in targeted for gene,_ in pairs)
        assert bool(row.fixed95pct_weight_gate)==(fraction>=.95)
    assert [key for key in keys if 'log' in key.lower() and 'mean' in key.lower()]==report['mean_log_summary_keys']
    audit=dict(status='PASS_ALL58_FIXED_GENE_WEIGHT_CACHE_COVERAGE_AND_UNEXTRACTED_NOT_ABSENT',
        finished_UTC=datetime.now(timezone.utc).isoformat(),coverage_sha256=report['coverage_sha256'],
        no_outcomes_or_effects=True,identity_mapping_or_patient_independence_not_validated=True)
    target.write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit))


if __name__=='__main__':main()

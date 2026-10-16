"""All58 frozen programme slots against complete native NSCLC gene axis.

Gene coverage does not validate matching fine identities or clinical endpoints.
No scores, learned weights or per-case records leave the private server.
"""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);N=R/'NSCLC_full_native_gene_units_20261016_v1'
    O=R/'NSCLC_full_native_programme_coverage_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
    terminal=json.loads((N/'private_terminal.json').read_text());gp=N/'private_native_gene_axis.tsv'
    assert sha(gp)==terminal['outputs'][gp.name]['sha256']
    native=pd.read_csv(gp,sep='\t',keep_default_na=False);assert native.gene.is_unique and native.source_measured.all()
    assert native.native_column.tolist()==list(range(31831));available=set(native.gene)
    V=Path(json.loads((R/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    F=V/'state_after_DIALOGUE_frozen_projection_20261016_v1';fi=json.loads((F/'inventory.json').read_text())
    for name in ['public_slot_scope.tsv','public_fixed_gene_weights.tsv']:assert sha(F/name)==fi['files'][name]
    scope=pd.read_csv(F/'public_slot_scope.tsv',sep='\t');weights=pd.read_csv(F/'public_fixed_gene_weights.tsv',sep='\t')
    prior=pd.read_csv(R/'NSCLC_existing_frozen_programme_measurement_20261016_v1/all58_private_fixed_programme_cache_coverage.tsv',sep='\t',dtype={'context':str})
    (O/'plan_before_native_coverage.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
        source_native_axis_sha256=sha(gp),frozen_weight_sha256=sha(F/'public_fixed_gene_weights.tsv'),
        scope='Original2 fixed CRC endpoint programmes plus all56 contextual DIALOGUE slots including3 unsupported/failed slots.',
        matching='Unique exact symbols on the already audited author uppercase gene axis; no fuzzy matching, refit, threshold reduction or imputation.',
        criterion='Complete measurement and original >=95% absolute weight criterion recorded separately.',
        audit_scope='Independent Python scalar sum and NumPy vector sum for every supported slot; full expression cache raw audit separately pending.',
        limitations='A measured gene programme is not a fine-cell identity match, functional state, independent patient linkage or clinical validation.'),indent=2)+'\n')
    rows=[];gene_rows=[]
    def add(w,fields):
        assert w.gene.is_unique and np.isfinite(w.weight).all()
        ws=w.weight.to_numpy(float);provided=w.gene.isin(available).to_numpy();total=np.abs(ws).sum();assert total>0
        f=float(np.abs(ws[provided]).sum()/total)
        scalar=sum(abs(float(value)) for gene,value in zip(w.gene,w.weight) if gene in available)/sum(abs(float(value)) for value in w.weight)
        assert abs(f-scalar)<1e-12
        rows.append(dict(**fields,programme_genes=len(w),native_provided_genes=int(provided.sum()),absolute_weight_native=f,
            complete_native_measurement=bool(provided.all()),native95_weight_gate=bool(f>=.95),fine_identity_match_validated=False))
        gene_rows.extend(dict(**fields,gene=g,weight=float(weight),native_provided=g in available) for g,weight in zip(w.gene,w.weight))
    for endpoint in ['Tnode9','Macnode18']:
        p=R/'CRC_fine_pair_DIALOGUE_20261014_v1/private_transfer_fixed_weights_v1'/(endpoint+'_MCP1_genes.tsv')
        add(pd.read_csv(p,sep='\t'),dict(model='original_CRC_fixed_pair_MCP1',context='original_pair',endpoint=endpoint,program='MCP1',supported=True))
    for row in scope.itertuples():
        fields=dict(model='2026-10-16_CONTEXT_DIALOGUE',context=str(row.context),endpoint=row.type,program=row.program,supported=row.status=='FROZEN_RECONSTRUCTION_PASS')
        w=weights[weights.context.eq(row.context)&weights.type.eq(row.type)&weights.program.eq(row.program)]
        if fields['supported']:add(w[['gene','weight']],fields)
        else:rows.append(dict(**fields,programme_genes=0,native_provided_genes=0,absolute_weight_native=np.nan,complete_native_measurement=False,native95_weight_gate=False,fine_identity_match_validated=False))
    result=pd.DataFrame(rows);assert len(result)==58
    keys=['model','context','endpoint','program'];assert set(result.set_index(keys).index)==set(prior.set_index(keys).index)
    result=result.merge(prior[keys+['absolute_weight_in_existing_cache','fixed95pct_weight_gate']],on=keys,validate='one_to_one')
    assert (result.loc[result.supported,'absolute_weight_native']+1e-12>=result.loc[result.supported,'absolute_weight_in_existing_cache']).all()
    result.to_csv(O/'all58_private_full_native_coverage.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(gene_rows).to_csv(O/'private_all_supported_gene_measurement.tsv',sep='\t',index=False)
    report=dict(status='ALL58_FULL_NATIVE_SYMBOL_COVERAGE_SCALAR_REPLAY_COMPLETE_NO_SCORES_OR_CLINICAL_TESTS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        slots=len(result),supported_slots=int(result.supported.sum()),native95_eligible=int(result.native95_weight_gate.sum()),complete_native=int(result.complete_native_measurement.sum()),
        original_two_endpoint_coverage=result[result.model.eq('original_CRC_fixed_pair_MCP1')].to_dict('records'),
        coverage_sha256=sha(O/'all58_private_full_native_coverage.tsv'),gene_axis_sha256=sha(gp),full_expression_independent_audit_required=True,
        fine_identity_mapping_validated=False,scores_computed=False,clinical_tests_executed=False)
    (O/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()

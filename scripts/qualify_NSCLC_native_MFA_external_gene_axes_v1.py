"""Every frozen native factor/view on every existing source gene axis; no effects."""
import json
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import load,sha,VIEWS
from NSCLC_native_MFA_clinical_inputs_v1 import matched_quota
from external_functional_coverage_inputs_v1 import read_inputs


def main():
    data=load();root=data['root'];clinical=root/'NSCLC_native_MFA_clinical_results_20261016_v1/private_independent_clinical_result_audit.json'
    assert json.loads(clinical.read_text())['status'].startswith('PASS_ALL40_')
    out=root/'NSCLC_native_MFA_external_gene_qualification_20261016_v1';assert not out.exists();matched_quota(root);out.mkdir(mode=0o700)
    provider=Path(__file__).with_name('external_functional_coverage_inputs_v1.py')
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),provider_sha256=sha(provider),
        prior_clinical_audit_sha256=sha(clinical),primary_model_sha256=sha(data['model_path']),
        question='After the fixed native clinical family, qualify all10 factors x4 frozen views against ALL existing source native physical gene axes, before any external score or association. No favorable factor, cancer or outcome selection.',
        model_coefficients='For each view A=solve(I+W.T*diag(Tau)*W, W.T*diag(Tau)); each factor row is the fixed own-view linear projection on selected expression inputs. No weight relearning or top-gene trimming.',
        gene_mapping='Use prior accepted source measured/native symbol/canonical/stable-alias physical-column lookup. Require unique gene->column and unique column->gene within each view. Ambiguous/colliding symbols unavailable, never expression-imputed.',
        physical_gate='>=500 uniquely mapped features; >=95% absolute coefficientL1, >=90% positive and negative coefficient mass where nonzero. All missing reasons and all source/factor slots retained.',
        important_limits='Physical coverage only. Does NOT establish matching CD8/Treg/macrophage/B-Plasma cells, count normalization, cross-platform score calibration, source/patient independence, or clinical sample eligibility. Exact original GSE243013 marked discovery-only; all other study overlaps require review.',
        existing_provider='Reuse established physical-axis loader; its old programme/reference input is read for provenance but NOT used for new coefficient or source selection.',
        new_expression_score_or_P=False,new_download=False,old_identity_or_model_changed=False)
    (out/'plan_before_external_native_gene_qualification.json').write_text(json.dumps(plan,indent=2)+'\n')
    existing=read_inputs(root);sources=existing['source_axes'];assert len(sources)==61
    rows=[]
    for vi,view in enumerate(VIEWS):
        w=data['model']['weights_'+str(vi)];tau=data['model']['precision_'+str(vi)]
        coefficient=np.linalg.solve(np.eye(10)+w.T@(tau[:,None]*w),w.T*tau)
        genes=data['genes'][data['prepared'][vi]['selected_native_columns']]
        assert coefficient.shape==(10,len(genes)) and np.isfinite(coefficient).all()
        for (arm,dataset),axis in sorted(sources.items()):
            candidates=[axis['lookup'].get(g,set()) for g in genes]
            physical=[next(iter(v)) if len(v)==1 else -1 for v in candidates]
            counts=pd.Series([p for p in physical if p>=0]).value_counts().to_dict()
            unique=np.array([p>=0 and counts.get(p,0)==1 for p in physical],bool)
            for k in range(10):
                v=coefficient[k];absolute=np.abs(v);positive=np.maximum(v,0);negative=np.maximum(-v,0)
                total=absolute.sum();p=positive.sum();n=negative.sum()
                coverage=float(absolute[unique].sum()/total) if total>1e-12 else np.nan
                pc=float(positive[unique].sum()/p) if p>1e-12 else np.nan;nc=float(negative[unique].sum()/n) if n>1e-12 else np.nan
                passed=bool(unique.sum()>=500 and coverage>=.95 and (p<=1e-12 or pc>=.9) and (n<=1e-12 or nc>=.9))
                rows.append(dict(arm=arm,dataset=dataset,view=view,factor='NSCLC_MFA_20261016_F%02d'%(k+1),selected_features=len(genes),
                    uniquely_mapped_features=int(unique.sum()),missing_symbols=sum(len(v)==0 for v in candidates),ambiguous_symbols=sum(len(v)>1 for v in candidates),
                    mapped_column_collisions=int(sum(len(v)==1 for v in candidates)-unique.sum()),absolute_coefficient_coverage=coverage,positive_mass_coverage=pc,negative_mass_coverage=nc,
                    physical_coverage_pass=passed,source_role='DISCOVERY_NOT_EXTERNAL' if dataset=='GSE243013' else 'OVERLAP_AND_PATIENT_INDEPENDENCE_UNREVIEWED',
                    compatible_lineage_capture_or_expression_calibration_accepted=False))
    table=pd.DataFrame(rows);assert len(table)==2440 and not table.duplicated(['arm','dataset','view','factor']).any()
    table.to_csv(out/'private_all2440_native_factor_source_gene_coverage.tsv',sep='\t',index=False)
    terminal=dict(status='COMPLETE_ALL2440_NATIVE_FACTOR_VIEW_SOURCE_PHYSICAL_GENE_COVERAGE_NO_EXPRESSION_OR_EXTERNAL_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        source_sha256=sha(__file__),files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
        old_learning_genes_used_for_selection=False,external_identity_or_clinical_validation_claim=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

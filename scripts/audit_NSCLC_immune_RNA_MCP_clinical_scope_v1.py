"""Independently verify all243 source/clinical keys and strict scope rules."""
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from preflight_NSCLC_native_MFA_clinical_coverage_v1 import COMMON, FIELDS


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    source=root/'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1'
    assert sha(source/'private_terminal.json')=='09d3fc4d485277429c0c88fa9a15b82a042d29e2f58883f8adbae971b887786a'
    t=json.loads((source/'private_terminal.json').read_text())
    for name,item in t['files'].items():assert sha(source/name)==item['sha256']
    out=root/'NSCLC_immune_RNA_MCP_clinical_scope_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_all243_key_scope_replay.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,
        scope='All author case keys/raw clinical fields, strict full-string therapy rule, original coding integrity, all eligibility and source summary counts. No RNA score or effect read.'),indent=2)+'\n')
    raw=[]
    for block in pd.read_csv(root/'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz',usecols=FIELDS,dtype=str,keep_default_na=False,chunksize=150000):raw.append(block.drop_duplicates())
    raw=pd.concat(raw).drop_duplicates().set_index('sampleID');assert len(raw)==243 and raw.index.is_unique
    original=pd.read_csv(root/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv',sep='\t',keep_default_na=False,dtype=str).set_index('sampleID')
    saved=pd.read_csv(source/'private_all243_current_RNA_clinical_scope.tsv',sep='\t',keep_default_na=False,dtype=str).set_index('sampleID')
    assert set(raw.index)==set(original.index)==set(saved.index)
    for name in FIELDS[1:]:assert saved.loc[raw.index,name].tolist()==raw[name].tolist()
    for name in COMMON+['response']:assert np.allclose(pd.to_numeric(saved.loc[raw.index,name],errors='coerce'),pd.to_numeric(original.loc[raw.index,name],errors='coerce'),equal_nan=True)
    agents='(?:Tislelizumab|Pembrolizumab|Camrelizumab|Sintilimab|Nivolumab|Toripalimab)'
    token=agents+r'(?:\(C\d+(?:-C\d+)?\))?';pattern=re.compile('^'+token+'(?: \\+ '+token+')*$')
    for case,row in raw.iterrows():
        s=saved.loc[case];drug=row['anti-PD1_therapy']
        role='EXPLICIT_ICI_AUTHOR_CASE' if drug=='Yes' or pattern.fullmatch(drug) else 'NO_ICI' if drug=='No' else 'AMBIGUOUS_OR_UNKNOWN_ICI'
        assert s.strict_treatment_role==role
        continuous=pd.to_numeric(s[COMMON],errors='coerce').to_numpy(float);complete=bool(np.isfinite(continuous).all())
        response={'pCR':1.,'MPR':1.,'non-MPR':0.}.get(row.pathological_response,np.nan)
        value=pd.to_numeric(pd.Series([s.response]),errors='coerce').iloc[0];assert (np.isnan(value) and np.isnan(response)) or value==response
        assert s.nonoutcome_controls_complete==str(complete)
        assert s.author_pathology_evaluable==str(bool(np.isfinite(response)))
        assert s.prospective_RNA_coordination_metadata_eligible==str(role=='EXPLICIT_ICI_AUTHOR_CASE' and complete)
        assert s.prospective_RNA_clinical_eligible==str(role=='EXPLICIT_ICI_AUTHOR_CASE' and complete and np.isfinite(response))
        assert s.reviewed_independent_patient=='False'
    coverage=pd.read_csv(source/'private_all_source_clinical_metadata_coverage.tsv',sep='\t');assert coverage.author_materials.sum()==243
    for r in coverage.itertuples(index=False):
        g=saved[saved.strict_treatment_role.eq(r.treatment_role)&saved.cancer_type.eq(r.cancer_type)]
        assert len(g)==r.author_materials
        for col in ['complete_nonoutcome_controls','pathology_evaluable','prospective_RNA_clinical_eligible','prospective_RNA_coordination_metadata_eligible']:
            field={'complete_nonoutcome_controls':'nonoutcome_controls_complete','pathology_evaluable':'author_pathology_evaluable'}.get(col,col)
            assert int(g[field].eq('True').sum())==getattr(r,col)
        response=pd.to_numeric(g.response,errors='coerce');assert response.eq(0).sum()==r.response0 and response.eq(1).sum()==r.response1 and r.reviewed_independent_patients==0
    accepted=dict(status='ACCEPTED_ALL243_ORIGINAL_CLINICAL_KEYS_STRICT_ROLES_AND_SOURCE_COVERAGE',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),
        all_author_records_checked=243,reviewed_independent_patients=0,no_RNA_or_clinical_effects=True)
    (out/'private_terminal.json').write_text(json.dumps(accepted,indent=2)+'\n')
    print(json.dumps(dict(status=accepted['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

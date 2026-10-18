"""All243 original author clinical records and strict therapy roles, no scores."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from preflight_NSCLC_native_MFA_clinical_coverage_v1 import COMMON, FIELDS


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    source=root/'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz'
    previous=root/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv'
    rolespath=root/'NSCLC_native_multicell_discovery_preflight_20261016_v1/private_all243_author_case_therapy_partition.tsv'
    assert sha(source)=='24155466454df53d01eb6cbbfcd56abeb45adb68f74e721af22a9dce41040b2f'
    assert sha(previous)=='29cb244219d628cd3509a73c30c59957964f9579fb4bd73e79aac9bf1aaccb7d'
    assert sha(rolespath)=='18e2bd76b2c87cee9965108e84f19f55b8b2643e73664df3520ae17446fa33cc'
    out=root/'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,
        input_sha256={p.name:sha(p) for p in [source,previous,rolespath]},
        scope='All243 native author materials. Exact original clinical columns and previously frozen codings, strict explicit-ICI role distinct from substring PD1_recorded. No old208/four-view selection; no RNA score or association read.',
        outcome='Author pCR/MPR versus non-MPR, no invented rate threshold or alternate outcome selection.',
        sampling_context='Official GSE/BioProject describes post-neoadjuvant surgical material. Source-level context, not new per-case stage review; not baseline prediction.',
        independence='Unique sampleID and source-level patient counts do not prove reviewed independent patient aliases. Keep author-case units and unverified status.',
        source_history='This dataset and its outcomes have previously been explored in project. Current RNA programme learned in other sources, but any new association is exploratory source transfer, not untouched project-wide validation.')
    (out/'plan_before_all243_clinical_scope.json').write_text(json.dumps(plan,indent=2)+'\n')
    records=[]
    for chunk in pd.read_csv(source,usecols=FIELDS,dtype=str,keep_default_na=False,chunksize=200000):records.append(chunk.drop_duplicates())
    raw=pd.concat(records).drop_duplicates().set_index('sampleID',drop=False)
    assert len(raw)==243 and raw.index.is_unique
    old=pd.read_csv(previous,sep='\t',dtype=str,keep_default_na=False).set_index('sampleID',drop=False)
    roles=pd.read_csv(rolespath,sep='\t',dtype=str,keep_default_na=False).set_index('sampleID')
    assert set(raw.index)==set(old.index)==set(roles.index)
    for field in FIELDS:assert raw[field].tolist()==old.loc[raw.index,field].tolist()
    assert raw['anti-PD1_therapy'].tolist()==roles.loc[raw.index,'therapy_original'].tolist()
    assert raw.cancer_type.tolist()==roles.loc[raw.index,'cancer_type'].tolist()
    expected=raw.pathological_response.map({'pCR':1.,'MPR':1.,'non-MPR':0.})
    assert np.allclose(expected,pd.to_numeric(old.loc[raw.index,'response'],errors='coerce'),equal_nan=True)
    joined=old.loc[sorted(old.index)].copy()
    for field in COMMON+['response']:joined[field]=pd.to_numeric(joined[field],errors='coerce')
    joined['strict_treatment_role']=roles.loc[joined.index,'treatment_role']
    joined['reviewed_independent_patient']=False
    joined['nonoutcome_controls_complete']=np.isfinite(joined[COMMON].to_numpy(float)).all(axis=1)
    joined['author_pathology_evaluable']=np.isfinite(joined.response)
    joined['prospective_RNA_clinical_eligible']=joined.strict_treatment_role.eq('EXPLICIT_ICI_AUTHOR_CASE')&joined.nonoutcome_controls_complete&joined.author_pathology_evaluable
    joined['prospective_RNA_coordination_metadata_eligible']=joined.strict_treatment_role.eq('EXPLICIT_ICI_AUTHOR_CASE')&joined.nonoutcome_controls_complete
    # Keep all exclusions. Actual RNA-view coverage and source transport remain
    # a separate future gate, not an assumption at this metadata-only stage.
    joined.to_csv(out/'private_all243_current_RNA_clinical_scope.tsv',sep='\t',index=False)
    summary=[]
    for (role,cancer),g in joined.groupby(['strict_treatment_role','cancer_type'],sort=True):
        summary.append(dict(dataset='GSE243013',source_sampling_context='post_neoadjuvant_surgery',treatment_role=role,cancer_type=cancer,
            author_materials=len(g),complete_nonoutcome_controls=int(g.nonoutcome_controls_complete.sum()),
            pathology_evaluable=int(g.author_pathology_evaluable.sum()),response0=int(g.response.eq(0).sum()),response1=int(g.response.eq(1).sum()),
            prospective_RNA_clinical_eligible=int(g.prospective_RNA_clinical_eligible.sum()),
            prospective_RNA_coordination_metadata_eligible=int(g.prospective_RNA_coordination_metadata_eligible.sum()),reviewed_independent_patients=0))
    pd.DataFrame(summary).to_csv(out/'private_all_source_clinical_metadata_coverage.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()};assert sum(v['bytes'] for v in files.values())<4*1024**2
    terminal=dict(status='COMPLETE_ALL243_SOURCE_CLINICAL_KEYS_ROLES_AND_COVERAGE_NO_RNA_EFFECTS',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,
        author_materials=243,reviewed_independent_patients=0,RNA_or_clinical_association_read=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

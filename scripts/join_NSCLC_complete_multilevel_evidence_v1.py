"""Connect every original pair/function across immutable evidence stages, no new tests."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);Broot=r/'NSCLC_third_lineage_background_descriptions_20261016_v1';Mroot=r/'NSCLC_matched_function_coordination_20261016_v1';Hroot=r/'NSCLC_paired_clinical_complete_evidence_review_20261016_v1'
    proofs=dict(third=Broot/'private_independent_third_lineage_audit.json',matched=Mroot/'private_independent_matched_function_audit.json',historical=Hroot/'private_independent_display_audit.json')
    for p in proofs.values():assert json.loads(p.read_text())['status'].startswith('PASS_')
    paths=dict(third=Broot/'private_all28600_third_lineage_descriptions.tsv',matched=Mroot/'private_all28600_matched_function_descriptions.tsv',historical=Hroot/'private_all572_complete_same_context_evidence.tsv')
    out=r/'NSCLC_complete_multilevel_functional_evidence_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),input_sha256={k:sha(p) for k,p in paths.items()},input_audit_sha256={k:sha(p) for k,p in proofs.items()},scope='Every original572 pair x50 function with all historical single/paired clinical and coordination fields, original matched-function fields and third-lineage subset fields; no significance or effect filtering.',
              limits='Different historical single-end and B-captured case sets remain separate. All P/q keep their original names/families and never apply to new functional descriptions. No new P/FDR/candidate-confirmation/independent-patient claim.')
    (out/'private_connection_plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    tables={k:pd.read_csv(p,sep='\t') for k,p in paths.items()};keys=['record','context','program','author_subtype_A','author_subtype_B'];fullkeys=keys+['reference'];base=tables['third'];assert len(base)==len(tables['matched'])==28600 and len(tables['historical'])==572
    assert not base.duplicated(fullkeys).any() and not tables['matched'].duplicated(fullkeys).any() and not tables['historical'].duplicated(keys).any()
    assert set(map(tuple,base[fullkeys].to_numpy()))==set(map(tuple,tables['matched'][fullkeys].to_numpy()))
    assert set(map(tuple,base[keys].to_numpy()))==set(map(tuple,tables['historical'][keys].to_numpy()))
    result=base.copy();orders=base[fullkeys].copy()
    for name,source_keys in [('matched',fullkeys),('historical',keys)]:
        f=tables[name];new=f.rename(columns={col:name+'_'+col for col in f.columns if col not in source_keys});result=result.merge(new,on=source_keys,how='left',sort=False,validate='one_to_one' if name=='matched' else 'many_to_one')
    pd.testing.assert_frame_equal(result[fullkeys],orders);assert len(result)==28600
    result['new_functional_P_or_BH_available']=False;result['source_case_identity_fully_reviewed']=False;result['old_and_new_case_sets_interchangeable']=False
    path=out/'private_all28600_complete_multilevel_evidence.tsv';result.to_csv(path,sep='\t',index=False,na_rep='NA')
    # Reconstruct each source's closed key lookup separately from pandas merge.
    disk=pd.read_csv(path,sep='\t');assert disk[fullkeys].equals(orders)
    for name,source_keys in [('third',fullkeys),('matched',fullkeys),('historical',keys)]:
        f=tables[name].set_index(source_keys);idx=pd.MultiIndex.from_frame(disk[source_keys]);expected=f.reindex(idx).reset_index(drop=True)
        fields=[col for col in f.columns];observed=disk[[col if name=='third' else name+'_'+col for col in fields]].copy();observed.columns=fields
        pd.testing.assert_frame_equal(observed.reset_index(drop=True),expected,check_dtype=False,check_exact=False,rtol=1e-10,atol=1e-12)
    assert not disk.new_functional_P_or_BH_available.any() and not disk.source_case_identity_fully_reviewed.any() and not disk.old_and_new_case_sets_interchangeable.any()
    text=['# Complete multilevel programme/function evidence (private)','',
          'Every572 original programme/author-subtype pair has all50 functional positions retained. Exact context, MCP and both endpoint labels are immutable keys.',
          'Current fields: explicit B-captured cases and all before/after third-lineage descriptions. matched_ fields: previous whole matched-case functional descriptions. historical_ fields: original clinical/coordination results, with their original family names and case counts.',
          'Historical single-end case sets can differ from paired sets. Third-lineage subset counts can differ from original paired sets. New functions have no new P/BH/CI and cannot inherit historical q.',
          'No table row establishes direct signalling, independent patient replication, protein activity, clinical prediction or causal mediation. Preserve original failed/unsupported families. B-cell conditioning can reflect real shared ecology as well as sample environment.','',
          'The TSV was checked against all fields of all three immutable source tables by separate closed-key reindexing after CSV roundtrip. Upstream OLS/measurement audits remain the scientific verification; this is a full data-connection consistency check.']
    (out/'PRIVATE_COMPLETE_MULTILEVEL_REVIEW.md').write_text('\n'.join(text)+'\n')
    terminal=dict(status='PASS_COMPLETE_ALL28600_EXACT_CONTEXT_MCP_ENDPOINT_FUNCTION_CONNECTION_AND_ALL_SOURCE_FIELD_LOOKUPS_NO_NEW_INFERENCE',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),rows=len(disk),columns=len(disk.columns),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_P_BH=False,goal_complete=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()

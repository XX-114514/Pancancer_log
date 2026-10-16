"""Independent exact-key and evidence-boundary audit of the complete clinical review."""
import hashlib,json,os,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from PIL import Image


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def truth(x):return str(x).lower()=='true'


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=r/'NSCLC_paired_clinical_complete_evidence_review_20261016_v1';rc=r/'HNSCC_receiver_followup_job_control_20261016_v2/review_NSCLC_paired_clinical_evidence_v1'
    receipt=json.loads((rc/'launch_receipt.json').read_text());deadline=time.monotonic()+21600
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split()
        assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(10)
    terminal=json.loads((out/'private_terminal.json').read_text())
    for n,h in terminal['files'].items():assert sha(out/n)==h
    clinical_root=r/'NSCLC_paired_programme_clinical_increment_results_20261016_v1';assert json.loads((clinical_root/'private_independent_audit.json').read_text())['status'].startswith('PASS_ALL572')
    t=pd.read_csv(clinical_root/'private_all572_paired_clinical_results.tsv',sep='\t').set_index('record');c=pd.read_csv(r/'NSCLC_two_end_coordination_results_20261016_v1/private_all572_two_end_results.tsv',sep='\t').set_index('record')
    single=pd.read_csv(r/'NSCLC_frozen_signature_clinical_results_20261016_v1/private_all_planned_signature_clinical_results.tsv',sep='\t').set_index(['context','program','type','sub_cell_type']);assert single.index.is_unique
    merged=pd.read_csv(out/'private_all572_complete_same_context_evidence.tsv',sep='\t').set_index('record');assert len(merged)==len(t)==len(c)==572 and merged.index.equals(t.index)
    for rec,z in merged.iterrows():
        a=t.loc[rec];b=c.loc[rec]
        for key in ['context','program','author_subtype_A','author_subtype_B']:assert a[key]==b[key]==z[key]
        for key in t.columns:
            if pd.api.types.is_numeric_dtype(t[key]):assert np.isclose(a[key],z[key],atol=1e-10,rtol=1e-10,equal_nan=True)
            elif pd.isna(a[key]):assert pd.isna(z[key])
            else:assert str(a[key])==str(z[key])
        for typ,end in [('typeA','A'),('typeB','B')]:
            old=single.loc[(a.context,a.program,typ,a['author_subtype_'+end])]
            assert old.eligible_cases==z['old_'+end+'_single_end_cases'];assert np.isclose(old.BH_q_full_planned_family,z['old_'+end+'_single_end_BH_q'],equal_nan=True)
        assert np.isclose(b.BH_q_full572,z.old_coordination_BH_q,equal_nan=True)
        assert np.isclose(b.partial_rank_rho,z.old_pathology_controlled_coordination_rho,equal_nan=True)
        assert np.isclose(b.nonshared_gene_partial_rho,z.old_nonshared_coordination_rho,equal_nan=True)
        estimable=pd.notna(a.pair_IUT_P)
        ci=all(a[k+'_CI_low']>0 or a[k+'_CI_high']<0 for k in ['A','B']) if estimable else False
        both=bool(estimable and a.BH_q_pair572<.05 and ci and truth(a.A_LOO_same_sign) and truth(a.B_LOO_same_sign))
        coord=bool(pd.notna(b.BH_q_full572) and b.BH_q_full572<.05)
        assert truth(z.both_conditional_clinical_supported)==both and truth(z.coordination_supported_same_cases)==coord
        expected='INPUT_INSUFFICIENT' if not estimable else 'TWO_END_CONDITIONAL_CLINICAL_AND_COORDINATION' if both and coord else 'TWO_END_CONDITIONAL_CLINICAL_ONLY' if both else 'COORDINATED_PROGRAMME_WITHOUT_VERIFIED_BOTH_END_INCREMENT' if coord else 'NO_VERIFIED_PAIRED_CLINICAL_INCREMENT'
        assert z.evidence_class==expected
        if estimable:
            new=np.load(r/'NSCLC_paired_programme_clinical_increment_preflight_20261016_v1/private_designs'/(rec+'.npz'),allow_pickle=False)
            old=np.load(r/'NSCLC_two_end_coordination_preflight_20261016_v1/private_pair_designs'/(rec+'.npz'),allow_pickle=False)
            assert np.array_equal(new['cases'],old['cases']) and np.array_equal(new['Y'],old['response']) and truth(z.exact_case_join_verified)
    display=pd.read_csv(out/'private_top12_display_only.tsv',sep='\t');e=merged[merged.pair_IUT_P.notna()].reset_index().sort_values(['BH_q_pair572','pair_IUT_P','record']).head(12)
    assert display.record.tolist()==e.record.tolist()
    files=[]
    for context,program in merged[merged.pair_IUT_P.notna()][['context','program']].drop_duplicates().itertuples(index=False):
        p=out/f'private_complete_forest_C{context}_{program}.png';assert p.exists()
        with Image.open(p) as im:im.verify()
        files.append(p.name)
    assert len(files)==8
    proof=dict(status='PASS_ALL572_EXACT_CASE_CONTEXT_ENDPOINT_CLINICAL_COORDINATION_JOINS_CLASSIFICATIONS_AND_ALL_FOREST_FILES_VISUAL_REVIEW_SEPARATE',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_code_sha256=sha(__file__),
               producer_terminal_sha256=sha(out/'private_terminal.json'),full_rows_checked=len(merged),complete_forest_files=files,new_P_or_independent_cohort_claim=False)
    p=out/'private_independent_display_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()

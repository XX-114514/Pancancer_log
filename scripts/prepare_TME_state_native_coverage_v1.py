"""Register fixed state-resource coverage on accepted native gene axes only.

No expression matrix, per-cell state score or clinical outcome is read. Passing
coverage is not programme applicability, function validation or CM admission.
"""
import argparse,collections,datetime,hashlib,json,sys
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def main(R):
    cfg=load(R/'config/spec_v1.json');H=Path(cfg['model']);W=Path(cfg['wilcoxon']);wc=load(W/'config/spec_v1.json');M=Path(wc['original_marker_stage']);A=Path(wc['native_profile_auditor'])
    sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(R,cfg);ap=A/'outputs/all_native_profiles_acceptance_v1.json';assert sha(ap)==wc['native_profile_phase_acceptance_SHA256'];phase=load(ap);assert phase['all_native_source_gene_axes_and_all_group_sufficient_statistics_replayed']
    # This phase was accepted before the distinct official Wilcoxon phase.
    # A later rank/auditor failure cannot erase its already checked statistics.
    consumer=Path(cfg['accepted_profile_consumer']);assert load(consumer/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0;ip=consumer/'outputs/private_matched_reference_profile_index_v1.json';assert sha(ip)==load(consumer/'provenance/terminal_v1.json')['index_SHA256'];assert load(ip)['native_profile_acceptance_SHA256']==sha(ap)
    curpath=Path(cfg['curation']);cur=load(curpath);assert cur['scored_catalog_SHA256']==cfg['catalog_SHA256'];matches=[Path(p) for p in cfg['catalog_candidates'] if Path(p).exists() and sha(p)==cfg['catalog_SHA256']];assert matches,'No exact frozen catalogue bytes available';catpath=matches[0];cat=load(catpath);assert cur['scoring_genes_unchanged'];assert not cur['all_module_admission'];scope=load(H/'outputs/scope_and_common_genes_v1.json');common=set(scope['genes'])
    cmeta={r['program_id']:r for r in cur['program_metadata']};assert len(cmeta)==len(cat['programs']);assert len({r['program_id'] for r in cat['programs']})==len(cat['programs'])
    verified={r['source']:r['source_profile_SHA256'] for r in load(ap)['source_checks']};index=load(M/'outputs/native_profiles_index_v1.json');rows=[];source_rows=[]
    for r in index:
        p=M/'outputs'/r['path'];assert sha(p)==r['SHA256']==verified[r['source']]
        with np.load(p,allow_pickle=False) as f:genes=f['native_gene_names'].astype(str).tolist();n=f['n_cells'].copy()
        assert len(genes)==len(set(genes)) and common<=set(genes);native=set(genes);status=collections.Counter()
        for pr in cat['programs']:
            m=cmeta[pr['program_id']];original=set(pr['genes']);assert original;measured=sorted(original&native);shared=sorted(original&common);eligible=len(measured)>=cfg['minimum_genes'] and len(measured)/len(original)>=cfg['minimum_fraction'];common_ok=len(shared)>=cfg['minimum_genes'] and len(shared)/len(original)>=cfg['minimum_fraction'];status['native_coverage_pass' if eligible else 'native_coverage_fail']+=1
            rows.append(dict(source=r['source'],program_id=pr['program_id'],n_original_raw=len(pr['genes']),n_original_unique=len(original),n_native_measured=len(measured),native_gene_fraction=len(measured)/len(original),native_coverage_pass=eligible,n_common_measured=len(shared),common_coverage_pass=common_ok,native_measured_genes=measured,missing_exact_symbols=sorted(original-native),reference_parents=m['reference_parents'],exact_Cytopus_gene_match=m['exact_Cytopus_gene_match'],citation_status=m['citation_status'],module_admitted=False,identity_applicability_pending=True))
        source_rows.append(dict(source=r['source'],n_fit=int(n.sum()),n_unique_native_measured_genes=len(native),native_profile_SHA256=sha(p),**status))
    assert [r['source'] for r in index]==scope['sources'];assert sum(r['n_fit'] for r in source_rows)==scope['n_fit_records']==phase['n_fit'];assert len(index)==phase['n_sources'];assert len(rows)==len(index)*len(cat['programs'])
    save(R/'outputs/private_state_program_source_native_coverage_v1.json',rows);save(R/'outputs/private_source_state_measurement_coverage_v1.json',source_rows)
    manifest=dict(schema='TME_state_resource_native_coverage_preparation_v1',n_sources=len(index),n_fit=sum(r['n_fit'] for r in source_rows),n_programs=len(cat['programs']),n_source_program_records=len(rows),catalog_SHA256=sha(catpath),curation_SHA256=sha(curpath),native_profile_acceptance_SHA256=sha(ap),accepted_profile_consumer_index_SHA256=sha(ip),native_statistics_phase_not_whole_Wilcoxon_accepted=True,minimum_genes=cfg['minimum_genes'],minimum_fraction=cfg['minimum_fraction'],exact_symbols_only_no_fuzzy_aliases=True,native_and_common_ranking_universes_not_conflated=True,expression_read=False,state_scores_computed=False,coverage_not_biological_state_or_module_admission=True,original_citation_review_still_separate=True,cell_state_module_fitting_not_started=True)
    save(R/'outputs/state_native_coverage_preparation_v1.json',manifest);storage_gate(R,cfg)
    save(R/'provenance/terminal_v1.json',dict(status='FIXED_STATE_RESOURCE_NATIVE_GENE_COVERAGE_REGISTERED_NO_STATE_SCORES',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),manifest_SHA256=sha(R/'outputs/state_native_coverage_preparation_v1.json')))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise

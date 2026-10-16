"""Independent cell-code recount and set intersections for all capture opportunities."""
import hashlib,json,os,time
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd,zarr
from anndata._io.specs import read_elem
from external_lineage_qualification_inputs_v2 import load,STRATA


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=r/'complete_external_lineage_opportunities_20261016_v2';job=r/'HNSCC_receiver_followup_job_control_20261016_v2/preflight_complete_external_lineage_opportunities_v2';receipt=json.loads((job/'launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    terminal=json.loads((out/'private_terminal.json').read_text());plan=json.loads((out/'plan_before_complete_capture_counts.json').read_text());assert plan['code_sha256']==receipt['source_sha256']
    for name,h in terminal['files'].items():assert sha(out/name)==h
    data=load(r);assert data['input_sha256']==plan['input_sha256'];assert sha(Path(__file__).with_name('external_lineage_qualification_inputs_v2.py'))==plan['input_helper_sha256']
    lib=pd.read_csv(out/'private_all_library_latest_metadata.tsv',sep='\t',keep_default_na=False,dtype={'dataset':str,'library_id':str,'patient_block':str,'resampling_block':str});pd.testing.assert_frame_equal(lib,data['libraries'],check_dtype=False)
    got=pd.read_csv(out/'private_all_original_node_library_counts.tsv.gz',sep='\t',dtype={'dataset':str,'library_id':str,'node_id':str});keys=['arm','dataset','library_id','node_id','major'];pd.testing.assert_frame_equal(got.sort_values(keys).reset_index(drop=True),data['counts'].sort_values(keys).reset_index(drop=True),check_dtype=False)
    # Fresh categorical mixed-radix coding of every atlas cell, independent of groupby.
    atlas=r.parents[2]/'Project_v5/runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr';z=zarr.open_group(str(atlas),mode='r');cols=['dataset','library_id','node_id','major'];cats=[pd.Categorical(np.asarray(read_elem(z['obs'][k])).astype(str)) for k in cols];encoded=np.zeros(len(cats[0]),np.int64)
    for v in cats:assert (v.codes>=0).all();encoded=encoded*len(v.categories)+v.codes
    _,first,n=np.unique(encoded,return_index=True,return_counts=True);fresh=pd.DataFrame({col:np.asarray(v)[first].astype(str) for col,v in zip(cols,cats)});fresh.insert(0,'arm','archived_immune_atlas');fresh['n_cells']=n
    observed=got[got.arm.eq('archived_immune_atlas')];pd.testing.assert_frame_equal(observed.sort_values(keys).reset_index(drop=True),fresh.sort_values(keys).reset_index(drop=True),check_dtype=False);assert n.sum()==2670057
    definitions=pd.read_csv(out/'private_all_pool_node_definitions.tsv',sep='\t',dtype=str);pools=data['pool_nodes']
    assert {(a,p,n) for a,pp in pools.items() for p,ns in pp.items() for n in ns}==set(map(tuple,definitions[['arm','pool','node_id']].to_numpy()))
    cells=defaultdict(dict)
    for row in got.itertuples():cells[(row.arm,row.dataset,row.library_id)][row.node_id]=int(row.n_cells)
    counts=pd.read_csv(out/'private_all_library_pool_counts.tsv.gz',sep='\t',dtype={'dataset':str,'library_id':str});count_lookup={};eligible=defaultdict(set);meta={};strata=defaultdict(set);sources=defaultdict(set)
    for row in lib.itertuples():
        k=(row.arm,row.dataset,row.library_id);assert k in cells;meta[k]=row;sources[(row.arm,row.dataset)].add(k);strata[(row.arm,row.dataset,*[getattr(row,s) for s in STRATA])].add(k)
    for row in counts.itertuples():
        k=(row.arm,row.dataset,row.library_id);wanted=sum(value for node,value in cells[k].items() if node in pools[row.arm][row.pool]);assert row.n_cells==wanted
        kk=(*k,row.pool);assert kk not in count_lookup;count_lookup[kk]=wanted
        if wanted>=20:eligible[(row.arm,row.pool)].add(k)
    assert len(count_lookup)==sum(len(pools[x.arm]) for x in lib.itertuples())
    def check(summary,keys):
        verified={meta[k].patient_block for k in keys if meta[k].reviewed_patient};fallback=sum(not meta[k].reviewed_patient for k in keys)
        assert summary.libraries==len(keys) and summary.reviewed_patient_blocks==len(verified) and summary.unreviewed_library_fallback==fallback
    raw=pd.read_csv(out/'private_all_source_raw_major_pairs.tsv',sep='\t');ss=pd.read_csv(out/'private_all_source_stratum_pairs.tsv.gz',sep='\t',keep_default_na=False);candidate=pd.read_csv(out/'private_all_frozen_candidate_pair_opportunities.tsv',sep='\t');paired=pd.read_csv(out/'private_all_HNSCC_reviewed_paired_opportunities.tsv',sep='\t');zero=False
    for table in [raw,candidate]:
        assert not table.duplicated(['arm','dataset','lineage_A','lineage_B']).any()
        for row in table.itertuples():
            universe=sources[(row.arm,row.dataset)];members=universe&eligible[(row.arm,row.lineage_A)]&eligible[(row.arm,row.lineage_B)];check(row,members);assert row.all_source_libraries==len(universe) and not row.new_model_or_independence_accepted;zero|=len(members)==0
            if row.scope=='frozen_candidate_extension':
                triple=members&eligible[(row.arm,'extension::B_plus_Plasma')];assert row.B_plus_Plasma_libraries==len(triple) and row.B_plus_Plasma_reviewed_patient_blocks==len({meta[k].patient_block for k in triple if meta[k].reviewed_patient})
    for row in ss.itertuples():
        universe=strata[(row.arm,row.dataset,*[getattr(row,s) for s in STRATA])];members=universe&eligible[(row.arm,row.lineage_A)]&eligible[(row.arm,row.lineage_B)];check(row,members);assert row.all_stratum_libraries==len(universe) and not row.new_model_or_independence_accepted
    for row in paired.itertuples():
        members=sources[('archived_immune_atlas',row.dataset)]&eligible[('archived_immune_atlas',row.lineage_A)]&eligible[('archived_immune_atlas',row.lineage_B)]
        if row.scope=='two_end_plus_B_Plasma':members &= eligible[('archived_immune_atlas','extension::B_plus_Plasma')]
        phases=defaultdict(set);pre=post=other=0
        for key in members:
            m=meta[key]
            if m.timepoint not in ['pre','post']:other+=1
            if not m.reviewed_patient or m.timepoint not in ['pre','post']:continue
            phases[m.patient_block].add(m.timepoint);pre+=m.timepoint=='pre';post+=m.timepoint=='post'
        assert row.reviewed_pre_post_patient_pairs==sum({'pre','post'}<=p for p in phases.values()) and row.pre_libraries==pre and row.post_libraries==post and row.other_phase_libraries==other
    assert len(raw)==61*153 and len(candidate)==18 and len(paired)==4 and zero
    assert not lib.untouched_independent_validation_accepted.any() and not lib.new_clinical_model_accepted.any() and not lib.exact_NSCLC_endpoint_bridge.any()
    proof=dict(status='PASS_ALL61_SOURCE_ARM_ALL1690_LIBRARY_ALL_INHERITED_NODE_AND_POOL_COUNTS_ALL_MAJOR_CANDIDATE_STRATA_AND_PATIENT_PAIR_OPPORTUNITIES',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),input_helper_sha256=plan['input_helper_sha256'],producer_terminal_sha256=sha(out/'private_terminal.json'),fresh_atlas_cells_recounted=2670057,raw_major_pair_rows_checked=len(raw),stratum_rows_checked=len(ss),candidate_rows_checked=len(candidate),zero_capture_and_patient_repeat_phase_cases_checked=True,new_expression_scores_P_or_clinical_acceptance=False)
    p=out/'private_independent_opportunity_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()

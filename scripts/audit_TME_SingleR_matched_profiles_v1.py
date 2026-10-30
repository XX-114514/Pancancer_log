"""Check complete group-reference axes and replay full source query batches.

Full source batches preserve SingleR's query-dependent pruning distribution;
subsampling query groups would not be an equivalent pruning replay.
"""
import argparse,datetime,hashlib,json,os,signal,subprocess,sys
from pathlib import Path
import numpy as np
import pandas as pd

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,value):
    with Path(p).open('x') as f:json.dump(value,f,indent=2,allow_nan=False)
def main(R):
    cfg=load(R/'config/spec_v1.json');S=Path(cfg['annotation']);P=Path(cfg['profiles'])
    sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(R,cfg)
    for root in [S,P]:assert load(root/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    manifest=S/'outputs/matched_SingleR_proposal_manifest_v1.json';assert sha(manifest)==load(S/'provenance/terminal_v1.json')['manifest_SHA256']
    m=load(manifest);idxp=P/'outputs/private_matched_reference_profile_index_v1.json';assert sha(idxp)==m['profile_index_SHA256'];idx=load(idxp)
    refs=load(S/'provenance/frozen_public_reference_inputs_v1.json');t=pd.read_csv(S/'outputs/private_all_matched_SingleR_proposals_v1.tsv',sep='\t');expected={};sizes={};n_profiles=0
    for s in idx['metadata']:
        p=P/'outputs'/s['private_profile_file'];assert sha(p)==s['profile_SHA256'];sizes[s['source']]=sum(g['n_cells'] for g in s['groups'])
        for g in s['groups']:expected[(s['source'],g['cluster'])]=g['n_cells'];n_profiles+=1
    assert n_profiles==idx['n_profiles']==m['n_profiles'];assert sum(expected.values())==idx['n_fit'];assert len(t)==4*n_profiles
    assert not t[['source','cluster','reference','label_level']].duplicated().any()
    coverage=[]
    for ref in refs:
        assert sha(ref['matrix_path'])==ref['matrix_SHA256'] and sha(ref['labels_path'])==ref['labels_SHA256']
        vocab=pd.read_csv(S/'outputs'/('private_reference_'+ref['reference']+'_vocabulary_v1.tsv'),sep='\t')
        for level in ['main','fine']:
            q=t[(t.reference==ref['reference'])&(t.label_level==level)];assert set(zip(q.source,q.cluster))==set(expected)
            assert set(q.status)=={'official_SingleR_group_reference_proposal'};assert (q.n_shared_reference_genes>=cfg['minimum_shared_genes']).all();assert np.isfinite(q.delta_next).all()
            assert set(q.raw_label)<=set(vocab['label_'+level]);assert set(q.pruned_label.dropna())<=set(vocab['label_'+level]);assert (q.pruned_label.dropna()==q.loc[q.pruned_label.notna(),'raw_label']).all()
            coverage.append(dict(reference=ref['reference'],level=level,n_profiles=len(q),n_pruned=int(q.pruned_label.isna().sum()),n_cells_in_pruned_groups=sum(expected[(r.source,r.cluster)] for r in q.itertuples() if pd.isna(r.pruned_label))))
    pruned=t[t.pruned_label.isna()].groupby('source').size()
    probes=sorted({sorted(sizes)[0],max(sizes,key=sizes.get),pruned.idxmax() if len(pruned) else sorted(sizes)[-1]})
    pd.DataFrame([dict(source=s,profile_path=str(P/'outputs'/next(r['private_profile_file'] for r in idx['metadata'] if r['source']==s))) for s in probes]).to_csv(R/'config/private_replay_sources_v1.tsv',sep='\t',index=False)
    pd.DataFrame(refs).to_csv(R/'config/private_replay_references_v1.tsv',sep='\t',index=False)
    target=R/'outputs/private_official_query_replay_v1.tsv'
    cmd=[cfg['Rscript'],'--vanilla',str(R/'scripts/audit_TME_SingleR_matched_profiles_v1.R'),str(R/'config/private_replay_references_v1.tsv'),str(R/'config/private_replay_sources_v1.tsv'),str(S/'outputs'),str(target),m['SingleR_version']]
    with (R/'provenance/private_R_stdout_v1.txt').open('x') as o,(R/'provenance/private_R_stderr_v1.txt').open('x') as e:
        child=subprocess.Popen(cmd,stdout=o,stderr=e,start_new_session=True);stat=Path('/proc')/str(child.pid)/'stat';ticks=stat.read_text().rsplit(')',1)[1].split()[19];save(R/'provenance/R_child_launch_v1.json',dict(PID=child.pid,startticks=ticks))
        try:code=child.wait(timeout=cfg['R_wall_seconds'])
        except subprocess.TimeoutExpired:
            assert stat.exists() and stat.read_text().rsplit(')',1)[1].split()[19]==ticks;os.killpg(child.pid,signal.SIGTERM);code=child.wait(timeout=30)
        save(R/'provenance/R_child_terminal_v1.json',dict(actual_exit_code=code))
    assert code==0,'Full source query replay failed'
    replay=pd.read_csv(target,sep='\t');assert len(replay)==4*len(probes);assert replay[['labels_equal','pruned_equal','delta_equal','query_axis_reversed']].to_numpy().all()
    storage_gate(R,cfg)
    accept=dict(schema='accepted_official_SingleR_group_reference_execution_v1',n_profiles=n_profiles,n_sources=len(sizes),n_fit=idx['n_fit'],coverage=coverage,annotation_manifest_SHA256=sha(manifest),table_SHA256=sha(S/'outputs/private_all_matched_SingleR_proposals_v1.tsv'),reference_inputs_SHA256=sha(S/'provenance/frozen_public_reference_inputs_v1.json'),independent_full_source_replay_records=len(replay),original_pruning_query_batches_preserved=True,query_axis_reverse_invariance_checked=True,profile_annotations_are_not_per_cell_mixture_resolution=True,final_CM_input_qualification_not_granted=True)
    save(R/'outputs/accepted_SingleR_group_execution_v1.json',accept)
    save(R/'provenance/terminal_v1.json',dict(status='SINGLER_GROUP_EXECUTION_FULL_AXES_AND_REVERSED_QUERY_REPLAY_ACCEPTED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(R/'outputs/accepted_SingleR_group_execution_v1.json')))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise

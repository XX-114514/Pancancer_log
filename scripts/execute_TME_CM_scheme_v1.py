"""Run one accepted CoVarNet scheme, with a separate synthetic API mode."""
import argparse,datetime,hashlib,json,os,subprocess,sys
from pathlib import Path
import numpy as np
import pandas as pd

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def synthetic_inputs(R):
    rng=np.random.default_rng(20261030);n=36;u=np.linspace(-2,2,n);values=[]
    for j in range(3):
        logits=np.stack([u+.1*j,.7*u,-u,-.5*u])+rng.normal(0,.15,(4,n));p=np.exp(logits-logits.max(axis=0));values.extend(p/p.sum(axis=0))
    samples=['ARTIFICIAL_UNIT_'+str(j) for j in range(n)];fq=pd.DataFrame(np.stack(values),index=['ARTIFICIAL_NODE_'+str(j) for j in range(12)],columns=samples);fq.insert(0,'majorCluster',['ARTIFICIAL_PARENT_'+str(j//4) for j in range(12)]);fp=R/'config/private_synthetic_frequency_v1.tsv';fq.to_csv(fp,sep='\t',index_label='subCluster')
    source=np.repeat(['ARTIFICIAL_A','ARTIFICIAL_B','ARTIFICIAL_C'],[8,12,16]);seen={};blocks=[]
    for s in source:
        pos=seen.get(s,0);seen[s]=pos+1;size={'ARTIFICIAL_A':1,'ARTIFICIAL_B':2,'ARTIFICIAL_C':3}[s];blocks.append(s+'::block_'+str(pos//size))
    mp=R/'config/private_synthetic_metadata_v1.tsv';pd.DataFrame(dict(sampleID=samples,source=source,patient_candidate_block=blocks,analysis_unit_kind='artificial')).to_csv(mp,sep='\t',index=False)
    return fp,mp
def main(R):
    cfg=load(R/'config/spec_v1.json');sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(R,cfg);synthetic=cfg['analysis_role']=='synthetic_validation'
    if synthetic:frequency,metadata=synthetic_inputs(R)
    else:
        a=load(cfg['input_acceptance']['path']);assert sha(cfg['input_acceptance']['path'])==cfg['input_acceptance']['SHA256'];assert a['schema']=='accepted_TME_CM_scheme_input_v1' and a['accepted_purpose']==cfg['analysis_role'];assert a['complete_measured_denominators_no_missing_to_zero'] and a['annotation_and_unit_masks_bound']
        for name in ['frequency','metadata']:
            b=a[name];assert sha(b['path'])==b['SHA256']
        frequency=Path(a['frequency']['path']);metadata=Path(a['metadata']['path'])
    request=dict(cfg['R_parameters'],frequency=str(frequency),metadata=str(metadata),author_code=str(R/'scripts/frozen_author_CoVarNet_v1.R'),output=str(R/'outputs'),analysis_role=cfg['analysis_role'],synthetic_counterexamples=synthetic)
    assert sha(request['author_code'])==cfg['executed_author_code_SHA256'];save(R/'config/frozen_R_request_v1.json',request);input_hash={str(p):sha(p) for p in [frequency,metadata]};cmd=[cfg['Rscript'],'--vanilla',str(R/'scripts/run_TME_CoVarNet_scheme_v1.R'),str(R/'config/frozen_R_request_v1.json')]
    with (R/'provenance/private_R_stdout_v1.txt').open('x') as stdout,(R/'provenance/private_R_stderr_v1.txt').open('x') as stderr:
        child=subprocess.Popen(cmd,stdout=stdout,stderr=stderr,start_new_session=True);p=Path('/proc')/str(child.pid)/'stat';ticks=p.read_text().rsplit(')',1)[1].split()[19];save(R/'provenance/R_child_launch_v1.json',dict(PID=child.pid,startticks=ticks))
        try:code=child.wait(timeout=cfg['R_wall_seconds'])
        except subprocess.TimeoutExpired:
            import signal
            assert p.exists() and p.read_text().rsplit(')',1)[1].split()[19]==ticks;os.killpg(child.pid,signal.SIGTERM);code=child.wait(timeout=30)
        save(R/'provenance/R_child_terminal_v1.json',dict(actual_exit_code=code))
    assert code==0,'CoVarNet R execution failed; preserve actual error history'
    report=load(R/'outputs/CoVarNet_execution_report_v1.json');assert report['NMF_version']==cfg['expected_NMF_version'];fq=pd.read_csv(frequency,sep='\t',index_col=0);meta=pd.read_csv(metadata,sep='\t');assert report['n_nodes']==len(fq) and report['n_analysis_observations']==len(meta);weights=pd.read_csv(R/'outputs/private_analysis_weights_v1.tsv',sep='\t');assert weights.sampleID.tolist()==fq.columns[1:].tolist();assert np.isclose(weights.weight.mean(),1)
    for seed in request['final_seeds']:
        root=R/'outputs'/('seed_'+str(seed));frames={name:pd.read_csv(root/('private_'+name+'_v1.tsv'),sep='\t',index_col=0) for name in ['W','S','A','H','H_weighted']};assert frames['W'].index.tolist()==fq.index.tolist() and frames['A'].index.tolist()==fq.index.tolist();assert frames['H'].columns.tolist()==fq.columns[1:].tolist() and frames['H_weighted'].columns.tolist()==fq.columns[1:].tolist();parts={name:f.to_numpy(float) for name,f in frames.items()}
        np.testing.assert_allclose(parts['W']@parts['S'],parts['A'],rtol=1e-9,atol=1e-10);np.testing.assert_allclose(parts['H']*weights.weight.to_numpy()[None,:],parts['H_weighted'],rtol=1e-9,atol=1e-10);assert parts['A'].shape==(len(fq),report['selected_rank']) and parts['H'].shape==(report['selected_rank'],len(meta))
        for x in parts.values():assert np.isfinite(x).all() and (x>=-1e-12).all()
    assert all(sha(p)==digest for p,digest in input_hash.items());storage_gate(R,cfg)
    save(R/'outputs/numeric_readback_acceptance_v1.json',dict(n_final_seeds=len(request['final_seeds']),all_saved_W_S_effective_A_and_weighted_unweighted_H_checked=True,all_input_files_unchanged=True,synthetic_not_human_discovery=synthetic,report_SHA256=sha(R/'outputs/CoVarNet_execution_report_v1.json')))
    save(R/'provenance/terminal_v1.json',dict(status='AUTHOR_COVARNET_REAL_NSNMF_SYNTHETIC_API_CHECK_PASS' if synthetic else 'ACCEPTED_SCHEME_COVARNET_NMF_AND_NETWORK_NUMERIC_EXECUTION_COMPLETE',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_nodes=report['n_nodes'],n_observations=report['n_analysis_observations'],selected_rank=report['selected_rank'],synthetic_only=synthetic,acceptance_SHA256=sha(R/'outputs/numeric_readback_acceptance_v1.json')))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise

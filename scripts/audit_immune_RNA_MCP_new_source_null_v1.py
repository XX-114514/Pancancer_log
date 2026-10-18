"""Independent QR replay of every finite null dataset and wild stream."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from immune_RNA_MCP_new_source_io_v1 import read,unique_designs,wait
from immune_RNA_MCP_new_source_math_v1 import design,independent_design,independent_wild


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2);out=root/'immune_RNA_MCP_new_source_family_null_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700);source=root/'immune_RNA_MCP_new_source_family_null_check_20261019_v1'
    terminal=wait(root,out,'calibrate_immune_RNA_MCP_new_source_family_v1',source/'private_terminal.json',{'COMPLETE_ALL_UNIQUE_NULL_DESIGNS_WITH_PREDECLARED_PER_DESIGN_FLAGS_NO_TRUE_EFFECTS'})
    for name,v in terminal['files'].items():assert sha(source/name)==v['sha256']
    _,_,_,registry,arrays=read(root);unique,mapping=unique_designs(registry,arrays);assert len(unique)==terminal['unique_designs'];stats=pd.read_csv(source/'private_all_null_simulation_statistics.tsv',sep='\t').set_index(['design_key','error_model','simulation']);summary=pd.read_csv(source/'private_all_unique_design_null_scenarios.tsv',sep='\t').set_index(['design_key','error_model']);savedmap=pd.read_csv(source/'private_complete_original_model_design_mapping.tsv',sep='\t').set_index('model_id')
    rng=np.random.default_rng(2026101905);datasets=0;wilds=0;maxerror=0.;flags={}
    for number,(key,row) in enumerate(unique):
        a=arrays[row.model_id];m,status=design(a['x'],a['y'],a['C'],a['names'],row.kind);ind,st=independent_design(a['x'],a['y'],a['C'],a['names'],row.kind);assert status==st=='ESTIMABLE';D=m['D'];n=len(D);mu=D[:,:-1]@np.linspace(-.2,.2,D.shape[1]-1);q=a['C'][:,0];q=(q-q.mean())/q.std(ddof=0) if q.std(ddof=0)>0 else np.zeros(n);passes=True
        for ei,error in enumerate(['normal','t3_unit_variance','RNA_quality_heteroskedastic_normal']):
            rejected=0
            for simulation in range(256):
                e=rng.standard_t(3,n)/np.sqrt(3) if error=='t3_unit_variance' else rng.normal(size=n)
                if error=='RNA_quality_heteroskedastic_normal':e*=np.exp(.6*q)
                y=mu+e;s=2026101905+number*100000+ei*10000+simulation;g,values=independent_wild(ind,999,s,y=y);count=int((np.abs(values)>=abs(g['t'])-1e-12).sum());saved=stats.loc[(key,error,simulation)];assert saved.wild_seed==s and saved.wild_exceedances==count and saved.representative_model==row.model_id;err=max(abs(saved.beta-g['beta']),abs(saved.HC3_t-g['t']));assert err<1e-8;maxerror=max(maxerror,err);rejected+=(1+count)/1000<.05;datasets+=1;wilds+=999
            r=summary.loc[(key,error)];assert r.rejections==rejected and r.simulations==256 and abs(r.rejection_rate-rejected/256)<1e-12;passes=passes and rejected/256<=.10
            (out/'private_progress.json').write_text(json.dumps(dict(phase='INDEPENDENT_NULL_STREAM_REPLAY',datasets_replayed=datasets,wild_statistics_regenerated=wilds,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
        flags[key]=passes
    for model,key in mapping.items():assert savedmap.loc[model,'design_key']==key and savedmap.loc[model,'pass_all_null_scenarios']==flags[key]
    assert datasets==terminal['simulated_datasets']==len(stats) and sum(not x for x in flags.values())==terminal['failed_unique_designs']
    receipt=dict(status='ACCEPTED_ALL_UNIQUE_NULL_SIMULATIONS_RANDOM_STREAMS_AND_PER_DESIGN_FLAGS_NO_TRUE_EFFECTS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),resource_guard=guard,unique_designs=len(unique),all_scenarios=len(unique)*3,simulated_datasets_checked=datasets,complete_wild_statistics_regenerated=wilds,maximum_error=maxerror,failed_unique_designs=sum(not x for x in flags.values()),actual_association_effects_computed=False)
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(dict(status=receipt['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

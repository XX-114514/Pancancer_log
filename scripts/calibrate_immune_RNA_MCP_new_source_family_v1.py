"""Finite per-design null check; no observed association fits."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from HNSCC_author_state_clinical_kernel_v2 import fit
from immune_RNA_MCP_new_source_io_v1 import read,unique_designs
from immune_RNA_MCP_new_source_math_v1 import design,wild_values

ERRORS=['normal','t3_unit_variance','RNA_quality_heteroskedastic_normal']


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2);_,audit,_,registry,arrays=read(root)
    math=root/'immune_RNA_MCP_new_source_math_check_20261019_v1/private_terminal.json';assert sha(math)=='8066936bd66e2b65fb811204dbbb15cb55a93c183e9ddf6005c1228072777bb7'
    unique,mapping=unique_designs(registry,arrays);out=root/'immune_RNA_MCP_new_source_family_null_check_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,input_acceptance_sha256=sha(audit),math_acceptance_sha256=sha(math),unique_designs=len(unique),scope='Every distinct eligible adjusted-source or clinical design, no actual association estimates. Marginal rank tests use direct case permutations and do not require this working-model screen.',simulations=256,wild_per_simulation=999,errors=ERRORS,seed=2026101905,gate='Per unique design every3 working-error rejection rates at nominal.05 <=.10. Failing designs stay in preeligible BH denominator with P1 and no adjusted inference. No threshold/covariate retry. Other passed designs and exact marginal questions continue.',limits='Finite working-model anomaly diagnostic, not a guarantee under all misspecification. Synthetic nuisance mean only; target coefficient0.')
    (out/'plan_before_all_unique_design_null_checks.json').write_text(json.dumps(plan,indent=2)+'\n');rng=np.random.default_rng(2026101905);rows=[];summary=[]
    for number,(key,row) in enumerate(unique):
        a=arrays[row.model_id];m,status=design(a['x'],a['y'],a['C'],a['names'],row.kind);assert status=='ESTIMABLE';D=m['D'];n=len(D);mu=D[:,:-1]@np.linspace(-.2,.2,D.shape[1]-1);q=a['C'][:,0];q=(q-q.mean())/q.std(ddof=0) if q.std(ddof=0)>0 else np.zeros(n)
        for ei,error in enumerate(ERRORS):
            rejected=0
            for simulation in range(256):
                e=rng.standard_t(3,n)/np.sqrt(3) if error=='t3_unit_variance' else rng.normal(size=n)
                if error=='RNA_quality_heteroskedastic_normal':e*=np.exp(.6*q)
                g=fit(mu+e,D);s=2026101905+number*100000+ei*10000+simulation;values=wild_values(g,999,s);count=int((np.abs(values)>=abs(g['t'])-1e-12).sum());rejected+=(1+count)/1000<.05
                rows.append(dict(design_key=key,representative_model=row.model_id,error_model=error,simulation=simulation,wild_seed=s,beta=g['beta'],HC3_t=g['t'],wild_exceedances=count))
            summary.append(dict(design_key=key,representative_model=row.model_id,error_model=error,simulations=256,rejections=int(rejected),rejection_rate=rejected/256))
            (out/'private_progress.json').write_text(json.dumps(dict(phase='NULL_WORKING_MODELS_ONLY',completed_scenarios=len(summary),total_scenarios=len(unique)*3,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    pd.DataFrame(rows).to_csv(out/'private_all_null_simulation_statistics.tsv',sep='\t',index=False);table=pd.DataFrame(summary);table.to_csv(out/'private_all_unique_design_null_scenarios.tsv',sep='\t',index=False)
    passes=table.groupby('design_key').rejection_rate.max().le(.10).to_dict();pd.DataFrame([dict(model_id=model,design_key=key,pass_all_null_scenarios=passes[key]) for model,key in mapping.items()]).to_csv(out/'private_complete_original_model_design_mapping.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<8*1024**2
    terminal=dict(status='COMPLETE_ALL_UNIQUE_NULL_DESIGNS_WITH_PREDECLARED_PER_DESIGN_FLAGS_NO_TRUE_EFFECTS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,unique_designs=len(unique),simulated_datasets=len(rows),failed_unique_designs=sum(not x for x in passes.values()),actual_association_effects_computed=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(dict(status=terminal['status'],terminal_sha256=sha(out/'private_terminal.json'),unique_designs=len(unique),failed_designs=terminal['failed_unique_designs'])))


if __name__=='__main__':main()

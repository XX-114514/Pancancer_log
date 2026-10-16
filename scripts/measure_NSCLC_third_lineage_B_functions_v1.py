"""All frozen-mask functions in the complete author B lineage, existing cache only."""
import hashlib,json,os,subprocess
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from disjoint_hallmark_ULM_kernel_v1 import score


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);pre=r/'NSCLC_disjoint_hallmark_preflight_20261016_v1';native=r/'NSCLC_full_native_gene_units_20261016_v1';qcroot=r/'NSCLC_author_RNA_composition_20261016_v1'
    question=r/'NSCLC_third_lineage_background_plan_20261016_v1/PRIVATE_NEXT_QUESTION_AND_SCOPE.json';frozen=json.loads(question.read_text())
    proofs=dict(matched_function_audit=r/'NSCLC_matched_function_coordination_20261016_v1/private_independent_matched_function_audit.json',RNA_QC_audit=r/'NSCLC_RNA_composition_sensitivity_20261016_v1/private_independent_sensitivity_audit.json',native_unit_audit=r/'NSCLC_full_native_gene_unit_audit_20261016_v1/private_terminal.json')
    for k,p in proofs.items():assert sha(p)==frozen['input_sha256'][k] and json.loads(p.read_text())['status'].startswith('PASS_')
    paths=dict(question=question,units=native/'private_all_native_units_and_quality.tsv',matrix=native/'mean_logCP10K.npy',QC=qcroot/'private_all_native_unit_RNA_composition.tsv',masks=pre/'private_physical_column_masks.json',references=pre/'private_reference_axis.tsv',coverage=pre/'private_all1200_functional_scope.tsv',official=r/'fullgene_receiver_target_code_v1/collectri_reference_v1/official_decoupler192_ULM_source.txt')
    hashes={k:sha(p) for k,p in paths.items()};old=json.loads((native/'private_terminal.json').read_text());original=json.loads((pre/'plan_before_functional_RNA_measurements.json').read_text())
    for k,name in [('matrix','mean_logCP10K.npy'),('units','private_all_native_units_and_quality.tsv')]:assert hashes[k]==old['outputs'][name]['sha256']
    assert hashes['masks']==frozen['input_sha256']['frozen_masks'] and hashes['official']==original['input_sha256']['official_ULM']
    assert hashes['QC']==json.loads((qcroot/'private_terminal.json').read_text())['files'][paths['QC'].name]
    kernel=Path(__file__).with_name('disjoint_hallmark_ULM_kernel_v1.py');assert sha(kernel)==original['kernel_sha256']
    df=subprocess.run(['df','-P',str(r.resolve())],capture_output=True,text=True,timeout=20);quota=subprocess.run(['quota','-g','-v','-w','USER002'],capture_output=True,text=True,timeout=20);assert df.returncode==quota.returncode==0
    fs=df.stdout.splitlines()[-1].split()[0];qr=[x.split() for x in quota.stdout.splitlines() if x.split() and x.split()[0]==fs];assert len(qr)==1;remaining=(int(qr[0][3])-int(qr[0][1].rstrip('*')))*1024
    assert remaining>500_000_000_000+1_073_741_824+268_435_456
    ram=int(next(x for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')).split()[1])*1024;assert ram>2*1024**3
    out=r/'NSCLC_third_lineage_B_functions_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),kernel_sha256=sha(kernel),input_sha256=hashes,input_audit_sha256={k:sha(p) for k,p in proofs.items()},
              scope='Every author B cell subtype pooled within every original source case; all16 original context/MCP physical masks x all50 references. No clinical or association filter.',
              pooling='Original per-cell log1pCP10K means averaged with native unit cell counts; preserve full original RNA denominator. RAM only pooled expression, no expression copy on disk.',
              quality='Cell-count-weighted original unit mean logRNA/loggenes and author mito/ribo percentages; B fraction of ALL original source cells per author case and log1p B cell count.',
              measurement='Same pinned ULM1.9.2 and BOTH-partner learned-gene exclusion in targets/background; preserve measured zero and RNA target support. Not protein activity or state label.',
              future_association='All572 pairs x50functions; compare before/after B function on identical explicitly B-captured cases; no new P/BH or original family changes.',
              reserve_decimal_bytes=500_000_000_000,quota_free_bytes=remaining,RAM_available_bytes=ram,output_cap_bytes=268_435_456)
    (out/'plan_before_any_B_function_measurement.json').write_text(json.dumps(plan,indent=2)+'\n')
    u=pd.read_csv(paths['units'],sep='\t').sort_values('unit_index');qc=pd.read_csv(paths['QC'],sep='\t').set_index('unit_index');b=u[u.major_cell_type.eq('B cell')].copy();assert len(b)>0 and b.unit_index.is_unique
    assert b.sub_cell_type.nunique()==10 and set(b.sampleID)<=set(u.sampleID)
    totals=u.groupby('sampleID').n_cells.sum();q=qc.reindex(b.unit_index);assert np.array_equal(q.n_cells.to_numpy(int),b.n_cells.to_numpy(int))
    b['mean_author_pct_mito']=q.mean_author_pct_mito.to_numpy(float);b['mean_author_pct_ribo']=q.mean_author_pct_ribo.to_numpy(float)
    cases=sorted(b.sampleID.unique());m=np.load(paths['matrix'],mmap_mode='r');assert m.shape[0]==len(u);pooled=np.empty((len(cases),m.shape[1]),float);axis=[];mapping=[]
    for i,case in enumerate(cases):
        f=b[b.sampleID.eq(case)];n=int(f.n_cells.sum());weights=f.n_cells.to_numpy(float)/n;assert n>0
        pooled[i]=weights@np.asarray(m[f.unit_index.to_numpy(int),:]);row=dict(pool_index=i,sampleID=case,n_cells=n,author_subtypes=int(f.sub_cell_type.nunique()),fraction=n/int(totals.loc[case]),logN=float(np.log1p(n)))
        for field in ['mean_logRNA','mean_loggenes','mean_author_pct_mito','mean_author_pct_ribo']:row[field]=float(weights@f[field].to_numpy(float))
        assert np.isfinite(list(row[k] for k in ['fraction','logN','mean_logRNA','mean_loggenes','mean_author_pct_mito','mean_author_pct_ribo'])).all();axis.append(row)
        for unit,weight in zip(f.unit_index,weights):mapping.append(dict(pool_index=i,unit_index=int(unit),cell_weight=float(weight)))
    assert np.isfinite(pooled).all() and (pooled>=0).all()
    pd.DataFrame(axis).to_csv(out/'private_all_B_case_axis_and_quality.tsv',sep='\t',index=False);pd.DataFrame(mapping).to_csv(out/'private_all_B_native_unit_pool_mapping.tsv',sep='\t',index=False)
    masks=json.loads(paths['masks'].read_text());refs=pd.read_csv(paths['references'],sep='\t').sort_values('reference_column');cov=pd.read_csv(paths['coverage'],sep='\t');mask_names=sorted(masks);assert len(mask_names)==16 and len(refs)==50
    arrays={k:np.full((len(mask_names),len(cases),50),np.nan,float) for k in ['ULM_score','target_mean','nontarget_mean']};arrays['expressed_target_genes']=np.full((len(mask_names),len(cases),50),-1,np.int16)
    for i,name in enumerate(mask_names):
        mask=masks[name];keep=np.ones(m.shape[1],bool);keep[mask['forbidden_native_columns']]=False;bg=np.flatnonzero(keep);reverse=np.full(m.shape[1],-1,int);reverse[bg]=np.arange(len(bg))
        input_rows=cov[cov.context.eq(mask['context'])&cov.program.eq(mask['program'])];assert input_rows.status.eq('MEASURABLE_INPUT_NO_RNA_SCORE').all() and input_rows.reference.nunique()==50
        targets=[reverse[np.asarray(mask['targets'][ref],int)] for ref in refs.reference];assert all((x>=0).all() and len(x)>=20 and len(bg)-len(x)>=500 for x in targets)
        for start in range(0,len(cases),64):
            stop=min(start+64,len(cases));result=score(pooled[start:stop,bg],targets)
            for key,value in result.items():arrays[key][i,start:stop,:]=value
        (out/'private_measurement_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),masks_completed=i+1,total_masks=len(mask_names)))+'\n')
    pd.DataFrame(dict(mask_index=np.arange(len(mask_names)),mask_id=mask_names)).to_csv(out/'private_all16_mask_axis.tsv',sep='\t',index=False)
    np.savez_compressed(out/'private_all_B_function_scores_and_QA.npz',**arrays);assert sum(p.stat().st_size for p in out.iterdir())<plan['output_cap_bytes']
    terminal=dict(status='COMPLETE_ALL_B_CASES_ALL16_MASKS50_FUNCTIONS_AUDIT_PENDING_NO_ASSOCIATIONS',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),native_B_units=len(b),native_B_cells=int(b.n_cells.sum()),author_B_cases=len(cases),mask_count=len(mask_names),functional_score_values=arrays['ULM_score'].size,files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_P_BH=False,new_download=False,new_expression_copy=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()

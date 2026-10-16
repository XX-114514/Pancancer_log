"""Measure every held-out-gene function in all frozen longitudinal HNSCC pools."""
import hashlib,json,os,subprocess
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.sparse import csr_matrix
from HNSCC_native_multicell_function_inputs_v1 import inputs
from disjoint_hallmark_ULM_kernel_v1 import score


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);data=inputs(r);out=r/'HNSCC_native_multicell_functions_20261016_v1';assert not out.exists()
    filesystem=subprocess.check_output(['df','-P',str(r.resolve())],text=True).splitlines()[-1].split()[0];quota=subprocess.check_output(['quota','-g','-v','-w','USER002'],text=True);qr=[s.split() for s in quota.splitlines() if s.split() and s.split()[0]==filesystem];assert len(qr)==1;free=(int(qr[0][3])-int(qr[0][1].rstrip('*')))*1024;assert free>500_000_000_000+1_073_741_824+268_435_456
    out.mkdir(mode=0o700);axis=data['axis'];mapping=data['mapping'];counts=data['zarr']['layers']['counts'];N,G=tuple(counts.attrs['shape']);U=len(axis)
    kernel=Path(__file__).with_name('disjoint_hallmark_ULM_kernel_v1.py');old=json.loads((r/'NSCLC_disjoint_hallmark_preflight_20261016_v1/plan_before_functional_RNA_measurements.json').read_text());assert sha(kernel)==old['kernel_sha256'] and data['input_sha256']['official']==old['input_sha256']['official_ULM']
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),input_helper_sha256=sha(Path(__file__).with_name('HNSCC_native_multicell_function_inputs_v1.py')),kernel_sha256=sha(kernel),input_sha256=data['input_sha256'],coverage_input_sha256=data['coverage_input_sha256'],selection='HNSCC selected by reviewed longitudinal patient mapping, >=20 complete paired-capture patients and full native functional coverage, without new effects. Other61-source qualification retained. All57 libraries and all4 disjoint frozen CD8/Macro/Treg/B+Plasma pools, including0 and low capture.',question='Do reviewed pre/post patients show coordinated functions across immune pools, and which changes are common versus endpoint limited? Measure all16 original learned-gene exclusion masks and all50 functions before effects. The eight masks from CD8-Macro/Macro-Treg contexts are relevant primary contexts; others remain full sensitivity metadata, not extra independent observations.',scope='All16 masks x all228 pool-library positions x50 references; no cell subsampling. Original native counts, full source RNA denominator, per-cell log1pCP10K mean pooled within original library/pool. No complete expression matrix persisted.',quality='Per-cell mean logRNA/loggenes and symbol-prefix-defined MT-/RPL/RPS RNA percentages, plus pool cell count and fraction of all source library cells; not author-reported percentages.',limits='Raw RNA functional pattern, not protein activity, exact old172 subtype, direct communication or baseline prediction. Recurrence remains post2 and is not ordinary post. No new P/BH or clinical model in measurement stage; original HNSCC programme clinical family closed.',quota_free_bytes=free,output_cap_bytes=268_435_456,new_download=False)
    (out/'plan_before_any_multicell_function_measurement.json').write_text(json.dumps(plan,indent=2)+'\n')
    means=np.zeros((U,G),float);quality=np.zeros((U,7),float);mito=np.char.startswith(data['symbols'],'MT-');ribo=np.char.startswith(data['symbols'],'RPL')|np.char.startswith(data['symbols'],'RPS')
    for start in range(0,N,4096):
        stop=min(start+4096,N);keep=mapping[start:stop]>=0
        if not keep.any():continue
        ptr=np.asarray(counts['indptr'][start:stop+1],np.int64);x=csr_matrix((np.asarray(counts['data'][int(ptr[0]):int(ptr[-1])]),np.asarray(counts['indices'][int(ptr[0]):int(ptr[-1])]),ptr-ptr[0]),shape=(stop-start,G))[keep].copy();ids=mapping[start:stop][keep]
        x.sum_duplicates();x.eliminate_zeros();x.sort_indices();assert (x.data>0).all() and np.isfinite(x.data).all() and np.equal(x.data,np.floor(x.data)).all() and data['native_mask'][x.indices].all()
        total=np.asarray(x.sum(axis=1)).ravel().astype(float);detected=np.diff(x.indptr);assert (total>0).all();logged=x.astype(float);logged.data=np.log1p(logged.data*np.repeat(10000/total,detected));P=csr_matrix((np.ones(len(ids)),(ids,np.arange(len(ids)))),shape=(U,len(ids)));summed=(P@logged).tocoo();means[summed.row,summed.col]+=summed.data
        values=[np.ones(len(ids)),total,np.log1p(total),np.log1p(detected),100*np.asarray(x[:,mito].sum(axis=1)).ravel()/total,100*np.asarray(x[:,ribo].sum(axis=1)).ravel()/total,detected]
        for j,value in enumerate(values):quality[:,j]+=np.bincount(ids,weights=value,minlength=U)
        (out/'private_measurement_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),last_original_row=stop,selected_cells_accumulated=int(quality[:,0].sum())))+'\n')
    assert np.array_equal(quality[:,0],axis.n_cells);valid=np.flatnonzero(axis.n_cells.to_numpy()>0);means[valid]/=quality[valid,:1]
    for j,name in enumerate(['mean_full_RNA','mean_logRNA','mean_loggenes','mean_symbol_pct_mito','mean_symbol_pct_ribo','mean_detected_genes'],start=1):axis[name]=np.divide(quality[:,j],quality[:,0],out=np.full(U,np.nan),where=quality[:,0]>0)
    names=data['masks'].sort_values('mask_id').reset_index(drop=True);references=sorted(data['references']);arrays={k:np.full((16,U,50),np.nan) for k in ['ULM_score','target_mean','nontarget_mean']};arrays['expressed_target_genes']=np.full((16,U,50),-1,np.int16);coverage=data['coverage'].set_index(['mask_id','reference']);axis_native=data['physical'];params=data['parameters'];physical_masks={}
    for i,mask in enumerate(names.itertuples()):
        learned=set(params.loc[params.context.eq(mask.context)&params.program.eq(mask.program),'gene']);forbidden=set()
        for g in learned:forbidden.update(axis_native['lookup'].get(g,set()))
        bg=np.asarray(sorted(axis_native['native']-forbidden),int);reverse={c:j for j,c in enumerate(bg)};targets=[]
        for reference in references:
            remaining=data['references'][reference]-learned;gene_to_col={g:next(iter(axis_native['lookup'][g])) for g in remaining if len(axis_native['lookup'].get(g,set()))==1};hist=Counter(gene_to_col.values());physical=sorted(c for c in gene_to_col.values() if hist[c]==1 and c not in forbidden);row=coverage.loc[(mask.mask_id,reference)];assert len(physical)==row.native_unique_targets and len(bg)==row.native_background_genes;targets.append(np.asarray([reverse[c] for c in physical],int))
        result=score(means[np.ix_(valid,bg)],targets)
        for key,value in result.items():arrays[key][i,valid,:]=value
        physical_masks[str(mask.mask_id)]={'background':bg.tolist(),'targets':{ref:bg[ix].tolist() for ref,ix in zip(references,targets)}}
    axis.to_csv(out/'private_all228_pool_library_axis_and_quality.tsv',sep='\t',index=False,na_rep='NA');names.to_csv(out/'private_all16_mask_axis.tsv',sep='\t',index=False);pd.DataFrame(dict(reference_index=np.arange(50),reference=references)).to_csv(out/'private_all50_reference_axis.tsv',sep='\t',index=False)
    np.savez_compressed(out/'private_all_multicell_function_scores_and_QA.npz',**arrays);(out/'private_physical_function_masks.json').write_text(json.dumps(physical_masks,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in out.iterdir())<plan['output_cap_bytes']
    terminal=dict(status='COMPLETE_ALL_NATIVE_HNSCC_POOL_LIBRARY_MASK_FUNCTION_MEASUREMENTS_AUDIT_PENDING_NO_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),pool_library_slots=U,selected_cells=int(quality[:,0].sum()),function_slots=16*U*50,zero_cell_pool_slots=int((axis.n_cells==0).sum()),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},full_expression_persisted=False,new_P_BH=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()

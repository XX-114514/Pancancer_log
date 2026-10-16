"""Independent count regrouping, Python TMM, QR regression and full feature audit."""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as stream:
        for b in iter(lambda:stream.read(16777216),b''):h.update(b)
    return h.hexdigest()


def tmm(obs,ref):
    no=float(obs.sum());nr=float(ref.sum());ok=(obs>0)&(ref>0);o=obs[ok].astype(float);r=ref[ok].astype(float)
    logratio=np.log2((o/no)/(r/nr));average=.5*(np.log2(o/no)+np.log2(r/nr));variance=(no-o)/no/o+(nr-r)/nr/r;finite=np.isfinite(logratio)&np.isfinite(average)&(average>-1e10);logratio=logratio[finite];average=average[finite];variance=variance[finite]
    if not len(logratio) or np.max(abs(logratio))<1e-6:return 1.
    n=len(logratio);lo_l=np.floor(n*.3)+1;hi_l=n+1-lo_l;lo_a=np.floor(n*.05)+1;hi_a=n+1-lo_a;rl=rankdata(logratio);ra=rankdata(average);keep=(rl>=lo_l)&(rl<=hi_l)&(ra>=lo_a)&(ra<=hi_a);weight=1/variance[keep];value=float(np.sum(logratio[keep]*weight)/np.sum(weight));return 2**value if np.isfinite(value) else 1.


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=root/'NSCLC_native_MFA_inputs_20261016_v1';job=root/'HNSCC_receiver_followup_job_control_20261016_v2/prepare_NSCLC_native_MFA_inputs_v1';receipt=json.loads((job/'launch_receipt.json').read_text());deadline=time.monotonic()+3600
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split();assert s[0]!='Z' and s[19]==str(receipt['proc_start_ticks']) and time.monotonic()<deadline;time.sleep(3)
    terminal=json.loads((out/'private_terminal.json').read_text());assert terminal['code_sha256']==receipt['source_sha256']
    for name,h in terminal['files'].items():assert sha(out/name)==h
    plan=json.loads((out/'plan_before_native_expression_and_MFA_fit.json').read_text());qual=root/'NSCLC_native_multicell_discovery_preflight_20261016_v1';previous=json.loads((qual/'plan_before_any_native_discovery_effect.json').read_text());roles=pd.read_csv(qual/'private_all243_author_case_therapy_partition.tsv',sep='\t',keep_default_na=False,dtype={'sampleID':str});axis=pd.read_csv(out/'private_case_axis.tsv',sep='\t',keep_default_na=False,dtype={'sampleID':str});expected=roles[roles.discovery_role.isin(['discovery','internal_holdout'])].sort_values('sampleID').rename(columns={'discovery_role':'partition'});assert axis.sampleID.tolist()==expected.sampleID.tolist() and axis.partition.tolist()==expected.partition.tolist() and axis.therapy_original.tolist()==expected.therapy_original.tolist() and axis.cancer_type.tolist()==expected.cancer_type.tolist()
    drugs=['Tislelizumab','Pembrolizumab','Camrelizumab','Sintilimab','Nivolumab','Toripalimab']
    for row in axis.itertuples():
        matches={x for x in drugs if x in row.therapy_original};group=next(iter(matches)) if len(matches)==1 else 'MULTI_OR_UNSPECIFIED';assert row.therapy_group==group and row.holdout_block==row.cancer_type+'::'+group
    train=axis.partition.eq('discovery').to_numpy();cases=axis.sampleID.tolist();assert len(cases)==208 and train.sum()==145
    native=root/'NSCLC_full_native_gene_units_20261016_v1';units=pd.read_csv(native/'private_all_native_units_and_quality.tsv',sep='\t',dtype={'sampleID':str});genes=pd.read_csv(native/'private_native_gene_axis.tsv',sep='\t');raw=np.load(native/'sum_raw_counts.npy',mmap_mode='r');mt=genes.gene.str.startswith('MT-').to_numpy();rb=genes.gene.str.match(r'^(RPL|RPS)').to_numpy();qa=pd.read_csv(out/'private_all832_case_view_quality.tsv',sep='\t',dtype={'sampleID':str}).set_index(['sampleID','view']);summaries=pd.read_csv(out/'private_all_four_view_input_summary.tsv',sep='\t').set_index('view');maximum=0.
    for view in plan['primary_views']:
        saved=pd.read_csv(out/('private_'+view+'_pseudobulk.tsv.gz'),sep='\t',index_col=0);assert saved.columns.tolist()==cases and saved.index.astype(str).tolist()==genes.gene.astype(str).tolist();counts=saved.to_numpy().T;assert np.array_equal(counts,counts.astype(np.int64));covariates=[];parent='T/NK cell' if view in ['author_CD8','author_Treg'] else 'Myeloid cell' if view=='author_Macrophage' else 'B cell'
        for i,case in enumerate(cases):
            rows=units[units.sampleID.eq(case)&units.sub_cell_type.isin(previous['view_author_labels'][view])];n=int(rows.n_cells.sum());expected_counts=np.zeros(len(genes),np.int64)
            for index in rows.unit_index.to_numpy(int):expected_counts+=raw[index]
            assert np.array_equal(counts[i],expected_counts);total=int(expected_counts.sum());parent_n=int(units.loc[units.sampleID.eq(case)&units.major_cell_type.eq(parent),'n_cells'].sum());q=[np.log1p(n),float((rows.mean_logRNA*rows.n_cells).sum()/n),float((rows.mean_loggenes*rows.n_cells).sum()/n),100*expected_counts[mt].sum()/total,100*expected_counts[rb].sum()/total]
            if view!='author_B_Plasma':q.append(n/parent_n)
            covariates.append(q);one=qa.loc[(case,view)];assert one.n_cells==n and one.parent_cells==parent_n and one.full_native_RNA==total
        norm=pd.read_csv(out/('private_'+view+'_TMM.tsv'),sep='\t',dtype={'sampleID':str});assert norm.sampleID.tolist()==cases;reference=counts[train].sum(axis=0);factors=np.array([tmm(x,reference) for x in counts]);maximum=max(maximum,float(np.max(abs(factors-norm.raw_factor))));assert np.allclose(factors,norm.raw_factor,rtol=1e-10,atol=1e-12);factor=factors/np.exp(np.log(factors[train]).mean());lib=counts.sum(axis=1);assert np.array_equal(lib,norm.library_size) and np.allclose(lib*factor,norm.effective_library_size,rtol=1e-10,atol=1e-6)
        y=np.log2(1+counts*1e6/(lib*factor)[:,None]);cov=np.array(covariates);mean=cov[train].mean(axis=0);sd=cov[train].std(axis=0);d=np.column_stack([np.ones(len(cov)),(cov-mean)/sd]);q,rr=np.linalg.qr(d[train],mode='reduced');coef=np.linalg.solve(rr,q.T@y[train]);residual=y-d@coef;var=residual[train].var(axis=0);detected=(counts[train]>=100).sum(axis=0);eligible=(detected>=37)&~(mt|rb)&(var>1e-12);candidate=np.flatnonzero(eligible);selected=candidate[np.lexsort((genes.gene.to_numpy()[candidate],-var[candidate]))[:2000]];a=np.load(out/('private_'+view+'_prepared.npz'),allow_pickle=False);assert a['cases'].astype(str).tolist()==cases and np.array_equal(a['train'],train);assert np.array_equal(a['selected_native_columns'],selected) and np.array_equal(a['selected_genes'],genes.gene.to_numpy(str)[selected]);assert np.allclose(a['quality'],cov,rtol=1e-12,atol=1e-12) and np.allclose(a['design'],d,rtol=1e-12,atol=1e-12) and np.allclose(a['regression_coefficients'],coef[:,selected],rtol=1e-9,atol=1e-10)
        rms=np.sqrt(np.mean(residual[train][:,selected]**2));assert np.isclose(a['training_view_RMS'],rms,rtol=1e-10) and np.allclose(a['X'],residual[:,selected]/rms,rtol=1e-9,atol=1e-9)
        feature=pd.read_csv(out/('private_'+view+'_all_native_feature_selection.tsv.gz'),sep='\t');assert feature.native_column.tolist()==list(range(31831)) and np.array_equal(feature.discovery_ge100_cases,detected) and np.array_equal(feature.technical_excluded,mt|rb) and np.array_equal(feature.eligible,eligible) and np.array_equal(feature.selected,np.isin(np.arange(len(genes)),selected));assert np.allclose(feature.technical_residual_variance,var,rtol=1e-8,atol=1e-10)
        summary=summaries.loc[view];assert summary.selected_features==len(selected) and summary.eligible_features==eligible.sum() and summary.technical_covariates==len(cov[0]) and not summary.holdout_used_for_selection
    proof=dict(status='PASS_ALL832_NATIVE_COUNTS_PYTHON_EDGE_R_TMM_PARITY_FULL_QC_QR_HVG_AND_FROZEN_PARTITION_NO_LEAKAGE',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),all_native_features_each_view=31831,case_view_rows_checked=832,maximum_raw_TMM_factor_difference=maximum,native_raw_audit_inherited_not_replayed_all_source_cells=True,real_model_or_clinical_effects=False);p=out/'private_independent_MFA_input_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()

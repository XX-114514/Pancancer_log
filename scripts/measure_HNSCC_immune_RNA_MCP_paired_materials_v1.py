"""Full frozen RNA assay in all qualified paired materials; no paired effects."""
import gzip
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy.io import mmread
from scipy.sparse import csr_matrix, save_npz, vstack

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

VIEWS = ['R_FOXP3_capture', 'E_effector_capture', 'Myeloid']
STRESS = ['FOS','FOSB','JUN','JUNB','JUND','ATF3','EGR1','DUSP1','HSPA1A','HSPA1B','DNAJB1','IER2','IER3','NR4A1','BTG2']
FIELDS = ['cellQ','geneQ','fraction','logN']


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);cap=128*1024**2;guard=resource(root,cap)
    paired=root/'HNSCC_immune_RNA_MCP_paired_materials_20261019_v1'
    assert sha(paired/'private_independent_metadata_acceptance_20261019_v1.json')=='8dc0f5b74ea8cf41ff10a62eb5c804064a0dc558ec90e787f4e9cbac5b4461c9'
    scope=root/'immune_RNA_MCP_discovery_scope_20261018_v1'
    assert sha(scope/'private_terminal.json')=='0f4833a07bd5cdf5e8f6806d06805777dddf78adad60a6151eca991dfdda6196'
    for n,v in json.loads((scope/'private_terminal.json').read_text())['files'].items():assert sha(scope/n)==v['sha256']
    frozen=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1'
    assert sha(root/'immune_RNA_MCP_native_expression_input_audit_20261018_v1/private_terminal.json')=='58385bca98d751470d7f7026b0e6ff4ca7f1c584dacbfa6b7eb5ac644f7b5d4b'
    assert sha(frozen/'private_terminal.json')=='2474355a003d96232b752ad2c77ba92a506f4c68829c7e0adf434045cf03d656'
    for n,v in json.loads((frozen/'private_terminal.json').read_text())['files'].items():assert sha(frozen/n)==v['sha256']
    old=root/'immune_CM_fixed_panel_same_cell_RNA_20261018_v1'
    proof=json.loads((root/'immune_CM_fixed_panel_same_cell_RNA_audit_20261018_v1/private_terminal.json').read_text())
    assert proof['status'].startswith('ACCEPTED_ALL_ORIGINAL_T_CELLS_') and sha(old/'private_terminal.json')==proof['producer_terminal_sha256']
    for n,v in json.loads((old/'private_terminal.json').read_text())['files'].items():assert sha(old/n)==v['sha256']
    oldscope=root/'existing_immune_fine_native_measurement_scope_20261016_v1'
    for n,v in json.loads((oldscope/'private_terminal.json').read_text())['files'].items():assert sha(oldscope/n)==v
    pairs=pd.read_csv(paired/'private_all_exact_pre_post_patient_material_pairs.tsv',sep='\t');pairs=pairs[~pairs.post_used_in_original_discovery]
    assert len(pairs)==26 and pairs.patient_registry.is_unique
    materials=[]
    for r in pairs.itertuples(index=False):
        for phase in ['pre','post']:materials.append(dict(patient_registry=r.patient_registry,phase=phase,sampleID=getattr(r,phase+'_sampleID')))
    materials=pd.DataFrame(materials);assert materials.sampleID.is_unique and len(materials)==52
    overlay=pd.read_csv(scope/'private_all236_current_case_phase_overlay.tsv',sep='\t');overlay=overlay[overlay.threshold.eq(1)].set_index('sampleID')
    materials['all_material_cells']=overlay.loc[materials.sampleID,'annotated_scope_cells'].to_numpy(int)
    ready=root/'immune_CM_RNA_supported_multicell_views_preflight_20261018_v1'
    assert sha(ready/'private_independent_count_readiness_acceptance_20261018_v2.json')=='3f1261b9083f1cefcf0e072daecdef7c3db5d1fbf7f44b4ac1ff237efb126c15'
    capture=pd.read_csv(ready/'private_all_original_libraries_multicell_view_counts.tsv',sep='\t');capture=capture[capture.threshold.eq(1)].set_index('sampleID')
    lookup=pd.Index(materials.sampleID)
    cache=np.load(old/'private_fixed34_gene_T_counts_and_depth.npz');raw=cache['raw_gene_counts'];markers=cache['genes'].tolist()
    T=pd.read_csv(old/'private_original_T_unit_scope.tsv',sep='\t');tcase=lookup.get_indexer(T.dataset.astype(str)+'::'+T.library_id.astype(str))[cache['original_T_unit_index']]
    fixed=root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json';panels=json.loads(fixed.read_text())['panels']
    core=(raw[:,[markers.index(g) for g in panels['pan_T']]]>0).sum(axis=1)>=2
    reg=((raw[:,[markers.index(g) for g in panels['regulatory_associated']]]>0).sum(axis=1)>=2)&(raw[:,markers.index('FOXP3')]>0)
    eff=(raw[:,[markers.index(g) for g in ['GZMB','PRF1','CTSW','NKG7']]]>0).sum(axis=1)>=2
    units=pd.read_csv(oldscope/'private_all_frozen_T_myeloid_node_units.tsv',sep='\t')
    mapping=np.load(oldscope/'private_original_cell_to_fine_node_unit.npy',mmap_mode='r')
    unitcase=lookup.get_indexer(units.dataset.astype(str)+'::'+units.library_id.astype(str));unitcase[~units.major.eq('Myeloid')]=-1
    mr=np.flatnonzero(mapping>=0);mc=unitcase[mapping[mr]]
    pieces=[];coverage=[]
    specs=[(cache['global_original_T_row_index'][core&reg&~eff&(tcase>=0)],tcase[core&reg&~eff&(tcase>=0)]),
           (cache['global_original_T_row_index'][core&eff&~reg&(tcase>=0)],tcase[core&eff&~reg&(tcase>=0)]),(mr[mc>=0],mc[mc>=0])]
    count_fields=['FOXP3_supported_R_only','effector_only','original_myeloid_cells']
    for vi,(rows,cases) in enumerate(specs):
        for m in range(len(materials)):
            rows_m=rows[cases==m];n=len(rows_m);assert n==int(capture.loc[materials.iloc[m].sampleID,count_fields[vi]]) and n>=20
            chosen=set(sorted(rows_m,key=lambda r:hashlib.sha256(('2026101802|'+VIEWS[vi]+'|'+str(int(r))).encode()).hexdigest())[:100])
            pieces.append(pd.DataFrame(dict(global_row=rows_m,material_index=m,view_index=vi,selected=[int(r) in chosen for r in rows_m])))
            coverage.append(dict(material_index=m,view=VIEWS[vi],full_view_cells=n,selected_cells=len(chosen),all_material_cells=int(materials.iloc[m].all_material_cells)))
    cells=pd.concat(pieces,ignore_index=True).sort_values('global_row').reset_index(drop=True);assert cells.global_row.is_unique
    cov=pd.DataFrame(coverage);assert len(cov)==156
    out=root/'HNSCC_immune_RNA_MCP_paired_native_measurement_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,paired_acceptance_sha256=sha(paired/'private_independent_metadata_acceptance_20261019_v1.json'),
        scope='All26 exact reviewed pre/post pairs; all52 materials/3views/2frozen programmes. No outcomes, phase differences, correlations or P read/computed.',
        selection='Original complete RNA gate and whole originalMyeloid; allfull capture retained; primary smallest100 SHA256(2026101802|view|global_row), all-capture mean sensitivity.',
        normalisation='Full original RNA and detected genes from native counts; original discovery mean/SD/weights and nuisance coefficients. HNSCC sourceCC=0. No gene/quality refit using new post cells.',
        fixed_stress_complement=STRESS,stress_meaning='Remove direct fixed-gene contribution and its original training nuisance component, not a newly learned or stress-free programme.',
        measurements_only=True,output_cap_bytes=cap,no_new_download=True)
    (out/'plan_before_paired_native_expression.json').write_text(json.dumps(plan,indent=2)+'\n')
    weights_path=root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv'
    assert sha(weights_path)=='fd7316c67beb687c403a752c36b5989644e432f59a83119b21eade5efb0d89a7'
    weights=pd.read_csv(weights_path,sep='\t');weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')]
    gm=pd.read_csv(scope/'private_marker_excluded_native_gene_map.tsv',sep='\t').set_index('gene')
    contract=next(r for r in json.loads((old/'private_source_gene_measurement_contract.json').read_text()) if r['dataset']=='GSE286827')
    markercols=np.array(contract['physical_columns'],int);assert contract['requested_genes']==markers and (markercols>=0).all()
    active=gm.loc[weights.gene,'GSE286827'].to_numpy(int);assert (active>=0).all()
    physical=np.array(sorted(set(active)|set(markercols)),int);colmap={int(c):i for i,c in enumerate(physical)}
    models={};model_records=[];traininputs=frozen/'private_inputs'
    for vi,view in enumerate(VIEWS):
        md=pd.read_csv(traininputs/(view+'_discovery_metadata.tsv'),sep='\t');D=np.c_[np.ones(len(md)),md[['cellQ','geneQ','sourceCC','fraction','logN']].to_numpy(float)]
        genes=Path(str(traininputs/(view+'_discovery'))+'_genes.tsv').read_text().splitlines()
        with gzip.open(str(traininputs/(view+'_discovery'))+'_tpm.mtx.gz','rb') as h:tr=mmread(h).tocsr()
        B=np.zeros((len(physical),4));center=[];coef=[]
        for pi,program in enumerate(['MCP1','MCP2']):
            w=weights[weights.view.eq(view)&weights.program.eq(program)];assert w.gene.is_unique
            for si,variant in enumerate(['full','fixed_stress_complement']):
                wv=w if si==0 else w[~w.gene.str.upper().isin(STRESS)]
                beta=wv.weight.to_numpy()/wv.discovery_sd.to_numpy();cols=gm.loc[wv.gene,'GSE286827'].to_numpy(int);assert len(set(cols))==len(cols)
                k=2*pi+si;B[[colmap[int(c)] for c in cols],k]=beta;mu=float(beta@wv.discovery_mean.to_numpy());center.append(mu)
                tidx=pd.Index(genes).get_indexer(wv.gene);assert (tidx>=0).all()
                rawtrain=np.asarray(tr[tidx].T.dot(beta)).ravel()-mu
                co=np.linalg.lstsq(D,rawtrain,rcond=None)[0];coef.append(co)
                model_records.append(dict(view=view,program=program,variant=variant,active_genes=len(wv),center=mu,**{name:float(v) for name,v in zip(['intercept','cellQ','geneQ','sourceCC','fraction','logN'],co)}))
        f=np.load(traininputs/(view+'_frozen_projection.npz'))
        models[vi]=dict(B=B,center=np.array(center),coef=np.array(coef),mu=f['quality_mean'],sd=f['quality_sd'])
    atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1';z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r');counts=z['layers/counts'];nc,nf=map(int,counts.attrs['shape']);assert len(mapping)==nc
    candidates=cells.global_row.to_numpy(int);n=len(cells);RNA=np.zeros(n,np.int64);detected=np.zeros(n,np.int32);quality=np.empty((n,4));scores=np.empty((n,4));selected_parts=[]
    totalmaterial=materials.all_material_cells.to_numpy(int);capture_n=cov.pivot(index='material_index',columns='view',values='full_view_cells')
    for start in range(0,nc,8192):
        stop=min(nc,start+8192);lo,hi=np.searchsorted(candidates,[start,stop])
        if hi>lo:
            ptr=np.asarray(counts['indptr'][start:stop+1],np.int64)
            full=csr_matrix((np.asarray(counts['data'][int(ptr[0]):int(ptr[-1])]),np.asarray(counts['indices'][int(ptr[0]):int(ptr[-1])]),ptr-ptr[0]),shape=(stop-start,nf))[candidates[lo:hi]-start].copy()
            full.sum_duplicates();full.eliminate_zeros();full.sort_indices();assert (full.data>0).all() and np.array_equal(full.data,np.floor(full.data)) and np.isfinite(full.data).all()
            RNA[lo:hi]=np.asarray(full.sum(axis=1)).ravel().astype(np.int64);detected[lo:hi]=np.diff(full.indptr);assert (RNA[lo:hi]>0).all()
            rawunion=full[:,physical].tocsr();take=cells.iloc[lo:hi].selected.to_numpy(bool)
            if take.any():selected_parts.append(rawunion[take])
            logged=rawunion.astype(float);logged.data=np.log1p(logged.data*10000/np.repeat(RNA[lo:hi],np.diff(logged.indptr)))
            for vi,model in models.items():
                local=np.flatnonzero(cells.iloc[lo:hi].view_index.to_numpy()==vi)
                if not len(local):continue
                ids=lo+local;mi=cells.iloc[ids].material_index.to_numpy(int);cnt=capture_n.loc[mi,VIEWS[vi]].to_numpy(int)
                q=np.c_[np.log1p(RNA[ids]),np.log1p(detected[ids]),cnt/totalmaterial[mi],np.log1p(cnt)]
                q=(q-model['mu'])/model['sd'];quality[ids]=q
                scores[ids]=np.asarray(logged[local]@model['B'])-model['center']-np.c_[np.ones(len(ids)),q]@model['coef'][:,[0,1,2,4,5]].T
            tind=np.flatnonzero(cells.iloc[lo:hi].view_index.to_numpy()!=2)
            if len(tind):
                ix=np.searchsorted(cache['global_original_T_row_index'],candidates[lo:hi][tind]);assert np.array_equal(cache['global_original_T_row_index'][ix],candidates[lo:hi][tind])
                assert np.array_equal(rawunion[tind][:,[colmap[int(c)] for c in markercols]].toarray(),raw[ix])
                assert np.array_equal(RNA[lo+tind],cache['full_native_total_counts'][ix]) and np.array_equal(detected[lo+tind],cache['full_native_detected_genes'][ix])
        if start//8192%16==0 or stop==nc:(out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),native_rows_scanned=stop,native_rows=nc,captured_cells_read=int(hi),captured_cells=n))+'\n')
    assert np.isfinite(scores).all() and np.isfinite(quality).all()
    chosen=cells.selected.to_numpy(bool);sel=cells.loc[chosen].copy();sel['fullRNA']=RNA[chosen];sel['full_detected_genes']=detected[chosen]
    for j,f in enumerate(FIELDS):sel[f]=quality[chosen,j]
    for k,f in enumerate(['MCP1_full','MCP1_fixed_stress_complement','MCP2_full','MCP2_fixed_stress_complement']):sel[f]=scores[chosen,k]
    save_npz(out/'private_selected_raw_union_counts.npz',vstack(selected_parts,format='csr'),compressed=True);np.save(out/'private_union_native_columns.npy',physical,allow_pickle=False)
    sel.to_csv(out/'private_selected_cells_and_all_scores.tsv.gz',sep='\t',index=False);materials.to_csv(out/'private_all52_paired_materials.tsv',sep='\t',index=False);cov.to_csv(out/'private_all156_capture_scope.tsv',sep='\t',index=False)
    pd.DataFrame(model_records).to_csv(out/'private_all12_original_frozen_models.tsv',sep='\t',index=False)
    case_rows=[];domain=[]
    for vi,view in enumerate(VIEWS):
        train=pd.read_csv(traininputs/(view+'_discovery_metadata.tsv'),sep='\t')
        for j,f in enumerate(FIELDS):
            v=quality[cells.view_index.eq(vi),j];a=float(train[f].min());b=float(train[f].max())
            domain.append(dict(view=view,covariate=f,training_min=a,training_max=b,captured_cells=len(v),new_min=float(v.min()),new_max=float(v.max()),outside_training_cells=int(((v<a)|(v>b)).sum())))
        for m in range(52):
            ix=np.flatnonzero(cells.view_index.eq(vi)&cells.material_index.eq(m));sx=ix[chosen[ix]];assert len(sx)==min(100,len(ix))
            for pi,program in enumerate(['MCP1','MCP2']):
                for si,variant in enumerate(['full','fixed_stress_complement']):
                    k=2*pi+si;case_rows.append(dict(sampleID=materials.iloc[m].sampleID,patient_registry=materials.iloc[m].patient_registry,phase=materials.iloc[m].phase,material_index=m,view=view,program=program,variant=variant,score=float(scores[sx,k].mean()),all_captured_mean_sensitivity=float(scores[ix,k].mean()),selected_cells=len(sx),full_view_cells=len(ix),all_material_cells=int(totalmaterial[m]),**{f:float(quality[sx,j].mean()) for j,f in enumerate(FIELDS)}))
    case=pd.DataFrame(case_rows);assert len(case)==624
    baseline=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_all_case_mean_scores.tsv',sep='\t')
    baseline=baseline[baseline.fit_id.eq('K2_seed20261018_source_adjusted')&baseline.partition.eq('discovery')]
    original=pd.read_csv(scope/'private_metadata_only_case_source_partition.tsv',sep='\t').set_index('fit_sample').sampleID
    baseline['sampleID']=baseline['sample'].map(original);assert baseline.sampleID.notna().all()
    matches=case[case.variant.eq('full')].merge(baseline[['sampleID','view','program','score']],on=['sampleID','view','program'],suffixes=('_new','_accepted'),validate='one_to_one')
    anchor=float(np.max(np.abs(matches.score_new-matches.score_accepted)));assert len(matches)>0 and anchor<1e-8
    scales=baseline.groupby(['view','program']).score.std(ddof=0).reset_index(name='discovery_case_sd');assert len(scales)==6 and scales.discovery_case_sd.gt(0).all()
    case.to_csv(out/'private_all624_material_programme_measurements.tsv',sep='\t',index=False);scales.to_csv(out/'private_all6_discovery_case_scales.tsv',sep='\t',index=False);pd.DataFrame(domain).to_csv(out/'private_all12_training_domain_descriptions.tsv',sep='\t',index=False)
    matches.to_csv(out/'private_all_original_discovery_anchor_checks.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()};assert sum(v['bytes'] for v in files.values())<cap
    terminal=dict(status='COMPLETE_ALL52_PAIRED_NATIVE_FULL_FROZEN_RNA_MEASUREMENTS_NO_EFFECTS_AUDIT_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,all_captured_cells=n,selected_cells=int(chosen.sum()),measurements=624,original_anchor_measurements=len(matches),maximum_original_anchor_error=anchor,no_paired_or_clinical_effects=True)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(dict(status=terminal['status'],terminal_sha256=sha(out/'private_terminal.json'),all_captured_cells=n,selected_cells=int(chosen.sum()))))


if __name__=='__main__':main()

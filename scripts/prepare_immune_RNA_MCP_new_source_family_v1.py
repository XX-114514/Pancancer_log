"""Freeze and assemble all new-context/source and author clinical questions."""
import hashlib
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from immune_RNA_MCP_new_source_math_v1 import design

VIEWS=['R_FOXP3_capture','E_effector_capture','Myeloid']
PAIRS=[(0,1),(0,2),(1,2)]
QFIELDS=['cellQ','geneQ','fraction','logN']
CONTEXT=['dataset','cancer_label','tissue','material_class','phase_class','capture_scope','ex_vivo']
CLOSED={'GSE286827','GSE294109','GSE305078'}
TARGET='CXG_5b9b89d2_2e3e_4ad1_8965_666f64540fea'


def verify(directory,expected):
    assert sha(directory/'private_terminal.json')==expected
    t=json.loads((directory/'private_terminal.json').read_text())
    for name,v in t.get('files',{}).items():assert sha(directory/name)==v['sha256']
    return t


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,16*1024**2)
    native=root/'immune_RNA_MCP_all_atlas_native_measurement_20261019_v2';verify(native,'2463d66299d7b1f1be96d99a19554a834735e688440716edfa00cd1aa01575e9')
    na=root/'immune_RNA_MCP_all_atlas_native_measurement_audit_20261019_v1/private_terminal.json';assert sha(na)=='f7a6db048cfccac2a98de84023ce5048d4faa595353d2589767718d0e661cfc0'
    clinical=root/'CAR_product_RNA_MCP_clinical_metadata_qualification_20261019_v2';verify(clinical,'00a13975d441239e4844b17359124e96b73d127e2f04d747fbdfc6bf1929cdd9')
    ca=root/'CAR_product_RNA_MCP_clinical_metadata_audit_20261019_v3/private_terminal.json';assert sha(ca)=='dde9f19ebf0f933ac9ab883c883d45e32eee4db008075423403f20d274b9adc3'
    out=root/'immune_RNA_MCP_new_source_family_inputs_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),math_sha256=sha(Path(__file__).with_name('immune_RNA_MCP_new_source_math_v1.py')),resource_guard=guard,native_acceptance_sha256=sha(na),clinical_acceptance_sha256=sha(ca),
        source_scope='All atlas sources excluding closed HNSCC,ccRCC,GSE305078; all existing source/context strata and2programmes/3pairs. Same-material pairs then equal-material reviewed-case mean within context; unknown original sample units retained.',context_fields=CONTEXT,closed_sources=sorted(CLOSED),
        source_families='Marginal Spearman n>=10; adjusted partial-rank8quality/proportion controls n>=20 and10residualdf. Separate BH/BY over ALL metadata-preeligible slots. Constant/failed-design slots remain P1 in those universes. All ineligible slots retained.',
        clinical_scope='All3pairs/2programmes/2fixed contrasts/4predefined endpoints retained. Metadata-preeligible n>=20, binary each>=10. Primary response3m family separate from complete secondary toxicity family; no inference-dependent endpoint choice.',
        clinical_contrasts='0.5*(A/SD_A+ B/SD_B) shared_mean and0.5*(A/SD_A-B/SD_B) A_minus_B. Original38discovery-case SDs only; no clinical weight fit. A clinical contrast is not both endpoints independently clinical or direct communication.',
        clinical_controls='Both8original quality/proportion means,sex_female,parent coreT authorCD4/CD8 fractions. Whole stable SVD nuisance span, no automatic covariate deletion. Extra4within R/E authorCD4/CD8 fractions sensitivity only.',
        endpoints={'response3m':'R1 NR0, original CRvsnonCR, unknown remains missing; primary','ICANS_group':'Grade3-4=1,Grade0-2=0; secondary','CRS_max_grade':'numeric0..4 author ordinal severity, partial-rank effect; no threshold; secondary','prolonged_cytopenia':'Yes1 No0 only; all other tokens unreviewed; time cutoff unverified, secondary author-label association only'},
        inference='19999case permutations marginal/Wilcoxon and HC3 restricted-null wild adjusted;2000casebootstrap, binary endpoint-group stratified, others wholecase;CI>=1900valid,LOO,MCupper. Main corrected clinical effects use fullbackground. Allpoint sensitivities and raw Wilcoxon/Spearman separate fullfamily correction.',
        sensitivities=['all_captured_own_quality','fixed15_stress_complement','RNA_quality_only','abundance_only','plus_view_lineage_mix','component_A','component_B'],
        limits='Single CAR manufacturedproduct cohort, priorprojectsource use; no untouched-project,wholeTME,prediction,causalinteraction or formalTreg/CD8 claim. Marginalsource coactivity may reflect quality/abundance/drugmix; sources/contexts sharing cases not independent validation. No new expression downloads.',
        effect_values_or_P_computed=False)
    (out/'plan_before_any_new_source_or_clinical_effect.json').write_text(json.dumps(plan,indent=2)+'\n')
    lib=pd.read_csv(native/'private_all400_materials.tsv',sep='\t',keep_default_na=False).set_index('sampleID',drop=False)
    opp=root/'immune_RNA_MCP_all_atlas_opportunities_20261019_v1';pairs=pd.read_csv(opp/'private_all2400_material_programme_pair_qualification.tsv.gz',sep='\t')
    meta=pd.read_csv(clinical/'private_all59_exact_author_endpoint_and_patient_bridge.tsv',sep='\t',keep_default_na=False).set_index('library_id',drop=False)
    collection=root.parents[1]/'runs/20261010_clinical_scrna_expansion_v1';conversion=json.loads((collection/'manifests'/f'{TARGET}_native_conversion.json').read_text());accepted=json.loads((collection/'manifests'/f'accepted_{TARGET}.json').read_text())
    fields=['sample_id','donor_id','response3m','3mo PET/CT','cell_type','sex','ICANS group','CRS max grade','prolonged cytopenia']
    def original_obs(path):
        with h5py.File(path,'r') as h:
            o=h['obs'];index=np.asarray(read_elem(o[str(o.attrs['_index'])])).astype(str);d=pd.DataFrame({f:np.asarray(read_elem(o[f])).astype(str) for f in fields},index=pd.Index(index,name='original_cell_id'))
        return d
    original=original_obs(conversion['source_path']);qc=zarr.open_group(accepted['output'],mode='r');converted=original_obs(str(read_elem(qc['uns/raw_counts_source'])));pd.testing.assert_frame_equal(original,converted)
    assert len(original)==417167 and original.index.is_unique
    for f in fields:assert original.groupby('donor_id')[f].nunique(dropna=False).eq(1).all()
    bydonor=original.groupby('donor_id',sort=True).first();endpoints=[]
    for key,r in meta.iterrows():
        d=bydonor.loc[r.original_donor_id];assert d.response3m==r.author_response3m and d['3mo PET/CT']==r.author_PET_CT_3month and d.sample_id==r.original_sample_id
        response={'R':1.,'NR':0.}.get(d.response3m,np.nan);icans={'Grade 3-4':1.,'Grade 0-2':0.}.get(d['ICANS group'],np.nan);cyto={'Yes':1.,'No':0.}.get(d['prolonged cytopenia'],np.nan)
        try:crs=float(d['CRS max grade'])
        except ValueError:crs=np.nan
        assert not np.isfinite(crs) or crs in [0,1,2,3,4]
        assert d.sex in ['male','female']
        endpoints.append(dict(sampleID=TARGET+'::'+key,library_id=key,case_block=TARGET+'::'+r.patient_block,original_donor_id=r.original_donor_id,response3m=response,ICANS_group=icans,CRS_max_grade=crs,prolonged_cytopenia=cyto,sex_female=float(d.sex=='female'),**{'raw_'+f:d[f] for f in ['sex','ICANS group','CRS max grade','prolonged cytopenia']}))
    endpoint=pd.DataFrame(endpoints).set_index('sampleID',drop=False);assert len(endpoint)==59 and endpoint.case_block.is_unique
    atlas=zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),mode='r');labels=np.asarray(read_elem(atlas['obs/source_author_label'])).astype(str);lower=np.char.lower(labels);isCD4=np.char.find(lower,'cd4')>=0;isCD8=np.char.find(lower,'cd8')>=0;atlas_dataset=np.asarray(read_elem(atlas['obs/dataset'])).astype(str);assert not (isCD4[atlas_dataset==TARGET]&isCD8[atlas_dataset==TARGET]).any()
    old=root/'immune_CM_fixed_panel_same_cell_RNA_20261018_v1';cache=np.load(old/'private_fixed34_gene_T_counts_and_depth.npz');T=pd.read_csv(old/'private_original_T_unit_scope.tsv',sep='\t');cm=cache['genes'].tolist();fixed=root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json';panel=json.loads(fixed.read_text())['panels'];core=np.count_nonzero(cache['raw_gene_counts'][:,[cm.index(g) for g in panel['pan_T']]]>0,axis=1)>=2
    unitkeys=(T.dataset.astype(str)+'::'+T.library_id.astype(str)).to_numpy(str);tkeys=unitkeys[cache['original_T_unit_index']];tr=cache['global_original_T_row_index'];pools=pd.read_csv(native/'private_all_capture_pools.tsv',sep='\t');membership=np.load(native/'private_all_captured_membership.npz');poolids=membership['pool_index'];gr=membership['global_row'];profile=[]
    for key in endpoint.index:
        rr=tr[core&(tkeys==key)];assert len(rr)>0
        row=dict(sampleID=key,parent_core_T=len(rr),parent_CD4=float(isCD4[rr].mean()),parent_CD8=float(isCD8[rr].mean()))
        for vi,view in enumerate(VIEWS[:2]):
            p=pools[pools.sampleID.eq(key)&pools['view'].eq(view)]
            if len(p):
                ids=gr[poolids==int(p.pool_index.iloc[0])];row['view'+str(vi)+'_CD4']=float(isCD4[ids].mean());row['view'+str(vi)+'_CD8']=float(isCD8[ids].mean())
            else:row['view'+str(vi)+'_CD4']=np.nan;row['view'+str(vi)+'_CD8']=np.nan
        profile.append(row)
    profile=pd.DataFrame(profile).set_index('sampleID');endpoint=endpoint.join(profile.drop(columns=[],errors='ignore'),validate='one_to_one')
    endpoint.to_csv(out/'private_all59_original_endpoints_and_lineage_controls.tsv',sep='\t',index=False)
    measures=pd.read_csv(native/'private_all3488_programme_measurements.tsv',sep='\t');index=measures.set_index(['sampleID','view','program','variant']);assert index.index.is_unique
    official=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_all_case_mean_scores.tsv',sep='\t');official=official[official.fit_id.eq('K2_seed20261018_source_adjusted')&official.partition.eq('discovery')];scales=official.groupby(['view','program']).score.std(ddof=0);assert len(scales)==6 and scales.gt(0).all()
    scales.rename('discovery_case_sd').to_csv(out/'private_all6_original_discovery_case_scales.tsv',sep='\t')
    registry=[];model_no=0
    def save_model(record,arrays):
        nonlocal model_no
        model_no+=1;record['model_id']='N%04d'%model_no;record['array_file']='private_'+record['model_id']+'_inputs.npz';np.savez_compressed(out/record['array_file'],**arrays);registry.append(record)
    quality_names=np.array(['quality::'+end+'::'+q for end in ['A','B'] for q in QFIELDS])
    def material_arrays(keys,va,vb,program):
        values=[]
        for key in keys:
            a=index.loc[(key,va,program,'full')];b=index.loc[(key,vb,program,'full')];sa=index.loc[(key,va,program,'fixed_stress_complement')];sb=index.loc[(key,vb,program,'fixed_stress_complement')]
            values.append([a.score,b.score,a.all_captured_mean_sensitivity,b.all_captured_mean_sensitivity,sa.score,sb.score,*a[QFIELDS].to_numpy(float),*b[QFIELDS].to_numpy(float),*a[['all_'+q for q in QFIELDS]].to_numpy(float),*b[['all_'+q for q in QFIELDS]].to_numpy(float)])
        return np.array(values,float).reshape(len(keys),22)
    source_contexts=[]
    for context,g in lib[~lib.dataset.isin(CLOSED)].groupby(CONTEXT,sort=True,dropna=False):
        descriptor=dict(zip(CONTEXT,map(str,context)));context_id=hashlib.sha256(json.dumps(descriptor,sort_keys=True).encode()).hexdigest()[:16];source_contexts.append(dict(context_id=context_id,**descriptor,materials=len(g)))
        for ai,bi in PAIRS:
            va,vb=VIEWS[ai],VIEWS[bi];pair=va+'|'+vb
            for program in ['MCP1','MCP2']:
                qualified=pairs[pairs.sampleID.isin(g.index)&pairs.pair.eq(pair)&pairs.program.eq(program)&pairs.eligible];keys=sorted(qualified.sampleID);m=material_arrays(keys,va,vb,program)
                if keys:
                    case_labels=lib.loc[keys,'case_block'].to_numpy(str);group=pd.DataFrame(m);group['case']=case_labels;averaged=group.groupby('case',sort=True).mean();cases=averaged.index.to_numpy(str);a=averaged.to_numpy(float)
                else:cases=np.array([],str);a=np.empty((0,22))
                for mode,minimum in [('source_marginal',10),('source_adjusted',20)]:
                    C=np.empty((len(a),0)) if mode=='source_marginal' else a[:,6:14];names=np.array([],str) if mode=='source_marginal' else quality_names;n=len(a);preeligible=n>=minimum
                    status='METADATA_PREELIGIBLE' if preeligible else 'FEWER_THAN_REQUIRED_CASE_UNITS'
                    if mode=='source_adjusted' and preeligible:
                        _,status=design(a[:,0],a[:,1],C,names,'coordination')
                    elif mode=='source_marginal' and preeligible:status='ESTIMABLE' if min(np.ptp(a[:,0]),np.ptp(a[:,1]))>0 else 'CONSTANT_TARGET_OR_READOUT'
                    record=dict(family=mode,dataset=descriptor['dataset'],context_id=context_id,pair=pair,program=program,contrast='coactivity',endpoint='none',n_units=n,materials=len(keys),reviewed_units=int(lib.loc[keys].groupby('case_block').reviewed_patient.all().sum()) if keys else 0,metadata_preeligible=preeligible,status=status,kind='coordination',binary_endpoint=False)
                    save_model(record,dict(cases=cases,material_keys=np.array(keys),material_case_labels=lib.loc[keys,'case_block'].to_numpy(str) if keys else np.array([],str),x=a[:,0],y=a[:,1],C=C,names=names,allcapture_x=a[:,2],allcapture_y=a[:,3],allcapture_C=np.empty((len(a),0)) if mode=='source_marginal' else a[:,14:22],stress_x=a[:,4],stress_y=a[:,5],strata=np.zeros(len(a),int)))
    pd.DataFrame(source_contexts).to_csv(out/'private_all_new_source_contexts.tsv',sep='\t',index=False)
    for ai,bi in PAIRS:
        va,vb=VIEWS[ai],VIEWS[bi];pair=va+'|'+vb
        for program in ['MCP1','MCP2']:
            qualified=pairs[pairs.dataset.eq(TARGET)&pairs.pair.eq(pair)&pairs.program.eq(program)&pairs.eligible]
            for ep in ['response3m','ICANS_group','CRS_max_grade','prolonged_cytopenia']:
                keys=sorted(set(qualified.sampleID)&set(endpoint.index[np.isfinite(endpoint[ep])]))
                a=material_arrays(keys,va,vb,program);md=endpoint.loc[keys];n=len(keys);binary=ep!='CRS_max_grade';x=md[ep].to_numpy(float);n0=int((x==0).sum());n1=int((x==1).sum());preeligible=bool(n>=20 and (not binary or min(n0,n1)>=10));background=md[['sex_female','parent_CD4','parent_CD8']].to_numpy(float);C=np.c_[a[:,6:14],background];allC=np.c_[a[:,14:22],background];names=np.r_[quality_names,['sex_female','parent_CD4','parent_CD8']];fineC=np.c_[C,md[['view0_CD4','view0_CD8','view1_CD4','view1_CD8']].to_numpy(float)];fine_names=np.r_[names,['view0_CD4','view0_CD8','view1_CD4','view1_CD8']]
                scale=np.array([scales.loc[(va,program)],scales.loc[(vb,program)]])
                for contrast,sign in [('shared_mean',1),('A_minus_B',-1)]:
                    y=.5*(a[:,0]/scale[0]+sign*a[:,1]/scale[1]);ay=.5*(a[:,2]/scale[0]+sign*a[:,3]/scale[1]);sy=.5*(a[:,4]/scale[0]+sign*a[:,5]/scale[1]);kind='clinical' if binary else 'coordination';status='FEWER_THAN_REQUIRED_CASE_UNITS_OR_BINARY_GROUPS'
                    if preeligible:_,status=design(x,y,C,names,kind)
                    record=dict(family='clinical_primary' if ep=='response3m' else 'clinical_secondary',dataset=TARGET,context_id='manufactured_infusion_product',pair=pair,program=program,contrast=contrast,endpoint=ep,n_units=n,materials=n,reviewed_units=n,metadata_preeligible=preeligible,status=status,kind=kind,binary_endpoint=binary,n0=n0,n1=n1)
                    save_model(record,dict(cases=md.case_block.to_numpy(str),material_keys=np.array(keys),x=x,y=y,C=C,names=names,allcapture_x=x,allcapture_y=ay,allcapture_C=allC,stress_y=sy,component_A=a[:,0]/scale[0],component_B=a[:,1]/scale[1],fineC=fineC,fine_names=fine_names,strata=x.astype(int) if binary else np.zeros(n,int)))
    reg=pd.DataFrame(registry);assert len(reg[reg.family.str.startswith('clinical')])==48
    reg.to_csv(out/'private_all_registered_models_and_pre_effect_eligibility.tsv',sep='\t',index=False)
    summary=reg.groupby('family').agg(registry_slots=('model_id','size'),metadata_preeligible=('metadata_preeligible','sum'),estimable=('status',lambda x:int(x.eq('ESTIMABLE').sum()))).reset_index();summary.to_csv(out/'private_family_technical_dimensions.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<16*1024**2
    terminal=dict(status='COMPLETE_ALL_NEW_SOURCE_CONTEXT_AND_CAR_PRIMARY_SECONDARY_INPUTS_NO_EFFECTS_AUDIT_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,registered_models=len(reg),effect_values_or_P_computed=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(dict(status=terminal['status'],terminal_sha256=sha(out/'private_terminal.json'),technical=summary.to_dict(orient='records'))))


if __name__=='__main__':main()

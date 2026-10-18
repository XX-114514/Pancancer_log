"""Every effect-free input reconstructed from original metadata and assay pivots."""
import json
import os
import hashlib
from datetime import datetime,timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from immune_RNA_MCP_new_source_math_v1 import independent_design

VIEWS=['R_FOXP3_capture','E_effector_capture','Myeloid']
Q=['cellQ','geneQ','fraction','logN']
CONTEXT=['dataset','cancer_label','tissue','material_class','phase_class','capture_scope','ex_vivo']
TARGET='CXG_5b9b89d2_2e3e_4ad1_8965_666f64540fea'


def close(a,b):assert np.allclose(np.asarray(a,float),np.asarray(b,float),rtol=0,atol=1e-8,equal_nan=True)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2);inp=root/'immune_RNA_MCP_new_source_family_inputs_20261019_v2';assert sha(inp/'private_terminal.json')=='9796f0ba9fc23a189814765cf7b9d6774600f17f2fd3fecbf83beac2b45faf11'
    term=json.loads((inp/'private_terminal.json').read_text());assert not term['effect_values_or_P_computed']
    for n,v in term['files'].items():assert sha(inp/n)==v['sha256']
    plan=json.loads((inp/'plan_before_any_new_source_or_clinical_effect.json').read_text());assert plan['closed_sources']==['GSE286827','GSE294109','GSE305078']
    native=root/'immune_RNA_MCP_all_atlas_native_measurement_20261019_v2';lib=pd.read_csv(native/'private_all400_materials.tsv',sep='\t',keep_default_na=False).set_index('sampleID',drop=False);pools=pd.read_csv(native/'private_all_capture_pools.tsv',sep='\t');member=np.load(native/'private_all_captured_membership.npz')
    collection=root.parents[1]/'runs/20261010_clinical_scrna_expansion_v1';cv=json.loads((collection/'manifests'/f'{TARGET}_native_conversion.json').read_text());accepted=json.loads((collection/'manifests'/f'accepted_{TARGET}.json').read_text());qc=zarr.open_group(accepted['output'],mode='r');fields=['sample_id','donor_id','response3m','3mo PET/CT','cell_type','sex','ICANS group','CRS max grade','prolonged cytopenia']
    def readobs(path):
        with h5py.File(path,'r') as h:
            g=h['obs'];return pd.DataFrame({f:np.asarray(read_elem(g[f])).astype(str) for f in fields},index=pd.Index(np.asarray(read_elem(g[str(g.attrs['_index'])])).astype(str),name='original_cell_id'))
    original=readobs(cv['source_path']);pd.testing.assert_frame_equal(original,readobs(str(read_elem(qc['uns/raw_counts_source']))));assert len(original)==417167
    for f in set(fields)-{'cell_type'}:assert original.groupby('donor_id')[f].nunique(dropna=False).eq(1).all()
    donor=original.groupby('donor_id').first();prior=pd.read_csv(root/'CAR_product_RNA_MCP_clinical_metadata_qualification_20261019_v2/private_all59_exact_author_endpoint_and_patient_bridge.tsv',sep='\t',keep_default_na=False).set_index('library_id');md=pd.read_csv(inp/'private_all59_original_endpoints_and_lineage_controls.tsv',sep='\t',keep_default_na=False).set_index('sampleID',drop=False)
    assert len(md)==59 and md.case_block.is_unique
    atlas=zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),mode='r');get=lambda f:np.asarray(read_elem(atlas['obs'][f])).astype(str);source=get('dataset');carrows=np.flatnonzero(source==TARGET);ids=get('original_cell_id')[carrows];pos=original.index.get_indexer(ids);assert (pos>=0).all()
    labels=original.cell_type.to_numpy(str)[pos];lower=np.char.lower(labels);cd4=np.zeros(len(source),bool);cd8=cd4.copy();cd4[carrows]=np.char.find(lower,'cd4')>=0;cd8[carrows]=np.char.find(lower,'cd8')>=0;assert not (cd4&cd8).any()
    old=root/'immune_CM_fixed_panel_same_cell_RNA_20261018_v1';cache=np.load(old/'private_fixed34_gene_T_counts_and_depth.npz');units=pd.read_csv(old/'private_original_T_unit_scope.tsv',sep='\t');panel=json.loads((root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json').read_text())['panels'];counts=cache['raw_gene_counts'];gn=cache['genes'].tolist();core=(counts[:,[gn.index(g) for g in panel['pan_T']]]>=1).sum(axis=1)>=2;keys=(units.dataset.astype(str)+'::'+units.library_id.astype(str)).to_numpy(str)[cache['original_T_unit_index']];rr=cache['global_original_T_row_index'][core];kk=keys[core]
    parent=pd.DataFrame(dict(sampleID=kk,CD4=cd4[rr].astype(int),CD8=cd8[rr].astype(int))).groupby('sampleID').agg(parent_core_T=('CD4','size'),parent_CD4=('CD4','mean'),parent_CD8=('CD8','mean'))
    for key,row in md.iterrows():
        p=prior.loc[row.library_id];d=donor.loc[p.original_donor_id];assert row.case_block==TARGET+'::'+p.patient_block and row.original_donor_id==p.original_donor_id
        expected={'response3m':{'R':1.,'NR':0.}.get(d.response3m,np.nan),'ICANS_group':{'Grade 0-2':0.,'Grade 3-4':1.}.get(d['ICANS group'],np.nan),'prolonged_cytopenia':{'No':0.,'Yes':1.}.get(d['prolonged cytopenia'],np.nan),'sex_female':float(d.sex=='female')}
        try:expected['CRS_max_grade']=float(d['CRS max grade'])
        except ValueError:expected['CRS_max_grade']=np.nan
        for f,v in expected.items():close(float(row[f]) if row[f]!='' else np.nan,v)
        for f in ['sex','ICANS group','CRS max grade','prolonged cytopenia']:assert row['raw_'+f]==d[f]
        close(row[['parent_core_T','parent_CD4','parent_CD8']],parent.loc[key]);assert row.parent_core_T==lib.loc[key,'core_T']
        for vi,view in enumerate(VIEWS[:2]):
            pool=pools.loc[pools.sampleID.eq(key)&pools['view'].eq(view),'pool_index'];rows=member['global_row'][member['pool_index']==int(pool.iloc[0])] if len(pool) else np.array([],int)
            for name,flags in [('CD4',cd4),('CD8',cd8)]:
                value=float(row['view'+str(vi)+'_'+name]) if row['view'+str(vi)+'_'+name]!='' else np.nan;close(value,float(flags[rows].mean()) if len(rows) else np.nan)
    # After metadata has been checked, numeric columns can safely use NaN.
    md=pd.read_csv(inp/'private_all59_original_endpoints_and_lineage_controls.tsv',sep='\t').set_index('sampleID',drop=False)
    measures=pd.read_csv(native/'private_all3488_programme_measurements.tsv',sep='\t');wide={p:g.pivot(index='sampleID',columns=['view','variant'],values=['score','all_captured_mean_sensitivity',*Q,*['all_'+f for f in Q]]) for p,g in measures.groupby('program')}
    official=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_all_case_mean_scores.tsv',sep='\t');official=official[official.fit_id.eq('K2_seed20261018_source_adjusted')&official.partition.eq('discovery')];scales={(v,p):float(np.std(g.score,ddof=0)) for (v,p),g in official.groupby(['view','program'])};savedsc=pd.read_csv(inp/'private_all6_original_discovery_case_scales.tsv',sep='\t').set_index(['view','program'])
    for key,sd in scales.items():close(savedsc.loc[key,'discovery_case_sd'],sd)
    def vector(keys,va,vb,program):
        w=wide[program];cols=[('score',va,'full'),('score',vb,'full'),('all_captured_mean_sensitivity',va,'full'),('all_captured_mean_sensitivity',vb,'full'),('score',va,'fixed_stress_complement'),('score',vb,'fixed_stress_complement')]+[(f,v,'full') for fset in [Q,['all_'+q for q in Q]] for v in [va,vb] for f in fset]
        return w.reindex(keys).loc[:,cols].to_numpy(float).reshape(len(keys),22)
    pairs=pd.read_csv(root/'immune_RNA_MCP_all_atlas_opportunities_20261019_v1/private_all2400_material_programme_pair_qualification.tsv.gz',sep='\t');contexts=pd.read_csv(inp/'private_all_new_source_contexts.tsv',sep='\t',keep_default_na=False).set_index('context_id');registry=pd.read_csv(inp/'private_all_registered_models_and_pre_effect_eligibility.tsv',sep='\t',keep_default_na=False);assert len(registry)==420 and registry.model_id.is_unique
    expected_contexts={}
    for key,g in lib[~lib.dataset.isin(plan['closed_sources'])].groupby(CONTEXT,sort=True,dropna=False):
        desc=dict(zip(CONTEXT,map(str,key)));cid=hashlib.sha256(json.dumps(desc,sort_keys=True).encode()).hexdigest()[:16];expected_contexts[cid]=(desc,g.index.tolist())
    assert set(contexts.index)==set(expected_contexts)
    qnames=np.array(['quality::'+v+'::'+f for v in ['A','B'] for f in Q]);checked=0;expected_set=set()
    for cid,(desc,contextkeys) in expected_contexts.items():
        assert all(str(contexts.loc[cid,f])==desc[f] for f in CONTEXT) and contexts.loc[cid,'materials']==len(contextkeys)
        for va,vb in [(VIEWS[0],VIEWS[1]),(VIEWS[0],VIEWS[2]),(VIEWS[1],VIEWS[2])]:
            for p in ['MCP1','MCP2']:
                for family in ['source_marginal','source_adjusted']:expected_set.add((family,cid,va+'|'+vb,p))
    actual_set={tuple(r) for r in registry[registry.family.str.startswith('source')][['family','context_id','pair','program']].to_numpy()};assert actual_set==expected_set
    clinical_set={(va+'|'+vb,p,ep,c) for va,vb in [(VIEWS[0],VIEWS[1]),(VIEWS[0],VIEWS[2]),(VIEWS[1],VIEWS[2])] for p in ['MCP1','MCP2'] for ep in ['response3m','ICANS_group','CRS_max_grade','prolonged_cytopenia'] for c in ['shared_mean','A_minus_B']};assert {tuple(r) for r in registry[registry.family.str.startswith('clinical')][['pair','program','endpoint','contrast']].to_numpy()}==clinical_set
    for r in registry.itertuples(index=False):
        a=np.load(inp/r.array_file);va,vb=r.pair.split('|');qual=pairs[pairs.pair.eq(r.pair)&pairs.program.eq(r.program)&pairs.eligible]
        if r.family.startswith('source'):
            keys=sorted(set(qual.sampleID)&set(expected_contexts[r.context_id][1]));m=vector(keys,va,vb,r.program);frame=pd.DataFrame(m);frame['case']=lib.loc[keys,'case_block'].to_numpy(str);agg=frame.groupby('case',sort=True).mean();cases=agg.index.to_numpy(str);v=agg.to_numpy(float).reshape(len(cases),22);x=v[:,0];y=v[:,1];C=np.empty((len(cases),0)) if r.family=='source_marginal' else v[:,6:14];names=np.array([],str) if r.family=='source_marginal' else qnames
            assert a['material_case_labels'].tolist()==lib.loc[keys,'case_block'].tolist();close(a['stress_x'],v[:,4]);close(a['stress_y'],v[:,5]);close(a['allcapture_x'],v[:,2]);close(a['allcapture_y'],v[:,3]);close(a['allcapture_C'],np.empty((len(cases),0)) if r.family=='source_marginal' else v[:,14:22]);preeligible=len(cases)>=(10 if r.family=='source_marginal' else 20)
            status='FEWER_THAN_REQUIRED_CASE_UNITS' if not preeligible else 'ESTIMABLE'
            if preeligible and r.family=='source_marginal' and min(np.ptp(x),np.ptp(y))==0:status='CONSTANT_TARGET_OR_READOUT'
        else:
            keys=sorted(set(qual.loc[qual.dataset.eq(TARGET),'sampleID'])&set(md.index[np.isfinite(md[r.endpoint])]))
            v=vector(keys,va,vb,r.program);mm=md.loc[keys];cases=mm.case_block.to_numpy(str);x=mm[r.endpoint].to_numpy(float);sign=1 if r.contrast=='shared_mean' else -1;s=np.array([scales[(va,r.program)],scales[(vb,r.program)]]);y=(v[:,0]/s[0]+sign*v[:,1]/s[1])/2;bg=mm[['sex_female','parent_CD4','parent_CD8']].to_numpy(float);C=np.c_[v[:,6:14],bg];names=np.r_[qnames,['sex_female','parent_CD4','parent_CD8']];preeligible=bool(len(cases)>=20 and (not r.binary_endpoint or min((x==0).sum(),(x==1).sum())>=10));status='FEWER_THAN_REQUIRED_CASE_UNITS_OR_BINARY_GROUPS'
            close(a['allcapture_x'],x);close(a['allcapture_y'],(v[:,2]/s[0]+sign*v[:,3]/s[1])/2);close(a['allcapture_C'],np.c_[v[:,14:22],bg]);close(a['stress_y'],(v[:,4]/s[0]+sign*v[:,5]/s[1])/2);close(a['component_A'],v[:,0]/s[0]);close(a['component_B'],v[:,1]/s[1]);close(a['fineC'],np.c_[C,mm[['view0_CD4','view0_CD8','view1_CD4','view1_CD8']].to_numpy(float)]);assert a['fine_names'].tolist()==np.r_[names,['view0_CD4','view0_CD8','view1_CD4','view1_CD8']].tolist();assert a['strata'].tolist()==(x.astype(int) if r.binary_endpoint else np.zeros(len(x),int)).tolist()
        assert a['cases'].tolist()==cases.tolist() and a['material_keys'].tolist()==keys and a['names'].tolist()==names.tolist();close(a['x'],x);close(a['y'],y);close(a['C'],C);assert r.n_units==len(cases) and r.materials==len(keys) and r.metadata_preeligible==preeligible
        if preeligible and r.family!='source_marginal':_,status=independent_design(x,y,C,names,r.kind)
        assert r.status==status;checked+=1
    out=root/'immune_RNA_MCP_new_source_family_input_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700);receipt=dict(status='ACCEPTED_ALL420_FROZEN_NEW_SOURCE_AND_CLINICAL_INPUTS_ORIGINAL_METADATA_AND_FULL_NUISANCE_DESIGNS_NO_EFFECTS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(inp/'private_terminal.json'),resource_guard=guard,all_models_checked=checked,all_original_metadata_identities=417167,all_original_metadata_fields=9,all_products_lineage_and_endpoint_records=59,new_association_or_clinical_effects=False)
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(dict(status=receipt['status'],terminal_sha256=sha(out/'private_terminal.json'),all_models_checked=checked)))


if __name__=='__main__':main()

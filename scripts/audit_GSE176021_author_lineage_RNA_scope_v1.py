"""Independently reconstruct all reference-lineage membership and paired scope."""
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
from immune_RNA_MCP_new_source_io_v1 import wait


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2);src=root/'GSE176021_author_lineage_RNA_scope_20261019_v2';out=root/'GSE176021_author_lineage_RNA_scope_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    t=wait(root,out,'prepare_GSE176021_author_lineage_RNA_scope_v2',src/'private_terminal.json',{'COMPLETE_AUTHOR_LINEAGE_ALL_MATERIAL_SCOPE_AND_PAIRED_FAMILY_FROZEN_NO_NEW_RNA_EFFECT'})
    for name,v in t['files'].items():assert sha(src/name)==v['sha256']
    side=root/'GSE176021_exact_original_author_annotation_sidecar_20261019_v2';assert sha(root/'GSE176021_exact_original_author_annotation_audit_20261019_v2/private_terminal.json')=='d8d21a1c27465d9b8248857c0028e5106509e43426cbdfef9f44ccb1d71c4169'
    native=root/'immune_RNA_MCP_all_atlas_native_measurement_20261019_v2';expected=pd.read_csv(native/'private_all400_materials.tsv',sep='\t',keep_default_na=False);expected=expected[expected.dataset.eq('GSE176021')].reset_index(drop=True);expected['new_material_index']=np.arange(len(expected));lib=pd.read_csv(src/'private_all_source_materials.tsv',sep='\t',keep_default_na=False);pd.testing.assert_frame_equal(lib,expected)
    a=zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),mode='r');ds=np.asarray(read_elem(a['obs/dataset'])).astype(str);gr=np.flatnonzero(ds=='GSE176021');library=np.asarray(read_elem(a['obs/library_id'])).astype(str)[gr];keys={s:i for i,s in enumerate(lib.sampleID)};mi=np.array([keys['GSE176021::'+x] for x in library]);assert np.array_equal(np.bincount(mi,minlength=len(lib)),lib.annotated_scope_cells)
    bridge=np.load(side/'private_all_atlas_cells_exact_author_bridge.npz');assert np.array_equal(gr,bridge['global_atlas_row']);names=json.loads((side/'private_author_label_lookup_and_coverage.json').read_text())['CD3']['labels'];c=bridge['CD3_label_code'];R=np.array([x>=0 and names[x]=='CD4-Treg' for x in c]);E=np.array([x>=0 and names[x][:4]=='CD8-' for x in c]);assert not (R&E).any()
    old=root/'immune_CM_fixed_panel_same_cell_RNA_20261018_v1';cache=np.load(old/'private_fixed34_gene_T_counts_and_depth.npz');T=cache['global_original_T_row_index'];native_counts=cache['raw_gene_counts'];marker={g:i for i,g in enumerate(cache['genes'])};p=root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json';panels=json.loads(p.read_text())['panels']
    npos=lambda genes:np.sum(native_counts[:,[marker[x] for x in genes]]>=1,axis=1);core=npos(panels['pan_T'])>=2;reg=(npos(panels['regulatory_associated'])>=2)&(npos(['FOXP3'])==1);eff=npos(['GZMB','PRF1','CTSW','NKG7'])>=2;capR=np.isin(gr,T[core&reg&~eff]);capE=np.isin(gr,T[core&eff&~reg]);oldp=pd.read_csv(native/'private_all_capture_pools.tsv',sep='\t');oldkeys=set(oldp[['sampleID','view']].itertuples(index=False,name=None))
    scope=pd.read_csv(src/'private_all_scope_material_view_positions.tsv',sep='\t');pools=pd.read_csv(src/'private_all_eligible_pools.tsv',sep='\t');saved=np.load(src/'private_all_scope_cell_memberships.npz');all_rows=[];all_p=[];all_s=[];p=0
    for si,sc in enumerate(['legacy_capture_anchor','author_confirmed_capture','all_author_lineage']):
        for m,material in lib.iterrows():
            for vi,(view,legacy,author) in enumerate(zip(['R_FOXP3_capture','E_effector_capture'],[capR,capE],[R,E])):
                mask=legacy if si==0 else legacy&author if si==1 else author;ids=gr[(mi==m)&mask];oldcount=int(material[view]);eligible=len(ids)>=20 and oldcount>=20 and (si!=0 or (material.sampleID,view) in oldkeys);row=scope[(scope.scope==sc)&scope.sampleID.eq(material.sampleID)&scope['view'].eq(view)];assert len(row)==1;r=row.iloc[0]
                assert r.cells==len(ids) and r.original_capture_cells==oldcount and bool(r.eligible)==eligible and bool(r.old_anchor_available)==((material.sampleID,view) in oldkeys) and r.new_material_index==m and r.annotated_scope_cells==material.annotated_scope_cells
                if si==0:assert len(ids)==oldcount
                if not eligible:assert r.pool_index==-1;continue
                z=pools.iloc[p];assert r.pool_index==z.pool_index==p and z.scope==sc and z.sampleID==material.sampleID and z['view']==view
                hashes=np.array([hashlib.sha256(('2026101802|'+view+'|'+str(x)).encode()).hexdigest() for x in ids]);chosen=np.zeros(len(ids),bool);chosen[np.argsort(hashes)[:100]]=True;assert z.selected_cells==chosen.sum();all_rows.extend(ids);all_p.extend([p]*len(ids));all_s.extend(chosen);p+=1
    assert p==len(pools);order=np.lexsort((all_p,all_rows));assert np.array_equal(saved['global_row'],np.array(all_rows)[order]) and np.array_equal(saved['pool_index'],np.array(all_p)[order]) and np.array_equal(saved['selected'],np.array(all_s)[order])
    pairs=pd.read_csv(src/'private_all_scope_material_pair_positions.tsv',sep='\t',keep_default_na=False);ps=pd.read_csv(src/'private_all_paired_case_positions.tsv',sep='\t',keep_default_na=False);coverage=pd.read_csv(src/'private_all_source_context_coverage.tsv',sep='\t')
    for sc in ['legacy_capture_anchor','author_confirmed_capture','all_author_lineage']:
        use=[]
        for r in lib.itertuples(index=False):
            part=scope[scope.scope.eq(sc)&scope.sampleID.eq(r.sampleID)];ok=bool(part.eligible.all());rr=pairs[pairs.scope.eq(sc)&pairs.sampleID.eq(r.sampleID)].iloc[0];assert rr.case_block==r.case_block and rr.tissue==r.tissue and str(rr.both_views_eligible).lower()==str(ok).lower()
            if ok:use.append(r.sampleID)
        x=lib[lib.sampleID.isin(use)];sets=x.groupby('case_block').tissue.apply(set);paired=sorted(k for k,v in sets.items() if 'tumor' in v and 'adjacent normal lung' in v);assert ps[ps.scope.eq(sc)].case_block.tolist()==paired
        for r in coverage[coverage.scope.eq(sc)].itertuples(index=False):
            if r.tissue=='PAIRED_TUMOUR_ADJACENT_NORMAL':assert r.case_units==len(paired) and r.complete_materials==len(x[x.case_block.isin(paired)&x.tissue.isin(['tumor','adjacent normal lung'])])
            else:assert r.case_units==x[x.tissue.eq(r.tissue)].case_block.nunique() and r.complete_materials==int(x.tissue.eq(r.tissue).sum())
    # Check actual original source-case identity, not just a curated flag.
    rawpath=root.parents[2]/'Project_v2/Downstream/0716_Data/raw_geo_supplements/GSE176021/GSE176021_rawdata.h5ad'
    with h5py.File(rawpath,'r') as h:rp=np.asarray(read_elem(h['obs/patient_id'])).astype(str)[bridge['original_raw_row']]
    links=pd.DataFrame({'raw_patient':rp,'sampleID':['GSE176021::'+x for x in library]}).merge(lib[['sampleID','case_block','reviewed_patient']],on='sampleID',validate='many_to_one');assert links.groupby('sampleID').raw_patient.nunique().eq(1).all()
    reviewed=links[links.reviewed_patient.astype(str).str.lower().isin(['true','1'])];assert reviewed.groupby('case_block').raw_patient.nunique().eq(1).all() and reviewed.groupby('raw_patient').case_block.nunique().eq(1).all()
    gm=pd.read_csv(src/'private_all904_native_gene_roles_and_mapping.tsv',sep='\t');orig=pd.read_csv(native/'private_all_source_auxiliary_native_gene_mapping.tsv.gz',sep='\t');orig=orig[orig.dataset.eq('GSE176021')].set_index('gene').loc[gm.gene].reset_index();pd.testing.assert_frame_equal(gm,orig)
    registered=pd.read_csv(src/'private_all6_paired_tissue_registered_tests.tsv',sep='\t');n=int(ps.scope.eq('all_author_lineage').sum());assert len(registered)==6 and registered.paired_case_units.eq(n).all() and registered.metadata_preeligible.eq(n>=10).all();assert t['scope_effect_or_clinical_read'] is False
    d={'status':'ACCEPTED_ALL_AUTHOR_LINEAGE_SCOPES_MEMBERSHIPS_CASE_IDENTITIES_AND_PROSPECTIVE6_NO_NEW_EFFECT','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'producer_terminal_sha256':sha(src/'private_terminal.json'),'positions':len(scope),'eligible_pools':len(pools),'membership_rows':len(saved['global_row']),'unique_cells':len(np.unique(saved['global_row'])),'primary_paired_case_units':n,'reviewed_raw_source_patient_keys':reviewed.raw_patient.nunique(),'no_new_expression_scores_or_association':True}
    (out/'private_terminal.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'status':d['status'],'terminal_sha256':sha(out/'private_terminal.json'),'primary_paired_case_units':n}))


if __name__=='__main__':main()

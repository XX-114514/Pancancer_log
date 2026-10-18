"""Freeze complete author-lineage assay and prospective paired-tissue questions."""
import hashlib
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

SCOPES=['legacy_capture_anchor','author_confirmed_capture','all_author_lineage']
VIEWS=['R_FOXP3_capture','E_effector_capture']


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,16*1024**2)
    side=root/'GSE176021_exact_original_author_annotation_sidecar_20261019_v2';audit=root/'GSE176021_exact_original_author_annotation_audit_20261019_v2/private_terminal.json';assert sha(audit)=='d8d21a1c27465d9b8248857c0028e5106509e43426cbdfef9f44ccb1d71c4169'
    native=root/'immune_RNA_MCP_all_atlas_native_measurement_20261019_v2';na=root/'immune_RNA_MCP_all_atlas_native_measurement_audit_20261019_v1/private_terminal.json';assert sha(na)=='f7a6db048cfccac2a98de84023ce5048d4faa595353d2589767718d0e661cfc0'
    for directory,expected in [(side,'b4960a3c2d69256c669ca094540466f513480081a403c6499043c2c4bebaf80c'),(native,'2463d66299d7b1f1be96d99a19554a834735e688440716edfa00cd1aa01575e9')]:
        assert sha(directory/'private_terminal.json')==expected
        for f,v in json.loads((directory/'private_terminal.json').read_text())['files'].items():assert sha(directory/f)==v['sha256']
    out=root/'GSE176021_author_lineage_RNA_scope_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan={'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'resource_guard':guard,'annotation_acceptance_sha256':sha(audit),'native_acceptance_sha256':sha(na),'scopes':SCOPES,'primary':'all_author_lineage; original author CD3 CD4-Treg vs union of all original CD8-prefixed labels. These unions are unchanged by the published figure swaps. This is a reference-qualified assay, not authority/state replacement.','sensitivity':'author_confirmed_capture: intersection of original frozen exclusive marker captures with same author lineages. legacy_capture_anchor exact original primary/hash100/allcapture anchors only; no retesting closed family.','all_source_scope':'Every existing GSE176021 material and all tissues; both frozen programmes/full and fixed15stress-complement. Missing author/low cell/gene failures remain explicit. No outcome or new programme score read here.','selection':'Original 2026101802|original_view|global_row SHA order, min100, each scope/material/view >=20 author/capture cells and original capture >=20; all-captured mean sensitivity. No new gene/model selection.','background':'Primary preserve original capture-defined fraction/logN per material to retain meanings of learned nuisances; per-cell log RNA/detected genes use original training normalization. Own-population fraction/logN is a separately labelled constant-offset sensitivity; never substitute silently. Unknown source offset0/intercept cancels within same-view paired differences; no calibrated cross-source level claim.','prospective_family':'All6 primary all_author_lineage paired tumour-minus-adjacentnormal slots: 2programmes x (R shift,E shift,R/E cochange), on one shared set of case units complete for both views in both tissues, >=10. Collapse materials equally within source/case/tissue; source-author case identities/processing limitations retained. No MPR effect.','paired_tests':'Paired shifts: two-sided exact Wilcoxon signed-rank sign enumeration when <=20 nonzero differences, else19999 sign randomizations; effect median paired delta in original discovery-case SD units. Symmetry assumption explicit; exact direction-only binomial sign sensitivity separately BH4. Cochange: Spearman case permutation19999. Main6 BH/BY includes every preeligible slot, constant/unavailable numerical designsP1;2000 whole-pair bootstrap,CI>=1900valid,LOO,MC upper. No covariate deletion or fully8QC-adjusted claim at small N.','predefined_point_sensitivities':['author_confirmed_capture_on_same_cases','all_capture_mean','fixed15_stress_complement','own_population_background','unresidualized_weighted_RNA','RNA_quality_only_background','drop_one_case'],'coactivity_followup':'All scopes/tissues/both programmes R/E coactivity point/LOO and composition descriptions only; no new P or significance re-screen for the closed original-source family.','limits':'Source-context biological association, not direct communication, clinical prediction or independent replication. Author labels/mapping do not measure antigen specificity. No new expression download.'}
    (out/'plan_before_author_lineage_RNA_or_tissue_effects.json').write_text(json.dumps(plan,indent=2)+'\n')
    lib=pd.read_csv(native/'private_all400_materials.tsv',sep='\t',keep_default_na=False);lib=lib[lib.dataset.eq('GSE176021')].reset_index(drop=True);lookup=pd.Index(lib.sampleID);assert lookup.is_unique;lib['new_material_index']=np.arange(len(lib))
    z=zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),mode='r');read=lambda f:np.asarray(read_elem(z['obs/'+f])).astype(str);source=read('dataset');global_rows=np.flatnonzero(source=='GSE176021');sids=np.char.add('GSE176021::',read('library_id')[global_rows]);mi=lookup.get_indexer(sids);assert (mi>=0).all() and np.array_equal(np.bincount(mi,minlength=len(lib)),lib.annotated_scope_cells.to_numpy(int))
    bridge=np.load(side/'private_all_atlas_cells_exact_author_bridge.npz');assert np.array_equal(global_rows,bridge['global_atlas_row']);labels=json.loads((side/'private_author_label_lookup_and_coverage.json').read_text())['CD3']['labels'];codes=bridge['CD3_label_code'];author=[codes==labels.index('CD4-Treg'),np.isin(codes,[i for i,n in enumerate(labels) if n.startswith('CD8-')])];assert not (author[0]&author[1]).any()
    old=root/'immune_CM_fixed_panel_same_cell_RNA_20261018_v1';cache=np.load(old/'private_fixed34_gene_T_counts_and_depth.npz');mark=cache['genes'].tolist();counts=cache['raw_gene_counts'];fixed=root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json';panels=json.loads(fixed.read_text())['panels'];positive=lambda names:np.count_nonzero(counts[:,[mark.index(g) for g in names]]>0,axis=1);core=positive(panels['pan_T'])>=2;reg=(positive(panels['regulatory_associated'])>=2)&(counts[:,mark.index('FOXP3')]>0);eff=positive(['GZMB','PRF1','CTSW','NKG7'])>=2;rawT=cache['global_original_T_row_index'];legacy=[]
    for mask in [core&reg&~eff,core&eff&~reg]:legacy.append(np.isin(global_rows,rawT[mask]))
    rows=[];membership=[];scope=[];pairs=[]
    oldp=pd.read_csv(native/'private_all_capture_pools.tsv',sep='\t');oldmem=np.load(native/'private_all_captured_membership.npz')
    for si,name in enumerate(SCOPES):
        for m,r in lib.iterrows():
            for vi,view in enumerate(VIEWS):
                take=(mi==m)&(legacy[vi] if si==0 else legacy[vi]&author[vi] if si==1 else author[vi]);gr=global_rows[take];oldcount=int(r[view]);n=len(gr);eligible=n>=20 and oldcount>=20
                if si==0:assert n==oldcount
                record={'scope':name,'scope_index':si,'new_material_index':m,'sampleID':r.sampleID,'view':view,'view_index':vi,'cells':n,'original_capture_cells':oldcount,'annotated_scope_cells':int(r.annotated_scope_cells),'eligible':eligible,'status':'ELIGIBLE_METADATA' if eligible else 'NEW_SCOPE_FEWER_THAN20' if n<20 else 'ORIGINAL_CAPTURE_FEWER_THAN20_BACKGROUND_NOT_QUALIFIED'}
                if eligible:
                    p=len(rows);record['pool_index']=p;selection=set(sorted(gr,key=lambda x:hashlib.sha256(('2026101802|'+view+'|'+str(int(x))).encode()).hexdigest())[:100]);rows.append({**record,'selected_cells':len(selection)});membership.extend((int(x),p,int(x) in selection) for x in gr)
                    if si==0:
                        op=oldp[(oldp.sampleID==r.sampleID)&(oldp['view']==view)];assert len(op)==1;opid=int(op.pool_index.iloc[0]);assert np.array_equal(gr,oldmem['global_row'][oldmem['pool_index']==opid])
                        assert set(oldmem['global_row'][(oldmem['pool_index']==opid)&oldmem['selected']])==selection
                else:record['pool_index']=-1
                scope.append(record)
            aa=scope[-2:];pairs.append({'scope':name,'sampleID':r.sampleID,'case_block':r.case_block,'reviewed_patient':r.reviewed_patient,'tissue':r.tissue,'phase_class':r.phase_class,'capture_scope':r.capture_scope,'ex_vivo':r.ex_vivo,'both_views_eligible':bool(all(v['eligible'] for v in aa)),'R_cells':aa[0]['cells'],'E_cells':aa[1]['cells']})
    pools=pd.DataFrame(rows);members=pd.DataFrame(membership,columns=['global_row','pool_index','selected']).sort_values(['global_row','pool_index']).reset_index(drop=True);assert not members.duplicated(['global_row','pool_index']).any();scope=pd.DataFrame(scope);pairs=pd.DataFrame(pairs)
    gene_roles=pd.read_csv(root/'immune_RNA_MCP_all_atlas_projection_contract_20261019_v1/private_frozen_all_programme_LR_and_diagnostic_gene_roles.tsv',sep='\t');gm=pd.read_csv(native/'private_all_source_auxiliary_native_gene_mapping.tsv.gz',sep='\t');gm=gm[gm.dataset.eq('GSE176021')].set_index('gene').loc[gene_roles.gene].reset_index();weights=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv',sep='\t');w=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')&weights.view.isin(VIEWS)];assert gm.set_index('gene').loc[w.gene,'native_column'].ge(0).all()
    counts_out=[];paired=[]
    for name in SCOPES:
        d=pairs[(pairs.scope==name)&pairs.both_views_eligible]
        for t,g in d.groupby('tissue'):counts_out.append({'scope':name,'tissue':t,'complete_materials':len(g),'case_units':g.case_block.nunique(),'reviewed_case_units':g[g.reviewed_patient.astype(str).str.lower().isin(['true','1'])].case_block.nunique()})
        by=d.groupby('case_block').tissue.apply(set);shared=sorted(c for c,s in by.items() if {'tumor','adjacent normal lung'}<=s)
        for c in shared:paired.append({'scope':name,'case_block':c,'both_views_both_tissues':True,'tumour_materials':int(((d.case_block==c)&d.tissue.eq('tumor')).sum()),'normal_materials':int(((d.case_block==c)&d.tissue.eq('adjacent normal lung')).sum())})
        counts_out.append({'scope':name,'tissue':'PAIRED_TUMOUR_ADJACENT_NORMAL','complete_materials':int(d[d.case_block.isin(shared)&d.tissue.isin(['tumor','adjacent normal lung'])].shape[0]),'case_units':len(shared),'reviewed_case_units':int(d[d.case_block.isin(shared)&d.reviewed_patient.astype(str).str.lower().isin(['true','1'])].case_block.nunique())})
    unique_cells=members.global_row.nunique();selected=int(members.selected.sum());npools=len(pools);bound=int(len(members)*11+selected*112+npools*len(gm)*48+npools*16384+unique_cells*12+8*1024**2);assert bound<64*1024**2
    registry=[];npaired=sum(x['scope']=='all_author_lineage' for x in paired)
    for program in ['MCP1','MCP2']:
        for outcome in ['R_tissue_shift','E_tissue_shift','RE_tissue_cochange']:registry.append({'model_id':'AP'+str(len(registry)+1).zfill(2),'scope':'all_author_lineage','program':program,'outcome':outcome,'paired_case_units':npaired,'metadata_preeligible':npaired>=10,'family':'paired_tissue_main6'})
    lib.to_csv(out/'private_all_source_materials.tsv',sep='\t',index=False);scope.to_csv(out/'private_all_scope_material_view_positions.tsv',sep='\t',index=False);pools.to_csv(out/'private_all_eligible_pools.tsv',sep='\t',index=False);pairs.to_csv(out/'private_all_scope_material_pair_positions.tsv',sep='\t',index=False);pd.DataFrame(paired).to_csv(out/'private_all_paired_case_positions.tsv',sep='\t',index=False);pd.DataFrame(counts_out).to_csv(out/'private_all_source_context_coverage.tsv',sep='\t',index=False);pd.DataFrame(registry).to_csv(out/'private_all6_paired_tissue_registered_tests.tsv',sep='\t',index=False);gm.to_csv(out/'private_all904_native_gene_roles_and_mapping.tsv',sep='\t',index=False);np.savez_compressed(out/'private_all_scope_cell_memberships.npz',global_row=members.global_row.to_numpy(np.int32),pool_index=members.pool_index.to_numpy(np.int16),selected=members.selected.to_numpy(bool))
    files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<16*1024**2
    receipt={'status':'COMPLETE_AUTHOR_LINEAGE_ALL_MATERIAL_SCOPE_AND_PAIRED_FAMILY_FROZEN_NO_NEW_RNA_EFFECT','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'files':files,'positions':len(scope),'eligible_pools':npools,'membership_rows':len(members),'unique_native_cells':unique_cells,'selected_memberships':selected,'native_prospective_cap_bytes':64*1024**2,'conservative_native_uncompressed_bound_bytes':bound,'paired_primary_case_units':npaired,'scope_effect_or_clinical_read':False}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json'),'eligible_pools':npools,'paired_primary_case_units':npaired,'uncompressed_bound_bytes':bound}))


if __name__=='__main__':main()

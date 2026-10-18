"""All existing atlas materials, frozen pairwise RNA captures and full gene axes."""
import hashlib
import json
import os
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

VIEWS=['R_FOXP3_capture','E_effector_capture','Myeloid']
PAIRS=[(VIEWS[0],VIEWS[1]),(VIEWS[0],VIEWS[2]),(VIEWS[1],VIEWS[2])]


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);cap=16*1024**2;guard=resource(root,cap)
    scope=root/'immune_RNA_MCP_discovery_scope_20261018_v1';assert sha(scope/'private_terminal.json')=='0f4833a07bd5cdf5e8f6806d06805777dddf78adad60a6151eca991dfdda6196'
    for n,v in json.loads((scope/'private_terminal.json').read_text())['files'].items():assert sha(scope/n)==v['sha256']
    previous_plan=json.loads((scope/'plan_before_expression_or_MCP_fitting.json').read_text())
    libpath=root/'NSCLC_native_MFA_external_lineage_qualification_20261016_v2/private_all1690_library_metadata_and_label_diagnostics.tsv';assert sha(libpath)==previous_plan['current_library_metadata_sha256']
    fields=['arm','dataset','library_id','reviewed_patient','patient_block','resampling_block','patient_evidence','phase_class','cancer_label','tissue','material_class','capture_scope','ex_vivo','annotated_scope_cells','latest_sidecar']
    lib=pd.read_csv(libpath,sep='\t',usecols=fields,keep_default_na=False);lib=lib[lib.arm.eq('archived_immune_atlas')].copy().sort_values(['dataset','library_id']).reset_index(drop=True)
    lib['sampleID']=lib.dataset.astype(str)+'::'+lib.library_id.astype(str);lib['reviewed_patient']=lib.reviewed_patient.astype(str).str.lower().eq('true');assert len(lib)==400 and lib.sampleID.is_unique and lib.dataset.nunique()==9
    old=root/'immune_CM_fixed_panel_same_cell_RNA_20261018_v1';audit=root/'immune_CM_fixed_panel_same_cell_RNA_audit_20261018_v1/private_terminal.json';proof=json.loads(audit.read_text());assert proof['status'].startswith('ACCEPTED_ALL_ORIGINAL_T_CELLS_') and sha(old/'private_terminal.json')==proof['producer_terminal_sha256']
    for n,v in json.loads((old/'private_terminal.json').read_text())['files'].items():assert sha(old/n)==v['sha256']
    oldscope=root/'existing_immune_fine_native_measurement_scope_20261016_v1'
    for n,v in json.loads((oldscope/'private_terminal.json').read_text())['files'].items():assert sha(oldscope/n)==v
    ready=root/'immune_CM_RNA_supported_multicell_views_preflight_20261018_v1';assert sha(ready/'private_independent_count_readiness_acceptance_20261018_v2.json')=='3f1261b9083f1cefcf0e072daecdef7c3db5d1fbf7f44b4ac1ff237efb126c15'
    original=pd.read_csv(ready/'private_all_original_libraries_multicell_view_counts.tsv',sep='\t');original=original[original.threshold.eq(1)].set_index('sampleID');assert len(original)==118
    chosen=pd.read_csv(scope/'private_metadata_only_case_source_partition.tsv',sep='\t',keep_default_na=False)
    discovery=chosen[chosen.partition.eq('reviewed_case_discovery')]
    lib['case_block']=lib.dataset.astype(str)+'::'+lib.resampling_block.astype(str)
    lib['in_original_CM118']=lib.sampleID.isin(original.index);lib['in_original_RNA_material_selection']=lib.sampleID.isin(chosen.sampleID)
    lib['original_discovery_material']=lib.sampleID.isin(discovery.sampleID);lib['reviewed_discovery_case_overlap']=lib.reviewed_patient&lib.case_block.isin(discovery.case_block)
    paired=root/'HNSCC_immune_RNA_MCP_paired_native_measurement_20261019_v1/private_all52_paired_materials.tsv'
    assert sha(root/'HNSCC_immune_RNA_MCP_paired_native_measurement_audit_20261019_v1/private_terminal.json')=='f3c916ecb8f4d40b535ca0c4486a0ae8a825dbff1d362ae6b64e721b5c981640'
    prior=set(chosen.loc[~chosen.partition.eq('small_source_projection_only'),'sampleID'])|set(pd.read_csv(paired,sep='\t',usecols=['sampleID']).sampleID)
    lib['previously_RNA_measured_material']=lib.sampleID.isin(prior)
    out=root/'immune_RNA_MCP_all_atlas_opportunities_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    fixed=root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json';panels=json.loads(fixed.read_text())['panels']
    (out/'plan_before_all400_capture_and_gene_opportunities.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,metadata_sha256=sha(libpath),old_all_T_RNA_acceptance_sha256=sha(audit),original400_scope='All archived immune atlas materials, including non-primary/T-enriched/unknown contexts. No new scores, clinical outcomes, patient promotions or fitting.',
        capture='Same threshold1/core>=100; per pair both>=20 and source capture genes measured. Missing third view does not invalidate a measurable pair. Group by source/cancer/tissue/material/phase/capture/exvivo.',
        native_axes='All6 fixed programmes/all9sources, exact unique native symbol/canonical/stableID alias; no missing zero fill/gene trimming, forbid within-axis physical collisions or capture-marker alias overlaps.',
        prior_data='CM118, original RNA material/case overlap and all already measured original66/paired52 marked. A new material is not a new patient.',
        future='After independent entire opportunity audit, freeze full usable source/programme/pair slots and designs before new effects. No previously closed family reopened.'),indent=2)+'\n')
    T=pd.read_csv(old/'private_original_T_unit_scope.tsv',sep='\t',keep_default_na=False);cache=np.load(old/'private_fixed34_gene_T_counts_and_depth.npz');raw=cache['raw_gene_counts'];genes=cache['genes'].tolist();tu=cache['original_T_unit_index']
    keys=pd.Index(lib.sampleID);t_material=keys.get_indexer(T.dataset.astype(str)+'::'+T.library_id.astype(str))[tu];assert (t_material>=0).all()
    contracts={r['dataset']:r for r in json.loads((old/'private_source_gene_measurement_contract.json').read_text())}
    complete=np.array([all(contracts[s]['complete_panels'][k] for k in ['pan_T','regulatory_associated','cytotoxic']) for s in T.dataset])[tu]
    core=(raw[:,[genes.index(g) for g in panels['pan_T']]]>=1).sum(axis=1)>=2
    reg=((raw[:,[genes.index(g) for g in panels['regulatory_associated']]]>=1).sum(axis=1)>=2)&(raw[:,genes.index('FOXP3')]>=1)
    eff=(raw[:,[genes.index(g) for g in ['GZMB','PRF1','CTSW','NKG7']]]>=1).sum(axis=1)>=2
    lib['capture_gene_panels_complete']=[all(contracts[s]['complete_panels'][k] for k in ['pan_T','regulatory_associated','cytotoxic']) for s in lib.dataset]
    lib['core_T']=np.bincount(t_material[core&complete],minlength=400)
    lib[VIEWS[0]]=np.bincount(t_material[core&reg&~eff&complete],minlength=400)
    lib[VIEWS[1]]=np.bincount(t_material[core&eff&~reg&complete],minlength=400)
    units=pd.read_csv(oldscope/'private_all_frozen_T_myeloid_node_units.tsv',sep='\t',keep_default_na=False)
    my=units[units.major.eq('Myeloid')].copy();my['sampleID']=my.dataset.astype(str)+'::'+my.library_id.astype(str);lib[VIEWS[2]]=my.groupby('sampleID').n_cells.sum().reindex(lib.sampleID,fill_value=0).to_numpy(int)
    assert (lib[VIEWS].sum(axis=1)<=lib.annotated_scope_cells).all()
    checks=lib[lib.in_original_CM118].set_index('sampleID').loc[original.index]
    for new,oldfield in [('core_T','core_T'),(VIEWS[0],'FOXP3_supported_R_only'),(VIEWS[1],'effector_only'),(VIEWS[2],'original_myeloid_cells')]:assert np.array_equal(checks[new],original[oldfield])
    atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1';z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r');assert sha(atlas/'manifests/immune_authority_release_public.json')==previous_plan['atlas_release_sha256']
    order=np.asarray(read_elem(z['uns/source_order'])).astype(str).tolist();measured=np.asarray(read_elem(z['varm/measured_in_source']),bool);symbol=np.asarray(read_elem(z['var/gene_symbol'])).astype(str);canonical=np.asarray(read_elem(z['var/canonical_gene_identity'])).astype(str);identity=np.asarray(read_elem(z['var/feature_identity_status'])).astype(str);aliases=np.asarray(read_elem(z['var/source_gene_aliases'])).astype(str)
    wp=root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv';assert sha(wp)=='fd7316c67beb687c403a752c36b5989644e432f59a83119b21eade5efb0d89a7';weights=pd.read_csv(wp,sep='\t');weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')]
    oldmap=pd.read_csv(scope/'private_marker_excluded_native_gene_map.tsv',sep='\t').set_index('gene');maps=[];axes=[]
    for source in sorted(lib.dataset.unique()):
        si=order.index(source);names=defaultdict(set)
        for c in np.flatnonzero(measured[:,si]):
            if canonical[c].startswith('AMBIGUOUS::') or identity[c]=='ambiguous_source_feature':continue
            labels={symbol[c],canonical[c]}
            if identity[c]=='stable_Ensembl_ID':labels.update(aliases[c].split('|'))
            for name in labels:
                if re.fullmatch(r'[A-Za-z][A-Za-z0-9._-]*',name) and not name.startswith('ENSG') and name not in ['nan','None','__SOURCE_VALUE_MISSING__']:names[name].add(int(c))
        markercolumns={c for c in contracts[source]['physical_columns'] if c>=0}
        for (view,program),w in weights.groupby(['view','program'],sort=True):
            columns=[next(iter(names[g])) if len(names[g])==1 else -1 for g in w.gene];collisions={c for c in columns if c>=0 and columns.count(c)>1};marker_hits={c for c in columns if c in markercolumns}
            for gene,col in zip(w.gene,columns):
                if source in oldmap.columns:assert col==int(oldmap.loc[gene,source])
                maps.append(dict(dataset=source,view=view,program=program,gene=gene,native_column=col,physical_collision=col in collisions,group_marker_physical_overlap=col in marker_hits))
            axes.append(dict(dataset=source,view=view,program=program,active_positions=len(w),uniquely_measured_positions=sum(c>=0 for c in columns),physical_collision_positions=sum(c in collisions for c in columns),group_marker_overlap_positions=sum(c in marker_hits for c in columns),complete_axis=all(c>=0 for c in columns) and not collisions and not marker_hits))
    axis=pd.DataFrame(axes);assert len(axis)==54;axis.to_csv(out/'private_all54_source_axis_qualification.tsv',sep='\t',index=False);pd.DataFrame(maps).to_csv(out/'private_all_active_gene_source_native_mapping.tsv.gz',sep='\t',index=False,compression='gzip')
    lib.to_csv(out/'private_all400_material_RNA_capture_scope.tsv',sep='\t',index=False)
    byaxis=axis.set_index(['dataset','view','program']).complete_axis.to_dict();pairrows=[]
    for r in lib.itertuples(index=False):
        for va,vb in PAIRS:
            eligible=bool(r.capture_gene_panels_complete and r.core_T>=100 and getattr(r,va)>=20 and getattr(r,vb)>=20)
            for p in ['MCP1','MCP2']:
                pairrows.append(dict(sampleID=r.sampleID,dataset=r.dataset,pair=va+'|'+vb,program=p,capture_eligible=eligible,complete_both_axes=bool(byaxis[(r.dataset,va,p)] and byaxis[(r.dataset,vb,p)]),eligible=bool(eligible and byaxis[(r.dataset,va,p)] and byaxis[(r.dataset,vb,p)])))
    pair=pd.DataFrame(pairrows);assert len(pair)==2400;pair.to_csv(out/'private_all2400_material_programme_pair_qualification.tsv.gz',sep='\t',index=False,compression='gzip')
    full=pair.merge(lib,on=['sampleID','dataset'],validate='many_to_one');summary=[];context=[]
    def summarize(g):
        q=g[g.eligible];reviewed=q[q.reviewed_patient]
        return dict(all_materials=len(g),capture_eligible_materials=int(g.capture_eligible.sum()),full_gene_eligible_materials=len(q),eligible_reviewed_case_blocks=int(reviewed.case_block.nunique()),eligible_reviewed_materials=len(reviewed),eligible_unreviewed_sample_units=int((~q.reviewed_patient).sum()),eligible_outside_CM118=int((~q.in_original_CM118).sum()),eligible_previously_unmeasured_materials=int((~q.previously_RNA_measured_material).sum()),eligible_reviewed_non_discovery_cases=int(reviewed.loc[~reviewed.reviewed_discovery_case_overlap,'case_block'].nunique()),coverage_opportunity_not_model_or_independence_acceptance=True)
    for (source,pairname,p),g in full.groupby(['dataset','pair','program'],sort=True):summary.append(dict(dataset=source,pair=pairname,program=p,**summarize(g)))
    contextfields=['dataset','cancer_label','tissue','material_class','phase_class','capture_scope','ex_vivo','pair','program']
    for key,g in full.groupby(contextfields,sort=True,dropna=False):context.append(dict(zip(contextfields,key),**summarize(g)))
    assert len(summary)==54;pd.DataFrame(summary).to_csv(out/'private_all54_source_programme_pair_opportunities.tsv',sep='\t',index=False);pd.DataFrame(context).to_csv(out/'private_all_context_programme_pair_opportunities.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()};assert sum(v['bytes'] for v in files.values())<cap
    t=dict(status='COMPLETE_ALL400_MATERIAL_CAPTURE_ALL54_FULL_AXES_AND_PAIR_OPPORTUNITIES_NO_EFFECTS_AUDIT_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,all_materials=400,source_axes=54,material_programme_pairs=2400,original118_capture_counts_exact=True,no_new_expression_or_scores_or_effects=True)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

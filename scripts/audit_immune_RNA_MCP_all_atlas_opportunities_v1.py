"""All400 capture counts,54 native axes and2400 pair/context slots replayed."""
import json
import os
from collections import Counter,defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

VIEWS=['R_FOXP3_capture','E_effector_capture','Myeloid']


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    source=root/'immune_RNA_MCP_all_atlas_opportunities_20261019_v1';assert sha(source/'private_terminal.json')=='2a37f7af2f73ac5ceddcb065626df411ca6f016271fce7f0830a1bd658f83cd1'
    term=json.loads((source/'private_terminal.json').read_text());assert term['source_sha256']=='0a0abfac3bb38a730b3b6a9e60cc3430bf51122cd33aa12d3e07abee2682de56'
    for n,v in term['files'].items():assert sha(source/n)==v['sha256']
    out=root/'immune_RNA_MCP_all_atlas_opportunities_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_all400_capture_axis_context_replay.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,producer_terminal_sha256=sha(source/'private_terminal.json'),scope='All400 materials via per-T-node four-state contingency then library sums; raw original-cell map independently recounts My units; all54 native axes rebuilt using required-name-only maps; all2400 eligibility records/source/context sets checked. No scores or new effects.'),indent=2)+'\n')
    lib=pd.read_csv(source/'private_all400_material_RNA_capture_scope.tsv',sep='\t',keep_default_na=False);assert len(lib)==400 and lib.sampleID.is_unique
    scope=root/'immune_RNA_MCP_discovery_scope_20261018_v1';priorplan=json.loads((scope/'plan_before_expression_or_MCP_fitting.json').read_text())
    lp=root/'NSCLC_native_MFA_external_lineage_qualification_20261016_v2/private_all1690_library_metadata_and_label_diagnostics.tsv';assert sha(lp)==priorplan['current_library_metadata_sha256']
    m=pd.read_csv(lp,sep='\t',keep_default_na=False);m=m[m.arm.eq('archived_immune_atlas')].copy();m['sampleID']=m.dataset.astype(str)+'::'+m.library_id.astype(str);m=m.set_index('sampleID').loc[lib.sampleID]
    fields=['dataset','library_id','patient_block','resampling_block','patient_evidence','phase_class','cancer_label','tissue','material_class','capture_scope','ex_vivo','annotated_scope_cells','latest_sidecar']
    for f in fields:assert np.array_equal(lib[f].astype(str),m[f].astype(str))
    assert np.array_equal(lib.reviewed_patient,m.reviewed_patient.astype(str).str.lower().eq('true'))
    chosen=pd.read_csv(scope/'private_metadata_only_case_source_partition.tsv',sep='\t',keep_default_na=False);discovery=chosen[chosen.partition.eq('reviewed_case_discovery')];cb=lib.dataset.astype(str)+'::'+lib.resampling_block.astype(str);assert np.array_equal(lib.case_block,cb)
    original=pd.read_csv(root/'immune_CM_RNA_supported_multicell_views_preflight_20261018_v1/private_all_original_libraries_multicell_view_counts.tsv',sep='\t');original=original[original.threshold.eq(1)].set_index('sampleID')
    previous=set(chosen.loc[~chosen.partition.eq('small_source_projection_only'),'sampleID'])|set(pd.read_csv(root/'HNSCC_immune_RNA_MCP_paired_native_measurement_20261019_v1/private_all52_paired_materials.tsv',sep='\t').sampleID)
    expected_flags={'in_original_CM118':lib.sampleID.isin(original.index),'in_original_RNA_material_selection':lib.sampleID.isin(chosen.sampleID),'original_discovery_material':lib.sampleID.isin(discovery.sampleID),'reviewed_discovery_case_overlap':lib.reviewed_patient&cb.isin(discovery.case_block),'previously_RNA_measured_material':lib.sampleID.isin(previous)}
    for f,v in expected_flags.items():assert np.array_equal(lib[f],v)
    old=root/'immune_CM_fixed_panel_same_cell_RNA_20261018_v1';cache=np.load(old/'private_fixed34_gene_T_counts_and_depth.npz');raw=cache['raw_gene_counts'];genes=cache['genes'].tolist();tu=cache['original_T_unit_index'];T=pd.read_csv(old/'private_original_T_unit_scope.tsv',sep='\t',keep_default_na=False)
    contracts={r['dataset']:r for r in json.loads((old/'private_source_gene_measurement_contract.json').read_text())};fixed=root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json';panels=json.loads(fixed.read_text())['panels']
    complete=np.array([all(contracts[s]['complete_panels'][k] for k in ['pan_T','regulatory_associated','cytotoxic']) for s in T.dataset]);cc=complete[tu]
    core=np.count_nonzero(raw[:,[genes.index(g) for g in panels['pan_T']]]>0,axis=1)>=2
    reg=(raw[:,genes.index('FOXP3')]>0)&(np.count_nonzero(raw[:,[genes.index(g) for g in panels['regulatory_associated']]]>0,axis=1)>=2)
    eff=np.count_nonzero(raw[:,[genes.index(g) for g in ['GZMB','PRF1','CTSW','NKG7']]]>0,axis=1)>=2
    state=reg.astype(int)+2*eff.astype(int);take=core&cc
    totals=np.bincount(tu[take]*4+state[take],minlength=len(T)*4).reshape(len(T),4);T['sampleID']=T.dataset.astype(str)+'::'+T.library_id.astype(str)
    for i in range(4):T['state'+str(i)]=totals[:,i]
    grouped=T.groupby('sampleID')[['state0','state1','state2','state3']].sum().reindex(lib.sampleID,fill_value=0)
    assert np.array_equal(lib.core_T,grouped.sum(axis=1)) and np.array_equal(lib[VIEWS[0]],grouped.state1) and np.array_equal(lib[VIEWS[1]],grouped.state2)
    assert np.array_equal(lib.capture_gene_panels_complete,[all(contracts[s]['complete_panels'][k] for k in ['pan_T','regulatory_associated','cytotoxic']) for s in lib.dataset])
    oldscope=root/'existing_immune_fine_native_measurement_scope_20261016_v1';units=pd.read_csv(oldscope/'private_all_frozen_T_myeloid_node_units.tsv',sep='\t',keep_default_na=False);mapping=np.load(oldscope/'private_original_cell_to_fine_node_unit.npy',mmap_mode='r');counts=np.bincount(mapping[mapping>=0],minlength=len(units));assert np.array_equal(counts,units.set_index('unit_index').n_cells.reindex(np.arange(len(units))))
    my=Counter()
    for r in units[units.major.eq('Myeloid')].itertuples(index=False):my[str(r.dataset)+'::'+str(r.library_id)]+=int(counts[r.unit_index])
    assert np.array_equal(lib[VIEWS[2]],[my[s] for s in lib.sampleID])
    overlap=lib[lib.in_original_CM118].set_index('sampleID').loc[original.index]
    for f,g in [('core_T','core_T'),(VIEWS[0],'FOXP3_supported_R_only'),(VIEWS[1],'effector_only'),(VIEWS[2],'original_myeloid_cells')]:assert np.array_equal(overlap[f],original[g])
    atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1';z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r')
    order=np.asarray(read_elem(z['uns/source_order'])).astype(str).tolist();measured=np.asarray(read_elem(z['varm/measured_in_source']),bool);symbol=np.asarray(read_elem(z['var/gene_symbol'])).astype(str);canonical=np.asarray(read_elem(z['var/canonical_gene_identity'])).astype(str);identity=np.asarray(read_elem(z['var/feature_identity_status'])).astype(str);aliases=np.asarray(read_elem(z['var/source_gene_aliases'])).astype(str)
    w=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv',sep='\t');w=w[w.fit_id.eq('K2_seed20261018_source_adjusted')];requested=set(w.gene)
    maps=pd.read_csv(source/'private_all_active_gene_source_native_mapping.tsv.gz',sep='\t').set_index(['dataset','view','program','gene']);axis=pd.read_csv(source/'private_all54_source_axis_qualification.tsv',sep='\t').set_index(['dataset','view','program']);assert len(axis)==54
    complete_axes={};map_records=0
    for dataset in sorted(lib.dataset.unique()):
        matches={g:set() for g in requested};si=order.index(dataset)
        for col in np.flatnonzero(measured[:,si]):
            if identity[col]=='ambiguous_source_feature' or canonical[col].startswith('AMBIGUOUS::'):continue
            labels={symbol[col],canonical[col]}
            if identity[col]=='stable_Ensembl_ID':labels|=set(aliases[col].split('|'))
            for g in labels&requested:matches[g].add(int(col))
        forbidden={c for c in contracts[dataset]['physical_columns'] if c>=0}
        for (v,p),gw in w.groupby(['view','program'],sort=True):
            cols={g:next(iter(matches[g])) if len(matches[g])==1 else -1 for g in gw.gene};freq=Counter(c for c in cols.values() if c>=0);collision=sum(c>=0 and freq[c]>1 for c in cols.values());marker=sum(c in forbidden for c in cols.values());qualified=all(c>=0 for c in cols.values()) and collision==0 and marker==0
            r=axis.loc[(dataset,v,p)];assert r.active_positions==len(gw) and r.uniquely_measured_positions==sum(c>=0 for c in cols.values()) and r.physical_collision_positions==collision and r.group_marker_overlap_positions==marker and r.complete_axis==qualified
            complete_axes[(dataset,v,p)]=qualified
            for g,c in cols.items():
                q=maps.loc[(dataset,v,p,g)];assert q.native_column==c and q.physical_collision==bool(c>=0 and freq[c]>1) and q.group_marker_physical_overlap==bool(c in forbidden);map_records+=1
    assert map_records==len(maps)
    pair=pd.read_csv(source/'private_all2400_material_programme_pair_qualification.tsv.gz',sep='\t');assert len(pair)==2400 and not pair.duplicated(['sampleID','pair','program']).any();lm=lib.set_index('sampleID')
    for r in pair.itertuples(index=False):
        a,b=r.pair.split('|');m=lm.loc[r.sampleID];eligible=bool(m.capture_gene_panels_complete and m.core_T>=100 and m[a]>=20 and m[b]>=20);genesok=bool(complete_axes[(r.dataset,a,r.program)] and complete_axes[(r.dataset,b,r.program)])
        assert r.capture_eligible==eligible and r.complete_both_axes==genesok and r.eligible==bool(eligible and genesok)
    full=pair.merge(lib,on=['sampleID','dataset'],validate='many_to_one')
    def expected(g):
        valid=g[g.eligible];rv=valid[valid.reviewed_patient]
        return dict(all_materials=len(g),capture_eligible_materials=int(g.capture_eligible.sum()),full_gene_eligible_materials=len(valid),eligible_reviewed_case_blocks=len(set(rv.case_block)),eligible_reviewed_materials=len(rv),eligible_unreviewed_sample_units=int((~valid.reviewed_patient).sum()),eligible_outside_CM118=int((~valid.in_original_CM118).sum()),eligible_previously_unmeasured_materials=int((~valid.previously_RNA_measured_material).sum()),eligible_reviewed_non_discovery_cases=len(set(rv.loc[~rv.reviewed_discovery_case_overlap,'case_block'])))
    nchecked=0
    for filename,fields in [('private_all54_source_programme_pair_opportunities.tsv',['dataset','pair','program']),('private_all_context_programme_pair_opportunities.tsv',['dataset','cancer_label','tissue','material_class','phase_class','capture_scope','ex_vivo','pair','program'])]:
        table=pd.read_csv(source/filename,sep='\t',keep_default_na=False).set_index(fields);groups=full.groupby(fields,sort=True,dropna=False);assert len(table)==len(groups)
        for key,g in groups:
            r=table.loc[key]
            for f,value in expected(g).items():assert r[f]==value
            assert r.coverage_opportunity_not_model_or_independence_acceptance;nchecked+=1
    t=dict(status='ACCEPTED_ALL400_CAPTURE_COUNTS_54_NATIVE_AXES_2400_PAIRS_AND_ALL_SOURCE_CONTEXT_SUMMARIES',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),materials_replayed=400,native_axes_replayed=54,native_mapping_records_replayed=map_records,material_programme_pairs_replayed=2400,source_context_summary_rows_replayed=nchecked,original118_capture_counts_exact=True,no_new_scores_or_effects=True)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

"""Finite same-case common immune background sensitivity; arrays RAM only.

Background uses agreement of frozen current coarse lineage and original
author broad annotation. Missing reference is not absence; the intersection
is an analysis view, never an annotation update. No new P/BH or clinical test.
"""
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import json,os,resource,traceback
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
import analyze_CM02_RNA_multigene_in_RAM_v1 as baseline
import CM02_shared_background_rank_v1 as sensitivity

ONTOLOGY={'B':['B cells','GCB','plasma B'],
          'Myeloid':['TAMs','mDC','pDC','Neutrophils','prolif Neu']}
sha=baseline.sha


def broad_intersection(tab,cur):
    parent=Path(os.environ['PANCANCER_PARENT_OUTPUT'])
    terminal=json.loads((parent/'private_terminal.json').read_text())
    for name,item in terminal['files'].items():assert sha(parent/name)==item['SHA256']
    reader=baseline.load(parent.parent/'probe_CM02_author_library_barcode_bridge_v1.py',terminal['source_SHA256'],'accepted_broad_reader')
    name,field,expected=reader.SOURCES[0];assert field=='manual_celltype_broad'
    op=next(r for r in json.loads((parent/'private_all18_library_barcode_opportunities.json').read_text())
        if r['annotation_source']==name and r['barcode_mode']=='within_reviewed_library_16nt' and r['current_cell_key_field']=='original_cell_id')
    lookup={}
    for row in tab.to_dict('records'):
        key=(reader.norm(row[op['patient_alias_field']]),reader.norm(row[op['time_alias_field']]))
        assert key not in lookup;lookup[key]=row['library_id']
    author=reader.read_author(name,field,expected)
    author['library_id']=[lookup.get((reader.norm(p),reader.norm(t))) for p,t in zip(author.patient,author.phase)]
    author['barcode']=[reader.barcode(v,'within_reviewed_library_16nt') for v in author.cell_key]
    author=author.loc[author.library_id.notna()&author.barcode.notna(),['library_id','barcode','label']]
    assert not author.duplicated(['library_id','barcode']).any() and not cur.duplicated(['library_id','barcode']).any()
    joined=cur.merge(author,on=['library_id','barcode'],how='inner',validate='one_to_one')
    got=Counter({tuple(map(str,k)):int(n) for k,n in joined.groupby(['timepoint','major','node_id','label']).size().items()})
    old=json.loads((parent/'private_all_source_node_state_crosscounts.json').read_text())
    expected_counts=Counter({(r['phase'],r['current_major'],r['current_node'],r['author_label']):r['matched_cells'] for r in old
        if r['annotation_source']==name and r['barcode_mode']=='within_reviewed_library_16nt' and r['current_cell_key_field']=='original_cell_id'})
    assert got==expected_counts
    coverage=[];selected=[];registry=tab.set_index('library_id')
    for phase in ['pre','post']:
        for major,labels in ONTOLOGY.items():
            original=cur.loc[cur.major.eq(major)&cur.timepoint.eq(phase)]
            matched=joined.loc[joined.major.eq(major)&joined.timepoint.eq(phase)]
            agreed=matched.loc[matched.label.isin(labels)].copy();selected.append(agreed)
            counts=agreed.library_id.map(registry.patient_registry).value_counts()
            coverage.append(dict(phase=phase,current_major=major,current_cells=len(original),broad_matched_cells=len(matched),
                concordant_cells=len(agreed),concordant_fraction_among_matched=len(agreed)/len(matched) if len(matched) else None,
                author_T_fraction_among_matched=float(matched.label.isin(['CD4+T','CD8+T','Tstr','prolif T']).mean()) if len(matched) else None,
                reference_coverage_of_current=len(matched)/len(original) if len(original) else None,
                reviewed_phase_cases=int(tab.timepoint_registry.eq(phase).sum()),cases_atleast20_concordant_cells=int((counts>=20).sum()),
                declared_author_lineage_labels=labels,intersection_not_annotation_update=True))
    background=pd.concat(selected,ignore_index=True).sort_values('position').reset_index(drop=True)
    assert background.position.is_unique and not background.major.eq('T').any()
    return background,coverage,dict(original_broad_crosscounts_replayed=True,public_source_SHA256=expected)


def main():
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert out.parent.parent==Path('/dev/shm') and not out.exists();out.mkdir(mode=0o700)
    def save(name,v):(out/name).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    try:
        root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1'
        previous=Path(os.environ['PANCANCER_BASELINE_OUTPUT']);original_input=Path(os.environ['PANCANCER_RNA_INPUT_OUTPUT'])
        for folder,expected in [(previous,os.environ['PANCANCER_BASELINE_TERMINAL_SHA256']),(original_input,os.environ['PANCANCER_RNA_INPUT_TERMINAL_SHA256'])]:
            assert sha(folder/'private_terminal.json')==expected
            for name,item in json.loads((folder/'private_terminal.json').read_text())['files'].items():assert sha(folder/name)==item['SHA256']
        old=json.loads((previous/'private_all32_RNA_programme_associations.json').read_text())['results']
        for name,digest in json.loads((previous/'plan_before_native_RNA_effects.json').read_text())['script_SHA256'].items():assert sha(Path(__file__).parent/name)==digest
        save('plan_before_shared_background_scores.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source='GSE286827',scope='2026-10-14_CM02_v1_T_whole_TME_primary_tumor',
            script_SHA256={p.name:sha(p) for p in Path(__file__).parent.glob('*.py')},original_terminal_SHA256=sha(previous/'private_terminal.json'),
            background_ontology=ONTOLOGY,minimum_concordant_background_cells_per_lineage=20,minimum_cases=20,minimum_residual_df=10,maximum_condition=1e8,
            slots='original eight G2M/G2M opportunities only',joint_background='equal arithmetic mean of B and Myeloid raw case mean G2M scores',
            controls='all original eight continuous controls, arm and HPV, plus one joint background; no automatic dropping',
            comparisons='refit original and expanded model on identical joint-background-eligible cases; B-only/Myeloid-only points on same cases',
            bootstrap=2000,minimum_valid_bootstrap=1800,seed=[2026102206,'original zero-based slot'],resampling='paired cases within arm x HPV strata',
            no_new_P_or_BH=True,no_clinical_tests=True,no_individual_or_learned_arrays_written=True,
            not_independent_validation=True,background_may_capture_real_shared_biology_not_only_artifact=True))
        fine=Path(os.environ['PANCANCER_FINE_LOADER']);ft=Path(os.environ['PANCANCER_FINE_SCOPE_OUTPUT'])/'private_terminal.json'
        loader=baseline.load(fine,json.loads(ft.read_text())['source_SHA256'],'accepted_fine');tab,cur,unique,shared,matrices,table=loader.load_inputs()
        net=baseline.load(os.environ['PANCANCER_NETWORK_RESOLVER'],os.environ['PANCANCER_NETWORK_RESOLVER_SHA256'],'accepted_network')
        contexts,membership=net.resolve(atlas,['all','reg','CD8'],['39697d077b6cb8bf','c3f761a9da61f578']);nodes=contexts['CM02_network_member_union']
        eligibility=baseline.load(os.environ['PANCANCER_ELIGIBILITY_SOURCE'],'a05c1d06b3c873b8abc0e78da20e420b20a0768e241487a43b17439fbfa2d3bd','accepted_eligibility')
        assert eligibility.reference_coverage(tab,unique,nodes)==json.loads((original_input/'private_author_parent_phase_coverage.json').read_text())
        assert membership==json.loads((original_input/'private_original_network_context_qualification.json').read_text())
        background,background_coverage,broad_audit=broad_intersection(tab,cur)
        save('private_background_lineage_coverage.json',dict(source='GSE286827',coverage=background_coverage,**broad_audit))
        unique=unique.loc[unique.major.eq('T')].sort_values('position').reset_index(drop=True);registry=tab.set_index('library_id')
        assert cur.library_id.map(registry.timepoint_registry).eq(cur.timepoint).all()
        assert not set(unique.position)&set(background.position)
        union=np.sort(np.r_[unique.position.to_numpy(int),background.position.to_numpy(int)]);assert len(np.unique(union))==len(union)
        z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r')
        selected=np.flatnonzero(np.asarray(read_elem(z['obs/dataset'])).astype(str)=='GSE286827')
        assert len(selected)==len(cur) and np.array_equal(cur.position,np.arange(len(cur)))
        var=read_elem(z['var']);order=list(map(str,read_elem(z['uns/source_order'])));measured=np.asarray(read_elem(z['varm/measured_in_source']),bool)[:,order.index('GSE286827')]
        human=measured&var.feature_identity_status.astype(str).eq('stable_Ensembl_ID').to_numpy();assert human.sum()>3000
        signature=Path(os.environ['PANCANCER_PROGRAMME_DEFINITION']);assert sha(signature)==os.environ['PANCANCER_PROGRAMME_SHA256']
        programmes=json.loads(signature.read_text())['programs'];columns=baseline.programme_columns(var,human,programmes)
        gene=dict(source='GSE286827',measured_stable_human_features=int(human.sum()),programmes=eligibility.feature_coverage(var,human,programmes))
        assert gene==json.loads((original_input/'private_source_gene_programme_qualification.json').read_text())
        scores,total,detected,score_audit=baseline.read_scores(z,selected[union],human,measured,columns,
            lambda v:save('private_progress.json',dict(observed_UTC=datetime.now(timezone.utc).isoformat(),**v)))
        tloc=np.searchsorted(union,unique.position);bloc=np.searchsorted(union,background.position)
        units,coverage=baseline.inputs.aggregate(tab,cur,unique,nodes,scores[tloc],total[tloc],detected[tloc])
        prev_audit=json.loads((previous/'private_RNA_input_and_score_audit.json').read_text());assert coverage==prev_audit['case_coverage']
        assert len(unique)==prev_audit['score_audit']['source_cells_scored']
        background['case']=background.library_id.map(registry.patient_registry);background['G2M']=scores[bloc,1]
        bg=background.groupby(['timepoint','case','major']).agg(n=('position','size'),score=('G2M','mean'))
        results=[];replayed=0
        for slot,definition in enumerate(baseline.statistics.fixed_slots()):
            x,y,c,cat,cases=baseline.inputs.pair_arrays(units,definition)
            fit=baseline.statistics.fit(x,y,c,cat,baseline.statistics.categorical_spec(cat))
            np.testing.assert_allclose([fit['rho'],fit['unadjusted'],fit['condition']],[old[slot]['rho'],old[slot]['unadjusted_rho'],old[slot]['design_condition']],atol=1e-11,rtol=1e-11);replayed+=1
            if definition['left_programme']!='g2m.genes' or definition['right_programme']!='g2m.genes':continue
            accepted=[];bs=[];ms=[];reason=Counter()
            for i,case in enumerate(cases):
                values=[];missing=[]
                for lineage in ['B','Myeloid']:
                    key=(definition['phase'],case,lineage)
                    if key not in bg.index or bg.loc[key,'n']<20:missing.append(lineage)
                    else:values.append(float(bg.loc[key,'score']))
                if missing:reason['AND'.join(missing)+'_CONCORDANT_REFERENCE_BELOW20_OR_UNCAPTURED']+=1
                else:accepted.append(i);bs.append(values[0]);ms.append(values[1]);reason['JOINT_ELIGIBLE']+=1
            ix=np.array(accepted,int)
            result=sensitivity.analyze(x[ix],y[ix],c[ix],cat[ix],bs,ms,slot)
            independent=sensitivity.analyze(x[ix],y[ix],c[ix],cat[ix],bs,ms,slot,independent=True);sensitivity.assert_replay(result,independent)
            results.append(dict(original_slot=slot+1,**definition,original_T_eligible_cases=len(cases),background_eligibility_reasons=dict(reason),
                independent_all_paired_bootstrap_replayed=True,**result))
            save('private_progress.json',dict(observed_UTC=datetime.now(timezone.utc).isoformat(),stage='paired_background_bootstrap_and_independent_replay',completed_slots=len(results)))
        assert replayed==32 and len(results)==8
        save('private_all8_shared_background_sensitivities.json',dict(source='GSE286827',results=results,no_new_P_or_BH=True))
        save('private_score_and_scope_audit.json',dict(source='GSE286827',score_audit=score_audit,
            original_T_subset_scored_cells=len(unique),background_scored_cells=len(background),original_all32_scalars_and_case_coverage_replayed=True,
            background_intersection_current_T_disjoint=True,no_individual_or_learned_arrays_written=True))
        assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024<=8*1024**3
        save('private_terminal.json',dict(status='COMPLETE_SAME_CASE_COMMON_IMMUNE_BACKGROUND_SENSITIVITY',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,no_new_P_or_BH=True,no_individual_or_learned_arrays_written=True,
            files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()}))
        print(json.dumps(dict(status='COMPLETE_SAME_CASE_COMMON_IMMUNE_BACKGROUND_SENSITIVITY',terminal_SHA256=sha(out/'private_terminal.json'))))
    except Exception as error:
        save('private_failed_terminal.json',dict(status='COMMON_IMMUNE_BACKGROUND_NOT_ACCEPTED',observed_UTC=datetime.now(timezone.utc).isoformat(),error_class=type(error).__name__,
            frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)],no_full_error_or_individual_values_written=True))
        raise SystemExit(1)


if __name__=='__main__':main()

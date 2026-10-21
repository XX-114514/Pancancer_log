"""Fixed multigene source associations with RAM-only individual measurements.

Only source-level coverage, fixed-slot effects, and hashes may be saved. Exact
public author-reference bytes enter through the preverified RAM cache driver;
no new individual expression/clinical/score tables are persisted or exported.
"""
from collections import defaultdict
from datetime import datetime,timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import traceback
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
import RNA_programme_case_inputs_v1 as inputs
import RNA_programme_partial_rank_v1 as statistics
import ucell_native_sparse_scores_v1 as scoring
import project_original_CM02_activity_in_RAM_v1 as projection
import audit_RNA_programme_partial_rank_v1 as audit

SOURCE='GSE286827'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path,expected,name):
    assert sha(path)==expected
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def programme_columns(var,human,programmes):
    lookup=defaultdict(set)
    for column in np.flatnonzero(human):
        row=var.iloc[column];names={str(var.index[column]),str(row.gene_symbol)}|set(str(row.source_gene_aliases).split('|'))
        for name in names:lookup[name].add(int(column))
    result=[]
    for name in ['s.genes','g2m.genes']:
        genes=programmes[name];mapped=[next(iter(lookup[g])) for g in genes if len(lookup[g])==1]
        counts={column:mapped.count(column) for column in mapped};mapped=[column for column in mapped if counts[column]==1]
        assert len(mapped)>=np.ceil(.8*len(genes)) and len(mapped)==len(set(mapped));result.append(np.array(mapped,np.int64))
    return result


def read_scores(z,selected,human,measured,programmes,checkpoint):
    raw=z['layers/counts'];shape=tuple(raw.attrs['shape']);ptr=raw['indptr'][:]
    assert len(ptr)==shape[0]+1 and ptr[0]==0 and np.all(np.diff(ptr)>=0)
    assert ptr[-1]==raw['data'].shape[0]==raw['indices'].shape[0] and np.all(np.diff(selected)>0)
    scores=np.zeros((len(selected),2));total=np.zeros(len(selected),np.int64);detected=total.copy();replayed=0;completed=0
    for start in range(0,shape[0],5000):
        stop=min(start+5000,shape[0]);lo,hi=np.searchsorted(selected,[start,stop])
        if lo==hi:continue
        a,b=int(ptr[start]),int(ptr[stop]);indices=raw['indices'][a:b];values=raw['data'][a:b];localptr=ptr[start:stop+1]-a
        rows=selected[lo:hi]-start
        ss,tt,dd=scoring.score_selected_CSR_rows(localptr,indices,values,rows,human,measured,programmes)
        # Independent dense pandas ranking on both boundary rows of each block.
        for position in sorted({0,len(rows)-1}):
            row=int(rows[position]);begin,end=map(int,localptr[row:row+2]);dense=np.zeros(shape[1],np.int64)
            np.add.at(dense,np.asarray(indices[begin:end],int),np.asarray(values[begin:end],np.int64))
            ranks=np.full(shape[1],np.nan);ranks[human]=pd.Series(-dense[human]).rank(method='average').to_numpy()
            expected=[]
            for columns in programmes:
                n=len(columns);minimum=n*(n+1)/2;maximum=min(1500,int(human.sum()))
                expected.append(1-(np.minimum(ranks[columns],maximum).sum()-minimum)/(n*maximum-minimum))
            np.testing.assert_allclose(ss[position],expected,atol=1e-12,rtol=1e-12)
            assert dense[human].sum()==tt[position] and np.count_nonzero(dense[human])==dd[position];replayed+=1
        scores[lo:hi]=ss;total[lo:hi]=tt;detected[lo:hi]=dd;completed+=hi-lo
        checkpoint(dict(stage='native_CSR_scores',completed_source_cells=int(completed),independent_dense_score_rows=replayed))
        assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024<=8*1024**3
    assert completed==len(selected)
    return scores,total,detected,dict(source_cells_scored=int(completed),independent_dense_boundary_rows=replayed,
        all_selected_rows_canonicalized_before_rank=True,no_rank_or_individual_array_written=True)


def independent_base_replay(x,y,c,cat,result):
    if result['status']!='ESTIMABLE':return False
    # Build ranks/contrasts independently of the production rank/QR helpers.
    columns=[np.ones(len(x))]
    for j in range(c.shape[1]):
        r=pd.Series(c[:,j]).rank(method='average').to_numpy();columns.append((r-r.mean())/r.std())
    for j in range(cat.shape[1]):
        for label in sorted(set(cat[:,j]))[1:]:columns.append((cat[:,j]==label).astype(float))
    z=np.column_stack(columns);rx=pd.Series(x).rank().to_numpy();ry=pd.Series(y).rank().to_numpy()
    ex=rx-z@np.linalg.lstsq(z,rx,rcond=None)[0];ey=ry-z@np.linalg.lstsq(z,ry,rcond=None)[0]
    np.testing.assert_allclose(result['rho'],np.corrcoef(ex,ey)[0,1],atol=1e-11,rtol=1e-11)
    np.testing.assert_allclose(result['unadjusted_rho'],np.corrcoef(rx,ry)[0,1],atol=1e-12)
    return True


def main():
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert out.parent.parent==Path('/dev/shm') and not out.exists();out.mkdir(mode=0o700)
    def save(name,value):(out/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    try:
        root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1'
        accepted_input=Path(os.environ['PANCANCER_RNA_INPUT_OUTPUT']);terminal=accepted_input/'private_terminal.json'
        assert sha(terminal)==os.environ['PANCANCER_RNA_INPUT_TERMINAL_SHA256']
        old=json.loads(terminal.read_text());assert old['status']=='COMPLETE_OUTCOME_FREE_RNA_PROGRAMME_INPUT_QUALIFICATION'
        for name,v in old['files'].items():assert sha(accepted_input/name)==v['SHA256']
        eligibility=load(os.environ['PANCANCER_ELIGIBILITY_SOURCE'],old['source_SHA256'],'accepted_RNA_eligibility')
        signature=Path(os.environ['PANCANCER_PROGRAMME_DEFINITION']);assert sha(signature)==os.environ['PANCANCER_PROGRAMME_SHA256']
        definition=json.loads(signature.read_text());assert definition['excluded_prior_marker_genes']==['MKI67','TOP2A','UBE2C']
        fine=Path(os.environ['PANCANCER_FINE_LOADER']);ft=Path(os.environ['PANCANCER_FINE_SCOPE_OUTPUT'])/'private_terminal.json'
        loader=load(fine,json.loads(ft.read_text())['source_SHA256'],'accepted_fine_reference')
        network=load(os.environ['PANCANCER_NETWORK_RESOLVER'],os.environ['PANCANCER_NETWORK_RESOLVER_SHA256'],'accepted_CM02_network')
        scripts=Path(__file__).parent
        save('plan_before_native_RNA_effects.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),
            source=SOURCE,scope='2026-10-14_CM02_v1_T_whole_TME_primary_tumor',slots=statistics.fixed_slots(),
            script_SHA256={p.name:sha(p) for p in scripts.glob('*.py')},input_terminal_SHA256=sha(terminal),
            programme_SHA256=sha(signature),prior_markers_excluded=True,minimum_parent=20,minimum_cases=20,
            maxRank=1500,minimum_measured_coverage=.8,bootstrap=2000,minimum_valid_bootstrap=1800,
            permutations=2000,seed=2026102204,min_residual_df=10,maximum_design_condition=1e8,
            controls=['left_fraction_allT','right_fraction_allT','left_log1p_cells','right_log1p_cells','left_mean_log1p_UMI','right_mean_log1p_UMI','left_mean_log1p_detected','right_mean_log1p_detected'],
            categorical_controls=['arm','HPV'],no_continuous_control_dropped=True,
            residual_permutation='rankY Freedman-Lane within arm x HPV; pseudoY not reranked; conditional approximation',
            bootstrap_design='paired reviewed cases within arm x HPV, rerank/refit',BH_family=32,
            unestimable_P1_in_family_not_negative=True,post2_coverage_only=True,
            independent_actual_score_replay='first/last selected raw row per CSR block, dense pandas ranks',
            independent_observed_statistics='all slots and all random replicates; independent pandas ranks and ordinary least squares',
            activity_correspondence='all 3 parents x2 programmes x2 phases x2 contexts, rho only on >=20 cases, no P/CI or clinical family',
            no_individual_or_learned_arrays_written=True,no_outcome_or_sequence_selected=True))
        tab,cur,unique,shared,matrices,table=loader.load_inputs()
        contexts,membership=network.resolve(atlas,['all','regulatory_anchor','cytotoxic_anchor'],['39697d077b6cb8bf','c3f761a9da61f578']);nodes=contexts['CM02_network_member_union']
        coverage=eligibility.reference_coverage(tab,unique,nodes)
        assert coverage==json.loads((accepted_input/'private_author_parent_phase_coverage.json').read_text())
        assert membership==json.loads((accepted_input/'private_original_network_context_qualification.json').read_text())
        z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r')
        selected=np.flatnonzero(np.asarray(read_elem(z['obs/dataset'])).astype(str)==SOURCE)
        assert len(selected)==len(cur) and np.array_equal(cur.position,np.arange(len(cur)))
        var=read_elem(z['var']);source_order=list(map(str,read_elem(z['uns/source_order'])))
        measured=np.asarray(read_elem(z['varm/measured_in_source']),bool)[:,source_order.index(SOURCE)]
        human=measured&var.feature_identity_status.astype(str).eq('stable_Ensembl_ID').to_numpy()
        gene=dict(source=SOURCE,measured_stable_human_features=int(human.sum()),programmes=eligibility.feature_coverage(var,human,definition['programs']))
        assert gene==json.loads((accepted_input/'private_source_gene_programme_qualification.json').read_text())
        assert all(r['source_coverage_accepted'] for r in gene['programmes'])
        columns=programme_columns(var,human,definition['programs'])
        for name,col in zip(['s.genes','g2m.genes'],columns):assert len(col)==next(r for r in gene['programmes'] if r['program']==name)['uniquely_mapped_measured_genes']
        unique=unique.loc[unique.major.eq('T')].sort_values('position').reset_index(drop=True)
        scores,total,detected,score_receipt=read_scores(z,selected[unique.position.to_numpy(int)],human,measured,columns,lambda v:save('private_progress.json',dict(observed_UTC=datetime.now(timezone.utc).isoformat(),**v)))
        units,coverage2=inputs.aggregate(tab,cur,unique,nodes,scores,total,detected)
        for r in coverage2:
            previous=next(p for p in coverage['pair_coverage'] if all(r[k]==p[k] for k in ['context','phase','left_parent','right_parent']))
            assert r['eligible_cases']==previous['cases_with_concordant_known_arm_and_HPV']
        save('private_RNA_input_and_score_audit.json',dict(gene_coverage=gene,case_coverage=coverage2,score_audit=score_receipt,
            original_coverage_and_network_exactly_replayed=True,new_individual_arrays_RAM_only=True))
        results=[];audited=0
        for slot,definition_slot in enumerate(statistics.fixed_slots()):
            x,y,c,cat,cases=inputs.pair_arrays(units,definition_slot)
            result=statistics.analyze(x,y,c,cat,slot)
            audited+=independent_base_replay(x,y,c,cat,result)
            results.append(dict(slot=slot+1,**definition_slot,**result))
            save('private_progress.json',dict(stage='fixed32_associations',completed_slots=slot+1,observed_UTC=datetime.now(timezone.utc).isoformat()))
        assert len(results)==32;statistics.bh_fixed(results)
        save('private_all32_RNA_programme_associations.json',dict(source=SOURCE,results=results,independent_base_replayed_slots=audited,
            pre_post_not_independent_replication=True,reference_subset_not_lineage_truth=True,
            intervals_missing_not_negative=True,residual_exchange_approximation=True))
        audit_rows=[]
        for slot,definition_slot in enumerate(statistics.fixed_slots()):
            x,y,c,cat,cases=inputs.pair_arrays(units,definition_slot)
            audit_rows.append(dict(slot=slot+1,**audit.replay(x,y,c,cat,slot,results[slot])))
            save('private_progress.json',dict(stage='independent_all32_random_replay',completed_slots=slot+1,observed_UTC=datetime.now(timezone.utc).isoformat()))
        save('private_all32_independent_inference_replay.json',dict(source=SOURCE,slots=audit_rows,
            numerical_replay_not_new_biological_replication=True))
        activity,projection_receipt=projection.project(root,atlas,tab,Path(os.environ['PANCANCER_ORIGINAL_PHASE_SOURCE']),os.environ['PANCANCER_ORIGINAL_PHASE_SOURCE_SHA256'])
        save('private_original_activity_projection_parity.json',projection_receipt)
        descriptions=[]
        for context in inputs.CONTEXTS:
            for phase in ['pre','post']:
                for role in inputs.PARENTS:
                    records=[(case,r) for (co,ph,case,parent),r in sorted(units.items()) if co==context and ph==phase and parent==role and r['n']>=20 and r['arm'] in ['D','D+T'] and r['HPV'] in ['positive','negative']]
                    for programme,field in [('s.genes','S'),('g2m.genes','G2M')]:
                        a=np.array([activity[case,phase] for case,r in records]);s=np.array([r[field] for case,r in records])
                        value=statistics.correlation(statistics.average_rank(a),statistics.average_rank(s)) if len(records)>=20 else None
                        descriptions.append(dict(context=context,phase=phase,parent=role,programme=programme,n_cases=len(records),
                            rho=value,status='DESCRIPTIVE_SAME_SOURCE_ONLY' if value is not None else 'NOT_ESTIMABLE',P_not_computed=True))
        assert len(descriptions)==24
        save('private_all24_programme_activity_descriptions.json',dict(source=SOURCE,descriptions=descriptions,
            circular_same_source_not_independent_validation=True,no_new_clinical_model_or_fitted_weights=True))
        assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024<=8*1024**3
        save('private_terminal.json',dict(status='COMPLETE_FIXED32_NATIVE_RNA_PROGRAMME_ASSOCIATIONS',
            completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),
            all_fixed_slots_retained=True,no_individual_or_learned_arrays_written=True,
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
            files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()}))
        print(json.dumps(dict(status='COMPLETE_FIXED32_NATIVE_RNA_PROGRAMME_ASSOCIATIONS',terminal_SHA256=sha(out/'private_terminal.json'))))
    except Exception as error:
        save('private_failed_terminal.json',dict(status='FIXED32_RNA_PROGRAMME_ANALYSIS_NOT_ACCEPTED',
            observed_UTC=datetime.now(timezone.utc).isoformat(),error_class=type(error).__name__,
            frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)],
            individual_values_and_full_errors_not_returned=True))
        raise SystemExit(1)


if __name__=='__main__':main()

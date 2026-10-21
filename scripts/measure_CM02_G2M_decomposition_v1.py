"""Finite raw-score diagnostics with individual arrays only in RAM.

Preserves all closed rank tests. Returns additive source covariance terms,
measurement concentration summaries, author/node counts and coarse non-T
coverage only. No new P, CI, clinical test or learned/individual export.
"""
from datetime import datetime,timezone
import hashlib,json,os,resource,traceback
from pathlib import Path
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
import analyze_CM02_RNA_multigene_in_RAM_v1 as baseline
import CM02_G2M_additive_diagnostics_v1 as diagnostic
import CM02_author_node_bridge_v1 as bridge


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_extended(z,selected,human,measured,programmes,checkpoint):
    raw=z['layers/counts'];shape=tuple(raw.attrs['shape']);ptr=raw['indptr'][:]
    assert len(ptr)==shape[0]+1 and ptr[0]==0 and np.all(np.diff(ptr)>=0) and ptr[-1]==raw['data'].shape[0]==raw['indices'].shape[0]
    scores=np.zeros((len(selected),2));total=np.zeros(len(selected),np.int64);detected=total.copy();diag=np.zeros((len(selected),2,4));completed=0;replayed=0
    for start in range(0,shape[0],5000):
        stop=min(start+5000,shape[0]);lo,hi=np.searchsorted(selected,[start,stop])
        if lo==hi:continue
        a,b=int(ptr[start]),int(ptr[stop]);index=raw['indices'][a:b];value=raw['data'][a:b];pp=ptr[start:stop+1]-a;rows=selected[lo:hi]-start
        ss,tt,dd=baseline.scoring.score_selected_CSR_rows(pp,index,value,rows,human,measured,programmes)
        for j,row in enumerate(rows):
            begin,end=map(int,pp[row:row+2]);idx=np.asarray(index[begin:end],int);val=np.asarray(value[begin:end],np.int64)
            take=human[idx]&(val>0);unique,inverse=np.unique(idx[take],return_inverse=True);counts=np.zeros(len(unique),np.int64);np.add.at(counts,inverse,val[take]);query=[]
            for columns in programmes:
                found=np.searchsorted(unique,columns);ok=found<len(unique);good=np.zeros(len(columns),bool);good[ok]=unique[found[ok]]==columns[ok]
                values=np.zeros(len(columns));values[good]=counts[found[good]];query.append(values)
            check,current=diagnostic.signature_diagnostics(counts,query,int(human.sum()));np.testing.assert_allclose(check,ss[j],atol=1e-12,rtol=1e-12);diag[lo+j]=current
            if j in [0,len(rows)-1]:
                dense=np.zeros(shape[1],np.int64);np.add.at(dense,idx,val);rank=np.full(shape[1],np.nan);rank[human]=pd.Series(-dense[human]).rank(method='average').to_numpy()
                for k,columns in enumerate(programmes):
                    deficit=1500-np.minimum(rank[columns],1500);mass=deficit.sum()
                    expected=[np.count_nonzero(dense[columns]),np.count_nonzero(deficit),mass**2/(deficit@deficit) if mass else 0,deficit.max()/mass if mass else 0]
                    np.testing.assert_allclose(current[k],expected,atol=1e-12,rtol=1e-12)
                assert dense[human].sum()==tt[j] and np.count_nonzero(dense[human])==dd[j];replayed+=1
        scores[lo:hi]=ss;total[lo:hi]=tt;detected[lo:hi]=dd;completed+=hi-lo
        checkpoint(dict(stage='native_signature_support',completed_source_cells=int(completed),independent_dense_diagnostic_rows=replayed))
        assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024<=8*1024**3
    assert completed==len(selected)
    return scores,total,detected,diag,dict(source_cells_scored=int(completed),independent_dense_boundary_rows=replayed,
        all_selected_rows_canonicalized_before_rank=True,no_rank_or_individual_array_written=True)


def main():
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert out.parent.parent==Path('/dev/shm') and not out.exists();out.mkdir(mode=0o700)
    def save(name,v):(out/name).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    try:
        root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1'
        previous=Path(os.environ['PANCANCER_BASELINE_OUTPUT']);stability=Path(os.environ['PANCANCER_STABILITY_OUTPUT'])
        for folder,expected in [(previous,os.environ['PANCANCER_BASELINE_TERMINAL_SHA256']),(stability,os.environ['PANCANCER_STABILITY_TERMINAL_SHA256'])]:
            assert sha(folder/'private_terminal.json')==expected
            for name,v in json.loads((folder/'private_terminal.json').read_text())['files'].items():assert sha(folder/name)==v['SHA256']
        old=json.loads((previous/'private_all32_RNA_programme_associations.json').read_text())['results']
        for name,v in json.loads((previous/'plan_before_native_RNA_effects.json').read_text())['script_SHA256'].items():assert sha(Path(__file__).parent/name)==v
        save('plan_before_raw_score_diagnostics.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source='GSE286827',scope='2026-10-14_CM02_v1_T_whole_TME_primary_tumor',
            code_SHA256={p.name:sha(p) for p in Path(__file__).parent.glob('*.py')},original_terminal_SHA256=sha(previous/'private_terminal.json'),stability_terminal_SHA256=sha(stability/'private_terminal.json'),
            decomposition='all8 fixed G2M/G2M case sets; raw means residualized on original ranked controls; two fixed reference mean schemes; all168 author fine-pair covariance terms',
            reference_means=['pooled_cell_weighted_reference','equal_observed_case_reference'],absent_observed_term_not_biological_absence=True,
            multigene_support='all3 parents and13 fine labels xpre/post x2 contexts; >=20 cells per case; >=5 cases to report distribution',
            metrics=['detected_genes','rank_contributing_genes','effective_gene_number','maximum_mass_fraction_positive_cells','fractions_active_genes_ge1_ge2_ge5_ge10'],
            coarse_nonT='all observed current nonT majors, coverage only; no nonT expression scores or effects',
            no_new_P_CI_or_clinical_tests=True,not_a_decomposition_of_primary_Spearman=True,no_individual_or_learned_arrays_written=True))
        fine=Path(os.environ['PANCANCER_FINE_LOADER']);ft=Path(os.environ['PANCANCER_FINE_SCOPE_OUTPUT'])/'private_terminal.json';loader=baseline.load(fine,json.loads(ft.read_text())['source_SHA256'],'accepted_fine')
        tab,cur,unique,shared,matrices,table=loader.load_inputs()
        net=baseline.load(os.environ['PANCANCER_NETWORK_RESOLVER'],os.environ['PANCANCER_NETWORK_RESOLVER_SHA256'],'accepted_network');contexts,membership=net.resolve(atlas,['all','reg','CD8'],['39697d077b6cb8bf','c3f761a9da61f578']);nodes=contexts['CM02_network_member_union']
        original_input=Path(os.environ['PANCANCER_RNA_INPUT_OUTPUT']);assert sha(original_input/'private_terminal.json')==os.environ['PANCANCER_RNA_INPUT_TERMINAL_SHA256']
        for name,v in json.loads((original_input/'private_terminal.json').read_text())['files'].items():assert sha(original_input/name)==v['SHA256']
        eligibility=baseline.load(os.environ['PANCANCER_ELIGIBILITY_SOURCE'],'a05c1d06b3c873b8abc0e78da20e420b20a0768e241487a43b17439fbfa2d3bd','accepted_eligibility')
        assert eligibility.reference_coverage(tab,unique,nodes)==json.loads((original_input/'private_author_parent_phase_coverage.json').read_text())
        assert membership==json.loads((original_input/'private_original_network_context_qualification.json').read_text())
        inputs=baseline.inputs;registry=tab.set_index('library_id');unique=unique.loc[unique.major.eq('T')].sort_values('position').reset_index(drop=True)
        assert cur.library_id.map(registry.timepoint_registry).eq(cur.timepoint).all()
        old_fine=json.loads((stability/'private_all52_author_fine_profiles.json').read_text())['profiles']
        cross=[dict(phase=p,current_major=g,current_node=n,author_label=l,matched_cells=int(v)) for (p,g,n,l),v in unique.groupby(['timepoint','major','node_id','label']).size().items()]
        node_bridge=bridge.bridge(cross,membership,old_fine)
        z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r');selected=np.flatnonzero(np.asarray(read_elem(z['obs/dataset'])).astype(str)=='GSE286827')
        assert len(selected)==len(cur) and np.array_equal(cur.position,np.arange(len(cur)))
        lookups=[]
        for field in ['subCluster','subcluster']:
            values=np.asarray(read_elem(z['obs/'+field])).astype(str)[selected];lookup={}
            for node,g in cur.loc[cur.major.eq('T')].groupby('node_id'):
                names=set(values[g.position]);assert len(names)==1;lookup[node]=next(iter(names))
            lookups.append(lookup)
        assert lookups[0]==lookups[1]
        for row in node_bridge:
            for v in row['node_counts']:v['current_frozen_label']=lookups[0][v['current_node']]
        save('private_all26_author_node_bridge.json',dict(source='GSE286827',rows=node_bridge,authority_two_label_fields_agree=True))
        nonT=[]
        current=cur.copy();current['case']=current.library_id.map(registry.patient_registry)
        for phase in ['pre','post']:
            phase_tab=tab.loc[tab.timepoint_registry.eq(phase)]
            for major in sorted(set(cur.major)-{'T'}):
                g=current.loc[current.major.eq(major)&current.timepoint.eq(phase)];counts=g.groupby('case').size();qual=counts[counts>=20].index
                nonT.append(dict(phase=phase,current_major=major,reviewed_phase_cases=len(phase_tab),cases_observed=len(counts),cases_atleast20_current_cells=len(qual),
                    current_label_not_independent_lineage_validation=True,expression_not_scored=True))
        save('private_coarse_nonT_coverage_only.json',dict(source='GSE286827',coverage=nonT,no_new_nonT_effects=True))
        var=read_elem(z['var']);order=list(map(str,read_elem(z['uns/source_order'])));measured=np.asarray(read_elem(z['varm/measured_in_source']),bool)[:,order.index('GSE286827')];human=measured&var.feature_identity_status.astype(str).eq('stable_Ensembl_ID').to_numpy();assert human.sum()>3000
        signature=Path(os.environ['PANCANCER_PROGRAMME_DEFINITION']);assert sha(signature)==os.environ['PANCANCER_PROGRAMME_SHA256'];programmes=json.loads(signature.read_text())['programs'];columns=baseline.programme_columns(var,human,programmes)
        gene=dict(source='GSE286827',measured_stable_human_features=int(human.sum()),programmes=eligibility.feature_coverage(var,human,programmes))
        assert gene==json.loads((original_input/'private_source_gene_programme_qualification.json').read_text())
        scores,total,detected,diag,score_audit=read_extended(z,selected[unique.position.to_numpy(int)],human,measured,columns,lambda v:save('private_progress.json',dict(observed_UTC=datetime.now(timezone.utc).isoformat(),**v)))
        prev_audit=json.loads((previous/'private_RNA_input_and_score_audit.json').read_text());assert score_audit==prev_audit['score_audit']
        units,coverage=inputs.aggregate(tab,cur,unique,nodes,scores,total,detected);assert coverage==prev_audit['case_coverage']
        f=unique.copy();f['case']=f.library_id.map(registry.patient_registry);f['G2M']=scores[:,1]
        all_decompositions=[];all_fine=[];replayed=0
        for slot,definition in enumerate(baseline.statistics.fixed_slots()):
            x,y,c,cat,cases=inputs.pair_arrays(units,definition);fit=baseline.statistics.fit(x,y,c,cat,baseline.statistics.categorical_spec(cat))
            np.testing.assert_allclose([fit['rho'],fit['unadjusted'],fit['condition']],[old[slot]['rho'],old[slot]['unadjusted_rho'],old[slot]['design_condition']],atol=1e-11,rtol=1e-11);replayed+=1
            if definition['left_programme']!=definition['right_programme'] or definition['left_programme']!='g2m.genes':continue
            current=f.loc[f.timepoint.eq(definition['phase'])]
            if definition['context']==inputs.CONTEXTS[1]:current=current.loc[current.node_id.isin(nodes)]
            grouped=current.groupby(['case','label']).agg(n=('position','size'),score_sum=('G2M','sum'))
            arrays=[];labelsets=[inputs.REGULATORY if definition['left_parent']==inputs.PARENTS[0] else inputs.CONVENTIONAL,inputs.CD8]
            for labels in labelsets:
                counts=np.zeros((len(cases),len(labels)),np.int64);sums=np.zeros_like(counts,dtype=float)
                for i,case in enumerate(cases):
                    for j,label in enumerate(labels):
                        if (case,label) in grouped.index:counts[i,j]=int(grouped.loc[(case,label),'n']);sums[i,j]=grouped.loc[(case,label),'score_sum']
                arrays.append((counts,sums))
            np.testing.assert_allclose(arrays[0][1].sum(axis=1)/arrays[0][0].sum(axis=1),x,atol=1e-13)
            np.testing.assert_allclose(arrays[1][1].sum(axis=1)/arrays[1][0].sum(axis=1),y,atol=1e-13)
            result,fine_pairs=diagnostic.additive_decomposition(arrays[0][0],arrays[1][0],arrays[0][1],arrays[1][1],fit['z'],*labelsets)
            all_decompositions.append(dict(original_slot=slot+1,**definition,original_partial_Spearman=old[slot]['rho'],**result))
            all_fine.extend(dict(original_slot=slot+1,**definition,**v) for v in fine_pairs)
        assert replayed==32 and len(all_decompositions)==8 and len(all_fine)==168
        save('private_all8_raw_score_decompositions.json',dict(source='GSE286827',rows=all_decompositions,no_new_primary_rank_tests=True))
        save('private_all168_fine_pair_covariance_terms.json',dict(source='GSE286827',rows=all_fine,terms_not_independently_validated_communications=True))
        profiles=diagnostic.support_profiles(tab,unique,nodes,diag,inputs)
        for r in profiles:
            if r['reference_kind']=='fine':
                previous_row=next(p for p in old_fine if p['context']==r['context'] and p['phase']==r['phase'] and p['author_fine_label']==r['reference_label'])
                assert r['cases_atleast20_observed_cells']==previous_row['cases_atleast20_fine_cells'] and r['source_cell_observations']==previous_row['unique_reference_cell_observations']
        save('private_all64_multigene_support_profiles.json',dict(source='GSE286827',profiles=profiles,positive_rank_genes_not_cell_cycle_identity=True))
        save('private_score_and_scope_audit.json',dict(source='GSE286827',score_audit=score_audit,all32_observed_scalars_and_52_fine_counts_replayed=True,no_new_individual_or_learned_arrays_written=True))
        assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024<=8*1024**3
        save('private_terminal.json',dict(status='COMPLETE_RAW_SCORE_COMPOSITION_AND_MULTIGENE_DIAGNOSTICS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,no_new_P_CI_or_clinical_tests=True,no_individual_or_learned_arrays_written=True,
            files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()}))
        print(json.dumps(dict(status='COMPLETE_RAW_SCORE_COMPOSITION_AND_MULTIGENE_DIAGNOSTICS',terminal_SHA256=sha(out/'private_terminal.json'))))
    except Exception as error:
        save('private_failed_terminal.json',dict(status='RAW_SCORE_DIAGNOSTICS_NOT_ACCEPTED',observed_UTC=datetime.now(timezone.utc).isoformat(),error_class=type(error).__name__,
            frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)],no_full_error_or_individual_values_written=True))
        raise SystemExit(1)


if __name__=='__main__':main()

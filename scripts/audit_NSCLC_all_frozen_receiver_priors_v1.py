"""Independently verify all receiver backgrounds, grouped-tie AP and shortlists.

Uses manual precision-recall sums/AUC ranks rather than sklearn metrics.
This audits prior/input ranking, not native single-cell receptor complexes.
"""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd,pyreadr
from scipy.stats import rankdata

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def same(a,b):assert np.allclose(a,b,rtol=1e-9,atol=1e-10,equal_nan=True),(a,b)

def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=R/'NSCLC_all_frozen_receiver_prior_ranking_20261016_v1'
    launch=json.loads((R/'HNSCC_receiver_followup_job_control_20261016_v2/rank_NSCLC_all_frozen_receiver_priors_v1/launch_receipt.json').read_text());deadline=time.monotonic()+7200
    while not (O/'private_terminal.json').exists():
        stat=Path('/proc')/str(launch['pid'])/'stat';assert stat.exists()
        s=stat.read_text().rsplit(')',1)[1].split();assert s[19]==str(launch['proc_start_ticks']) and s[0]!='Z'
        assert time.monotonic()<deadline;time.sleep(10)
    terminal=json.loads((O/'private_terminal.json').read_text());output=O/'private_independent_audit.json';assert not output.exists()
    plan=json.loads((O/'plan_before_any_ranking_or_sender_effect.json').read_text())
    for f,h in terminal['files'].items():assert sha(O/f)==h
    prior=R/'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1/verified_prior_inputs_v2'
    for f,h in plan['prior_sha256'].items():assert sha(prior/f)==h
    T=next(iter(pyreadr.read_r(str(prior/'ligand_target_matrix_nsga2r_final.rds')).values()));T.index=T.index.astype(str);T.columns=T.columns.astype(str)
    LR=next(iter(pyreadr.read_r(str(prior/'lr_network_human_21122021.rds')).values()))
    N=R/'NSCLC_full_native_gene_units_20261016_v1';units=pd.read_csv(N/'private_all_native_units_and_quality.tsv',sep='\t')
    native=pd.Index(pd.read_csv(N/'private_native_gene_axis.tsv',sep='\t').gene)
    hits=np.load(N/'positive_cell_counts.npy',mmap_mode='r');detection=pd.read_csv(O/'private_all_lineage_native_detection.tsv.gz',sep='\t')
    expressed={}
    for family,prefix in {'Macrophage':'M\u03c6_','CD8':'CD8T_','Treg':'CD4T_Treg_'}.items():
        us=units[units.sub_cell_type.str.startswith(prefix)];counts=np.zeros(len(native),np.int64)
        for i in us.unit_index:counts+=hits[int(i)]
        f=detection[detection.source_family.eq(family)].set_index('gene').loc[native]
        assert f.native_cells.eq(int(us.n_cells.sum())).all() and np.array_equal(f.positive_cells,counts)
        gate=counts*10>=int(us.n_cells.sum());assert np.array_equal(gate,f.expressed10pct)
        expressed[family]=set(native[gate])
    ppath=R/'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv';assert sha(ppath)==terminal['parameters_sha256']
    params=pd.read_csv(ppath,sep='\t');scope=pd.read_csv(R/'NSCLC_context_native_programme_projection_20261016_v1/private_all56_projection_scope.tsv',sep='\t');scope=scope[scope.status.str.startswith('PARTIAL_FROZEN')]
    V=Path(json.loads((R/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    backgrounds=pd.read_csv(O/'private_all24_native_prior_backgrounds.tsv.gz',sep='\t')
    coverage=pd.read_csv(O/'private_all48_signed_programme_coverage.tsv',sep='\t')
    ranks=pd.read_csv(O/'private_all_ligand_rankings.tsv.gz',sep='\t');slots=pd.read_csv(O/'private_all960_prospective_prior_slots.tsv',sep='\t',keep_default_na=False)
    assert len(scope)==24 and len(coverage)==48 and len(ranks)==58848 and len(slots)==960
    assert not ranks.duplicated(['signature','programme_sign','ligand']).any()
    assert not slots.duplicated(['signature','programme_sign','prospective_position']).any()
    assert set(backgrounds.signature)==set(ranks.signature)==set(slots.signature)
    checked=0
    for row in scope.itertuples():
        signature=f'{row.context}_{row.type}_{row.program}';family=row.source_family
        discovery=set((V/'state_after_DIALOGUE_inputs_20261016_v3'/f'context_{row.context:04d}'/(row.type+'_discovery_genes.tsv')).read_text().splitlines())
        bg=sorted(discovery&expressed[family]&set(T.index))
        f=backgrounds[backgrounds.signature.eq(signature)].set_index('gene');assert set(f.index)==set(bg) and f.index.is_unique
        pp=params[params.context.eq(row.context)&params.type.eq(row.type)&params.program.eq(row.program)];weights=dict(zip(pp.gene,pp.weight))
        assert f.loc[bg,'programme_sign'].tolist()==['positive' if weights.get(g,0)>0 else 'negative' if weights.get(g,0)<0 else 'background' for g in bg]
        B=T.loc[bg]
        for sign in ['positive','negative']:
            y=np.array([(weights.get(g,0)>0 if sign=='positive' else weights.get(g,0)<0) for g in bg],int);n=int(y.sum());ok=n>=10 and n<len(bg)
            cv=coverage[coverage.signature.eq(signature)&coverage.programme_sign.eq(sign)];assert len(cv)==1
            assert int(cv.expressed_prior_genes.iloc[0])==n and int(cv.background_genes.iloc[0])==len(bg) and bool(cv.ranking_available.iloc[0])==ok
            table=ranks[ranks.signature.eq(signature)&ranks.programme_sign.eq(sign)].set_index('ligand');assert set(table.index)==set(T.columns)
            for ligand in T.columns:
                r=table.loc[ligand];x=B[ligand].to_numpy(float)
                receptors=sorted({str(r) for r in LR.loc[LR['from'].astype(str).eq(ligand),'to']}&expressed[family])
                actual='' if pd.isna(r.receiver_flat_receptors) else str(r.receiver_flat_receptors)
                assert actual=='|'.join(receptors) and bool(r.receiver_compatible)==bool(receptors)
                if ok:
                    order=np.argsort(-x,kind='stable');xx=x[order];yy=y[order];end=np.r_[np.flatnonzero(np.diff(xx)!=0),len(xx)-1]
                    tp=np.cumsum(yy)[end];delta=np.diff(np.r_[0,tp]);ap=float(np.sum(delta/n*tp/(end+1)))
                    auc=float((rankdata(x)[y.astype(bool)].sum()-n*(n+1)/2)/(n*(len(y)-n)))
                    same(ap,r.average_precision);same(ap-n/len(bg),r.AP_excess_over_prevalence);same(auc,r.AUROC)
                else:assert pd.isna(r.average_precision) and pd.isna(r.AP_excess_over_prevalence) and pd.isna(r.AUROC)
                checked+=1
            expected=table[table.receiver_compatible&table.AP_excess_over_prevalence.gt(0)].reset_index().sort_values(['AP_excess_over_prevalence','ligand'],ascending=[False,True]).head(20)
            ss=slots[slots.signature.eq(signature)&slots.programme_sign.eq(sign)].sort_values('prospective_position')
            assert ss.prospective_position.tolist()==list(range(1,21))
            assert ss.ligand.tolist()==expected.ligand.tolist()+['']*(20-len(expected))
            assert set(table[table.prospective_top20].index)==set(expected.ligand)
        (O/'private_audit_progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),signature=signature,all_prior_rankings_checked=checked)))
    refpath=R/'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv';assert sha(refpath)==terminal['curated_reference_sha256']
    ref=pd.read_csv(refpath,sep='\t',keep_default_na=False);edges=pd.read_csv(O/'private_all_shortlist_reference_interfaces.tsv',sep='\t')
    assert not edges.duplicated(['signature','ligand','receptor']).any()
    acceptable={'simplex_reference_supported','all_constituents_supported_complex_assembly_unverified'};expected_edges=0
    for slot in slots[slots.ligand.ne('')][['signature','source_family','ligand']].drop_duplicates().itertuples():
        er=ref[ref.ligand.eq(slot.ligand)];ee=edges[edges.signature.eq(slot.signature)&edges.ligand.eq(slot.ligand)].set_index('receptor')
        assert set(ee.index)==set(er.receptor);expected_edges+=len(er)
        for edge in er.itertuples():
            rr=ee.loc[edge.receptor];parts=str(edge.receptor).split('_');measured=all(g in native for g in parts);marginal=measured and all(g in expressed[slot.source_family] for g in parts)
            assert rr.reference_status==edge.reference_status and bool(rr.receptor_components_measured)==measured and bool(rr.all_component_pooled_marginals10pct)==marginal
            assert bool(rr.reference_and_marginal_compatible)==(edge.reference_status in acceptable and marginal)
            assert not rr.same_cell_RNA_verified and not rr.complex_assembly_verified
    assert expected_edges==len(edges)
    result=dict(status='PASS_ALL24_BACKGROUNDS_ALL58848_GROUPED_TIE_AP_AUC_ALL960_SLOTS_AND_LITERAL_REFERENCE_INTERFACES',
        finished_UTC=datetime.now(timezone.utc).isoformat(),audit_code_sha256=sha(__file__),producer_terminal_sha256=sha(O/'private_terminal.json'),
        ranked_records=checked,fixed_slots=len(slots),reference_interfaces=len(edges),no_sender_response_P_or_BH=True,same_cell_receptor_RNA_verified=False,causal_or_clinical_axis_validated=False)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()

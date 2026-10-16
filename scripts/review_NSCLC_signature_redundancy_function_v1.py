"""All24 measured signatures: redundancy and reference gene-set descriptions.

No clinical selection or new clinical tests. Enrichment describes gene lists,
not independent functional RNA responses, cell states or mechanisms.
"""
import hashlib,json,os,itertools,re
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata,fisher_exact
from scipy.cluster.hierarchy import linkage,fcluster
from scipy.spatial.distance import squareform
from statsmodels.stats.multitest import multipletests

COMMON=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary','response']
UNIT=['fraction','logN','mean_logRNA','mean_loggenes']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def corr(x,y,Z=None):
    a=rankdata(x).astype(float);b=rankdata(y).astype(float);a-=a.mean();b-=b.mean()
    if Z is not None:a-=Z@np.linalg.lstsq(Z,a,rcond=None)[0];b-=Z@np.linalg.lstsq(Z,b,rcond=None)[0]
    d=np.linalg.norm(a)*np.linalg.norm(b)
    return float(a@b/d) if d>1e-12 else np.nan


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);V=Path(json.loads((R/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    P=R/'NSCLC_context_native_programme_projection_20261016_v1'
    assert json.loads((P/'private_independent_audit.json').read_text())['status'].startswith('PASS_ALL56_SCOPE_')
    O=R/'NSCLC_signature_redundancy_function_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
        scope='All24 measured frozen signatures, all276 signature pairs, all same-family author-subtype diagnostics; no previous clinical q filter.',
        genotype='Report native signed raw-weight cosine, actual linear coefficient (weight/frozenSD) cosine, native gene-set overlap and original direction.',
        clusters='Within native source families only, complete-linkage on1-abs(effective-coefficient cosine), cut.05. Descriptive near-duplicate labels, no biological states or BH-family reduction.',
        score_redundancy='Same author subtype and cases, raw signature Spearman and rank residual correlation with same9 clinical and4 subtype controls; >=50 cases, complete rank14. Descriptive only, no P.',
        function='All24 x2 signed gene lists x all50 existing Hallmark sets; measured native intersect original model input gene universe; >=80%/10 Hallmark genes measured, hypergeometric Fisher upper tail, full2400 BH including unavailable p1.',
        category_diagnostics='Original abs-weight fractions of HSP/DNAJ-prefix, exact ribosomal-symbol pattern, MT-prefix and MALAT1. Gene families are not automatic technical-artifact calls.',
        boundary='Reference-list enrichment is not a gene-disjoint functional response measurement, new clinical validation, directed communication or formal state assignment.')
    (O/'plan_before_redundancy_and_reference_annotation.json').write_text(json.dumps(plan,indent=2)+'\n')
    s=pd.read_csv(P/'private_native_author_subtype_programme_scores.tsv.gz',sep='\t')
    slot=pd.read_csv(P/'private_all56_projection_scope.tsv',sep='\t');slot=slot[slot.status.str.startswith('PARTIAL_FROZEN')].copy();assert len(slot)==24
    params=pd.read_csv(R/'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv',sep='\t')
    native=pd.Index(pd.read_csv(R/'NSCLC_full_native_gene_units_20261016_v1/private_native_gene_axis.tsv',sep='\t').gene)
    c=pd.read_csv(R/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv',sep='\t').set_index('sampleID')
    slot['signature']=[f'{r.context}_{r.type}_{r.program}' for r in slot.itertuples()];slot=slot.sort_values('signature').reset_index(drop=True)
    vectors=np.zeros((24,len(native)));rawv=vectors.copy();sets=[];categories=[];lookups={}
    for i,r in enumerate(slot.itertuples()):
        p=params[params.context.eq(r.context)&params.type.eq(r.type)&params.program.eq(r.program)]
        ix=native.get_indexer(p.gene);keep=ix>=0;vectors[i,ix[keep]]=p.weight.to_numpy()[keep]/p.gene_sd.to_numpy()[keep];rawv[i,ix[keep]]=p.weight.to_numpy()[keep]
        sets.append(set(p.loc[keep,'gene']));lookups[r.signature]=p
        total=np.abs(p.weight).sum();g=p.gene.astype(str)
        masks={'HSP_DNAJ_prefix':g.str.startswith(('HSP','DNAJ')),'ribosomal_symbol_pattern':g.map(lambda x:bool(re.fullmatch(r'(RPS[0-9]{1,2}[A-Z]?|RPL[0-9]{1,2}[A-Z]?|RPLP[0-9]|RPSA)',x))),
            'MT_prefix':g.str.startswith('MT-'),'MALAT1':g.eq('MALAT1')}
        for name,mask in masks.items():categories.append(dict(signature=r.signature,source_family=r.source_family,category=name,genes=int(mask.sum()),original_absolute_weight_fraction=float(np.abs(p.weight[mask]).sum()/total)))
    norm=vectors/np.linalg.norm(vectors,axis=1)[:,None];nr=rawv/np.linalg.norm(rawv,axis=1)[:,None];cosine=norm@norm.T;rawcos=nr@nr.T;pair=[]
    for i,j in itertools.combinations(range(24),2):
        pair.append(dict(signature_A=slot.iloc[i].signature,signature_B=slot.iloc[j].signature,family_A=slot.iloc[i].source_family,family_B=slot.iloc[j].source_family,
            effective_coefficient_cosine=float(cosine[i,j]),original_weight_cosine=float(rawcos[i,j]),shared_native_genes=len(sets[i]&sets[j]),native_gene_Jaccard=len(sets[i]&sets[j])/len(sets[i]|sets[j])))
    groups=[]
    for family,f in slot.groupby('source_family',sort=True):
        ix=f.index.to_numpy();distance=np.clip(1-np.abs(cosine[np.ix_(ix,ix)]),0,1);np.fill_diagonal(distance,0)
        label=fcluster(linkage(squareform(distance,checks=False),method='complete'),t=.05,criterion='distance') if len(ix)>1 else [1]
        for index,lab in zip(ix,label):groups.append(dict(signature=slot.loc[index,'signature'],source_family=family,descriptive_weight_cluster=family+'_W'+str(lab)))
    pd.DataFrame(pair).to_csv(O/'private_all276_signature_weight_pairs.tsv',sep='\t',index=False)
    pd.DataFrame(groups).to_csv(O/'private_all24_descriptive_weight_clusters.tsv',sep='\t',index=False)
    pd.DataFrame(categories).to_csv(O/'private_all96_gene_family_weight_descriptions.tsv',sep='\t',index=False)
    diagnostic=[]
    for family,f in slot.groupby('source_family',sort=True):
        for a,b in itertools.combinations(f.itertuples(),2):
            aa=s[s.context.eq(a.context)&s.type.eq(a.type)&s.program.eq(a.program)]
            bb=s[s.context.eq(b.context)&s.type.eq(b.type)&s.program.eq(b.program)]
            for subtype,af in aa.groupby('sub_cell_type',sort=True):
                bf=bb[bb.sub_cell_type.eq(subtype)];af=af.set_index('sampleID');bf=bf.set_index('sampleID');keys=sorted(set(af.index)&set(bf.index));af=af.loc[keys];bf=bf.loc[keys];cc=c.loc[keys]
                use=af.n_cells.ge(20).to_numpy()&bf.n_cells.ge(20).to_numpy()&cc.clinical_eligible.to_numpy(bool);af=af.loc[use];bf=bf.loc[use];cc=cc.loc[use];n=len(af)
                assert np.array_equal(af.unit_index,bf.unit_index)
                row=dict(signature_A=a.signature,signature_B=b.signature,source_family=family,author_subtype=subtype,cases=n,raw_score_spearman=np.nan,full_control_partial_rank_rho=np.nan,status='INSUFFICIENT_CASES')
                if n>=50:
                    cov=np.column_stack([cc[COMMON].to_numpy(float),af[UNIT].to_numpy(float)]);sd=cov.std(0);D=np.column_stack([np.ones(n),(cov-cov.mean(0))/np.where(sd>1e-12,sd,1)])
                    row['raw_score_spearman']=corr(af.raw_weighted_score,bf.raw_weighted_score)
                    if (sd>1e-12).all() and np.linalg.matrix_rank(D)==14:
                        row.update(status='ESTIMABLE_DESCRIPTIVE_NO_P',full_control_partial_rank_rho=corr(af.raw_weighted_score,bf.raw_weighted_score,D))
                    else:row['status']='FULL_CONTROLS_NOT_IDENTIFIABLE'
                diagnostic.append(row)
    pd.DataFrame(diagnostic).to_csv(O/'private_all_same_subtype_score_redundancy.tsv',sep='\t',index=False,na_rep='NA')
    refs=R/'CRC_fine_pair_DIALOGUE_code_20261014_v1/hallmark50_public_reference_v3';prov=json.loads((refs/'reference_provenance.json').read_text());assert len(prov['records'])==50
    hallmarks={}
    for r in prov['records']:
        path=refs/(r['name']+'.json');assert sha(path)==r['sha256'];hallmarks[r['name']]=set(json.loads(path.read_text())[r['name']]['geneSymbols'])
    enrich=[]
    for r in slot.itertuples():
        p=lookups[r.signature];background=set((V/'state_after_DIALOGUE_inputs_20261016_v3'/f'context_{r.context:04d}'/(r.type+'_discovery_genes.tsv')).read_text().splitlines())&set(native)
        assert set(p.gene)&set(native)<=background
        for direction in ['positive','negative']:
            fg=set(p.loc[p.weight.gt(0) if direction=='positive' else p.weight.lt(0),'gene'])&background
            for name,reference in hallmarks.items():
                genes=reference&background;eligible=len(genes)>=10 and len(genes)/len(reference)>=.8
                hit=len(fg&genes);pv=1.;status='REFERENCE_UNIVERSE_COVERAGE_INSUFFICIENT'
                if eligible and fg:
                    contingency=[[hit,len(fg)-hit],[len(genes)-hit,len(background-fg-genes)]]
                    pv=float(fisher_exact(contingency,alternative='greater').pvalue);status='REFERENCE_ENRICHMENT_NOT_FUNCTIONAL_RESPONSE'
                elif eligible:status='EMPTY_SIGNED_FOREGROUND'
                enrich.append(dict(signature=r.signature,source_family=r.source_family,direction=direction,reference=name,universe_genes=len(background),foreground_genes=len(fg),
                    reference_genes_measured=len(genes),reference_genes_total=len(reference),overlap=hit,Fisher_p=pv,status=status))
    e=pd.DataFrame(enrich);assert len(e)==2400;e['BH_q_all2400']=multipletests(e.Fisher_p,method='fdr_bh')[1]
    e.to_csv(O/'private_all2400_signed_reference_enrichments.tsv',sep='\t',index=False)
    report=dict(status='ALL24_SIGNATURE_REDUNDANCY_AND2400_SIGNED_REFERENCE_DESCRIPTIONS_COMPLETE',finished_UTC=datetime.now(timezone.utc).isoformat(),
        signature_pairs=len(pair),same_subtype_diagnostics=len(diagnostic),signed_reference_records=len(e),no_new_clinical_tests=True,
        formal_cell_states_assigned=False,independent_functional_responses_measured=False,clinical_BH_family_changed=False,
        files={p.name:sha(p) for p in O.glob('*.tsv')})
    (O/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()

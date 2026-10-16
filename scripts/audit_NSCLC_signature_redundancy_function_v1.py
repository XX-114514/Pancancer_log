"""Independent scalar cosine, same-cell OLS and hypergeometric reference audit."""
import hashlib,json,os,itertools,re,math
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata,hypergeom
import statsmodels.api as sm

COMMON=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary','response']
UNIT=['fraction','logN','mean_logRNA','mean_loggenes']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def same(a,b):assert np.allclose(a,b,atol=1e-8,rtol=1e-8,equal_nan=True),(a,b)


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=R/'NSCLC_signature_redundancy_function_20261016_v1';out=O/'private_independent_audit.json';assert not out.exists()
    terminal=json.loads((O/'private_terminal.json').read_text())
    for name,h in terminal['files'].items():assert sha(O/name)==h
    V=Path(json.loads((R/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    params=pd.read_csv(R/'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv',sep='\t')
    native=set(pd.read_csv(R/'NSCLC_full_native_gene_units_20261016_v1/private_native_gene_axis.tsv',sep='\t').gene)
    P=R/'NSCLC_context_native_programme_projection_20261016_v1';scores=pd.read_csv(P/'private_native_author_subtype_programme_scores.tsv.gz',sep='\t')
    slot=pd.read_csv(P/'private_all56_projection_scope.tsv',sep='\t');slot=slot[slot.status.str.startswith('PARTIAL_FROZEN')]
    look={f'{r.context}_{r.type}_{r.program}':params[params.context.eq(r.context)&params.type.eq(r.type)&params.program.eq(r.program)] for r in slot.itertuples()};assert len(look)==24
    pair=pd.read_csv(O/'private_all276_signature_weight_pairs.tsv',sep='\t');assert len(pair)==276
    cosine={}
    for row in pair.itertuples():
        a=look[row.signature_A];b=look[row.signature_B]
        aa={r.gene:(float(r.weight),float(r.weight/r.gene_sd)) for r in a.itertuples() if r.gene in native}
        bb={r.gene:(float(r.weight),float(r.weight/r.gene_sd)) for r in b.itertuples() if r.gene in native}
        for j,col in [(0,'original_weight_cosine'),(1,'effective_coefficient_cosine')]:
            value=math.fsum(aa[g][j]*bb[g][j] for g in set(aa)&set(bb))/math.sqrt(math.fsum(x[j]**2 for x in aa.values())*math.fsum(x[j]**2 for x in bb.values()))
            same(value,getattr(row,col))
        assert row.shared_native_genes==len(set(aa)&set(bb));same(row.native_gene_Jaccard,len(set(aa)&set(bb))/len(set(aa)|set(bb)))
        cosine[tuple(sorted([row.signature_A,row.signature_B]))]=row.effective_coefficient_cosine
    cluster=pd.read_csv(O/'private_all24_descriptive_weight_clusters.tsv',sep='\t');assert cluster.signature.is_unique and set(cluster.signature)==set(look)
    for _,f in cluster.groupby('descriptive_weight_cluster'):
        assert f.source_family.nunique()==1
        for a,b in itertools.combinations(f.signature,2):assert abs(cosine[tuple(sorted([a,b]))])>=.95-1e-12
    for _,f in cluster.groupby('source_family'):
        groups=list(f.groupby('descriptive_weight_cluster'))
        for (_,a),(_,b) in itertools.combinations(groups,2):assert min(abs(cosine[tuple(sorted([x,y]))]) for x in a.signature for y in b.signature)<.95+1e-12
    categories=pd.read_csv(O/'private_all96_gene_family_weight_descriptions.tsv',sep='\t');assert len(categories)==96
    for row in categories.itertuples():
        p=look[row.signature];pairs=list(zip(p.gene,p.weight));predicate={
            'HSP_DNAJ_prefix':lambda g:g.startswith(('HSP','DNAJ')),
            'ribosomal_symbol_pattern':lambda g:re.fullmatch(r'(RPS[0-9]{1,2}[A-Z]?|RPL[0-9]{1,2}[A-Z]?|RPLP[0-9]|RPSA)',g) is not None,
            'MT_prefix':lambda g:g.startswith('MT-'),'MALAT1':lambda g:g=='MALAT1'}[row.category]
        assert row.genes==sum(predicate(g) for g,_ in pairs)
        value=math.fsum(abs(w) for g,w in pairs if predicate(g))/math.fsum(abs(w) for _,w in pairs);same(value,row.original_absolute_weight_fraction)
    c=pd.read_csv(R/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv',sep='\t').set_index('sampleID')
    d=pd.read_csv(O/'private_all_same_subtype_score_redundancy.tsv',sep='\t');assert len(d)==1106
    for row in d.itertuples():
        parts=[]
        for name in [row.signature_A,row.signature_B]:
            context,typ,program=name.split('_');f=scores[scores.context.eq(int(context))&scores.type.eq(typ)&scores.program.eq(program)&scores.sub_cell_type.eq(row.author_subtype)].set_index('sampleID');parts.append(f)
        a,b=parts;ids=sorted(set(a.index)&set(b.index));a=a.loc[ids];b=b.loc[ids];cc=c.loc[ids]
        use=a.n_cells.ge(20).to_numpy()&b.n_cells.ge(20).to_numpy()&cc.clinical_eligible.to_numpy(bool);a=a.loc[use];b=b.loc[use];cc=cc.loc[use]
        assert len(a)==row.cases and np.array_equal(a.unit_index,b.unit_index)
        if len(a)<50:assert row.status=='INSUFFICIENT_CASES';continue
        x=rankdata(a.raw_weighted_score);y=rankdata(b.raw_weighted_score);same(np.corrcoef(x,y)[0,1],row.raw_score_spearman)
        if row.status=='ESTIMABLE_DESCRIPTIVE_NO_P':
            cov=np.column_stack([cc[COMMON].to_numpy(float),a[UNIT].to_numpy(float)]);Z=np.column_stack([np.ones(len(a)),(cov-cov.mean(0))/cov.std(0)])
            ex=sm.OLS(x,Z).fit().resid;ey=sm.OLS(y,Z).fit().resid;same(np.corrcoef(ex,ey)[0,1],row.full_control_partial_rank_rho)
    refs=R/'CRC_fine_pair_DIALOGUE_code_20261014_v1/hallmark50_public_reference_v3';prov=json.loads((refs/'reference_provenance.json').read_text());reference={}
    for r in prov['records']:
        p=refs/(r['name']+'.json');assert sha(p)==r['sha256'];reference[r['name']]=set(json.loads(p.read_text())[r['name']]['geneSymbols'])
    t=pd.read_csv(O/'private_all2400_signed_reference_enrichments.tsv',sep='\t');assert len(t)==2400
    for row in t.itertuples():
        p=look[row.signature];context,typ,program=row.signature.split('_')
        bg=set((V/'state_after_DIALOGUE_inputs_20261016_v3'/f'context_{int(context):04d}'/(typ+'_discovery_genes.tsv')).read_text().splitlines())&native
        fg={g for g,w in zip(p.gene,p.weight) if (w>0 if row.direction=='positive' else w<0)}&bg;rr=reference[row.reference]&bg
        assert row.universe_genes==len(bg) and row.foreground_genes==len(fg) and row.reference_genes_measured==len(rr) and row.overlap==len(fg&rr)
        eligible=len(rr)>=10 and len(rr)/len(reference[row.reference])>=.8
        pv=float(hypergeom.sf(len(fg&rr)-1,len(bg),len(rr),len(fg))) if eligible and fg else 1.;same(pv,row.Fisher_p)
    order=np.argsort(t.Fisher_p.to_numpy(),kind='stable');expected=np.minimum.accumulate((t.Fisher_p.to_numpy()[order]*len(t)/np.arange(1,len(t)+1))[::-1])[::-1];expected=np.minimum(expected,1)
    same(t.BH_q_all2400.to_numpy()[order],expected)
    result=dict(status='PASS_ALL276_SCALAR_COSINES24_CLUSTER_RULES96_WEIGHT_DESCRIPTIONS1106_SCORE_AND2400_HYPERGEOMETRIC_RECORDS',
        finished_UTC=datetime.now(timezone.utc).isoformat(),producer_terminal_sha256=sha(O/'private_terminal.json'),audit_code_sha256=sha(__file__),
        new_clinical_or_coordination_tests=False,independent_functional_response_measured=False,formal_states_assigned=False)
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()

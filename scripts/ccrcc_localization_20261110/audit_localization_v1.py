"""Independent regional sum, exact-null and frozen-identity checks."""
from pathlib import Path
from datetime import datetime,timezone
import sys,json,hashlib,math,itertools
import numpy as np
import pandas as pd


def main(r,o):
    original=r/'current_identity_forward_v1'
    d=pd.read_csv(r/'spatial_pilot_v1/identity_candidates_PRIVATE.tsv.gz',sep='\t',index_col=0)
    sha=hashlib.sha256((r/'spatial_pilot_v1/identity_candidates_PRIVATE.tsv.gz').read_bytes()).hexdigest()
    assert sha=='d65cbcb6f8508f2598934d9de3c95a072fde4a53ba2871a13689728a8e262bde'
    assert d.QC_keep.sum()==358841
    lr=pd.read_csv(o/'LR108_exact_expectation_block_robustness.tsv',sep='\t')
    cm=pd.read_csv(o/'CM12_exact_expectation_block_robustness.tsv',sep='\t')
    assert len(lr)==108 and len(cm)==12
    assert lr.status.eq('MISSING_PANEL_COMPONENT').sum()==36
    assert lr.loc[lr.status.eq('MISSING_PANEL_COMPONENT'),'exact_expected'].isna().all()
    maxerr=0.; maxz=0.
    for table,regfile,keys in [(lr,'LR_regional_PRIVATE.tsv.gz',['sender','receiver','ligand','receptor','radius_um','tile_um']),
                               (cm,'CM_regional_PRIVATE.tsv.gz',['radius_um','tile_um','relation'])]:
        regional=pd.read_csv(o/regfile,sep='\t')
        for _,line in table.iterrows():
            q=np.ones(len(regional),bool)
            for k in keys:q &= regional[k].eq(line[k])
            a=regional.loc[q]
            if not len(a):assert line.get('status')=='MISSING_PANEL_COMPONENT';continue
            for expected,actual in [('observed','observed'),('exact_expected','expected'),('difference','excess')]:
                error=abs(a[actual].sum()-line[expected]);maxerr=max(maxerr,error);assert error<1e-6
            assert line.positive_blocks_at_least20==((a.opportunities>=20)&(a.excess>0)).sum()
            rest=(a.excess.sum()-a.excess)/(a.opportunities.sum()-a.opportunities)*100
            assert np.isclose(rest.min(),line.fixed_null_leave_one_block_out_min_pp,rtol=0,atol=1e-6)
            assert np.isclose(rest.max(),line.fixed_null_leave_one_block_out_max_pp,rtol=0,atol=1e-6)
    # Independent pandas group means; no reuse of producer's strata routine.
    z=np.load(original/'communication_inputs_PRIVATE.npz');inc=z['included'];es=z['es'];er=z['er'];values=z['expression']
    genes=['CCL19','CXCL12','CCR7','CXCR4'];p=d.iloc[inc]
    null=np.load(original/'communication_null_PRIVATE.npz')
    for tile in [250,500]:
        group=pd.DataFrame({'x':np.floor(p.x_centroid.to_numpy()/tile).astype(int),'y':np.floor(p.y_centroid.to_numpy()/tile).astype(int),'depth':np.floor(np.log2(p.rna_total.to_numpy())).astype(int),'identity':p.identity_candidate.to_numpy()})
        for j,g in enumerate(genes):group[g]=values[:,j].astype(int)
        means=group.groupby(['x','y','depth','identity'],sort=False)[genes].transform('mean').to_numpy()
        a=lr[(lr.tile_um==tile)&lr.status.ne('MISSING_PANEL_COMPONENT')]
        for _,line in a.iterrows():
            si=['Blood_EC','Mural','Lymphatic'].index(line.sender);ti=['B','T','NK'].index(line.receiver)
            take=(z['distance']<=line.radius_um)&(z['pair_code']==si*3+ti)
            expected=float(np.dot(means[es[take],genes.index(line.ligand)],means[er[take],genes.index(line.receptor)]))
            err=abs(expected-line.exact_expected);maxerr=max(maxerr,err);assert err<1e-6
            j=(([50,100].index(line.radius_um)*2+['CCL19','CXCL12'].index(line.ligand))*9+si*3+ti)
            sem=null[str(tile)][:,j].std(ddof=1)/math.sqrt(999)
            zz=abs(expected-null[str(tile)][:,j].mean())/(sem or 1)
            maxz=max(maxz,zz);assert zz<6
    graph=np.load(original/'graph_PRIVATE.npz');imm=graph['immune_indices']
    labels=d.identity_candidate.to_numpy();is_ec=labels[graph['sender']]=='Blood_EC'
    # Check independently 20 edge-bearing ECs per radius/grid with integer combinations.
    eclist=np.unique(graph['sender'][is_ec]);picks=np.unique(np.linspace(0,len(eclist)-1,20).astype(int));checks=0
    for tile in [250,500]:
        immune=d.iloc[imm];keys=list(zip(np.floor(immune.x_centroid/tile).astype(int),np.floor(immune.y_centroid/tile).astype(int),np.floor(np.log2(immune.rna_total)).astype(int)))
        from collections import Counter,defaultdict
        counts=defaultdict(Counter)
        for k,name in zip(keys,labels[imm]):counts[k][name]+=1
        from localize_existing_identity_v1 import strata,shared_expectation
        st=strata(immune[['x_centroid','y_centroid']].to_numpy(),immune.rna_total.to_numpy(),tile)
        ec=np.flatnonzero(d.QC_keep.to_numpy()&(labels=='Blood_EC'));ecmap={x:i for i,x in enumerate(ec)}
        codes=pd.Categorical(labels[imm],categories=['B','T','NK']).codes
        for rad in [50,100]:
            good=is_ec&(graph['distance']<=rad)
            owner=np.array([ecmap[x] for x in graph['sender'][good]])
            fast=shared_expectation(owner,st[graph['receiver_row'][good]],codes,st,len(ec))
            for e in eclist[picks]:
                take=good&(graph['sender']==e);g=Counter(keys[j] for j in graph['receiver_row'][take])
                products=[]
                for forbidden in [('B',),('NK',),('B','NK')]:
                    probability=1.
                    for key,k in g.items():
                        c=counts[key];n=sum(c.values());allowed=n-sum(c[a] for a in forbidden)
                        probability*=math.comb(allowed,k)/math.comb(n,k) if allowed>=k else 0.
                    products.append(probability)
                ref=1-products[0]-products[1]+products[2]
                err=abs(ref-fast[ecmap[e]]);maxerr=max(maxerr,err);assert err<1e-7;checks+=1
    # Exhaustive counterexample: B/NK identity draws without replacement are dependent.
    from localize_existing_identity_v1 import shared_expectation
    toy=np.array([0,0,1,2]);allst=np.zeros(4,dtype=int)
    exact=shared_expectation(np.zeros(2,dtype=int),np.zeros(2,dtype=int),toy,allst,1)[0]
    truth=np.mean([0 in toy[list(q)] and 2 in toy[list(q)] for q in itertools.combinations(range(4),2)])
    assert abs(exact-truth)<1e-12
    receipt={'status':'PASS_frozen_identity_region_sums_exact_LR_and_hypergeometric_checks','UTC':datetime.now(timezone.utc).isoformat(),'LR108_slots':len(lr),'CM12_slots':len(cm),'CCL21_missing_slots':36,'independent_LR_expectations_checked':72,'shared_EC_integer_combination_checks':checks,'shared_EC_exhaustive_toy_PASS':True,'maximum_arithmetic_error':maxerr,'max_LR_exact_expectation_difference_from_prior_MC_mean_in_MC_SE':maxz,'identity_sha256':sha,'new_P_values':0,'no_new_identity_gates':True,'not_histology_or_patient_replication':True}
    (o/'independent_audit.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)


if __name__=='__main__':main(Path(sys.argv[1]),Path(sys.argv[2]))

"""Fresh native-count per-gene readback and all aggregate/tail checks."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd,anndata as ad

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def frame(p):
    with np.load(p,allow_pickle=True) as z:return pd.DataFrame(z['data'],index=z['index'],columns=z['columns'])
def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);D=R/'state_lineage_coexpression_20261014_v2';plan=json.loads(Path(os.environ['PANCANCER_LINEAGE_DIAGNOSTIC_PLAN']).read_text());inv=json.loads((D/'inventory.json').read_text());src=json.loads((R/'state_singlecell_inputs_acceptance_20261014_v2/accepted_input_inventory.json').read_text());summary=pd.read_csv(D/'all110_source_major_panel_diagnostics.tsv',sep='\t');nodes=pd.read_csv(D/'all_source_node_foreign_RNA_diagnostics.tsv',sep='\t');tail=pd.read_csv(D/'all990_discovery_module_tail_diagnostics.tsv',sep='\t');coverage=pd.read_csv(D/'all55_gene_panel_coverage.tsv',sep='\t').fillna('');checks=0;aggs=0;actual_node_keys=set();training=[]
    for n,h in inv['output_sha256'].items():assert sha(D/n)==h
    for sid,item in src['verified_dataset_files'].items():
        assert sha(Path(item['path']))==item['sha256'];a=ad.read_h5ad(item['path']);total=a.obs.full_RNA_total.to_numpy();z=np.load(D/('private_'+sid+'_panel_quantities.npz'),allow_pickle=False);assert np.array_equal(z['cell_ids'],a.obs_names.to_numpy(str));q=pd.DataFrame(z['values'],index=a.obs_names,columns=z['columns']);union=sorted({g for genes in plan['panels'].values() for g in genes if g in a.var_names});x=a.X[:,a.var_names.get_indexer(union)].tocsc();pos={g:i for i,g in enumerate(union)};fresh={}
        for panel,genes in plan['panels'].items():
            count=np.zeros(a.n_obs,dtype=np.int64);det=np.zeros(a.n_obs,dtype=np.int64);missing=[]
            for gene in genes:
                if gene not in pos:missing.append(gene);continue
                j=pos[gene];lo,hi=x.indptr[j:j+2];rows=x.indices[lo:hi];vals=x.data[lo:hi];count[rows]+=vals.astype(np.int64);det[rows]+=vals>0
            frac=count/total;assert np.array_equal(count,q[panel+'__umis']) and np.array_equal(det,q[panel+'__det']) and np.array_equal(frac,q[panel+'__frac']);checks+=3*a.n_obs;fresh[panel]=dict(det=det,frac=frac)
            cv=coverage[(coverage.gse_id==sid)&(coverage.panel==panel)].iloc[0];assert cv.measured_genes==len(genes)-len(missing) and cv.missing_genes==';'.join(missing)
            for major in ['T cell','Macrophage']:
                keep=a.obs.major.eq(major).to_numpy();r=summary[(summary.gse_id==sid)&(summary.major==major)&(summary.panel==panel)].iloc[0];assert r.cells==keep.sum();vv=frac[keep];assert np.allclose([np.median(vv),np.quantile(vv,.95),np.quantile(vv,.99),max(vv)],[r.fraction_median,r.fraction_q95,r.fraction_q99,r.fraction_max],atol=1e-12);assert r.multi_marker_cells==int(((det>=3)&keep).sum());aggs+=1
        for major in ['T cell','Macrophage']:
            keep=a.obs.major.eq(major).to_numpy();core=fresh['T_core' if major=='T cell' else 'Myeloid_core']['det']>=2
            for panel in plan['panels']:
                r=summary[(summary.gse_id==sid)&(summary.major==major)&(summary.panel==panel)].iloc[0]
                for t in plan['RNA_fraction_thresholds']:
                    tag=str(t).replace('.','_');good=(fresh[panel]['det']>=3)&(fresh[panel]['frac']>=t);assert r['multi_at_fraction_'+tag]==int((good&keep).sum());assert r['samecell_owncore_and_multi_at_fraction_'+tag]==int((good&keep&core).sum())
            # Unique strings avoids the exact categorical empty-group bug.
            labels=a.obs.subCluster.astype(str).to_numpy()
            for node in np.unique(labels[keep]):
                mask=keep&(labels==node)
                for panel in plan['foreign_panels']:
                    r=nodes[(nodes.gse_id==sid)&(nodes.major==major)&(nodes.subCluster==node)&(nodes.panel==panel)].iloc[0];good=(fresh[panel]['det']>=3)&(fresh[panel]['frac']>=.01);assert r.cells==mask.sum() and r.multi_marker_fraction1pct_cells==int((good&mask).sum()) and r.samecell_owncore_cells==int((good&mask&core).sum());actual_node_keys.add((sid,major,node,panel))
        training.append(q.loc[a.obs.state_training_selected.to_numpy(bool)])
    assert actual_node_keys==set(nodes[['gse_id','major','subCluster','panel']].itertuples(index=False,name=None));tr=pd.concat(training);assert len(tr)==12000
    for major,name in [('T cell','T_cell'),('Macrophage','Macrophage')]:
        for k in [10,15,20]:
            u=frame(R/'cNMF_state_contrast_20261014_v1'/name/'cnmf/cnmf_tmp'/f'cnmf.usages.k_{k}.dt_0_5.consensus.df.npz');u=u.div(u.sum(axis=1),axis=0);v=tr.loc[u.index]
            for component in u.columns:
                # Python stable sort is independent of producer numpy argsort.
                ix=sorted(range(len(u)),key=lambda i:(-u.iloc[i][component],i))[:300]
                for panel in plan['panels']:
                    r=tail[(tail.major==major)&(tail.K==k)&(tail.component==component)&(tail.panel==panel)].iloc[0];det=v[panel+'__det'].to_numpy();frac=v[panel+'__frac'].to_numpy();assert np.isclose(np.median(frac[ix]),r.tail_fraction_median,atol=1e-12);assert r.tail_multi_marker_fraction1pct_cells==int(((det[ix]>=3)&(frac[ix]>=.01)).sum());assert r.all6000_multi_marker_fraction1pct_cells==int(((det>=3)&(frac>=.01)).sum())
    old=R/'state_lineage_coexpression_20261014_v1/private_GSE148071_panel_quantities.npz';new=D/'private_GSE148071_panel_quantities.npz';za=np.load(old,allow_pickle=False);zb=np.load(new,allow_pickle=False)
    for k in za.files:assert np.array_equal(za[k],zb[k])
    O=R/'state_lineage_coexpression_audit_20261014_v1';assert not O.exists();O.mkdir();out=dict(status='PASS_ALL240533_NATIVE_PANEL_VALUES_AGGREGATES_AND990_TAILS',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(D/'inventory.json'),native_cell_panel_values_checked=checks,source_major_aggregate_rows=aggs,source_node_rows=len(nodes),module_tail_rows=len(tail),all_actual_node_keys_exact=True,failed_v1_first_source_private_quantities_unchanged=True,original_identities_changed=False,technical_or_biological_origin_determined=False,formal_state_labels=False);(O/'independent_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()

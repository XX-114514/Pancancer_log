"""Full native-cell descriptive lineage RNA check; no reassignment/diagnosis."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc
import numpy as np,pandas as pd,anndata as ad

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()
def frame(p):
    with np.load(p,allow_pickle=True) as z:return pd.DataFrame(z['data'],index=z['index'],columns=z['columns'])
def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=Path(os.environ['PANCANCER_LINEAGE_DIAGNOSTIC_PLAN']);plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256'];I=R/'state_singlecell_inputs_acceptance_20261014_v2/accepted_input_inventory.json';assert sha(I)==plan['input_inventory_sha256'];inv=json.loads(I.read_text());panels=plan['panels'];O=R/'state_lineage_coexpression_20261014_v2';assert not O.exists();O.mkdir()
    summaries=[];nodes=[];coverage=[];tail=[];train={};total=0
    for sid,item in inv['verified_dataset_files'].items():
        p=Path(item['path']);assert sha(p)==item['sha256'];a=ad.read_h5ad(p);assert a.obs_names.is_unique and a.var_names.is_unique;assert np.isfinite(a.X.data).all() and (a.X.data>=0).all();assert np.array_equal(np.asarray(a.X.sum(axis=1)).ravel(),a.obs.full_RNA_total.to_numpy());total+=a.n_obs;expr={};full={}
        for name,genes in panels.items():
            ix=a.var_names.get_indexer(genes);matched=ix>=0;full[name]=bool(matched.all());x=a.X[:,ix[matched]].toarray();det=(x>0).sum(axis=1);umis=x.sum(axis=1);frac=umis/a.obs.full_RNA_total.to_numpy();expr[name]=dict(det=det,umis=umis,frac=frac)
            coverage.append(dict(gse_id=sid,panel=name,required_genes=len(genes),measured_genes=int(matched.sum()),missing_genes=';'.join(np.array(genes)[~matched]),full_panel_measured=full[name]))
        # Per-cell panel quantities stay private; IDs never in public aggregate tables.
        private=pd.DataFrame(index=a.obs_names)
        for name,e in expr.items():
            for key in ['det','umis','frac']:private[name+'__'+key]=e[key]
        selected=a.obs.state_training_selected.to_numpy(bool);train[sid]=private.loc[selected].copy();np.savez_compressed(O/('private_'+sid+'_panel_quantities.npz'),cell_ids=a.obs_names.to_numpy(str),columns=private.columns.to_numpy(str),values=private.to_numpy())
        for major in ['T cell','Macrophage']:
            keep=a.obs.major.eq(major).to_numpy();identity='T_core' if major=='T cell' else 'Myeloid_core';core=expr[identity]['det']>=2
            for name,e in expr.items():
                vals=e['frac'][keep];multi=e['det']>=3
                row=dict(gse_id=sid,major=major,panel=name,cells=int(keep.sum()),libraries=int(a.obs.loc[keep].sample_uid.nunique()),inherited_analysis_blocks=int(a.obs.loc[keep].analysis_block.nunique()),full_panel_measured=full[name],fraction_median=float(np.median(vals)),fraction_q95=float(np.quantile(vals,.95)),fraction_q99=float(np.quantile(vals,.99)),fraction_max=float(vals.max()),multi_marker_cells=int((multi&keep).sum()),own_core_atleast2_cells=int((core&keep).sum()))
                for threshold in plan['RNA_fraction_thresholds']:
                    good=multi&(e['frac']>=threshold);tag=str(threshold).replace('.','_');row['multi_at_fraction_'+tag]=int((good&keep).sum());row['samecell_owncore_and_multi_at_fraction_'+tag]=int((good&core&keep).sum())
                summaries.append(row)
            for node,meta in a.obs.loc[keep].groupby('subCluster',observed=True):
                ii=a.obs_names.get_indexer(meta.index);assert len(ii)>=1
                for name in plan['foreign_panels']:
                    e=expr[name];good=(e['det'][ii]>=3)&(e['frac'][ii]>=.01);nodes.append(dict(gse_id=sid,major=major,subCluster=node,panel=name,cells=len(ii),full_panel_measured=full[name],multi_marker_fraction1pct_cells=int(good.sum()),samecell_owncore_cells=int((good&core[ii]).sum()),fraction_q99=float(np.quantile(e['frac'][ii],.99))))
        print(json.dumps(dict(gse_id=sid,cells=a.n_obs,phase='ALL_NATIVE_CELLS_PANEL_COUNTS_COMPLETE')),flush=True);del a,private;gc.collect()
    assert total==240533;training=pd.concat(train.values());assert len(training)==12000 and training.index.is_unique
    for major,name in [('T cell','T_cell'),('Macrophage','Macrophage')]:
        for k in [10,15,20]:
            u=frame(R/'cNMF_state_contrast_20261014_v1'/name/'cnmf/cnmf_tmp'/f'cnmf.usages.k_{k}.dt_0_5.consensus.df.npz');u=u.div(u.sum(axis=1),axis=0);v=training.loc[u.index];assert len(v)==6000
            for component in u.columns:
                # Exactly300 largest normalized usage; tie broken by frozen order, diagnostic only.
                ii=np.argsort(-u[component].to_numpy(),kind='stable')[:300]
                for panel in panels:
                    det=v[panel+'__det'].to_numpy();frac=v[panel+'__frac'].to_numpy();tail.append(dict(major=major,K=k,component=int(component),panel=panel,tail_cells=300,tail_definition='top5pct_normalized_usage_diagnostic_not_state_threshold',tail_fraction_median=float(np.median(frac[ii])),tail_multi_marker_fraction1pct_cells=int(((det[ii]>=3)&(frac[ii]>=.01)).sum()),all6000_multi_marker_fraction1pct_cells=int(((det>=3)&(frac>=.01)).sum())))
    for n,x in [('all110_source_major_panel_diagnostics.tsv',pd.DataFrame(summaries)),('all_source_node_foreign_RNA_diagnostics.tsv',pd.DataFrame(nodes)),('all55_gene_panel_coverage.tsv',pd.DataFrame(coverage)),('all990_discovery_module_tail_diagnostics.tsv',pd.DataFrame(tail))]:x.to_csv(O/n,sep='\t',index=False)
    out=dict(status='COMPLETE_NATIVE_LINEAGE_RNA_DIAGNOSTICS_NOT_IDENTITY_OR_STATE_DECISION',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),cells=total,training_cells=12000,source_major_panels=len(summaries),gene_coverage_records=len(coverage),module_tail_records=len(tail),identity_labels_changed=False,technical_vs_biological_origin_determined=False,formal_state_labels=False,output_sha256={p.name:sha(p) for p in O.iterdir() if p.is_file() and not p.name.startswith('private_')});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()

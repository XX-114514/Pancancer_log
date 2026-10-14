"""Check prior selected-panel major denominators against full source units."""
import json,os
from datetime import datetime,timezone
from pathlib import Path
import h5py,pandas as pd
from anndata._io.specs import read_elem


def main():
    project=Path(os.environ['PANCANCER_PROJECT_ROOT']);root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    annotation=project/'Project_v3/annotation_rebuild/20261010_uniform_expression_annotation_v1'
    out=root/'old_MCP_parent_denominator_scope_audit_20261014_v2';assert not out.exists();out.mkdir()
    catalogue=pd.read_csv(annotation/'all_node_number_catalog.tsv',sep='\t')
    required={'node_id','majorCluster'}
    if not required.issubset(catalogue):
        print(json.dumps(dict(status='CATALOGUE_SCHEMA_REVIEW_REQUIRED_NO_NUMERICAL_CHANGE',fields=list(catalogue.columns))))
        return
    records=[];detail=[]
    for source in ['GSE161529','GSE182109','GSE183904']:
        path=root/'state_singlecell_inputs_20261014_v1/private_native_sources'/(source+'.h5ad')
        with h5py.File(path) as f:obs=read_elem(f['obs'])
        source_units=pd.read_csv(annotation/'by_gse'/source/'units.tsv',sep='\t')
        if not {'sample_uid','node_index','n_cells'}.issubset(source_units):
            print(json.dumps(dict(status='SOURCE_UNIT_SCHEMA_REVIEW_REQUIRED_NO_NUMERICAL_CHANGE',fields=list(source_units.columns))))
            return
        source_units['major']=catalogue.majorCluster.to_numpy()[source_units.node_index]
        complete=source_units.groupby(['sample_uid','major']).n_cells.sum()
        selected=obs.groupby(['sample_uid','major'],observed=True).size()
        for key,count in selected.items():
            assert key in complete.index and complete.loc[key]>=count
            detail.append(dict(source=source,sample_uid=key[0],major=key[1],selected_panel_cells=int(count),full_source_major_cells=int(complete.loc[key])))
        for major in ['T cell','Macrophage']:
            rows=[r for r in detail if r['source']==source and r['major']==major]
            records.append(dict(source=source,major=major,source_sample_units=len(rows),
              samples_with_panel_below_full_major=sum(r['selected_panel_cells']<r['full_source_major_cells'] for r in rows),
              previous_fraction_scope='selected13node_panel_major_denominator',
              original_unadjusted_MCP_P_BH_changed=False))
    pd.DataFrame(detail).to_csv(out/'private_selected_vs_full_major_counts.tsv',sep='\t',index=False)
    terminal=dict(status='PREVIOUS_SELECTED_PANEL_PARENT_DENOMINATOR_SCOPE_AUDITED',finished_UTC=datetime.now(timezone.utc).isoformat(),summary=records)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()

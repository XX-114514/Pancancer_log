"""Reconcile exact stored aliases on unique, measured stable features only.

Keep the v1 check, original weights, and 95% coverage criterion unchanged.
"""
import hashlib
import json
import os
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem


def main():
    project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    contract_path=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1/config/source_contract_v2_private.json'
    contract=json.loads(contract_path.read_text())
    rec=next(r for r in contract['records'] if r['dataset']=='GSE286827')
    out=root/'HNSCC_frozen_MCP_gene_axis_20261014_v2'
    assert not out.exists()
    g=h5py.File(rec['path'],'r') if str(rec['path']).endswith('.h5ad') else zarr.open_group(rec['path'],mode='r')
    var=read_elem(g[rec['var_key']]);columns=np.arange(len(var))
    assert {'source_gene_aliases','feature_identity_status'}.issubset(var.columns)
    if rec.get('new_parent_dataset'):
        source_order=list(read_elem(g['uns'])['source_order'])
        measured=np.asarray(g['varm/measured_in_source'])[:,source_order.index(rec['new_parent_dataset'])]
        columns=np.flatnonzero(measured);var=var.iloc[columns].copy()
    names=defaultdict(set);old_names=defaultdict(set)
    for column,(index,row) in zip(columns,var.iterrows()):
        if str(index).startswith('AMBIGUOUS::') or row.feature_identity_status=='ambiguous_source_feature':
            continue
        known={str(index)}
        for field in ['gene_symbol','feature_name','gene_name_original']:
            if field in var and str(row[field]) not in ['', 'nan', 'None', '__SOURCE_VALUE_MISSING__']:
                known.add(str(row[field]))
        for name in known:old_names[name].add(int(column))
        # Stored aliases were collected by stable source feature identity.
        # No fuzzy, case-insensitive, or expression-driven matching.
        if row.feature_identity_status=='stable_Ensembl_ID':
            known.update(n for n in str(row.source_gene_aliases).split('|') if n not in ['', 'nan', 'None'])
        for name in known:names[name].add(int(column))
    out.mkdir()
    summaries=[]
    for endpoint in ['Tnode9','Macnode18']:
        path=root/'CRC_fine_pair_DIALOGUE_20261014_v1/private_transfer_fixed_weights_v1'/(endpoint+'_MCP1_genes.tsv')
        w=pd.read_csv(path,sep='\t')
        assert w.gene.is_unique and np.isfinite(w.weight).all()
        w['native_column']=[next(iter(names[gene])) if len(names[gene])==1 else -1 for gene in w.gene]
        w['match_kind']=['previous_exact_name' if len(old_names[gene])==1 else
                         'stored_stable_ID_alias' if len(names[gene])==1 else
                         'ambiguous' if len(names[gene])>1 else 'unmeasured' for gene in w.gene]
        # Distinct frozen genes cannot silently share a physical column.
        duplicated=w.native_column.ge(0)&w.native_column.duplicated(keep=False)
        w.loc[duplicated,'match_kind']='frozen_gene_physical_column_collision'
        w.loc[duplicated,'native_column']=-1
        measured=w.native_column.ge(0)
        fraction=float(np.abs(w.loc[measured,'weight']).sum()/np.abs(w.weight).sum())
        summary=dict(endpoint=endpoint,frozen_genes=len(w),measured_genes=int(measured.sum()),
          stored_alias_recovered=int(w.match_kind.eq('stored_stable_ID_alias').sum()),
          ambiguous_genes=int(w.match_kind.eq('ambiguous').sum()),column_collision_genes=int(duplicated.sum()),
          absolute_weight_retained=fraction,complete_score_measurable=bool(measured.all()),partial95_measurable=bool(fraction>=.95))
        summaries.append(summary)
        w.to_csv(out/(endpoint+'_private_frozen_gene_columns.tsv'),sep='\t',index=False)
    pd.DataFrame(summaries).to_csv(out/'aggregate_frozen_program_gene_coverage.tsv',sep='\t',index=False)
    terminal=dict(status='EXACT_STORED_ALIAS_MEASUREMENT_AUDIT_COMPLETE_NO_SCORES',
      finished_UTC=datetime.now(timezone.utc).isoformat(),source_contract_sha256=hashlib.sha256(contract_path.read_bytes()).hexdigest(),
      producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),coverage=summaries,
      weights_refitted=False,threshold_changed=False,clinical_tests_executed=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(terminal))


if __name__=='__main__':main()

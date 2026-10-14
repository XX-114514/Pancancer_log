"""Inspect existing supplemental ICB programme-validation opportunities.
Selected exports cannot validate population composition. No new downloads.
"""
import json,os
from datetime import datetime,timezone
from pathlib import Path
import h5py,numpy as np,pandas as pd
from anndata._io.specs import read_elem


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    atlas=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    contract=json.loads((atlas/'config/source_contract_v2_private.json').read_text())
    records=[r for r in contract['records'] if r['dataset'].startswith('ICB_')]
    assert records
    paths={r['path'] for r in records};assert len(paths)==1
    out=root/'ICB_fixed_programme_metadata_review_20261014_v1';assert not out.exists();out.mkdir()
    with h5py.File(next(iter(paths))) as handle:
        obs=read_elem(handle['obs']);var=read_elem(handle['raw/var'])
    rows=[];schemas=[];coverage=[]
    for record in records:
        md=obs.iloc[np.load(record['row_positions_ref'])].copy()
        assert len(md)==record['n_cells']
        columns=[str(c) for c in md.columns]
        schemas.append(dict(source=record['dataset'],role=record['role'],fields=columns))
        # Report only field-level categories, never donor or cell keys.
        categories={c:md[c].astype(str).value_counts().to_dict() for c in columns
                    if c in ['pre_post','Resp','Response','response','response_group','author_cell_type_update','author_cell_type','Cancer_type_update','Treatment','treatment']}
        rows.append(dict(source=record['dataset'],role=record['role'],n_cells=len(md),
          patient_tokens=md['PMID_donor_id'].astype(str).nunique(),categories=categories,
          population_composition_eligible=False,independent_validation_accepted=False))
    aliases={}
    for i,(name,row) in enumerate(var.iterrows()):
        names={str(name)}
        for col in ['gene_symbol','feature_name','gene_name_original']:
            if col in var and str(row[col]) not in ['','nan','None']:names.add(str(row[col]))
        for name in names:aliases.setdefault(name,set()).add(i)
    for endpoint in ['Tnode9','Macnode18']:
        w=pd.read_csv(root/'CRC_fine_pair_DIALOGUE_20261014_v1/private_transfer_fixed_weights_v1'/(endpoint+'_MCP1_genes.tsv'),sep='\t')
        take=np.array([len(aliases.get(g,set()))==1 for g in w.gene])
        coverage.append(dict(endpoint=endpoint,frozen_genes=len(w),measured_genes=int(take.sum()),
          absolute_weight_measurable=float(np.abs(w.weight[take]).sum()/np.abs(w.weight).sum())))
    (out/'private_native_schema.json').write_text(json.dumps(schemas,indent=2)+'\n')
    (out/'aggregate_opportunities.json').write_text(json.dumps(dict(sources=rows,coverage=coverage),indent=2)+'\n')
    print(json.dumps(dict(status='EXISTING_ICB_METADATA_AND_GENE_COVERAGE_AUDITED_NO_SCORES',
      finished_UTC=datetime.now(timezone.utc).isoformat(),schemas=schemas,sources=rows,coverage=coverage)))


if __name__=='__main__':main()

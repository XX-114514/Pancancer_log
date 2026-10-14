"""Reconcile exact stable IDs and audit clinical/type scope in existing ICB.
No programme scores or clinical effect estimates; all patient keys private.
"""
import json,os,re
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
import h5py,numpy as np,pandas as pd
from anndata._io.specs import read_elem


ID=re.compile(r'ENSG[0-9]{11}(?:\.[0-9]+)?')


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    atlas=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    contract=json.loads((atlas/'config/source_contract_v2_private.json').read_text())
    records=[r for r in contract['records'] if r['dataset'].startswith('ICB_')]
    paths={r['path'] for r in records};assert len(paths)==1
    out=root/'ICB_clinical_programme_scope_20261014_v2';assert not out.exists();out.mkdir()
    with h5py.File(next(iter(paths))) as handle:
        obs=read_elem(handle['obs']);var=read_elem(handle['raw/var'])
    rec=next(r for r in contract['records'] if r['dataset']=='GSE286827')
    import zarr
    parent=zarr.open_group(rec['path'],mode='r');parent_var=read_elem(parent[rec['var_key']])
    by_ID=defaultdict(set)
    for i,(name,row) in enumerate(var.iterrows()):
        values=[str(name)]+[str(row[c]) for c in var if c in ['gene_id','feature_id','gene_id_original','ensembl_id']]
        for identifier in values:
            if ID.fullmatch(identifier):by_ID[identifier.split('.')[0]].add(i)
    names=defaultdict(set)
    for name,row in parent_var.iterrows():
        if not ID.fullmatch(str(name)) or row.get('feature_identity_status')!='stable_Ensembl_ID':continue
        known={str(row.get('gene_symbol',''))}|set(str(row.get('source_gene_aliases','')).split('|'))
        for alias in known:
            if alias:names[alias].add(str(name).split('.')[0])
    exact=defaultdict(set)
    for i,(name,row) in enumerate(var.iterrows()):
        aliases={str(name)}|{str(row[c]) for c in ['gene_symbol','feature_name','gene_name_original'] if c in var}
        for alias in aliases:
            if alias not in ['','nan','None']:exact[alias].add(i)
    coverage=[]
    for endpoint in ['Tnode9','Macnode18']:
        w=pd.read_csv(root/'CRC_fine_pair_DIALOGUE_20261014_v1/private_transfer_fixed_weights_v1'/(endpoint+'_MCP1_genes.tsv'),sep='\t')
        columns=[];modes=[]
        for gene in w.gene:
            candidates=set(exact[gene]);mode='native_exact_name'
            if not candidates and len(names[gene])==1:
                candidates=set(by_ID.get(next(iter(names[gene])),set()));mode='exact_stored_stable_ID_bridge'
            columns.append(next(iter(candidates)) if len(candidates)==1 else -1)
            modes.append(mode if len(candidates)==1 else 'unmeasurable_or_ambiguous')
        w['native_column']=columns;w['match_mode']=modes
        duplicates=w.native_column.ge(0)&w.native_column.duplicated(keep=False)
        w.loc[duplicates,'native_column']=-1;w.loc[duplicates,'match_mode']='physical_column_collision'
        valid=w.native_column.ge(0)
        coverage.append(dict(endpoint=endpoint,native_features=len(var),native_var_fields=list(var.columns),
          stable_ID_features=len(by_ID),frozen_genes=len(w),measured_genes=int(valid.sum()),
          absolute_weight_measurable=float(np.abs(w.weight[valid]).sum()/np.abs(w.weight).sum()),
          exact_ID_bridge_recovered=int(w.match_mode.eq('exact_stored_stable_ID_bridge').sum())))
        w.to_csv(out/(endpoint+'_private_frozen_gene_columns.tsv'),sep='\t',index=False)
    aggregates=[];units=[]
    for rec in records:
        md=obs.iloc[np.load(rec['row_positions_ref'])]
        clinical_fields=['outcome','Combined_outcome','Treatment.or.Mode.of.Action','Cancer_type_update','Primary_or_met']
        source_units=[]
        for (patient,time),g in md.groupby(['PMID_donor_id','pre_post'],observed=True):
            counts=g.author_cell_type_update.astype(str).value_counts()
            unit=dict(source=rec['dataset'],patient_token=str(patient),phase=str(time),n_cells=len(g),
              author_macrophage_cells=int(counts.get('Macrophage',0)+counts.get('TAM',0)),
              author_unsplit_myeloid_cells=int(counts.get('Myeloid_cell',0)),
              author_CD8_cells=int(sum(n for name,n in counts.items() if name.startswith('CD8'))),
              author_unsplit_T_cells=int(counts.get('T_cell',0)))
            for field in clinical_fields:
                values=sorted(g[field].astype(str).unique())
                unit[field]='|'.join(values);unit[field+'_conflict']=len(values)>1
            source_units.append(unit);units.append(unit)
        frame=pd.DataFrame(source_units);pre=frame[frame.phase.eq('Pre')]
        aggregates.append(dict(source=rec['dataset'],role=rec['role'],baseline_patients=len(pre),
          all_phase_outcomes=frame.outcome.value_counts().to_dict(),
          baseline_combined_outcomes=pre.Combined_outcome.value_counts().to_dict(),
          baseline_treatment_categories=pre['Treatment.or.Mode.of.Action'].value_counts().to_dict(),
          baseline_any_clinical_field_conflicts=int(pre[[c+'_conflict' for c in clinical_fields]].any(axis=1).sum()),
          baseline_atleast20_author_macrophage=int(pre.author_macrophage_cells.ge(20).sum()),
          baseline_atleast20_unsplit_myeloid=int(pre.author_unsplit_myeloid_cells.ge(20).sum()),
          baseline_atleast20_both_macrophage_CD8=int((pre.author_macrophage_cells.ge(20)&pre.author_CD8_cells.ge(20)).sum()),
          baseline_atleast20_both_unsplit_myeloid_T=int((pre.author_unsplit_myeloid_cells.ge(20)&pre.author_unsplit_T_cells.ge(20)).sum()),
          population_composition_eligible=False,clinical_endpoint_definition_accepted=False))
    pd.DataFrame(units).to_csv(out/'private_source_patient_phase_scope.tsv',sep='\t',index=False)
    report=dict(status='ALL8_ICB_SCOPE_AND_STABLE_ID_GENE_AUDIT_COMPLETE_NO_TESTS',finished_UTC=datetime.now(timezone.utc).isoformat(),
      coverage=coverage,sources=aggregates,no_new_expression_download=True,weights_changed=False)
    (out/'aggregate_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()

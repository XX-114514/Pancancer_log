"""Independent sparse cell-set intersections for all frozen receptor edges."""
import hashlib,json,os
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
import anndata as ad
import numpy as np,pandas as pd
from diagnose_frozen_MCP_two_end_functions_v1 import hnscc_inputs


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    out=root/'HNSCC_prior_receptor_qualification_20261016_v1';path=out/'private_independent_audit.json';assert not path.exists()
    plan=json.loads((out/'plan_before_qualification.json').read_text())
    reference_path=root/'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv'
    assert hashlib.sha256(reference_path.read_bytes()).hexdigest()==plan['curated_reference_sha256']
    reference=pd.read_csv(reference_path,sep='\t').set_index(['ligand','receptor'])
    edges=pd.read_csv(out/'private_all_shortlist_receptor_edges.tsv',sep='\t')
    summaries=pd.read_csv(out/'private_all_frozen_ligand_qualification.tsv',sep='\t')
    g,source_rows,cells,var,eligible,measured=hnscc_inputs(root,project)
    original=pd.read_csv(root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1/private_patient_subtype_scores.tsv',sep='\t')
    patients=set(original.loc[original.subtype_candidate.eq('Macro_C1QC'),'patient_id'])
    mask=cells.patient_id.isin(patients)&cells.timepoint.eq('pre')&cells.library_id.isin(eligible.library_id)&cells.major.eq('Myeloid')&cells.subtype.eq('Macro_C1QC')&cells.evidence.eq('multigene_subtype_candidate')
    positions=np.flatnonzero(mask.to_numpy());receiver=cells.iloc[positions]
    aliases=defaultdict(set)
    for i,(name,row) in enumerate(var.iterrows()):
        if not measured[i] or str(name).startswith('AMBIGUOUS::'):continue
        names={str(name),str(row.get('gene_symbol',''))}
        if row.get('feature_identity_status')=='stable_Ensembl_ID':names|=set(str(row.get('source_gene_aliases','')).split('|'))
        for name in names:
            if name not in ['','nan','None']:aliases[name].add(i)
    columns={}
    for row in edges.itertuples():
        genes=row.receptor.split('_');candidate=[next(iter(aliases[g])) for g in genes] if all(len(aliases[g])==1 for g in genes) else None
        columns[row.ligand,row.receptor]=candidate if candidate is not None and len(candidate)==len(set(candidate)) else None
    needed=sorted({c for values in columns.values() if values is not None for c in values});positive_cells={c:set() for c in needed}
    matrix=ad.io.sparse_dataset(g['layers/counts'])
    for start in range(0,len(positions),1000):
        raw=matrix[source_rows[positions[start:start+1000]],:].tocsr()[:,needed].tocoo()
        assert np.isfinite(raw.data).all() and (raw.data>=0).all()
        for local_col,physical in enumerate(needed):
            positive_cells[physical].update((raw.row[(raw.col==local_col)&(raw.data>0)]+start).tolist())
    patient_names=receiver.patient_id.to_numpy(str);patient_sizes=receiver.groupby('patient_id').size()
    coverage=[]
    accepted=set(plan['accepted_reference_status'])
    for row in edges.itertuples():
        exact=reference.loc[(row.ligand,row.receptor)];assert exact.reference_status==row.reference_status
        values=columns[row.ligand,row.receptor]
        assert (values is not None)==bool(row.receptor_components_measured)
        joint=set.intersection(*(positive_cells[c] for c in values)) if values is not None else None
        if joint is None:
            assert pd.isna(row.same_cell_receptor_positive_cells);valid=False
            any_patients=onepct_patients=None
        else:
            assert len(joint)==row.same_cell_receptor_positive_cells
            assert np.isclose(len(joint)/len(positions),row.same_cell_receptor_positive_fraction,atol=1e-12)
            valid=len(joint)*100>=len(positions)
            counts=pd.Series(patient_names[sorted(joint)]).value_counts().reindex(patient_sizes.index,fill_value=0)
            any_patients=int(counts.gt(0).sum());onepct_patients=int((counts*100>=patient_sizes).sum())
        assert valid==bool(row.joint_RNA_gate)
        assert (valid and row.reference_status in accepted)==bool(row.prior_reference_and_RNA_qualified)
        coverage.append(dict(ligand=row.ligand,receptor=row.receptor,receiver_patients=len(patient_sizes),
            patients_with_any_same_cell_receptor_RNA=any_patients,patients_with_atleast1pct_same_cell_receptor_RNA=onepct_patients,
            original_qualification_unchanged=True))
    for row in summaries.itertuples():
        subset=edges[edges.ligand.eq(row.ligand)]
        good=subset[subset.prior_reference_and_RNA_qualified]
        assert len(subset)==row.literal_reference_edges and len(good)==row.reference_and_joint_RNA_qualified_edges
        assert bool(len(good))==bool(row.prospective_receiver_qualified)
    pd.DataFrame(coverage).to_csv(out/'private_all_edge_patient_RNA_coverage.tsv',sep='\t',index=False,na_rep='NA')
    report=dict(status='PASS_ALL149_CURATED_EDGES_SPARSE_RECEPTOR_CELLSET_INTERSECTIONS_AND_PATIENT_COVERAGE',
        finished_UTC=datetime.now(timezone.utc).isoformat(),literal_edges=len(edges),unique_ligands=len(summaries),
        qualified_table_sha256=hashlib.sha256((out/'private_all_frozen_ligand_qualification.tsv').read_bytes()).hexdigest(),
        receptor_RNA_not_protein_assembly=True,qualification_or_statistical_family_changed=False)
    path.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()

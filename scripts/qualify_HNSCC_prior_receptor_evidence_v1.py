"""Cross-check frozen ligand shortlist with curated edges and same-cell RNA.

No replacement of failed shortlist entries with lower-ranked ligands. No
sender-response effect tests, direct binding or protein-complex claims.
"""
import hashlib,json,os
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
import anndata as ad
import numpy as np,pandas as pd
from diagnose_frozen_MCP_two_end_functions_v1 import hnscc_inputs


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    prior=root/'HNSCC_receiver_prior_ligand_ranking_20261016_v1'
    audit=json.loads((prior/'private_independent_audit.json').read_text());assert audit['status'].startswith('PASS_NATIVE_RECEIVER')
    rank_path=prior/'private_all_ligand_prior_rankings.tsv'
    assert hashlib.sha256(rank_path.read_bytes()).hexdigest()==audit['prior_ranking_sha256']
    ranks=pd.read_csv(rank_path,sep='\t');selected=ranks[ranks.prospective_top20]
    ligands=sorted(set(selected.ligand));assert len(selected)==40
    reference_path=root/'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv'
    reference=pd.read_csv(reference_path,sep='\t',keep_default_na=False);assert len(reference)==4624
    out=root/'HNSCC_prior_receptor_qualification_20261016_v1';assert not out.exists();out.mkdir()
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        frozen_shortlist_sha256=audit['prior_ranking_sha256'],curated_reference_sha256=hashlib.sha256(reference_path.read_bytes()).hexdigest(),
        source='Pinned all4624 resource cross-reference, not new curated evidence or a receptor-specific NicheNet target prior.',
        accepted_reference_status=['simplex_reference_supported','all_constituents_supported_complex_assembly_unverified'],
        exact_ligand_matching='Literal whole ligand identifier only; component-only ligand hits are not promoted to assembled ligand.',
        receiver_RNA='All receptor components uniquely measured; >=1% of original full native C1QC baseline cells detect all components in the same cell.',
        binding='Constituent references and RNA compatibility only; assembled complex/protein binding/current-tumour signaling not verified.',
        failed_entries='Retain complete original40 sign slots, no replacement with rank21+ and no significance tests.')
    (out/'plan_before_qualification.json').write_text(json.dumps(plan,indent=2)+'\n')
    g,source_rows,cells,var,eligible,measured=hnscc_inputs(root,project)
    sub=pd.read_csv(root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1/private_patient_subtype_scores.tsv',sep='\t')
    patients=set(sub.loc[sub.subtype_candidate.eq('Macro_C1QC'),'patient_id'])
    mask=cells.timepoint.eq('pre')&cells.patient_id.isin(patients)&cells.library_id.isin(eligible.library_id)&cells.major.eq('Myeloid')&cells.subtype.eq('Macro_C1QC')&cells.evidence.eq('multigene_subtype_candidate')
    positions=np.flatnonzero(mask.to_numpy());assert len(positions)>0
    lookup=defaultdict(set)
    for c,(name,row) in enumerate(var.iterrows()):
        if not measured[c] or str(name).startswith('AMBIGUOUS::'):continue
        aliases={str(name),str(row.get('gene_symbol',''))}
        if row.get('feature_identity_status')=='stable_Ensembl_ID':aliases.update(str(row.get('source_gene_aliases','')).split('|'))
        for alias in aliases:
            if alias not in ['','nan','None']:lookup[alias].add(c)
    edges=reference[reference.ligand.isin(ligands)].copy().sort_values(['ligand','receptor'])
    assert not edges.duplicated(['ligand','receptor']).any()
    edge_columns=[];all_columns=set()
    for edge in edges.itertuples():
        genes=str(edge.receptor).split('_')
        columns=[next(iter(lookup[gene])) for gene in genes] if all(len(lookup.get(gene,set()))==1 for gene in genes) else None
        if columns is not None and len(set(columns))!=len(columns):columns=None
        edge_columns.append(columns)
        if columns is not None:all_columns.update(columns)
    all_columns=sorted(all_columns);local={c:i for i,c in enumerate(all_columns)}
    counts=np.zeros(len(edges),np.int64);single=np.zeros(len(all_columns),np.int64)
    matrix=ad.io.sparse_dataset(g['layers/counts'])
    for start in range(0,len(positions),2000):
        raw=matrix[source_rows[positions[start:start+2000]],:].tocsr();raw.eliminate_zeros()
        positive=raw[:,all_columns].toarray()>0;single+=positive.sum(axis=0)
        for i,columns in enumerate(edge_columns):
            if columns is not None:counts[i]+=np.all(positive[:,[local[c] for c in columns]],axis=1).sum()
    rows=[]
    accepted=set(plan['accepted_reference_status'])
    for i,edge in enumerate(edges.itertuples()):
        columns=edge_columns[i];joint=int(counts[i]) if columns is not None else None
        native_ok=columns is not None and joint*100>=len(positions)
        rows.append(dict(ligand=edge.ligand,receptor=edge.receptor,reference_status=edge.reference_status,
            n_supported=int(edge.n_supported),n_curators_rejected=int(edge.n_curators_rejected),
            receptor_components_measured=columns is not None,native_C1QC_cells=len(positions),
            same_cell_receptor_positive_cells=joint,same_cell_receptor_positive_fraction=joint/len(positions) if joint is not None else None,
            joint_RNA_gate=native_ok,prior_reference_and_RNA_qualified=edge.reference_status in accepted and native_ok,
            assembled_complex_verified=False,current_tumour_signaling_verified=False))
    details=pd.DataFrame(rows)
    summaries=[]
    for ligand in ligands:
        group=details[details.ligand.eq(ligand)] if len(details) else pd.DataFrame()
        good=group[group.prior_reference_and_RNA_qualified] if len(group) else pd.DataFrame()
        signs=sorted(set(selected.loc[selected.ligand.eq(ligand),'programme_sign']))
        summaries.append(dict(ligand=ligand,programme_signs='|'.join(signs),literal_reference_edges=len(group),
            reference_and_joint_RNA_qualified_edges=len(good),
            prospective_receiver_qualified=bool(len(good)),qualified_receptors='|'.join(good.receptor.astype(str)) if len(good) else '',
            no_lower_rank_replacement=True,protein_or_causal_axis_verified=False))
    details.to_csv(out/'private_all_shortlist_receptor_edges.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(summaries).to_csv(out/'private_all_frozen_ligand_qualification.tsv',sep='\t',index=False)
    report=dict(status='ALL_FROZEN_PRIOR_LIGANDS_CURATED_EDGE_AND_NATIVE_JOINT_RNA_CHECKED_NO_EFFECTS',
        finished_UTC=datetime.now(timezone.utc).isoformat(),original_sign_slots=len(selected),unique_ligands=len(ligands),
        literal_reference_edges=len(details),ligands_with_qualified_constituent_reference_and_RNA=sum(r['prospective_receiver_qualified'] for r in summaries),
        not_direct_binding_or_complex_assembly_validation=True,no_new_P_BH=True,
        qualified_table_sha256=hashlib.sha256((out/'private_all_frozen_ligand_qualification.tsv').read_bytes()).hexdigest())
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    print(pd.DataFrame(summaries).to_json(orient='records'))


if __name__=='__main__':main()

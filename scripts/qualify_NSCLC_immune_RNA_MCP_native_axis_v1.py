"""Outcome-free full native-axis identity qualification for the frozen RNA MCP."""
import json
import os
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata.io import read_elem

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    paths=dict(native=root/'GSE243013_genes.csv.gz',
        weights=root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv',
        mapping=root/'immune_RNA_MCP_discovery_scope_20261018_v1/private_marker_excluded_native_gene_map.tsv',
        units=root/'NSCLC_full_native_gene_units_20261016_v1/private_all_native_units_and_quality.tsv')
    expected=dict(native='ac325c9ce6bfaf2270eb931e7d947bf31d1fb151b14cb166ed9f6c8f4b34a104',
        weights='fd7316c67beb687c403a752c36b5989644e432f59a83119b21eade5efb0d89a7',
        mapping='17a214a298e743ca205bd57f7419e6b4484afd5e67d9141711dc99bd7d5ff50c',
        units='a032ef126f10ce20a327eefbddad43ea9fcb565cdb6b1587ef7917846d1728f9')
    for key,p in paths.items():assert sha(p)==expected[key]
    out=root/'NSCLC_immune_RNA_MCP_native_axis_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,input_sha256=expected,
        scope='All6 fixed K2 source-adjusted view/programme active gene axes and fixed34 capture markers; all original native units. No expression values, clinical outcomes or effects read.',
        identity='Direct unique uppercase native symbol first. Otherwise require target physical columns in all3 frozen sources share a stable Ensembl identity, unique authority alias identity, exactly one native column. Reject all within-programme physical-column collisions.',
        gate='All active genes AND all fixed capture markers measured for exact programme/view projection. Coverage fraction alone cannot qualify missing genes; no missing zero fill.',
        limits='Native deposition has symbols only, not stable IDs. Supported alias inference uses existing frozen authority metadata, not original native stable-ID proof. It is reported separately. No fuzzy/manual matching or outcome selection.',
        no_new_download=True,no_clinical_or_programme_association=True)
    (out/'plan_before_all_native_identity_qualification.json').write_text(json.dumps(plan,indent=2)+'\n')
    native=pd.read_csv(paths['native'],keep_default_na=False).geneSymbol.str.strip().str.upper().to_numpy(str)
    assert len(native)==31831 and len(set(native))==len(native)
    index={name:i for i,name in enumerate(native)}
    weights=pd.read_csv(paths['weights'],sep='\t');weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')&weights.weight.ne(0)].copy()
    gm=pd.read_csv(paths['mapping'],sep='\t').set_index('gene');assert gm.index.is_unique
    z=zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),mode='r')
    symbol=np.asarray(read_elem(z['var/gene_symbol'])).astype(str)
    identity=np.asarray(read_elem(z['var/feature_identity_status'])).astype(str)
    canonical=np.asarray(read_elem(z['var/canonical_gene_identity'])).astype(str)
    aliases=np.asarray(read_elem(z['var/source_gene_aliases'])).astype(str)
    alias_ids=defaultdict(set);id_names=defaultdict(set)
    for i in np.flatnonzero(identity=='stable_Ensembl_ID'):
        if not canonical[i].startswith('ENSG'):continue
        names={symbol[i].upper(),canonical[i].upper()}|set(a.upper() for a in aliases[i].split('|') if a)
        for name in names:alias_ids[name].add(canonical[i]);id_names[canonical[i]].add(name)
    records=[]
    for r in weights.itertuples(index=False):
        gene=str(r.gene);col=index.get(gene.upper(),-1);status='DIRECT_UNIQUE_NATIVE_SYMBOL' if col>=0 else 'UNMEASURED_NO_QUALIFIED_NATIVE_ALIAS'
        target_ids=set();proof=False
        physical=[int(gm.loc[gene,s]) for s in ['GSE286827','GSE294109','GSE305078']]
        if all(c>=0 and identity[c]=='stable_Ensembl_ID' and canonical[c].startswith('ENSG') for c in physical):
            target_ids={canonical[c] for c in physical};proof=len(target_ids)==1
        candidates=[]
        if col<0 and proof:
            target=next(iter(target_ids))
            candidates=sorted({index[name] for name in id_names[target] if name in index and alias_ids[name]=={target}})
            if len(candidates)==1:col=candidates[0];status='SUPPORTED_UNIQUE_STABLE_ID_ALIAS'
            elif len(candidates)>1:status='AMBIGUOUS_MULTIPLE_NATIVE_ALIAS_COLUMNS'
        records.append(dict(view=r.view,program=r.program,gene=gene,weight=float(r.weight),native_column=col,status=status,
            same_stable_ID_all3=proof,target_stable_ID='|'.join(sorted(target_ids)),native_alias_candidates=len(candidates)))
    result=pd.DataFrame(records)
    for _,g in result.groupby(['view','program']):
        duplicated=g[g.native_column.ge(0)].native_column.duplicated(keep=False)
        where=g[g.native_column.ge(0)].index[duplicated]
        result.loc[where,'status']='WITHIN_PROGRAMME_NATIVE_COLUMN_COLLISION';result.loc[where,'native_column']=-1
    fixed=root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    panels=json.loads(fixed.read_text())['panels'];markers=sorted(set(g for genes in panels.values() for g in genes));assert len(markers)==34
    marker=pd.DataFrame(dict(gene=markers,native_column=[index.get(g.upper(),-1) for g in markers]))
    summary=[]
    for (view,program),g in result.groupby(['view','program'],sort=True):
        valid=g.native_column.ge(0);total=float(g.weight.abs().sum())
        summary.append(dict(dataset='GSE243013',view=view,program=program,active_genes=len(g),
            direct=int(g.status.eq('DIRECT_UNIQUE_NATIVE_SYMBOL').sum()),supported_alias=int(g.status.eq('SUPPORTED_UNIQUE_STABLE_ID_ALIAS').sum()),
            missing_or_ambiguous=int((~valid).sum()),absolute_weight_retained=float(g.loc[valid,'weight'].abs().sum()/total),
            fixed_capture_markers=len(marker),fixed_capture_measured=int(marker.native_column.ge(0).sum()),
            full_gene_axis_qualified=bool(valid.all() and marker.native_column.ge(0).all())))
    result.to_csv(out/'private_all_active_gene_native_axis.tsv',sep='\t',index=False)
    marker.to_csv(out/'private_fixed34_capture_gene_axis.tsv',sep='\t',index=False)
    pd.DataFrame(summary).to_csv(out/'private_all6_native_axis_qualification.tsv',sep='\t',index=False)
    units=pd.read_csv(paths['units'],sep='\t',keep_default_na=False)
    assert len(units)==11701 and units.unit_index.tolist()==list(range(len(units)))
    author=[]
    for major,g in units.groupby('major_cell_type'):
        author.append(dict(dataset='GSE243013',author_major=str(major),native_cells=int(g.n_cells.sum()),native_materials=int(g.sampleID.nunique()),author_subtypes=int(g.sub_cell_type.nunique())))
    pd.DataFrame(author).to_csv(out/'private_all_native_author_major_coverage.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()}
    assert sum(v['bytes'] for v in files.values())<4*1024**2
    terminal=dict(status='COMPLETE_ALL6_NATIVE_IDENTITY_AND_FIXED_CAPTURE_AXIS_INDEPENDENT_ACCEPTANCE_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,all_native_cells=int(units.n_cells.sum()),
        all_native_materials=int(units.sampleID.nunique()),no_expression_or_outcome_values_read=True)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

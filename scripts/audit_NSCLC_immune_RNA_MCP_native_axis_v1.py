"""Independently replay every fixed RNA gene/native-axis mapping and coverage."""
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
    source=root/'NSCLC_immune_RNA_MCP_native_axis_20261019_v1'
    assert sha(source/'private_terminal.json')=='7f1ddaa54cfb2fd550c4892ce6e63363888641103bf57f6a986f974314a33e5d'
    terminal=json.loads((source/'private_terminal.json').read_text())
    for name,item in terminal['files'].items():assert sha(source/name)==item['sha256']
    out=root/'NSCLC_immune_RNA_MCP_native_axis_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_full_axis_replay.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,
        scope='All active learned gene positions, fixed34 markers, all6 source coverage descriptions and all original native-author major counts. Independent native-name scans and stable-ID alias relation reconstruction; no score or outcome read.'),indent=2)+'\n')
    names=pd.read_csv(root/'GSE243013_genes.csv.gz',keep_default_na=False).geneSymbol.str.strip().str.upper().to_numpy(str)
    assert len(names)==len(set(names))==31831
    gm=pd.read_csv(root/'immune_RNA_MCP_discovery_scope_20261018_v1/private_marker_excluded_native_gene_map.tsv',sep='\t').set_index('gene')
    weights=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv',sep='\t')
    weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')&weights.weight.ne(0)]
    saved=pd.read_csv(source/'private_all_active_gene_native_axis.tsv',sep='\t',keep_default_na=False)
    assert not saved.duplicated(['view','program','gene']).any()
    assert set(map(tuple,saved[['view','program','gene']].to_numpy()))==set(map(tuple,weights[['view','program','gene']].to_numpy()))
    z=zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),mode='r')
    canonical=np.asarray(read_elem(z['var/canonical_gene_identity'])).astype(str)
    identity=np.asarray(read_elem(z['var/feature_identity_status'])).astype(str)
    symbol=np.asarray(read_elem(z['var/gene_symbol'])).astype(str)
    aliases=np.asarray(read_elem(z['var/source_gene_aliases'])).astype(str)
    by_name=defaultdict(set)
    for i,kind in enumerate(identity):
        if kind!='stable_Ensembl_ID' or not canonical[i].startswith('ENSG'):continue
        labels=[symbol[i],canonical[i],*aliases[i].split('|')]
        for label in labels:
            if label:by_name[label.upper()].add(canonical[i])
    replay=[]
    for item in saved.itertuples(index=False):
        direct=np.flatnonzero(names==item.gene.upper());assert len(direct)<=1
        physical=gm.loc[item.gene,['GSE286827','GSE294109','GSE305078']].to_numpy(int)
        stable=all(i>=0 and identity[i]=='stable_Ensembl_ID' and canonical[i].startswith('ENSG') for i in physical)
        ids={canonical[i] for i in physical} if stable else set()
        proof=stable and len(ids)==1
        assert bool(item.same_stable_ID_all3)==proof and item.target_stable_ID=='|'.join(sorted(ids))
        candidates=[]
        if len(direct):col=int(direct[0]);status='DIRECT_UNIQUE_NATIVE_SYMBOL'
        else:
            col=-1;status='UNMEASURED_NO_QUALIFIED_NATIVE_ALIAS'
            if proof:
                target=next(iter(ids));candidates=[i for i,name in enumerate(names) if by_name.get(name,set())=={target}]
                if len(candidates)==1:col=candidates[0];status='SUPPORTED_UNIQUE_STABLE_ID_ALIAS'
                elif len(candidates)>1:status='AMBIGUOUS_MULTIPLE_NATIVE_ALIAS_COLUMNS'
        assert item.native_alias_candidates==len(candidates)
        replay.append(dict(view=item.view,program=item.program,gene=item.gene,native_column=col,status=status))
    replay=pd.DataFrame(replay)
    for _,g in replay.groupby(['view','program']):
        values=g[g.native_column.ge(0)];freq=values.native_column.value_counts();collision=set(freq[freq>1].index)
        take=g.index[g.native_column.isin(collision)];replay.loc[take,'native_column']=-1;replay.loc[take,'status']='WITHIN_PROGRAMME_NATIVE_COLUMN_COLLISION'
    assert replay.native_column.tolist()==saved.native_column.tolist() and replay.status.tolist()==saved.status.tolist()
    fixed=root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    genes=sorted({g for p in json.loads(fixed.read_text())['panels'].values() for g in p});assert len(genes)==34
    markers=pd.read_csv(source/'private_fixed34_capture_gene_axis.tsv',sep='\t');assert markers.gene.tolist()==genes
    marker_cols=[]
    for gene in genes:
        found=np.flatnonzero(names==gene.upper());assert len(found)<=1;marker_cols.append(int(found[0]) if len(found) else -1)
    assert marker_cols==markers.native_column.tolist()
    summary=pd.read_csv(source/'private_all6_native_axis_qualification.tsv',sep='\t');assert len(summary)==6
    for row in summary.itertuples(index=False):
        g=saved[saved.view.eq(row.view)&saved.program.eq(row.program)];w=weights[weights.view.eq(row.view)&weights.program.eq(row.program)].set_index('gene').weight.reindex(g.gene)
        assert np.allclose(w.to_numpy(),g.weight.to_numpy(),rtol=0,atol=1e-15)
        valid=g.native_column.ge(0);fraction=np.sum(np.abs(w.to_numpy()[valid]))/np.sum(np.abs(w.to_numpy()))
        assert row.active_genes==len(g) and row.direct==g.status.eq('DIRECT_UNIQUE_NATIVE_SYMBOL').sum()
        assert row.supported_alias==g.status.eq('SUPPORTED_UNIQUE_STABLE_ID_ALIAS').sum() and row.missing_or_ambiguous==(~valid).sum()
        assert abs(row.absolute_weight_retained-fraction)<1e-12 and row.fixed_capture_measured==sum(c>=0 for c in marker_cols)
        assert bool(row.full_gene_axis_qualified)==bool(valid.all() and min(marker_cols)>=0)
    units=pd.read_csv(root/'NSCLC_full_native_gene_units_20261016_v1/private_all_native_units_and_quality.tsv',sep='\t',keep_default_na=False)
    coverage=pd.read_csv(source/'private_all_native_author_major_coverage.tsv',sep='\t',keep_default_na=False)
    assert set(coverage.author_major)==set(units.major_cell_type)
    for row in coverage.itertuples(index=False):
        g=units[units.major_cell_type.eq(row.author_major)]
        assert row.native_cells==sum(map(int,g.n_cells)) and row.native_materials==len(set(g.sampleID)) and row.author_subtypes==len(set(g.sub_cell_type))
    assert terminal['all_native_cells']==sum(map(int,units.n_cells))==1254749 and terminal['all_native_materials']==len(set(units.sampleID))==243
    t=dict(status='ACCEPTED_ALL_ACTIVE_NATIVE_AXIS_IDENTITIES_FIXED34_MARKERS_AND_ALL_SOURCE_COVERAGE',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),
        active_gene_positions_checked=len(saved),marker_positions_checked=34,view_programmes_checked=6,all_native_units=11701,
        expression_values_or_clinical_outcomes_read=False)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n')
    print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

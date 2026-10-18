"""Replay external panel measurement flags against accepted native gene masks.

Only metadata and plans. No expression or effect values, no P or downloads.
"""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from describe_NSCLC_native_interface_mechanistic_followup_v1 import read_scope, SCOPE_SHA
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);resource(root,1024**2)
    out=root/'NSCLC_candidate_external_exact_RNA_preflight_20261018_v1'
    terminal=json.loads((out/'private_terminal.json').read_text())
    assert terminal['source_sha256']=='ba0c74250f765468447b9a95baf025860e0a057070e285aa111d72b7abb9acb1'
    assert terminal['scope_sha256']==SCOPE_SHA
    for name,entry in terminal['files'].items():assert sha(out/name)==entry['sha256']
    cal=root/'NSCLC_native_MFA_external_RNA_universe_calibration_20261016_v1'
    proof=cal/'private_independent_native_RNA_universe_audit.json'
    assert sha(proof)=='f9822656c7be74fa7953b04c8807b5da718bc2084bf3a946e7943504ed31d5b3'
    assert json.loads(proof.read_text())['producer_terminal_sha256']==sha(cal/'private_terminal.json')
    manifest=json.loads((cal/'private_terminal.json').read_text())
    for name in ['private_native_universe_masks.tsv.gz','private_all61_source_native_universe_map.tsv']:
        assert sha(cal/name)==manifest['files'][name]['sha256']
    masks=pd.read_csv(cal/'private_native_universe_masks.tsv.gz',sep='\t',index_col=0)
    assert masks.isin([0,1]).all().all();masks=masks.astype(bool)
    mapping=pd.read_csv(cal/'private_all61_source_native_universe_map.tsv',sep='\t')
    scope,candidates=read_scope(dict(root=root))
    desired=set(scope['ligand_specs'])|{g for s in scope['ligand_specs'].values() for c in s['complexes'] for g in c}
    assert desired<=set(masks.index)
    gene=pd.read_csv(out/'private_all_requested_gene_mappings.tsv',sep='\t')
    route=pd.read_csv(out/'private_all_candidate_source_routes.tsv',sep='\t')
    sources=pd.read_csv(out/'private_all61_exact_RNA_source_routes.tsv',sep='\t')
    assert len(gene)==61*len(desired) and len(route)==61*len(candidates) and len(sources)==61
    assert not gene.duplicated(['arm','dataset','gene']).any() and not route.duplicated(['record','arm','dataset']).any()
    for item in mapping.itertuples():
        g=gene[gene.arm.eq(item.arm)&gene.dataset.eq(item.dataset)]
        assert set(g.gene)==desired
        available={x.gene:bool(masks.loc[x.gene,'mask_%03d'%item.mask_index]) for x in g.itertuples()}
        assert all((x.source_column>=0)==available[x.gene] for x in g.itertuples())
        present=g[g.source_column.ge(0)];assert present.source_column.is_unique
        for c in candidates.itertuples():
            row=route[route.record.eq(c.record)&route.arm.eq(item.arm)&route.dataset.eq(item.dataset)]
            assert len(row)==1;row=row.iloc[0]
            complete=sum(all(available[g] for g in rec) for rec in scope['ligand_specs'][c.ligand]['complexes'])
            assert bool(row.ligand_unique_measured)==available[c.ligand]
            assert row.complete_functional_complexes_measured==complete
            assert bool(row.measurement_opportunity)==bool(available[c.ligand] and complete)
    result=dict(status='PASS_ALL_SOURCE_GENE_FLAGS_AND_CANONICAL_COMPLEX_ROUTES_AGAINST_ACCEPTED_NATIVE_MASKS',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        producer_terminal_sha256=sha(out/'private_terminal.json'),scope_sha256=SCOPE_SHA,
        limits='All requested gene flags and candidate/source routes checked; readable container and cache availability are producer metadata declarations. No new count values, expression, full raw audit, pool equivalence or model validity established.')
    path=out/'private_independent_exact_RNA_metadata_audit.json';assert not path.exists()
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],completed_UTC=result['completed_UTC'],proof_sha256=sha(path))))


if __name__=='__main__':main()

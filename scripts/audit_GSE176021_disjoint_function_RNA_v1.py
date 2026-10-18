"""Replay every native selected cell and every official score using scalar ranks."""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import zarr
from scipy.sparse import csr_matrix
from scipy.stats import rankdata
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from immune_RNA_MCP_new_source_io_v1 import wait


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 2*1024**2)
    out = root/'GSE176021_disjoint_function_native_RNA_audit_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    src = root/'GSE176021_disjoint_function_native_RNA_20261019_v1'
    terminal = wait(root,out,'measure_GSE176021_disjoint_function_RNA_v1',src/'private_terminal.json',{'COMPLETE_ALL50_NATIVE_FUNCTIONAL_SCORES_FULL_REPLAY_PENDING'})
    for name, info in terminal['files'].items():
        assert sha(src/name) == info['sha256']
    assert sha(root/'GSE176021_disjoint_function_reference_audit_20261019_v1/private_terminal.json') == 'ee8ab9e6cd44c0039f42da978d623a02ce2ac17cee62033851e3228f99444976'
    scope = root/'GSE176021_author_lineage_RNA_scope_20261019_v2'
    inp = root/'GSE176021_paired_tissue_inputs_20261019_v1'
    cases = np.load(inp/'private_AP01_inputs.npz')['cases']
    pairs = pd.read_csv(scope/'private_all_scope_material_pair_positions.tsv',sep='\t')
    pairs = pairs[(pairs.scope == 'all_author_lineage') & pairs.case_block.isin(cases) & pairs.tissue.isin(['tumor','adjacent normal lung']) & pairs.both_views_eligible]
    assert len(pairs) == 60
    pools = pd.read_csv(scope/'private_all_eligible_pools.tsv',sep='\t')
    pools = pools[(pools.scope == 'all_author_lineage') & pools.sampleID.isin(pairs.sampleID)]
    assert len(pools) == 120
    old = dict(np.load(root/'GSE176021_author_lineage_native_RNA_20261019_v1/private_selected_cell_quality_scores.npz'))
    ix = np.flatnonzero(np.isin(old['pool_index'],pools.pool_index))
    ix = ix[np.argsort(old['global_row'][ix])]
    expected_rows = old['global_row'][ix]
    result = dict(np.load(src/'private_all_selected_cell_functional_scores.npz'))
    assert len(ix) == 11939 and np.array_equal(result['global_row'],expected_rows) and np.array_equal(result['pool_index'],old['pool_index'][ix])
    refroot = root/'GSE176021_disjoint_function_reference_scope_20261019_v1'
    definitions = json.loads((refroot/'private_all50_full_and_disjoint_reference_gene_definitions.json').read_text())
    universe = np.load(refroot/'private_all_source_measured_rank_universe.npz')['native_columns']
    references = [d['reference'] for d in definitions]
    assert np.array_equal(result['references'],references)
    ids = [np.array(list(d['unique_source_mapping'].values()),int) for d in definitions]
    assert all(len(set(a)) == len(a) for a in ids)
    assert result['ucell'].shape == (11939,50,3) and result['mean_logCP10K'].shape == (11939,50)
    z = zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),mode='r')
    layer = z['layers/counts']
    nf = int(layer.attrs['shape'][1])
    measured = np.zeros(nf,bool);measured[universe] = True
    digest = {k:hashlib.sha256() for k in ['global_row_int64','native_indices_int64','native_data_int64','native_row_endpoints_int64']}
    offset = 0
    max_uc, max_ml = 0., 0.
    for block in np.unique(expected_rows//2048):
        at = np.flatnonzero(expected_rows//2048 == block)
        rr = expected_rows[at].astype(np.int64)
        lo, hi = int(block)*2048, min(int(layer.attrs['shape'][0]),(int(block)+1)*2048)
        ptr = np.asarray(layer['indptr'][lo:hi+1],np.int64)
        matrix = csr_matrix((np.asarray(layer['data'][int(ptr[0]):int(ptr[-1])]),np.asarray(layer['indices'][int(ptr[0]):int(ptr[-1])]),ptr-ptr[0]),shape=(hi-lo,nf))[rr-lo].copy()
        matrix.sum_duplicates();matrix.eliminate_zeros();matrix.sort_indices();matrix = matrix.astype(np.int64)
        assert np.all(matrix.data > 0) and measured[matrix.indices].all()
        for key, array in [('global_row_int64',rr),('native_indices_int64',matrix.indices),('native_data_int64',matrix.data),('native_row_endpoints_int64',matrix.indptr[1:]+offset)]:
            digest[key].update(np.asarray(array,dtype='<i8').tobytes())
        offset += matrix.nnz
        for j, pos in enumerate(at):
            start, stop = matrix.indptr[j:j+2]
            columns, values = matrix.indices[start:stop],matrix.data[start:stop]
            total = int(sum(values));detected = len(values)
            assert total == old['fullRNA'][ix[pos]] and detected == old['detected'][ix[pos]]
            ranks = rankdata(-values,method='average')
            gene_log = np.zeros(nf,float)
            gene_log[columns] = np.log1p(values*10000./total)
            log_expected = np.array([np.sum(gene_log[g])/len(g) for g in ids])
            error = float(np.max(np.abs(log_expected-result['mean_logCP10K'][pos])))
            max_ml = max(max_ml,error);assert error < 1e-10
            for k, rmax in enumerate([1500,1000,2000]):
                grank = np.full(nf,float(rmax));grank[columns] = np.minimum(ranks,rmax)
                scores = []
                for genes in ids:
                    n = len(genes);minimum = n*(n+1)/2
                    scores.append(1-(sum(grank[genes])-minimum)/(n*rmax-minimum))
                error = float(np.max(np.abs(np.array(scores)-result['ucell'][pos,:,k])))
                max_uc = max(max_uc,error);assert error < 1e-7
        (out/'private_progress.json').write_text(json.dumps({'phase':'FULL_NATIVE_SCALAR_REPLAY','selected_cells_replayed':int(at[-1])+1,'selected_cells':len(expected_rows),'observed_UTC':datetime.now(timezone.utc).isoformat()})+'\n')
    official = json.loads((src/'private_native_canonical_and_software_receipt.json').read_text())
    assert official['digests'] == {k:d.hexdigest() for k,d in digest.items()} and offset == official['canonical_nonzero_values']
    table = pd.read_csv(src/'private_all6000_functional_pool_values.tsv',sep='\t')
    assert len(table) == 6000 and not table.duplicated(['pool_index','reference']).any()
    pmeta = pools.set_index('pool_index')
    lookup = {name:j for j,name in enumerate(references)}
    for row in table.itertuples(index=False):
        ids_here = np.flatnonzero(result['pool_index'] == row.pool_index)
        p = pmeta.loc[row.pool_index]
        assert row.sampleID == p.sampleID and row.view == p['view'] and row.selected_cells == len(ids_here) == p.selected_cells
        r = lookup[row.reference]
        np.testing.assert_allclose([row.ucell1500,row.ucell1000,row.ucell2000,row.mean_logCP10K],[*[result['ucell'][ids_here,r,k].astype(float).sum()/len(ids_here) for k in range(3)],result['mean_logCP10K'][ids_here,r].sum()/len(ids_here)],rtol=0,atol=1e-10)
    receipt = {'status':'ACCEPTED_COMPLETE_NATIVE_ALL50_FUNCTIONAL_SCORES_AND6000_POOL_VALUES','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'producer_terminal_sha256':sha(src/'private_terminal.json'),'resource_guard':guard,'selected_cells':len(expected_rows),'canonical_nonzero_values':int(offset),'all_ucell_positions':int(result['ucell'].size),'all_mean_log_positions':int(result['mean_logCP10K'].size),'all_pool_reference_positions':len(table),'maximum_scalar_ucell_error':max_uc,'maximum_scalar_mean_log_error':max_ml,'no_functional_association_or_gene_fit':True}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__ == '__main__':
    main()

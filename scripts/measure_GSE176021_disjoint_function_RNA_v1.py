"""Frozen, full-reference native RNA assay; official UCell in the existing runtime."""
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path


def worker():
    import contextlib
    import hashlib
    import importlib.metadata
    import inspect
    import io
    import pickle
    import anndata as ad
    import numpy as np
    import pandas as pd
    import pyucell
    assert importlib.metadata.version('pyucell') == '0.8.0'
    source = Path(inspect.getfile(pyucell.compute_ucell_scores))
    config = pickle.load(sys.stdin.buffer)
    var = pd.DataFrame(index=config['features'])
    receipt = {'python':sys.version.split()[0], 'pyucell':'0.8.0', 'official_module_sha256':hashlib.sha256(source.read_bytes()).hexdigest(), 'api':str(inspect.signature(pyucell.compute_ucell_scores))}
    pickle.dump(receipt, sys.stdout.buffer, protocol=4)
    sys.stdout.buffer.flush()
    while True:
        x = pickle.load(sys.stdin.buffer)
        if x is None:
            break
        scores = []
        for rank in [1500, 1000, 2000]:
            a = ad.AnnData(X=x.copy(), var=var.copy())
            with contextlib.redirect_stdout(io.StringIO()):
                pyucell.compute_ucell_scores(a, signatures=config['signatures'], max_rank=rank, ties_method='average', missing_genes='skip', chunk_size=64, n_jobs=1, device='cpu')
            scores.append(a.obs[[name+'_UCell' for name in config['references']]].to_numpy(np.float32))
        pickle.dump(np.stack(scores, axis=2), sys.stdout.buffer, protocol=4)
        sys.stdout.buffer.flush()


def main():
    import hashlib
    import pickle
    import subprocess
    import numpy as np
    import pandas as pd
    import zarr
    from scipy.sparse import csr_matrix
    from NSCLC_native_MFA_interpretation_inputs_v1 import sha
    from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 32*1024**2)
    required = {
        'GSE176021_author_lineage_native_RNA_audit_20261019_v1':'6950e1c15fda41b93962cfd763348898802d8875f1eef6719ec395b25392730d',
        'GSE176021_disjoint_function_reference_audit_20261019_v1':'ee8ab9e6cd44c0039f42da978d623a02ce2ac17cee62033851e3228f99444976',
        'GSE176021_author_subtype_dynamics_audit_20261019_v1':'c4a8be3debb2e22b5b9282b34410d5567943cf918c8e8e7eb6498a72925b7fe6'}
    for name, digest in required.items():
        assert sha(root/name/'private_terminal.json') == digest
    runtime = root/'pyucell_runtime_20261014_v1/venv/bin/python'
    assert runtime.exists()
    scope = root/'GSE176021_author_lineage_RNA_scope_20261019_v2'
    inp = root/'GSE176021_paired_tissue_inputs_20261019_v1'
    cases = np.load(inp/'private_AP01_inputs.npz')['cases']
    pair = pd.read_csv(scope/'private_all_scope_material_pair_positions.tsv', sep='\t')
    pair = pair[pair.scope.eq('all_author_lineage') & pair.case_block.isin(cases) & pair.tissue.isin(['tumor','adjacent normal lung']) & pair.both_views_eligible]
    assert len(pair) == 60
    pool = pd.read_csv(scope/'private_all_eligible_pools.tsv', sep='\t')
    pool = pool[pool.scope.eq('all_author_lineage') & pool.sampleID.isin(pair.sampleID)]
    assert len(pool) == 120
    old = dict(np.load(root/'GSE176021_author_lineage_native_RNA_20261019_v1/private_selected_cell_quality_scores.npz'))
    take = np.flatnonzero(np.isin(old['pool_index'], pool.pool_index))
    take = take[np.argsort(old['global_row'][take])]
    rows, pids = old['global_row'][take], old['pool_index'][take]
    assert len(rows) == len(set(rows)) == 11939
    refs = root/'GSE176021_disjoint_function_reference_scope_20261019_v1'
    definitions = json.loads((refs/'private_all50_full_and_disjoint_reference_gene_definitions.json').read_text())
    universe = np.load(refs/'private_all_source_measured_rank_universe.npz')['native_columns']
    references = [r['reference'] for r in definitions]
    assert len(set(references)) == 50 and np.all(np.diff(universe) > 0)
    signatures = {}
    mapping = []
    for d in definitions:
        columns = np.array(list(d['unique_source_mapping'].values()), int)
        assert len(columns) == len(set(columns)) and len(columns) >= 20
        pos = np.searchsorted(universe, columns)
        assert np.array_equal(universe[pos], columns)
        signatures[d['reference']] = ['F'+str(x) for x in columns]
        mapping.append(pos)
    cap = len(rows)*50*20 + len(rows)*100 + 120*50*4*80 + 2*1024**2
    assert cap < 32*1024**2
    out = root/'GSE176021_disjoint_function_native_RNA_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = {'created_UTC':datetime.now(timezone.utc).isoformat(), 'source_sha256':sha(__file__), 'resource_guard':guard, 'conservative_output_bound_bytes':cap, 'prerequisite_acceptances':required, 'reference_definitions_sha256':sha(refs/'private_all50_full_and_disjoint_reference_gene_definitions.json'), 'scope':'All50 qualified complements,10 original shared pairs/60materials/120pools/11939 same fixed selected cells. No additional cell or pathway selection; no new formal state or gene fit.', 'rank_background':'Every original source-measured authoritative count feature, including zero entries. Unique physical feature required per signature; no missing gene imputation. Official existing pyucell0.8.0 CPU/average ties/skip, maxrank1500 primary and1000/2000 fixed sensitivities. Mean logCP10K uses full original RNA denominators.', 'native_contract':'512-global-row count windows; canonicalize duplicates/zeros/order, compare every selected fullRNA/detected to prior independently accepted measurement. Score in memory only; no new raw-count matrix saved. Full native independent replay and scalar score audit required before scientific review.', 'future_association_contract':'All50 references x2 frozen MCPs x2 predictor lineages x2 response lineages=400 paired-delta slots. Same10 cases/equal-material means. All within/cross lineages explicit. No functional association executed here; full mathematics/input/resource contract required before real tests.', 'limits':'Excluded both-programme/both-view active predictor union. Complements are not full canonical signatures, and independent genes do not imply direct interaction or eliminate shared quality/upstream effects.'}
    (out/'plan_before_complete_functional_RNA_measurement.json').write_text(json.dumps(plan, indent=2)+'\n')
    z = zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'), mode='r')
    counts = z['layers/counts']
    nf = int(counts.attrs['shape'][1])
    nr = int(counts.attrs['shape'][0])
    incidence = csr_matrix((np.concatenate([np.full(len(p),1/len(p)) for p in mapping]), (np.concatenate(mapping), np.concatenate([np.full(len(p),j) for j,p in enumerate(mapping)]))), shape=(len(universe),50))
    uc = np.empty((len(rows),50,3),np.float32)
    ml = np.empty((len(rows),50),np.float64)
    hashes = {name:hashlib.sha256() for name in ['global_row_int64','native_indices_int64','native_data_int64','native_row_endpoints_int64']}
    nnz = 0
    error_log = (out/'private_official_worker_stderr.log').open('xb')
    child = subprocess.Popen([str(runtime),'-u',str(Path(__file__).resolve()),'--worker'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=error_log)
    def send(value):
        pickle.dump(value,child.stdin,protocol=4)
        child.stdin.flush()
    send({'features':['F'+str(x) for x in universe], 'signatures':signatures, 'references':references})
    software = pickle.load(child.stdout)
    try:
        for window in np.unique(rows//512):
            loc = np.flatnonzero(rows//512 == window)
            rr = rows[loc].astype(np.int64)
            start, stop = int(window)*512, min(nr,(int(window)+1)*512)
            ptr = np.asarray(counts['indptr'][start:stop+1],np.int64)
            x = csr_matrix((np.asarray(counts['data'][int(ptr[0]):int(ptr[-1])]), np.asarray(counts['indices'][int(ptr[0]):int(ptr[-1])]), ptr-ptr[0]), shape=(stop-start,nf))[rr-start].copy()
            x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
            assert np.isfinite(x.data).all() and (x.data>0).all() and np.equal(x.data,np.rint(x.data)).all()
            x = x.astype(np.int64)
            rna, detected = np.asarray(x.sum(axis=1)).ravel(), np.diff(x.indptr)
            assert np.array_equal(rna,old['fullRNA'][take[loc]]) and np.array_equal(detected,old['detected'][take[loc]])
            for name, a in [('global_row_int64',rr),('native_indices_int64',x.indices),('native_data_int64',x.data),('native_row_endpoints_int64',x.indptr[1:]+nnz)]:
                hashes[name].update(np.asarray(a,dtype='<i8').tobytes())
            nnz += x.nnz
            y = x[:,universe].tocsr()
            assert np.array_equal(np.asarray(y.sum(axis=1)).ravel(),rna)
            lg = y.astype(float)
            lg.data = np.log1p(lg.data*10000/np.repeat(rna,np.diff(lg.indptr)))
            ml[loc] = (lg@incidence).toarray()
            for first in range(0,len(loc),64):
                last = min(len(loc),first+64)
                send(y[first:last])
                values = pickle.load(child.stdout)
                assert values.shape == (last-first,50,3) and np.isfinite(values).all() and (values >= 0).all() and (values <= 1).all()
                uc[loc[first:last]] = values
            if int(window)%32 == 0:
                (out/'private_progress.json').write_text(json.dumps({'phase':'OFFICIAL_ALL50_NATIVE_FUNCTIONAL_RNA','selected_cells_read':int(loc[-1])+1,'selected_cells':len(rows),'observed_UTC':datetime.now(timezone.utc).isoformat()})+'\n')
        send(None)
        assert child.wait(timeout=60) == 0
    finally:
        if child.poll() is None:
            child.terminate();child.wait(timeout=30)
        error_log.close()
    np.savez_compressed(out/'private_all_selected_cell_functional_scores.npz',global_row=rows,pool_index=pids,references=np.array(references),ucell=uc,mean_logCP10K=ml)
    records = []
    for p in pool.itertuples(index=False):
        ids = np.flatnonzero(pids == p.pool_index)
        assert len(ids) == p.selected_cells
        for j, ref in enumerate(references):
            records.append({'pool_index':p.pool_index,'sampleID':p.sampleID,'view':p.view,'reference':ref,'selected_cells':len(ids),'ucell1500':float(uc[ids,j,0].astype(float).mean()),'ucell1000':float(uc[ids,j,1].astype(float).mean()),'ucell2000':float(uc[ids,j,2].astype(float).mean()),'mean_logCP10K':float(ml[ids,j].mean())})
    pd.DataFrame(records).to_csv(out/'private_all6000_functional_pool_values.tsv',sep='\t',index=False)
    (out/'private_native_canonical_and_software_receipt.json').write_text(json.dumps({'digests':{k:h.hexdigest() for k,h in hashes.items()},'selected_cells':len(rows),'canonical_nonzero_values':int(nnz),'software':software},indent=2)+'\n')
    files = {p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()}
    actual = sum(v['bytes'] for v in files.values())
    assert actual < cap
    receipt = {'status':'COMPLETE_ALL50_NATIVE_FUNCTIONAL_SCORES_FULL_REPLAY_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'files':files,'selected_cells':len(rows),'functional_cell_score_positions':int(uc.size+ml.size),'pool_reference_positions':len(records),'actual_bytes':actual,'no_new_association_or_gene_fit':True}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json'),'bytes':actual}))


if __name__ == '__main__':
    if sys.argv[1:] == ['--worker']:
        worker()
    else:
        main()

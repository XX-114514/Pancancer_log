"""Natural-library reference diagnostic; all RNA/profiles/case predictions RAM only.

Full captured RNA and all original classes remain in the denominator. Unsupported
classes are measured as uncovered RNA, never silently declared absent. This is
an outcome-free measurement diagnostic, not full-atlas cell quantification.
"""
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.linalg import solve_triangular
from scipy.optimize import nnls
from scipy.sparse import csr_matrix
from scipy.stats import rankdata


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def targets(values, basis, marginal):
    names = list(basis); rows = [v for v in values]
    for name, members in marginal.items():
        names.append(name)
        rows.append(values[[basis.index(c) for c in members]].sum(axis=0))
    d = dict(zip(names, rows))
    names.append('paired_R_only_E_minimum')
    rows.append(np.minimum(d['T_anchored_reg_only'], d['T_effector_without_anchored_reg']))
    return names, np.asarray(rows)


def blocks(frame, libraries):
    lookup = {v: i for i, v in enumerate(libraries)}
    parent = list(range(len(libraries)))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    owners = {}
    missing = {'', 'na', 'nan', 'none', 'null', '<na>', 'unknown', 'n/a', 'no_data'}
    for source, library, token in frame[['source', 'library', 'patient_id']].drop_duplicates().itertuples(index=False, name=None):
        token = token.strip()
        if token.casefold() in missing:
            continue
        key = (source, token); i = lookup[library]
        if key in owners:
            parent[find(i)] = find(owners[key])
        else:
            owners[key] = i
    roots = [find(i) for i in range(len(libraries))]
    return pd.Categorical(roots).codes


def solve(reference, y, weights, problem, direct):
    c = problem['A'].T @ (y * weights)
    rhs = solve_triangular(problem['L'], c, lower=True)
    x, _ = nnls(problem['L'].T, rhs, maxiter=1000)
    grad = problem['G'] @ x - c
    scale = max(float(np.max(np.abs(c))), 1e-15)
    kkt = max(float(np.maximum(-grad, 0).max()), float(np.abs(x * grad).max())) / scale
    assert kkt < 1e-7 and (x >= 0).all() and np.isfinite(x).all()
    assert x.sum() > 1e-12
    error = 0.
    if direct:
        other, _ = nnls(problem['A'], y * weights, maxiter=1000)
        error = float(np.max(np.abs(x - other)))
        np.testing.assert_allclose(x, other, atol=1e-6, rtol=1e-6)
    return x, kkt, error


def moments(truth, estimate):
    a, b = rankdata(truth), rankdata(estimate)
    return dict(sum_truth=float(truth.sum()), sum_estimate=float(estimate.sum()),
                sum_truth2=float(truth @ truth), sum_estimate2=float(estimate @ estimate),
                sum_product=float(truth @ estimate), sum_absolute_error=float(np.abs(estimate-truth).sum()),
                sum_squared_error=float(np.sum((estimate-truth)**2)),
                sum_rank_truth=float(a.sum()), sum_rank_estimate=float(b.sum()),
                sum_rank_truth2=float(a @ a), sum_rank_estimate2=float(b @ b), sum_rank_product=float(a @ b))


def fixture(math):
    rng=np.random.default_rng(2026102101)
    reference=rng.uniform(.01,1,(70,19));reference/=reference.sum(axis=0)
    weights=1/np.sqrt(reference.mean(axis=1)+1e-6)
    problem,reason=math.geometry(reference,weights);assert not reason
    truth=rng.dirichlet(np.ones(19))*.8
    estimate,kkt,error=solve(reference,reference@truth,weights,problem,True)
    np.testing.assert_allclose(estimate,truth,atol=1e-7,rtol=1e-7)
    assert abs(estimate.sum()-.8)<1e-7 and abs((estimate/estimate.sum()).sum()-1)<1e-12
    libraries=['S1::L1','S1::L2','S2::L1','S1::L3']
    frame=pd.DataFrame(dict(source=['S1','S1','S2','S1'],library=libraries,patient_id=['001','001','001','NA']))
    b=blocks(frame,libraries)
    assert b[0]==b[1] and b[0]!=b[2] and b[0]!=b[3] and b[2]!=b[3]
    return dict(exact_known_coefficients_recovered=True,unnormalized_mass_distinct_from_normalized_share=True,same_source_claimed_identity_linked=True,cross_source_and_missing_tokens_not_linked=True,full_objective_KKT=kkt,direct_coefficient_error=error)


def main():
    cfg = json.load(sys.stdin); prior = cfg['prior_reference_plan']
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    old_control = Path(os.environ['PANCANCER_PRIOR_REFERENCE_CONTROL'])
    old_output = Path(os.environ['PANCANCER_PRIOR_REFERENCE_OUTPUT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent == Path('/dev/shm') and not out.exists()
    memory = int(next(v.split()[1] for v in Path('/proc/meminfo').read_text().splitlines() if v.startswith('MemAvailable:'))) * 1024
    assert memory > 8 * 1024**3
    assert sha(old_output/'private_terminal.json') == cfg['prior_reference_terminal_SHA256']
    receipt = json.loads((old_control/'private_launch_receipt.json').read_text())
    assert sha(old_control/'run_CM_marker_disjoint_reference_diagnostic_v1.py') == receipt['source_SHA256']
    assert sha(old_control/'qualify_CM_lineage_state_reference_inputs_v1.py') == receipt['helper_SHA256']
    solver = Path(os.environ['PANCANCER_MIXTURE_SOLVER'])
    assert sha(solver) == 'af2705f6e24d0e2e44218cb2fbdb1affe51d4e9cc40744efe0e4d3f6f7bcb5b9'
    spec = importlib.util.spec_from_file_location('accepted_mixture_math', solver)
    math = importlib.util.module_from_spec(spec); spec.loader.exec_module(math)
    sys.path.insert(0, str(old_control))
    from run_CM_marker_disjoint_reference_diagnostic_v1 import select_rows, digest_rows
    from qualify_CM_lineage_state_reference_inputs_v1 import label_T
    atlas = root.parents[1]/'runs/20261013_immunotherapy_atlas_v1'
    assert sha(atlas/'manifests/immune_authority_release_public.json') == prior['release_SHA256']
    out.mkdir(mode=0o700)
    save = lambda name, value: (out/name).write_text(json.dumps(value, indent=2)+'\n')
    save('private_plan_before_native_expression.json', dict(plan=cfg, created_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), memory_available_bytes=memory))
    save('private_numerical_fixture.json',fixture(math))
    z = zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'), mode='r')
    md = pd.DataFrame({k: np.asarray(read_elem(z['obs'][k])).astype(str) for k in ['dataset','library_id','major','patient_id']})
    md['source'] = md.dataset; md['library'] = md.dataset+'::'+md.library_id; md['row'] = np.arange(len(md))
    material_path = atlas/'CM/T_whole_TME_primary_tumor/inputs/sample_metadata_private.tsv'
    assert sha(material_path) == prior['CM_material_SHA256']
    materials = pd.read_csv(material_path, sep='\t', usecols=['sampleID'], keep_default_na=False)
    material = md.library.isin(materials.sampleID).to_numpy(); sources = prior['source_scope']
    assert set(md.loc[material,'source']) == set(sources)
    old = root/'immune_CM_fixed_panel_same_cell_RNA_20261018_v1'
    assert sha(old/'private_fixed34_gene_T_counts_and_depth.npz') == prior['T_cache_SHA256']
    cache = dict(np.load(old/'private_fixed34_gene_T_counts_and_depth.npz', allow_pickle=False))
    tr = cache['global_original_T_row_index']; np.testing.assert_array_equal(tr, np.flatnonzero(md.major.eq('T')))
    fixed = root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    assert sha(fixed) == prior['fixed_panel_contract_SHA256']
    panels = json.loads(fixed.read_text())['panels']
    md['T_capture_state'] = ''
    md.loc[tr,'T_capture_state'] = label_T(cache['raw_gene_counts'], cache['genes'].tolist(), panels, prior['anchor_gene'], prior['effector_genes'], 1)
    md['nonT_major'] = md.major; md['candidate'] = md.major
    md.loc[md.major.eq('T'),'candidate'] = md.loc[md.major.eq('T'),'T_capture_state']
    trows = select_rows(md.loc[material & md.major.eq('T').to_numpy()], 'T_capture_state', prior)
    nrows = select_rows(md.loc[material & md.major.ne('T').to_numpy()], 'nonT_major', prior)
    selected = np.unique(np.r_[trows,nrows]); qualification = cfg['coverage_qualification']
    assert digest_rows(trows) == qualification['T_selected_rows_SHA256']
    assert digest_rows(nrows) == qualification['nonT_selected_rows_SHA256']
    assert digest_rows(selected) == qualification['selected_union_native_rows_SHA256']
    allclasses = cfg['all_original_classes']; core = cfg['supported_reference_classes']
    assert len(allclasses) == 24 and len(core) == 19
    frame = md.loc[material].copy(); native_rows = frame.row.to_numpy(); libraries = sorted(set(frame.library))
    li = pd.Index(libraries).get_indexer(frame.library); lib_source = np.array([sources.index(v.split('::',1)[0]) for v in libraries])
    lib_block = blocks(frame, libraries); nb = int(lib_block.max())+1; nl = len(libraries)
    block_source = np.array([next(iter(set(lib_source[lib_block==b]))) for b in range(nb)])
    assert all(len(set(lib_source[lib_block==b]))==1 for b in range(nb))
    ci = pd.Index(allclasses).get_indexer(frame.candidate); assert (ci >= 0).all()
    ri = pd.Index(core).get_indexer(frame.candidate)
    reference_flag = ri>=0
    nref = np.zeros((nb,len(core)),np.int64); np.add.at(nref,(lib_block[li[reference_flag]],ri[reference_flag]),1)
    selected_frame = frame.loc[reference_flag,['source','library','candidate']].copy()
    selected_frame['block'] = lib_block[li[reference_flag]]
    libref = np.zeros_like(nref)
    for (b,c), group in selected_frame.groupby(['block','candidate']):
        libref[b,core.index(c)] = group.library.nunique()
    for row in qualification['source_class_coverage']:
        g = frame[frame.source.eq(row['source']) & frame.candidate.eq(row['reference_class'])]
        assert len(g)==row['original_cells']
        if row['reference_class'] in core:
            c=core.index(row['reference_class']); bs=block_source==sources.index(row['source'])
            assert nref[bs,c].sum()==row['original_cells'] and libref[bs,c].sum()==row['original_material_units']
    source_order = list(np.asarray(read_elem(z['uns/source_order'])).astype(str))
    measured = np.asarray(read_elem(z['varm/measured_in_source']),bool)
    canonical = np.asarray(read_elem(z['var/canonical_gene_identity'])).astype(str)
    status = np.asarray(read_elem(z['var/feature_identity_status'])).astype(str)
    mp = old/'private_source_gene_measurement_contract.json'; old_manifest=json.loads((old/'private_terminal.json').read_text())
    assert sha(mp)==old_manifest['files'][mp.name]['sha256']; maps=json.loads(mp.read_text())
    excluded={int(v) for m in maps for v in m['physical_columns'] if v>=0}
    good=measured[:,[source_order.index(s) for s in sources]].all(axis=1)
    good &= ~np.char.startswith(canonical,'AMBIGUOUS::') & (status!='ambiguous_source_feature') & ~pd.Index(canonical).duplicated(keep=False)
    good[list(excluded)]=False; common=np.flatnonzero(good)
    previous=json.loads((old_output/'private_plan_before_RNA_and_model_fitting.json').read_text())
    assert len(common)==25188 and digest_rows(common)==previous['common_physical_gene_indices_SHA256']
    save('private_metadata_and_selection_contract.json',dict(material_cells=len(frame),material_libraries=nl,conservative_claimed_identity_components=nb,all_original_classes=len(allclasses),supported_reference_classes=len(core),selection_SHA256=digest_rows(selected),common_axis_SHA256=digest_rows(common),individual_keys_exported=False))
    bulk=np.zeros((nl,len(common))); native_total=np.zeros(nl); truth_mass=np.zeros((nl,len(allclasses)))
    ref=np.zeros((nb*len(core),len(common))); ref_mass=np.zeros(nb*len(core)); ref_num=np.zeros(nb*len(core),np.int64)
    bulk_gene_proof=np.zeros(len(common)); ref_gene_proof=np.zeros(len(common)); cells_read=0
    counts=z['layers/counts']; nc,nf=map(int,counts.attrs['shape'])
    for start in range(0,nc,cfg['chunk_rows']):
        stop=min(nc,start+cfg['chunk_rows']); a=int(np.searchsorted(native_rows,start)); b=int(np.searchsorted(native_rows,stop))
        if a==b:continue
        ptr=np.asarray(counts['indptr'][start:stop+1],np.int64); first,last=int(ptr[0]),int(ptr[-1])
        x=csr_matrix((np.asarray(counts['data'][first:last]),np.asarray(counts['indices'][first:last]),ptr-ptr[0]),shape=(stop-start,nf))[native_rows[a:b]-start]
        x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
        assert np.isfinite(x.data).all() and (x.data>0).all() and np.array_equal(x.data,np.floor(x.data))
        src=frame.source.to_numpy()[a:b]
        for s in set(src):assert measured[x[src==s].indices,source_order.index(s)].all()
        total=np.asarray(x.sum(axis=1)).ravel().astype(np.int64); tt=frame.major.to_numpy()[a:b]=='T'
        if tt.any():
            pos=np.searchsorted(tr,native_rows[a:b][tt]);np.testing.assert_array_equal(tr[pos],native_rows[a:b][tt])
            np.testing.assert_array_equal(total[tt],cache['full_native_total_counts'][pos]);np.testing.assert_array_equal(np.diff(x.indptr)[tt],cache['full_native_detected_genes'][pos])
            for s in set(src[tt]):
                take=tt&(src==s);cols=next(m for m in maps if m['dataset']==s)['physical_columns']
                np.testing.assert_array_equal(x[take][:,cols].toarray(),cache['raw_gene_counts'][np.searchsorted(tr,native_rows[a:b][take])])
        sub=x[:,common]; local=li[a:b]
        local_libraries,inverse=np.unique(local,return_inverse=True)
        incidence=csr_matrix((np.ones(b-a),(inverse,np.arange(b-a))),shape=(len(local_libraries),b-a))
        bulk[local_libraries]+=np.asarray((incidence@sub).toarray()); native_total+=np.bincount(local,weights=total,minlength=nl)
        truth_mass+=np.bincount(local*len(allclasses)+ci[a:b],weights=total,minlength=nl*len(allclasses)).reshape(nl,len(allclasses))
        bulk_gene_proof+=np.asarray(sub.sum(axis=0)).ravel()
        use=reference_flag[a:b]
        if use.any():
            key=lib_block[local[use]]*len(core)+ri[a:b][use]
            local_keys,inverse=np.unique(key,return_inverse=True)
            inc=csr_matrix((np.ones(int(use.sum())),(inverse,np.arange(int(use.sum())))),shape=(len(local_keys),int(use.sum())))
            ref[local_keys]+=np.asarray((inc@sub[use]).toarray());ref_mass+=np.bincount(key,weights=total[use],minlength=len(ref));ref_num+=np.bincount(key,minlength=len(ref))
            ref_gene_proof+=np.asarray(sub[use].sum(axis=0)).ravel()
        cells_read+=b-a;save('private_progress.json',dict(phase='READING_ORIGINAL_LIBRARY_RNA',cells_read=cells_read,expected_cells=len(frame)))
    assert cells_read==len(frame) and (native_total>0).all()
    np.testing.assert_array_equal(bulk.sum(axis=0),bulk_gene_proof);np.testing.assert_array_equal(ref.sum(axis=0),ref_gene_proof)
    np.testing.assert_array_equal(truth_mass.sum(axis=1),native_total);np.testing.assert_array_equal(ref_num.reshape(nb,len(core)),nref)
    ref=ref.reshape(nb,len(core),len(common));ref_mass=ref_mass.reshape(nb,len(core))
    bulk/=native_total[:,None];truth=truth_mass/native_total[:,None]; core_truth=truth[:,[allclasses.index(c) for c in core]]
    coarse=prior['T_classes']+['NK','Other_supported_reference_pool','Boundary_supported_reference_pool']
    mapping=np.zeros((len(core),len(coarse)))
    for i,c in enumerate(core):mapping[i,coarse.index(c if c in prior['T_classes']+['NK'] else 'Boundary_supported_reference_pool' if c in prior['boundary_majors'] else 'Other_supported_reference_pool')]=1
    source_raw=[ref[block_source==s].sum(axis=0) for s in range(len(sources))]
    source_mass=[ref_mass[block_source==s].sum(axis=0) for s in range(len(sources))]
    runs=defaultdict(list); failures=Counter(); numeric=[]; checked_direct=set()
    for si,source in enumerate(sources):
        lib_indices=np.flatnonzero(lib_source==si)
        for domain in cfg['source_domains']:
            groups=[(-1,lib_indices)] if domain=='heldout_source' else [(int(b),lib_indices[lib_block[lib_indices]==b]) for b in sorted(set(lib_block[lib_indices]))]
            for block,queries in groups:
                if domain=='heldout_source':
                    train=[s for s in range(len(sources)) if s!=si]
                    rf=np.mean([source_raw[s]/source_mass[s][:,None] for s in train],axis=0).T
                    rc=np.mean([(mapping.T@source_raw[s])/(source_mass[s]@mapping)[:,None] for s in train],axis=0).T
                else:
                    count=nref[block_source==si].sum(axis=0)-nref[block]; libs=libref[block_source==si].sum(axis=0)-libref[block]
                    if np.any(count<20) or np.any(libs<2):
                        for architecture in cfg['reference_architectures']:
                            for scaling in cfg['scalings']:failures[(source,domain,architecture,scaling,'REQUIRED_TRAINING_CLASS_SUPPORT')]+=len(queries)
                        continue
                    raw=source_raw[si]-ref[block];mass=source_mass[si]-ref_mass[block]
                    rf=(raw/mass[:,None]).T;rc=((mapping.T@raw)/(mass@mapping)[:,None]).T
                for scaling in cfg['scalings']:
                    weights=1/np.sqrt(rf.mean(axis=1)+1e-6) if scaling=='train_mean_weighted' else np.ones(len(common))
                    for architecture,reference,basis,true_base in [('separate_original_major',rf,core,core_truth),('merge_supported_nonT_background',rc,coarse,core_truth@mapping)]:
                        problem,reason=math.geometry(reference,weights)
                        if problem is None:
                            failures[(source,domain,architecture,scaling,reason)]+=len(queries);continue
                        for qi,lib in enumerate(queries):
                            proof_key=(source,domain,architecture,scaling)
                            direct=proof_key not in checked_direct
                            x,kkt,error=solve(reference,bulk[lib],weights,problem,direct)
                            if direct:checked_direct.add(proof_key)
                            numeric.append(dict(kkt=kkt,direct=direct,direct_error=error))
                            names,t=targets(true_base[lib,:,None],basis,cfg['marginal_targets'])
                            for variant,estimate in [('raw_nonnegative_coefficients',x),('normalized_supported_reference_sum',x/x.sum())]:
                                nn,pred=targets(estimate[:,None],basis,cfg['marginal_targets']);assert nn==names
                                runs[(source,domain,architecture,scaling,variant)].append((int(lib),int(lib_block[lib]),names,t[:,0],pred[:,0],float(x.sum())))
                save('private_progress.json',dict(phase='FITTING_FIXED_REFERENCE_CONTEXTS',source=source,domain=domain,completed_model_fits=len(numeric)))
    metrics=[];contexts=[]
    for key,rows in sorted(runs.items()):
        source,domain,architecture,scaling,variant=key;names=rows[0][2]
        actual=np.array([r[3] for r in rows]);predicted=np.array([r[4] for r in rows])
        context=dict(source=source,domain=domain,architecture=architecture,scaling=scaling,coefficient_view=variant)
        contexts.append(dict(context,evaluated_libraries=len(rows),conservative_blocks=len(set(r[1] for r in rows)),expected_libraries=int(np.sum(lib_source==sources.index(source))),coefficient_sum_min=min(r[5] for r in rows),coefficient_sum_max=max(r[5] for r in rows)))
        for j,name in enumerate(names):
            s=math.summarize(actual[:,j],predicted[:,j]);s['libraries']=s.pop('mixtures')
            metrics.append(dict(context,target=name,**s,**moments(actual[:,j],predicted[:,j])))
    existing=set(runs)
    for source in sources:
        for domain in cfg['source_domains']:
            for architecture,basis in [('separate_original_major',core),('merge_supported_nonT_background',coarse)]:
                for scaling in cfg['scalings']:
                    for variant in cfg['coefficient_views']:
                        key=(source,domain,architecture,scaling,variant)
                        if key in existing:continue
                        context=dict(source=source,domain=domain,architecture=architecture,scaling=scaling,coefficient_view=variant)
                        contexts.append(dict(context,evaluated_libraries=0,conservative_blocks=0,expected_libraries=int(np.sum(lib_source==sources.index(source))),coefficient_sum_min=None,coefficient_sum_max=None,status='UNESTIMABLE_ALL_LIBRARIES'))
                        names,_=targets(np.zeros((len(basis),1)),basis,cfg['marginal_targets'])
                        for name in names:metrics.append(dict(context,target=name,libraries=0,status='UNESTIMABLE_ALL_LIBRARIES',MAE=None,RMSE=None,bias=None,Pearson=None,Spearman=None))
    assert len(contexts)==64 and len(metrics)==1376
    paired=[]
    for source in sources:
        for domain in cfg['source_domains']:
            for scaling in cfg['scalings']:
                for variant in cfg['coefficient_views']:
                    a=runs.get((source,domain,'separate_original_major',scaling,variant),[]);b=runs.get((source,domain,'merge_supported_nonT_background',scaling,variant),[])
                    aa={r[0]:r for r in a};bb={r[0]:r for r in b};shared=sorted(set(aa)&set(bb))
                    if not shared:continue
                    for target in sorted(set(a[0][2])&set(b[0][2])):
                        av=[];bv=[]
                        for lib in shared:
                            x,y=aa[lib],bb[lib];i=x[2].index(target);j=y[2].index(target);assert abs(x[3][i]-y[3][j])<1e-12
                            av.append(abs(x[4][i]-x[3][i]));bv.append(abs(y[4][j]-y[3][j]))
                        paired.append(dict(source=source,domain=domain,scaling=scaling,coefficient_view=variant,target=target,matched_libraries=len(shared),fine_MAE=float(np.mean(av)),merged_MAE=float(np.mean(bv)),fine_minus_merged_MAE=float(np.mean(np.asarray(av)-bv)),inferential_P_or_FDR=False))
    coverage=[]
    missing=[allclasses.index(c) for c in cfg['unsupported_reference_classes']]
    for si,source in enumerate(sources):
        take=lib_source==si;u=truth[take][:,missing].sum(axis=1)
        coverage.append(dict(source=source,material_cells=int(frame.source.eq(source).sum()),material_libraries=int(take.sum()),conservative_blocks=int(len(set(lib_block[take]))),unsupported_RNA_fraction_mean=float(u.mean()),unsupported_RNA_fraction_max=float(u.max()),common_gene_RNA_fraction_mean=float(bulk[take].sum(axis=1).mean()),all24_truth_mass_conserved=True,patient_independence_not_newly_verified=True))
    save('private_all_context_coverage.json',contexts);save('private_all_target_recovery_metrics.json',metrics);save('private_paired_background_comparisons.json',paired);save('private_all_source_uncovered_RNA.json',coverage)
    save('private_all_context_failures.json',[dict(source=k[0],domain=k[1],architecture=k[2],scaling=k[3],reason=k[4],libraries=v) for k,v in sorted(failures.items())])
    save('private_full_numeric_proofs.json',dict(model_fits=len(numeric),KKT_max=max(r['kkt'] for r in numeric),direct_checks=sum(r['direct'] for r in numeric),direct_max_coefficient_error=max(r['direct_error'] for r in numeric),all_material_T_cache_counts_depth_detected_exact=True,full_gene_aggregation_sum_proofs=True,all24_native_RNA_mass_conserved=True,all_original_selection_hashes_exact=True))
    files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()}
    assert sum(v['bytes'] for v in files.values())<cfg['output_cap_bytes']
    terminal=dict(status='COMPLETE_NATURAL_LIBRARY_REFERENCE_DIAGNOSTIC_SOURCE_REVIEW_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),all_raw_RNA_profiles_predictions_RAM_only=True,no_clinical_tests_or_formal_identity_adoption=True,files=files)
    save('private_terminal.json',terminal)
    print(json.dumps(dict(status=terminal['status'],terminal_SHA256=sha(out/'private_terminal.json'))))


if __name__=='__main__':
    main()

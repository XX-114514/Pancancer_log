"""Private I/O and outcome-free scope for complete native MFA interpretation.

No learned feature names, patient rows or weights may be printed or published.
"""
import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

VIEWS = ['author_CD8', 'author_Treg', 'author_Macrophage', 'author_B_Plasma']
PANELS = {
    'SEX_LINKED_RNA': 'XIST TSIX KDM5D DDX3Y EIF1AY UTY ZFY USP9Y TMSB4Y RPS4Y1'.split(),
    'IEG_STRESS_RNA': 'FOS JUN JUNB JUND ATF3 EGR1 HSPA1A HSPA1B HSP90AA1 DNAJB1 DDIT3'.split(),
    'T_IDENTITY': 'CD3D CD3E TRAC CD247 LCK CD2 IL7R'.split(),
    'CYTOTOXIC_TNK': 'NKG7 GNLY PRF1 GZMB GZMH CCL5 CST7'.split(),
    'TREG_IDENTITY': 'FOXP3 IL2RA CTLA4 IKZF2 TIGIT TNFRSF18'.split(),
    'MYELOID_IDENTITY': 'LST1 TYROBP FCER1G CTSS AIF1 CSF1R LYZ'.split(),
    'B_IDENTITY': 'CD79A CD79B MS4A1 CD37 CD74 CD19'.split(),
    'PLASMA_IDENTITY': 'MZB1 JCHAIN SDC1 DERL3 XBP1 IGKC IGHG1 IGHA1'.split(),
}
BASIS = {
    'NSCLC_native_MFA_inputs_20261016_v1/private_terminal.json': '4ee0c9fc89eabe08c41545d36cd9ef4191dbe05294efb8038fa5031ef88ecd1c',
    'NSCLC_native_MFA_inputs_20261016_v1/private_independent_MFA_input_audit.json': '3b6b8086841151de638b8bb6f7f1373fc856debeae2742c392562c1c00607b8d',
    'NSCLC_native_MFA_discovery_20261016_v3/private_independent_MFA_discovery_audit.json': '3ebc1087f7d3030c2d6b36caa96f2869a4bcdfb44e5deebbe03651276cb9ad3f',
    'NSCLC_native_MFA_holdout_20261016_v2/private_independent_MFA_holdout_audit.json': 'fcc4a9144139b2ac247b93443dd557fab9381d5f75733b1a7340aba9ed72a081',
    'NSCLC_native_MFA_holdout_complete_review_20261016_v1/PRIVATE_NEXT_ALL_FACTOR_FUNCTION_AND_MIXTURE_SCOPE.json': '720e394d8a272dd998b88e734b37aeafc929498aab93880f6c968fd7d470bf2c',
    'NSCLC_full_native_gene_units_20261016_v1/private_terminal.json': '8e8d81678227a506842315b1b2a5e28b604c0f350e475c62d589ea7e0d9719bb',
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(16777216), b''):
            h.update(block)
    return h.hexdigest()


def quota_gate(root, cap=256 * 1024**2):
    result = subprocess.run(['quota', '-g', '-s', 'USER002'], capture_output=True, text=True, timeout=20)
    assert result.returncode == 0
    m = re.search(r'^\s*([0-9.]+)([KMGT])\s+\S+\s+([0-9.]+)([KMGT])', result.stdout, re.M)
    assert m, 'Cannot resolve group quota'
    scale = {'K': 1024, 'M': 1024**2, 'G': 1024**3, 'T': 1024**4}
    remaining = float(m[3]) * scale[m[4]] - float(m[1]) * scale[m[2]]
    free = shutil.disk_usage(root).free
    assert min(remaining, free) > 500_000_000_000 + 1024**3 + cap, '500GB reserve plus rounding buffer and output cap'
    return dict(rounded_group_remaining_bytes=remaining, filesystem_free_bytes=free,
                reserve_bytes=500_000_000_000, rounding_buffer_bytes=1024**3, output_cap_bytes=cap)


def load():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    for name, digest in BASIS.items():
        assert sha(root / name) == digest, 'Frozen upstream evidence changed'
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    terminal = json.loads((native / 'private_terminal.json').read_text())
    for name in ['private_all_native_units_and_quality.tsv', 'private_native_gene_axis.tsv', 'mean_logCP10K.npy']:
        assert sha(native / name) == terminal['outputs'][name]['sha256']
    units = pd.read_csv(native / 'private_all_native_units_and_quality.tsv', sep='\t', keep_default_na=False)
    genes = pd.read_csv(native / 'private_native_gene_axis.tsv', sep='\t', keep_default_na=False).gene.to_numpy(str)
    assert len(set(genes)) == len(genes) == 31831
    inp = root / 'NSCLC_native_MFA_inputs_20261016_v1'
    cases = pd.read_csv(inp / 'private_case_axis.tsv', sep='\t', keep_default_na=False)
    prepared = []
    for view in VIEWS:
        with np.load(inp / ('private_' + view + '_prepared.npz'), allow_pickle=False) as z:
            prepared.append({k: z[k] for k in z.files})
    model_path = root / 'NSCLC_native_MFA_discovery_20261016_v3/seed_20261017/private_frozen_parameters_and_own_view_scores.npz'
    with np.load(model_path, allow_pickle=False) as z:
        model = {k: z[k] for k in z.files}
    assert np.array_equal(model['cases'].astype(str), cases.sampleID.to_numpy(str))
    assert model['views'].astype(str).tolist() == VIEWS
    train = model['train'].astype(bool)
    assert train.sum() == 145 and (~train).sum() == 63
    assert np.array_equal(train, cases.partition.eq('discovery').to_numpy())
    for j, p in enumerate(prepared):
        assert np.array_equal(p['cases'].astype(str), cases.sampleID.to_numpy(str))
        assert np.array_equal(p['train'].astype(bool), train)
        assert np.array_equal(genes[p['selected_native_columns']], p['selected_genes'].astype(str))
        assert p['design'].shape[1] == p['quality'].shape[1] + 1
        assert model['weights_' + str(j)].shape == (2000, 10)
    union = np.unique(np.concatenate([p['selected_native_columns'] for p in prepared]))
    remaining = np.setdiff1d(np.arange(len(genes)), union)
    refroot = Path(__file__).parent / 'hallmark50_public_reference_v3'
    provenance = json.loads((refroot / 'reference_provenance.json').read_text())
    refs = {}
    for record in provenance['records']:
        path = refroot / (record['name'] + '.json')
        assert sha(path) == record['sha256']
        values = json.loads(path.read_text())[record['name']]['geneSymbols']
        # Public reference symbols are normalized to the native uppercase axis.
        refs[record['name']] = set(str(g).strip().upper() for g in values)
    assert len(refs) == 50
    plan = json.loads((root / 'NSCLC_native_multicell_discovery_preflight_20261016_v1/plan_before_any_native_discovery_effect.json').read_text())
    labels = {v: plan['view_author_labels'][v] for v in VIEWS}
    meta = root / 'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz'
    assert sha(meta) == '24155466454df53d01eb6cbbfcd56abeb45adb68f74e721af22a9dce41040b2f'
    tiny = []
    for chunk in pd.read_csv(meta, usecols=['sampleID', 'gender', 'age'], chunksize=200000, keep_default_na=False):
        tiny.append(chunk.drop_duplicates())
    tiny = pd.concat(tiny).drop_duplicates()
    assert tiny.groupby('sampleID').size().eq(1).all()
    demographic = cases[['sampleID']].merge(tiny, on='sampleID', how='left', validate='one_to_one')
    assert demographic.gender.notna().all()
    demographic['female'] = demographic.gender.map({'F': 1., 'M': 0.})
    demographic['age_numeric'] = pd.to_numeric(demographic.age, errors='coerce')
    demographic.loc[~demographic.age_numeric.between(1, 130), 'age_numeric'] = np.nan
    return dict(root=root, native=native, units=units, genes=genes, cases=cases, prepared=prepared,
                model=model, train=train, union=union, remaining=remaining, refs=refs,
                labels=labels, demographic=demographic, refroot=refroot, model_path=model_path)


def pool_native(data):
    """Pool all four views in gene blocks, with true cell-count denominators."""
    cases = data['cases'].sampleID.to_numpy(str)
    ci = {s: i for i, s in enumerate(cases)}
    u = data['units']
    row, col, count = [], [], []
    subtype_counts = []
    for vi, view in enumerate(VIEWS):
        labels = data['labels'][view]
        table = u[u.sub_cell_type.isin(labels)].pivot(index='sampleID', columns='sub_cell_type', values='n_cells')
        table = table.reindex(index=cases, columns=labels).fillna(0).astype(int)
        assert (table.sum(axis=1) >= 20).all()
        subtype_counts.append(table.to_numpy())
        for item in u[u.sampleID.isin(ci) & u.sub_cell_type.isin(labels)].itertuples():
            row.append(vi * len(cases) + ci[item.sampleID]); col.append(item.unit_index); count.append(item.n_cells)
    total = np.bincount(row, weights=count, minlength=4 * len(cases))
    weights = np.asarray(count, float) / total[np.asarray(row)]
    mixing = csr_matrix((weights, (row, col)), shape=(4 * len(cases), len(u)))
    source = np.load(data['native'] / 'mean_logCP10K.npy', mmap_mode='r')
    pooled = np.empty((4 * len(cases), len(data['genes'])), float)
    for start in range(0, source.shape[1], 512):
        pooled[:, start:start+512] = mixing @ source[:, start:start+512]
    assert np.isfinite(pooled).all() and (pooled >= 0).all()
    return pooled.reshape(4, len(cases), -1), subtype_counts

# Ingest the 13 accepted supplementary sources (inventories/source_search/candidates_detailed_check_v1.tsv,
# verdict=pass) from their downloaded GEO files into one raw-count h5ad per source.
# Only tumour samples are kept (adjacent normal, normal lymph node and blood samples are dropped, as
# in the existing analysis set); sample-level patient / sampling site / treatment come from GEO
# characteristics and the papers (rules below, one block per source). Genes are keyed by symbol and
# duplicate symbols are summed. QC is a separate step (qc_supplement_v1.py).
# Usage: ingest_supplement_v1.py RAW_DIR OUT_DIR [SOURCE ...]
import sys, os, re, io, tarfile, gzip, json
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.io, anndata as ad, h5py

raw, out = sys.argv[1:3]
only = set(sys.argv[3:])
os.makedirs(out, exist_ok=True)


def collapse(X, genes):
    """sum columns with the same symbol; X cells x genes"""
    genes = np.asarray(genes).astype(str)
    u, inv = np.unique(genes, return_inverse=True)
    if len(u) == len(genes):
        order = np.argsort(genes)
        return X.tocsc()[:, order].tocsr(), genes[order]
    agg = sp.csr_matrix((np.ones(len(genes)), (np.arange(len(genes)), inv)), shape=(len(genes), len(u)))
    return (X @ agg).tocsr(), u


def as_counts(X):
    X = sp.csr_matrix(X)
    if X.nnz and not np.allclose(X.data[:200000], np.round(X.data[:200000])):
        raise ValueError('non-integer counts')
    X.data = np.rint(X.data)
    return X.astype(np.int32)


def tar_members(path):
    t = tarfile.open(path)
    return t, {m.name: m for m in t.getmembers()}


def read_mm(fobj):
    """MatrixMarket coordinate reader tolerant of blank '%' comment lines (Singleron); returns cells x genes"""
    head = None
    while True:
        line = fobj.readline().decode().strip()
        if line.startswith('%') or not line:
            continue
        head = [int(x) for x in line.split()]; break
    nr, nc, nnz = head
    t = pd.read_csv(fobj, sep=r'\s+', header=None, dtype={0: np.int32, 1: np.int32, 2: np.float64}, engine='c')
    assert len(t) == nnz, (len(t), nnz)
    return sp.csr_matrix((t[2].values, (t[1].values - 1, t[0].values - 1)), shape=(nc, nr))


def floor_call(X, bc):
    """cell calling for unfiltered matrices: keep barcodes passing the QC floor (>=200 genes, >=500 UMI)"""
    X = sp.csr_matrix(X)
    ok = (np.diff(X.indptr) >= 200) & (np.asarray(X.sum(1)).ravel() >= 500)
    return X[ok], bc[ok]


def read_mtx_triplet(open_fn, prefix):
    # file objects from tar members / open() are not decompressed automatically
    gz = lambda name: gzip.open(open_fn(name))
    mtx = read_mm(gz(prefix + 'matrix.mtx.gz'))
    bc = pd.read_csv(gz(prefix + 'barcodes.tsv.gz'), header=None, sep='\t')[0].astype(str).values
    ft = None
    for fname in ('features.tsv.gz', 'genes.tsv.gz'):
        try:
            ft = pd.read_csv(gz(prefix + fname), header=None, sep='\t'); break
        except (KeyError, FileNotFoundError):
            continue
    if ft.shape[1] >= 3:
        ft = ft[ft[2].astype(str).str.contains('Gene Expression')] if (ft[2] != 'Gene Expression').any() else ft
        mtx = mtx[:, ft.index.values]
    sym = ft[1] if ft.shape[1] >= 2 else ft[0]
    return mtx, bc, sym.astype(str).values


def read_10x_h5(fobj):
    with h5py.File(fobj, 'r') as h:
        g = h['matrix']
        X = sp.csc_matrix((g['data'][:], g['indices'][:], g['indptr'][:]), shape=tuple(g['shape'][:])).T.tocsr()
        bc = g['barcodes'][:].astype(str)
        f = g['features']
        keep = f['feature_type'][:].astype(str) == 'Gene Expression'
        return X[:, keep], bc, f['name'][:].astype(str)[keep]


def dense_to_sparse(path_or_buf, sep='\t', index_col=0, chunksize=1000):
    rows, genes = [], []
    cells = None
    for ch in pd.read_csv(path_or_buf, sep=sep, index_col=index_col, chunksize=chunksize):
        cells = ch.columns.astype(str).values
        genes.extend(ch.index.astype(str))
        rows.append(sp.csr_matrix(ch.values))
    M = sp.vstack(rows).T.tocsr()
    return M, cells, np.asarray(genes)


def build(src, parts, cancer_main, extra_uns):
    """parts: list of (X, barcodes, genes, obs-dict)"""
    allg = sorted(set().union(*[set(p[2]) for p in parts]))
    gi = {g: i for i, g in enumerate(allg)}
    mats, obs = [], []
    for X, bc, genes, meta in parts:
        X, genes = collapse(as_counts(X), genes)
        cols = np.array([gi[g] for g in genes])
        X = sp.csr_matrix((X.data, cols[X.indices], X.indptr), shape=(X.shape[0], len(allg)))
        mats.append(X)
        o = pd.DataFrame({'origin_barcode': bc})
        for k, v in meta.items():
            o[k] = v
        obs.append(o)
    X = sp.vstack(mats).tocsr()
    o = pd.concat(obs, ignore_index=True)
    o['source'] = src
    o['cell_id'] = src + '|' + o['sample_id'].astype(str) + '|' + o['origin_barcode'].astype(str)
    if o.cell_id.duplicated().any():
        raise ValueError(f'{src}: duplicated cell ids')
    o.index = o.cell_id.values
    if 'analysis_main' not in o:
        o['analysis_main'] = cancer_main
    a = ad.AnnData(X=X, obs=o, var=pd.DataFrame(index=allg))
    a.uns['ingest'] = json.dumps(extra_uns)
    a.write_h5ad(os.path.join(out, f'{src}.h5ad'), compression='gzip')
    s = o.groupby('sample_id').size()
    print(json.dumps({'source': src, 'cells': int(a.n_obs), 'genes': int(a.n_vars), 'samples': int(len(s)),
                      'patients': int(o.patient_id.nunique())}), flush=True)


def common(**kw):
    base = {'treatment_status': 'unknown', 'treatment_detail': '', 'fresh_frozen': 'unknown',
            'capture_scope': 'whole_tissue', 'suspension': 'cell', 'cell_calling': 'author_filtered_matrix'}
    base.update(kw)
    return base


R = lambda s, *p: os.path.join(raw, s, *p)
run = lambda s: not only or s in only

# ---- Papillary thyroid carcinoma
if run('GSE163203'):  # Microwell-seq dge per sample; PTC samples only (Adj_ dropped)
    t, mem = tar_members(R('GSE163203', 'GSE163203_RAW.tar'))
    parts = []
    for n in sorted(mem):
        if 'Adj' in n or not n.endswith('_dge.txt.gz'):
            continue
        sample = re.sub(r'^GSM\d+_|_dge\.txt\.gz$', '', n)
        pid = re.search(r'_(\d+)(_\d+)?$', sample).group(1)
        X, cells, genes = dense_to_sparse(gzip.open(t.extractfile(mem[n])), sep=r'\s+')
        parts.append((X, cells, genes, common(sample_id=sample, gsm=n.split('_')[0], patient_id=f'GSE163203:{pid}',
                     sampling_site='primary_organ_site', material='PTC primary tumour', fresh_frozen='fresh',
                     technology='Microwell-seq', hashimoto='withHT' in sample)))
    build('GSE163203', parts, 'Papillary thyroid carcinoma',
          {'evidence': 'GEO GSM titles (PTC with/without Hashimoto; adjacent samples dropped); PMID 34805166 fresh tissue'})

if run('GSE191288'):  # 10x h5 per lobe; non-tumour sample of patient 4 dropped
    t, mem = tar_members(R('GSE191288', 'GSE191288_RAW.tar'))
    parts = []
    for n in sorted(mem):
        m = re.match(r'(GSM\d+)_T(\d)([LR])\.h5$', n)
        if not m:
            continue
        X, bc, genes = read_10x_h5(io.BytesIO(t.extractfile(mem[n]).read()))
        parts.append((X, bc, genes, common(sample_id=f'T{m.group(2)}{m.group(3)}', gsm=m.group(1), patient_id=f'GSE191288:P{m.group(2)}',
                     sampling_site='primary_organ_site', material=f'PTC {"left" if m.group(3) == "L" else "right"} lobe',
                     fresh_frozen='fresh', technology='10x')))
    build('GSE191288', parts, 'Papillary thyroid carcinoma', {'evidence': 'GEO titles (left/right lobe per patient); PMID 35515000 freshly resected'})

# ---- Osteosarcoma
if run('GSE162454'):
    t, mem = tar_members(R('GSE162454', 'GSE162454_RAW.tar'))
    parts = []
    for pre in sorted({re.sub(r'(barcodes|features|matrix)\.(tsv|mtx)\.gz$', '', n) for n in mem}):
        X, bc, genes = read_mtx_triplet(lambda f: t.extractfile(mem[f]), pre)
        sample = re.sub(r'^GSM\d+_|_$', '', pre)
        parts.append((X, bc, genes, common(sample_id=sample, gsm=pre.split('_')[0], patient_id=f'GSE162454:{sample}',
                     sampling_site='unspecified', material='osteosarcoma tumour tissue', treatment_status='untreated_or_pretreatment',
                     treatment_detail='treatment-naive (PMID 34367994 title)', fresh_frozen='fresh', technology='10x')))
    build('GSE162454', parts, 'Osteosarcoma', {'evidence': 'PMID 34367994: treatment-naive OS, fresh specimens'})

if run('GSE299023'):
    t, mem = tar_members(R('GSE299023', 'GSE299023_RAW.tar'))
    inst = {'GSM9030790': 'CCMC', 'GSM9030791': 'CCMC', 'GSM9030792': 'CCMC', 'GSM9030793': 'CHCO', 'GSM9030794': 'CHCO', 'GSM9030795': 'CHCO'}
    parts = []
    for pre in sorted({re.sub(r'(barcodes|features|matrix)\.(tsv|mtx)\.gz$', '', n) for n in mem}):
        X, bc, genes = read_mtx_triplet(lambda f: t.extractfile(mem[f]), pre)
        gsm = pre.split('_')[0]; sample = re.sub(r'^GSM\d+_|_$', '', pre)
        parts.append((X, bc, genes, common(sample_id=sample, gsm=gsm, patient_id=f'GSE299023:{sample}', sampling_site='primary_organ_site',
                     material='osteosarcoma primary tumour biopsy', treatment_status='untreated_or_pretreatment',
                     treatment_detail='treatment-naive biopsy (PMID 40647416)',
                     fresh_frozen='viably_cryopreserved' if inst[gsm] == 'CHCO' else 'unknown', technology='10x',
                     capture_scope='whole_tissue_viable_cell_FACS', institution=inst[gsm])))
    build('GSE299023', parts, 'Osteosarcoma', {'evidence': 'PMID 40647416: untreated primary tumour biopsies; CHCO viably cryopreserved; viable-cell FACS'})

# ---- PDAC
if run('GSE205013'):
    gsm = pd.read_json(R('..', '..', '..', 'analysis', 'source_search_v1', 'detailed_v1', 'GSE205013.json'), typ='series')['samples']
    info = {s['gsm']: s for s in gsm}
    t, mem = tar_members(R('GSE205013', 'GSE205013_RAW.tar'))
    parts = []
    for pre in sorted({re.sub(r'(barcodes|features|matrix)\.(tsv|mtx)\.gz$', '', n) for n in mem}):
        X, bc, genes = read_mtx_triplet(lambda f: t.extractfile(mem[f]), pre)
        g = pre.split('_')[0]; s = info[g]; ch = s['chars']
        met = 'liver met' in ch.get('tissue', '')
        parts.append((X, bc, genes, common(sample_id=s['title'], gsm=g, patient_id=f'GSE205013:{s["title"]}',
                     sampling_site='metastatic_site' if met else 'primary_organ_site', material=ch.get('tissue', ''),
                     treatment_status='on_or_after_treatment' if ch.get('treatment') == 'Treated' else 'untreated_or_pretreatment',
                     treatment_detail=f'{ch.get("treatment")}; {ch.get("procedure")}; {ch.get("stage")}', fresh_frozen='fresh', technology='10x')))
    build('GSE205013', parts, 'Pancreatic ductal adenocarcinoma', {'evidence': 'GEO characteristics tissue/treatment/procedure; PMID 36781852 fresh'})

if run('GSE242230'):
    ann = pd.read_csv(R('GSE242230', 'GSE242230_annotations.txt.gz'), sep='\t')
    t, mem = tar_members(R('GSE242230', 'GSE242230_RAW.tar'))
    parts = []
    for pre in sorted({re.sub(r'(barcodes|features|matrix)\.(tsv|mtx)\.gz$', '', n) for n in mem}):
        X, bc, genes = read_mtx_triplet(lambda f: t.extractfile(mem[f]), pre)
        sample = re.sub(r'^GSM\d+_|_$', '', pre)
        a = ann[ann.sample_id == sample].set_index('cell_id')
        key = pd.Index([f'{sample}_{b}' for b in bc])
        # matrices are unfiltered: keep only cells listed in the author annotation file
        inann = key.isin(a.index)
        X, bc, key = X[inann], bc[inann], key[inann]
        parts.append((X, bc, genes, common(sample_id=sample, gsm=pre.split('_')[0], patient_id=f'GSE242230:{sample}', cell_calling='author_annotation_list',
                     sampling_site='primary_organ_site', material='EUS-FNB primary PDAC', treatment_status='untreated_or_pretreatment',
                     treatment_detail='treatment-naive at diagnosis (PMID 37857854)', fresh_frozen='fresh_or_cryopreserved_per_sample_unknown',
                     technology='10x', author_cell_type=a.reindex(key)['cell_type'].astype(str).values,
                     author_filtered=a.reindex(key)['filtered'].astype(str).values)))
    build('GSE242230', parts, 'Pancreatic ductal adenocarcinoma', {'evidence': 'PMID 37857854 treatment-naive EUS-FNB; annotations file author cell types / filtered flag'})

# ---- Hepatocellular carcinoma / ICC
if run('GSE149614'):
    meta = pd.read_csv(R('GSE149614', 'GSE149614_HCC.metadata.updated.txt.gz'), sep='\t').set_index('Cell')
    X, cells, genes = dense_to_sparse(R('GSE149614', 'GSE149614_HCC.scRNAseq.S71915.count.txt.gz'))
    m = meta.reindex(cells)
    keep = m['site'].isin(['Tumor', 'PVTT', 'Lymph']).values | m['site'].astype(str).str.contains('Lymph|PVTT|Tumor').values
    keep &= ~m['site'].astype(str).str.contains('Normal').values
    site = m['site'].astype(str).values
    parts = [(X[keep], cells[keep], genes, common(sample_id=m['sample'].values[keep], gsm='', patient_id='GSE149614:' + m['patient'].astype(str).values[keep],
              sampling_site=np.where(pd.Series(site[keep]).str.contains('Tumor'), 'primary_organ_site',
                                     np.where(pd.Series(site[keep]).str.contains('PVTT'), 'tumour_thrombus', 'metastatic_site')),
              material=site[keep], fresh_frozen='fresh', technology='10x', author_cell_type=m['celltype'].astype(str).values[keep],
              virus=m['virus'].astype(str).values[keep], tumour_stage=m['stage'].astype(str).values[keep]))]
    build('GSE149614', parts, 'Hepatocellular carcinoma', {'evidence': 'GEO metadata.updated site/patient; PMID 35933472 fresh; author-filtered count table',
                                                           'note': 'non-tumour liver dropped; PVTT kept as tumour_thrombus; LN metastasis kept'})

def nih_liver(src, sets, sample_col, cell_col, type_by_sample, keep_sample, patient_of):
    parts = []
    for pre in sets:
        X, bc, genes = read_mtx_triplet(lambda f: open(R(src, f), 'rb'), pre)
        info = pd.read_csv(R(src, pre + ('samples.txt.gz' if src == 'GSE125449' else 'Info.txt.gz')), sep='\t')
        smap = info.set_index(cell_col)[sample_col].reindex(bc).astype(str).values
        ok = np.array([keep_sample(s) for s in smap])
        parts.append((X[ok], bc[ok], genes, common(sample_id=smap[ok], gsm='', patient_id=[f'{src}:{patient_of(s)}' for s in smap[ok]],
                     sampling_site='primary_organ_site', material='primary liver cancer', fresh_frozen='fresh', technology='10x',
                     analysis_main=[type_by_sample(s) for s in smap[ok]])))
    return parts

if run('GSE125449'):
    j = pd.read_json(R('..', '..', '..', 'analysis', 'source_search_v1', 'detailed_v1', 'GSE125449.json'), typ='series')['samples']
    ctype = {s['title']: s['chars']['cancer type'] for s in j}
    tomain = lambda t: 'Hepatocellular carcinoma' if 'Hepato' in t else 'Intrahepatic cholangiocarcinoma'
    parts = nih_liver('GSE125449', ['GSE125449_Set1_', 'GSE125449_Set2_'], 'Sample', 'Cell Barcode',
                      lambda s: tomain(ctype.get(s, '')), lambda s: s in ctype, lambda s: re.search(r'_(P\d+)_', s).group(1))
    build('GSE125449', parts, 'Hepatocellular carcinoma', {'evidence': 'GEO samples.txt cell->sample; GSM cancer type; PMID 31588021 fresh',
                                                           'note': 'barcodes are not unique across Set1/Set2; cell_id includes sample'})

if run('GSE189903'):
    j = pd.read_json(R('..', '..', '..', 'analysis', 'source_search_v1', 'detailed_v1', 'GSE189903.json'), typ='series')['samples']
    ctype = {re.search(r'S\d+', s['title']).group(0): (s['chars']['cancer type'], s['chars']['tissue type']) for s in j}
    info = pd.read_csv(R('GSE189903', 'GSE189903_Info.txt.gz'), sep='\t')
    s2id = dict(zip(info['Sample'], info['S_ID']))
    tomain = lambda t: 'Hepatocellular carcinoma' if 'Hepato' in t else 'Intrahepatic cholangiocarcinoma'
    parts = nih_liver('GSE189903', ['GSE189903_'], 'Sample', 'Cell', lambda s: tomain(ctype[s2id[s]][0]),
                      lambda s: s in s2id and 'Adjacent' not in ctype[s2id[s]][1], lambda s: re.match(r'(\d+[CH])', s).group(1))
    build('GSE189903', parts, 'Hepatocellular carcinoma', {'evidence': 'GEO Info.txt cell->sample, GSM cancer/tissue type; adjacent non-tumour dropped; PMID 36476645',
                                                           'note': 'patient = case number + C (ICC) / H (HCC) of sample name: 3 ICC + 4 HCC = 7 cases as in the paper (multiregional)'})

# ---- Bladder
if run('GSE277524'):
    t, mem = tar_members(R('GSE277524', 'GSE277524_RAW.tar'))
    j = pd.read_json(R('..', '..', '..', 'analysis', 'source_search_v1', 'detailed_v1', 'GSE277524.json'), typ='series')['samples']
    info = {s['gsm']: s for s in j}
    parts = []
    for pre in sorted({re.sub(r'(barcodes|genes|features|matrix)\.(tsv|mtx)\.gz$', '', n) for n in mem}):
        X, bc, genes = read_mtx_triplet(lambda f: t.extractfile(mem[f]), pre)
        X, bc = floor_call(X, bc)  # GEO holds unfiltered matrices only
        g = pre.split('_')[0]; ch = info[g]['chars']; sample = re.sub(r'^GSM\d+_|_$', '', pre)
        rec = ch.get('tumor stage') == 'Recurrence'
        parts.append((X, bc, genes, common(sample_id=sample, gsm=g, patient_id=f'GSE277524:{sample.split("_")[0]}',
                     sampling_site='primary_organ_site', material=f'{ch.get("tissue")} urothelial tumour ({ch.get("tumor stage")})',
                     treatment_status='on_or_after_treatment' if rec else 'untreated_or_pretreatment',
                     treatment_detail=('recurrent tumour after prior resection' if rec else 'primary, treatment-naive (PMID 39702554)') + f'; procedure {ch.get("treatment")}',
                     fresh_frozen='fresh', technology='Singleron GEXSCOPE', cell_calling='QC_floor_on_unfiltered_matrix')))
    build('GSE277524', parts, 'Bladder cancer', {'evidence': 'GEO characteristics; PMID 39702554', 'note': 'P03_R1 is a ureter tumour (upper tract)'})

# ---- Esophageal
if run('GSE160269'):
    parts = []
    for frac in ('CD45pos', 'CD45neg'):
        cells_meta = pd.read_csv(R('GSE160269', f'GSE160269_{frac}_cells.txt.gz'), sep=r'\s+')
        X, cells, genes = dense_to_sparse(R('GSE160269', f'GSE160269_{frac}_UMIs.txt.gz'), sep=r'\s+', chunksize=500)
        smp = cells_meta.set_index('cell')['sample'].reindex(cells).astype(str).values
        ok = pd.Series(smp).str.endswith('T').values
        parts.append((X[ok], cells[ok], genes, common(sample_id=[f'{s}_{frac}' for s in smp[ok]], gsm='', patient_id=[f'GSE160269:{s[:-1]}' for s in smp[ok]],
                     sampling_site='primary_organ_site', material='ESCC primary tumour', fresh_frozen='fresh', technology="10x 5'" if frac == 'CD45pos' else '10x',
                     capture_scope=f'{frac}_sorted_fraction_recombined',
                     author_cell_type=cells_meta.set_index('cell')['annotated_type'].reindex(cells).astype(str).values[ok])))
    build('GSE160269', parts, 'Esophageal cancer', {'evidence': 'GEO CD45pos/CD45neg UMI tables + cell sample sheet; PMID 34489433 fresh',
                                                    'note': 'immune and non-immune fractions are separate sorted libraries: immune/non-immune ratio is not a tissue composition'})

if run('GSE203067'):
    X, bc, genes = read_mtx_triplet(lambda f: open(R('GSE203067', f), 'rb'), 'GSE203067_')
    code = {1: 'L1', 2: 'L2', 3: 'L3', 4: 'L4', 5: 'N1', 6: 'N2', 7: 'N3', 8: 'N4', 9: 'T1', 10: 'T2', 11: 'T3', 12: 'T4'}
    smp = np.array([code[int(b.rsplit('-', 1)[1])] for b in bc])
    ok = np.isin(smp, ['L1', 'L2', 'L3', 'T1', 'T2', 'T3', 'T4'])
    parts = [(X[ok], bc[ok], genes, common(sample_id=smp[ok], gsm='', patient_id=[f'GSE203067:P{s[1]}' for s in smp[ok]],
              sampling_site=np.where(pd.Series(smp[ok]).str.startswith('L'), 'metastatic_site', 'primary_organ_site'),
              material=np.where(pd.Series(smp[ok]).str.startswith('L'), 'metastatic lymph node', 'ESCC primary tumour'),
              treatment_status='untreated_or_pretreatment', treatment_detail='no preoperative treatment (PMC9982558)', technology='10x'))]
    build('GSE203067', parts, 'Esophageal cancer', {'evidence': 'GEO README barcode suffix -> sample; GSM characteristics; PMID 36709495',
                                                    'note': 'L4 is a normal lymph node and N* adjacent normal: dropped'})

# ---- Acral melanoma
if run('GSE215120'):
    t, mem = tar_members(R('GSE215120', 'GSE215120_RAW.tar'))
    parts = []
    for n in sorted(mem):
        m = re.match(r'(GSM\d+)_AM(\d)(pre|post)?_filtered_feature_bc_matrix\.h5$', n)
        if not m:
            continue
        X, bc, genes = read_10x_h5(io.BytesIO(t.extractfile(mem[n]).read()))
        post = m.group(3) == 'post'
        parts.append((X, bc, genes, common(sample_id=f'AM{m.group(2)}{m.group(3) or ""}', gsm=m.group(1), patient_id=f'GSE215120:AM{m.group(2)}',
                     sampling_site='unspecified', material='acral melanoma tumour',
                     treatment_status='on_or_after_treatment' if post else 'untreated_or_pretreatment',
                     treatment_detail='with anti-PD1 treatment' if post else 'no-treatment (GEO)', fresh_frozen='fresh', technology='10x')))
    build('GSE215120', parts, 'Acral melanoma', {'evidence': 'GEO characteristics treatment; PMID 36433984 fresh', 'note': 'cutaneous melanoma samples dropped (outside target)'})

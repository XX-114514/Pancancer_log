# Ingest the v2 supplementary sources (inventories/source_search/v2/candidates_detailed_check_v2.tsv,
# verdict=pass; TCGA-coverage expansion 2026-11-01) into one raw-count h5ad per source.
# Same rules as v1: tumour samples only (adjacent normal, normal lymph node, blood, premalignant,
# organoid and spatial samples dropped), no lineage-enriched libraries, sample-level patient /
# sampling site / treatment from GEO records and papers. QC runs afterwards (qc_supplement_v1.py).
# Usage: ingest_supplement_v2.py RAW_DIR OUT_DIR DETAIL_JSON_DIR RSCRIPT [SOURCE ...]
import sys, os, re, io, gzip, json, subprocess, tempfile
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.io
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ingest_common import (tar_members, read_mtx_triplet, read_10x_h5, dense_to_sparse, build, common)

raw, out, djson, rscript = sys.argv[1:5]
only = set(sys.argv[5:])
os.makedirs(out, exist_ok=True)
R = lambda s, *p: os.path.join(raw, s, *p)
run = lambda s: not only or s in only
gsm_info = lambda s: {x['gsm']: x for x in json.load(open(f'{djson}/{s}.json'))['samples']}
prefixes = lambda mem, rx=r'(barcodes|features|genes|matrix)\.(tsv|mtx)\.gz$': sorted({re.sub(rx, '', n) for n in mem if re.search(rx, n)})
B = lambda src, parts, main, ev: build(src, parts, main, ev, out)

# ---- Head and neck squamous cell carcinoma
if run('GSE181919'):  # dense UMI table + barcode metadata in the same column order
    md = pd.read_csv(R('GSE181919', 'GSE181919_Barcode_metadata.txt.gz'), sep='\t')
    X, cells, genes = dense_to_sparse(R('GSE181919', 'GSE181919_UMI_counts.txt.gz'), index_col=None)
    assert len(cells) == len(md)
    t = md['tissue.type'].astype(str).values
    ok = np.isin(t, ['CA', 'LN'])
    parts = [(X[ok], md.index.astype(str).values[ok], genes, common(sample_id=md['sample.id'].astype(str).values[ok], gsm='',
              patient_id='GSE181919:' + md['patient.id'].astype(str).values[ok],
              sampling_site=np.where(t[ok] == 'LN', 'metastatic_site', 'primary_organ_site'),
              material=(pd.Series(np.where(t[ok] == 'LN', 'HNSCC lymph node metastasis', 'HNSCC primary')) + ' (' + md['subsite'].astype(str).values[ok] + ', ' + md['hpv'].astype(str).values[ok] + ')').values,
              technology='10x', author_cell_type=md['cell.type'].astype(str).values[ok]))]
    B('GSE181919', parts, 'Head and neck squamous cell carcinoma', {'evidence': 'GEO barcode metadata tissue.type; CA and LN kept, NL and LP dropped'})

if run('GSE234933'):  # per-sample rds dgCMatrix; converted with R to MatrixMarket
    import tarfile
    sa = pd.read_csv(R('GSE234933', 'GSE234933_MGH_HNSCC_sample_annotation.txt.gz'), sep='\t').set_index('Sample')
    ca = pd.read_csv(R('GSE234933', 'GSE234933_MGH_HNSCC_cell_annotation.txt.gz'), sep='\t').set_index('sample_barcode')
    t = tarfile.open(R('GSE234933', 'GSE234933_MGH_HNSCC_gex_raw_counts.tar.gz'))
    parts = []
    with tempfile.TemporaryDirectory() as td:
        for m in t.getmembers():
            if not m.name.endswith('_gex_raw_counts.rds'):
                continue
            sample = os.path.basename(m.name).replace('_gex_raw_counts.rds', '')
            p = os.path.join(td, 'x.rds'); open(p, 'wb').write(t.extractfile(m).read())
            subprocess.run([rscript, '-e', f'm <- readRDS("{p}"); m <- as(m, "CsparseMatrix"); Matrix::writeMM(m, "{td}/m.mtx"); '
                            f'writeLines(rownames(m), "{td}/g.txt"); writeLines(colnames(m), "{td}/c.txt")'], check=True, capture_output=True)
            X = scipy.io.mmread(f'{td}/m.mtx').T.tocsr()
            genes = open(f'{td}/g.txt').read().split('\n')[:X.shape[1]]
            bc = np.array(open(f'{td}/c.txt').read().split('\n')[:X.shape[0]])
            key = pd.Index([b if b.startswith(sample + '_') else f'{sample}_{b}' for b in bc])
            inann = key.isin(ca.index)  # matrices may include unannotated barcodes: keep author-annotated cells
            loc = str(sa.loc[sample, 'Anatomic location of scRNA-seq specimen'])
            site = 'primary_organ_site' if loc == 'Primary' else ('metastatic_site' if 'metastasis' in loc.lower() else 'unspecified')
            parts.append((X[inann], bc[inann], genes, common(sample_id=sample, gsm='', patient_id=f'GSE234933:{sample}', sampling_site=site,
                         material=f'HNSCC {loc}; origin {sa.loc[sample, "Original anatomic site"]}; HPV {sa.loc[sample, "HPV-Status"]}',
                         fresh_frozen='fresh', technology='10x', cell_calling='author_annotation_list',
                         author_cell_type=ca.reindex(key[inann])['major_state'].astype(str).values)))
    B('GSE234933', parts, 'Head and neck squamous cell carcinoma', {'evidence': 'GEO: no cell-type enrichment; 52 fresh tumours (24 primary, 16 recurrence, 12 distant metastasis)',
                                                                    'note': 'patient = sample (paper: 51 patients for 52 samples; the shared patient is not identified in GEO)'})

if run('GSE173468'):
    info = gsm_info('GSE173468'); t, mem = tar_members(R('GSE173468', 'GSE173468_RAW.tar'))
    parts = []
    for pre in prefixes(mem):
        g = pre.split('_')[0]; ch = info[g]['chars']
        if ch.get('tissue type') != 'tumor':
            continue
        X, bc, genes = read_mtx_triplet(lambda f: t.extractfile(mem[f]), pre)
        title = info[g]['title']
        parts.append((X, bc, genes, common(sample_id=re.sub(r'^GSM\d+_|_$', '', pre), gsm=g, patient_id='GSE173468:' + ch['subject id'].replace('Patient ', ''),
                     sampling_site='metastatic_site' if 'metastatic' in title else 'primary_organ_site', material='HPV+ HNSCC', technology='10x')))
    B('GSE173468', parts, 'Head and neck squamous cell carcinoma', {'evidence': 'GEO tissue type; adjacent normal dropped'})

if run('GSE215403'):
    info = gsm_info('GSE215403'); t, mem = tar_members(R('GSE215403', 'GSE215403_RAW.tar'))
    parts = []
    for pre in prefixes(mem):
        g = pre.split('_')[0]
        X, bc, genes = read_mtx_triplet(lambda f: t.extractfile(mem[f]), pre)
        pid = re.search(r'\((P\d+)\)', info[g]['title']).group(1)
        parts.append((X, bc, genes, common(sample_id=re.sub(r'^GSM\d+_|_$', '', pre), gsm=g, patient_id=f'GSE215403:{pid}',
                     sampling_site='primary_organ_site', material='OSCC gingivo-buccal primary', technology='10x')))
    B('GSE215403', parts, 'Head and neck squamous cell carcinoma', {'evidence': 'GEO GSM titles/tissue'})

if run('GSE310797'):
    info = gsm_info('GSE310797'); t, mem = tar_members(R('GSE310797', 'GSE310797_RAW.tar'))
    parts = []
    for n in sorted(mem):
        m = re.match(r'(GSM\d+)_(.+)_count\.txt\.gz$', n)
        if not m or info[m.group(1)]['title'].startswith('SP_'):
            continue
        X, cells, genes = dense_to_sparse(gzip.open(t.extractfile(mem[n])))
        parts.append((X, cells, genes, common(sample_id=m.group(2), gsm=m.group(1), patient_id=f'GSE310797:{info[m.group(1)]["title"]}',
                     sampling_site='primary_organ_site', material='OSCC tumour', technology='10x')))
    B('GSE310797', parts, 'Head and neck squamous cell carcinoma', {'evidence': 'GEO; spatial arrays SP_* dropped', 'note': 'one sample per patient (OSCC_1-16)'})

# ---- Cutaneous melanoma
if run('GSE215120_CM'):  # cutaneous samples of the acral study (already downloaded in supplement_v1)
    v1raw = os.path.join(os.path.dirname(os.path.dirname(raw.rstrip('/'))), 'supplement_v1', 'raw')
    t, mem = tar_members(os.path.join(v1raw, 'GSE215120', 'GSE215120_RAW.tar'))
    parts = []
    for n in sorted(mem):
        m = re.match(r'(GSM\d+)_CM(\d)(_lym)?_filtered_feature_bc_matrix\.h5$', n)
        if not m:
            continue
        X, bc, genes = read_10x_h5(io.BytesIO(t.extractfile(mem[n]).read()))
        parts.append((X, bc, genes, common(sample_id=f'CM{m.group(2)}{m.group(3) or ""}', gsm=m.group(1), patient_id=f'GSE215120:CM{m.group(2)}',
                     sampling_site='metastatic_site' if m.group(3) else 'unspecified', material='cutaneous melanoma' + (' lymph node' if m.group(3) else ''),
                     treatment_status='untreated_or_pretreatment', treatment_detail='no-treatment (GEO)', fresh_frozen='fresh', technology='10x')))
    B('GSE215120_CM', parts, 'Cutaneous melanoma', {'evidence': 'GEO genotype cutaneous melanoma; PMID 36433984', 'note': 'source key GSE215120_CM keeps it apart from the acral samples of GSE215120'})

if run('GSE269936'):
    info = gsm_info('GSE269936'); t, mem = tar_members(R('GSE269936', 'GSE269936_RAW.tar'))
    parts = []
    for pre in prefixes(mem):
        g = pre.split('_')[0]; title = info[g]['title']
        X, bc, genes = read_mtx_triplet(lambda f: t.extractfile(mem[f]), pre)
        parts.append((X, bc, genes, common(sample_id=title, gsm=g, patient_id=f'GSE269936:{title.split("_")[0]}',
                     sampling_site='metastatic_site', material='metastatic melanoma', technology='10x')))
    B('GSE269936', parts, 'Cutaneous melanoma', {'evidence': 'GEO tissue metastatic melanoma', 'note': 'patient = number before "_" in sample title (several lesions per patient)'})

if run('GSE222446'):
    info = gsm_info('GSE222446'); t, mem = tar_members(R('GSE222446', 'GSE222446_RAW.tar'))
    parts = []
    for n in sorted(mem):
        g = n.split('_')[0]; ch = info[g]['chars']
        X, cells, genes = dense_to_sparse(gzip.open(t.extractfile(mem[n])))
        parts.append((X, cells, genes, common(sample_id=info[g]['title'].split(' ')[0], gsm=g, patient_id=f'GSE222446:{ch["patient"]}',
                     sampling_site='metastatic_site', material=f'melanoma biopsy site {ch.get("biopsy site")}', treatment_status='other',
                     treatment_detail='pre-ACT (TIL therapy) biopsy; prior lines not recorded in GEO', technology='10x',
                     capture_scope='whole_tissue_viable_cells')))
    B('GSE222446', parts, 'Cutaneous melanoma', {'evidence': 'GEO patient/time/biopsy site; viable-cell counts', 'note': 'counts genes x cells per sample'})

if run('GSE244983'):
    an = pd.read_csv(R('GSE244983', 'GSE244983_SingleCellAnnotations.txt.gz'), sep='\t', index_col=0)
    X, cells, genes = dense_to_sparse(R('GSE244983', 'GSE244983_RawCounts_scRNAseq.txt.gz'))
    smp = an.reindex(cells)['Sample'].astype(str).values
    parts = [(X, cells, genes, common(sample_id=smp, gsm='', patient_id=['GSE244983:' + s for s in smp], sampling_site='metastatic_site',
              material='melanoma biopsy', treatment_status=np.where(pd.Series(smp).str.contains('naive', case=False), 'untreated_or_pretreatment', 'on_or_after_treatment'),
              treatment_detail=[f'sample {s} (ICB-resistance cohort)' for s in smp], technology='10x'))]
    B('GSE244983', parts, 'Cutaneous melanoma', {'evidence': 'GEO raw counts + annotations; treatment from sample name'})

# ---- Ovarian / endometrial
if run('GSE165897'):
    ci = pd.read_csv(R('GSE165897', 'GSE165897_cellInfo_HGSOC.tsv.gz'), sep='\t').set_index('cell')
    X, cells, genes = dense_to_sparse(R('GSE165897', 'GSE165897_UMIcounts_HGSOC.tsv.gz'))
    m = ci.reindex(cells)
    post = m['treatment_phase'].astype(str).str.contains('post', case=False).values
    parts = [(X, cells, genes, common(sample_id=m['sample'].astype(str).values, gsm='', patient_id='GSE165897:' + m['patient_id'].astype(str).values,
              sampling_site='unspecified', material='HGSOC ' + m['anatomical_location'].astype(str).values,
              treatment_status=np.where(post, 'on_or_after_treatment', 'untreated_or_pretreatment'),
              treatment_detail=m['treatment_phase'].astype(str).values, technology='10x', author_cell_type=m['cell_type'].astype(str).values))]
    B('GSE165897', parts, 'Ovarian cancer', {'evidence': 'GEO cellInfo treatment_phase / anatomical_location',
                                             'note': 'sampling site per anatomical location (peritoneum, omentum, ovary...) not mapped to primary/metastatic'})

if run('GSE286057'):
    info = gsm_info('GSE286057'); t, mem = tar_members(R('GSE286057', 'GSE286057_RAW.tar'))
    parts = []
    for pre in prefixes(mem):
        g = pre.split('_')[0]
        if not info[g]['title'].endswith('_tumor'):
            continue
        X, bc, genes = read_mtx_triplet(lambda f: t.extractfile(mem[f]), pre)
        parts.append((X, bc, genes, common(sample_id=info[g]['title'], gsm=g, patient_id=f'GSE286057:{info[g]["title"].split("_")[0]}',
                     sampling_site='metastatic_site', material='HGSOC omental metastasis', treatment_status='untreated_or_pretreatment',
                     treatment_detail='treatment-naive cohort (paper)', fresh_frozen='fresh', technology='10x')))
    B('GSE286057', parts, 'Ovarian cancer', {'evidence': 'GEO tissue; only tumour (omentum majus) samples kept; paper: fresh, treatment-naive'})

if run('GSE173682'):
    info = gsm_info('GSE173682'); t, mem = tar_members(R('GSE173682', 'GSE173682_RAW.tar'))
    parts = []
    for n in sorted(mem):
        m = re.match(r'(GSM\d+)_matrix-(.+)\.mtx\.gz$', n)
        if not m or not info[m.group(1)]['title'].endswith('_RNA'):
            continue
        g, lib = m.group(1), m.group(2); ch = info[g]['chars']
        if ch['histology'] == 'GIST':
            continue
        X = scipy.io.mmread(gzip.open(t.extractfile(mem[n]))).T.tocsr()
        bc = pd.read_csv(gzip.open(t.extractfile(mem[f'{g}_barcodes-{lib}.tsv.gz'])), header=None)[0].astype(str).values
        ft = pd.read_csv(gzip.open(t.extractfile(mem[f'{g}_features-{lib}.tsv.gz'])), header=None, sep='\t')
        if ft.shape[1] >= 3:
            keepf = ft[2].astype(str).str.contains('Gene Expression').values
            X, ft = X[:, keepf], ft[keepf]
        endo = ch['tumor_site'] == 'Endometrium'
        parts.append((X, bc, ft[1].astype(str).values, common(sample_id=lib, gsm=g, patient_id=f'GSE173682:{ch["patient"]}',
                     sampling_site='primary_organ_site' if 'Primary' in ch['metastatic status'] else 'metastatic_site',
                     material=f'{ch["tumor_site"]} {ch["histology"]}', treatment_status='untreated_or_pretreatment',
                     treatment_detail='treatment-naive debulking surgery (paper)', technology='10x',
                     analysis_main='Endometrial cancer' if endo else 'Ovarian cancer')))
    B('GSE173682', parts, 'Endometrial cancer', {'evidence': 'GEO histology/tumor_site; RNA libraries only (ATAC dropped); GIST sample dropped'})

# ---- Sarcoma, lymphoma, mesothelioma
if run('GSE212527'):
    info = gsm_info('GSE212527'); t, mem = tar_members(R('GSE212527', 'GSE212527_RAW.tar'))
    ann = pd.read_csv(R('GSE212527', 'GSE212527_Cell.Annotations.txt.gz'), sep='\t')
    parts = []
    for pre in prefixes(mem):
        g = pre.split('_')[0]; title = info[g]['title']
        X, bc, genes = read_mtx_triplet(lambda f: t.extractfile(mem[f]), pre)
        cells = set(ann.loc[ann['Sample'].astype(str) == title, 'Cell_Barcode'].astype(str))
        ok = np.isin(bc, list(cells))  # matrices are unfiltered: keep author-annotated cells
        X, bc = X[ok], bc[ok]
        parts.append((X, bc, genes, common(sample_id=title, gsm=g, patient_id=f'GSE212527:{title.split("_")[0]}', cell_calling='author_annotation_list',
                     sampling_site='unspecified', material=info[g]['chars']['tissue'], technology='10x')))
    B('GSE212527', parts, 'Soft tissue sarcoma', {'evidence': 'GEO tissue (LMS / UPS)', 'note': 'SRC242 biceps + pectoralis = one patient'})

if run('GSE182434'):
    an = pd.read_csv(R('GSE182434', 'GSE182434_cell_annotation.txt.gz'), sep='\t').set_index('ID')
    X, cells, genes = dense_to_sparse(R('GSE182434', 'GSE182434_raw_count_matrix.txt.gz'))
    m = an.reindex(cells)
    ok = (m['Tissue'].astype(str) == 'DLBCL').values
    parts = [(X[ok], cells[ok], genes, common(sample_id=m['Sample'].astype(str).values[ok], gsm='', patient_id='GSE182434:' + m['Patient'].astype(str).values[ok],
              sampling_site='primary_organ_site', material='DLBCL biopsy (COO ' + m['COO'].astype(str).values[ok] + ')', technology='10x',
              author_cell_type=m['CellType'].astype(str).values[ok]))]
    B('GSE182434', parts, 'Diffuse large B-cell lymphoma', {'evidence': 'GEO cell annotation Tissue == DLBCL (follicular lymphoma and tonsil dropped)'})

if run('GSE253973'):
    info = gsm_info('GSE253973'); t, mem = tar_members(R('GSE253973', 'GSE253973_RAW.tar'))
    parts = []
    for pre in prefixes(mem):
        g = pre.split('_')[0]
        if info[g]['chars'].get('tissue') != 'Mesothelioma':
            continue
        X, bc, genes = read_mtx_triplet(lambda f: t.extractfile(mem[f]), pre)
        sid = info[g]['chars']['tissue source']
        parts.append((X, bc, genes, common(sample_id=sid, gsm=g, patient_id=f'GSE253973:{sid}', sampling_site='primary_organ_site',
                     material='pleural mesothelioma biopsy', treatment_status='untreated_or_pretreatment', treatment_detail='treatment-naive (paper)', technology='10x')))
    B('GSE253973', parts, 'Mesothelioma', {'evidence': 'GEO tissue; pleura and PDO samples dropped; paper: treatment-naive biopsies'})

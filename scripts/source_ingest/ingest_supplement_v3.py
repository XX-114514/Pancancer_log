# Ingest the v3 supplementary sources (inventories/source_search/v3/candidates_detailed_check_v3.tsv,
# verdict=pass; common-cancer sources 2026-11-01). Same rules as v1/v2: tumour samples only, no
# lineage-enriched libraries, sample-level patient / sampling site / treatment from GEO and papers.
# Usage: ingest_supplement_v3.py RAW_DIR OUT_DIR DETAIL_JSON_DIR [SOURCE ...]
import sys, os, re, io, gzip, json, tarfile, zipfile
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.io
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ingest_common import (tar_members, read_mtx_triplet, read_mm, dense_to_sparse, build, common)

raw, out, djson = sys.argv[1:4]
only = set(sys.argv[4:])
os.makedirs(out, exist_ok=True)
R = lambda s, *p: os.path.join(raw, s, *p)
run = lambda s: not only or s in only
gsm_info = lambda s: {x['gsm']: x for x in json.load(open(f'{djson}/{s}.json'))['samples']}
prefixes = lambda mem, rx=r'(barcodes|features|genes|matrix)\.(tsv|mtx)\.gz$': sorted({re.sub(rx, '', n) for n in mem if re.search(rx, n)})
B = lambda src, parts, main, ev: build(src, parts, main, ev, out)

if run('GSE117570'):  # dense per-sample tables (integer despite "processed" name); tumour samples only
    info = gsm_info('GSE117570'); t, mem = tar_members(R('GSE117570', 'GSE117570_RAW.tar'))
    parts = []
    for n in sorted(mem):
        g = n.split('_')[0]; ch = info[g]['chars']
        if 'umor' not in ch.get('tissue', ''):
            continue
        X, cells, genes = dense_to_sparse(gzip.open(t.extractfile(mem[n])))
        title = info[g]['title']
        parts.append((X, cells, genes, common(sample_id=title, gsm=g, patient_id=f'GSE117570:{title.split("_")[0]}', sampling_site='primary_organ_site',
                     material=ch.get('diagnosis', ''), technology='10x', analysis_fine='LUAD' if 'adeno' in ch.get('diagnosis', '').lower() else 'LUSC')))
    B('GSE117570', parts, 'NSCLC', {'evidence': 'GEO tissue/diagnosis; adjacent normal dropped'})

if run('GSE161277'):
    info = gsm_info('GSE161277'); t, mem = tar_members(R('GSE161277', 'GSE161277_RAW.tar'))
    parts = []
    for pre in prefixes(mem):
        g = pre.split('_')[0]
        if info[g]['chars'].get('tumor stage') != 'carcinoma':
            continue
        X, bc, genes = read_mtx_triplet(lambda f: t.extractfile(mem[f]), pre)
        title = info[g]['title']
        parts.append((X, bc, genes, common(sample_id=title.replace(' ', '_'), gsm=g, patient_id=f'GSE161277:{title.split(" ")[0]}', sampling_site='primary_organ_site',
                     material='CRC carcinoma (MSS)', treatment_status='untreated_or_pretreatment', treatment_detail='untreated (paper)', technology='10x')))
    B('GSE161277', parts, 'Colorectal cancer', {'evidence': 'GEO tumor stage; normal/adenoma/para-cancer/blood dropped; paper: four untreated MSS patients'})

if run('GSE167297'):
    info = gsm_info('GSE167297'); t, mem = tar_members(R('GSE167297', 'GSE167297_RAW.tar'))
    parts = []
    for n in sorted(mem):
        g = n.split('_')[0]; title = info[g]['title']
        if 'cancer' not in title:
            continue
        X, cells, genes = dense_to_sparse(gzip.open(t.extractfile(mem[n])))
        parts.append((X, cells, genes, common(sample_id=title, gsm=g, patient_id=f'GSE167297:{title.split("_")[0]}', sampling_site='primary_organ_site',
                     material=f'diffuse-type gastric cancer, {title.split("_")[1]} layer', technology='10x')))
    B('GSE167297', parts, 'Gastric cancer', {'evidence': 'GEO tissue; normal gastric tissue dropped'})

if run('GSE173278'):
    md = pd.read_csv(R('GSE173278', 'GSE173278_scRNAseq_all_cells_metadata.csv.gz'))
    X, bc, genes = read_mtx_triplet(lambda f: open(R('GSE173278', f), 'rb'), 'GSE173278_scRNAseq_all_cells_')
    assert len(bc) == len(md)
    ok = (md['source'].astype(str) == 'tis').values
    parts = [(X[ok], bc[ok], genes, common(sample_id=md['sample'].astype(str).values[ok], gsm='', patient_id='GSE173278:' + md['patient'].astype(str).values[ok],
              sampling_site='primary_organ_site', material='primary glioblastoma tissue (' + md['reg_stg'].astype(str).values[ok] + ')',
              technology='10x'))]
    B('GSE173278', parts, 'Diffuse glioma', {'evidence': 'GEO metadata source == tis (organoid and cell-line samples dropped)'})

if run('GSE176078'):
    t = tarfile.open(R('GSE176078', 'GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz'))
    gt = lambda n: t.extractfile(f'Wu_etal_2021_BRCA_scRNASeq/{n}')
    X = read_mm(gt('count_matrix_sparse.mtx'))
    genes = pd.read_csv(gt('count_matrix_genes.tsv'), header=None)[0].astype(str).values
    bc = pd.read_csv(gt('count_matrix_barcodes.tsv'), header=None)[0].astype(str).values
    md = pd.read_csv(gt('metadata.csv'), index_col=0).reindex(bc)
    sub = md['subtype'].astype(str).values if 'subtype' in md else ''
    parts = [(X, bc, genes, common(sample_id=md['orig.ident'].astype(str).values, gsm='', patient_id='GSE176078:' + md['orig.ident'].astype(str).values,
              sampling_site='primary_organ_site', material=['primary breast tumour ' + s for s in sub], technology='10x',
              author_cell_type=md['celltype_major'].astype(str).values if 'celltype_major' in md else ''))]
    B('GSE176078', parts, 'Breast cancer', {'evidence': 'Wu et al. 2021 count matrix + metadata (26 primary tumours)'})

if run('GSE180286'):
    info = gsm_info('GSE180286'); t, mem = tar_members(R('GSE180286', 'GSE180286_RAW.tar'))
    parts = []
    for n in sorted(mem):
        g = n.split('_')[0]; ch = info[g]['chars']
        X, cells, genes = dense_to_sparse(gzip.open(t.extractfile(mem[n])))
        ln = ch.get('tissue') == 'lymph node'
        parts.append((X, cells, genes, common(sample_id=info[g]['title'], gsm=g, patient_id=f'GSE180286:P{ch["patient id"]}',
                     sampling_site='metastatic_site' if ln else 'primary_organ_site', material=f'breast cancer {ch.get("tissue")} ({ch.get("molecular typing")})',
                     treatment_status='untreated_or_pretreatment', treatment_detail='untreated (GEO)', technology='Singleron GEXSCOPE')))
    B('GSE180286', parts, 'Breast cancer', {'evidence': 'GEO characteristics; lymph nodes are metastatic (paper title)'})

if run('GSE190888'):
    info = gsm_info('GSE190888'); t, mem = tar_members(R('GSE190888', 'GSE190888_RAW.tar'))
    parts = []
    for n in sorted(mem):
        g = n.split('_')[0]; ch = info[g]['chars']
        if 'bladder cancer' not in ch.get('disease state', ''):
            continue
        X, cells, genes = dense_to_sparse(gzip.open(t.extractfile(mem[n])), sep=',')
        rec = 'recurrent' in ch['disease state']
        parts.append((X, cells, genes, common(sample_id=info[g]['title'].replace(' ', '_'), gsm=g, patient_id=f'GSE190888:{ch["patient id"]}',
                     sampling_site='primary_organ_site', material=f'{ch["disease state"]} (grade {ch.get("grade")})',
                     treatment_status='on_or_after_treatment' if rec else 'unknown', treatment_detail='recurrent after prior resection' if rec else '',
                     technology='Singleron GEXSCOPE')))
    B('GSE190888', parts, 'Bladder cancer', {'evidence': 'GEO disease state; cystitis glandularis dropped'})

if run('GSE205506'):
    info = gsm_info('GSE205506'); t, mem = tar_members(R('GSE205506', 'GSE205506_RAW.tar'))
    parts = []
    for pre in prefixes(mem):
        g = pre.split('_')[0]; ch = info[g]['chars']
        if ch.get('genotype') != 'tumor':
            continue
        X, bc, genes = read_mtx_triplet(lambda f: t.extractfile(mem[f]), pre)
        tr = ch.get('treatment', '')
        parts.append((X, bc, genes, common(sample_id=info[g]['title'].split(':')[0], gsm=g, patient_id=f'GSE205506:{ch["subject"]}',
                     sampling_site='primary_organ_site', material=f'dMMR CRC tumour ({ch.get("tissue")})',
                     treatment_status='untreated_or_pretreatment' if tr == 'untreated' else 'on_or_after_treatment', treatment_detail=tr, technology='10x')))
    B('GSE205506', parts, 'Colorectal cancer', {'evidence': 'GEO genotype tumor/normal and treatment fields; normal mucosa dropped'})

if run('GSE207422'):
    md = pd.read_excel(R('GSE207422', 'GSE207422_NSCLC_scRNAseq_metadata.xlsx')).set_index('Sample')
    X, cells, genes = dense_to_sparse(R('GSE207422', 'GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz'))
    smp = np.array([c.rsplit('_', 1)[0] for c in cells])
    m = md.reindex(smp)
    pre = m['Resource'].astype(str).str.contains('Pre').values
    parts = [(X, cells, genes, common(sample_id=smp, gsm='', patient_id='GSE207422:' + m['Patient'].astype(str).values, sampling_site='primary_organ_site',
              material=('NSCLC ' + m['Pathology'].astype(str) + ' ' + m['Resource'].astype(str)).values,
              treatment_status=np.where(pre, 'untreated_or_pretreatment', 'on_or_after_treatment'),
              treatment_detail=(m['PD1 Antibody'].astype(str) + ' + ' + m['Chemotherapy'].astype(str) + '; response ' + m['Pathologic Response'].astype(str)).values,
              technology='BD Rhapsody WTA'))]
    B('GSE207422', parts, 'NSCLC', {'evidence': 'GEO metadata xlsx (resource pre-treatment biopsy / post-treatment surgery)'})

if run('GSE207493'):
    t = tarfile.open(R('GSE207493', 'GSE207493_scRNA-seq.tar.gz'))
    z = zipfile.ZipFile(io.BytesIO(t.extractfile('scRNA-seq.zip').read()))
    parts = []
    for s in sorted({n.split('/')[0] for n in z.namelist() if n.endswith('matrix.mtx.gz')}):
        X, bc, genes = read_mtx_triplet(lambda f: io.BytesIO(z.read(f)), f'{s}/')
        parts.append((X, bc, genes, common(sample_id=s, gsm='', patient_id=f'GSE207493:{s}', sampling_site='primary_organ_site',
                     material='ccRCC total tumour homogenate', technology='10x')))
    B('GSE207493', parts, 'Renal cell carcinoma', {'evidence': 'GEO: total ccRCC homogenate, 19 patients; scATAC dropped'})

"""Fetch selected-gene mRNA and clinical data for iAtlas-harmonised ICI cohorts from the public cBioPortal API."""
import json, urllib.request, pandas as pd, sys, os, time
API = 'https://www.cbioportal.org/api'
O = os.environ['TFINE_ROOT'] + '/bulk_v2/data'; os.makedirs(O, exist_ok=True)
studies = ['blca_iatlas_imvigor210_2017', 'mel_iatlas_liu_2019', 'mel_iatlas_gide_2019', 'mel_iatlas_riaz_nivolumab_2017', 'mel_iatlas_hugo_ucla_2016']
panels = json.load(open(O + '/../frozen_arm_panels.json'))
gep = ['CCL5','CD27','CD274','CD276','CD8A','CMKLR1','CXCL9','CXCR6','HLA-DQA1','HLA-DRB1','HLA-E','IDO1','LAG3','NKG7','PDCD1LG2','PSMB10','STAT1','TIGIT']
v1 = ['FOXP3','IL2RA','CTLA4','IKZF2','TNFRSF18','GZMB','PRF1','CTSW','NKG7','CD8A','CD8B']
genes = sorted(set(panels['Treg'] + panels['Tex'] + ['CD3D','CD3E','TRAC','GZMA','PRF1'] + gep + v1))
def req(path, body=None):
    for i in range(4):
        try:
            r = urllib.request.Request(API + path, data=None if body is None else json.dumps(body).encode(),
                                       headers={'Content-Type': 'application/json', 'Accept': 'application/json', 'User-Agent': 'generic-research-client'})
            return json.load(urllib.request.urlopen(r, timeout=120))
        except Exception as e:
            print('retry', path, e, flush=True); time.sleep(5 * (i + 1))
    raise RuntimeError(path)
gmap = {g['hugoGeneSymbol']: g['entrezGeneId'] for g in req('/genes/fetch?geneIdType=HUGO_GENE_SYMBOL', genes)}
print('genes resolved', len(gmap), 'missing', sorted(set(genes) - set(gmap)))
for s in studies:
    d = req(f'/molecular-profiles/{s}_rna_seq_mrna/molecular-data/fetch?projection=SUMMARY',
            {'entrezGeneIds': list(gmap.values()), 'sampleListId': f'{s}_all'})
    inv = {v: k for k, v in gmap.items()}
    E = pd.DataFrame([(x['sampleId'], inv[x['entrezGeneId']], x['value']) for x in d], columns=['sample', 'gene', 'value'])
    E.pivot_table(index='gene', columns='sample', values='value').to_csv(f'{O}/{s}_expr.tsv', sep='\t')
    for t in ['SAMPLE', 'PATIENT']:
        c = req(f'/studies/{s}/clinical-data?clinicalDataType={t}&projection=SUMMARY')
        C = pd.DataFrame([(x.get('sampleId', x['patientId']), x['patientId'], x['clinicalAttributeId'], x['value']) for x in c], columns=['id', 'patient', 'attr', 'value'])
        C.pivot_table(index=['id', 'patient'], columns='attr', values='value', aggfunc='first').to_csv(f'{O}/{s}_clin_{t}.tsv', sep='\t')
    print(s, 'samples', E['sample'].nunique(), 'genes', E.gene.nunique(), flush=True)

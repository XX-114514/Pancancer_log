"""Build CoVarNet inputs (frequency within major cell type, as freq_calculate) on curated T states.
Designs:
  A_original_units : same 118 units / independent_unit blocks as original T_whole_TME_primary_tumor (pre+post averaged per patient)
  B_baseline       : as A but GSE286827 post/post2 removed (pre-treatment / untreated only)
  R_replication    : GSE169246 tumour Pre (non-blood) + GSE176021 tumour; blocks = patient
Nodes: CD4 {Treg, CD4_CXCL13_Tph, CD4_memory, CD4_naive, CD4_ISG}; CD8 {Tex, CD8_GZMK_Tem, CD8_Tem_mixed, CD8_TRM, CD8_Temra, CD8_naive, CD8_ISG}.
Units need >=20 CD4 and >=20 CD8 curated cells (original rule: >=20 per measured major).
"""
import pandas as pd, numpy as np, os
W = os.environ['TFINE_ROOT']; O = W + '/covarnet_rerun'
ORIG = os.environ['PROJECT_V5_ROOT'] + '/runs/20261013_immunotherapy_atlas_v1/CM/T_whole_TME_primary_tumor/inputs/sample_metadata_private.tsv'
st = pd.read_csv(W + '/cells_phaseB_states.tsv.gz', sep='\t', index_col=0)
ob = pd.read_csv(W + '/cells_phaseA.tsv.gz', sep='\t', index_col=0, usecols=['row', 'dataset', 'sampleID', 'patient_id', 'tissue', 'timepoint', 'cohort'])
d = st.join(ob)
nodes = {'CD4': ['Treg', 'CD4_CXCL13_Tph', 'CD4_memory', 'CD4_naive', 'CD4_ISG'],
         'CD8': ['Tex', 'CD8_GZMK_Tem', 'CD8_Tem_mixed', 'CD8_TRM', 'CD8_Temra', 'CD8_naive', 'CD8_ISG']}
node_major = {n: m for m, v in nodes.items() for n in v}
d = d[d.state.isin(node_major)]
orig = pd.read_csv(ORIG, sep='\t')
print('original units', len(orig), 'matched in extracted data', orig.sampleID.isin(d.sampleID.astype(str)).sum())

def build(name, cells, meta):
    ct = pd.crosstab(cells.state, cells.sampleID.astype(str)).reindex(list(node_major)).fillna(0)
    maj = pd.Series(node_major)
    ok = [u for u in ct.columns if all(ct.loc[maj.index[maj == m], u].sum() >= 20 for m in nodes)]
    ct = ct[ok]
    fq = ct.div(ct.groupby(maj).transform('sum'))
    fq.insert(0, 'majorCluster', maj.loc[fq.index])
    fq.index.name = 'subCluster'
    meta = meta[meta.sampleID.isin(ok)].copy()
    os.makedirs(f'{O}/{name}/inputs', exist_ok=True)
    fq[['majorCluster'] + list(meta.sampleID)].to_csv(f'{O}/{name}/inputs/frequency_private.tsv', sep='\t')
    meta.to_csv(f'{O}/{name}/inputs/sample_metadata_private.tsv', sep='\t', index=False)
    print(name, 'units', len(meta), 'blocks', meta.independent_unit.nunique(), meta.dataset.value_counts().to_dict())

disc = d[d.cohort == 'discovery']
mA = orig[['sampleID', 'dataset', 'patient_block', 'independent_unit', 'patient_verified', 'timepoint']].copy()
build('A_original_units', disc, mA)
mB = mA[~mA.timepoint.astype(str).isin(['post', 'post2'])]
build('B_baseline', disc, mB)
rep = d[(d.cohort == 'replication') & (((d.dataset == 'GSE169246') & (d.timepoint == 'Pre') & (d.tissue != 'blood')) |
                                       ((d.dataset == 'GSE176021') & (d.tissue == 'tumor')))]
mR = rep.groupby('sampleID').agg(dataset=('dataset', 'first'), patient_id=('patient_id', 'first'), timepoint=('timepoint', 'first')).reset_index()
mR['patient_block'] = mR.dataset + '::' + mR.patient_id.astype(str)
mR['independent_unit'] = mR.patient_block; mR['patient_verified'] = True
build('R_replication', rep, mR[['sampleID', 'dataset', 'patient_block', 'independent_unit', 'patient_verified', 'timepoint']])

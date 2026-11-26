import zarr, numpy as np, pandas as pd, sys
Z=sys.argv[1]
g=zarr.open(Z,mode='r')
def col(name):
    a=g['obs'][name]
    if isinstance(a,zarr.hierarchy.Group) and 'categories' in a:
        cats=a['categories'][:]; codes=a['codes'][:]
        out=np.array(cats,dtype=object)[codes]; out[codes<0]=None; return out
    return a[:]
cols=['dataset','major','subcluster','capture_scope','tissue','timepoint','sampleID','patient_block','CM_composition_eligible','capture_bias_status','identity_scope','patient_identity_verified','annotation_status']
df=pd.DataFrame({c:col(c) for c in cols})
df.to_parquet('obs_small.parquet')
S=['GSE294109','GSE256326','GSE305078','GSE286827']
d=df[df.dataset.isin(S)]
pd.set_option('display.width',250); pd.set_option('display.max_rows',200)
print(d.groupby(['dataset','capture_scope','tissue']).size())
print(d.groupby(['dataset','timepoint']).size())
print(df.groupby(['dataset','capture_scope']).size())
print(d[d.major=='T'].groupby(['CM_composition_eligible']).size())
print(d.annotation_status.value_counts().head(10))
print(d.identity_scope.value_counts().head(10))

from pathlib import Path
import zipfile,pandas as pd,json
r=Path(__file__).resolve().parents[1]
with zipfile.ZipFile(r/'sources/cellphonedb_v4.1.0.zip') as z:
 t={n:pd.read_csv(z.open(n)) for n in z.namelist()}
g=t['gene_table.csv'];p=t['protein_table.csv'];co=t['complex_composition_table.csv'];inter=t['interaction_table.csv'];multi=t['multidata_table.csv'];join=g.merge(p,left_on='protein_id',right_on='id_protein');mp=join.groupby('protein_multidata_id').hgnc_symbol.apply(lambda x:set(x.dropna())).to_dict()
for key,a in co.groupby('complex_multidata_id'):mp[key]=set().union(*(mp.get(v,set()) for v in a.protein_multidata_id))
meas=pd.read_csv(r/'audit/LR_gene_measurability.tsv',sep='\t');genes=set(meas.loc[meas.measured,'gene']);rows=[]
for row in inter.itertuples():
 a=mp.get(row.multidata_1_id,set());b=mp.get(row.multidata_2_id,set());allg=a|b;rows.append({'id_cp_interaction':row.id_cp_interaction,'partnerA_genes':';'.join(sorted(a)),'partnerB_genes':';'.join(sorted(b)),'all_subunits_measured':bool(a and b and allg<=genes),'missing_subunits':';'.join(sorted(allg-genes))})
pd.DataFrame(rows).to_csv(r/'audit/LR_database_background.tsv',sep='\t',index=False);print('database rows',len(rows),'measurable',sum(x['all_subunits_measured'] for x in rows))

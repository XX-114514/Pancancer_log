from pathlib import Path
import json, ast, hashlib
import pandas as pd
from panel_review_notes import review, MAJOR_NOTES
import argparse
parser=argparse.ArgumentParser()
parser.add_argument('--input-dir',type=Path,required=True)
parser.add_argument('--result-dir',type=Path,required=True)
args=parser.parse_args()
R=args.result_dir
P=args.input_dir/'inputs'
cat=json.loads((R/'panel_catalogue.json').read_text())
cross=pd.read_csv(P/'v7_label_legacy_cluster_crosswalk.tsv',sep='\t')
initial=pd.read_csv(P/'legacy_initial_cluster_gates.tsv',sep='\t')
expanded=pd.read_csv(P/'legacy_expanded_cluster_gates.tsv',sep='\t')
sources={
'Tcell':'https://github.com/Coolgenome/TCM;https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf',
'NK':'https://pmc.ncbi.nlm.nih.gov/articles/PMC6718415/;https://pmc.ncbi.nlm.nih.gov/articles/PMC4612277/',
'Dendritic':'https://pmc.ncbi.nlm.nih.gov/articles/PMC12485386/',
'Granulocyte':'https://pmc.ncbi.nlm.nih.gov/articles/PMC4935811/',
}
rows=[]
for c in cat:
    role,action,note=review(c['major'],c['old_panel'])
    a=initial[(initial.major_celltype==c['major'])&(initial.top_subtype==c['old_panel'])]
    b=expanded[(expanded.major_celltype==c['major'])&(expanded.top_subtype==c['old_panel'])]
    ids=set(a.cluster_uid)|set(b.cluster_uid)
    mapping=cross[cross.cluster_uid.isin(ids)]
    rows.append(dict(major=c['major'],old_panel=c['old_panel'],old_genes=','.join(c['genes']),
       old_canonical=','.join(c['canonical']),proposed_role=role,action=action,review=note,
       historical_top_clusters=','.join(sorted(ids)),historical_top_v7_cells=int(mapping.n_cells.sum()),
       initial_labels=';'.join(a.final_cluster_label.dropna().astype(str)),
       expanded_rescued_clusters=','.join(b.loc[b.rescued,'cluster_uid']),
       label_count_meaning='Historical top-panel exposure; not proven misannotations. Clusters may appear under multiple panels across stages; do not sum rows.',
       evidence_source=sources.get(c['major'],'https://github.com/lambrechtslab/Pan-cancer_Lodietal'),
       evidence_boundary='Manual audit proposal informed by cited sources; not a verbatim paper panel or validated replacement classifier.',
       missing_legacy_genes=','.join(c['missing_from_legacy_object'])))
df=pd.DataFrame(rows).replace('',pd.NA)
df.to_csv(R/'all_162_panel_reviews.tsv',sep='\t',index=False,na_rep='NA')
maj=json.loads((R/'major_panel_catalogue.json').read_text())
mdf=pd.DataFrame([dict(major=k,positive_genes=','.join(v),exclusion_genes=','.join(maj['MAJOR_EXCLUDE_MARKERS'][k]),review=MAJOR_NOTES[k]) for k,v in maj['MAJOR_MARKERS'].items()])
mdf.to_csv(R/'all_17_major_panel_reviews.tsv',sep='\t',index=False,na_rep='NA')
lines=['# V7 全部旧 panel 审查表','', '17 个父类正向/排除组合与 162 个亚型候选 panel。旧定义是历史事实；修订意见是本次审计建议，不是新的细胞注释。统计列只统计历史第一候选出现范围，不代表这些细胞均错注；不同阶段同群可进入多个 panel 行，禁止跨行相加。','', '## 父类面板','', '| 父类 | 旧正向 marker | 旧排除 marker | 审查意见 |','|---|---|---|---|']
for r in mdf.to_dict('records'):lines.append('| '+' | '.join(str(r[k]) for k in ['major','positive_genes','exclusion_genes','review'])+' |')
for major,g in df.groupby('major',sort=False):
    lines+=['',f'## {major}：{len(g)} 个候选','', '| 旧 panel | 全部 marker（旧 canonical 另列） | 旧 canonical | 建议角色 | 审查意见 |','|---|---|---|---|---|']
    for r in g.to_dict('records'):lines.append('| '+' | '.join(str(r[k]) for k in ['old_panel','old_genes','old_canonical','proposed_role','review'])+' |')
(R/'all_panel_reviews.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
assert len(df)==162 and len(mdf)==17
print({'reviewed_subtype_panels':len(df),'reviewed_major_pairs':len(mdf),'roles':df.proposed_role.value_counts().to_dict()})

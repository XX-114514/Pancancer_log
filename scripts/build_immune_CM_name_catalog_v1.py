#!/usr/bin/env python3
"""Build auditable descriptive aliases; never edit original CM identities or fit outcomes."""
import argparse,csv,json,hashlib
from pathlib import Path
from datetime import datetime,timezone

OVERRIDES={
 ('whole_TME','CM01'):('GZMK–IL7R T/NK共变','GZMK–IL7R T/NK covariance'),
 ('whole_TME','CM02'):('IGHM–LTB B候选/T/NK共变','IGHM–LTB B-proposal/T/NK covariance'),
 ('T_whole_TME_primary_tumor','CM01'):('DNAJB1–CCL4 T细胞应激/趋化相关共变','DNAJB1–CCL4 stress/chemokine-associated T covariance'),
 ('T_whole_TME_primary_tumor','CM02'):('CTLA4–GZMB T细胞共变','CTLA4–GZMB T covariance'),
 ('T_whole_TME_primary_tumor','CM03'):('XIST–LTB/LEF1 T细胞混合共变','XIST–LTB/LEF1 mixed T covariance'),
 ('T_whole_TME_primary_tumor','CM04'):('GZMK–IL7R–CTSW T细胞共变','GZMK–IL7R–CTSW T covariance'),
 ('whole_TME_primary_tumor','CM01'):('IGHG3浆细胞候选共变','IGHG3 plasma-proposal covariance'),
 ('whole_TME_primary_tumor','CM02'):('免疫球蛋白–CTLA4/GZMB多谱系共变','Immunoglobulin–CTLA4/GZMB multilineage covariance'),
 ('whole_TME_primary_tumor','CM03'):('C1QB/LTB髓系–B候选/NK共变','C1QB/LTB myeloid–B-proposal/NK covariance'),
 ('whole_TME_primary_tumor','CM04'):('LTB/LEF1 T–IGHM B候选共变','LTB/LEF1 T–IGHM B-proposal covariance'),
 ('whole_TME_primary_tumor','CM05'):('CBLB/CD2 B候选–APOE髓系/NK共变','CBLB/CD2 B-proposal–APOE myeloid/NK covariance'),
 ('whole_TME_primary_tumor','CM06'):('STMN1/TYMS增殖相关淋巴细胞共变','STMN1/TYMS proliferation-associated lymphoid covariance'),
 ('whole_TME_primary_tumor','CM07'):('NKG7 T–B候选/浆细胞/髓系混合共变','NKG7 T–B-proposal/plasma/myeloid mixed covariance'),
 ('whole_TME_primary_tumor','CM08'):('IGHG4浆细胞候选–DNAJB1共变','IGHG4 plasma-proposal–DNAJB1 covariance'),
 ('whole_TME_primary_tumor','CM09'):('CREM NK–CXCL8髓系/浆细胞候选共变','CREM NK–CXCL8 myeloid/plasma-proposal covariance'),
 ('whole_TME_primary_tumor','CM10'):('NR4A1 B候选–JCHAIN浆细胞候选共变','NR4A1 B-proposal–JCHAIN plasma-proposal covariance'),
 ('whole_TME_primary_tumor','CM11'):('LAMP3 B候选–CTLA4 T/浆细胞候选共变','LAMP3 B-proposal–CTLA4 T/plasma-proposal covariance'),
 ('whole_TME_primary_tumor','CM12'):('CD2浆细胞候选–C1QB/CDKN1A髓系共变','CD2 plasma-proposal–C1QB/CDKN1A myeloid covariance'),
}

def main(a):
 m=json.loads(a.evidence.read_text());rows=[];pairs=set()
 for x in m['modules']:
  key=(x['scope'],x['factor']);assert key not in pairs;pairs.add(key)
  single=x['maximum_node_weight']>=.9
  core=[n for n in x['network_members'] if n['weight']>=.01]
  nodes=x['top_loading_nodes'][:1] if single else (core[:3] or x['top_loading_nodes'][:3])
  short=' / '.join(n['node'] for n in nodes)
  name_cn,name_en=OVERRIDES.get(key,(short+' 描述性候选',short+' descriptive candidate'))
  if x['n_network_members']==0:
   status='loading_factor_without_network_support';name_cn+='（仅载荷，未获网络支持）'
  elif single:
   status='raw_W_single_node_dominant';name_cn=short+' 原始W单节点主导因子';name_en=short+' raw-W single-node-dominant factor'
  elif not x['n_internal_edges_both_endpoints_weight_ge_1pct']:
   status='network_members_present_but_material_loading_edge_not_demonstrated';name_cn+='（载荷与网络需复核）'
  else:status='descriptive_covariance_candidate'
  flags=[]
  if x['rank_at_survey_boundary']:flags.append('rank_survey_boundary')
  if x['n_unverified_library_units']:flags.append('unverified_library_units_present')
  if x['LOSO_unestimable_count']:flags.append('LOSO_unestimable_present')
  if any(n['major'] in ['B','Plasma'] for n in nodes):flags.append('B_or_plasma_identity_proposals_require_review')
  if any(n['marker']=='XIST' for n in x['top_loading_nodes'][:3]):flags.append('sex_related_expression_requires_review')
  if x['scope']=='CAR_product':flags.append('infusion_product_not_tissue_TME')
  rows.append(dict(module_key='ICM20261014v1::'+x['scope']+'::'+x['factor'],scope=x['scope'],original_factor=x['factor'],
    proposed_name_cn=name_cn,proposed_name_en=name_en,evidence_class=status,
    naming_anchor_nodes=short,top_loading_nodes=';'.join(n['node'] for n in x['top_loading_nodes'][:5]),
    original_network_members=';'.join(n['node'] for n in x['network_members']),
    diagnostic_flags=';'.join(flags),naming_basis='raw_W_and_frozen_network;not_effective_WS_or_outcome',
    functional_name_confirmed='false',clinical_significance_confirmed='false',original_model_changed='false'))
 assert len(rows)==79 and len(pairs)==79
 a.output.mkdir(parents=True,exist_ok=True)
 with (a.output/'CM_name_catalog_v1.tsv').open('x') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
 md=['# 2026-10-14免疫图谱CM（v1）完整暂定名称','',
 '每行均保留原范围与编号。范围重叠，79行为因子实例，不能计作79个独立生态模块。名称为描述性别名，不覆盖原模型；细胞身份沿用候选标签。未获网络支持与原始W单节点主导因子明确标注。','',
 '命名使用原始W与冻结网络；实际nsNMF重建为F=W S H。原始W的集中度不等于有效基底W S只由一个节点组成。1%权重仅是本轮可读性诊断，不改变旧CM门槛或成员。','']
 for scope in dict.fromkeys(x['scope'] for x in rows):
  md+=['## '+scope,'','| 原编号 | 暂定名称 | 解释等级 |','|---|---|---|']
  for x in rows:
   if x['scope']==scope:md.append('| '+x['original_factor']+' | '+x['proposed_name_cn']+' | '+x['evidence_class']+' |')
  md.append('')
 (a.output/'CM_NAMES_CN.md').write_text('\n'.join(md))
 receipt=dict(observed_utc=datetime.now(timezone.utc).isoformat(),status='ALL79_SCOPE_QUALIFIED_DESCRIPTIVE_ALIASES_GENERATED',
  source_evidence_sha256=hashlib.sha256(a.evidence.read_bytes()).hexdigest(),
  n_factor_instances=79,n_with_original_network_members=62,n_without_original_network_members=17,
  n_raw_W_peak_ge90pct=sum(x['maximum_node_weight']>=.9 for x in m['modules']),
  n_supported_but_no_edge_with_both_raw_W_endpoints_ge1pct=sum(x['n_network_members']>0 and not x['n_internal_edges_both_endpoints_weight_ge_1pct'] for x in m['modules']),
  source_files_sha256=m['input_files'],naming_producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  learned_weights_and_H_excluded_from_public_catalog=True,original_model_or_annotation_changed=False,no_new_clinical_tests=True,
  diagnostic_thresholds_not_original_CM_acceptance_gates=True,raw_W_not_effective_WS=True,
  catalogue_sha256=hashlib.sha256((a.output/'CM_name_catalog_v1.tsv').read_bytes()).hexdigest())
 (a.output/'naming_receipt_public.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({k:v for k,v in receipt.items() if k!='source_files_sha256'}))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--evidence',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

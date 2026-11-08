import pathlib,sys,pandas as pd,json
h=pathlib.Path(sys.argv[1]);r=pathlib.Path(sys.argv[2]);assert json.loads((r/'input_audit_v1/acceptance.json').read_text())['status']=='passed'
q=pd.read_csv(r/'input_audit_v1/r3_sample_concentration.tsv',sep='\t',dtype={'cluster':str}).set_index(['pool','cluster'])
spec={
'Myeloid':[
('0','Mono_Mac_inflammatory_like','keep','inflammatory','LYZ .96 CSF1R .77 CD14 .82 C1QC .70; EREG IL1B OLR1; mixed mono/mac identity kept coarse'),
('1','cDC2_like','keep','','FCER1A .58 CD1C .62; CD1A CD1E CLEC10A DE'),
('2','Mono_classical_like','keep','','FCN1 VCAN S100A8/9/12 DE; LYZ .97 CD14 .73; C1QC .09'),
('3','Mono_FCGR3A_like','keep','','FCGR3A .79 CD14 .54 C1QC .21; LST1 CDKN1C LILRB2 DE'),
('4','Mac_metallothionein_state','keep_state','metallothionein_inflammatory','MT1G MT1H MT1X MT2A CXCL5 DE; SPP1 .81; state not identity subtype'),
('5','Myeloid_stress_state','keep_state','stress','HSPA6 BAG3 HSPB1 DNAJB1 HSPH1 DE'),
('6','Mac_FOLR2_CCL3_like','keep','inflammatory','C1QC .89 FOLR2 .49 CD14 .79; RNASE1 STAB1 SLC40A1 CCL3/4 DE'),
('7','Mac_MARCO_GPNMB_like','keep','','MARCO .65 TREM2 .84 SPP1 .72; CCL18 CHIT1 GPNMB APOC1 DE; FABP4 .13 insufficient for alveolar naming'),
('8','Myeloid_proliferating_state','keep_state','proliferating','MKI67 .68 TOP2A CENPF TYMS DE; doublet z1.43 alone not removal'),
('9','Myeloid_low_depth','keep_flag_low_depth','','729 genes normal MT3.99%; no positive DE marker; identity assigned only by flagged transfer'),
('10','Mac_CCL18_LGMN_like','keep','','C1QC .96 FOLR2 .66 TREM2 .76; CCL18 LGMN GPNMB APOE DE'),
('11','Mac_C1Q_FOLR2_like','keep','','C1QC .98 FOLR2 .56; C1QA/B/C SLC40A1 RNASE1 STAB1 DE'),
('12','Myeloid_identity_pending','set_aside_review','','293 cells doublet z1.83; SLC8A1 FMNL2 PLXDC2; no clear alternate-lineage marker combination; do not infer doublets from score alone'),
('13','Neutrophil','keep','','CSF3R .60; FCGR3B S100A8 CXCR2 FPR2 DE; low gene depth consistent with identity not automatic low quality'),
('14','Mac_glioma_TREM2_like','keep_flag_source_restricted','','Top source .757; C1QC1.0 TREM2 .95 CD14 .96 P2RY12 .05; retain contextual macrophage identity'),
('15','cDC1_like','keep','','CLEC9A .79 XCR1 .48 WDFY4 IRF8 DE; doublet z1.63 without alternate lineage insufficient for removal'),
('16','mregDC_like','keep','','LAMP3 .87 CCL19/22 CCR7 FSCN1 DE; retain established DC activation-associated phenotype'),
('17','pDC','keep','','LILRA4 .83 GZMB TCF4 IRF7 PLD4 DE'),
('18','Mast','keep','','TPSAB1 .95 KIT .75 CPA3 HPGDS DE'),
('19','Osteoclast_like','keep','','CTSK ACP5 MMP9 TCIRG1 SIGLEC15 DE; osteosarcoma enriched; COL1A1 .54 warrants context sensitivity; doublet z .94'),
('20','Non_myeloid_identity_pending','set_aside_review','','95 cells; LYZ .07 CSF1R .01 KRT8 .38; NDUFA4L2 KRT18 KRT8 DE; possible residual renal epithelial/tumour identity, not confirmed malignancy'),
('21','Myeloid_IFN_state','keep_state','IFN','ISG15 IFIT1/2/3 MX1 RSAD2; 53 cells all one source; sample flag evaluated separately')],
'Endothelial':[
('0','EC_capillary_FLT1_like','keep','','PECAM1 .92 CDH5 .74 EMCN .75 KDR .81; FLT1 RBP7 RGCC DE'),
('1','EC_tip_ESM1_like','keep','','ESM1 .77 KDR .90 PLVAP .98; ESM1 ANGPT2 KCNE3 DE'),
('2','EC_capillary_CA4_like','keep','','CA4 .41 EMCN .83 BTNL9 EDNRB RBP7 DE'),
('3','EC_venous_COL15A1_like','keep','','COL15A1 IGFBP5 PRCP DE; ACKR1 .43 PLVAP .96'),
('4','EC_arterial_like','keep','','GJA5 .58 EFNB2 .80; SEMA3G HEY1 CXCL12 DE'),
('5','EC_immediate_early_state','keep_state','immediate_early','FOSB EGR1 GADD45B SOCS3 NR4A1 DE; endothelial identity retained'),
('6','EC_venous_ACKR1_like','keep','','ACKR1 .72 PECAM1 .68 VWF .67; CLEC4G/CLEC4M approximately zero; correct stale cluster-6 sinusoidal note'),
('7','EC_venous_ACKR1_like','keep','','ACKR1 .90 VWF .97; NR2F2 SELP CPE CLU DE'),
('8','EC_inflammatory_state','keep_state','inflammatory','SELE ICAM1 VCAM1 CCL2 IL6 DE; ACKR1 .87; state separate from subtype'),
('9','EC_activated_state','keep_state','activated','HMOX1 SERPINE1 THBS1 PGF DE; avoid calling independent hypoxic identity'),
('10','EC_low_depth','keep_flag_low_depth','','792 genes MT6.23%; no positive DE marker; PECAM1 .64 VWF .58; transfer only'),
('11','EC_PGF_angiogenic_like','keep','','PGF CXCR4 PXDN MCAM DE; ESM1 .56 PECAM1 .94 CDH5 .78'),
('12','EC_sinusoidal_like_context_pending','set_aside_review','','CLEC4G .58 CLEC4M .40 MRC1 .64 CDH5 .58; mixed tissue context needs source-level evidence before liver-specific identity; preserve cells pending review'),
('13','EC_BBB_like','keep_tissue','','MFSD2A ABCG2 ABCB1 SLC2A1 DE; glioma .562; tissue-context label retained outside pan-cancer reference'),
('14','EC_proliferating_state','keep_state','proliferating','MKI67 .67 TOP2A CENPF DE; PECAM1 .97 CDH5 .87; doublet z2.34 without alternate identity insufficient for removal'),
('15','EC_venous_MHCII_like','keep','MHCII','ACKR1 .94 VWF .93; HLA-DRA DPA1 DPB1 DRB1 DE'),
('16','EC_glioma_ESM1_like','keep_flag_source_restricted','','Glioma .962 top source .755; ESM1 .62 ANGPT2 CA2 DE'),
('17','Lymphatic_EC','keep','','CCL21 .93 PDPN .69 PROX1 .64; doublet z2.46 but coherent lymphatic programme and weak alternate markers'),
('18','EC_source_restricted_state','keep_flag_source_restricted','stress_IFN','Top source .975; GADD45B DNAJB1 SOCS3 IFIT2; sample concentration checked separately')]
}
rows=[]
for pool,rs in spec.items():
 c=pd.read_csv(h/f'p0/r3_pools/{pool}/clusters.tsv',sep='\t',dtype={'leiden_r1':str}).set_index('leiden_r1')
 assert {x[0] for x in rs}==set(c.index)
 for cl,lab,act,state,reason in rs:
  qc=q.loc[(pool,cl)];flag=qc.max_source_sample_share>=.70
  if flag and act.startswith('keep') and act!='keep_flag_low_depth':
   act='keep_flag_sample_specific';reason+='; source-qualified max sample share >=0.70'
  rows.append(dict(pool=pool,cluster=cl,r3_label=lab,action=act,state=state,reason=reason,cells=int(c.loc[cl,'cells']),max_sample_share=float(qc.max_source_sample_share),review_status='working_evidence_review_not_final_publication'))
o=r/'review_v1';o.mkdir(exist_ok=False);pd.DataFrame(rows).to_csv(o/'lineage_r3_reviewed.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).groupby(['pool','action']).cells.sum().to_string())

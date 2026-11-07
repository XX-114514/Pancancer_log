"""Direct Scanpy paper-style DE, measured-universe sensitivity and response networks."""
import json, re, hashlib
from datetime import datetime, timezone
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

r=PROJECT_ROOT/'runs/20261107_V4_CoVarNet_communication_v1';out=r/'results/cytokine';out.mkdir(exist_ok=False)
a=ad.read_h5ad(r/'inputs/full_common_PRIVATE.h5ad');members=pd.read_csv(r/'audit/network_members.tsv',sep='\t')
signatures=pd.read_csv(r/'sources/signatures_one_to_one.tsv.gz',sep='\t');allgenes=set(a.var_names)
def mapping(node,major):
    if major=='Bcell':return 'B_cell'
    if major in ['Macrophage','Monocyte']:return major
    if major=='NK':return 'NK_cell'
    if major=='Mast':return 'Mast_cell'
    if major=='ILC':return 'ILC'
    for token,kind in [('Treg','Treg'),('CD8','T_cell_CD8'),('CD4','T_cell_CD4'),('gamma_delta','T_cell_gd'),('cDC1','cDC1'),('cDC2','cDC2'),('pDC','pDC'),('Neutrophil','Neutrophil')]:
        if token in node:return kind
    return None
def bh(v):return multipletests(np.nan_to_num(np.asarray(v,float),nan=1),method='fdr_bh')[1]
mappingrows=[];results=[];degproof=[];producermeans={}
for cm,g in members.groupby('cm',sort=True):
    for node in g.subCluster:
        major=str(a.obs.loc[a.obs.subCluster.eq(node),'original_major'].iloc[0]);kind=mapping(node,major)
        mask=a.obs.subCluster.eq(node);b=a[mask].copy();b.obs['comparison']=np.where(b.obs.CMT.eq(cm),'High','Low')
        nh=int(b.obs.comparison.eq('High').sum());nl=b.n_obs-nh
        producermeans[(cm,node)]=np.asarray(b[b.obs.comparison.eq('High')].X.mean(0)).ravel() if nh else None
        status='READY_CELL_LEVEL_ONLY' if kind and min(nh,nl)>=2 else 'NO_MATCHED_SIGNATURE_TYPE_OR_TOO_FEW_CELLS'
        mappingrows.append({'cm':cm,'node':node,'major':major,'signature_type':kind,'n_high_cells':nh,'n_low_cells':nl,'status':status})
        if status!='READY_CELL_LEVEL_ONLY':continue
        sc.tl.rank_genes_groups(b,groupby='comparison',groups=['High'],reference='Low',method='t-test',use_raw=False,corr_method='benjamini-hochberg',pts=False)
        de=sc.get.rank_genes_groups_df(b,group='High');de.to_csv(out/(cm+'_'+str(len(mappingrows))+'_DE_PRIVATE.tsv.gz'),sep='\t',index=False)
        detected=np.asarray((b.X>0).sum(0)).ravel()>0;measured=set(b.var_names[detected])
        up=set(de.loc[(de.logfoldchanges>np.log2(1.2))&(de.pvals_adj<.05),'names'])
        down=set(de.loc[(de.logfoldchanges< -np.log2(1.2))&(de.pvals_adj<.05),'names'])
        lib=signatures[signatures.sheet.eq(kind)];annotated=set(lib.human_gene)
        degproof.append({'cm':cm,'node':node,'n_high_cells':nh,'n_low_cells':nl,'n_detected_common_genes':len(measured),'n_up_DE':len(up),'n_down_DE':len(down)})
        for branch in ['author_unsigned_annotation_universe','directional_measured_universe']:
            universe=annotated if branch.startswith('author') else measured
            for (cytokine,symbol),t in lib.groupby(['cytokine','cytokine_symbol_mouse'],sort=True):
                for direction in ['unsigned'] if branch.startswith('author') else ['up','down']:
                    signature=set(t.human_gene if direction=='unsigned' else t.loc[t.direction.eq(direction),'human_gene'])&universe
                    query=(up if direction in ['unsigned','up'] else down)&universe;overlap=len(signature&query)
                    eligible=10<=len(signature)<=500 and len(query)>0 and len(universe)>0
                    p=float(hypergeom.sf(overlap-1,len(universe),len(signature),len(query))) if eligible else np.nan
                    results.append({'branch':branch,'cm':cm,'node':node,'signature_type':kind,'cytokine':cytokine,'cytokine_symbol_mouse':symbol,'direction':direction,
                        'universe_size':len(universe),'signature_size':len(signature),'DE_size':len(query),'overlap':overlap,'p':p,'status':'tested' if eligible else 'signature_or_query_unassessable'})
        (out/'progress.json').write_text(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),'n_CM_member_records_processed':len(mappingrows)}))
pd.DataFrame(mappingrows).to_csv(out/'CM_member_signature_eligibility.tsv',sep='\t',index=False)
pd.DataFrame(degproof).to_csv(out/'CM_member_DE_summary.tsv',sep='\t',index=False)
table=pd.DataFrame(results);table['q_within_CM_member']=np.nan;table['q_full_family']=np.nan
if len(table):
    for _,g in table.groupby(['branch','cm','node']):table.loc[g.index,'q_within_CM_member']=bh(g.p)
    for _,g in table.groupby('branch'):table.loc[g.index,'q_full_family']=bh(g.p)
    table.loc[table.status.ne('tested'),['q_within_CM_member','q_full_family']]=np.nan
table.to_csv(out/'response_full_family.tsv',sep='\t',index=False)
table.groupby(['branch','cm','cytokine']).q_within_CM_member.min().rename('minimum_member_q_descriptive_not_CM_adjusted_q').reset_index().to_csv(out/'CM_minimum_member_q_descriptive.tsv',sep='\t',index=False)
# Only unambiguous symbols are mapped. Family/engineered/uncertain cytokines remain explicit NA.
alias={'41BBL':['TNFSF9'],'APRIL':['TNFSF13'],'Adiponectin':['ADIPOQ'],'BAFF':['TNFSF13B'],'C3a':['C3'],'C5a':['C5'],'CD27L':['CD70'],'CD30L':['TNFSF8'],'CD40L':['CD40LG'],'Cardiotrophin-1':['CTF1'],'Decorin':['DCN'],'EGF':['EGF'],'FGF-basic':['FGF2'],'FasL':['FASLG'],'Flt3l':['FLT3LG'],'G-CSF':['CSF3'],'GDNF':['GDNF'],'GITRL':['TNFSF18'],'GM-CSF':['CSF2'],'HGF':['HGF'],'IFNa1':['IFNA1'],'IFNb':['IFNB1'],'IFNe':['IFNE'],'IFNg':['IFNG'],'IFNk':['IFNK'],'IFNl2':['IFNL2'],'IGF-I':['IGF1'],'IL1ra':['IL1RN'],'IL36RA':['IL36RN'],'IL12':['IL12A','IL12B'],'IL23':['IL23A','IL12B'],'IL27':['IL27','EBI3'],'IL30':['IL27'],'LIF':['LIF'],'LIGHT':['TNFSF14'],'LTA1-B2':['LTA','LTB'],'LTA2-B1':['LTA','LTB'],'Leptin':['LEP'],'M-CSF':['CSF1'],'Noggin':['NOG'],'OSM':['OSM'],'OX40L':['TNFSF4'],'Persephin':['PSPN'],'Prolactin':['PRL'],'RANKL':['TNFSF11'],'Resistin':['RETN'],'SCF':['KITLG'],'TGF-beta-1':['TGFB1'],'TL1A':['TNFSF15'],'TNFa':['TNF'],'TPO':['THPO'],'TRAIL':['TNFSF10'],'TSLP':['TSLP'],'TWEAK':['TNFSF12']}
for name in signatures.cytokine_symbol_mouse.unique():
    if re.fullmatch(r'IL\d+[a-fA-F]?',name):alias.setdefault(name,[name.upper()])
pd.DataFrame([{'cytokine_symbol_mouse':x,'human_subunits':';'.join(alias.get(x,[])),'status':'explicit_mapping_requires_RNA_not_protein_interpretation' if x in alias else 'UNRESOLVED_IDENTITY_NO_PRODUCER_EDGE'} for x in sorted(signatures.cytokine_symbol_mouse.unique())]).to_csv(out/'cytokine_producer_symbol_mapping.tsv',sep='\t',index=False)
edges=[]
for branch,qkey in [('author_unsigned_annotation_universe','q_within_CM_member'),('directional_measured_universe','q_full_family')]:
    valid=table[table.branch.eq(branch)&table[qkey].lt(.05)]
    for row in valid.itertuples():
        subunits=alias.get(row.cytokine_symbol_mouse,[])
        if not subunits:continue
        responders=set(valid.loc[valid.cm.eq(row.cm)&valid.cytokine.eq(row.cytokine),'node'])
        for prod in members.loc[members.cm.eq(row.cm),'subCluster']:
            means=producermeans.get((row.cm,prod))
            if prod in responders or means is None:continue
            available=all(s in a.var_names for s in subunits)
            expressions=[float(means[a.var_names.get_loc(s)]) if s in a.var_names else np.nan for s in subunits]
            allpass=available and all(v>.1 for v in expressions)
            for subunit,mean in zip(subunits,expressions):
                if np.isfinite(mean) and mean>.1:edges.append({'branch':branch,'cm':row.cm,'producer':prod,'cytokine':row.cytokine,'subunit':subunit,'receiver':row.node,'direction':row.direction,'mean_expression':mean,'response_q':getattr(row,qkey),'all_required_subunits_mean_gt_0_1':allpass,'meaning':'RNA production-response association; not protein secretion or intervention'})
pd.DataFrame(edges,columns=['branch','cm','producer','cytokine','subunit','receiver','direction','mean_expression','response_q','all_required_subunits_mean_gt_0_1','meaning']).to_csv(out/'production_response_network.tsv',sep='\t',index=False)
receipt={'status':'DIRECT_SCANPY_CELL_LEVEL_RESPONSE_FINISHED_PENDING_INDEPENDENT_AUDIT','UTC':datetime.now(timezone.utc).isoformat(),'scanpy':sc.__version__,
    'n_tests':len(table),'n_eligible_CM_members':len(degproof),'n_production_response_subunit_edges':len(edges),'patient_or_intervention_claim':False,
    'primary_DE':'direct Scanpy method=t-test implementation Welch; positive log2FC>log2(1.2), gene BH<.05',
    'enrichment':'scipy exact hypergeometric reconstructed, clusterProfiler not installed; 10<=signature<=500',
    'orthology':'reuse audited reciprocal one-to-one current MGI mapping, not authors historical snapshot'}
(out/'complete.json').open('x').write(json.dumps(receipt,indent=2));print(json.dumps(receipt))

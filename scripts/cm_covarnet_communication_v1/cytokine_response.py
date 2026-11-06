"""Cell-level paper-style response inference from exact sufficient statistics; not patient inference."""
from pathlib import Path
import numpy as np,pandas as pd,json
from scipy.stats import ttest_ind_from_stats,hypergeom
from statsmodels.stats.multitest import multipletests
r=Path(__file__).resolve().parents[1];o=r/'results/cytokine';o.mkdir(exist_ok=True)
z=np.load(r/'private/CMT_identity_stats.npz');ix=pd.read_csv(r/'private/CMT_identity_stats_index.tsv',sep='\t');genes=z['genes'];n=z['n'];s=z['sums'];ss=z['sq'];detect=z['det'];members=pd.read_csv(r/'audit/network_members.tsv',sep='\t');sig=pd.read_csv(r/'sources/signatures_one_to_one.tsv.gz',sep='\t');eligible=ix.CMT!='NO_H'
def bh(p):
 p=np.asarray(p);return multipletests(np.nan_to_num(p,nan=1),method='fdr_bh')[1]
def mapping(node,major):
 if major=='Bcell':return 'B_cell'
 if major=='NK':return 'NK_cell'
 if major in ['Monocyte','Macrophage']:return major
 if major=='Mast':return 'Mast_cell'
 if major=='ILC':return 'ILC'
 if 'pDC' in node:return 'pDC'
 if 'CD8' in node:return 'T_cell_CD8'
 return None
results=[];quals=[];degrows=[]
for row in members.itertuples():
 typ=mapping(row.subCluster,row.majorCluster);high=(ix.CMT==row.cm)&(ix.identity==row.subCluster);low=(ix.CMT!=row.cm)&eligible&(ix.identity==row.subCluster)
 nh=int(n[high].sum());nl=int(n[low].sum());quals.append(dict(cm=row.cm,identity=row.subCluster,signature_type=typ,n_high_cells=nh,n_low_cells=nl,status='cell_level_only' if typ and min(nh,nl)>=2 else 'identity_signature_or_cells_unqualified'))
 if not typ or min(nh,nl)<2:continue
 sh=s[high].sum(0);sl=s[low].sum(0);mh=sh/nh;ml=sl/nl;vh=np.maximum((ss[high].sum(0)-sh*sh/nh)/(nh-1),0);vl=np.maximum((ss[low].sum(0)-sl*sl/nl)/(nl-1),0)
 # Scanpy t-test uses Welch; keep that code behavior, document prose Student discrepancy.
 _,pv=ttest_ind_from_stats(mh,np.sqrt(vh),nh,ml,np.sqrt(vl),nl,equal_var=False);qv=bh(pv);lfc=np.log2((np.expm1(mh)+1e-9)/(np.expm1(ml)+1e-9));measured=(detect[high|low].sum(0)>0);up=set(genes[measured&(lfc>np.log2(1.2))&(qv<.05)]);down=set(genes[measured&(lfc<-np.log2(1.2))&(qv<.05)])
 d=pd.DataFrame({'gene':genes,'log2FC':lfc,'p':pv,'q':qv,'measurable':measured,'mean_high':mh,'mean_low':ml});d.to_csv(o/f'{row.cm}_{row.Index}_DE.tsv.gz',sep='\t',index=False)
 lib=sig[sig.sheet==typ];background=set(genes[measured]);annotated=set(lib.human_gene)
 for branch in ['author_unsigned_annotation_universe','project_directional_measured_universe']:
  U=annotated if branch.startswith('author') else background
  for (cyt,cys),t in lib.groupby(['cytokine','cytokine_symbol_mouse']):
   for direction in (['unsigned'] if branch.startswith('author') else ['up','down']):
    S=set(t.human_gene if direction=='unsigned' else t.loc[t.direction==direction,'human_gene'])&U;G=(up if direction in ['unsigned','up'] else down)&U;k=len(S&G)
    status='tested' if 10<=len(S)<=500 and len(U)>0 else 'signature_size_unqualified'
    p=float(hypergeom.sf(k-1,len(U),len(S),len(G))) if status=='tested' else np.nan
    results.append(dict(branch=branch,cm=row.cm,identity=row.subCluster,signature_type=typ,cytokine=cyt,cytokine_symbol_mouse=cys,direction=direction,universe=len(U),signature_size=len(S),DE_size=len(G),overlap=k,p=p,status=status))
a=pd.DataFrame(results);a['q_within_identity']=np.nan;a['q_full_family']=np.nan
for key,g in a.groupby(['branch','cm','identity']):a.loc[g.index,'q_within_identity']=bh(g.p)
for key,g in a.groupby('branch'):a.loc[g.index,'q_full_family']=bh(g.p)
a.loc[a.status!='tested',['q_within_identity','q_full_family']]=np.nan
a.to_csv(o/'response_full_family.tsv',sep='\t',index=False);pd.DataFrame(quals).to_csv(o/'identity_eligibility.tsv',sep='\t',index=False)
# Minimum within-subset q is the author's descriptive CM display, not a valid CM-level q.
a.groupby(['branch','cm','cytokine']).q_within_identity.min().reset_index(name='minimum_subset_q_descriptive').to_csv(o/'CM_response_minimum_q.tsv',sep='\t',index=False)
alias={'IFNg':['IFNG'],'TNFa':['TNF'],'IL12':['IL12A','IL12B'],'IL23':['IL23A','IL12B'],'IL27':['IL27','EBI3'],'TGFb1':['TGFB1']};edges=[]
for row in a[(a.branch=='project_directional_measured_universe')&(a.q_full_family<.05)].itertuples():
 cygenes=alias.get(row.cytokine_symbol_mouse,[row.cytokine_symbol_mouse.upper()]);responders=set(a[(a.cm==row.cm)&(a.cytokine==row.cytokine)&(a.branch==row.branch)&(a.q_full_family<.05)].identity)
 for prod in members[members.cm==row.cm].subCluster:
  if prod in responders:continue
  mask=(ix.CMT==row.cm)&(ix.identity==prod);nn=n[mask].sum()
  if nn<2:continue
  for cg in cygenes:
   if cg not in genes:continue
   mean=float(s[mask][:,np.where(genes==cg)[0][0]].sum()/nn)
   if mean>.1:edges.append(dict(cm=row.cm,producer=prod,cytokine=row.cytokine,subunit=cg,receiver=row.identity,mean_expression=mean,response_q=row.q_full_family,direction=row.direction))
pd.DataFrame(edges,columns=['cm','producer','cytokine','subunit','receiver','mean_expression','response_q','direction']).to_csv(o/'production_response_network.tsv',sep='\t',index=False)
(o/'complete.json').write_text(json.dumps({'status':'cell_level_response_inference_only','n_tests':len(a),'n_network_edges':len(edges),'patient_pseudobulk':'not_estimable_under_frozen_eligibility','scanpy_equivalent':'Welch_from_exact_sums_and_squares; implementation_deviation_from_direct_scanpy_call'},indent=2));print('response completed',len(a),len(edges))

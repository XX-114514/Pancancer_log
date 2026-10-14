"""Private pathway/nuisance diagnostics; no result or gene exports.
All files remain on the authorized analysis server pending publication scope.
"""
import os,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import hypergeom
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=R/'CRC_fine_pair_DIALOGUE_20261014_v1';C=R/'fullgene_receiver_target_code_v1';D=O/'private_program_interpretation_v1';assert not D.exists();D.mkdir();ref=json.loads((C/'program_reference_coverage_v1.json').read_text());sets={p['program']:set(p['shared_primary_score_genes']) for p in ref['programs']};sets['sex_linked_nuisance_check']={'XIST','TSIX','DDX3Y','EIF1AY','KDM5D','UTY','ZFY','TMSB4Y','USP9Y','PRKY','AMELY','NLGN4Y','TBL1Y','RPS4Y1'};sets['immediate_early_stress_nuisance_check']={'FOS','FOSB','JUN','JUNB','JUND','EGR1','DUSP1','IER2','IER3','ATF3','NR4A1','NR4A2','HSPA1A','HSPA1B','HSP90AA1','DDIT3','PPP1R15A'};sig=pd.read_csv(O/'discovery_active_gene_signatures.tsv',sep='\t');rows=[]
 for typ,g in sig.groupby('cell_type'):
  bg=set(pd.read_csv(O/'private_inputs'/(typ+'_discovery_genes.tsv'),header=None)[0]);allgenes=set(g.genes)
  for direction in ['all','up','down']:
   genes=allgenes if direction=='all' else set(g.loc[g.direction.eq(direction),'genes']);assert genes.issubset(bg)
   for name,s in sets.items():
    measurable=s&bg;overlap=genes&measurable;p=hypergeom.sf(len(overlap)-1,len(bg),len(measurable),len(genes));rows.append(dict(cell_type=typ,direction=direction,reference=name,background_genes=len(bg),signature_genes=len(genes),reference_measurable_genes=len(measurable),overlap_genes=len(overlap),exploratory_enrichment_p=float(p),overlap_gene_names=';'.join(sorted(overlap))))
 out=pd.DataFrame(rows);p=out.exploratory_enrichment_p.to_numpy();order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);out['exploratory_BH']=q;out.to_csv(D/'all36_private_pathway_nuisance_diagnostics.tsv',sep='\t',index=False)
 (D/'scope.json').write_text(json.dumps(dict(status='PRIVATE_PROGRAM_DIAGNOSTICS_COMPLETE_NOT_MECHANISM_VALIDATION',unweighted_signature_overlap=True,reference_sets='four previously frozen Hallmark/common dictionaries plus two explicit nuisance sets',cell_type_labels_frozen=True,genes_not_selected_using_clinical_or_holdout=True,results_retained_server_only=True,publication_scope_pending_user_authorization=True),indent=2)+'\n');print('PRIVATE_PROGRAM_DIAGNOSTICS_COMPLETE_RESULTS_REMAIN_SERVER')
if __name__=='__main__':main()

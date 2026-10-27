"""Additional T-lineage RNA proposals, with measured competing-core gates.

These rules define a working RNA evidence tier. They do not establish purity,
ontogeny, function, alpha/beta versus gamma/delta identity, or malignancy.
"""
import numpy as np

FAMILIES=['T_NK','B_Plasma','Myeloid','Mast','Endothelial','Fibroblast','Mural',
          'Epithelial','Glial','Astrocyte','Neuron','Erythroid']
CORE_GENES={
 'T_NK':['CD3D','CD3E','TRAC','TRBC1','TRBC2','TRDC','TRGC1','TRGC2','LCK','LAT','ZAP70','BCL11B'],
 'B_Plasma':['CD79A','CD79B','MS4A1','CD19','CD22','PAX5','MZB1','JCHAIN','DERL3','TNFRSF17'],
 'Myeloid':['CSF1R','SPI1','LST1','LYZ','AIF1','CLEC9A','XCR1','CD1C','FCER1A','LILRA4','CLEC4C','IL3RA','TCF4','CSF3R','FCGR3B','CEACAM8','CXCR2','MMP8','CAMP'],
 'Mast':['KIT','TPSAB1','TPSB2','CPA3','HDC','MS4A2'],
 'Endothelial':['PECAM1','CDH5','CLDN5','VWF','ESAM','PROX1','FLT4','CCL21','PDPN'],
 'Fibroblast':['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA'],
 'Mural':['RGS5','CSPG4','PDGFRB','MCAM','MYH11','CNN1'],
 'Epithelial':['EPCAM','CDH1','KRT7','KRT8','KRT18','KRT19','KRT5','KRT14','SFN'],
 'Glial':['OLIG1','OLIG2','SOX10','MBP','PLP1','MOG','MOBP','MAG'],
 'Astrocyte':['GFAP','AQP4','SLC1A3'],
 'Neuron':['SNAP25','SYT1','RBFOX3','STMN2'],
 # Preserve the complete previously audited anchored erythroid support;
 # missing new quantitative anchors never become a negative result.
 'Erythroid':[],
}
REASONS=['not_baseline_unresolved','additional_T_RNA_core_reference_local_concordance',
 'TCR_complex_or_chain_or_signaling_not_complete','competing_lineage_core_detected',
 'old_competing_family_core_measurement_incomplete','reference_family_not_unique_T_NK',
 'low_reference_not_T_or_unassigned','source_local_group_not_T_or_T_NK',
 'without_fitted_source_group','zero_native_RNA']

def evaluate(detected,genes,measured,old_native_mask,reference_mask,low_reference_T,
             local_T_context,fit,positive,baseline_unresolved):
    d=np.asarray(detected,bool);measured=np.asarray(measured,bool);genes=list(genes)
    if d.ndim!=2 or d.shape[1]!=len(genes) or measured.shape!=(len(genes),) or len(set(genes))!=len(genes):raise ValueError('Gene axes differ')
    n=len(d);old=np.asarray(old_native_mask);ref=np.asarray(reference_mask)
    arrays=[old,ref,*[np.asarray(x) for x in [low_reference_T,local_T_context,fit,positive,baseline_unresolved]]]
    if any(a.shape!=(n,) for a in arrays):raise ValueError('Cell axes differ')
    if old.dtype.kind not in 'iu' or ref.dtype.kind not in 'iu' or (old>=1<<len(FAMILIES)).any() or (ref>=1<<len(FAMILIES)).any():raise ValueError('Unknown family bit')
    if d[:,~measured].any():raise ValueError('Unmeasured detections fabricated')
    by={g:j for j,g in enumerate(genes)}
    def g(name):return d[:,by[name]] if name in by and measured[by[name]] else np.zeros(n,bool)
    def count(names):return sum((g(v).astype(np.uint8) for v in names),np.zeros(n,np.uint8))
    alpha_chain=g('TRAC')|g('TRBC1')|g('TRBC2')
    gamma_chain=g('TRDC')&(g('TRGC1')|g('TRGC2'))
    tcore=g('CD3D')&g('CD3E')&(alpha_chain|gamma_chain)&(count(['LCK','LAT','ZAP70','BCL11B'])>=1)
    bcore=((g('CD79A')|g('CD79B'))&(count(['MS4A1','CD19','CD22','PAX5'])>=1))|(g('CD79A')&g('CD79B'))
    plasma=g('MZB1')&(count(['JCHAIN','DERL3','TNFRSF17'])>=1)
    myeloid=(g('CSF1R')|g('SPI1'))&(count(['CSF1R','SPI1','LST1','LYZ','AIF1'])>=3)
    dc=(g('CLEC9A')&g('XCR1'))|(g('CD1C')&g('FCER1A'))|((g('LILRA4')|g('CLEC4C'))&(count(['LILRA4','CLEC4C','IL3RA','TCF4'])>=3))
    granulocyte=(g('CSF3R')|g('FCGR3B'))&(count(['CEACAM8','CXCR2','MMP8','CAMP'])>=1)
    cores={
     'B_Plasma':bcore|plasma,'Myeloid':myeloid|dc|granulocyte,
     'Mast':(count(['TPSAB1','TPSB2','CPA3'])>=1)&(count(CORE_GENES['Mast'])>=2),
     'Endothelial':(count(['PECAM1','CDH5','CLDN5','VWF','ESAM'])>=3)|((g('PROX1')|g('FLT4'))&(count(['PROX1','FLT4','CCL21','PDPN'])>=3)),
     'Fibroblast':count(CORE_GENES['Fibroblast'])>=3,
     'Mural':(g('RGS5')|g('CSPG4')|g('MYH11'))&(count(CORE_GENES['Mural'])>=3),
     'Epithelial':count(CORE_GENES['Epithelial'])>=3,
     'Glial':(count(['OLIG1','OLIG2','SOX10'])>=1)&(count(CORE_GENES['Glial'])>=3),
     'Astrocyte':count(CORE_GENES['Astrocyte'])>=2,
     'Neuron':count(CORE_GENES['Neuron'])>=2,
     'Erythroid':(old&(1<<FAMILIES.index('Erythroid')))!=0,
    }
    competitor=np.zeros(n,np.uint16);unknown=np.zeros(n,np.uint16);availability={}
    for name,values in cores.items():
        bit=1<<FAMILIES.index(name);competitor|=values.astype(np.uint16)*bit
        complete=all(v in by and measured[by[v]] for v in CORE_GENES[name])
        availability[name]=complete
        if not complete:unknown|=np.asarray((old&bit)!=0,np.uint16)*bit
    reason=np.zeros(n,np.uint8);target=np.asarray(baseline_unresolved,bool)
    # Apply all gates without letting later ones overwrite earlier reasons.
    pending=target.copy()
    for code,ok in [(9,positive),(2,tcore),(3,competitor==0),(4,unknown==0),
                    (5,ref==1),(6,low_reference_T),(8,fit),(7,local_T_context)]:
        failed=pending&~np.asarray(ok,bool);reason[failed]=code;pending[failed]=False
    reason[pending]=1
    if np.any(pending&((old&1)==0)):raise ValueError('TCR core unexpectedly lacks frozen T_NK evidence')
    return {'review_code':reason,'proposed_T':pending,'TCR_core':tcore,'alpha_chain_detected':alpha_chain,
        'gamma_chain_detected':gamma_chain,'competing_core_mask':competitor,'unassessed_old_competitor_mask':unknown,
        'competing_core_measurement_complete':availability}

"""Native quantities for context review and exact-group descriptive markers.

Panel membership is not a classifier. Full RNA totals remain the denominator;
unmeasured/ambiguous features are missing, not biological zeros.
"""
import re
import numpy as np
import scipy.sparse as sp

CONTEXT_PANELS={
 'T_receptor':['CD3D','CD3E','CD247','TRAC','TRBC1','TRBC2','LCK','LAT','ZAP70','BCL11B'],
 'T_identity_and_state':['CD4','CD8A','CD8B','FOXP3','IL2RA','CTLA4','TIGIT','CCR7','IL7R','TCF7','LEF1','GZMK','GZMB','PRF1','PDCD1','LAG3','HAVCR2','TOX','TRDC','TRGC1','TRGC2'],
 'NK_context':['NKG7','GNLY','KLRD1','NCR1','KLRF1','FCGR3A'],
 'B_receptor':['CD19','CD22','MS4A1','CD79A','CD79B','CD37','CD74','PAX5','BANK1'],
 'Plasma_secretory':['MZB1','JCHAIN','SDC1','DERL3','XBP1','TNFRSF17','IGHG1','IGHM','IGKC','IGLC2'],
 'Myeloid_context':['LST1','LYZ','FCER1G','CTSS','AIF1','TYROBP','SPI1','CSF1R','CD68','FCN1','C1QA','C1QB','C1QC','APOE','FOLR2','TREM2','MRC1'],
 'DC_context':['FCER1A','CD1C','CLEC10A','CLEC9A','XCR1','BATF3','LILRA4','IL3RA','CLEC4C','TCF4','PLD4','TSPAN13','CCDC50','AXL','SIGLEC6'],
 'Granulocyte_context':['CSF3R','FCGR3B','S100A8','S100A9','CEACAM8','CXCR2','MMP8','CAMP','LCN2','S100A12'],
 'Mast_basophil_context':['KIT','TPSAB1','TPSB2','CPA3','HDC','MS4A2','CLC','GATA2'],
 'Vascular_context':['PECAM1','CDH5','CLDN5','VWF','EMCN','ESAM','ERG','KDR','FLT1','PLVAP','ACKR1','PROX1','FLT4','PDPN','CCL21','LYVE1','MMRN1'],
 'Fibroblast_context':['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA','CXCL12'],
 'Mural_contractile_context':['RGS5','CSPG4','PDGFRB','MCAM','ACTA2','TAGLN','MYH11','CNN1'],
 'Epithelial_context':['EPCAM','CDH1','KRT7','KRT8','KRT18','KRT19','KRT5','KRT14','SFN'],
 'Neural_glial_context':['OLIG1','OLIG2','SOX10','MBP','PLP1','MOG','MOBP','MAG','GFAP','AQP4','SLC1A3','ALDOC','SOX9','SNAP25','SYT1','RBFOX3','STMN2'],
 'Erythroid_progenitor_context':['HBB','HBA1','HBA2','AHSP','SLC4A1','GYPA','CD34','SPINK2','MEIS1'],
 'Generic_RNA_state_context':['IFIT1','ISG15','MX1','IFIT3','MKI67','TOP2A','STMN1','UBE2C'],
}
CONTEXT_GENES=sorted({g for panel in CONTEXT_PANELS.values() for g in panel})

def canonical_native(raw):
    if not sp.issparse(raw):raise ValueError('Sparse native counts required')
    if (not np.isfinite(raw.data).all() or (raw.data<0).any() or
        (raw.data>2**53).any() or not np.equal(raw.data,np.floor(raw.data)).all()):
        raise ValueError('Native finite nonnegative exact-count domain required')
    x=raw.astype(np.float64).tocsr();x.sum_duplicates();x.eliminate_zeros()
    if (x.data>2**53).any():raise ValueError('Coalesced slot exceeds exact-count domain')
    return x

def context_values(raw,physical_columns):
    columns=list(physical_columns);selected=[int(j) for j in columns if j is not None]
    if len(selected)!=len(set(selected)):raise ValueError('Context symbols collide on physical slots')
    out=np.full((raw.shape[0],len(columns)),np.nan,np.float64);measured=np.asarray([j is not None for j in columns],bool)
    if selected:out[:,measured]=raw[:,selected].toarray()
    return out,measured

def update_group_stats(stat,raw,totals,codes):
    totals=np.asarray(totals,np.float64);codes=np.asarray(codes)
    if raw.shape[0]!=len(totals) or codes.shape!=totals.shape or raw.shape[1]!=stat['cp_sum'].shape[1]:raise ValueError('Group/gene/native axes differ')
    if codes.dtype.kind not in 'iu' or (codes<-1).any() or (codes>=len(stat['n'])).any():raise ValueError('Unregistered group code')
    if not np.isfinite(totals).all() or (totals<0).any() or (totals>2**53).any():raise ValueError('Invalid complete native RNA totals')
    fit=codes>=0
    if (totals[fit]<=0).any():raise ValueError('Fitted group cannot include zero-RNA rows')
    if np.any(np.asarray(raw.sum(axis=1)).ravel()>totals+1e-6):raise ValueError('Measured gene subset exceeds full native RNA')
    if not fit.any():return
    groups,local=np.unique(codes[fit],return_inverse=True);x=raw[fit].astype(np.float64)
    cp=x.multiply((10000/totals[fit])[:,None]).tocsr();log=cp.copy();log.data=np.log1p(log.data)
    indicator=sp.csr_matrix((np.ones(len(local)),(local,np.arange(len(local)))),shape=(len(groups),len(local)))
    stat['n'][groups]+=np.bincount(local,minlength=len(groups))
    for name,values in [('cp_sum',cp),('log_sum',log),('detected',x>0)]:stat[name][groups]+=(indicator@values).toarray()

def expression_candidates(genes,own_cp,own_det,other_cp,other_det,minimum_detection=.1,limit=20):
    """Same frozen descriptive ranking, plus explicit no-comparator fallback."""
    genes=list(genes);own_cp=np.asarray(own_cp);own_det=np.asarray(own_det)
    eligible=np.asarray([not re.match(r'^(MT-|RPL|RPS|ENSG|IG[HKL][VDJC])',g) for g in genes])&(own_det>=minimum_detection)&(own_cp>0)
    effect=None
    if other_cp is not None:
        other_cp=np.asarray(other_cp);other_det=np.asarray(other_det);effect=np.log2((own_cp+.1)/(other_cp+.1))
        keep=eligible&(own_det>other_det)&(effect>0);score=effect*np.maximum(own_det-other_det,0)
        order=sorted(np.flatnonzero(keep),key=lambda j:(-score[j],genes[j]))[:limit]
        kind='own_subset_expression_specificity_not_inferential_DEG'
    else:order=[];kind='own_subset_expression_without_comparator_not_inferential_DEG'
    # Absence of a positive specificity contrast does not remove native own
    # expression. Mark the fallback tier rather than pretending it is a DEG.
    if not order:
        order=sorted(np.flatnonzero(eligible),key=lambda j:(-own_cp[j]*own_det[j],genes[j]))[:limit]
        kind='own_subset_expression_without_positive_specificity_not_inferential_DEG' if other_cp is not None else kind
    return [{'gene':genes[j],'native_mean_CP10K':float(own_cp[j]),'native_detection_fraction':float(own_det[j]),
             'comparator_detection_fraction':float(other_det[j]) if other_det is not None else None,
             'descriptive_log2_mean_ratio':float(effect[j]) if effect is not None else None,'marker_evidence_type':kind} for j in order]

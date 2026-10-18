"""Frozen axis designs; no outcomes, identities or learned names in source."""
import numpy as np
from NSCLC_immune_RNA_MCP_cohort_math_v1 import design
from NSCLC_frozen_axis_block_math_v1 import geometry

SENSITIVITIES=['unadjusted','clinical_RNA_counts_no_fractions','main_plus_total_communication']


def controls(a,j):
    s=int(a['sender_index'][j]);r=int(a['receiver_index'][j]);C=np.c_[a['clinical_controls'],a['global_quality'],np.log1p(a['unit_cells'][:,s]),np.log1p(a['unit_cells'][:,r]),a['unit_quality'][:,s],a['unit_quality'][:,r],a['within_major_fraction'][:,s],a['within_major_fraction'][:,r]]
    names=a['clinical_control_names'].tolist()+['logN_all_immune','all_immune_mean_logRNA','all_immune_mean_loggenes','logN_sender','logN_receiver','sender_meanlogRNA','sender_meanloggenes','receiver_meanlogRNA','receiver_meanloggenes','sender_within_major_fraction','receiver_within_major_fraction']
    return C,names


def eligibility(a,j):
    C,names=controls(a,j);y=a['standardized_axis_score'][:,j];valid=a['endpoint_measurable'][:,j]&np.isfinite(np.c_[a['response'],y,C]).all(1)
    return valid,C,names


def prepare(a,j,mode='main'):
    valid,C,names=eligibility(a,j);x=a['response'][valid];y=a['standardized_axis_score'][valid,j];C=C[valid];blocks=a['profile_blocks'][valid];why='ESTIMABLE'
    if mode=='unadjusted':C=np.empty((len(x),0));names=[]
    elif mode=='clinical_RNA_counts_no_fractions':C=C[:,:-2];names=names[:-2]
    elif mode=='main_plus_total_communication':C=np.c_[C,a['log1p_mean_observed_normalized_LRscore'][valid],a['captured_author_types'][valid]];names+=['total_observed_LR_background','captured_author_types']
    elif mode!='main':raise ValueError('UNDECLARED_SENSITIVITY')
    if min(int((x==0).sum()),int((x==1).sum()))<10:return None,None,'FEWER_THAN10_PER_RESPONSE',valid
    model,why=design(x,y,C,names,'clinical')
    if model is None:return None,None,why,valid
    if model['df']<20:return None,None,'FEWER_THAN20_RESIDUAL_DF',valid
    try:g=geometry(model,blocks)
    except ValueError as e:return model,None,str(e),valid
    return model,g,'ESTIMABLE',valid

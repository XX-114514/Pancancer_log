"""Frozen complete function-family designs; no scientific outputs embedded."""
import hashlib
import numpy as np
from NSCLC_immune_RNA_MCP_cohort_math_v1 import design
from NSCLC_immune_RNA_MCP_cohort_audit_math_v1 import prepare as independent_design
from NSCLC_frozen_axis_block_math_v2 import geometry
from NSCLC_frozen_axis_block_audit_math_v1 import geometry as independent_geometry

PARTITIONS = ['DISCOVERY', 'HELDOUT']
QUESTIONS = ['clinical_response', 'axis_given_response', 'sender_given_receiver_and_response']
SENSITIVITIES = ['clinical_agent_only', 'full_without_parent_fractions']


def seed(slot, stream):
    return int.from_bytes(hashlib.sha256(('2026101951|' + slot + '|' + stream).encode()).digest()[:8], 'little')


def prepare(a, k, f, partition, question, mode='main', independent=False):
    C = a['controls'][:, k].copy()
    names = a['control_names'].tolist()
    x = (a['response'] if question == 'clinical_response' else a['sender_expression_score'][:, k] if question == 'sender_given_receiver_and_response' else a['axis_score'][:, k])
    if mode == 'clinical_agent_only':
        C = C[:, :int(a['base_control_count'])]
        names = names[:int(a['base_control_count'])]
    elif mode == 'full_without_parent_fractions':
        C = C[:, :-2]
        names = names[:-2]
    elif mode != 'main':
        raise ValueError('UNDECLARED_SENSITIVITY')
    if question in ['axis_given_response', 'sender_given_receiver_and_response']:
        C = np.column_stack([C, a['response']])
        names.append('original_pCR_or_MPR')
        if question == 'sender_given_receiver_and_response':
            C = np.column_stack([C, a['receiver_expression_score'][:, k]])
            names.append('receiver_receptor_expression')
    elif question != 'clinical_response':
        raise ValueError('UNDECLARED_QUESTION')
    y = a['function_score'][k, :, f]
    valid = ((a['partition'] == partition) & a['endpoint_measurable'][:, k]
             & np.isfinite(np.column_stack([x, y, C, a['response']])).all(1))
    if not a['function_eligible'][k, f]:
        return None, None, 'FUNCTION_GENE_OR_VARIANCE_COVERAGE_INSUFFICIENT', valid
    response = a['response'][valid]
    if min(int((response == 0).sum()), int((response == 1).sum())) < 10:
        return None, None, 'FEWER_THAN10_PER_RESPONSE', valid
    m, status = (independent_design if independent else design)(x[valid], y[valid], C[valid], names, 'clinical')
    if m is None:
        return None, None, status, valid
    if m['df'] < 20:
        return m, None, 'FEWER_THAN20_RESIDUAL_DF', valid
    try:
        g = (independent_geometry if independent else geometry)(m, a['profile_blocks'][valid])
    except ValueError as e:
        return m, None, str(e), valid
    return m, g, 'ESTIMABLE', valid


def slots(a):
    for k in range(len(a['scope_ids'])):
        for f, reference in enumerate(a['references']):
            for partition in PARTITIONS:
                for question in QUESTIONS:
                    name = '%s|%s|%s|%s' % (a['scope_ids'][k], reference, partition, question)
                    yield name, k, f, partition, question

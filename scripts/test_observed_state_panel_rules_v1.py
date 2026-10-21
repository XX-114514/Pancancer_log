"""Synthetic assay mapping checks, without human inputs."""
import qualify_observed_state_panel_rules_v1 as q


def main():
    full=set().union(*(q.leaves(v) for v in list(q.RNA_R.values())+list(q.RNA_E.values())))
    assert q.describe(q.AND(q.RNA_R['R_LAG3'],q.RNA_E['E_cycling']),full)['complete_formula_measurement']
    partial=full-{'TRAC'}
    row=q.describe(q.AND(q.RNA_R['R_LAG3'],q.RNA_E['E_cycling']),partial)
    assert row['positive_evidence_route_measurable'] and not row['complete_formula_measurement']
    assert not q.support(q.RNA_R['R_LAG3'],full-{'IL2RA'})
    assert q.support(q.RNA_RA['R_LAG3'],(full-{'IL2RA'})|{'CTLA4'})
    assert not q.support(q.RNA_E['E_cycling'],full-{'TOP2A','UBE2C'})
    assert 'LAG3' in q.canonical_proteins(['CD223'])
    assert 'MKI67' in q.canonical_proteins(['Ki-67'])
    assert 'PDCD1' in q.canonical_proteins(['CD279'])
    assert 'IL7R' in q.canonical_proteins(['CD127'])
    assert 'TNFRSF9' not in q.canonical_proteins(['CD134','OX40','GITR'])
    assert 'TCF7' not in q.canonical_proteins(['TCF1','CD45RO'])
    assert 'LAG3' not in q.canonical_proteins(['prob.LAG3','LAG3.prob'])
    assert 'TCF7' in q.canonical_proteins(['TCF1/TCF7'])
    # Complementary Ki67 protein evidence never becomes exact RNA measurement.
    result=q.qualify(dict(source='synthetic',panel_key='P',assay='protein',measured_marker_names=['CD3','CD4','FoxP3','CD8a','CD223','Ki67'],input_data_column_authority='synthetic',proof_SHA256='synthetic'))
    chosen=next(x for x in result['pair_opportunities'] if x['R_state']=='R_LAG3' and x['E_state']=='E_cycling')
    assert chosen['positive_evidence_route_measurable'] and not chosen['RNA_formula_measurement_equivalent']
    print('PASS_14_SYNTHETIC_FULL_PARTIAL_ALTERNATE_PARENT_ALIAS_AND_PROTEIN_BOUNDARY_CHECKS')


if __name__=='__main__':main()

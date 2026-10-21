"""Assay-only observability of fixed RNA rules and complementary protein proxies.

No expression values, cell records, weights, coordinates or clinical rows.
An incomplete OR rule can support observed positives without supporting a
complete negative classification. Protein proxies never become RNA identities.
"""
import argparse
from itertools import product
import json
from pathlib import Path
import re


def AND(*items): return {'all': list(items)}
def OR(*items): return {'any': list(items)}


T = OR('CD3D','CD3E','TRAC')
R = AND(T,'FOXP3','IL2RA')
RA = AND(T,'FOXP3',OR('IL2RA','CTLA4'))
E = AND(T,'CD8A')
CYCLE = AND('MKI67',OR('TOP2A','UBE2C'))
RNA_R = {'R_TNFRSF9': AND(R,'TNFRSF9'), 'R_LAG3': AND(R,'LAG3'), 'R_cycling': AND(R,CYCLE)}
RNA_RA = {'R_TNFRSF9': AND(RA,'TNFRSF9'), 'R_LAG3': AND(RA,'LAG3'), 'R_cycling': AND(RA,CYCLE)}
RNA_E = {'E_dysfunction': AND(E,'PDCD1',OR('TOX','HAVCR2')),
         'E_cycling': AND(E,CYCLE), 'E_memory': AND(E,OR('TCF7','IL7R'))}
PR = AND('CD3','CD4','FOXP3')
PE = AND('CD3','CD8')
PROTEIN_R = {'R_TNFRSF9': AND(PR,'TNFRSF9'), 'R_LAG3': AND(PR,'LAG3'), 'R_cycling': AND(PR,'MKI67')}
PROTEIN_E = {'E_dysfunction': AND(PE,'PDCD1',OR('TOX','HAVCR2')),
             'E_cycling': AND(PE,'MKI67'), 'E_memory': AND(PE,OR('TCF7','IL7R'))}

# Exact antigen names only. OX40/CD134, GITR and CD45RO are deliberately not
# substituted for TNFRSF9/CD137 or TCF7/IL7R. TCF1 alone needs a source binding.
ALIASES = {
    'CD3':['CD3'], 'CD4':['CD4'], 'CD8':['CD8','CD8a','CD8A'], 'FOXP3':['FOXP3','FoxP3'],
    'IL2RA':['IL2RA','CD25'], 'TNFRSF9':['TNFRSF9','CD137','4-1BB'],
    'LAG3':['LAG3','LAG-3','CD223'], 'MKI67':['MKI67','Ki67','Ki-67'],
    'TOP2A':['TOP2A'], 'UBE2C':['UBE2C'], 'PDCD1':['PDCD1','PD1','PD-1','CD279'],
    'TOX':['TOX'], 'HAVCR2':['HAVCR2','TIM3','TIM-3','CD366'],
    'TCF7':['TCF7','TCF1/TCF7','TCF-1 (TCF7)'], 'IL7R':['IL7R','CD127']}


def normalize(name):
    return re.sub('[^a-z0-9]', '', name.casefold())


def canonical_proteins(names):
    lookup = {}
    for target, aliases in ALIASES.items():
        for alias in aliases:
            key = normalize(alias)
            assert key not in lookup or lookup[key] == target
            lookup[key] = target
    return {lookup[normalize(n)] for n in names if normalize(n) in lookup}


def leaves(rule):
    if isinstance(rule,str): return {rule}
    assert set(rule) in [{'all'},{'any'}]
    return set().union(*(leaves(x) for x in next(iter(rule.values()))))


def support(rule, measured):
    if isinstance(rule,str): return rule in measured
    if 'all' in rule: return all(support(x,measured) for x in rule['all'])
    return any(support(x,measured) for x in rule['any'])


def describe(rule, measured):
    required = leaves(rule)
    return {'complete_formula_measurement': required <= measured,
            'positive_evidence_route_measurable': support(rule,measured),
            'unmeasured_terms': sorted(required-measured),
            'absence_of_cells_or_zero_expression_inferred': False}


def qualify(panel):
    assay = panel['assay']
    assert assay in ['RNA','protein']
    raw = panel['measured_marker_names']
    measured = set(raw) if assay == 'RNA' else canonical_proteins(raw)
    rr, ee = (RNA_R,RNA_E) if assay == 'RNA' else (PROTEIN_R,PROTEIN_E)
    rows = []
    for a,b in product(rr,ee):
        result = describe(AND(rr[a],ee[b]),measured)
        row = {'source':panel['source'],'panel_key':panel['panel_key'],'assay':assay,
            'R_state':a,'E_state':b,**result,'input_data_column_authority':panel['input_data_column_authority'],
            'new_effect_or_test':False,'clinical_validation_qualified':False}
        if assay == 'RNA':
            row['alternate_regulatory_parent'] = describe(AND(RNA_RA[a],ee[b]),measured)
        else:
            row['RNA_formula_measurement_equivalent'] = False
            row['regulatory_CD25_antigen_measured'] = 'IL2RA' in measured
            row['cycling_proxy_is_single_Ki67_not_full_original_multigene_rule'] = a=='R_cycling' or b=='E_cycling'
            row['complete_formula_measurement_refers_to_protein_proxy_only'] = True
        rows.append(row)
    return {'source':panel['source'],'panel_key':panel['panel_key'],'assay':assay,
            'canonical_measured_markers':sorted(measured),'pair_opportunities':rows,
            'context':panel.get('context',{}),'proof_SHA256':panel['proof_SHA256']}


def main():
    p=argparse.ArgumentParser();p.add_argument('--catalog',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    assert not args.output.exists()
    panels=json.loads(args.catalog.read_text())
    assert len({(x['source'],x['panel_key']) for x in panels})==len(panels)
    results=[qualify(x) for x in panels]
    args.output.write_text(json.dumps(results,indent=2)+'\n')
    print('COMPLETE_SOURCE_PANEL_OBSERVABILITY_ONLY_NO_NEW_EFFECTS')


if __name__=='__main__':main()

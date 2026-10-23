"""Read accepted metadata only; emit source/context coverage, never unit rows.

No expression, scores, gene lists, clinical outcomes, files or network writes.
This audits a frozen older capture definition; it does not confer cell identity.
"""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import io
import itertools
import json
from pathlib import Path
import sys

VIEWS = ['R_FOXP3_capture', 'E_effector_capture', 'Myeloid']
PROGRAMS = ['MCP1', 'MCP2']
THRESHOLDS = [20, 50, 100, 200, 500, 1000]
CONTEXT = ['dataset', 'cancer_label', 'tissue', 'material_class', 'phase_class', 'capture_scope', 'ex_vivo']
BASE = 'immune_RNA_MCP_all_atlas_opportunities_20261019_v1'
NATIVE = 'immune_RNA_MCP_all_atlas_native_measurement_20261019_v2'
TERMINALS = {
    BASE + '/private_terminal.json': '2a37f7af2f73ac5ceddcb065626df411ca6f016271fce7f0830a1bd658f83cd1',
    'immune_RNA_MCP_all_atlas_opportunities_audit_20261019_v1/private_terminal.json': '7f181b10ce21de7b84b0ebb0f34ca88a34946afcaabba38a0c48db2b24157e0e',
    NATIVE + '/private_terminal.json': '2463d66299d7b1f1be96d99a19554a834735e688440716edfa00cd1aa01575e9',
    'immune_RNA_MCP_all_atlas_native_measurement_audit_20261019_v1/private_terminal.json': 'f7a6db048cfccac2a98de84023ce5048d4faa595353d2589767718d0e661cfc0',
}


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def boolean(value):
    require(str(value).lower() in ['true', 'false'], 'boolean_not_explicit')
    return str(value).lower() == 'true'


def load_bound(root, name, expected, evidence):
    path = root / name
    require(path.stat().st_size <= 1024**2, 'metadata_size_limit')
    raw = path.read_bytes()
    got = hashlib.sha256(raw).hexdigest()
    require(got == expected, 'input_SHA_mismatch')
    evidence[name] = dict(SHA256=got, bytes=len(raw))
    return raw


def parse_tsv(raw):
    return list(csv.DictReader(io.StringIO(raw.decode()), delimiter='\t'))


def validate_materials(rows):
    require(len({r['sampleID'] for r in rows}) == len(rows), 'duplicate_material')
    for r in rows:
        for k in ['reviewed_patient', 'capture_gene_panels_complete', 'reviewed_discovery_case_overlap']:
            r[k] = boolean(r[k])
        for k in ['core_T', *VIEWS, 'annotated_scope_cells']:
            r[k] = int(r[k])
            require(r[k] >= 0, 'negative_cells')
        require(r['R_FOXP3_capture']+r['E_effector_capture'] <= r['core_T'], 'capture_exceeds_core')
        require(sum(r[v] for v in VIEWS) <= r['annotated_scope_cells'], 'capture_exceeds_material')
        require(all(k in r for k in CONTEXT), 'missing_context')
        require(bool(r['case_block']) or not r['reviewed_patient'], 'reviewed_case_missing')
    return rows


def count_units(rows):
    reviewed = [r for r in rows if r['reviewed_patient']]
    units = {(r['dataset'], r['case_block']) for r in reviewed}
    return dict(materials=len(rows), reviewed_case_blocks=len(units),
                reviewed_materials=len(reviewed),
                extra_reviewed_materials=len(reviewed)-len(units),
                unreviewed_library_units=sum(not r['reviewed_patient'] for r in rows),
                reviewed_cases_outside_original_MCP_discovery=len({
                    (r['dataset'], r['case_block']) for r in reviewed
                    if not r['reviewed_discovery_case_overlap']}))


def aggregate(rows, axes):
    groups = defaultdict(list)
    for r in rows:
        groups[tuple(r[k] for k in CONTEXT)].append(r)
    out = []
    for key, group in sorted(groups.items()):
        source = key[0]
        for left, right in itertools.combinations(VIEWS, 2):
            usable = [r for r in group if r['capture_gene_panels_complete'] and r['core_T'] >= 100]
            counts = []
            for n in THRESHOLDS:
                subset = [r for r in usable if min(r[left], r[right]) >= n]
                counts.append(dict(min_cells_each=n, **count_units(subset)))
            flags = {p: bool(axes[(source, left, p)] and axes[(source, right, p)]) for p in PROGRAMS}
            out.append(dict(**dict(zip(CONTEXT, key)), left=left, right=right,
                            available=count_units(group), complete_program_axes=flags,
                            thresholds=counts))
    return out


def check_artificial():
    def row(sample, case, reviewed, n, phase='pre'):
        return dict(dataset='ARTIFICIAL', cancer_label='toy', tissue='toy', material_class='toy',
                    phase_class=phase, capture_scope='toy', ex_vivo='false', sampleID=sample,
                    case_block=case, reviewed_patient=str(reviewed), capture_gene_panels_complete='true',
                    reviewed_discovery_case_overlap='false', core_T=str(2*n), R_FOXP3_capture=str(n),
                    E_effector_capture=str(n), Myeloid=str(n), annotated_scope_cells=str(3*n))
    source = [row('a','caseA',True,500), row('b','caseA',True,100),
              row('c','unknownC',False,500), row('d','caseA',True,100,'post')]
    rows = validate_materials([r.copy() for r in source])
    axes = {('ARTIFICIAL',v,p): p == 'MCP1' for v in VIEWS for p in PROGRAMS}
    summary = aggregate(rows, axes)
    pre = next(r for r in summary if r['phase_class']=='pre' and r['left']==VIEWS[0] and r['right']==VIEWS[1])
    t = next(r for r in pre['thresholds'] if r['min_cells_each']==100)
    require((t['materials'],t['reviewed_case_blocks'],t['extra_reviewed_materials'],t['unreviewed_library_units'])==(3,1,1,1), 'artificial_units_wrong')
    require(pre['complete_program_axes']=={'MCP1':True,'MCP2':False}, 'artificial_missing_axis')
    require(not any(x['materials'] for x in pre['thresholds'] if x['min_cells_each']==1000), 'empty_scope_wrong')
    require(len(summary)==6, 'phase_merge')
    for field,value in [('reviewed_patient','unknown'),('R_FOXP3_capture','-1')]:
        invalid=[r.copy() for r in source];invalid[0][field]=value
        try: validate_materials(invalid)
        except ValueError: pass
        else: raise AssertionError('invalid_metadata_accepted')
    try: validate_materials([source[0].copy(),source[0].copy()])
    except ValueError: pass
    else: raise AssertionError('duplicate_accepted')
    return dict(status='PASS_ARTIFICIAL_REPEATS_PHASES_UNKNOWN_PATIENT_MISSING_AXES_INVALID_METADATA')


def main(root):
    evidence = {}
    terms = {name:json.loads(load_bound(root,name,expected,evidence)) for name,expected in TERMINALS.items()}
    opp=terms[BASE+'/private_terminal.json']; native=terms[NATIVE+'/private_terminal.json']
    require(terms['immune_RNA_MCP_all_atlas_opportunities_audit_20261019_v1/private_terminal.json']['producer_terminal_sha256']==TERMINALS[BASE+'/private_terminal.json'], 'audit_producer_mismatch')
    def table(folder, term, name):
        return parse_tsv(load_bound(root,folder+'/'+name,term['files'][name]['sha256'],evidence))
    rows=validate_materials(table(BASE,opp,'private_all400_material_RNA_capture_scope.tsv'))
    axis=table(BASE,opp,'private_all54_source_axis_qualification.tsv')
    require(len(rows)==400 and len(axis)==54 and len({r['dataset'] for r in rows})==9,'frozen_scope_changed')
    axes={(r['dataset'],r['view'],r['program']):boolean(r['complete_axis']) for r in axis}
    require(len(axes)==54,'duplicate_axis')
    require(set(axes)=={(s,v,p) for s in {r['dataset'] for r in rows} for v in VIEWS for p in PROGRAMS},'incomplete_axis_grid')
    pools=table(NATIVE,native,'private_all_capture_pools.tsv')
    require(len(pools)==877,'pool_scope_changed')
    byid={r['sampleID']:r for r in rows}
    for p in pools:
        r=byid[p['sampleID']]
        require(int(p['full_view_cells'])==r[p['view']], 'native_full_count_mismatch')
        require(int(p['selected_cells'])==min(100,int(p['full_view_cells'])), 'selected_cap_changed')
    require(len({(p['sampleID'],p['view']) for p in pools})==877,'duplicate_native_pool')
    summary=aggregate(rows,axes)
    # Independent scalar reconstruction through counters, not subset extraction.
    for r in summary:
        for t in r['thresholds']:
            case_counts=Counter(); materials=0;unknown=0
            for raw in rows:
                same=all(raw[k]==r[k] for k in CONTEXT)
                enough=raw[r['left']]>=t['min_cells_each'] and raw[r['right']]>=t['min_cells_each']
                if same and enough and raw['capture_gene_panels_complete'] and raw['core_T']>=100:
                    materials+=1
                    if raw['reviewed_patient']:case_counts[(raw['dataset'],raw['case_block'])]+=1
                    else:unknown+=1
            require((materials,len(case_counts),unknown)==(t['materials'],t['reviewed_case_blocks'],t['unreviewed_library_units']),'independent_counts_disagree')
    for name, info in evidence.items():
        require(hashlib.sha256((root/name).read_bytes()).hexdigest()==info['SHA256'],'input_changed_during_read')
    return dict(status='COMPLETE_METADATA_ONLY_DISTRIBUTION_COVERAGE', observed_UTC=datetime.now(timezone.utc).isoformat(),
                code_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest() if '__file__' in globals() and Path(__file__).is_file() else None,
                evidence=evidence, input_materials=400, sources=9, native_capture_pools=877,
                whole_dataset_units=count_units(rows),
                native_selected_score_cap=100, full_score_distribution_available_from_saved_sums=False,
                context_pair_rows=len(summary), thresholds=THRESHOLDS, results=summary,
                no_expression_or_score_or_clinical_outcome_read=True, no_individual_rows_saved_or_returned=True,
                coverage_not_precision_or_identity_acceptance=True,checks=check_artificial())


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path);parser.add_argument('--check-artificial',action='store_true');args=parser.parse_args()
    try:
        result=check_artificial() if args.check_artificial else main(args.root)
        text=json.dumps(result,allow_nan=False)
        require(len(text.encode())<=512*1024,'aggregate_output_limit')
        print(text)
    except Exception as error:
        print(json.dumps(dict(status='COVERAGE_AUDIT_FAILED',error_class=type(error).__name__,individual_values_not_returned=True)))
        sys.exit(2)

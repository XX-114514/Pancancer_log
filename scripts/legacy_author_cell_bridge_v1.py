"""Exact author-metadata join for a declared historical cell-ID encoding.

This preserves the author sample prefix and complete barcode. It cannot
certify sample design, patient independence, QC, or an analysis view.
"""
from collections import Counter, defaultdict


def encoded_author_cell(source, namespace, cell_id):
    """Reproduce the explicitly reviewed colon-to-underscore encoding only."""
    if not all(isinstance(x,str) and x for x in [source,namespace,cell_id]):
        return None
    if any('|' in x for x in [source,namespace,cell_id]):
        raise ValueError('Reserved historical encoding separator')
    return source+'__'+namespace+'__'+cell_id.replace(':','_')


def decode_legacy_global_cell(source, value):
    prefix=source+'|project_v3_sample_id='
    if not isinstance(value,str) or not value.startswith(prefix) or value.count('|cell=')!=1:
        return None
    unit,cell=value[len(prefix):].split('|cell=')
    if not unit.startswith(source+'__') or not cell.startswith(unit+'__'):
        return None
    return cell


def exact_author_bridge(source, global_ids, author_cell_ids, author_namespaces):
    if len(author_cell_ids)!=len(author_namespaces):
        raise ValueError('Author metadata axes differ')
    lookup=defaultdict(list)
    for i,(cell,namespace) in enumerate(zip(author_cell_ids,author_namespaces)):
        key=encoded_author_cell(source,namespace,cell)
        if key is not None:lookup[key].append(i)
    decoded=[decode_legacy_global_cell(source,x) for x in global_ids]
    positions=[];status=[]
    for key in decoded:
        hits=lookup.get(key,[]) if key is not None else []
        if key is None:positions.append(-3);status.append('invalid_declared_namespace')
        elif not hits:positions.append(-1);status.append('unmatched')
        elif len(hits)>1:positions.append(-2);status.append('ambiguous_author_identity')
        else:positions.append(hits[0]);status.append('exact_declared_encoding')
    repeated_author_rows={k:v for k,v in Counter(i for i in positions if i>=0).items() if v>1}
    return {'author_row_position':positions,'status':status,
            'status_counts':dict(Counter(status)),
            'duplicate_author_encoded_keys':sum(len(v)>1 for v in lookup.values()),
            'repeated_author_rows_in_archive':repeated_author_rows,
            'complete_one_to_one':bool(positions) and all(i>=0 for i in positions) and not repeated_author_rows,
            'bare_barcode_matching_used':False,'source_design_or_QC_qualified':False}

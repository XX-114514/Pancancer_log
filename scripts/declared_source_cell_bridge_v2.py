"""Closed forward ID encodings with complete author-cell and library keys.

No barcode suffix is removed; no expression or phenotype chooses a match.
These mappings do not establish source design, clinical or QC eligibility.
"""
from collections import defaultdict,Counter

RECIPES=('new_cohort_source_prefixed_full_index','legacy_wrapped_h5_full_index')

def encoded_identity(source,raw_index,library,recipe):
    if recipe not in RECIPES:raise ValueError('Unreviewed source ID recipe')
    if not all(isinstance(v,str) and v for v in [source,raw_index,library]):return None
    if any('|' in v for v in [source,raw_index,library]):raise ValueError('Reserved ID delimiter')
    if recipe==RECIPES[0]:return source+'__'+raw_index,source+'|'+library
    unit=source+'__'+library
    return source+'|project_v3_sample_id='+unit+'|cell='+unit+'__'+raw_index,source+'|'+unit


def recover_rows(source,archive_global_ids,archive_units,raw_indices,raw_libraries,recipe):
    if len(archive_global_ids)!=len(archive_units) or len(raw_indices)!=len(raw_libraries):raise ValueError('ID and library axes differ')
    index=defaultdict(list)
    for i,(raw,lib) in enumerate(zip(raw_indices,raw_libraries)):
        key=encoded_identity(source,raw,lib,recipe)
        if key is not None:index[key].append(i)
    positions=[];status=[]
    for cell,lib in zip(archive_global_ids,archive_units):
        hits=index.get((cell,lib),[])
        if not hits:positions.append(-1);status.append('unmatched')
        elif len(hits)>1:positions.append(-2);status.append('ambiguous_author_key')
        else:positions.append(hits[0]);status.append('exact_complete_cell_and_library')
    repeated={k:v for k,v in Counter(i for i in positions if i>=0).items() if v>1}
    return {'raw_row_position':positions,'status_counts':dict(Counter(status)),
            'duplicate_author_composite_keys':sum(len(v)>1 for v in index.values()),
            'repeated_author_rows_in_archive':repeated,
            'complete_one_to_one':bool(positions) and all(x>=0 for x in positions) and not repeated,
            'barcode_suffix_removed':False,'biological_or_clinical_eligibility_assigned':False}

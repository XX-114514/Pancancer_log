"""Current native gene-axis qualification only; no cell rows or matrix reads.

Source entries arrive in RAM. Only source/cancer aggregates leave the reader.
Gene presence is not a qualification of capture, QC, geometry or clinical units.
"""
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

PARENT = ['CD3D', 'CD3E', 'TRAC', 'FOXP3', 'IL2RA', 'CD8A']
DEFINITIONS = {
    'R_TNFRSF9': PARENT + ['TNFRSF9'],
    'R_LAG3': PARENT + ['LAG3'],
    'R_cycling': PARENT + ['MKI67', 'TOP2A', 'UBE2C'],
    'E_dysfunction': PARENT + ['PDCD1', 'TOX', 'HAVCR2'],
    'E_cycling': PARENT + ['MKI67', 'TOP2A', 'UBE2C'],
    'E_memory': PARENT + ['TCF7', 'IL7R'],
    'alternative_R_parent': PARENT + ['CTLA4'],
}


def strings(values):
    return [v.decode('utf-8') if isinstance(v, bytes) else str(v) for v in values]


def qualify(names):
    names = [s.strip().upper() for s in names]
    counts = Counter(names)
    required = set(sum(DEFINITIONS.values(), []))
    return dict(
        complete={k: all(counts[g] == 1 for g in gs) for k, gs in DEFINITIONS.items()},
        missing=sorted(g for g in required if counts[g] == 0),
        ambiguous=sorted(g for g in required if counts[g] > 1),
        axis_SHA256=hashlib.sha256(json.dumps(names, separators=(',', ':')).encode()).hexdigest(),
    )


def read_axis(group, name, h5py):
    obj = group[name]
    if isinstance(obj, h5py.Dataset):
        return strings(obj[:])
    import numpy as np
    assert str(obj.attrs.get('encoding-type', '')) in ['categorical', "b'categorical'"]
    categories = strings(obj['categories'][:])
    codes = np.asarray(obj['codes'], dtype=int)
    assert (codes >= -1).all() and (codes < len(categories)).all()
    return [categories[int(i)] if i >= 0 else '' for i in codes]


def run(rows):
    import h5py
    groups = defaultdict(list)
    for row in rows:
        result = {}
        try:
            path = Path(row['path'])
            before = path.stat()
            with h5py.File(path, 'r') as h:
                var = h['var']
                index = var.attrs.get('_index', '_index')
                if isinstance(index, bytes):
                    index = index.decode()
                selected = 'gene_symbol' if 'gene_symbol' in var else str(index)
                names = read_axis(var, selected, h5py)
                result = qualify(names)
                assert len(names) == int(row['n_native_genes'])
                result['selected_gene_column'] = selected
            after = path.stat()
            assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
            result['status'] = 'COMPLETE_CURRENT_GENE_AXIS_ONLY'
        except (OSError, KeyError, ValueError, AssertionError) as error:
            result = dict(status='UNQUALIFIED_AXIS', error_type=type(error).__name__)
        groups[(row['dataset_id'], row['cancer'])].append(result)
    output = []
    for (source, cancer), records in sorted(groups.items()):
        good = [r for r in records if r['status'] == 'COMPLETE_CURRENT_GENE_AXIS_ONLY']
        output.append(dict(source=source, cancer=cancer, input_files=len(records),
            accepted_current_axis_files=len(good),
            failed_axis_reasons=dict(Counter(r.get('error_type') for r in records if r not in good)),
            number_complete_by_definition={k: sum(r['complete'][k] for r in good) for k in DEFINITIONS},
            absent_required_genes_union=sorted({g for r in good for g in r['missing']}),
            ambiguous_required_genes_union=sorted({g for r in good for g in r['ambiguous']}),
            distinct_native_gene_axes=len({r['axis_SHA256'] for r in good}),
            selected_gene_columns=sorted({r['selected_gene_column'] for r in good}),
            molecular_capture_QC_geometry_patient_phase_endpoint_not_qualified=True))
    return dict(status='COMPLETE_ALL_DECLARED_GENE_AXIS_SCOPE_NO_NEW_STATE_EFFECT',
                source_cancer_aggregates=output, no_X_obs_or_obsm_read=True,
                input_scope_SHA256=hashlib.sha256(json.dumps(rows, sort_keys=True).encode()).hexdigest())

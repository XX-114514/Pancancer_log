"""All cached source annotation fields, metadata only; no reference promotion."""
import hashlib
import json
import os
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem

from gamma_author_text_semantics_v2 import KINDS, label_kind, fixture


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def candidate_field(name):
    if re.search(r'barcode|sample|patient|donor|subject|library|response|treatment|cancer|tumou?r|tissue|time|batch|clinical', name, re.I):
        return False
    return bool(re.search(r'cell.?type|annotation|cell.?annot|sub.?type|(?:^|_)(?:cluster|leiden|louvain|cellclass|lineage|label|phenotype)(?:$|_)', name, re.I))


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    assert shutil.disk_usage('/dev/shm').free > 512 * 1024**2
    contract = root.parents[2] / 'Project_v5/runs/20261013_immunotherapy_atlas_v1/config/source_contract_v2_private.json'
    assert sha(contract) == '8d7f54a385389d0ecbc825155df727ed065cae2ed1c3f3ed584c1a7896647e71'
    prior = root / 'gamma_group_lineage_orthogonal_description_20261019_v1'
    audit = root / 'gamma_group_lineage_orthogonal_description_audit_20261019_v1/private_terminal.json'
    assert sha(audit) == '32d0dba881c9029c61839ccef8b449a27ddb1142453c3d889551fcff3a1e3ec9'
    proof = json.loads(audit.read_text())
    assert proof['producer_terminal_sha256'] == sha(prior/'private_terminal.json')
    t = json.loads((prior/'private_terminal.json').read_text())
    schema = prior/'private_all21_source_orthogonal_schema_inventory.json'
    assert sha(schema) == t['files'][schema.name]['sha256']
    old = json.loads(schema.read_text())
    records = json.loads(contract.read_text())['records']
    assert len(records) == len(old) == 21
    out.mkdir(mode=0o700)
    save = lambda name, value: (out/name).write_text(json.dumps(value, indent=2)+'\n')
    save('plan_before_cached_fine_author_values.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(),
         source_sha256=sha(__file__), semantic_source_sha256=sha(Path(__file__).with_name('gamma_author_text_semantics_v2.py')),
         source_contract_sha256=sha(contract), prior_schema_audit_sha256=sha(audit),
         record_positions=21, all_annotation_fields_preserved=True, max_fields_per_record=64,
         maximum_distinct_values_per_field=20000, labels_are_not_gold_standard=True,
         no_RNA_or_clinical_effects=True, storage='Task-owned tmpfs then private local small-copy; project read-only', output_cap_bytes=2*1024**2))
    save('private_semantic_fixture.json', fixture())
    summaries, terms, proofs = [], [], []
    for index, record in enumerate(records):
        previous = [v for v in old if v['source_path']==record['path'] and v['source_role']==record['role']]
        assert len(previous)==1
        previous = previous[0]
        selector = Path(record['row_positions_ref'])
        assert sha(selector)==previous['row_selector_sha256']
        rows = np.load(selector, allow_pickle=False)
        assert len(rows)==record['n_cells'] and len(rows)==len(np.unique(rows))
        path=Path(record['path'])
        before=(path.stat().st_size,path.stat().st_mtime_ns)
        h=h5py.File(path,'r') if path.suffix=='.h5ad' else zarr.open_group(str(path),mode='r')
        try:
            obs=h['obs'];fields=sorted(obs.keys())
            assert fields==previous['obs_fields']
            selected=[v for v in fields if candidate_field(v)]
            assert len(selected)<=64
            if not selected:
                summaries.append(dict(record_index=index,dataset=record['dataset'],role=record['role'],field='',status='NO_MATCHING_ANNOTATION_FIELD',selected_rows=len(rows)))
            for name in selected:
                raw=np.asarray(read_elem(obs[name]))
                assert len(raw)==previous['raw_metadata_rows']
                values=pd.Series(raw[rows]).fillna('AUTHOR_LABEL_MISSING').astype(str)
                counts=values.value_counts(dropna=False,sort=False)
                row=dict(record_index=index,dataset=record['dataset'],role=record['role'],field=name,selected_rows=len(rows),distinct_values=len(counts),values_sha256=hashlib.sha256(json.dumps(values.tolist(),ensure_ascii=False,separators=(',',':')).encode()).hexdigest())
                if len(counts)>20000:
                    row.update(status='HIGH_CARDINALITY_FIELD_NOT_TREATED_AS_IDENTITY_REFERENCE',explicit_gamma_cells=np.nan)
                else:
                    # Alternate unique count checks all selected metadata terms.
                    names,numbers=np.unique(values.to_numpy(str),return_counts=True)
                    assert dict(zip(names,numbers))==counts.to_dict()
                    numbers_by_kind=np.zeros(len(KINDS),np.int64)
                    for label,number in counts.items():
                        kind=label_kind(label);numbers_by_kind[kind]+=int(number)
                        terms.append(dict(record_index=index,dataset=record['dataset'],role=record['role'],field=name,author_label=label,semantic_kind=KINDS[kind],n_cells=int(number)))
                    assert int(numbers_by_kind.sum())==len(rows)
                    row.update(status='ANNOTATION_VALUES_COUNTED_NOT_REFERENCE_ACCEPTANCE',explicit_gamma_cells=int(numbers_by_kind[0]),delta_subtype_text_cells=int(numbers_by_kind[1]),gamma_marker_text_cells=int(numbers_by_kind[2]),mixed_T_NK_text_cells=int(numbers_by_kind[8]),independent_TCR_or_patient_validation=False)
                summaries.append(row)
        finally:
            if isinstance(h,h5py.File):h.close()
        assert before==(path.stat().st_size,path.stat().st_mtime_ns)
        proofs.append(dict(record_index=index,dataset=record['dataset'],role=record['role'],row_selector_sha256=sha(selector),selected_rows=len(rows),all_obs_schema_equal_prior=True,source_stat_unchanged=True,field_count=len(selected)))
        save('private_progress.json',dict(stage='CACHED_FINE_ANNOTATION_VALUES',completed_records=index+1,total_records=21,observed_UTC=datetime.now(timezone.utc).isoformat()))
    pd.DataFrame(summaries).to_csv(out/'private_complete21_annotation_field_qualification.tsv',sep='\t',index=False)
    pd.DataFrame(terms).to_csv(out/'private_complete_source_field_author_label_counts.tsv.gz',sep='\t',index=False,compression='gzip')
    save('private_all21_canonical_selector_and_source_proofs.json',proofs)
    receipt=dict(status='COMPLETE_ALL21_EXISTING_CACHE_FINE_ANNOTATION_METADATA_NO_RNA_OR_REFERENCE_PROMOTION',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),record_positions=21,source_role_overlap_not_removed=True,metadata_count_check='Pandas value_counts versus NumPy unique on every selected field',full_expression_or_cell_scores_persisted=False,original_identity_changed=False,scientific_outputs_public=False,files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(v['bytes'] for v in receipt['files'].values())<2*1024**2
    save('private_terminal.json',receipt)
    print(json.dumps(dict(status=receipt['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':
    main()

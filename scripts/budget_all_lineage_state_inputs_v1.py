from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,re,subprocess
import numpy as np,pandas as pd,zarr
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT'])
D=R/'all_lineage_state_discovery_row_selection_20261014_v1'
assert json.loads((R/'all_lineage_state_discovery_row_audit_20261014_v1/independent_audit.json').read_text())['status'].startswith('PASS_ALL100766')
z=np.load(D/'private_discovery_row_indices.npz');rows=z['canonical_global_rows'];m=pd.read_csv(D/'private_discovery_metadata.tsv.gz',sep='\t',keep_default_na=False)
contract=json.loads((AR/'20261010_uniform_expression_annotation_v1/input_contract.json').read_text());g=zarr.open_group(contract['core_zarr'],mode='r');ptr=g['X/indptr'][:]
nnz=ptr[rows+1]-ptr[rows];assert len(nnz)==100766
quota=subprocess.check_output(['/usr/bin/quota','-g'],text=True)
hit=re.search(r'^\s*(\d+)\s+0\s+(\d+)\s+\d+\s+0\s+0\s*$',quota,re.M);assert hit
free=(int(hit[2])-int(hit[1]))*1024
# Counts are one sparse common-axis copy per major, int32 values/indices.
count_bound=int(nnz.sum())*8+(len(rows)+18)*8
budget=int(np.ceil((count_bound*1.5+2**30)/2**30))
out=dict(status='ROW_POINTER_STORAGE_BUDGET_NO_EXPRESSION_COPIED',observed_UTC=datetime.now(timezone.utc).isoformat(),selected_cells=len(rows),native_nonzero_upper_bound=int(nnz.sum()),raw_CSR_byte_upper_bound=count_bound,requested_stage_GiB=budget,reserve_GiB=500,free_GiB=free/2**30,capacity_gate_pass=free>=(500+budget)*2**30,scope='One common-axis int32 count matrix per original major, native total/detected metadata only; no full4.53m copy; models budget separately',source_selection_inventory_sha256=hashlib.sha256((D/'inventory.json').read_bytes()).hexdigest(),major_bounds=[dict(major=v,cells=int(mask.sum()),native_nnz_bound=int(nnz[mask].sum())) for v in sorted(m.major.unique()) for mask in [m.major.eq(v).to_numpy()]])
O=R/'all_lineage_state_input_budget_20261014_v1';assert not O.exists();O.mkdir();(O/'budget.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))

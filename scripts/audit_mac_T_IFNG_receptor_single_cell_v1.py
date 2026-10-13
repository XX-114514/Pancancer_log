"""All fine-unit/marginal aggregation audit and deterministic fresh row examples."""
from pathlib import Path
import os,json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd,zarr
from scipy.sparse import load_npz
from anndata._io.specs import read_elem

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);P=R/'immune_same_name_family_source_review_v1';O=P/'mac_T_IFNG_receptor_single_cell_v1';U=AR/'20261010_uniform_expression_annotation_v1';plan=json.loads(Path(os.environ['PANCANCER_SINGLE_CELL_PLAN']).read_text());assert sha(__file__)==plan['auditor_sha256'];inv=json.loads((O/'inventory.json').read_text())
 for file,h in inv['outputs'].items():assert sha(O/file)==h
 E=pd.read_csv(O/'private_all1832_sample_node_RNA_counts.tsv',sep='\t');M=pd.read_csv(O/'private_all424_sample_major_RNA_counts.tsv',sep='\t');a=np.load(O/'private_selected_cell_gene_bits.npz');assert len(a['global_rows'])==plan['expected_cells'] and len(E)==plan['expected_units'];unit=a['unit_index'];b=a['gene_bits'];assert np.all(np.diff(a['global_rows'])>0) and (b<=7).all();assert np.array_equal(E.n_cells,np.bincount(unit,minlength=len(E)))
 for j,gene in enumerate(plan['genes']):
  expected=np.bincount(unit[b&(1<<j)>0],minlength=len(E));assert np.array_equal(E['positive_'+gene],expected)
 expected=np.bincount(unit[(b&6)==6],minlength=len(E));assert np.array_equal(E.positive_IFNGR1_IFNGR2,expected);assert (E.positive_IFNGR1_IFNGR2<=E[['positive_IFNGR1','positive_IFNGR2']].min(axis=1)).all()
 for sid,sub in E.groupby('gse_id'):
  units=pd.read_csv(U/'by_gse'/sid/'units.tsv',sep='\t');original=pd.read_csv(U/'all_node_number_catalog.tsv',sep='\t');units['node_id']=original.node_id.to_numpy()[units.node_index];units=units.set_index(['sample_uid','node_id']);genes=pd.read_csv(U/'by_gse'/sid/'genes.tsv',sep='\t');positive=load_npz(U/'by_gse'/sid/'positive_counts.npz');rows=units.loc[pd.MultiIndex.from_frame(sub[['sample_uid','node_id']]),'unit_row'].to_numpy()
  for gene in plan['genes']:col=int(np.flatnonzero(genes.gene==gene)[0]);assert np.array_equal(np.asarray(positive[rows,col].todense()).ravel(),sub['positive_'+gene])
 cols=['n_cells','positive_IFNG','positive_IFNGR1','positive_IFNGR2','positive_IFNGR1_IFNGR2'];rebuilt=E.groupby(['gse_id','sample_uid','major'],as_index=False)[cols].sum().sort_values(['gse_id','sample_uid','major']).reset_index(drop=True);pd.testing.assert_frame_equal(rebuilt,M.sort_values(['gse_id','sample_uid','major']).reset_index(drop=True),check_dtype=False)
 contract=json.loads((U/'input_contract.json').read_text());g=zarr.open_group(contract['core_zarr'],mode='r');X=g['X'];columns=json.loads(str(a['source_column_maps_json']));examples=[]
 for (sid,major),sub in E.groupby(['gse_id','major']):
  positions=np.flatnonzero(np.isin(unit,sub.unit_index));chosen=[]
  for code in [0,1,2,4,6,7]:
   eligible=positions[b[positions]==code]
   if len(eligible):chosen.extend(eligible[np.linspace(0,len(eligible)-1,min(8,len(eligible)),dtype=int)].tolist())
  for pos in sorted(set(chosen)):
   row=int(a['global_rows'][pos]);ip=X['indptr'][row:row+2].astype(int);ix=X['indices'][int(ip[0]):int(ip[1])];values=X['data'][int(ip[0]):int(ip[1])];code=0
   for j,gene in enumerate(plan['genes']):
    if (values[ix==columns[sid][gene]]>0).any():code|=1<<j
   assert code==int(b[pos]);examples.append(dict(gse_id=sid,major=major,bit_pattern=code))
 counts=pd.DataFrame(examples).groupby(['gse_id','major','bit_pattern']).size().reset_index(name='examples');counts.to_csv(O/'independent_raw_example_counts.tsv',sep='\t',index=False);j=dict(status='PASS_ALL1832_UNITS_ALL3_MARGINALS_AND_FRESH_ROW_EXAMPLES',finished_UTC=datetime.now(timezone.utc).isoformat(),inventory_sha256=sha(O/'inventory.json'),cells=plan['expected_cells'],all_unit_gene_checks=3*len(E),same_cell_bit_aggregation_all=True,fresh_raw_row_examples=len(examples),all_rows_independently_replayed=False,not_protein_assembly_or_independent_communication_validation=True,auditor_sha256=sha(__file__));(O/'independent_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j),flush=True)
if __name__=='__main__':main()

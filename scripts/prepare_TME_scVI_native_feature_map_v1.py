"""Freeze native column maps for scVI without copying any expression matrix."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main(P,out):
 assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 cfg=load(P/'config/spec_v1.json');J=Path(cfg['joined']);sys.path.insert(0,str(J/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection
 from native_annotation_review_v1 import native_measured_genes
 native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest'])
 selected=load(P/'outputs/selected_HVG_names_v1.json');scope=load(P/'outputs/scope_and_common_genes_v1.json');lineage=load(P/'outputs/source_input_lineage_v1.json');records=[]
 assert [r['source'] for r in lineage]==scope['sources']
 for r in lineage:
  names,cols=native_measured_genes(native,r['source']);by=dict(zip(names,map(int,cols)));columns=[by[g] for g in selected]
  assert len(columns)==len(set(columns))==cfg['n_HVGs']
  records.append(dict(r,selected_columns=columns));native._fields.clear()
 result=dict(producer=str(P),producer_manifest_SHA256=sha(P/'outputs/standard_TME_Harmony_manifest_v1.json'),axes_path=str(P/'outputs/private_cell_axes_v1.npz'),axes_SHA256=sha(P/'outputs/private_cell_axes_v1.npz'),genes=selected,genes_SHA256=sha(P/'outputs/selected_HVG_names_v1.json'),n_fit=scope['n_fit_records'],records=records,expression_matrix_copied=False)
 with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--producer',required=True);p.add_argument('--output',required=True);a=p.parse_args();main(Path(a.producer),Path(a.output))

"""Check plotted axes, all source/cluster denominators and bounded raw markers.

This is visualization provenance/readback, not proof that a UMAP or a cell
identity is biologically correct. PNG visual review remains a separate step.
"""
import argparse,collections,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from PIL import Image

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')

def load_required_display_axes(path):
 required=('UMAP_Harmony','major','fine','source','raw','corrected','genes')
 with np.load(path,allow_pickle=False) as f:return {k:f[k] for k in required}

def main(root):
 cfg=load(root/'config/spec_v1.json');R=Path(cfg['producer']);EA=Path(cfg['independent_native_auditor'])
 while not (R/'provenance/supervisor_terminal_v1.json').exists():time.sleep(3)
 assert load(R/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0;rc=load(R/'config/spec_v1.json');P=Path(rc['producer']);pc=load(P/'config/spec_v1.json');J=Path(rc['joined']);jm=load(J/'outputs/joined_native_atlas_delivery_v1.json');qs=load(Path(jm['projection_qualifier'])/'provenance/terminal_v1.json');assert sha(jm['qualified_projection'])==qs['output_sha256'];proj=load(jm['qualified_projection']);pi={(r['source'],r['model_id']):r['base_projection_record'] for r in proj['records']};baseline=load(Path(pc['working_delivery'])/'outputs/accepted_TME_scoped_working_annotation_v1.json');w=load(baseline['producer_manifest']['path']);wb=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};catalog=load(P/'outputs/cluster_hierarchy_named_catalog_v1.json')['nodes'];nodes={(r['source'],r['source_node_code']):r for r in catalog};sys.path.insert(0,str(R/'scripts'));sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(EA/'scripts'))
 from TME_visualization_sampling_v1 import balanced_indices
 from render_TME_UMAP_marker_review_v1 import MARKERS
 from native_archive_collection_v1 import NativeArchiveCollection
 from audit_full_TME_native_evidence_v1 import gather
 native=NativeArchiveCollection(jm['archive_manifest']);scope=load(R/'outputs/UMAP_display_scope_v1.json');coverage={r['source']:r for r in scope['source_coverage']};manifest=load(R/'outputs/visualization_manifest_v1.json');rt=load(R/'provenance/terminal_v1.json');assert sha(R/'outputs/visualization_manifest_v1.json')==rt['manifest_SHA256']
 # The legacy baseline dictionary was serialized as an object array. It is not
 # used by this numeric/current-label audit; never deserialize unused pickle
 # payloads merely because they share an NPZ with the required primitive axes.
 display=load_required_display_axes(R/'private_display_arrays_v1.npz')
 assert np.isfinite(display['UMAP_Harmony']).all() and display['UMAP_Harmony'].shape==(len(display['major']),2);n_fit=0;n_all=0;n_display=0;n_covered=0;n_checks=0
 for entry in load(P/'outputs/cluster_hierarchy_index_v1.json'):
  source=entry['source'];m=load(P/'outputs'/entry['metadata']);wr=wi[source]
  with np.load(m['arrays']['path'],allow_pickle=False) as f:pool=f['source_scope_position'];positions=f['source_native_position'];rows=f['native_row'];leaf=f['final_leaf_node']
  with np.load(wb/wr['arrays'],allow_pickle=False) as f:total=f['native_total_RNA'][pool];n_qualified=len(f['native_row']);np.testing.assert_array_equal(pool,np.flatnonzero(f['fit_pool_candidate_not_final_QC']))
  record=pi[source,'global_overview__shared'];p=Path(proj['base_projection'])/'outputs'/record['arrays'];assert sha(p)==record['arrays_sha256']
  with np.load(p,allow_pickle=False) as f:
   dom=f['source_native_position'];idx=np.searchsorted(dom,positions);hit=idx<len(dom);hit[hit]&=dom[idx[hit]]==positions[hit];raw=np.full((len(pool),50),np.nan,np.float32);corrected=raw.copy();ok=np.zeros(len(pool),bool);raw[hit]=f['raw_PCA'][idx[hit]];corrected[hit]=f['candidate_Harmony'][idx[hit]];ok[hit]=f['projection_numerical_supported'][idx[hit]];eligible=hit&ok&np.isfinite(raw).all(axis=1)&np.isfinite(corrected).all(axis=1)
  seed=rc['seed']+int(hashlib.sha256(source.encode()).hexdigest()[:6],16);take=balanced_indices(leaf,eligible,rc['global_per_source'],rc['minimum_per_leaf'],seed);where=np.flatnonzero(display['source']==source);assert len(where)==len(take)==coverage[source]['n_display_sample'];np.testing.assert_array_equal(display['raw'][where],raw[take]);np.testing.assert_array_equal(display['corrected'][where],corrected[take]);np.testing.assert_array_equal(display['major'][where],[nodes[source,int(k)]['cluster_RNA_proposal']['major'] for k in leaf[take]]);np.testing.assert_array_equal(display['fine'][where],[nodes[source,int(k)]['cluster_RNA_proposal']['celltype'] for k in leaf[take]]);assert set(leaf[take])==set(leaf[eligible]);n_covered+=len(set(leaf[take]));n_fit+=len(pool);n_all+=n_qualified;n_display+=len(take)
  check=np.unique(np.linspace(0,len(take)-1,min(len(take),cfg['native_checks_per_source'])).astype(int));selected=take[check];csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']];x=gather(csr,rows[selected]).astype(float).tocsr();np.testing.assert_array_equal(np.asarray(x.sum(axis=1)).ravel(),total[selected]);genes={g:j for j,g in enumerate(m['native_measured_genes'])}
  with np.load(m['statistics']['path'],allow_pickle=False) as f:cols=f['native_columns']
  for j,g in enumerate(MARKERS):
   actual=display['genes'][where[check],j]
   if g not in genes:assert np.isnan(actual).all()
   else:expected=np.log1p(np.asarray(x[:,[int(cols[genes[g]])]].toarray()).ravel()*10000/total[selected]);np.testing.assert_allclose(actual,expected,rtol=2e-6,atol=1e-6)
  n_checks+=len(selected);native._fields.clear()
 assert n_display==scope['n_display_records']==rt['n_display_records'];assert n_covered==scope['n_leaf_groups_represented']==rt['n_leaf_groups_represented'];assert n_fit==scope['full_fit_records'] and n_all==scope['full_qualified_records'];figs=[]
 for f in manifest['files']:
  p=R/'outputs'/f['path'];assert sha(p)==f['SHA256'] and p.stat().st_size==f['bytes']
  if p.suffix=='.png':
   with Image.open(p) as im:im.verify()
   with Image.open(p) as im:rgb=np.asarray(im.convert('RGB'));assert rgb.shape[0]>=400 and rgb.shape[1]>=400 and float(rgb.std())>5;figs.append(dict(name=f['path'],width=rgb.shape[1],height=rgb.shape[0]))
 assert len(figs)==manifest['n_figures']==rt['n_figures'];pdf=(R/'outputs/TME_UMAP_marker_review.pdf').read_bytes();assert pdf.startswith(b'%PDF') and b'%%EOF' in pdf[-100:];save(root/'outputs/visualization_readback_acceptance_v1.json',dict(n_sources=len(coverage),n_fit_records=n_fit,n_qualified_records=n_all,n_display_records=n_display,n_leaf_groups_represented=n_covered,n_native_marker_checks=n_checks,all_display_axes_and_masks_and_source_balanced_memberships_checked=True,marker_check_scope='bounded_per_source_independent_native_CSR_not_all_points',all_PNG_files_decode_and_are_nonblank=True,figures=figs,UMAP_is_visualization_not_identity_or_integration_validation=True));save(root/'provenance/terminal_v1.json',dict(status='UMAP_MARKER_VISUALIZATION_SCOPE_COORDINATES_BOUNDED_NATIVE_MARKERS_AND_IMAGE_READBACK_ACCEPTED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),producer_terminal_SHA256=sha(R/'provenance/terminal_v1.json'),all_no_new_identity_or_CM=True,n_figures=len(figs)))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException:
  import traceback
  save(root/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise

"""Add a cluster tree proposal view to the accepted TME working annotation.

The baseline consensus/author/state fields are preserved. Cells outside the
frozen first-round graph remain readable with an explicit outside-pool flag.
"""
import hashlib,json
from pathlib import Path
import numpy as np
from TME_scoped_annotation_view_v1 import TMEScopedAnnotationView

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

class TMEHierarchicalClusterView:
 def __init__(self,acceptance):
  self.acceptance=load(acceptance);c=self.acceptance;assert sha(c['catalog']['path'])==c['catalog']['SHA256'];assert sha(c['index']['path'])==c['index']['SHA256'];self.catalog=load(c['catalog']['path'])['nodes'];self.index={r['source']:r for r in load(c['index']['path'])};self.root=Path(c['index']['path']).parent;self.baseline=TMEScopedAnnotationView(c['baseline_working_acceptance']['path']);assert sha(c['baseline_working_acceptance']['path'])==c['baseline_working_acceptance']['SHA256'];self.sources=self.baseline.sources;self.cache={}
 def source(self,source):
  if source not in self.cache:
   item=self.index[source];mp=self.root/item['metadata'];assert sha(mp)==item['metadata_SHA256'];m=load(mp);ap=m['arrays']['path'];assert sha(ap)==m['arrays']['SHA256']
   with np.load(ap,allow_pickle=False) as f:a={k:f[k] for k in f.files}
   base=self.baseline.source(source)['arrays'];np.testing.assert_array_equal(a['native_row'],base['native_row'][a['source_scope_position']]);assert (np.diff(a['source_scope_position'])>0).all();self.cache[source]=dict(metadata=m,arrays=a,nodes={r['source_node_code']:r for r in self.catalog if r['source']==source})
  return self.cache[source]
 def metadata(self,source,positions):
  result=self.baseline.metadata(source,positions);p=np.asarray(positions,np.int64);s=self.source(source);a=s['arrays'];pool=a['source_scope_position'];where=np.searchsorted(pool,p);hit=where<len(pool);hit[hit]&=pool[where[hit]]==p[hit];extra={}
  names=['cluster_root_display','cluster_leaf_display','cluster_marker','cluster_major_proposal','cluster_celltype_proposal','cluster_annotation_status','cluster_membership_SHA256','cluster_state_program_names','cluster_view_status','cluster_T_evidence_mode']
  for name in names:
   v=np.empty(len(p),object);v[:]=None;extra[name]=v
  outside='outside_frozen_qualified_fit_pool' if self.acceptance['scope']=='source_qualified_fit_pool_all_lineages' else 'outside_frozen_first_round_pool'
  extra['cluster_view_status'][:]=outside;extra['cluster_leaf_number']=np.full(len(p),-1,np.int64)
  for j in np.flatnonzero(hit):
   pos=where[j];r=s['nodes'][int(a['final_leaf_node'][pos])];parent=s['nodes'][int(a['root_cluster_code'][pos])];call=r['cluster_RNA_proposal'];extra['cluster_root_display'][j]=parent['display_label'];extra['cluster_leaf_display'][j]=r['display_label'];extra['cluster_leaf_number'][j]=r['global_display_number'];extra['cluster_marker'][j]=r['marker'];extra['cluster_major_proposal'][j]=call['major'];extra['cluster_celltype_proposal'][j]=call['celltype'];extra['cluster_annotation_status'][j]=call['status'];extra['cluster_membership_SHA256'][j]=r['stable_membership_SHA256'];extra['cluster_state_program_names'][j]=tuple(k for k,v in r['descriptive_state_programs'].items() if v['descriptively_enriched']);extra['cluster_view_status'][j]='audited_numeric_cluster_RNA_proposal_not_final_taxonomy'
   extra['cluster_T_evidence_mode'][j]=r.get('assay_aware_T_evidence_mode','original_native_module_assay_coverage_in_proposal')
  assert not set(extra)&set(result['columns']);result['columns'].update(extra);return result
 def read_counts(self,source,positions,genes):return self.baseline.read_counts(source,positions,genes)

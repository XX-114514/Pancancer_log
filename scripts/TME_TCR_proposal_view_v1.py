"""Read accepted additional TCR evidence without replacing current consensus."""
import hashlib,json
from pathlib import Path
import numpy as np
from TME_scoped_annotation_view_v1 import TMEScopedAnnotationView
from joined_native_atlas_delivery_v1 import positions_checked

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

class TMETCRProposalView:
    def __init__(self,path):
        a=load(path)
        if a['schema']!='accepted_TME_TCR_boundary_proposal_v1':raise ValueError('Accepted proposal contract required')
        for binding in a['bindings'].values():
            if sha(binding['path'])!=binding['SHA256']:raise ValueError('Bound proposal evidence changed')
        self.contract=a;self.base=TMEScopedAnnotationView(a['bindings']['working_acceptance']['path'])
        self.root=Path(a['bindings']['proposal_index']['path']).parent;self.index={r['source']:r for r in load(a['bindings']['proposal_index']['path'])}
        if set(self.index)!=set(self.base.sources):raise ValueError('Proposal scope differs from accepted working annotation')
        self.sources=self.base.sources;self.cache={}
    def source(self,source):
        if source not in self.index:raise KeyError('Source outside accepted TME proposal scope')
        if self.cache.get('source')==source:return self.cache
        r=self.index[source];p=self.root/r['metadata']
        if sha(p)!=r['metadata_SHA256']:raise ValueError('Proposal source metadata changed')
        m=load(p);ap=self.root/m['arrays']
        if sha(ap)!=m['arrays_SHA256']:raise ValueError('Proposal arrays changed')
        with np.load(ap,allow_pickle=False) as f:a={k:f[k] for k in f.files}
        base=self.base.source(source)
        for k in ['native_row','source_native_position']:np.testing.assert_array_equal(a[k],base['arrays'][k])
        if m['inputs']['working_arrays_SHA256']!=self.base.index[source]['arrays_SHA256']:raise ValueError('Proposal baseline differs')
        self.cache={'source':source,'arrays':a,'metadata':m};return self.cache
    def metadata(self,source,positions):
        s=self.source(source);a=s['arrays'];p=positions_checked(positions,len(a['native_row']));out=self.base.metadata(source,p);c=out['columns'];proposed=a['proposed_T'][p]
        c['additional_TCR_coarse_proposal']=np.where(proposed,'T',None)
        c['additional_TCR_review']=np.asarray(self.contract['review_dictionary'],object)[a['review_code'][p]]
        for key in ['TCR_core','alpha_chain_detected','gamma_chain_detected','competing_core_mask','unassessed_old_competitor_mask']:
            c['additional_TCR__'+key]=a[key][p]
        c['additional_TCR_evidence_version']=np.repeat(self.contract['evidence_version'],len(p))
        out['provenance'].update(additional_TCR_proposal_not_applied_to_consensus=True,changed_membership_markers_required_before_release=True,new_CM_started=False)
        return out
    def read_counts(self,source,positions,genes):return self.base.read_counts(source,positions,genes)

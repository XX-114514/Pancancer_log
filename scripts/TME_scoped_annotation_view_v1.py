"""Read an independently audited TME working annotation over original counts."""
import hashlib,json
from pathlib import Path
import numpy as np
from joined_native_atlas_delivery_v1 import JoinedNativeAtlas,positions_checked

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

class TMEScopedAnnotationView:
    def __init__(self,acceptance):
        a=load(acceptance)
        if a['schema']!='accepted_TME_scoped_working_annotation_v1':raise ValueError('An accepted working-view contract is required')
        for k in ['producer_manifest','independent_audit_terminal']:
            b=a[k]
            if sha(b['path'])!=b['SHA256']:raise ValueError('Working view evidence changed')
        if a['final_module_input_ready'] is not False:raise ValueError('This reader does not certify final module eligibility')
        self.acceptance=a;self.manifest=load(a['producer_manifest']['path'])
        for b in self.manifest['bindings'].values():
            if sha(b['path'])!=b['SHA256']:raise ValueError('Bound source or annotation input changed')
        i=self.manifest['index'];g=self.manifest['group_catalog']
        if sha(i['path'])!=i['SHA256'] or sha(g['path'])!=g['SHA256']:raise ValueError('Annotation index or groups changed')
        self.root=Path(i['path']).parent;self.index={r['source']:r for r in load(i['path'])}
        self.sources=tuple(sorted(self.index));self.groups=load(g['path'])['groups'];self.atlas=JoinedNativeAtlas(self.manifest['joined_manifest']);self.cache={}

    def source(self,source):
        if source not in self.index:raise KeyError('Source is not in the qualified TME working scope')
        if self.cache.get('source')==source:return self.cache
        r=self.index[source];mp=self.root/r['metadata'];npz=self.root/r['arrays']
        if sha(mp)!=r['metadata_SHA256'] or sha(npz)!=r['arrays_SHA256']:raise ValueError('Source annotation changed')
        m=load(mp)
        with np.load(npz,allow_pickle=False) as f:a={k:f[k] for k in f.files}
        groups=sorted([g for g in self.groups if g['source']==source],key=lambda g:g['source_annotation_group_code'])
        if [g['source_annotation_group_code'] for g in groups]!=list(range(len(groups))):raise ValueError('Group lookup axis differs')
        self.cache={'source':source,'metadata':m,'arrays':a,'groups':groups};return self.cache

    def metadata(self,source,positions):
        s=self.source(source);a=s['arrays'];m=s['metadata'];p=positions_checked(positions,len(a['native_row']));native=a['source_native_position'][p]
        original=self.atlas.metadata(source,native);out=original['columns']
        for field in ['consensus_major','consensus_celltype','coarse_evidence','fine_evidence','confidence','subcluster','global_subcluster_id','global_subcluster_number','stable_group_sha256','marker','marker_origin','marker_evidence_type','label_version']:
            if field in out:out['baseline__'+field]=out[field]
        out['source_scope_position']=p
        for name in ['consensus_major','consensus_celltype']:
            out[name]=np.asarray(m['candidate_dictionaries'][name],object)[a[name+'_code'][p]]
        review=np.asarray(m['candidate_dictionaries']['annotation_review'],object)[a['annotation_review_code'][p]]
        out['annotation_review']=review;out['label_version']=np.repeat(m['label_version'],len(p))
        # Old conflict reasons describe the baseline and must not masquerade
        # as the evidence of a newly reviewed working proposal.
        code=a['annotation_review_code'][p]
        out['coarse_evidence']=out['baseline__coarse_evidence'].astype(object).copy()
        out['fine_evidence']=out['baseline__fine_evidence'].astype(object).copy()
        out['coarse_evidence'][code==1]='additional_native_context_pDC_related_working_proposal'
        out['fine_evidence'][code==1]='additional_native_context_not_orthogonal_identity_validation'
        out['fine_evidence'][code==5]='baseline_pDC_RNA_evidence_retained_related_alias'
        out['coarse_evidence'][code==6]='additional_TCR_core_reference_and_local_context_supported'
        out['fine_evidence'][code==6]='coarse_T_RNA_context_only_no_new_fine_identity'
        out['confidence']=np.asarray([str(r)+';'+str(c)+';'+str(f) for r,c,f in zip(review,out['coarse_evidence'],out['fine_evidence'])],object)
        out['author_label']=np.asarray(m['author_label_dictionary'],object)[a['author_label_code'][p]]
        out['fit_pool_candidate_not_final_QC']=a['fit_pool_candidate_not_final_QC'][p]
        out['source_qualification_code']=a['source_qualification_code'][p]
        group=[s['groups'][int(k)] if k>=0 else None for k in a['annotation_group_code'][p]]
        for name in ['global_subcluster_id','global_subcluster_number','stable_group_sha256','subcluster','marker','marker_origin','marker_evidence_type','parent_marker_borrowed','borrow_reason']:
            out[name]=np.asarray([g.get(name) if g is not None else (-1 if name=='global_subcluster_number' else None) for g in group],object)
        out['identity_release_status']=np.repeat('audited_scoped_working_annotation_not_final_taxonomy',len(p))
        original['provenance'].update(label_version=m['label_version'],scope='source_qualified_TME_not_final_QC',
            marker_parent_borrowing_is_explicit=True,subcluster_numbers_are_exact_groups_not_shared_types=True,
            biological_accuracy_not_certified_by_numeric_audit=True,final_module_input_ready=False,
            clinical_or_patient_independence_not_newly_qualified=True)
        return original

    def read_counts(self,source,positions,genes):
        a=self.source(source)['arrays'];p=positions_checked(positions,len(a['native_row']))
        return self.atlas.read_counts(source,a['source_native_position'][p],genes)

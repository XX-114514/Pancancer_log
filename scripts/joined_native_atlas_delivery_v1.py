"""Read audited native counts, guarded labels, states and projections together.

All inputs remain immutable. Exact structural group IDs are not shared cell
types. Scientific eligibility and candidate identity flags stay visible; this
interface does not convert technical acceptance into biological qualification.
"""
import json
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from audit_full_source_local_graphs_v1 import index_at
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection, text_axis_sha256
from summarize_full_native_annotation_evidence_v1 import read_source
from extract_full_source_local_marker_evidence_v1 import gather_csr_rows


def positions_checked(positions, n, limit=4096):
    p = np.asarray(positions)
    if p.ndim != 1 or len(p) > limit or (p.size and p.dtype.kind not in 'iu'):
        raise ValueError('Use a one-dimensional integer block of at most4096rows')
    p = p.astype(np.int64)
    if p.size and (p.min() < 0 or p.max() >= n):
        raise IndexError('Archive position outside source')
    return p


def domain_lookup(domain, positions, n):
    """Preserve query order/duplicates; absent model rows are explicitly masked."""
    q = np.asarray(domain)
    if q.ndim != 1 or q.dtype.kind not in 'iu' or (q.size and
            (q.min() < 0 or q.max() >= n or np.any(np.diff(q) <= 0))):
        raise ValueError('Model domain must be a sorted unique archive subaxis')
    p = positions_checked(positions, n)
    idx = np.searchsorted(q, p)
    available = idx < len(q)
    available[available] &= q[idx[available]] == p[available]
    return idx, available


def completed(stage):
    stage = Path(stage)
    t = json.loads((stage/'provenance/terminal_v1.json').read_text())
    s = json.loads((stage/'provenance/supervisor_terminal_v1.json').read_text())
    if s['actual_child_exit_code'] != 0 or s['wall_timeout']:
        raise ValueError('Successful real stage termination required')
    return t


class JoinedNativeAtlas:
    def __init__(self, manifest):
        self.spec = json.loads(Path(manifest).read_text())
        if self.spec['schema'] != 'joined_native_atlas_delivery_v1':
            raise ValueError('Unsupported joined delivery')
        for b in self.spec['bindings'].values():
            if sha256_file(b['path']) != b['sha256']:
                raise ValueError('Bound delivery evidence changed')
        self.native = NativeArchiveCollection(self.spec['archive_manifest'])
        self.catalog = Path(self.spec['catalog'])
        self.units = Path(self.spec['unit_views'])
        self.evidence = Path(self.spec['native_evidence'])
        self.review = Path(self.spec['local_review'])
        self.ci = index_at(self.catalog, 'full_guarded_identity_catalog_index_v1.json')
        self.ui = index_at(self.units, 'archive_unit_view_index_v1.json')
        self.ei = index_at(self.evidence, 'full_archive_reference_index_v1.json')
        self.li = index_at(self.review, 'full_local_identity_review_index_v1.json')
        self.sources = tuple(sorted(self.native.archive_partitions))
        if any(set(i) != set(self.sources) for i in [self.ci,self.ui,self.ei,self.li]):
            raise ValueError('Joined delivery source scope differs')
        fp = self.catalog/'outputs/guarded_subset_identity_catalog_v1.json'
        if sha256_file(fp) != completed(self.catalog)['catalog_sha256']:
            raise ValueError('Guarded group catalog changed')
        self.groups = {r['global_subcluster_number']:r for r in json.loads(fp.read_text())['groups']}
        op = Path(self.spec['qualified_projection'])
        if sha256_file(op) != completed(self.spec['projection_qualifier'])['output_sha256']:
            raise ValueError('Qualified projection manifest changed')
        self.projections = json.loads(op.read_text())
        self.pi = {(r['source'],r['model_id']):r for r in self.projections['records']}
        if len(self.pi) != len(self.projections['records']):
            raise ValueError('Repeated projection key')
        self._cache = {};self._verified = set()

    def _verify(self,path,digest):
        key=(str(path),digest)
        if key not in self._verified:
            if sha256_file(path)!=digest:raise ValueError('Bound sidecar changed')
            self._verified.add(key)

    def source(self,ds):
        if ds not in self.sources:raise KeyError(ds)
        if self._cache.get('source')==ds:return self._cache
        self._cache.clear();self.native._fields.clear()
        cm,ca=read_source(self.catalog,self.ci[ds]);em,ea=read_source(self.evidence,self.ei[ds])
        umpath=self.units/'outputs'/self.ui[ds]['metadata'];uapath=self.units/'outputs'/self.ui[ds]['arrays']
        self._verify(umpath,self.ui[ds]['metadata_sha256']);self._verify(uapath,self.ui[ds]['arrays_sha256'])
        um=json.loads(umpath.read_text())
        with np.load(uapath,allow_pickle=False) as f:ua={k:f[k] for k in f.files}
        lp=self.review/'outputs'/self.li[ds]['output'];self._verify(lp,self.li[ds]['output_sha256']);lm=json.loads(lp.read_text())
        rows=self.native.rows(ds);ids=self.native.native_ids(ds)
        for side in [ca,ea,ua]:np.testing.assert_array_equal(side['native_row'],rows)
        for m in [cm,em,um,lm]:
            if m['native_ids_sha256']!=text_axis_sha256(ids):raise ValueError('Joined native identity axis differs')
        self._cache.update(source=ds,cm=cm,ca=ca,em=em,ea=ea,um=um,ua=ua,lm=lm,rows=rows,ids=ids)
        return self._cache

    def metadata(self,ds,positions):
        s=self.source(ds);n=len(s['rows']);p=positions_checked(positions,n)
        ca,ea,ua,um=s['ca'],s['ea'],s['ua'],s['um']
        author=self.native.author_labels(ds)
        out={'source':np.repeat(ds,len(p)),'native_cell_id':s['ids'][p],
             'source_native_position':p,'native_row':s['rows'][p],
             'author_label':author['values'][p],
             'local_subcluster_id':np.asarray([ds+':parent:'+str(k) if k>=0 else '' for k in ca['parent_source_local_group_code'][p]],object),
             'fit_pool_candidate':ca['fit_pool_candidate'][p],
             'global_subcluster_number':ca['global_subcluster_number'][p],
             'guarded_source_subset_code':ca['guarded_source_subset_code'][p],
             'label_version':np.repeat(s['cm']['label_version'],len(p)),
             'identity_release_status':np.repeat('audited_guarded_candidate_not_authoritative_taxonomy',len(p))}
        dictionaries=s['cm']['candidate_dictionaries']
        for field in ['consensus_major','consensus_celltype','coarse_evidence','fine_evidence']:
            out[field]=np.asarray(dictionaries[field],object)[ca[field+'_code'][p]]
        out['confidence']=np.asarray([str(a)+';'+str(b) for a,b in zip(out['coarse_evidence'],out['fine_evidence'])],object)
        local={g['source_local_group_code']:g['local_major_proposal'] for g in s['lm']['groups']}
        out['local_major']=np.asarray([local.get(int(k),'unresolved_no_fitted_local_group') for k in ca['parent_source_local_group_code'][p]],object)
        groups=[self.groups.get(int(k)) for k in out['global_subcluster_number']]
        for field in ['global_subcluster_id','subcluster','marker','marker_evidence_type','marker_origin','stable_group_sha256','group_expression_evidence_status']:
            out[field]=np.asarray([g.get(field) if g else None for g in groups],object)
        unit=ua['native_expression_unit_code'][p];block=ua['analysis_block_code'][p]
        out['native_expression_unit']=ua['expression_unit_dictionary'][unit]
        out['analysis_block_kind']=np.asarray([um['analysis_block_dictionary'][k]['kind'] for k in block],object)
        out['analysis_block_original_key']=np.asarray([um['analysis_block_dictionary'][k]['original_key'] for k in block],object)
        out['patient_review_status']=np.asarray([(um['unit_records'][k]['inherited_frozen_unit_review'] or {}).get('patient_review_status','unknown') for k in unit],object)
        for j,view in enumerate(um['view_order']):
            out['view_status__'+view]=np.asarray([um['view_status_decoder'][str(k)] for k in ua['unit_view_status_code'][unit,j]],object)
        for j,state in enumerate(s['em']['state_order']):
            out['RNA_state__'+state]=ea['state_native_logCP10K'][p,j]
        # Preserve original fields verbatim, including CNV evaluability, prior
        # QC and treatment flags. Names/values alone never grant eligibility.
        from audit_native_malignancy_evidence_registry_v1 import PAIRS,STATUS,SCORES
        from audit_source_analysis_roles_v1 import STATUS_FIELDS
        carried=list(dict.fromkeys([f for pair in PAIRS for f in pair]+STATUS+SCORES+STATUS_FIELDS+
            ['patient_id','patient_uid','donor_id','sample_uid','sampleID','sample_id','library_id','treatment','treatment_stage','treatment_phase','response','timepoint','sample_class','cancer_type']))
        fields=set(self.native.field_names(ds));present=[f for f in carried if f in fields]
        for field in present:out['native__'+field]=self.native.field(ds,field)[p]
        return {'columns':out,'provenance':{'author_label_status':author['status'],
            'confidence_scale':'categorical_evidence_not_probability',
            'state_measurement_coverage':s['em']['state_measurement_coverage'],
            'RNA_states_are_descriptive_not_functional_identity':True,
            'original_carried_fields':present,'absent_carried_fields':[f for f in carried if f not in fields],
            'malignancy_status':'original_evidence_separate_from_lineage_no_new_call',
            'group_numbers_are_not_shared_biological_celltypes':True,
            'new_clinical_or_composition_eligibility_assigned':False}}

    def read_counts(self,ds,positions,genes):
        rows=self.native.rows(ds);p=positions_checked(positions,len(rows))
        if len(genes)!=len(set(genes)):raise ValueError('Repeated requested feature')
        lookup=self.native.feature_lookup(ds)
        cols=[next(iter(lookup[g])) if len(lookup.get(g,set()))==1 else None for g in genes]
        if len({c for c in cols if c is not None})!=sum(c is not None for c in cols):raise ValueError('Requested feature aliases collide')
        key=self.native.legacy_config['counts_key'] if self.native.is_legacy(ds) else self.native.primary.partitions[ds]['counts_key']
        full=gather_csr_rows(self.native.store(ds)[key],rows[p])
        panel=sp.hstack([full[:,c] if c is not None else sp.csr_matrix((len(p),1),dtype=full.dtype) for c in cols],format='csr') if cols else sp.csr_matrix((len(p),0),dtype=full.dtype)
        return {'counts':panel,'measured':np.asarray([c is not None for c in cols]),'genes':list(genes),
            'native_total_counts':np.asarray(full.sum(axis=1,dtype=np.float64)).ravel(),
            'native_row':rows[p],'count_origin':key,'missing_features_are_not_biological_zeros':True}

    def projection(self,ds,model_id,positions):
        n=len(self.native.rows(ds));p=positions_checked(positions,n)
        entry=self.pi[(ds,model_id)];record=entry['base_projection_record']
        if record['n_model_domain_rows']==0:return {'model_domain_available':np.zeros(len(p),bool),'arrays':{},'empty_model_domain':True}
        base=Path(self.projections['base_projection'])/'outputs';ap=base/record['arrays'];self._verify(ap,record['arrays_sha256'])
        with np.load(ap,allow_pickle=False) as f:values={k:f[k] for k in f.files}
        domain=values['source_native_position'];idx,available=domain_lookup(domain,p,n)
        np.testing.assert_array_equal(values['native_row'],self.native.rows(ds)[domain])
        exact=np.zeros(len(domain),bool)
        for patch in entry['exact_neighbor_overlays']:
            path=Path(patch['repair_arrays']);self._verify(path,patch['repair_arrays_sha256'])
            with np.load(path,allow_pickle=False) as f:v={k:f[k] for k in f.files}
            q,found=domain_lookup(domain,v['source_native_position'],n)
            if not found.all():raise ValueError('Exact overlay outside model domain')
            basis=patch['basis'];original_available=values[basis+'_reference_distance_available']
            np.testing.assert_array_equal(domain[q],domain[original_available])
            values[basis+'_reference_distance']=values[basis+'_reference_distance'].astype(np.float64)
            values[basis+'_reference_distance'][q]=v['exact_nearest_other_reference_distance']
            values[basis+'_outside_reference_radius_q99'][q]=v['outside_original_reference_radius_q99']
            exact[q]=True
        result={}
        for name,v in values.items():
            if v.ndim<1 or len(v)!=len(domain):raise ValueError('Unexpected projection subaxis')
            result[name]=np.ma.masked_all((len(p),)+v.shape[1:],dtype=v.dtype)
            result[name][available]=v[idx[available]]
        marker=np.zeros(len(p),bool);marker[available]=exact[idx[available]]
        return {'model_domain_available':available,'arrays':result,'exact_neighbor_overlay_applied':marker,
            'radius_flags_are_geometric_not_identity_confidence':True,'clinical_or_biological_transfer_qualified':False}

"""Outcome-free native-RNA marker/author-state correspondence; arrays RAM only.

All source cells, all author states and unresolved references retained. No
classifier learning, new clinical inference, or cross-assay QC equivalence.
"""
from collections import defaultdict
from datetime import datetime,timezone
from itertools import product
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import urllib.request

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.sparse import csr_matrix
import observed_T_marker_evidence_v1 as evidence
import frozen_CM_network_contexts_v1 as network_context

PARENTS=['dual_regulatory_evidence','costate_regulatory_sensitivity']
CONTEXTS=['full_current_source','CM02_regulatory_naming_anchor','CM02_cytotoxic_naming_anchor']
ANCHORS=['39697d077b6cb8bf','c3f761a9da61f578']
NODE_CONTEXTS={CONTEXTS[0]:None,CONTEXTS[1]:{ANCHORS[0]},CONTEXTS[2]:{ANCHORS[1]}}
PHASES=['pre','post','post2']
DEPTH=['all_native_filtered','low_native_capture','mid_native_capture','high_native_capture']


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def module_at(name,path,expected):
    assert sha(path)==expected
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m


def native_counts(z,selected,columns,human,measured,checkpoint):
    raw=z['layers/counts'];shape=tuple(raw.attrs['shape'])
    ptr=raw['indptr'][:];assert len(ptr)==shape[0]+1 and ptr[0]==0 and np.all(np.diff(ptr)>=0)
    assert ptr[-1]==raw['data'].shape[0]==raw['indices'].shape[0]
    counts=np.zeros((len(selected),len(columns)),np.int64);totals=np.zeros(len(selected),np.int64)
    native_detected=np.zeros(len(selected),np.int64)
    slot=np.full(shape[1],-1,np.int64);slot[columns]=np.arange(len(columns))
    completed=0
    for start in range(0,shape[0],5000):
        stop=min(shape[0],start+5000);lo,hi=np.searchsorted(selected,[start,stop])
        if lo==hi:continue
        a,b=int(ptr[start]),int(ptr[stop]);ind=np.asarray(raw['indices'][a:b],np.int64)
        values=np.asarray(raw['data'][a:b]);assert np.isfinite(values).all() and (values>=0).all() and (values==np.floor(values)).all()
        assert values.max(initial=0)<2**63;values=values.astype(np.int64)
        assert (ind>=0).all() and (ind<shape[1]).all() and values.sum(dtype=np.float64)<2**53
        pp=ptr[start:stop+1]-a
        block=csr_matrix((values,ind,pp),shape=(stop-start,shape[1]),copy=True)
        block.sum_duplicates();block.eliminate_zeros();block.sort_indices()
        rows=selected[lo:hi]-start
        current=block[rows]
        observed=current[:,columns].toarray().astype(np.int64)
        full=np.asarray(current[:,human].sum(axis=1)).ravel().astype(np.int64)
        native_detected[lo:hi]=np.diff(current[:,human].indptr)
        # Raw CSR entry arithmetic, separate from sparse library selection.
        ri=np.repeat(np.arange(stop-start,dtype=np.int64),np.diff(pp))
        wanted=slot[ind]>=0
        direct=np.bincount(ri[wanted]*len(columns)+slot[ind[wanted]],weights=values[wanted],minlength=(stop-start)*len(columns)).reshape(stop-start,len(columns))
        valid=human[ind]
        direct_total=np.bincount(ri[valid],weights=values[valid],minlength=stop-start)
        source_rows=np.zeros(stop-start,bool);source_rows[rows]=True
        assert not np.any((~measured[ind])&(values>0)&source_rows[ri])
        np.testing.assert_array_equal(observed,direct[rows].astype(np.int64))
        np.testing.assert_array_equal(full,direct_total[rows].astype(np.int64))
        counts[lo:hi]=observed;totals[lo:hi]=full;completed+=int(hi-lo)
        checkpoint(dict(stage='native_marker_counts_and_independent_raw_CSR',completed_source_cells=completed))
    assert completed==len(selected)
    return counts,totals,native_detected


def source_aggregate(frame,reference,depth,flags,parent,labels):
    patient_names=sorted(frame.patient.unique());patient=pd.Index(patient_names).get_indexer(frame.patient)
    phase=pd.Index(PHASES).get_indexer(frame.timepoint);label=pd.Index(labels).get_indexer(reference)
    assert (patient>=0).all() and (phase>=0).all() and (label>=0).all()
    n=len(patient_names);size=len(PHASES)*len(labels)*n
    group=(phase*len(labels)+label)*n+patient
    output=[]
    for context,di in product(CONTEXTS,range(4)):
        chosen=np.ones(len(frame),bool)
        if NODE_CONTEXTS[context] is not None:chosen&=frame.node_id.isin(NODE_CONTEXTS[context]).to_numpy()
        if di:chosen&=depth==(di-1)
        den=np.bincount(group[chosen],minlength=size).reshape(len(PHASES),len(labels),n)
        nums=np.column_stack([np.bincount(group[chosen],weights=flags[chosen,j],minlength=size) for j in range(flags.shape[1])]).astype(np.int64).reshape(len(PHASES),len(labels),n,-1)
        assert (nums>=0).all() and (nums<=den[...,None]).all()
        for pi,li in product(range(len(PHASES)),range(len(labels))):
            denominator=den[pi,li];good=denominator>=20
            metrics={}
            for j,name in enumerate(evidence.METRICS):
                numerator=nums[pi,li,:,j];total=int(denominator.sum());positive=int(numerator.sum())
                balanced=numerator[good]/denominator[good]
                metrics[name]=dict(positive_cell_observations=positive,
                    pooled_cell_observation_fraction=positive/total if total else None,
                    case_balanced_rate_min_q25_median_q75_max=np.quantile(balanced,[0,.25,.5,.75,1]).tolist() if good.sum()>=5 else None,
                    case_balanced_status='DESCRIPTIVE_REFERENCE_CORRESPONDENCE' if good.sum()>=5 else 'INSUFFICIENT_CASES_ATLEAST20_REFERENCE_CELLS')
            output.append(dict(source='GSE286827',phase=PHASES[pi],context=context,capture_scope=DEPTH[di],
                regulatory_parent=parent,reference_label=labels[li],cell_observations=int(denominator.sum()),
                reviewed_cases_with_observations=int((denominator>0).sum()),reviewed_cases_atleast20_reference_cells=int(good.sum()),
                measurements=metrics,author_reference_not_independent_biological_truth=True))
        assert int(den.sum())==int(chosen.sum())
    return output


def main():
    global CONTEXTS,NODE_CONTEXTS
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert out.parent.parent==Path('/dev/shm') and not out.exists();out.mkdir(mode=0o700)
    save=lambda n,v:(out/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    fine=Path(os.environ['PANCANCER_FINE_LOADER']);ft=Path(os.environ['PANCANCER_FINE_SCOPE_OUTPUT'])/'private_terminal.json'
    metadata=json.loads(ft.read_text())
    loader=module_at('pinned_fine_reference',fine,metadata['source_SHA256'])
    atlas=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']).parents[1]/'runs/20261013_immunotherapy_atlas_v1'
    NODE_CONTEXTS,membership=network_context.resolve(atlas,CONTEXTS,ANCHORS)
    CONTEXTS=list(NODE_CONTEXTS)
    save('private_original_network_member_identity.json',membership)
    save('plan_before_native_marker_reference_correspondence.json',dict(started_UTC=datetime.now(timezone.utc).isoformat(),
        source_SHA256=sha(__file__),recipe_SHA256=sha(evidence.__file__),network_context_SHA256=sha(network_context.__file__),fine_loader_SHA256=sha(fine),
        complete_original_network_contexts=True,network_membership_SHA256=membership['network_membership_SHA256'],
        posthoc_full_member_description_not_a_new_effect_family=True,
        fine_reference_terminal_SHA256=sha(ft),question='Same molecule-positive marker recipes versus every author state; no state learning or clinical tests.',
        native_QC='All source retained authority RNA cells; within-library native human RNA value-tertiles. Never treat full-native UMI20/50 as Xenium small-panel molecule QC.',
        all_genes=evidence.GENES,all_parents=PARENTS,all_contexts=CONTEXTS,all_phases=PHASES,all_capture_scopes=DEPTH,
        author_CD4_labels=loader.CD4,author_CD8_labels=loader.CD8,
        unresolved=['no_unique_fine_reference','overlapping_fine_reference'],
        case_balance='Reviewed patient/phase counts aggregate all libraries; >=20 reference cells and >=5 cases for descriptive rate quantiles; no P/FDR.',
        all_cell_scalar_checks=True,raw_CSR_independent_counts=True,no_cell_patient_expression_arrays_saved=True,
        source_transport='explicit official direct HTTPS, fixed original metadata payload SHA; no inherited proxy or fallback',
        no_clinical_or_TCR_values_selected=True,no_original_authority_modified=True))
    urllib.request.install_opener(urllib.request.build_opener(urllib.request.ProxyHandler({})))
    tab,frame,unique,shared,matrices,table=loader.load_inputs()
    atlas=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']).parents[1]/'runs/20261013_immunotherapy_atlas_v1'
    z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r')
    dataset=np.asarray(read_elem(z['obs/dataset'])).astype(str);selected=np.flatnonzero(dataset=='GSE286827')
    assert len(selected)==len(frame) and np.array_equal(frame.position,np.arange(len(frame)))
    var=read_elem(z['var']);sources=list(map(str,read_elem(z['uns/source_order'])))
    measured=np.asarray(read_elem(z['varm/measured_in_source']),bool)[:,sources.index('GSE286827')]
    human=measured & var.feature_identity_status.astype(str).eq('stable_Ensembl_ID').to_numpy()
    lookup=defaultdict(set)
    for col in np.flatnonzero(human):
        r=var.iloc[col]
        for value in {str(var.index[col]),str(r.gene_symbol)}|set(str(r.source_gene_aliases).split('|')):
            lookup[value].add(int(col))
    columns=[]
    for gene in evidence.GENES:
        assert len(lookup[gene])==1;columns.append(next(iter(lookup[gene])))
    assert len(set(columns))==len(columns)
    release=atlas/'manifests/immune_authority_release_public.json'
    save('private_native_input_identity.json',dict(source='GSE286827',authority_release_SHA256=sha(release),
        author_table_SHA256=sha(table),unique_measured_marker_columns=len(columns),native_source_feature_count=int(measured.sum()),
        native_stable_human_feature_count=int(human.sum()),all_required_markers_uniquely_mapped=True,
        native_column_map_SHA256=hashlib.sha256(json.dumps(dict(zip(evidence.GENES,columns)),sort_keys=True).encode()).hexdigest()))
    checkpoint=lambda v:save('private_checkpoint.json',v)
    counts,totals,detected=native_counts(z,selected,columns,human,measured,checkpoint)
    frame['patient']=frame.library_id.map(tab.set_index('library_id').patient_registry)
    assert frame.patient.notna().all() and set(frame.timepoint)<=set(PHASES)
    reference=np.full(len(frame),'no_unique_fine_reference',object)
    reference[list(shared)]='overlapping_fine_reference'
    reference[unique.position.to_numpy(int)]=unique.label.to_numpy(str)
    assert not set(unique.position)&shared
    labels=loader.CD4+loader.CD8+['no_unique_fine_reference','overlapping_fine_reference']
    assert set(reference)<=set(labels)
    depth=np.full(len(frame),-1,np.int8)
    for lib in frame.library_id.unique():
        positions=np.flatnonzero(frame.library_id.eq(lib).to_numpy())
        cuts=np.unique(np.quantile(totals[positions],[1/3,2/3]))
        depth[positions]=np.searchsorted(cuts,totals[positions],side='left')
    assert (depth>=0).all() and (depth<=2).all()
    summaries=[];background=[];parent_rates=[]
    for parent in PARENTS:
        flags=evidence.vector_evidence(counts,np.ones(len(counts),bool),parent)
        independent=evidence.scalar_evidence(counts,np.ones(len(counts),bool),parent)
        np.testing.assert_array_equal(flags,independent)
        assert (flags[:,1:5].sum(axis=1)==1).all()
        summaries.extend(source_aggregate(frame,reference,depth,flags,parent,labels))
        parent_rates.extend(network_context.parent_profiles(frame,depth,flags,parent,NODE_CONTEXTS,evidence.METRICS,PHASES,DEPTH))
        for (phase,major),idx in frame.groupby(['timepoint','major']).groups.items():
            ids=np.asarray(list(idx),int)
            background.append(dict(source='GSE286827',phase=phase,current_major=major,regulatory_parent=parent,
                cell_observations=len(ids),native_human_UMI_min_q25_median_q75_max=np.quantile(totals[ids],[0,.25,.5,.75,1]).tolist(),
                native_detected_human_genes_min_q25_median_q75_max=np.quantile(detected[ids],[0,.25,.5,.75,1]).tolist(),
                observed_marker_counts=dict(zip(evidence.METRICS,flags[ids].sum(axis=0).astype(int).tolist()))))
    assert len(summaries)==2*len(CONTEXTS)*4*3*15
    prior_root=Path(os.environ['PANCANCER_NATIVE_PARENT_OUTPUT'])
    prior_terminal=json.loads((prior_root/'private_terminal.json').read_text())
    prior_file=prior_root/'private_all1080_native_marker_reference_opportunities.json'
    assert sha(prior_file)==prior_terminal['files'][prior_file.name]['SHA256']
    base=['full_current_source','CM02_regulatory_naming_anchor','CM02_cytotoxic_naming_anchor']
    assert [r for r in summaries if r['context'] in base]==json.loads(prior_file.read_text())
    save('private_all_network_context_parent_state_rates.json',parent_rates)
    save('private_all_full_network_native_marker_reference_opportunities.json',summaries)
    save('private_all_source_major_capture_background.json',background)
    save('private_terminal.json',dict(status='COMPLETE_FULL_FROZEN_NETWORK_MEMBER_NATIVE_MARKER_REFERENCE_CORRESPONDENCE',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),recipe_SHA256=sha(evidence.__file__),
        network_context_SHA256=sha(network_context.__file__),all_original1080_anchor_opportunities_exactly_replayed=True,
        original_network_member_identity_SHA256=sha(out/'private_original_network_member_identity.json'),
        full_member_contexts=len(CONTEXTS),source_opportunities=len(summaries),parent_state_profiles=len(parent_rates),
        all_native_marker_and_human_totals_raw_CSR_replayed=True,all_source_cells_both_parents_scalar_checked=True,
        all_reference_labels_unmatched_overlap_depth_phase_and_anchor_opportunities_preserved=True,
        current_source_reference_join_reuses_prior_exact_library_barcode_and_unresolved_overlap=True,
        no_cross_assay_absolute_capture_equivalence_claim=True,no_new_clinical_pair_or_spatial_tests=True,
        no_cell_patient_or_native_arrays_saved=True,
        files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status='COMPLETE_FULL_FROZEN_NETWORK_MEMBER_NATIVE_MARKER_REFERENCE_CORRESPONDENCE',terminal_SHA256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

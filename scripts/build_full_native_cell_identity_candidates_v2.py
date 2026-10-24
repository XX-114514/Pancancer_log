"""All-archive cell-level guarded candidates using the accepted pilot v2 rule.

Categorical proposals retain native/group/reference conflict and model exposure.
They are not a final taxonomy, calibrated confidence or clinical eligibility.
"""
import argparse
import collections
import datetime
import gc
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
from annotation_representation_pilot_v1 import save
from audit_full_source_local_graphs_v1 import index_at
from cell_identity_confidence_v2 import FINE_TYPES,FINE_PANEL,decide_cell
from scoped_low_reference_identity_v1 import scoped_low_identity
from extract_full_source_local_marker_evidence_v1 import storage_gate
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
from native_lineage_routing_v1 import FAMILIES,reference_family
from summarize_full_native_annotation_evidence_v1 import read_source

MAJORS=sorted(set(['not_available','unresolved','T','NK','B','Plasma',*FAMILIES]))
FIELDS=['consensus_major','consensus_celltype','coarse_evidence','fine_evidence','direct_fine_candidate']
FLAGS=['group_major_disagreement','group_type_disagreement','fine_discovery_fit_available']


def save_arrays(path,arrays):
    with path.open('xb') as f:np.savez_compressed(f,**arrays)
    with np.load(path,allow_pickle=False) as f:
        if set(f.files)!=set(arrays):raise ValueError('Candidate serialization keys differ')
        for k,v in arrays.items():np.testing.assert_array_equal(f[k],v)


def unpack(code,names):return {name for j,name in enumerate(names) if int(code)&(1<<j)}


def evaluate_signature(row):
    native,reference,native_major,reference_major,group_major,group_type,own,ref,positive,finefit=map(int,row)
    native_major,reference_major=MAJORS[native_major],MAJORS[reference_major]
    return decide_cell(unpack(native,FAMILIES),unpack(reference,FAMILIES),
       None if native_major=='not_available' else native_major,None if reference_major=='not_available' else reference_major,
       MAJORS[group_major],CELLTYPES[group_type],unpack(own,FINE_TYPES),unpack(ref,FINE_TYPES),bool(positive),bool(finefit))


CELLTYPES=sorted(set(['unresolved',*MAJORS,*FINE_TYPES]))


def main(root):
    began=time.monotonic();cfg=json.loads((root/'config/full_cell_candidate_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path'])!=binding['sha256']:raise ValueError('Full cell candidate input changed')
    for filename,count in [('test_cell_identity_confidence_v2.py',10),('test_scoped_low_reference_identity_v1.py',3)]:
        test=subprocess.run([sys.executable,'-B',str(root/'tests'/filename)],capture_output=True,text=True,timeout=120)
        (root/'provenance'/(filename+'.stderr.txt')).write_text(test.stderr)
        if test.returncode or ('Ran %d tests'%count) not in test.stderr or 'skipped' in test.stderr.lower():raise ValueError('Actual full cell guard/scope counterexamples failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':13,'no_skips':True})
    review=Path(cfg['review']);graphs=Path(cfg['graphs']);reference=Path(cfg['reference'])
    index=index_at(review,'full_local_identity_review_index_v1.json');gi=index_at(graphs,'full_native_source_local_graph_index_v1.json');ri=index_at(reference,'full_archive_reference_index_v1.json')
    reader=NativeArchiveCollection(cfg['archive_manifest'])
    if any(set(v)!=set(reader.archive_partitions) for v in [index,gi,ri]):raise ValueError('Full candidate scope differs')
    models=json.loads((reference/'outputs/reference_models_v1.json').read_text())
    refmaps={m:[reference_family(m,label) for label in v['labels']] for m,v in models.items()}
    lowmaps=[scoped_low_identity(label) for label in models['immune_low']['labels']]
    major_index={s:j for j,s in enumerate(MAJORS)};type_index={s:j for j,s in enumerate(CELLTYPES)}
    popcount=np.array([bin(j).count('1') for j in range(1<<len(FAMILIES))],np.uint8)
    reports=[];total_replays=0
    for ds in sorted(index,key=lambda d:(index[d]['n_archive_rows'],d)):
        storage_gate(root,cfg);entry=index[ds];path=review/'outputs'/entry['output']
        if sha256_file(path)!=entry['output_sha256']:raise ValueError('Candidate group identity source changed')
        reviewed=json.loads(path.read_text());gm,ga=read_source(graphs,gi[ds]);rm,rv=read_source(reference,ri[ds]);rows=reader.rows(ds);n=len(rows)
        for side in [ga,rv]:np.testing.assert_array_equal(rows,side['native_row'])
        if any(m['native_ids_sha256']!=text_axis_sha256(reader.native_ids(ds)) for m in [reviewed,gm,rm]):raise ValueError('Candidate native identity differs')
        native=rv['native_family_bitmask'];ref=np.zeros(n,np.uint32);own=np.zeros(n,np.uint16);ref_fine=np.zeros(n,np.uint16)
        panel_order=rm['native_panel_bit_order']
        panel=lambda name:(rv['native_panel_bitmask']&(1<<panel_order.index(name)))!=0
        for j,t in enumerate(FINE_TYPES):
            if t=='CD4_T':
                k=rm['review_genes'].index('CD4')
                positive=((rv['native_gene_detection_packbits'][:,k//8]>>(k%8))&1).astype(bool)&panel('T') if rm['gene_measured'][k] else np.zeros(n,bool)
            else:positive=panel(FINE_PANEL[t])
            own|=positive.astype(np.uint16)<<j
        low_code=rv['immune_low__unique_code'];ref_major=np.full(n,major_index['not_available'],np.int16)
        for code,(major,fine) in enumerate(lowmaps):
            member=low_code==code
            if major is not None:ref_major[member]=major_index[major]
            value=sum(1<<j for j,t in enumerate(FINE_TYPES) if t in fine);ref_fine[member]=value
        for model in rm['source_reference_models']:
            codes=rv[model+'__unique_code'];npos=rv[model+'__n_positive_labels'];scope=rv['immune_reference_scope'] if model.startswith('immune_') else rv['organ_reference_scope']
            if not np.array_equal(codes>=0,npos==1) or np.any((npos>=0)&~scope):raise ValueError('Candidate reference scope changed')
            for j,family in enumerate(FAMILIES):ref|=np.isin(codes,[k for k,f in enumerate(refmaps[model]) if f==family]).astype(np.uint32)<<j
        if np.any(native>=(1<<len(FAMILIES))) or np.any(ref>=(1<<len(FAMILIES))):raise ValueError('Candidate evidence bit exceeds dictionary')
        native_major=np.full(n,major_index['not_available'],np.int16)
        for j,family in enumerate(FAMILIES):native_major[native==(1<<j)]=major_index[family]
        tn=native==(1<<FAMILIES.index('T_NK'));bp=native==(1<<FAMILIES.index('B_Plasma'))
        native_major[tn]=np.where(panel('T')[tn],major_index['T'],major_index['NK'])
        native_major[bp]=np.where(panel('Plasma')[bp],major_index['Plasma'],major_index['B'])
        parent=ga['source_local_group_code'];group_major=np.full(n,major_index['unresolved'],np.int16);group_type=np.full(n,type_index['unresolved'],np.int16)
        for group in reviewed['groups']:
            member=parent==group['source_local_group_code'];group_major[member]=major_index[group['local_major_proposal']];group_type[member]=type_index[group['local_celltype_proposal']]
        # Only equal complete evidence signatures share an evaluation. The rule
        # sees the same arguments as the accepted per-cell pilot implementation.
        signatures=np.column_stack([native,ref,native_major,ref_major,group_major,group_type,own,ref_fine,rv['positive_RNA'],ga['fine_graph_fitted']]).astype(np.int32)
        unique,inverse=np.unique(signatures,axis=0,return_inverse=True)
        decisions=[evaluate_signature(row) for row in unique]
        dictionaries={};output={'native_row':rows,'parent_source_local_group_code':parent,
          'fit_pool_candidate':ga['fit_pool_candidate'],'fine_graph_fitted':ga['fine_graph_fitted'],
          'native_family_bitmask':native,'reference_family_bitmask':ref,'own_fine_bitmask':own,'reference_fine_bitmask':ref_fine,
          'positive_native_RNA':rv['positive_RNA']}
        decoded={}
        for field in FIELDS:
            names,code=np.unique(np.asarray([d[field] for d in decisions]),return_inverse=True)
            if len(names)>32767:raise ValueError('Candidate categorical capacity exceeded')
            output[field+'_code']=code[inverse].astype(np.int16);dictionaries[field]=names.tolist();decoded[field]=names[code[inverse]]
        for field in FLAGS:output[field]=np.asarray([d[field] for d in decisions],bool)[inverse]
        hard=(popcount[native]>1)|(popcount[ref]>1)|((popcount[native]==1)&(popcount[ref]==1)&(native!=ref))|~rv['positive_RNA']
        if np.any(decoded['consensus_major'][hard]!='unresolved'):raise ValueError('Group label overwrote hard cell conflict')
        gamma=decoded['consensus_celltype']=='T_gamma_delta_RNA_candidate'
        if np.any(gamma&~panel('gamma_delta_T_related')):raise ValueError('Gamma candidate lacks own RNA evidence')
        np.testing.assert_array_equal(output['fine_discovery_fit_available'],ga['fine_graph_fitted'])
        selected=np.unique(np.r_[np.linspace(0,n-1,min(n,4096),dtype=np.int64),np.flatnonzero(hard)[:8],np.flatnonzero(gamma)[:8],np.flatnonzero(~ga['fit_pool_candidate'])[:8]])
        for i in selected:
            # Explicit individual rule calls validate signature expansion axes.
            a=signatures[i];nm=MAJORS[int(a[2])];rlm=MAJORS[int(a[3])]
            direct=decide_cell(unpack(a[0],FAMILIES),unpack(a[1],FAMILIES),None if nm=='not_available' else nm,None if rlm=='not_available' else rlm,
              MAJORS[int(a[4])],CELLTYPES[int(a[5])],unpack(a[6],FINE_TYPES),unpack(a[7],FINE_TYPES),bool(a[8]),bool(a[9]))
            for field in FIELDS:
                if direct[field]!=decoded[field][i]:raise ValueError('Individual versus compressed cell rule differs')
            for field in FLAGS:
                if direct[field]!=output[field][i]:raise ValueError('Individual versus compressed flag differs')
        total_replays+=len(selected)
        npz=root/'outputs'/(ds+'_guarded_identity_candidates_v1.npz');save_arrays(npz,output)
        old=reader.frozen_labels(ds);crosswalk=[]
        # A compact count crosswalk preserves archive exclusions and old labels.
        before=np.asarray(old['fine']).astype(str);after=decoded['consensus_celltype']
        bnames,bcode=np.unique(before,return_inverse=True);anames,acode=np.unique(after,return_inverse=True)
        packed,count=np.unique(bcode.astype(np.int64)*len(anames)+acode,return_counts=True)
        for value,number in zip(packed,count):crosswalk.append({'old_fine':str(bnames[value//len(anames)]),'guarded_celltype_candidate':str(anames[value%len(anames)]),'n_archive_records':int(number)})
        group_counts=[]
        for group in reviewed['groups']:
            code=group['source_local_group_code'];member=parent==code;hist=dict(collections.Counter(after[member]));group_counts.append({'source_local_group_code':code,'stable_parent_group_sha256':group['stable_local_group_sha256'],
             'guarded_celltype_counts':{str(k):int(v) for k,v in hist.items()},'parent_display_marker':group['display_marker_proposal'],
             'marker_origin':'parent_source_local_group_not_recomputed_for_guarded_subset','subgroup_marker_review_required':len(hist)>1})
        meta={'source':ds,'n_archive_rows':n,'native_ids_sha256':text_axis_sha256(reader.native_ids(ds)),
          'dictionaries':dictionaries,'native_family_order':FAMILIES,'fine_candidate_order':FINE_TYPES,
          'n_unique_evidence_signatures':len(unique),'n_individual_rule_calls_replayed':len(selected),
          'n_fit_pool_candidates':int(ga['fit_pool_candidate'].sum()),'n_archive_excluded_rows_retained':int((~ga['fit_pool_candidate']).sum()),
          'group_identity_review_record':entry,'graph_record':gi[ds],'native_reference_record':ri[ds],
          'label_version':cfg['label_version'],'fields_are_guarded_candidates_not_final_consensus_release':True,
          'confidence_is_evidence_provenance_not_calibrated_probability':True,'common_representation_and_final_taxonomy_pending':True,
          'author_labels_original_counts_and_RNA_states_retained_in_bound_native_inputs':True,
          'malignancy_is_separate_and_not_inferred_here':True,'clinical_and_composition_eligibility_assigned':False,
          'archive_excluded_rows_can_only_receive_direct_candidates_without_fitted_group_inference':True,
          'field_counts':{field:{str(k):int(v) for k,v in collections.Counter(decoded[field]).items()} for field in FIELDS},
          'old_fine_candidate_crosswalk':crosswalk,'parent_group_guarded_counts':group_counts,
          'global_subcluster_id_not_assigned':True,'own_marker_for_new_guarded_subsets_not_computed':True}
        mp=npz.with_suffix('.json');save(mp,meta)
        reports.append({'source':ds,'n_archive_rows':n,'arrays':npz.name,'arrays_sha256':sha256_file(npz),'metadata':mp.name,'metadata_sha256':sha256_file(mp)})
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'completed_sources':len(reports)})
        reader._fields.clear()
        del rv,ga,reviewed,output,decoded,unique,inverse,signatures,decisions,old,before,after,group_counts,crosswalk
        gc.collect()
    output=root/'outputs/full_native_cell_identity_candidate_index_v1.json';save(output,reports)
    save(root/'provenance/terminal_v1.json',{'status':'FULL_NATIVE_CELL_GUARDED_CANDIDATES_COMPLETE_INDEPENDENT_AUDIT_AND_TAXONOMY_PENDING',
      'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(reports),
      'n_archive_rows':sum(r['n_archive_rows'] for r in reports),'n_individual_rule_calls_replayed':total_replays,
      'index_sha256':sha256_file(output),'review_producer_terminal_sha256':sha256_file(review/'provenance/terminal_v1.json'),
      'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
      'new_final_identity':False,'common_model_complete':False,'native_counts_reread':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise

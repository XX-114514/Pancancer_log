"""All-row independent reconstruction of guarded cell candidate inputs/outputs.

Audits computation and scope. It neither calibrates identity confidence nor
turns parent-group markers into markers of newly guarded subsets.
"""
import argparse
import collections
import datetime
import gc
import json
from pathlib import Path
import resource
import time

import numpy as np
from annotation_representation_pilot_v1 import save
from audit_full_source_local_graphs_v1 import index_at
from cell_identity_confidence_v2 import FINE_TYPES,FINE_PANEL,low_reference_identity,decide_cell
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
from native_lineage_routing_v1 import FAMILIES,reference_family,local_family
from summarize_full_native_annotation_evidence_v1 import read_source

FIELDS=['consensus_major','consensus_celltype','coarse_evidence','fine_evidence','direct_fine_candidate']
FLAGS=['group_major_disagreement','group_type_disagreement','fine_discovery_fit_available']


def hist(values):return {str(k):int(v) for k,v in collections.Counter(values).items()}


def labels_from_bits(bits,names):return {name for j,name in enumerate(names) if int(bits)&(1<<j)}


def main(root):
    began=time.monotonic();spec=json.loads((root/'config/full_cell_candidate_audit_spec_v1.json').read_text());producer=Path(spec['producer'])
    cp=producer/'config/full_cell_candidate_spec_v1.json'
    if sha256_file(cp)!=spec['producer_config_sha256']:raise ValueError('Full cell producer configuration changed')
    cfg=json.loads(cp.read_text());request=json.loads((producer/'config/supervisor_request_v1.json').read_text())
    for name,digest in request['script_SHA256'].items():
        if sha256_file(producer/name)!=digest:raise ValueError('Full cell deployed code changed')
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Full cell input binding changed')
    index=index_at(producer,'full_native_cell_identity_candidate_index_v1.json');review=Path(cfg['review']);graphs=Path(cfg['graphs']);reference=Path(cfg['reference'])
    groups=index_at(review,'full_local_identity_review_index_v1.json');gi=index_at(graphs,'full_native_source_local_graph_index_v1.json');ri=index_at(reference,'full_archive_reference_index_v1.json')
    reader=NativeArchiveCollection(cfg['archive_manifest'])
    if any(set(v)!=set(reader.archive_partitions) for v in [index,groups,gi,ri]):raise ValueError('Full cell audit source scope differs')
    models=json.loads((reference/'outputs/reference_models_v1.json').read_text())
    maps={name:[reference_family(name,label) for label in model['labels']] for name,model in models.items()}
    low=[]
    for code,label in enumerate(models['immune_low']['labels']):
        expected_family=maps['immune_low'][code]
        if expected_family is None:major,fine=None,set()
        else:
            major,fine=low_reference_identity(label)
            if major is not None and local_family(major)!=expected_family:raise ValueError('Accepted low reference vocabulary inconsistent')
        low.append((major,fine))
    major_names=sorted(set(['not_available','unresolved','T','NK','B','Plasma',*FAMILIES]));type_names=sorted(set([*major_names,*FINE_TYPES]));midx={m:j for j,m in enumerate(major_names)};tidx={m:j for j,m in enumerate(type_names)}
    pop=np.array([bin(i).count('1') for i in range(1<<len(FAMILIES))],np.uint8);reports=[];total_cases=0
    for ds in sorted(index,key=lambda d:(index[d]['n_archive_rows'],d)):
        cm,ca=read_source(producer,index[ds]);gm,ga=read_source(graphs,gi[ds]);rm,rv=read_source(reference,ri[ds]);path=review/'outputs'/groups[ds]['output']
        if sha256_file(path)!=groups[ds]['output_sha256']:raise ValueError('Full cell local review changed')
        reviewed=json.loads(path.read_text());rows=reader.rows(ds);n=len(rows)
        for side in [ca,ga,rv]:np.testing.assert_array_equal(rows,side['native_row'])
        if any(m['native_ids_sha256']!=text_axis_sha256(reader.native_ids(ds)) for m in [cm,gm,rm,reviewed]):raise ValueError('Full cell native identity axis differs')
        for key in ['fit_pool_candidate','fine_graph_fitted']:np.testing.assert_array_equal(ca[key],ga[key])
        np.testing.assert_array_equal(ca['parent_source_local_group_code'],ga['source_local_group_code']);np.testing.assert_array_equal(ca['positive_native_RNA'],rv['positive_RNA'])
        np.testing.assert_array_equal(ca['native_family_bitmask'],rv['native_family_bitmask'])
        if cm['n_archive_rows']!=n or cm['n_fit_pool_candidates']!=int(ga['fit_pool_candidate'].sum()) or cm['n_archive_excluded_rows_retained']!=int((~ga['fit_pool_candidate']).sum()):raise ValueError('Full cell archive/model denominator differs')
        if cm['group_identity_review_record']!=groups[ds] or cm['graph_record']!=gi[ds] or cm['native_reference_record']!=ri[ds]:raise ValueError('Full cell source bindings differ')
        if not cm['fields_are_guarded_candidates_not_final_consensus_release'] or cm['clinical_and_composition_eligibility_assigned'] or not cm['global_subcluster_id_not_assigned']:raise ValueError('Full cell output scope promoted')
        decoded={}
        for field in FIELDS:
            dictionary=cm['dictionaries'][field];code=ca[field+'_code']
            if len(dictionary)!=len(set(dictionary)) or code.shape!=(n,) or not np.issubdtype(code.dtype,np.integer) or np.any(code<0) or np.any(code>=len(dictionary)):raise ValueError('Full cell categorical dictionary differs')
            decoded[field]=np.asarray(dictionary)[code]
            if hist(decoded[field])!=cm['field_counts'][field]:raise ValueError('Full cell source summary counts differ')
        expected_ref=np.zeros(n,np.uint32)
        for model in rm['source_reference_models']:
            code=rv[model+'__unique_code'];npos=rv[model+'__n_positive_labels'];scope=rv['immune_reference_scope'] if model.startswith('immune_') else rv['organ_reference_scope']
            if not np.array_equal(code>=0,npos==1) or np.any((npos>=0)&~scope):raise ValueError('Full cell reference scoring scope differs')
            for j,family in enumerate(FAMILIES):
                valid_codes=np.flatnonzero(np.asarray(maps[model],object)==family)
                expected_ref[np.isin(code,valid_codes)]|=1<<j
        np.testing.assert_array_equal(ca['reference_family_bitmask'],expected_ref)
        panel_order=rm['native_panel_bit_order'];bits=rv['native_panel_bitmask']
        has=lambda name:((bits>>panel_order.index(name))&1).astype(bool)
        own=np.zeros(n,np.uint16);ref_fine=np.zeros(n,np.uint16)
        for j,name in enumerate(FINE_TYPES):
            if name=='CD4_T':
                col=rm['review_genes'].index('CD4')
                positive=(((rv['native_gene_detection_packbits'][:,col//8]>>(col%8))&1).astype(bool)&has('T')) if rm['gene_measured'][col] else np.zeros(n,bool)
            else:positive=has(FINE_PANEL[name])
            own[positive]|=1<<j
        np.testing.assert_array_equal(ca['own_fine_bitmask'],own)
        low_major=np.full(n,midx['not_available'],np.int16);low_code=rv['immune_low__unique_code'];unscoped=np.zeros(n,bool)
        for code,(major,fine) in enumerate(low):
            where=low_code==code
            if major is not None:low_major[where]=midx[major]
            for j,name in enumerate(FINE_TYPES):
                if name in fine:ref_fine[where]|=1<<j
            if maps['immune_low'][code] is None:unscoped|=where
        np.testing.assert_array_equal(ca['reference_fine_bitmask'],ref_fine)
        native=rv['native_family_bitmask'];nm=np.full(n,midx['not_available'],np.int16)
        if np.any(native>=len(pop)):raise ValueError('Native family evidence code exceeds capacity')
        for j,family in enumerate(FAMILIES):nm[native==(1<<j)]=midx[family]
        tn=native==(1<<FAMILIES.index('T_NK'));bp=native==(1<<FAMILIES.index('B_Plasma'))
        nm[tn&has('T')]=midx['T'];nm[tn&~has('T')]=midx['NK'];nm[bp&has('Plasma')]=midx['Plasma'];nm[bp&~has('Plasma')]=midx['B']
        parent=ga['source_local_group_code'];gmajor=np.full(n,midx['unresolved'],np.int16);gtype=np.full(n,tidx['unresolved'],np.int16)
        for group in reviewed['groups']:
            where=parent==group['source_local_group_code'];gmajor[where]=midx[group['local_major_proposal']];gtype[where]=tidx[group['local_celltype_proposal']]
        signature=np.column_stack([native,expected_ref,nm,low_major,gmajor,gtype,own,ref_fine,rv['positive_RNA'],ga['fine_graph_fitted']]).astype(np.int32)
        unique,inverse=np.unique(signature,axis=0,return_inverse=True);decisions=[]
        for row in unique:
            a,b,c,d,e,f,g,h,i,j=map(int,row);native_major=major_names[c];reference_major=major_names[d]
            decisions.append(decide_cell(labels_from_bits(a,FAMILIES),labels_from_bits(b,FAMILIES),
              None if native_major=='not_available' else native_major,None if reference_major=='not_available' else reference_major,
              major_names[e],type_names[f],labels_from_bits(g,FINE_TYPES),labels_from_bits(h,FINE_TYPES),bool(i),bool(j)))
        if len(unique)!=cm['n_unique_evidence_signatures']:raise ValueError('Full cell evidence signature count differs')
        for field in FIELDS:np.testing.assert_array_equal(np.asarray([v[field] for v in decisions])[inverse],decoded[field])
        for field in FLAGS:np.testing.assert_array_equal(np.asarray([v[field] for v in decisions],bool)[inverse],ca[field])
        hard=(pop[native]>1)|(pop[expected_ref]>1)|((pop[native]==1)&(pop[expected_ref]==1)&(native!=expected_ref))|~rv['positive_RNA']
        if np.any(decoded['consensus_major'][hard]!='unresolved'):raise ValueError('Full cell hard conflict guard violated')
        gamma=decoded['consensus_celltype']=='T_gamma_delta_RNA_candidate'
        if np.any(gamma&~has('gamma_delta_T_related')):raise ValueError('Gamma candidate lacks own RNA support')
        old=reader.frozen_labels(ds);before=np.asarray(old['fine']).astype(str);after=decoded['consensus_celltype']
        bn,bc=np.unique(before,return_inverse=True);an,ac=np.unique(after,return_inverse=True);keys,count=np.unique(bc.astype(np.int64)*len(an)+ac,return_counts=True)
        expected_crosswalk={(str(bn[k//len(an)]),str(an[k%len(an)])):int(v) for k,v in zip(keys,count)}
        actual_crosswalk={(r['old_fine'],r['guarded_celltype_candidate']):r['n_archive_records'] for r in cm['old_fine_candidate_crosswalk']}
        if expected_crosswalk!=actual_crosswalk or len(actual_crosswalk)!=len(cm['old_fine_candidate_crosswalk']):raise ValueError('Full cell old-new crosswalk differs')
        bycode={g['source_local_group_code']:g for g in reviewed['groups']}
        if len(cm['parent_group_guarded_counts'])!=len(bycode):raise ValueError('Full cell parent group coverage differs')
        for g in cm['parent_group_guarded_counts']:
            code=g['source_local_group_code'];counts=hist(after[parent==code]);original=bycode[code]
            if g['guarded_celltype_counts']!=counts or g['stable_parent_group_sha256']!=original['stable_local_group_sha256'] or g['parent_display_marker']!=original['display_marker_proposal'] or g['marker_origin']!='parent_source_local_group_not_recomputed_for_guarded_subset' or g['subgroup_marker_review_required']!=(len(counts)>1):raise ValueError('Full cell parent marker provenance differs')
        total_cases+=len(unique);reports.append({'source':ds,'n_archive_rows':n,'n_evidence_signatures_replayed':len(unique),
          'n_unscoped_low_reference_rows':int(unscoped.sum()),'n_hard_conflict_rows':int(hard.sum()),'n_gamma_candidates_checked':int(gamma.sum()),
          'source_arrays_sha256':index[ds]['arrays_sha256'],'source_metadata_sha256':index[ds]['metadata_sha256']})
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'completed_sources':len(reports)})
        reader._fields.clear();del ca,ga,rv,reviewed,decoded,signature,unique,inverse,decisions,old,before,after;gc.collect()
    output=root/'outputs/full_native_cell_candidate_audit_v1.json';save(output,{'sources':reports,'all_archive_row_fields_reconstructed':True,
      'n_evidence_signatures_replayed':total_cases,'not_independent_biological_identity_truth':True,'native_counts_reread':False})
    save(root/'provenance/terminal_v1.json',{'status':'FULL_NATIVE_CELL_CANDIDATE_COMPUTATIONAL_AUDIT_COMPLETE_NOT_FINAL_LABELS','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'n_sources':len(reports),'n_archive_rows':sum(r['n_archive_rows'] for r in reports),'n_evidence_signatures_replayed':total_cases,
      'output_sha256':sha256_file(output),'producer_terminal_sha256':sha256_file(producer/'provenance/terminal_v1.json'),
      'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'new_final_identity':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise

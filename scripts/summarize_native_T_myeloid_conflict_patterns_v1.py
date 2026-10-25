"""Descriptive gene-detection signatures in native T/myeloid co-support rows.

Does not infer doublets, contamination or a new identity. Presence and missing
measurement are separate; stored reference labels remain candidate evidence.
"""
import argparse,collections,datetime,json,time,resource
from pathlib import Path
import numpy as np
from annotation_representation_pilot_v1 import save
from audit_full_source_local_graphs_v1 import index_at
from extract_full_source_local_marker_evidence_v1 import storage_gate
from full_context_atlas_v1 import sha256_file
from native_lineage_routing_v1 import FAMILIES
from summarize_full_native_annotation_evidence_v1 import read_source

GENES=['CD3D','CD3E','TRAC','CD247','LST1','LYZ','FCER1G','CTSS','AIF1','GNLY','KLRD1']


def main(root):
    start=time.monotonic();cfg=json.loads((root/'config/native_conflict_pattern_spec_v1.json').read_text())
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Conflict signature input binding changed')
    candidate=Path(cfg['candidate']);evidence=Path(cfg['native_evidence']);ci=index_at(candidate,'full_native_cell_identity_candidate_index_v1.json');ri=index_at(evidence,'full_archive_reference_index_v1.json')
    if set(ci)!=set(ri):raise ValueError('Conflict review source sets differ')
    report=[];total=collections.Counter();nrows=0
    for ds in sorted(ci):
        storage_gate(root,cfg);cm,c=read_source(candidate,ci[ds]);em,e=read_source(evidence,ri[ds]);np.testing.assert_array_equal(c['native_row'],e['native_row']);np.testing.assert_array_equal(c['native_family_bitmask'],e['native_family_bitmask'])
        both=((c['native_family_bitmask']&(1<<FAMILIES.index('T_NK')))!=0)&((c['native_family_bitmask']&(1<<FAMILIES.index('Myeloid')))!=0)
        covered={g:bool(em['gene_measured'][em['review_genes'].index(g)]) for g in GENES};codes=np.zeros(int(both.sum()),np.uint16)
        for j,g in enumerate(GENES):
            pos=em['review_genes'].index(g)
            if covered[g]:codes|=(((e['native_gene_detection_packbits'][both,pos//8]>>(pos%8))&1).astype(np.uint16)<<j)
        reference=c['reference_family_bitmask'][both].astype(np.int64);unique,counts=np.unique(np.column_stack([codes,reference]),axis=0,return_counts=True)
        patterns=[{'detection_pattern_bits':int(a),'reference_family_bits':int(b),'n_archive_records':int(n)} for (a,b),n in zip(unique,counts)]
        field=np.asarray(cm['dictionaries']['consensus_major'])[c['consensus_major_code'][both]]
        if len(field) and np.any(field!='unresolved'):raise ValueError('Hard native cross-family guard unexpectedly overridden')
        count={'n_archive_rows':len(both),'n_native_T_and_myeloid_cosupport':int(both.sum())}
        anchors=(1<<GENES.index('LST1'))|(1<<GENES.index('LYZ'));cd3=(1<<GENES.index('CD3D'))|(1<<GENES.index('CD3E'))
        if covered['LST1'] and covered['LYZ']:
            count['n_cosupport_with_detected_LST1_or_LYZ']=int(((codes&anchors)!=0).sum());count['n_cosupport_without_detected_LST1_or_LYZ_but_genes_measured']=int(((codes&anchors)==0).sum())
        else:count['LST1_LYZ_comparison']='measurement_incomplete_not_negative'
        if covered['CD3D'] and covered['CD3E']:count['n_cosupport_with_detected_CD3D_or_CD3E']=int(((codes&cd3)!=0).sum())
        total.update({k:v for k,v in count.items() if isinstance(v,int)});nrows+=len(both)
        report.append({'source':ds,**count,'gene_detection_bit_order':GENES,'gene_measured':covered,'reference_family_bit_order':FAMILIES,'patterns':patterns,'native_ids_sha256':cm['native_ids_sha256'],'no_quantitative_RNA_abundance_or_background_measured_here':True})
        save(root/'checkpoints'/(ds+'.json'),{'source':ds,'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':'native_T_myeloid_co_support_signature_counts_complete'})
    output=root/'outputs/native_T_myeloid_conflict_patterns_v1.json';save(output,{'sources':report,'aggregate_counts':dict(total),'status':'DESCRIPTIVE_NATIVE_DETECTION_REVIEW_NOT_NEW_CELL_IDENTITY','family_cosupport_is_not_doublet_or_contamination_diagnosis':True,'reference_labels_not_independent_truth':True,'new_clinical_or_composition_qualification':False})
    save(root/'provenance/terminal_v1.json',{'status':'FULL_NATIVE_T_MYELOID_DETECTION_SIGNATURE_REVIEW_COMPLETE','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(report),'n_archive_rows_checked':nrows,'aggregate_counts':dict(total),'output_sha256':sha256_file(output),'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'new_authoritative_labels':False,'native_counts_reread':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise

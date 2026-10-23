"""Distinguish group mixing from cell-level multi-lineage RNA support.

These observations do not diagnose doublets, ambient RNA or malignant status.
Legacy count mapping is checked against original per-source H5AD on fixed rows.
"""
import argparse
import collections
import datetime
import gzip
import json
from pathlib import Path
import resource
import time

import h5py
import numpy as np
from anndata._io.specs import read_elem

from annotation_representation_pilot_v1 import save
from full_context_atlas_v1 import FullContextAtlas,sha256_file,_h5_index
from native_panel_aggregate_v1 import gather_csr_rows
from native_annotation_review_v1 import REVIEW_GENES,marker_calls


def cell_families(calls):
    out={k:calls[p].copy() for k,p in {'T':'T','NK':'NK_related','B':'B','Plasma':'Plasma','Myeloid':'Myeloid',
         'Mast':'Mast','Endothelial':'Endothelial','Fibroblast':'Fibroblast','Mural':'Mural','Epithelial':'Epithelial',
         'Glial':'Glial_Oligodendrocyte_related','Astrocyte':'Astrocyte_related','Neuron':'Neuron_related','Erythroid':'Erythroid_related'}.items()}
    out['Myeloid']|=calls['DC_related']|calls['Neutrophil_related']
    out['Endothelial']|=calls['Lymphatic_related']
    out['NK']&=~out['T'];out['B']&=~out['Plasma']
    return out


def legacy_original_check(reader,ds,positions):
    if not reader.is_legacy(ds):return {'status':'not_applicable_to_primary_partition','original_primary_QC_evidence_not_repeated':True}
    source_path=Path(reader.legacy_config['native_source_root'])/(ds+'_symbol_counts.h5ad')
    canonical=reader.store(ds);key=reader.legacy_config['counts_key']
    gse=np.asarray(read_elem(canonical['obs/gse_id'])).astype(str)
    all_source=np.flatnonzero(gse==ds);selected=reader.rows(ds)[positions]
    fixed=np.unique(np.linspace(0,len(selected)-1,min(32,len(selected)),dtype=int))
    native_rows=np.searchsorted(all_source,selected[fixed])
    if not np.array_equal(all_source[native_rows],selected[fixed]):raise ValueError('Native source row bridge failed')
    genes=_h5_index(source_path,'var');lookup=reader.feature_lookup(ds)
    columns=np.array([next(iter(lookup[g])) for g in genes],np.int64)
    n_values=0
    with h5py.File(source_path,'r') as f:
        native=f['X']
        if native.attrs.get('encoding-type') not in ['csr_matrix',b'csr_matrix']:raise ValueError('Original source X is not qualified CSR')
        for local,row in zip(native_rows,selected[fixed]):
            a,b=map(int,native['indptr'][int(local):int(local)+2])
            values=np.zeros(len(genes),np.float64)
            np.add.at(values,np.asarray(native['indices'][a:b]),np.asarray(native['data'][a:b],np.float64))
            current=gather_csr_rows(canonical[key],np.array([row]))
            np.testing.assert_array_equal(np.asarray(current[:,columns].toarray()).ravel(),values)
            if float(current.sum())!=float(values.sum()):raise ValueError('Original versus canonical complete RNA denominator differs')
            n_values+=len(genes)
    return {'status':'original_H5AD_to_canonical_all_native_gene_readback_pass_on_fixed_rows',
            'n_rows':len(fixed),'n_native_gene_values':n_values,'row_selection':'evenly_spaced_pilot_positions_before_RNA_inspection',
            'all_cells_in_source_replayed':False}


def main(root,producer,audit):
    start=time.monotonic()
    if json.loads((audit/'provenance/supervisor_terminal_v1.json').read_text())['actual_child_exit_code']!=0:raise ValueError('Native review audit incomplete')
    cfg=json.loads((producer/'config/review_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path'])!=binding['sha256']:raise ValueError('Frozen annotation input changed')
    reader=FullContextAtlas(cfg['manifest']);pilot=Path(cfg['pilot']);discovery=Path(cfg['discovery'])
    inventory=json.loads((pilot/'outputs/pilot_input_inventory_v1.json').read_text())
    review_index=json.loads((producer/'outputs/source_annotation_review_index_v1.json').read_text())
    meta=json.loads((producer/'outputs/reference_models_v1.json').read_text())
    with np.load(discovery/'outputs'/('candidate_partitions_'+str(cfg['seed'])+'.npz'),allow_pickle=False) as f:
        parts={m:f[m] for m in ['A','B','C']}
    offset=0;summaries=[]
    for inv,entry in zip(inventory,review_index):
        ds=inv['dataset'];positions=np.load(pilot/inv['selector_file'],allow_pickle=False);n=len(positions)
        if ds!=entry['source'] or sha256_file(pilot/inv['selector_file'])!=inv['selector_sha256']:raise ValueError('Source or selector changed')
        review_path=producer/'outputs'/entry['output']
        if sha256_file(review_path)!=entry['output_sha256']:raise ValueError('Review table changed')
        with gzip.open(review_path,'rt') as f:review=json.load(f)
        original=legacy_original_check(reader,ds,positions)
        lookup=reader.feature_lookup(ds);measured=np.array([len(lookup.get(g,set()))==1 for g in REVIEW_GENES])
        cols=np.array([next(iter(lookup[g])) if measured[j] else 0 for j,g in enumerate(REVIEW_GENES)])
        z=reader.store(ds);key=reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
        pieces=[];core_pieces=[];dc_pieces=[];neutrophil_pieces=[]
        for lo in range(0,n,2048):
            raw=gather_csr_rows(z[key],reader.rows(ds)[positions[lo:lo+2048]])
            markers=raw[:,cols].toarray();markers[:,~measured]=0
            calls=marker_calls(markers>0,REVIEW_GENES,measured)
            core_pieces.append(calls['Myeloid']);dc_pieces.append(calls['DC_related']);neutrophil_pieces.append(calls['Neutrophil_related'])
            families=cell_families(calls)
            names=sorted(families);pieces.append(np.column_stack([families[k] for k in names]))
        support=np.vstack(pieces);del pieces
        core=np.concatenate(core_pieces);dc=np.concatenate(dc_pieces);neutrophil=np.concatenate(neutrophil_pieces)
        unique,patterns,counts=np.unique(support,axis=0,return_inverse=True,return_counts=True)
        with np.load(producer/'outputs'/(ds+'_reference_proposals_v1.npz'),allow_pickle=False) as f:
            low_unique=f['immune_low__unique_code'];low_positive=f['immune_low__n_positive_labels']
        pattern_rows=[]
        for j,pattern in enumerate(unique):
            member=patterns==j
            refs=collections.Counter(meta['immune_low']['labels'][c] for c in low_unique[member] if c>=0)
            pattern_rows.append({'native_supported_families':[name for name,v in zip(names,pattern) if v],
               'n_pilot_cells':int(counts[j]),'unique_immune_reference_counts_not_truth':dict(refs),
               'n_reference_unassigned':int((low_positive[member]==0).sum()),'n_reference_multiple':int((low_positive[member]>1).sum())})
        rows=[]
        for record in review['groups']:
            if record['annotation_evidence_status']!='cross_lineage_marker_conflict_review':continue
            member=parts[record['method']][offset:offset+n]==record['partition_code']
            family_columns=[names.index(g) for g in record['native_supported_families']]
            local=support[member][:,family_columns].sum(axis=1)
            rows.append({'source':ds,'method':record['method'],'source_local_group_id':record['source_local_group_id'],
                         'native_supported_families':record['native_supported_families'],'n_pilot_cells':int(member.sum()),
                         'n_cell_support_multiple_conflicting_families':int((local>=2).sum()),
                         'n_cell_support_exactly_one_conflicting_family':int((local==1).sum()),
                         'n_cell_support_no_conflicting_family':int((local==0).sum()),
                         'n_myeloid_core_panel':int(core[member].sum()),
                         'n_myeloid_family_from_neutrophil_panel_only':int((neutrophil[member]&~core[member]&~dc[member]).sum()),
                         'n_myeloid_family_from_DC_panel_only':int((dc[member]&~core[member]&~neutrophil[member]).sum())})
        payload={'source':ds,'original_count_mapping_check':original,'cell_support_patterns':pattern_rows,
                 'conflicting_group_cell_support':rows,'not_doublet_or_ambient_or_malignancy_diagnosis':True,'new_labels_assigned':False}
        path=root/'outputs'/(ds+'_lineage_conflicts_v1.json')
        save(path,payload)
        aggregate=[]
        for method in ['A','B','C']:
            group=collections.defaultdict(lambda:collections.Counter())
            for r in rows:
                if r['method']!=method:continue
                key='|'.join(r['native_supported_families']);group[key]['n_groups']+=1
                for name in ['n_pilot_cells','n_cell_support_multiple_conflicting_families','n_cell_support_exactly_one_conflicting_family','n_cell_support_no_conflicting_family','n_myeloid_core_panel','n_myeloid_family_from_neutrophil_panel_only','n_myeloid_family_from_DC_panel_only']:group[key][name]+=r[name]
            aggregate.append({'method':method,'family_pair_summaries':{k:dict(v) for k,v in group.items()}})
        summaries.append({'source':ds,'n_pilot_cells':n,'n_cell_multifamily':int((support.sum(axis=1)>1).sum()),
                          'top_cell_patterns':sorted(pattern_rows,key=lambda r:-r['n_pilot_cells'])[:8],
                          'original_count_mapping_check':original,'method_conflict_summary':aggregate,'output_sha256':sha256_file(path)})
        offset+=n
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds})
        if sum(p.stat().st_size for p in root.rglob('*') if p.is_file())>16*1024**2:raise ValueError('Lineage diagnostic output cap exceeded')
    save(root/'outputs/source_lineage_conflict_summary_v1.json',summaries)
    save(root/'provenance/terminal_v1.json',{'status':'SOURCE_LINEAGE_CONFLICT_DIAGNOSTIC_COMPLETE_NOT_NEW_IDENTITY_OR_QC_EXCLUSION',
         'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(summaries),
         'n_pilot_cells':offset,'output_sha256':sha256_file(root/'outputs/source_lineage_conflict_summary_v1.json'),
         'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['output','producer','audit']:p.add_argument('--'+name,required=True)
    a=p.parse_args();root=Path(a.output)
    try:main(root,Path(a.producer),Path(a.audit))
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)})
        raise

"""Independent vectorized guards, crosswalk audit and direct native fine replay.

No producer decision function is called. Native fine signals in the smallest
preselected source are reconstructed from direct CSR spans and gene booleans.
"""
import argparse
import collections
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import resource
import time
import numpy as np
import scipy.sparse as sp
from annotation_representation_pilot_v1 import save
from full_context_atlas_v1 import FullContextAtlas, sha256_file
from native_panel_aggregate_v1 import LINEAGE_PANELS
from native_annotation_review_v1 import TYPE_PANELS

FAMILIES=['T_NK','B_Plasma','Myeloid','Mast','Endothelial','Fibroblast','Mural','Epithelial','Glial','Astrocyte','Neuron','Erythroid']
FINE=['Treg','CD4_T','CD8_T','T_gamma_delta_RNA_candidate','Neutrophil','pDC','cDC1','cDC2_related','Macrophage','Monocyte']
FINE_PARENT=['T']*4+['Myeloid']*6


def family(major):
    if major in {'T','NK','T_NK'}:return 'T_NK'
    if major in {'B','Plasma','B_Plasma'}:return 'B_Plasma'
    return major if major in FAMILIES else 'unresolved'


def compatible(a,b):
    return a==b or (a in {'T_NK','B_Plasma'} and family(b)==a) or (b in {'T_NK','B_Plasma'} and family(a)==b)


def direct_rows(group,requested):
    unique,inverse=np.unique(requested,return_inverse=True);pieces=[];i=0
    while i<len(unique):
        j=i+1
        while j<len(unique) and unique[j]-unique[i]<12000:j+=1
        lo,hi=int(unique[i]),int(unique[j-1])+1;ptr=np.asarray(group['indptr'][lo:hi+1],np.int64)
        first,last=int(ptr[0]),int(ptr[-1])
        block=sp.csr_matrix((np.asarray(group['data'][first:last]),np.asarray(group['indices'][first:last]),ptr-first),shape=(hi-lo,int(group.attrs['shape'][1])))
        pieces.append(block[unique[i:j]-lo]);i=j
    return sp.vstack(pieces,format='csr')[inverse]


def main(root,producer,previous):
    started=time.monotonic();sup=json.loads((producer/'provenance/supervisor_terminal_v1.json').read_text())
    if sup['actual_child_exit_code']!=0 or sup['wall_timeout']:raise ValueError('Confidence producer incomplete')
    cfg=json.loads((producer/'config/cell_confidence_spec_v2.json').read_text())
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Frozen confidence input changed')
    term=json.loads((producer/'provenance/terminal_v1.json').read_text());idxp=producer/'outputs/cell_identity_confidence_index_v2.json'
    if sha256_file(idxp)!=term['output_sha256']:raise ValueError('Confidence index changed')
    index=json.loads(idxp.read_text());old_index={r['source']:r for r in json.loads((previous/'outputs/cell_identity_confidence_index_v1.json').read_text())}
    reader=FullContextAtlas(cfg['manifest']);pilot=Path(cfg['pilot']);disc=Path(cfg['discovery']);annot=Path(cfg['annotation'])
    inv=json.loads((pilot/'outputs/pilot_input_inventory_v1.json').read_text());audit_source=min(inv,key=lambda r:(r['n_selected_cells'],r['dataset']))['dataset']
    aindex={r['source']:r for r in json.loads((annot/'outputs/source_annotation_review_index_v1.json').read_text())}
    with np.load(disc/'outputs/candidate_partitions_2026102301.npz',allow_pickle=False) as f:parts={m:f[m] for m in ['A','B','C']}
    with np.load(disc/'outputs/derived_axis_and_model_mask.npz',allow_pickle=False) as f:
        source=f['source'];positive_all=f['modeled'];fit={'A':f['A_fine_discovery_eligible'],'B':positive_all,'C':f['C_fine_discovery_eligible']}
    reports=[];offset=0;native_values=0;direct_count=0;hard_count=0;corrected_count=0
    for entry,item in zip(index,inv):
        ds=entry['source'];n=item['n_selected_cells'];sl=slice(offset,offset+n)
        if ds!=item['dataset'] or not np.all(source[sl]==ds):raise ValueError('Confidence source alignment differs')
        positions=np.load(pilot/item['selector_file'],allow_pickle=False)
        if sha256_file(pilot/item['selector_file'])!=item['selector_sha256']:raise ValueError('Pilot selector changed')
        for name,digest in [('sidecar','sidecar_sha256'),('metadata','metadata_sha256')]:
            if sha256_file(producer/'outputs'/entry[name])!=entry[digest]:raise ValueError('Confidence output changed')
        with np.load(producer/'outputs'/entry['sidecar'],allow_pickle=False) as f:data={k:f[k] for k in f.files}
        with gzip.open(producer/'outputs'/entry['metadata'],'rt') as f:meta=json.load(f)
        old_entry=old_index[ds]
        with np.load(previous/'outputs'/old_entry['sidecar'],allow_pickle=False) as f:old_data={k:f[k] for k in f.files}
        with gzip.open(previous/'outputs'/old_entry['metadata'],'rt') as f:old_meta=json.load(f)
        def decode(field,old=False):
            d,md=(old_data,old_meta) if old else (data,meta)
            code=d[field+'_code'];labels=np.array(md['categorical_dictionaries'][field])
            if len(code)!=n or code.min()<0 or code.max()>=len(labels):raise ValueError('Categorical output bounds')
            return labels[code]
        for field in meta['categorical_dictionaries']:decode(field)
        if any(len(v)!=n for v in data.values()):raise ValueError('Sidecar row count differs')
        np.testing.assert_array_equal(data['native_row'],reader.rows(ds)[positions]);np.testing.assert_array_equal(data['positive_native_RNA'],positive_all[sl])
        z=reader.store(ds);ids=reader.field(ds,z['obs'].attrs.get('_index','_index'))[positions].astype(str)
        digest=hashlib.sha256('\0'.join(ds+'::'+x for x in ids).encode()).hexdigest()
        if digest!=entry['cell_axis_sha256'] or digest!=meta['cell_axis_sha256']:raise ValueError('Native cell identity axis changed')
        native=(data['native_family_bitmask'][:,None]&(1<<np.arange(12)))>0;ref=(data['reference_family_bitmask'][:,None]&(1<<np.arange(12)))>0
        own=(data['own_fine_bitmask'][:,None]&(1<<np.arange(10)))>0;ref_fine=(data['reference_fine_bitmask'][:,None]&(1<<np.arange(10)))>0
        nn,rn=native.sum(1),ref.sum(1);native_major=decode('native_major');ref_major=decode('low_reference_major')
        nf=np.array([FAMILIES[j] for j in native.argmax(1)]);rf=np.array([FAMILIES[j] for j in ref.argmax(1)])
        hard=~positive_all[sl]|(nn>1)|(rn>1)|((nn==1)&(rn==1)&(nf!=rf));hard_count+=int(hard.sum())
        with gzip.open(annot/'outputs'/aindex[ds]['output'],'rt') as f:groups=json.load(f)['groups']
        oldfine=reader.frozen_labels(ds)['fine'][positions].astype(str);corrections={}
        for method in ['A','B','C']:
            major=decode(method+'__consensus_major');ctype=decode(method+'__consensus_celltype');coarse_ev=decode(method+'__coarse_evidence');fine_ev=decode(method+'__fine_evidence')
            group={g['partition_code']:g for g in groups if g['method']==method};gm=np.array([group[int(c)]['local_major'] for c in parts[method][sl]]);gt=np.array([group[int(c)]['local_celltype'] for c in parts[method][sl]])
            gf=np.array([family(x) for x in gm]);expected=np.full(n,'unresolved',dtype='U32')
            ok_native=(nn==1)&~hard;expected[ok_native]=native_major[ok_native]
            within=ok_native&(ref_major!='not_available')&np.array([not compatible(a,b) for a,b in zip(native_major,ref_major)])
            expected[within]=nf[within]
            drop=(nn==0)&(rn==1)&~hard&(gf==rf)
            expected[drop]=np.where(ref_major[drop]!='not_available',ref_major[drop],gm[drop])
            drop_conflict=drop&(ref_major!='not_available')&np.array([not compatible(a,b) for a,b in zip(ref_major,gm)])
            expected[drop_conflict]=rf[drop_conflict]
            np.testing.assert_array_equal(major,expected)
            eligible=np.column_stack([np.array([compatible(parent,x) for x in expected]) for parent in FINE_PARENT])
            direct=own&ref_fine&eligible;direct[:,1]&=~direct[:,0]
            direct_ok=(direct.sum(1)==1)&ok_native&~within
            wanted=np.array(FINE)[direct.argmax(1)]
            np.testing.assert_array_equal(ctype[direct_ok],wanted[direct_ok]);direct_count+=int(direct_ok.sum())
            observed_direct=np.isin(fine_ev,['own_native_and_unique_reference_supported','own_gamma_RNA_and_reference_candidate_requires_orthogonal_review'])
            np.testing.assert_array_equal(observed_direct,direct_ok)
            for j,t in enumerate(FINE):
                assigned=ctype==t
                if t=='T_gamma_delta_RNA_candidate' and (assigned&~own[:,j]).any():raise ValueError('Gamma cell lacks direct RNA evidence')
                indirect=assigned&~direct_ok
                if (indirect&((gt!=t)|~fit[method][sl])).any():raise ValueError('Indirect fine label lacks fitted group context')
                if (indirect&~(own[:,j]|ref_fine[:,j])).any():raise ValueError('Indirect fine label lacks cell evidence')
            np.testing.assert_array_equal(data[method+'__fine_discovery_fit_available'],fit[method][sl]);np.testing.assert_array_equal(data[method+'__parent_discovery_group_code'],parts[method][sl])
            actual=collections.Counter(zip(oldfine,ctype));stored={(r['old_fine'],r['guarded_celltype']):r['n_cells'] for r in meta['old_fine_crosswalk'] if r['method']==method}
            if dict(actual)!=stored:raise ValueError('Old-new cell crosswalk differs')
            for r in meta['parent_group_records']:
                if r['method']!=method:continue
                keep=parts[method][sl]==r['partition_code']
                if r['n_cells']!=int(keep.sum()) or dict(collections.Counter(ctype[keep]))!=r['guarded_celltype_counts']:raise ValueError('Parent-group guarded denominator differs')
                if r['marker_origin']!='parent_discovery_group_not_recomputed_for_guarded_subset':raise ValueError('Guarded subset borrowed marker origin missing')
            old_major=decode(method+'__consensus_major',True);old_type=decode(method+'__consensus_celltype',True)
            np.testing.assert_array_equal(major[~drop_conflict],old_major[~drop_conflict]);np.testing.assert_array_equal(ctype[~drop_conflict],old_type[~drop_conflict])
            if (ctype[drop_conflict]!=rf[drop_conflict]).any():raise ValueError('Dropout correction does not retain broad family')
            corrections[method]=int(drop_conflict.sum());corrected_count+=int(drop_conflict.sum())
        if ds==audit_source:
            key=reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key'];blocks=[]
            for lo in range(0,n,2048):blocks.append(direct_rows(z[key],reader.rows(ds)[positions[lo:lo+2048]]))
            raw=sp.vstack(blocks,format='csr');del blocks;cols=raw.tocsc();lookup=reader.feature_lookup(ds);cache={}
            def det(g):
                if g not in cache:cache[g]=(cols[:,next(iter(lookup[g]))].toarray().ravel()>0) if len(lookup.get(g,set()))==1 else np.zeros(n,bool)
                return cache[g]
            def at_least_two(panel):return np.sum([det(g) for g in panel],axis=0)>=2
            reconstructed=np.column_stack([
                det('FOXP3')&(det('IL2RA')|det('CTLA4')),
                det('CD4')&at_least_two(LINEAGE_PANELS['T']),
                det('CD8A')&det('CD8B')&(det('CD3D')|det('CD3E')),
                det('TRDC')&(det('TRGC1')|det('TRGC2'))&det('CD3D'),
                at_least_two(LINEAGE_PANELS['Neutrophil_related'])&(det('CSF3R')|det('FCGR3B')),
                at_least_two(TYPE_PANELS['pDC_related']),at_least_two(TYPE_PANELS['cDC1_related']),at_least_two(TYPE_PANELS['cDC2_related']),
                at_least_two(TYPE_PANELS['Macrophage_related']),at_least_two(TYPE_PANELS['Monocyte_related'])])
            np.testing.assert_array_equal(reconstructed,own);native_values+=int(own.size);del raw,cols,cache
        reports.append({'source':ds,'n_pilot_cells':n,'all_cell_coarse_and_direct_fine_guards_independently_checked':True,'all_crosswalk_and_parent_denominators_checked':True,'v2_dropout_conflict_corrections_by_method':corrections})
        offset+=n;save(root/'checkpoints'/(ds+'.json'),{'source':ds,'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()})
    if offset!=len(source) or len(reports)!=len(inv):raise ValueError('Independent confidence scope incomplete')
    payload={'source_reports':reports,'n_pilot_cells':offset,'direct_native_fine_boolean_values_replayed':native_values,'independent_direct_fine_decision_cases_across_methods':direct_count,
        'independent_hard_conflict_cells':hard_count,'v2_correction_cases_across_methods':corrected_count,'numeric_audit_source_selection':'smallest_frozen_pilot_source',
        'native_replay_source':audit_source,'independent_identity_accuracy_claim':False,'v1_preserved':True}
    out=root/'outputs/cell_identity_confidence_audit_v2.json';save(out,payload)
    save(root/'provenance/terminal_v1.json',{'status':'COMMON_CELL_CONFIDENCE_V2_INDEPENDENT_GUARDS_NATIVE_FINE_AND_CROSSWALK_ACCEPTED',
        'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'output_sha256':sha256_file(out),
        'n_pilot_cells':offset,'direct_native_fine_boolean_values_replayed':native_values,'all_source_scope_complete':True,
        'wall_seconds':time.monotonic()-started,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        'new_authoritative_labels':False,'native_counts_copied':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--producer',required=True);p.add_argument('--previous',required=True);a=p.parse_args();root=Path(a.output)
    try:main(root,Path(a.producer),Path(a.previous))
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)})
        raise

"""All existing source/frozen-pool native count inputs; no external scores.

Retain only measured common-native gene columns per source. Full physical RNA
quality and44 fixed-marker controls are computed from the original cell counts.
"""
import hashlib
import json
import math
import shutil
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.sparse import csr_matrix

from NSCLC_external_lineage_inputs_v2 import load,cell_arms,sha,matched_quota,VIEWS
from external_functional_coverage_inputs_v1 import read_inputs


def mapped(axis,names):
    cols=[next(iter(axis['lookup'][g])) if g in axis['lookup'] and len(axis['lookup'][g])==1 else -1 for g in names]
    counts=Counter(c for c in cols if c>=0)
    return np.array([c if c>=0 and counts[c]==1 else -1 for c in cols],np.int32)


def main():
    data=load();root=data['root'];lib=data['lib'];n=len(lib);nu=n*4
    controls=root/'NSCLC_native_MFA_external_role_marker_profiles_20261016_v2'
    marker_proof=controls/'private_independent_all_role_marker_audit.json'
    assert sha(marker_proof)=='b40688206ea0abfa9ffd0bb795e85e871d59c745db13646416a3831f32d3c1b4'
    marker_terminal=json.loads((controls/'private_terminal.json').read_text())
    for name in ['private_fixed_marker_pool_sums.npz','private_all6760_frozen_role_units.tsv']:
        assert sha(controls/name)==marker_terminal['files'][name]['sha256']
    with np.load(controls/'private_fixed_marker_pool_sums.npz',allow_pickle=False) as z:
        marker={k:z[k] for k in z.files}
    markers=marker['markers'].astype(str);expected=marker['n_cells'];assert len(markers)==44 and len(expected)==nu
    axes=read_inputs(root)['source_axes']
    native_path=root/'NSCLC_full_native_gene_units_20261016_v1/private_native_gene_axis.tsv'
    genes=pd.read_csv(native_path,sep='\t',keep_default_na=False).gene.to_numpy(str)
    assert len(genes)==31831 and len(set(genes))==len(genes)
    calibration=root/'NSCLC_native_MFA_external_RNA_universe_calibration_20261016_v1'
    plan0=json.loads((calibration/'plan_before_native_count_universe_calibration.json').read_text())
    maskfile=calibration/'private_native_universe_masks.tsv.gz';assert sha(maskfile)==plan0['source_mask_sha256']
    masks=pd.read_csv(maskfile,sep='\t',index_col=0);assert masks.index.astype(str).tolist()==genes.tolist()
    mapping=pd.read_csv(calibration/'private_all61_source_native_universe_map.tsv',sep='\t').set_index(['arm','dataset'])
    infos={};by_library=np.full(n,-1,int);bounds=0
    for si,((arm,source),axis) in enumerate(sorted(axes.items())):
        libraries=lib.index[lib.arm.eq(arm)&lib.dataset.eq(source)].to_numpy(int)
        by_library[libraries]=si
        columns=mapped(axis,genes);native=np.flatnonzero(columns>=0)
        row=mapping.loc[(arm,source)];assert np.array_equal(columns>=0,masks['mask_%03d'%row.mask_index].to_numpy(bool))
        assert len(native)==int(row.measured_native_universe_genes)
        global_units=(libraries[:,None]*4+np.arange(4)).reshape(-1)
        inverse=np.full(nu,-1,int);inverse[global_units]=np.arange(len(global_units))
        infos[si]=dict(arm=arm,dataset=source,libraries=libraries,global_units=global_units,
            inverse_units=inverse,native_columns=native,source_columns=columns[native],marker_source_columns=mapped(axis,markers))
        bounds+=len(global_units)*len(native)*8
    assert len(infos)==61 and (by_library>=0).all()
    cap=int(math.ceil(bounds*1.002)+64*1024**2)
    quota=matched_quota(root)
    assert min(quota['group_remaining_bytes'],shutil.disk_usage(root).free)>500_000_000_000+1073741824+cap
    memory=int(next(s.split()[1] for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith('MemAvailable:')))*1024
    assert memory>2*bounds+3*1024**3
    out=root/'NSCLC_native_MFA_external_frozen_role_count_inputs_20261016_v1'
    assert not out.exists();out.mkdir(mode=0o700);(out/'private_source_counts').mkdir(mode=0o700)
    # Audit selection uses metadata and fixed hashing only, before any new RNA.
    audit_units=pd.read_csv(controls/'private_all6760_frozen_role_units.tsv',sep='\t',keep_default_na=False)
    audit_units['priority']=[hashlib.sha256(('20261016|'+r.arm+'|'+r.dataset+'|'+r.library_id+'|'+getattr(r,'view')).encode()).hexdigest() for r in audit_units.itertuples()]
    chosen=set()
    for _,g in audit_units.groupby(['arm','dataset'],sort=True):
        positive=g[g.n_cells>0];chosen.add(int((positive if len(positive) else g).sort_values('priority').index[0]))
    for _,g in audit_units.groupby('view',sort=True):
        positive=g[g.n_cells>0]
        if len(positive):chosen.add(int(positive.sort_values('priority').index[0]))
    for take in [audit_units.n_cells.eq(0),audit_units.n_cells.between(1,19),audit_units.n_cells.ge(20)]:
        chosen.update(audit_units[take].sort_values('priority').head(4).index.tolist())
    for i in audit_units.sort_values('priority').index:
        if len(chosen)>=96:break
        chosen.add(int(i))
    assert len(chosen)==96
    audit_units.loc[sorted(chosen)].to_csv(out/'private_prespecified96_raw_count_audit_units.tsv.gz',sep='\t',index=False,compression='gzip')
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),input_metadata_sha256=data['hashes'],
        marker_audit_sha256=sha(marker_proof),native_gene_axis_sha256=sha(native_path),native_RNA_universe_mask_sha256=sha(maskfile),
        predetermined_raw_audit_units_sha256=sha(out/'private_prespecified96_raw_count_audit_units.tsv.gz'),
        scope='All61 existing source-arm/1690 libraries x4 frozen candidate pools. All cells in each frozen pool retained, including small/zero library pools. No author-channel substitution or outcome filtering.',
        measurement='Sum original integer cell counts over frozen source/library/role. Source archives contain ONLY uniquely mapped measured common native columns, never expression-zero imputation for unmeasured genes. All common native genes retained, not merely learned features.',
        quality='Full physical SOURCE measured RNA and detected genes, common-native RNA/detection, symbol-prefix mito/ribo and zero-RNA cell counts. Original diagnostic parent counts retained separately; not declared equivalent to original lung broad-parent fraction.',
        normalization='No TMM or external factor scores in this stage. Original full native quality is not assumed platform-equivalent. Subsequent transfer definition and cohort model require separate outcome-free contract.',
        controls='Fresh raw44-marker logCP10K/positive-cell counts over ALL selected pools must match independently accepted prior footprints. This includes new original raw read, not just unit-cache reaggregation.',
        independent_raw_audit='96 metadata-selected units covering every source, all roles, low and zero capture; independently re-extract all measured common genes and quality after production.',
        resources=dict(quota=quota,full_uncompressed_count_bound_bytes=bounds,output_cap_bytes=cap,
            cap_reason='Uncompressed count arrays plus0.2percent archive bound and64MiB auxiliary/metadata allowance. New count-input phase, distinct from previous256MiB score-only phase.',memory_available_bytes=memory),
        original_labels_or_CM_changed=False,external_programme_scores_or_P=False,new_download=False)
    (out/'plan_before_all_external_native_count_extraction.json').write_text(json.dumps(plan,indent=2)+'\n')
    for info in infos.values():info['counts']=np.zeros((len(info['global_units']),len(info['native_columns'])),np.int64)
    quality=np.zeros((nu,12),float);mlog=np.zeros((nu,len(markers)));mpositive=np.zeros((nu,len(markers)),np.int64)
    oldbase=root.parents[2]/'Project_v3/annotation_rebuild/20261010_uniform_expression_annotation_v1'
    for arm in cell_arms(data):
        original=arm['arm']=='frozen_original_annotation'
        z=zarr.open_group(str(data['old_zarr'] if original else data['atlas']),mode='r')
        raw=z['X'] if original else z['layers']['counts'];nr,nf=raw.attrs['shape']
        unit=np.full(nr,-1,np.int32);roles=arm['roles'][0];take=roles>=0
        indices=np.arange(len(roles))
        if original:
            with np.load(data['old']/'private_canonical_row_selection.npz',allow_pickle=False) as f:indices=f['global_rows']
        assert len(indices)==len(roles) and len(np.unique(indices))==len(indices)
        unit[indices[take]]=arm['library_index'][take]*4+roles[take]
        if not original:
            source_order=np.asarray(read_elem(z['uns']['source_order'])).astype(str).tolist()
            measured=np.asarray(read_elem(z['varm']['measured_in_source']),bool)
            symbol=np.asarray(read_elem(z['var']['gene_symbol'])).astype(str)
        for si,info in infos.items():
            if info['arm']!=arm['arm']:continue
            common_map=np.full(nf,-1,np.int32);marker_map=np.full(nf,-1,np.int16)
            if original:
                g=pd.read_csv(oldbase/'by_gse'/info['dataset']/'genes.tsv',sep='\t',keep_default_na=False)
                sourcecols=g.core_column.to_numpy(int);assert len(sourcecols)==len(set(sourcecols))
                columns=sourcecols[info['source_columns']]
                present=np.zeros(nf,bool);present[sourcecols]=True
                mt=np.zeros(nf,bool);rb=np.zeros(nf,bool)
                mt[sourcecols]=g.gene.str.startswith('MT-').to_numpy()
                rb[sourcecols]=g.gene.str.match(r'^(RPL|RPS)').to_numpy()
                marker_columns=np.array([sourcecols[c] if c>=0 else -1 for c in info['marker_source_columns']],int)
            else:
                columns=info['source_columns'];present=measured[:,source_order.index(info['dataset'])]
                mt=pd.Series(symbol).str.startswith('MT-').to_numpy();rb=pd.Series(symbol).str.match(r'^(RPL|RPS)').to_numpy()
                marker_columns=info['marker_source_columns']
            common_map[columns]=np.arange(len(columns));valid=marker_columns>=0;marker_map[marker_columns[valid]]=np.flatnonzero(valid)
            info.update(common_map=common_map,marker_map=marker_map,present=present,mito=mt,ribo=rb)
        for start in range(0,nr,4096):
            stop=min(nr,start+4096);valid=unit[start:stop]>=0
            if valid.any():
                ptr=np.asarray(raw['indptr'][start:stop+1],np.int64);lo,hi=int(ptr[0]),int(ptr[-1])
                x=csr_matrix((np.asarray(raw['data'][lo:hi]),np.asarray(raw['indices'][lo:hi]),ptr-lo),shape=(stop-start,nf))[valid].copy()
                x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
                assert np.isfinite(x.data).all() and (x.data>0).all() and np.equal(x.data,np.floor(x.data)).all()
                ids=unit[start:stop][valid];source=by_library[ids//4]
                for si in np.unique(source):
                    mask=source==si;matrix=x[mask].copy();units=ids[mask];info=infos[int(si)];co=matrix.tocoo()
                    assert info['present'][co.col].all(),(info['dataset'],'Nonzero count outside full physical source measurement mask',int((~info['present'][co.col]).sum()))
                    total=np.asarray(matrix.sum(axis=1)).ravel();detected=np.diff(matrix.indptr)
                    common=info['common_map'][co.col];keep=common>=0
                    cc=csr_matrix((co.data[keep].astype(np.int64),(co.row[keep],common[keep])),shape=(len(units),len(info['native_columns'])))
                    cc.sum_duplicates();cc.eliminate_zeros();ct=np.asarray(cc.sum(axis=1)).ravel();cd=np.diff(cc.indptr)
                    local=info['inverse_units'][units];assert (local>=0).all()
                    group=csr_matrix((np.ones(len(units),np.int64),(local,np.arange(len(units)))),shape=(len(info['global_units']),len(units)))
                    pooled=(group@cc).tocoo();np.add.at(info['counts'],(pooled.row,pooled.col),pooled.data)
                    mito=np.bincount(co.row,weights=co.data*info['mito'][co.col],minlength=len(units))
                    ribo=np.bincount(co.row,weights=co.data*info['ribo'][co.col],minlength=len(units))
                    values=[np.ones(len(units)),total,ct,np.log1p(total),np.log1p(detected),np.log1p(ct),np.log1p(cd),mito,ribo,total==0,detected,cd]
                    for j,value in enumerate(values):quality[:,j]+=np.bincount(units,weights=value,minlength=nu)
                    target=info['marker_map'][co.col];keep=target>=0
                    rows=co.row[keep];positions=target[keep];scale=np.divide(10000.,total,out=np.zeros_like(total,dtype=float),where=total>0)
                    value=np.log1p(co.data[keep].astype(float)*scale[rows])
                    np.add.at(mlog,(units[rows],positions),value);np.add.at(mpositive,(units[rows],positions),1)
            if start//4096%32==0 or stop==nr:
                cp=dict(observed_UTC=datetime.now(timezone.utc).isoformat(),stage='native_count_stream',arm=arm['arm'],rows_read=stop,total_rows=nr,selected_cells_accumulated=int(quality[:,0].sum()))
                (out/'private_checkpoint.json').write_text(json.dumps(cp)+'\n')
        for info in infos.values():
            if info['arm']==arm['arm']:
                for field in ['common_map','marker_map','present','mito','ribo']:info.pop(field)
    assert np.array_equal(quality[:,0],expected)
    assert np.array_equal(mpositive,marker['positive_counts'])
    marker_error=float(np.max(abs(mlog-marker['sum_logCP10K'])))
    assert np.allclose(mlog,marker['sum_logCP10K'],atol=1e-6,rtol=1e-9),marker_error
    units=pd.read_csv(controls/'private_all6760_frozen_role_units.tsv',sep='\t',keep_default_na=False)
    with np.load(root/'NSCLC_native_MFA_external_lineage_qualification_20261016_v2/private_all1690_three_channel_role_counts.npz',allow_pickle=False) as f:parents=f['parent_counts'].reshape(-1)
    units['full_source_RNA_sum']=quality[:,1].astype(np.int64);units['common_native_RNA_sum']=quality[:,2].astype(np.int64)
    units['zero_full_RNA_cells']=quality[:,9].astype(np.int64);units['diagnostic_parent_cells']=parents
    for name,column in [('mean_log_full_source_RNA',3),('mean_log_full_source_detected_genes',4),('mean_log_common_native_RNA',5),('mean_log_common_native_detected_genes',6),('mean_full_source_detected_genes',10),('mean_common_native_detected_genes',11)]:
        units[name]=np.divide(quality[:,column],quality[:,0],out=np.full(nu,np.nan),where=quality[:,0]>0)
    for name,column in [('symbol_pct_mito_RNA',7),('symbol_pct_ribo_RNA',8)]:units[name]=np.divide(100*quality[:,column],quality[:,1],out=np.full(nu,np.nan),where=quality[:,1]>0)
    units['diagnostic_fraction_within_parent']=np.divide(quality[:,0],parents,out=np.full(nu,np.nan),where=parents>0)
    archives=[]
    for si,info in infos.items():
        rows=info['global_units'];assert (info['counts']>=0).all()
        assert np.array_equal(info['counts'].sum(axis=1),quality[rows,2])
        path=out/'private_source_counts'/('source_%03d.npz'%si)
        np.savez_compressed(path,counts=info['counts'],common_native_columns=info['native_columns'],global_unit_indices=rows)
        archives.append(dict(source_index=si,arm=info['arm'],dataset=info['dataset'],units=len(rows),common_native_genes=len(info['native_columns']),archive=str(path.relative_to(out)),sha256=sha(path),bytes=path.stat().st_size))
    units.to_csv(out/'private_all6760_native_count_units_and_quality.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(archives).to_csv(out/'private_all61_count_archive_manifest.tsv',sep='\t',index=False)
    pd.DataFrame(dict(native_column=np.arange(len(genes)),gene=genes)).to_csv(out/'private_complete_native_gene_axis.tsv',sep='\t',index=False)
    np.savez_compressed(out/'private_all44_fresh_raw_marker_controls.npz',n_cells=quality[:,0].astype(np.int64),sum_logCP10K=mlog,positive_counts=mpositive,markers=markers)
    terminal=dict(status='COMPLETE_ALL61_SOURCE6760_FROZEN_ROLE_NATIVE_COUNT_INPUTS_WITH_FRESH44_MARKER_CONTROL_AUDIT_PENDING',
        finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),fresh_marker_maximum_sum_log_error=marker_error,
        files={str(p.relative_to(out)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.rglob('*') if p.is_file()},
        original_annotations_or_CM_changed=False,external_programme_scores_or_P=False)
    assert sum(v['bytes'] for v in terminal['files'].values())<cap
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

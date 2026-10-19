"""All target summaries rebuilt with independent unique-value ranks and LS."""
import hashlib,json,os,shutil,tarfile,gzip
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import rankdata


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def integer_tokens(v):
    a=pd.to_numeric(pd.Series(v).astype(str),errors='raise').to_numpy()
    assert np.isfinite(a).all() and np.array_equal(a,np.floor(a))
    return a.astype(np.int64).astype(str)


def independent_ranks(x,cap):
    out=np.empty(x.shape,float)
    for i,row in enumerate(x):
        _,inverse,counts=np.unique(row,return_inverse=True,return_counts=True)
        average=1+len(row)-np.cumsum(counts)+(counts-1)/2
        out[i]=np.minimum(cap,average[inverse]);out[i,row<=0]=cap
    return out


def design(total,detected,area,labels):
    n=len(total)
    if n<30:return None,'FEWER_THAN30_AUTHOR_LINEAGE_CELLS'
    x=np.stack([rankdata(v,method='average') for v in [total,detected,area]],axis=1);s=x.std(0)
    if (s<=1e-12).any() or not np.isfinite(x).all():return None,'CONSTANT_REQUIRED_PANEL_QUALITY'
    x=(x-x.mean(0))/s;dummy=pd.get_dummies(pd.Series(labels),drop_first=True).to_numpy(float)
    X=np.column_stack([np.ones(n),x,dummy]);singular=np.linalg.svd(X,compute_uv=False)
    rank=int((singular>singular[0]*max(X.shape)*np.finfo(float).eps).sum())
    if rank<X.shape[1] or singular[0]/singular[-1]>=1e8:return None,'REQUIRED_PANEL_QUALITY_OR_AUTHOR_SUBTYPE_DESIGN_UNIDENTIFIABLE'
    if n-rank<20:return None,'FEWER_THAN20_RESIDUAL_DF'
    return X,'ESTIMABLE_CAPTURE_DESCRIPTION_ONLY'


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
    parent=Path('/dev/shm/pancancer-native-state-capture-20261020-v2/private_output')
    assert out.parent.parent==Path('/dev/shm') and not out.exists() and shutil.disk_usage('/dev/shm').free>512*1024**2
    assert sha(parent/'private_terminal.json')=='efa79e09ea535391339ec81d4ecdbb0088500cea23d9ef003814d8c020ea2a1d'
    terminal=json.loads((parent/'private_terminal.json').read_text())
    for name,d in terminal['files'].items():assert sha(parent/name)==d['sha256']
    cp=Path('/dev/shm/pancancer-native-state-probe-match-20261020-v1/private_output/private_frozen_matched_measurement_and_target_contract.json')
    c=json.loads(cp.read_text());f=pd.read_csv(parent/'private_all_field_proxy_capture_and_quality.tsv.gz',sep='\t',dtype={'sample_id':str,'FOV':str})
    labels_table=pd.read_csv(parent/'private_all_field_author_subtype_proxy_support.tsv.gz',sep='\t',dtype={'sample_id':str,'FOV':str})
    expected=f.set_index(['sample_id','FOV','major','proxy_group_id']);subtype=labels_table.set_index(['sample_id','FOV','major','author_label','proxy_group_id'])
    assert expected.index.is_unique and subtype.index.is_unique
    proofs={x['sample_id']:x for x in json.loads((parent/'private_all_native_raw_parity_proofs.json').read_text())['proofs']}
    old=json.loads((root/'immune_same_name_family_source_review_v1/spatial_contact_native_preflight_v3/inventory.json').read_text());original={x['sample_id']:x for x in old['inputs']}
    checked=0;subchecked=0;cells=0;archives={};seen=set();subseen=set();max_error=0.;ranks_count=0
    out.mkdir(mode=0o700)
    for ordinal,t in enumerate(c['targets']):
        assert sha(t['analysis_path'])==t['input_sha256'] and sha(t['cache_path'])==t['cache_sha256']
        a=ad.read_h5ad(t['analysis_path']);cells+=a.n_obs;assert a.n_obs==t['native_cells']
        with np.load(t['cache_path'],allow_pickle=False) as z:data={k:z[k] for k in ['author_labels','qc','fov','cell_id','coords','area','native_total','native_ngenes']}
        np.testing.assert_array_equal(np.asarray(a.X.sum(1)).ravel(),data['native_total'])
        np.testing.assert_array_equal(np.asarray((a.X>0).sum(1)).ravel(),data['native_ngenes'])
        proof=proofs[t['sample_id']];assert proof['all_native_RNA_positions_replayed']==a.n_obs*a.n_vars
        lookup={k:j for j,k in enumerate(zip(data['fov'].astype(str),data['cell_id'].astype(str)))};assert len(lookup)==a.n_obs
        archive=Path(original[t['sample_id']]['archive'])
        if str(archive) not in archives:
            archives[str(archive)]=sha(archive);assert archives[str(archive)]==original[t['sample_id']]['archive_sha256']
        mapping=dict(a.uns['archive_member_map'])
        with tarfile.open(archive,'r:*') as tar:
            stream=tar.extractfile(mapping['metadata']);stream=gzip.GzipFile(fileobj=stream) if mapping['metadata'].endswith('.gz') else stream
            meta=pd.read_csv(stream,usecols=['fov','cell_ID','Area','CenterX_global_px','CenterY_global_px'])
        ff=integer_tokens(meta.fov);cc=integer_tokens(meta.cell_ID);keep=cc!='0';meta=meta.loc[keep];ff=ff[keep];cc=cc[keep]
        keys=list(zip(ff,cc));assert len(keys)==len(set(keys)) and set(lookup)<=set(keys)
        select=np.array([k in lookup for k in keys]);positions=np.array([lookup[k] for k,take in zip(keys,select) if take])
        assert len(positions)==a.n_obs
        np.testing.assert_allclose(meta.loc[select,['CenterX_global_px','CenterY_global_px']].to_numpy(float)*.12028,data['coords'][positions],atol=1e-6,rtol=0)
        np.testing.assert_array_equal(meta.loc[select,'Area'].to_numpy(float),data['area'][positions])
        extra={field:int(((ff==field)&~select).sum()) for field in sorted(set(ff))}
        assert proof['metadata_only_counts_per_FOV']==extra and proof['raw_metadata_keys_outside_frozen_expression']==int((~select).sum())
        labels=data['author_labels'].astype(str);qc=data['qc'].astype(bool);fov=data['fov'].astype(str)
        target_rows=np.flatnonzero(qc&np.isin(labels,c['T_author_labels']+c['strict_macrophage_author_labels']))
        scores={};support={};bybg=defaultdict(list)
        for m in t['proxy_mappings']:bybg[tuple(m['background_indices'])].append(m)
        for bg,ms in bybg.items():
            cols=list(bg);local={j:k for k,j in enumerate(bg)}
            for m in ms:scores[m['proxy_group_id']]=np.full(a.n_obs,np.nan);support[m['proxy_group_id']]=np.zeros(a.n_obs,int)
            for start in range(0,len(target_rows),1024):
                ix=target_rows[start:start+1024];x=a.X[ix][:,cols].toarray();ranks=independent_ranks(x,len(cols));ranks_count+=ranks.size
                for m in ms:
                    sig=[local[j] for j in m['signature_indices']];n=len(sig);base=n*(n+1)/2;cap=m['maxrank'];assert cap==len(cols)
                    scores[m['proxy_group_id']][ix]=1-(ranks[:,sig].sum(1)-base)/(n*cap-base)
                    support[m['proxy_group_id']][ix]=(x[:,sig]>0).sum(1)
        for field in sorted(set(fov)):
            for major,author in [('T cell',c['T_author_labels']),('Macrophage',c['strict_macrophage_author_labels'])]:
                ix=np.flatnonzero(qc&(fov==field)&np.isin(labels,author));nn=len(ix)
                X,status=design(data['native_total'][ix],data['native_ngenes'][ix],data['area'][ix],labels[ix])
                for m in t['proxy_mappings']:
                    if m['major']!=major:continue
                    gid=m['proxy_group_id'];key=(t['sample_id'],field,major,gid);assert key not in seen;seen.add(key);r=expected.loc[key]
                    y=scores[gid][ix];multi=int((support[gid][ix]>=3).sum());var=float(np.var(y)) if nn else np.nan;resvar=np.nan;st=status
                    if X is not None:
                        residual=y-X@np.linalg.lstsq(X,y,rcond=None)[0];resvar=float(np.var(residual))
                        if var<=1e-12 or resvar<=var*1e-10:st='CONSTANT_OR_QUALITY_AUTHOR_EXPLAINED_PROXY'
                    if st=='ESTIMABLE_CAPTURE_DESCRIPTION_ONLY' and multi<10:st='FEWER_THAN10_CELLS_WITH3_SIGNATURE_PROBES'
                    assert r.author_lineage_cells==nn and r.signature_physical_probes==len(m['signature_indices'])
                    assert r.cells_with1_signature_probe==int((support[gid][ix]>=1).sum()) and r.cells_with3_signature_probes==multi
                    assert r.cells_with_atleast_half_signature==int((support[gid][ix]>=int(np.ceil(len(m['signature_indices'])/2))).sum())
                    assert r.score_positive_cells==int((y>0).sum()) and r.metadata_only_source_rows_in_FOV==extra.get(field,0)
                    assert st==r.status and (st=='ESTIMABLE_CAPTURE_DESCRIPTION_ONLY')==r.native_capture_qualified
                    values=[var,float(np.quantile(y,.75)-np.quantile(y,.25)) if nn else np.nan,resvar,resvar/var if var>1e-12 else np.nan]
                    actual=[r.score_variance,r.score_IQR,r.after_panel_quality_and_author_subtype_residual_variance,r.residual_variance_fraction]
                    np.testing.assert_allclose(values,actual,atol=2e-10,rtol=0,equal_nan=True)
                    err=np.abs(np.array(values)-np.array(actual));max_error=max(max_error,float(np.nanmax(err)) if np.isfinite(err).any() else 0);checked+=1
                    for label in sorted(set(labels[ix])):
                        at=ix[labels[ix]==label];v=scores[gid][at];sk=(*key[:3],label,gid);assert sk not in subseen;subseen.add(sk);sr=subtype.loc[sk]
                        assert sr.cells==len(at) and sr.cells_with3_signature_probes==int((support[gid][at]>=3).sum())
                        np.testing.assert_allclose(sr.score_variance,float(np.var(v)),atol=1e-14,rtol=0);subchecked+=1
        (out/'private_progress.json').write_text(json.dumps(dict(stage='INDEPENDENT_UNIQUE_VALUE_RANK_AND_ALL_FIELD_CHECKS',images_completed=ordinal+1,fields_checked=checked,subtype_positions_checked=subchecked,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
        del a,data,scores,support
    assert checked==len(f)==2426 and subchecked==len(labels_table)==4846 and cells==190038
    assert seen==set(expected.index) and subseen==set(subtype.index)
    s=pd.read_csv(parent/'private_complete_source_proxy_capture_summary.tsv',sep='\t').set_index(['dataset_id','cancer','major','proxy_group_id'])
    for key,g in f.groupby(['dataset_id','cancer','major','proxy_group_id'],sort=True):
        row=s.loc[key];assert row.target_image_positions==g.sample_id.nunique() and row.technical_fields==len(g)
        assert row.native_capture_qualified_fields==g.native_capture_qualified.sum() and row.total_author_cell_opportunities==g.author_lineage_cells.sum()
        assert row.cells_with3_signature_probes==g.cells_with3_signature_probes.sum()
        value=g.residual_variance_fraction.median();np.testing.assert_allclose(value,row.median_residual_variance_fraction,atol=2e-12,rtol=0,equal_nan=True)
    report=dict(status='ACCEPTED_ALL_NATIVE_TARGET_CAPTURE_FIELDS_SUBTYPES_SOURCE_SUMMARIES_AND_METADATA_DOMAIN',completed_UTC=datetime.now(timezone.utc).isoformat(),
                auditor_sha256=sha(__file__),producer_terminal_sha256=sha(parent/'private_terminal.json'),
                native_input_images=4,native_cells=190038,all_field_proxy_positions_checked=checked,all_author_subtype_proxy_positions_checked=subchecked,
                all_source_proxy_summary_positions_checked=len(s),independent_unique_value_rank_positions=ranks_count,max_field_numeric_error=max_error,
                actual_all_native_panel_totals_and_detected_checked=True,all_raw_metadata_unique_complete_joins_and_extra_counts_rebuilt=True,
                raw_expression_archive_rows_replayed_again=False,all_output_summary_fields_rebuilt_from_native_X=True,
                no_score_expression_graph_arrays_written=True,new_spatial_or_clinical_effects=False,scientific_outputs_public=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()

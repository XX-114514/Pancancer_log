"""Matched positive-signature UCell capture; arrays RAM only, no spatial effects."""
import hashlib,importlib.util,json,os,shutil,time
from datetime import datetime,timezone
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from anndata._io.specs import read_elem


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def ucell(x,signature,maxrank):
    """Accepted positive-signature formula, zero genes at maxrank, average ties."""
    assert x.ndim==2 and len(signature)>0 and maxrank>len(signature)
    ranks=np.minimum(maxrank,rankdata(-x,axis=1,method='average'))
    ranks[x<=0]=maxrank
    k=len(signature);base=k*(k+1)/2
    return 1-(ranks[:,signature].sum(1)-base)/(k*maxrank-base)


def direct_score(row,signature,cap):
    rr=[]
    for j in signature:
        v=row[j]
        rr.append(cap if v<=0 else min(cap,1+np.count_nonzero(row>v)+(np.count_nonzero(row==v)-1)/2))
    n=len(signature);base=n*(n+1)/2
    return 1-(sum(rr)-base)/(n*cap-base)


def quality_design(total,detected,area,labels):
    n=len(total);x=rankdata(np.column_stack([total,detected,area]),axis=0,method='average')
    sd=x.std(0)
    if not np.isfinite(x).all() or (sd<=1e-12).any():return None,'CONSTANT_REQUIRED_PANEL_QUALITY'
    x=(x-x.mean(0))/sd;levels=np.unique(labels)
    dummy=np.column_stack([(labels==k).astype(float) for k in levels[1:]]) if len(levels)>1 else np.empty((n,0))
    X=np.column_stack([np.ones(n),x,dummy]);U,s,_=np.linalg.svd(X,full_matrices=False)
    rank=int(np.sum(s>s[0]*max(X.shape)*np.finfo(float).eps));cond=float(s[0]/s[-1]) if s[-1]>0 else np.inf
    if rank<X.shape[1] or cond>=1e8:return None,'REQUIRED_PANEL_QUALITY_OR_AUTHOR_SUBTYPE_DESIGN_UNIDENTIFIABLE'
    if n-rank<20:return None,'FEWER_THAN20_RESIDUAL_DF'
    return U,'ESTIMABLE_CAPTURE_DESCRIPTION_ONLY'


def ref_check(group,V,panels,native,score_manifest,probe_module):
    gse=sorted(group['reference_GSEs'])[0];info=native[(gse,group['major'])]
    rel='private_cell_scores/'+Path(info['relative_file']).stem+'_'+group['axis_id'][:10]+'.npz'
    p=V/'state_spatial_partial_panel_calibration_20261016_v2'/rel;assert sha(p)==score_manifest[rel]['sha256']
    with np.load(p,allow_pickle=False) as z:
        j=z['programme_ids'].tolist().index(group['parent_programme_ids'][0]);expected=z['proxy'][:,j].astype(float)
        selected=np.unique([0,len(expected)-1,int(np.argmax(expected))]);expected=expected[selected]
    panel=panels[group['axis_id']];p=V/'state_spatial_proxy_native_cells_20261016_v1'/info['relative_file']
    with h5py.File(p,'r') as h:
        k=h['var'].attrs['_index'];k=k.decode() if isinstance(k,bytes) else k;genes=np.asarray(read_elem(h['var'][k]),str)
        names=probe_module.resolved(panel,set(genes));physical=probe_module.physical_keys(panel)
        key_to_gene={physical[panel['symbol_to_feature_index'][g]]:g for g in names};lookup={g:j for j,g in enumerate(genes)}
        cols=[lookup[key_to_gene[k]] for k in group['background_probe_keys']]
        X=h['X'];encoding=X.attrs['encoding-type'];encoding=encoding.decode() if isinstance(encoding,bytes) else encoding;assert encoding=='csr_matrix'
        values=np.zeros((len(selected),len(genes)))
        for i,row in enumerate(selected):
            lo,hi=int(X['indptr'][row]),int(X['indptr'][row+1]);values[i,X['indices'][lo:hi]]=X['data'][lo:hi]
    sig=[group['background_probe_keys'].index(k) for k in group['signature_probe_keys']]
    actual=ucell(values[:,cols],sig,group['maxrank'])
    np.testing.assert_allclose(actual,expected,atol=1e-6,rtol=0)
    return dict(source_GSE=gse,major=group['major'],proxy_group_id=group['proxy_group_id'],actual_cached_reference_cells=len(selected),
                maximum_official_cached_partial_UCell_error=float(np.max(np.abs(actual-expected))))


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
    control=out.parent;planp=Path(os.environ['PANCANCER_NATIVE_CAPTURE_PLAN']);plan=json.loads(planp.read_text())
    assert out.parent.parent==Path('/dev/shm') and not out.exists() and shutil.disk_usage('/dev/shm').free>512*1024**2
    assert sha(__file__)==plan['producer_sha256']
    parent=Path('/dev/shm/pancancer-native-state-probe-match-20261020-v1/private_output')
    audit=Path('/dev/shm/pancancer-native-state-probe-match-audit-20261020-v1/private_output/private_terminal.json')
    for _ in range(240):
        if audit.is_file():break
        time.sleep(5)
    else:raise TimeoutError('Same profile audit dependency; no producer restart')
    a=json.loads(audit.read_text());assert a['status'].startswith('ACCEPTED_ALL_PROBE_PROFILES') and a['producer_terminal_sha256']==sha(parent/'private_terminal.json')
    assert a['auditor_sha256']==plan['profile_auditor_sha256']
    assert sha(audit)==plan['profile_audit_terminal_sha256']
    t=json.loads((parent/'private_terminal.json').read_text())
    name='private_frozen_matched_measurement_and_target_contract.json';assert sha(parent/name)==t['files'][name]['sha256']
    c=json.loads((parent/name).read_text());assert len(c['targets'])==4
    helpers={}
    for key in ['author_spatial_native_reader_v3.py','qualify_native_spatial_state_probe_contract_v1.py']:
        p=control/key;assert sha(p)==plan['helper_source_sha256'][key]
        spec=importlib.util.spec_from_file_location(key[:-3],p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);helpers[key]=m
    reader=helpers['author_spatial_native_reader_v3.py'];probe=helpers['qualify_native_spatial_state_probe_contract_v1.py']
    V=Path(json.loads((root/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    panels0=json.loads((V/'state_spatial_proxy_calibration_inputs_20261016_v1/fixed_all14_panel_gene_axes.json').read_text());panels={};bysample={p['sample_id']:p for p in panels0}
    for p in panels0:panels.setdefault(p['gene_axis_sha256'],p)
    native={(x['gse_id'],x['major']):x for x in json.loads((V/'state_spatial_proxy_native_cells_20261016_v1/inventory.json').read_text())['private_native_files']}
    sm=json.loads((V/'state_spatial_partial_panel_calibration_20261016_v2/inventory.json').read_text());score_manifest={x['relative_file']:x for x in sm['private_outputs']}
    old=root/'immune_same_name_family_source_review_v1/spatial_contact_native_preflight_v3';old_inventory=json.loads((old/'inventory.json').read_text());original={x['sample_id']:x for x in old_inventory['inputs']}
    out.mkdir(mode=0o700);save=lambda n,x:(out/n).write_text(json.dumps(x,indent=2)+'\n')
    save('plan_before_native_RNA_capture.json',plan)
    used_groups={m['proxy_group_id'] for target in c['targets'] for m in target['proxy_mappings']}
    groupmap={g['proxy_group_id']:g for g in c['groups']};reference_checks=[]
    for gid in sorted(used_groups):reference_checks.append(ref_check(groupmap[gid],V,panels,native,score_manifest,probe))
    save('private_actual_cached_official_reference_score_checks.json',reference_checks)
    fields=[];labels_rows=[];raw_proofs=[];archive_checked={};numeric_checks=0;max_numeric_error=0.;all_cells=0
    for ordinal,target in enumerate(c['targets']):
        p=Path(target['analysis_path']);assert sha(p)==target['input_sha256'] and sha(target['cache_path'])==target['cache_sha256']
        archive=Path(original[target['sample_id']]['archive'])
        if str(archive) not in archive_checked:
            save('private_progress.json',dict(stage='CURRENT_SOURCE_ARCHIVE_SHA',completed_images=ordinal,observed_UTC=datetime.now(timezone.utc).isoformat()))
            archive_checked[str(archive)]=sha(archive);assert archive_checked[str(archive)]==original[target['sample_id']]['archive_sha256']
        panel=bysample[target['sample_id']];union=sorted({j for m in target['proxy_mappings'] for j in m['background_indices']})
        aliases={j:[] for j in union}
        for gene,j in panel['symbol_to_feature_index'].items():
            if j in aliases:aliases[j].append(gene)
        genes=[]
        for j in union:
            pref=panel['preferred_symbol_by_feature'][j];genes.append(pref if pref in aliases[j] else sorted(aliases[j])[0])
        assert len(genes)==len(set(genes)) and len(union)>0
        data,proof=reader.load_native(p,Path(target['cache_path']),archive,genes)
        raw_proofs.append(dict(dataset_id=target['dataset_id'],cancer=target['cancer'],sample_id=target['sample_id'],**proof))
        all_cells+=len(data['qc']);lookup={j:i for i,j in enumerate(union)}
        labels=data['author_labels'].astype(str);qc=data['qc'].astype(bool);fov=data['fov'].astype(str)
        total=data['native_total'].astype(float);det=data['native_ngenes'].astype(float);area=data['area'].astype(float)
        count=data['counts'];n=len(count);assert n==target['native_cells']
        background_groups={}
        for m in target['proxy_mappings']:background_groups.setdefault(tuple(m['background_indices']),[]).append(m)
        scores={};support={}
        for bg,mappings in background_groups.items():
            cols=[lookup[j] for j in bg];positions={j:i for i,j in enumerate(bg)}
            for m in mappings:scores[m['proxy_group_id']]=np.full(n,np.nan);support[m['proxy_group_id']]=np.zeros(n,np.uint8)
            chosen=qc&np.isin(labels,c['T_author_labels']+c['strict_macrophage_author_labels'])
            rows=np.flatnonzero(chosen)
            for start in range(0,len(rows),1024):
                ix=rows[start:start+1024];x=count[ix][:,cols].astype(float)
                ranks=rankdata(-x,axis=1,method='average');ranks[x<=0]=len(cols)
                for m in mappings:
                    sig=[positions[j] for j in m['signature_indices']];k=len(sig);cap=m['maxrank'];base=k*(k+1)/2
                    assert cap==min(1500,len(cols)) and cap>k
                    scores[m['proxy_group_id']][ix]=1-(np.minimum(cap,ranks[:,sig]).sum(1)-base)/(k*cap-base)
                    support[m['proxy_group_id']][ix]=(x[:,sig]>0).sum(1)
            for m in mappings:
                name=m['proxy_group_id'];own=np.flatnonzero(qc&np.isin(labels,c['T_author_labels'] if m['major']=='T cell' else c['strict_macrophage_author_labels']))
                sig=[positions[j] for j in m['signature_indices']]
                check=np.unique(np.r_[own[:2],own[np.argmax(scores[name][own])]]) if len(own) else np.array([],int)
                for row in check:
                    expected=direct_score(count[row,cols].astype(float),sig,m['maxrank']);err=abs(expected-scores[name][row]);assert err<1e-12;max_numeric_error=max(max_numeric_error,err);numeric_checks+=1
        for field in sorted(set(fov)):
            for major,author_set in [('T cell',c['T_author_labels']),('Macrophage',c['strict_macrophage_author_labels'])]:
                ix=np.flatnonzero((fov==field)&qc&np.isin(labels,author_set));nn=len(ix)
                U,status=(quality_design(total[ix],det[ix],area[ix],labels[ix]) if nn>=30 else (None,'FEWER_THAN30_AUTHOR_LINEAGE_CELLS'))
                for m in target['proxy_mappings']:
                    if m['major']!=major:continue
                    gid=m['proxy_group_id'];y=scores[gid][ix];k=len(m['signature_indices']);multi=int((support[gid][ix]>=3).sum());variance=float(np.var(y)) if nn else np.nan
                    residual_variance=np.nan;st=status
                    if U is not None:
                        residual=y-U@(U.T@y);residual_variance=float(np.var(residual))
                        if variance<=1e-12 or residual_variance<=variance*1e-10:st='CONSTANT_OR_QUALITY_AUTHOR_EXPLAINED_PROXY'
                    if st=='ESTIMABLE_CAPTURE_DESCRIPTION_ONLY' and multi<10:st='FEWER_THAN10_CELLS_WITH3_SIGNATURE_PROBES'
                    fields.append(dict(dataset_id=target['dataset_id'],cancer=target['cancer'],sample_id=target['sample_id'],FOV=field,
                        major=major,proxy_group_id=gid,author_lineage_cells=nn,signature_physical_probes=k,
                        metadata_only_source_rows_in_FOV=int(data['metadata_only_counts_per_fov'].get(field,0)),
                        cells_with1_signature_probe=int((support[gid][ix]>=1).sum()),cells_with3_signature_probes=multi,
                        cells_with_atleast_half_signature=int((support[gid][ix]>=int(np.ceil(k/2))).sum()),
                        score_positive_cells=int((y>0).sum()),score_variance=variance,score_IQR=float(np.quantile(y,.75)-np.quantile(y,.25)) if nn else np.nan,
                        after_panel_quality_and_author_subtype_residual_variance=residual_variance,
                        residual_variance_fraction=residual_variance/variance if variance>1e-12 else np.nan,status=st,
                        native_capture_qualified=st=='ESTIMABLE_CAPTURE_DESCRIPTION_ONLY',reviewed_independent_patients=0))
                    for label in sorted(set(labels[ix])):
                        at=ix[labels[ix]==label];v=scores[gid][at]
                        labels_rows.append(dict(dataset_id=target['dataset_id'],cancer=target['cancer'],sample_id=target['sample_id'],FOV=field,
                           major=major,author_label=label,proxy_group_id=gid,cells=len(at),cells_with3_signature_probes=int((support[gid][at]>=3).sum()),
                           score_variance=float(np.var(v)),author_subtype_not_new_formal_state=True))
        save('private_progress.json',dict(stage='NATIVE_ALL_PANEL_PARITY_AND_CAPTURE_COMPLETE',completed_images=ordinal+1,total_images=4,
            native_cells_completed=all_cells,field_proxy_rows=len(fields),observed_UTC=datetime.now(timezone.utc).isoformat()))
    assert all_cells==plan['expected_native_cells']
    f=pd.DataFrame(fields);f.to_csv(out/'private_all_field_proxy_capture_and_quality.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(labels_rows).to_csv(out/'private_all_field_author_subtype_proxy_support.tsv.gz',sep='\t',index=False,compression='gzip')
    save('private_all_native_raw_parity_proofs.json',dict(input_images=4,raw_source_archives_checked=archive_checked,proofs=raw_proofs,
          no_RNA_score_graph_arrays_written=True))
    summary=f.groupby(['dataset_id','cancer','major','proxy_group_id'],sort=True).agg(target_image_positions=('sample_id','nunique'),
             technical_fields=('FOV','size'),native_capture_qualified_fields=('native_capture_qualified','sum'),
             total_author_cell_opportunities=('author_lineage_cells','sum'),cells_with3_signature_probes=('cells_with3_signature_probes','sum'),
             median_residual_variance_fraction=('residual_variance_fraction','median')).reset_index()
    summary.to_csv(out/'private_complete_source_proxy_capture_summary.tsv',sep='\t',index=False)
    terminal=dict(status='COMPLETE_ALL4_NATIVE_PANEL_PARITY_MATCHED_PROXY_CAPTURE_AND_AUTHOR_QUALITY_DESCRIPTIONS_PENDING_AUDIT',
             completed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planp),
             matched_contract_audit_sha256=sha(audit),native_input_images=4,native_cells=all_cells,
             full_native_RNA_positions_checked=sum(x['all_native_RNA_positions_replayed'] for x in raw_proofs),
             reference_official_cache_proxy_checks=len(reference_checks),reference_actual_cells=sum(x['actual_cached_reference_cells'] for x in reference_checks),
             target_independent_scalar_checks=numeric_checks,max_target_scalar_error=max_numeric_error,
             field_proxy_positions=len(f),author_subtype_proxy_positions=len(labels_rows),source_proxy_positions=len(summary),
             no_RNA_score_graph_arrays_written=True,new_spatial_effects_P_or_clinical_outcomes=False,formal_states_or_interactions_validated=False,
             scientific_outputs_public=False,files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(x['bytes'] for x in terminal['files'].values())<8*1024**2
    save('private_terminal.json',terminal);print(json.dumps({k:terminal[k] for k in ['status','native_input_images','native_cells','full_native_RNA_positions_checked','field_proxy_positions','author_subtype_proxy_positions']}))


if __name__=='__main__':main()

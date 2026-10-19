"""Nested question: background-only controls or other-state conditional controls."""
import hashlib,importlib.util,json,os,shutil,traceback
from concurrent.futures import ProcessPoolExecutor,as_completed
from multiprocessing import get_context
from datetime import datetime,timezone
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata
from sklearn.cluster import KMeans

from spatial_blocked_multiview_math_v1 import seed,predict_views,field_metrics,aggregate


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def load_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def scores(count,mappings,rows):
    out=np.full((len(count),len(mappings)),np.nan)
    bybg={}
    for j,m in enumerate(mappings):bybg.setdefault(tuple(m['background_indices']),[]).append((j,m))
    for bg,ms in bybg.items():
        lookup={k:i for i,k in enumerate(bg)}
        for start in range(0,len(rows),1024):
            ix=rows[start:start+1024];x=count[ix][:,list(bg)].astype(float);cap=min(1500,len(bg))
            ranks=np.minimum(cap,rankdata(-x,axis=1,method='average'));ranks[x<=0]=cap
            for j,m in ms:
                sig=[lookup[k] for k in m['signature_indices']];n=len(sig);base=n*(n+1)/2
                assert cap==m['maxrank'] and n<cap
                out[ix,j]=1-(ranks[:,sig].sum(1)-base)/(n*cap-base)
    return out


def moments(x):return np.column_stack([x,x*x])


def frames_for_image(target,capture,contract,plan,graph_fn):
    assert sha(target['analysis_path'])==target['input_sha256'] and sha(target['cache_path'])==target['cache_sha256']
    a=ad.read_h5ad(target['analysis_path'])
    with np.load(target['cache_path'],allow_pickle=False) as z:
        data={k:z[k] for k in ['author_labels','qc','fov','coords','area','native_total','native_ngenes','negatives','background_q95']}
    assert a.n_obs==target['native_cells'] and np.isfinite(data['coords']).all()
    np.testing.assert_array_equal(np.asarray(a.X.sum(1)).ravel(),data['native_total'])
    own={'T cell':contract['T_author_labels'],'Macrophage':contract['strict_macrophage_author_labels']}
    labels=data['author_labels'].astype(str);fovs=data['fov'].astype(str);qc=data['qc'].astype(bool)
    mappings=target['proxy_mappings'];all_ids=[m['proxy_group_id'] for m in mappings]
    q=capture[capture.sample_id.eq(target['sample_id'])].set_index(['FOV','proxy_group_id']);assert q.index.is_unique
    proof_path=Path('/dev/shm/pancancer-native-state-capture-20261020-v2/private_output/private_all_native_raw_parity_proofs.json')
    extra=next(x for x in json.loads(proof_path.read_text())['proofs'] if x['sample_id']==target['sample_id'])['metadata_only_counts_per_FOV']
    center=data['coords'][qc].mean(0);scale=data['coords'][qc].std(0);assert (scale>0).all()
    output=[];scope=[];score_checks=0;max_score_error=0.
    for fov in sorted(set(fovs)):
        full=np.flatnonzero((fovs==fov)&qc)
        if len(full)<3:
            scope.append(dict(FOV=fov,status='FEWER_THAN3_QC_RNA_CELLS'));continue
        xy=data['coords'][full].astype(float);_,inv,mul=np.unique(xy,axis=0,return_inverse=True,return_counts=True);duplicate=mul[inv]>1
        ix=full[~duplicate];xy=data['coords'][ix].astype(float)
        if len(ix)<3:
            scope.append(dict(FOV=fov,status='FEWER_THAN3_UNIQUE_QC_RNA_CELLS'));continue
        lab=labels[ix];pop={k:np.isin(lab,v) for k,v in own.items()}
        reason='FIELD_INCLUDED'
        if min(int(v.sum()) for v in pop.values())<30:reason='FEWER_THAN30_PER_LINEAGE_AFTER_GEOMETRY'
        elif not all(bool(q.loc[(fov,g),'native_capture_qualified']) for g in all_ids):reason='NOT_ALL_FROZEN_PROXY_MEASUREMENTS_CAPTURE_QUALIFIED'
        scope.append(dict(FOV=fov,status=reason,QC_RNA_cells=len(ix),duplicate_centroids_removed=int(duplicate.sum()),T_cells=int(pop['T cell'].sum()),Mac_cells=int(pop['Macrophage'].sum())))
        if reason!='FIELD_INCLUDED':continue
        near,ng=graph_fn(xy,6,2.);far,fg=graph_fn(xy,12,5.)
        raw=a.X[ix].toarray();chosen=np.flatnonzero(pop['T cell']|pop['Macrophage'])
        total=np.log1p(data['native_total'][ix]);det=np.log1p(data['native_ngenes'][ix]);area=np.log1p(data['area'][ix]);negative=np.log1p(data['negatives'][ix].mean(1))
        quality=np.column_stack([total,det,area,negative]);normxy=(xy-center)/scale
        geometry=np.column_stack([normxy,normxy[:,0]**2,normxy[:,1]**2,normxy[:,0]*normxy[:,1]])
        missing_fraction=int(extra.get(fov,0))/(int((fovs==fov).sum())+int(extra.get(fov,0)))
        backgrounds={}
        for receiver in own:
            sender='Macrophage' if receiver=='T cell' else 'T cell'
            bits=[quality,geometry,np.full((len(ix),1),missing_fraction)]
            bits.append(np.column_stack([(lab==x).astype(float) for x in own[receiver]]))
            for W in [near,far]:
                bits.append(np.log1p(W.getnnz(axis=1))[:,None])
                for name in own:bits.append(np.asarray(W[:,pop[name]].sum(1)))
                for name in [sender,receiver]:bits.append(W[:,pop[name]]@moments(quality[pop[name]]))
                bits.append(W[:,pop[sender]]@np.column_stack([(lab[pop[sender]]==x).astype(float) for x in own[sender]]))
            backgrounds[receiver]=np.column_stack(bits)
            assert np.isfinite(backgrounds[receiver]).all()
        readouts={}
        for readout in ['RAW','NEGATIVE_EXCESS']:
            x=raw if readout=='RAW' else np.maximum(raw-data['background_q95'][ix,None],0.)
            s=scores(x,mappings,chosen)
            if readout=='RAW' and not duplicate.any():
                for j,m in enumerate(mappings):
                    y=s[pop[m['major']],j];reference=q.loc[(fov,m['proxy_group_id'])]
                    values=[np.var(y),np.quantile(y,.75)-np.quantile(y,.25)]
                    expected=[reference.score_variance,reference.score_IQR]
                    np.testing.assert_allclose(values,expected,atol=1e-12,rtol=0)
                    max_score_error=max(max_score_error,float(np.max(np.abs(np.array(values)-expected))));score_checks+=1
            views={}
            for receiver in own:
                sender='Macrophage' if receiver=='T cell' else 'T cell'
                cols=[j for j,m in enumerate(mappings) if m['major']==sender]
                values=s[pop[sender]][:,cols];assert np.isfinite(values).all()
                # Weighted observed sender-state mass: structural zero only for no observed sender.
                views[receiver]=(near[:,pop[sender]]@values,far[:,pop[sender]]@values)
            readouts[readout]=dict(scores=s,views=views)
        output.append(dict(FOV=fov,xy=xy,pop=pop,backgrounds=backgrounds,readouts=readouts,
                           center=xy.mean(0),wide_radius=fg['radius_in_source_units'],native_rows=ix,
                           nearest_spacing=ng['median_nearest_neighbor_distance']))
    del a,data
    return output,scope,dict(cached_field_raw_score_checks=score_checks,max_raw_score_reconstruction_error=max_score_error)


def run_image(ordinal):
    control=Path(os.environ['PANCANCER_MV_CONTROL']);out=Path(os.environ['PANCANCER_RAM_OUTPUT']);plan=json.loads((control/'plan.json').read_text())
    contract=json.loads(Path(plan['private_contract']).read_text());target=contract['targets'][ordinal]
    capture=pd.read_csv('/dev/shm/pancancer-native-state-capture-20261020-v2/private_output/private_all_field_proxy_capture_and_quality.tsv.gz',sep='\t',dtype={'sample_id':str,'FOV':str})
    graph_module=load_module(control/'author_spatial_inflow_math_v1.py','native_graph')
    fields,scope,checks=frames_for_image(target,capture,contract,plan,graph_module.graph)
    mappings=target['proxy_mappings'];case_rows=[];fold_rows=[];field_rows=[];importance_rows=[]
    okay=len(fields)>=6 and len(np.unique(np.stack([x['center'] for x in fields]),axis=0))>=3 if fields else False
    if okay:
        assignment=KMeans(n_clusters=3,n_init=10,random_state=seed(ordinal,'FOV_GEOMETRY')).fit_predict(np.stack([x['center'] for x in fields]))
        guard=max(x['wide_radius'] for x in fields)
    else:assignment=np.full(len(fields),-1);guard=np.nan
    for j,m in enumerate(mappings):
        receiver=m['major'];sender='Macrophage' if receiver=='T cell' else 'T cell'
        others=[k for k,x in enumerate(mappings) if x['major']==receiver and k!=j] if plan['include_other_receiver_state_proxies'] else []
        sender_ids=[x['proxy_group_id'] for x in mappings if x['major']==sender]
        for readout in ['RAW','NEGATIVE_EXCESS']:
            base=dict(image_ordinal=ordinal,dataset_id=target['dataset_id'],cancer=target['cancer'],sample_id=target['sample_id'],
                 response_lineage=receiver,predictor_lineage=sender,target_proxy_group_id=m['proxy_group_id'],readout=readout,
                 eligible_FOVs=len(fields),guard_radius_in_source_units=float(guard),reviewed_independent_patients=0)
            all_metrics=[];valid_folds=0
            if okay:
                for fold in range(3):
                    test=[x for k,x in enumerate(fields) if assignment[k]==fold];train=[x for k,x in enumerate(fields) if assignment[k]!=fold]
                    test_tree=cKDTree(np.vstack([x['xy'] for x in test]));train_parts=[]
                    for f in train:
                        rows=np.flatnonzero(f['pop'][receiver]);distance=test_tree.query(f['xy'][rows],k=1)[0];rows=rows[distance>guard]
                        if len(rows)>=30:train_parts.append((f,rows))
                    record=dict(base,fold=fold,train_FOVs=len(train_parts),test_FOVs=len(test),status='FEWER_THAN3_BUFFERED_TRAIN_FIELDS')
                    if len(train_parts)<3:
                        fold_rows.append(record);continue
                    per_field=min(80,min(len(rows) for _,rows in train_parts));rng=np.random.default_rng(seed(ordinal,receiver,fold,'CELL_SELECTION'))
                    train_sets=[(f,np.sort(rng.choice(rows,per_field,replace=False))) for f,rows in train_parts]
                    def arrays(parts):
                        B=[];N=[];F=[];Y=[];FV=[];indices=[]
                        for f,rows in parts:
                            s=f['readouts'][readout]['scores'];b=np.column_stack([f['backgrounds'][receiver][rows],s[rows][:,others]])
                            n,far=f['readouts'][readout]['views'][receiver]
                            B.append(b);N.append(n[rows]);F.append(far[rows]);Y.append(s[rows,j]);FV.extend([f['FOV']]*len(rows));indices.extend(f['native_rows'][rows].tolist())
                        return np.vstack(B),np.vstack(N),np.vstack(F),np.concatenate(Y),np.asarray(FV),np.asarray(indices,np.int64)
                    tr=arrays(train_sets);te=arrays([(f,np.flatnonzero(f['pop'][receiver])) for f in test])
                    assert not set(tr[4])&set(te[4]) and not np.intersect1d(tr[5],te[5]).size
                    if np.var(tr[3])<=1e-12:
                        record['status']='CONSTANT_TRAIN_TARGET';fold_rows.append(record);continue
                    record.update(train_cells=len(tr[3]),test_cells=len(te[3]),train_cells_per_field=per_field,
                                  train_selection_sha256=hashlib.sha256(tr[5].astype('<i8').tobytes()).hexdigest())
                    b,full,imp=predict_views(tr[0],tr[1],tr[2],tr[3],te[0],te[1],te[2],seed(ordinal,receiver,fold,'RF'),trees=64,leaf=5)
                    fm=field_metrics(te[3],b,full,te[4]);fold_metric=aggregate(fm);record.update(fold_metric)
                    if fold_metric['status']=='ESTIMABLE_BLOCKED_PREDICTION':valid_folds+=1;all_metrics.extend(fm)
                    fold_rows.append(record)
                    for rr in fm:field_rows.append(dict(base,fold=fold,**rr))
                    for view,values in zip(['near','far'],imp):
                        for sid,value in zip(sender_ids,values):importance_rows.append(dict(base,fold=fold,view=view,sender_proxy_group_id=sid,impurity_importance_not_effect=float(value)))
            metric=aggregate(all_metrics) if valid_folds==3 else dict(status='INCOMPLETE_OR_INSUFFICIENT_SPATIAL_FOLDS')
            case=dict(base,valid_folds=valid_folds,test_fields=0,test_cells=0,baseline_R2=np.nan,multiview_R2=np.nan,gain_R2=np.nan,
                      equal_field_target_variance=np.nan,equal_field_baseline_MSE=np.nan,equal_field_multiview_MSE=np.nan)
            case.update(metric);case_rows.append(case)
            (out/('private_image_%02d_progress.json'%ordinal)).write_text(json.dumps(dict(stage='REAL_BLOCKED_MULTIVIEW_PREDICTION',image_ordinal=ordinal,completed_cases=len(case_rows),total_cases=len(mappings)*2,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    for name,records in [('models',case_rows),('folds',fold_rows),('test_fields',field_rows),('importances',importance_rows),('field_scope',scope)]:
        table=pd.DataFrame(records)
        if len(table.columns)==0:
            table=pd.DataFrame(columns=['image_ordinal','target_proxy_group_id','readout','status'] if name!='field_scope' else ['FOV','status'])
        table.to_csv(out/('private_image_%02d_%s.tsv.gz'%(ordinal,name)),sep='\t',index=False,compression='gzip')
    (out/('private_image_%02d_checks.json'%ordinal)).write_text(json.dumps(checks,indent=2)+'\n')
    return dict(image_ordinal=ordinal,status='COMPLETE',model_cases=len(case_rows))


def safe_image(ordinal):
    try:return run_image(ordinal)
    except Exception as e:
        out=Path(os.environ['PANCANCER_RAM_OUTPUT']);d=dict(image_ordinal=ordinal,status='EXECUTION_FAILED',exception_class=type(e).__name__,source_lines=[x.lineno for x in traceback.extract_tb(e.__traceback__)])
        (out/('private_image_%02d_execution_failure.json'%ordinal)).write_text(json.dumps(d,indent=2)+'\n');return d


def main():
    control=Path(os.environ['PANCANCER_MV_CONTROL']);out=Path(os.environ['PANCANCER_RAM_OUTPUT']);plan=json.loads((control/'plan.json').read_text())
    assert not out.exists() and out.parent==control and control.parent==Path('/dev/shm') and shutil.disk_usage(control).free>512*1024**2
    assert sha(__file__)==plan['producer_sha256']
    for name,digest in plan['helper_sha256'].items():assert sha(control/name)==digest
    audit=Path('/dev/shm/pancancer-native-state-capture-audit-20261020-v2/private_output/private_terminal.json');assert sha(audit)=='95693e8e457d7042c93694563504ad56ca61b2e903c5dee8d20a31071510999c'
    cap=Path('/dev/shm/pancancer-native-state-capture-20261020-v2/private_output');ct=json.loads((cap/'private_terminal.json').read_text());assert sha(cap/'private_terminal.json')=='efa79e09ea535391339ec81d4ecdbb0088500cea23d9ef003814d8c020ea2a1d'
    for name in ['private_all_field_proxy_capture_and_quality.tsv.gz','private_all_native_raw_parity_proofs.json']:assert sha(cap/name)==ct['files'][name]['sha256']
    contract_path=Path(plan['private_contract']);assert sha(contract_path)==plan['private_contract_sha256'];contract=json.loads(contract_path.read_text());assert len(contract['targets'])==4
    out.mkdir(mode=0o700)
    source_scope=Path('/dev/shm/pancancer-native-state-probe-match-20261020-v1/private_output/private_complete14_native_target_scope.tsv')
    shutil.copyfile(source_scope,out/'private_complete14_input_coverage.tsv')
    family=[dict(image_ordinal=i,target_proxy_group_id=m['proxy_group_id'],response_lineage=m['major'],readout=r) for i,t in enumerate(contract['targets']) for m in t['proxy_mappings'] for r in ['RAW','NEGATIVE_EXCESS']]
    assert len(family)==126
    (out/'private_frozen_full_family_before_effects.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),plan_sha256=sha(control/'plan.json'),family=family),indent=2)+'\n')
    results=[]
    with ProcessPoolExecutor(max_workers=2,mp_context=get_context('spawn')) as executor:
        futures=[executor.submit(safe_image,i) for i in range(4)]
        for f in as_completed(futures):
            results.append(f.result());(out/'private_progress.json').write_text(json.dumps(dict(stage='REAL_IMAGE_MODEL_RESULTS',images_finished=len(results),images_total=4,execution_failures=sum(x['status']!='COMPLETE' for x in results),observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    failed=[x for x in results if x['status']!='COMPLETE']
    models=[]
    if not failed:
        for i in range(4):models.append(pd.read_csv(out/('private_image_%02d_models.tsv.gz'%i),sep='\t'))
        allmodels=pd.concat(models,ignore_index=True);assert len(allmodels)==126
        allmodels.to_csv(out/'private_complete126_blocked_model_cases.tsv',sep='\t',index=False)
    terminal=dict(status='COMPLETE_ALL126_REAL_BLOCKED_MULTIVIEW_CASES_PENDING_AUDIT' if not failed else 'INCOMPLETE_EXECUTION_FAILURE_NOT_SCIENTIFIC_NEGATIVE',
             completed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(control/'plan.json'),
             frozen_family_positions=126,images_completed=len(results)-len(failed),image_execution_failures=failed,
             all_real_effects_pending_independent_and_scientific_review=True,no_RNA_scores_graphs_models_or_cell_predictions_persisted=True,
             formal_state_or_clinical_or_direct_interaction_validation=False,scientific_outputs_public=False,
             files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(d['bytes'] for d in terminal['files'].values())<32*1024**2
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps({k:terminal[k] for k in ['status','frozen_family_positions','images_completed']}))


if __name__=='__main__':main()

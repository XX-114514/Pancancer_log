"""Independent native/programme, radius graph, full10000 null and QR audit.

No producer import. Uses HDF5 CSR extra genes, query-ball geometry, independently
encoded nuisance frames and QR projections. Only aggregate evidence is public.
"""
from pathlib import Path
from datetime import datetime, timezone
import os, json, hashlib
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree


def digest(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''): h.update(b)
    return h.hexdigest()


def stamp(): return datetime.now(timezone.utc).isoformat()


def independent_graph(z):
    rows, cols, distances = [], [], []
    for fv in sorted(set(z['fov']), key=int):
        idx = np.where(z['qc'] & (z['fov'] == fv))[0]
        if not len(idx): continue
        xy = z['coords'][idx] / .12028; radius = np.sqrt(z['area'][idx] / np.pi)
        domains = cKDTree(xy).query_ball_point(xy, 24 * max(radius))
        for i, domain in enumerate(domains):
            jj = np.array(domain, int); jj = jj[jj != i]
            d = np.sqrt(((xy[jj] - xy[i]) ** 2).sum(axis=1)) / (radius[i] + radius[jj])
            okay = d <= 12
            rows.extend([idx[i]] * int(okay.sum())); cols.extend(idx[jj[okay]]); distances.extend(d[okay])
    rows, cols, distances = np.asarray(rows, int), np.asarray(cols, int), np.asarray(distances)
    assert np.all(rows != cols) and np.array_equal(z['fov'][rows], z['fov'][cols])
    return {k: csr_matrix((np.ones(int((distances <= k).sum()), np.int64), (rows[distances <= k], cols[distances <= k])), shape=(len(z['qc']), len(z['qc']))) for k in [2, 6, 12]}


def append_source_genes(z, source, saved):
    ex = dict(np.load(saved, allow_pickle=False)); genes = list(ex['genes'])
    with h5py.File(source, 'r') as f:
        var = f['var']; alias = {}
        index = var['_index'].asstr()[:]
        symbols = var['gene_symbol'].asstr()[:] if 'gene_symbol' in var and isinstance(var['gene_symbol'],h5py.Dataset) else index
        for j, (g, sym) in enumerate(zip(index, symbols)):
            for name in set([g, sym]): alias.setdefault(name, set()).add(j)
        raw = csr_matrix((f['X/data'][:], f['X/indices'][:], f['X/indptr'][:]), shape=tuple(f['X'].attrs['shape']))
        present = np.array([len(alias.get(g, [])) == 1 for g in genes])
        counts = np.full((len(z['qc']), len(genes)), -1, np.int32)
        for j, g in enumerate(genes):
            if present[j]: counts[:, j] = raw[:, next(iter(alias[g]))].toarray().ravel()
        assert np.array_equal(counts, ex['counts']) and np.array_equal(present, ex['present'])
    z['counts'] = np.concatenate([z['counts'], counts], axis=1); z['genes'] = np.concatenate([z['genes'], ex['genes']]); z['present'] = np.r_[z['present'], present]
    return int(present.sum()) * len(z['qc'])


def masks(z, row, pp, sid):
    empty = np.zeros(len(z['qc']), bool); gi = dict(zip(z['genes'], range(len(z['genes']))))
    a = 'T' if row.sender_family.startswith('T_') else 'B'; b = 'T' if row.receiver_family.startswith('T_') else 'B'
    if sid.startswith('GSM8746298'): return (empty,) * 4 + ('AUTHOR_NOT_RETRIEVED',) * 2
    if sid.startswith(('GSM9046088', 'GSM9046089')) and 'B' in [a,b]: return (empty,) * 4 + ('AUTHOR_B_CATEGORY_UNAVAILABLE',) * 2
    sg = row.sender_family.rsplit('_',1)[1]; rg = row.receiver_family.rsplit('_',1)[1]
    if not z['present'][gi[sg]] or not z['present'][gi[rg]]: return (empty,) * 4 + ('MISSING_NAME_GENE',) * 2
    positive = z['counts'] > np.maximum(0, z['background_q95'][:,None])
    s = np.array([x in pp['author_'+a+'_labels'] for x in z['author_labels']]) & z['qc'] & positive[:,gi[sg]]
    r = np.array([x in pp['author_'+b+'_labels'] for x in z['author_labels']]) & z['qc'] & positive[:,gi[rg]]
    lg = row.ligand.split('_'); rec = row.receptor.split('_')
    missing_l = [g for g in lg if not z['present'][gi[g]]]; missing_r = [g for g in rec if not z['present'][gi[g]]]
    l = empty if missing_l else s & positive[:,[gi[g] for g in lg]].all(axis=1)
    q = empty if missing_r else r & positive[:,[gi[g] for g in rec]].all(axis=1)
    why = 'MISSING_LIGAND:'+ ';'.join(missing_l) if missing_l else 'MISSING_RECEPTOR:'+ ';'.join(missing_r) if missing_r else 'MEASURED'
    rwhy = 'MISSING_RECEPTOR:'+ ';'.join(missing_r) if missing_r else 'MEASURED'
    return s,r,l,q,why,rwhy


def score(z, genes):
    gi = dict(zip(z['genes'], range(len(z['genes'])))); okay = [g for g in genes if z['present'][gi[g]]]
    enough = len(okay) >= 4 and len(okay)/len(genes) >= .6
    if enough:
        columns = [np.log1p(z['counts'][:,gi[g]]*10000/np.maximum(z['native_total'],1)) for g in okay]
        y = np.vstack(columns).mean(axis=0)
    else: y = np.full(len(z['qc']), np.nan)
    return y,okay,enough


def design(z, idx, deg, ns, s):
    names = ['logRNA','lognGenes','logArea','logDegree','logSenderDegree']
    if len(idx)<40: return None,'INSUFFICIENT_ELIGIBLE_RECEIVERS',names
    raw = np.array([np.log1p(z['native_total'][idx]),np.log1p(z['native_ngenes'][idx]),np.log(z['area'][idx]),np.log1p(deg[idx]),np.log1p(ns[idx])]).T
    std = np.std(raw,axis=0)
    if (std<=1e-12).any():return None,'NONVARYING_MANDATORY_COVARIATE',names
    columns = [np.ones(len(idx))] + [((raw[:,j]-raw[:,j].mean())/std[j]) for j in range(5)];names=['intercept']+names
    for label,values in [('FOV',z['fov'][idx]),('author',z['author_labels'][idx]),('sender_membership',s[idx].astype(str))]:
        levels=sorted(set(values))
        for level in levels[:-1]:columns.append((values==level).astype(float));names.append(label+'='+level)
    x=np.asarray(columns).T
    if len(idx)-x.shape[1]<21:return None,'INSUFFICIENT_RESIDUAL_DF',names
    if np.linalg.matrix_rank(x)<x.shape[1] or np.linalg.cond(x)>1e8:return None,'COLLINEAR_NUISANCE_DESIGN',names
    return x,'READY',names


def fit(x,y,e):
    if x is None:return None,None,None,'NO_DESIGN'
    if not np.isfinite(y).all() or np.var(y)<=1e-12:return None,None,None,'NONVARYING_PROGRAMME'
    full=np.column_stack([x,e])
    if np.linalg.matrix_rank(full)<full.shape[1] or np.linalg.cond(full)>1e8:return None,None,None,'COLLINEAR_EXPOSURE_DESIGN'
    q,r=np.linalg.qr(x,mode='reduced');ry=y-q@(q.T@y);re=e-q@(q.T@e)
    if re@re<=1e-12:return None,None,None,'NONVARYING_RESIDUAL_EXPOSURE'
    # Full augmented QR gives coefficient independently of production residual beta.
    qq,rr=np.linalg.qr(full,mode='reduced');beta=float(np.linalg.solve(rr,qq.T@y)[-1])
    return beta,q,ry,'ESTIMABLE_CONDITIONAL'


def replay_labels(z,s,r,receptor,l,deg,fi,li,saved):
    sx=np.where(s)[0]; keys=[]
    med={}
    for fv in sorted(set(z['fov']),key=int):
        idx=z['qc']&(z['fov']==fv)
        med[fv]=(np.median(z['native_total'][idx]) if idx.any() else 0,np.median(deg[idx]) if idx.any() else 0)
    for i in sx:
        depth,density=med[z['fov'][i]]
        keys.append((int(z['fov'][i]),int(r[i]),int(receptor[i]),str(z['author_labels'][i]),int(z['native_total'][i]>depth),int(deg[i]>density)))
    strata=sorted(set(keys)); originals=l[sx].astype(np.uint8); labels=np.empty((10000,len(sx)),np.uint8);move=0
    for j,group in enumerate(strata):
        positions=np.where(np.array([value==group for value in keys]))[0];values=originals[positions]
        if np.unique(values).size==1:labels[:,positions]=values
        else:
            generator=np.random.default_rng(np.random.SeedSequence([20261012,fi,li,j]));uniform=generator.random((10000,len(positions)))
            order=np.argsort(uniform,axis=1,kind='stable');labels[:,positions]=np.take(values,order);move+=len(positions)
        assert np.all(np.count_nonzero(labels[:,positions],axis=1)==np.count_nonzero(values))
    archive=dict(np.load(saved,allow_pickle=False));assert np.array_equal(sx,archive['sender_indices']) and np.array_equal(labels,archive['draws']) and np.array_equal(originals,archive['original'])
    assert np.array_equal(np.asarray(['|'.join(map(str,key)) for key in keys]),archive['stratum_keys'])
    return sx,labels,move


def null_coefficients(graph,sx,idx,ns,draws,q,ry):
    a=graph[idx][:,sx];all_b=[];invariant=0
    for k in range(0,10000,80):
        e=(a@draws[k:k+80].T).astype(float)/ns[idx,None]
        re=e-q@(q.T@e);norm=(re*re).sum(axis=0);result=np.zeros(len(norm));okay=norm>1e-12
        result[okay]=np.einsum('i,ij->j',ry,re[:,okay])/norm[okay];invariant+=int((~okay).sum());all_b.extend(result)
    return np.asarray(all_b),invariant


def main():
    start=stamp();root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);p=root/'immune_same_name_family_source_review_v1';out=p/'spatial_receiver_programme_exposure_v1';old=p/'spatial_contact_native_preflight_v3'
    assert not (out/'independent_audit.json').exists();v=json.loads((out/'inventory.json').read_text());plan=json.loads((out/'predeclared_plan.json').read_text())
    assert digest(out/'predeclared_plan.json')==v['plan_sha256']
    for n,h in v['outputs'].items():assert digest(out/n)==h
    inv=json.loads((old/'inventory.json').read_text());pp=json.loads((old/'predeclared_plan.json').read_text())
    assert digest(old/'inventory.json')==plan['parent']['native_inventory_sha256']
    d=pd.read_csv(out/'all70_primary_programme_exposure.tsv',sep='\t',float_precision='round_trip');radius=pd.read_csv(out/'all140_radius_sensitivity.tsv',sep='\t',float_precision='round_trip');reduced=pd.read_csv(out/'all70_programme_name_sensitivity.tsv',sep='\t',float_precision='round_trip');phen=pd.read_csv(out/'all70_receptor_phenotype.tsv',sep='\t',float_precision='round_trip');mapping=pd.read_csv(out/'all75_record_mapping.tsv',sep='\t')
    candidates=pd.read_csv(p/'axis_validation_feasibility_v1/automatic_next_validation_records.tsv',sep='\t');fixed=pd.read_csv(p/'supported_LR_native_cell_panel_v1/fixed336_reference_supported_candidates.tsv',sep='\t');fixed['model_key']=fixed[['context_index','program','sender_family','receiver_family','ligand','receptor']].astype(str).agg('|'.join,axis=1)
    original=candidates.merge(fixed[['model_key','program_genes']],on='model_key',validate='one_to_one');keys=['sender_family','receiver_family','ligand','receptor','program_genes'];models=original[keys].drop_duplicates().sort_values(keys).reset_index(drop=True);assert len(models)==14 and len(mapping)==75 and len(d)==70
    assert set(mapping.model_key)==set(original.model_key) and mapping[['sample_id','model_key']].drop_duplicates().shape[0]==75
    interfaces=original[keys[:4]].drop_duplicates().sort_values(keys[:4]).reset_index(drop=True);interfaces['li']=range(len(interfaces));models=models.merge(interfaces,on=keys[:4],validate='many_to_one')
    total_cells=0;extra_values=0;checks=[];null_models=0;cached_draws_replayed=0
    for fi,inp in enumerate(sorted(inv['inputs'],key=lambda a:a['sample_id'])):
        sid=inp['sample_id'];cache=old/'private_native_cell_caches'/(sid+'.npz')
        assert digest(cache)==v['input_sha256'][sid]['native_cache_sha256']==inv['outputs'][str(cache.relative_to(old))] and digest(inp['path'])==v['input_sha256'][sid]['native_source_sha256']
        z=dict(np.load(cache,allow_pickle=False));total_cells+=len(z['qc']);extra_values+=append_source_genes(z,inp['path'],out/'private_extra_gene_counts'/(sid+'.npz'))
        graph=independent_graph(z);degrees={k:np.asarray(a.sum(axis=1)).ravel() for k,a in graph.items()};label_cache={}
        for j,row in models.iterrows():
            mid='SP%03d'%(j+1);target=d[(d.sample_id==sid)&(d.programme_interface_id==mid)].iloc[0];s,r,l,q,why,rwhy=masks(z,row,pp,sid);y,genes,pg=score(z,row.program_genes.split(';'));assert target.programme_coverage_ok==pg
            main_frame=None
            for rad in [2,6,12]:
                ns=np.asarray(graph[rad]@s.astype(np.int64)).ravel();nl=np.asarray(graph[rad]@l.astype(np.int64)).ravel();idx=np.where(q&(ns>0))[0];exp=nl[idx]/ns[idx]
                x,status,names=design(z,idx,degrees[rad],ns,s)
                if why!='MEASURED':status=why;x=None
                elif not pg:status='INSUFFICIENT_PROGRAMME_PANEL';x=None
                beta,qq,ry,fst=fit(x,y[idx],exp)
                if status=='READY':status=fst
                rr=target if rad==6 else radius[(radius.sample_id==sid)&(radius.programme_interface_id==mid)&(radius.radius_multiplier==rad)].iloc[0]
                if rad!=6:assert rr.status==status
                assert (pd.isna(rr.beta_per_unit_L_fraction) and beta is None) or (beta is not None and np.isclose(rr.beta_per_unit_L_fraction,beta,atol=1e-9,rtol=1e-7))
                if rad==6:
                    main_frame=idx,exp,x,qq,ry,ns
                    model=dict(np.load(out/'private_model_caches'/(sid+'__'+mid+'.npz'),allow_pickle=False));assert np.array_equal(idx,model['eligible']) and np.allclose(exp,model['exposure']) and np.allclose(y[idx],model['y'],equal_nan=True)
                    if x is not None:assert np.allclose(x,model['X'],atol=1e-10) and np.array_equal(np.asarray(names),model['nuisance_names'])
                    else:assert model['X'].shape==(0,0)
                    if status=='ESTIMABLE_CONDITIONAL':
                        li='LI%03d'%(row.li+1)
                        if li not in label_cache:
                            sx,draws,moving=replay_labels(z,s,r,q,l,degrees[6],fi,int(row.li),out/'private_ligand_draws'/(sid+'__'+li+'.npz'));label_cache[li]=sx,draws,moving;cached_draws_replayed+=10000
                        sx,draws,moving=label_cache[li];null,invariant=null_coefficients(graph[6],sx,idx,ns,draws,qq,ry)
                        assert np.allclose(null,model['null_betas'],atol=1e-8,rtol=1e-6) and target.invariant_null_exposure_draws==invariant and target.movable_sender_cells==moving
                        if np.var(null)<=1e-12:assert target.status=='ZERO_CONDITIONAL_NULL_VARIANCE' and pd.isna(target.p_two_sided)
                        else:
                            assert target.status=='ESTIMABLE_CONDITIONAL';observed=rr.beta_per_unit_L_fraction
                            # Compare saved producer statistic to independently recomputed
                            # draws; near numerical ties use the saved exact replay values.
                            exact=model['null_betas'];upper=(1+(exact>=observed).sum())/10001;lower=(1+(exact<=observed).sum())/10001;expect=min(1.,2*min(upper,lower));assert np.isclose(target.p_two_sided,expect,atol=1e-12)
                            assert np.isclose(target.null_mean,null.mean(),atol=1e-8);null_models+=1
                    else:assert target.status==status and pd.isna(target.p_two_sided) and model['null_betas'].size==0
            idx,exp,x,qq,ry,ns=main_frame
            names=set([row.sender_family.rsplit('_',1)[1],row.receiver_family.rsplit('_',1)[1]]);gg=[g for g in row.program_genes.split(';') if g not in names];yy,gm,okay=score(z,gg);beta,_,_,st=fit(x,yy[idx],exp) if okay else (None,None,None,'INSUFFICIENT_PROGRAMME_PANEL');rr=reduced[(reduced.sample_id==sid)&(reduced.programme_interface_id==mid)].iloc[0]
            assert (pd.isna(rr.beta_per_unit_L_fraction) and beta is None) or (beta is not None and np.isclose(rr.beta_per_unit_L_fraction,beta,atol=1e-9,rtol=1e-7))
            ridx=np.where(r)[0];xx,st,nn=design(z,ridx,degrees[6],ns,s);beta,_,_,fst=fit(xx,y[ridx],q[ridx].astype(float)) if pg and rwhy=='MEASURED' else (None,None,None,'NO_RECEPTOR_OR_PROGRAMME');rr=phen[(phen.sample_id==sid)&(phen.programme_interface_id==mid)].iloc[0]
            assert (pd.isna(rr.receptor_positive_minus_negative_programme) and beta is None) or (beta is not None and np.isclose(rr.receptor_positive_minus_negative_programme,beta,atol=1e-9,rtol=1e-7))
            checks.append(dict(programme_interface_id=mid,sample_id=sid,status='PASS_NATIVE_PROGRAMME_EXPOSURE_ALL_RADII_QR_AND_CONDITIONAL_NULL',primary_status=target.status))
            print(json.dumps(dict(audited=mid,sample_id=sid,at_UTC=stamp())),flush=True)
    assert total_cells==254458 and len(checks)==70
    raw=d.p_two_sided.fillna(1).to_numpy();order=np.argsort(raw);q=np.minimum.accumulate((raw[order]*70/np.arange(1,71))[::-1])[::-1].clip(0,1);expected=np.empty(70);expected[order]=q;good=d.p_two_sided.notna();assert np.allclose(expected[good],d.loc[good,'BH_full70']) and d.loc[~good,'BH_full70'].isna().all()
    pd.DataFrame(checks).to_csv(out/'independent_all70_checks.tsv',sep='\t',index=False)
    audit=dict(schema='spatial_receiver_programme_audit_v1',status='PASS_ALL70_NATIVE_PROGRAMME_GEOMETRY_EXPOSURE_QR_AND_FULL10000_CONDITIONAL_NULL',started_UTC=start,finished_UTC=stamp(),cells=total_cells,models=70,mapped_rows=75,radius_checks=140,name_sensitivity_checks=70,receptor_phenotype_checks=70,native_extra_values_checked=extra_values,estimable_null_models=null_models,unique_ligand_draws_replayed=cached_draws_replayed,joint_counts_preserved=True,self_or_cross_FOV_edges=0,producer_inventory_sha256=digest(out/'inventory.json'),auditor_sha256=digest(__file__),patient_level_inference=False,independence='No producer import; native HDF5 CSR extra values and accepted parent cache, query-ball graphs, independently encoded covariates and QR coefficients. Shared frozen null specification and random-key seed replayed; numerical ties checked against saved full draws. Not independent biological replication.',scientific_release=False,axis_validation_complete=False)
    (out/'independent_audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit),flush=True)


if __name__=='__main__':main()

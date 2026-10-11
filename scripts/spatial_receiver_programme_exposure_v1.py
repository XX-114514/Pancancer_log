"""All frozen receiver programmes vs sender-ligand spatial exposure.

PANCANCER_PREFLIGHT_ROOT=<run> python spatial_receiver_programme_exposure_v1.py
Individual cells, nuisance frames and all permutations remain private.
Conditional tissue randomization is not independent-patient inference.
"""
from pathlib import Path
from datetime import datetime, timezone
import os, json, hashlib
import numpy as np
import pandas as pd
import anndata as ad
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''): h.update(b)
    return h.hexdigest()


def now(): return datetime.now(timezone.utc).isoformat()


def records(p):
    d = pd.read_csv(p / 'axis_validation_feasibility_v1/automatic_next_validation_records.tsv', sep='\t')
    f = pd.read_csv(p / 'supported_LR_native_cell_panel_v1/fixed336_reference_supported_candidates.tsv', sep='\t')
    f['model_key'] = f[['context_index', 'program', 'sender_family', 'receiver_family', 'ligand', 'receptor']].astype(str).agg('|'.join, axis=1)
    d = d.merge(f[['model_key', 'program_genes']], on='model_key', validate='one_to_one')
    assert len(d) == 15 and d.model_key.is_unique
    for row in d.itertuples():
        assert set(row.program_genes.split(';')).isdisjoint(set(row.ligand.split('_') + row.receptor.split('_') + [row.receiver_family.rsplit('_', 1)[1]]))
    keys = ['sender_family', 'receiver_family', 'ligand', 'receptor', 'program_genes']
    definitions = d[keys].drop_duplicates().sort_values(keys).reset_index(drop=True)
    definitions['programme_interface_id'] = ['SP%03d' % (i + 1) for i in range(len(definitions))]
    assert len(definitions) == 14
    interfaces = d[keys[:4]].drop_duplicates().sort_values(keys[:4]).reset_index(drop=True)
    interfaces['ligand_interface_id'] = ['LI%03d' % (i + 1) for i in range(len(interfaces))]
    definitions = definitions.merge(interfaces, on=keys[:4], validate='many_to_one')
    return d.merge(definitions, on=keys, validate='many_to_one'), definitions


def extra_native(inp, z, extra, target):
    a = ad.read_h5ad(inp['path']); assert a.n_obs == len(z['qc'])
    assert np.array_equal(np.asarray(a.obs_names, str), z['cell_key'])
    lookup = {}
    for i, name in enumerate(a.var_names):
        aliases = {str(name)}
        if 'gene_symbol' in a.var: aliases.add(str(a.var.iloc[i].gene_symbol))
        for alias in aliases: lookup.setdefault(alias, set()).add(i)
    present = np.array([len(lookup.get(g, [])) == 1 for g in extra], bool)
    raw = np.full((a.n_obs, len(extra)), -1, np.int32)
    cols = [next(iter(lookup[g])) for g in extra if len(lookup.get(g, [])) == 1]
    raw[:, present] = a.X[:, cols].toarray()
    assert (raw[:, present] >= 0).all()
    np.savez_compressed(target, genes=np.asarray(extra), counts=raw, present=present)
    z['genes'] = np.r_[z['genes'], np.asarray(extra)]
    z['counts'] = np.column_stack([z['counts'], raw]); z['present'] = np.r_[z['present'], present]
    return [{'sample_id': inp['sample_id'], 'gene': g, 'aliases': len(lookup.get(g, [])), 'measurable': bool(present[j]), 'positive_cells': int((raw[:, j] > 0).sum()) if present[j] else None} for j, g in enumerate(extra)]


def graphs(z):
    all_i, all_j, all_d = [], [], []
    for fov in sorted(set(z['fov']), key=int):
        ix = np.flatnonzero(z['qc'] & (z['fov'] == fov))
        if not len(ix): continue
        xy = z['coords'][ix] / .12028; radius = np.sqrt(z['area'][ix] / np.pi)
        pairs = cKDTree(xy).query_pairs(24 * radius.max(), output_type='ndarray')
        if not len(pairs): continue
        ratio = np.linalg.norm(xy[pairs[:, 0]] - xy[pairs[:, 1]], axis=1) / (radius[pairs[:, 0]] + radius[pairs[:, 1]])
        mask = ratio <= 12; pairs, ratio = pairs[mask], ratio[mask]
        all_i.extend(np.r_[ix[pairs[:, 0]], ix[pairs[:, 1]]]); all_j.extend(np.r_[ix[pairs[:, 1]], ix[pairs[:, 0]]]); all_d.extend(np.r_[ratio, ratio])
    i, j, ratio = np.array(all_i, int), np.array(all_j, int), np.asarray(all_d)
    assert (i != j).all() and (z['fov'][i] == z['fov'][j]).all()
    n = len(z['qc'])
    return {k: csr_matrix((np.ones(int((ratio <= k).sum()), np.int64), (i[ratio <= k], j[ratio <= k])), shape=(n, n)) for k in [2, 6, 12]}


def endpoint_masks(z, row, pp, sid):
    n = len(z['qc']); no = np.zeros(n, bool); gi = {g: i for i, g in enumerate(z['genes'])}
    sk = 'T' if row.sender_family.startswith('T_') else 'B'; rk = 'T' if row.receiver_family.startswith('T_') else 'B'
    if sid.startswith('GSM8746298'): return no, no, no, no, 'AUTHOR_NOT_RETRIEVED', 'AUTHOR_NOT_RETRIEVED'
    if sid.startswith(('GSM9046088', 'GSM9046089')) and 'B' in [sk, rk]:
        return no, no, no, no, 'AUTHOR_B_CATEGORY_UNAVAILABLE', 'AUTHOR_B_CATEGORY_UNAVAILABLE'
    names = [row.sender_family.rsplit('_', 1)[1], row.receiver_family.rsplit('_', 1)[1]]
    if not all(z['present'][gi[g]] for g in names): return no, no, no, no, 'MISSING_NAME_GENE', 'MISSING_NAME_GENE'
    pos = (z['counts'] > 0) & (z['counts'] > z['background_q95'][:, None])
    s = z['qc'] & np.isin(z['author_labels'], pp['author_' + sk + '_labels']) & pos[:, gi[names[0]]]
    r = z['qc'] & np.isin(z['author_labels'], pp['author_' + rk + '_labels']) & pos[:, gi[names[1]]]
    lg = row.ligand.split('_'); rg = row.receptor.split('_')
    lm = [g for g in lg if not z['present'][gi[g]]]; rm = [g for g in rg if not z['present'][gi[g]]]
    l = no if lm else s & pos[:, [gi[g] for g in lg]].all(axis=1)
    q = no if rm else r & pos[:, [gi[g] for g in rg]].all(axis=1)
    status = 'MISSING_LIGAND:' + ';'.join(lm) if lm else 'MISSING_RECEPTOR:' + ';'.join(rm) if rm else 'MEASURED'
    rstatus = 'MISSING_RECEPTOR:' + ';'.join(rm) if rm else 'MEASURED'
    return s, r, l, q, status, rstatus


def programme(z, genes):
    gi = {g: i for i, g in enumerate(z['genes'])}; measured = [g for g in genes if z['present'][gi[g]]]
    good = len(measured) >= 4 and len(measured) / len(genes) >= .6
    y = np.mean(np.log1p(z['counts'][:, [gi[g] for g in measured]] / np.maximum(z['native_total'], 1)[:, None] * 10000), axis=1) if good else np.full(len(z['qc']), np.nan)
    return y, measured, good


def nuisance(z, ix, degree, sender_degree, smask):
    n = len(ix); fields = [np.log1p(z['native_total'][ix]), np.log1p(z['native_ngenes'][ix]), np.log(z['area'][ix]), np.log1p(degree[ix]), np.log1p(sender_degree[ix])]
    names = ['logRNA', 'lognGenes', 'logArea', 'logDegree', 'logSenderDegree']
    raw = np.column_stack(fields) if n else np.empty((0, 5))
    if n < 40: return None, 'INSUFFICIENT_ELIGIBLE_RECEIVERS', names
    sd = raw.std(axis=0)
    if (sd <= 1e-12).any(): return None, 'NONVARYING_MANDATORY_COVARIATE', names
    x = np.column_stack([np.ones(n), (raw - raw.mean(axis=0)) / sd]); names = ['intercept'] + names
    for name, values in [('FOV', z['fov'][ix]), ('author', z['author_labels'][ix]), ('sender_membership', smask[ix].astype(str))]:
        levels = sorted(set(values))
        for level in levels[:-1]: x = np.column_stack([x, values == level]); names.append(name + '=' + level)
    if n - x.shape[1] < 21: return None, 'INSUFFICIENT_RESIDUAL_DF', names
    if np.linalg.matrix_rank(x) < x.shape[1] or np.linalg.cond(x) > 1e8: return None, 'COLLINEAR_NUISANCE_DESIGN', names
    return x, 'READY', names


def adjusted_fit(x, y, exposure):
    if x is None: return None, None, None, 'NO_DESIGN'
    if not np.isfinite(y).all() or np.var(y) <= 1e-12: return None, None, None, 'NONVARYING_PROGRAMME'
    augmented = np.column_stack([x, exposure])
    if np.linalg.matrix_rank(augmented) < augmented.shape[1] or np.linalg.cond(augmented) > 1e8: return None, None, None, 'COLLINEAR_EXPOSURE_DESIGN'
    u, _, _ = np.linalg.svd(x, full_matrices=False)
    ey = y - u @ (u.T @ y); ee = exposure - u @ (u.T @ exposure)
    denom = ee @ ee
    if denom <= 1e-12: return None, None, None, 'NONVARYING_RESIDUAL_EXPOSURE'
    return float(ee @ ey / denom), u, ey, 'ESTIMABLE_CONDITIONAL'


def sender_draws(z, smask, rmask, receptor, ligand, degree, file_index, ligand_index):
    sx = np.flatnonzero(smask); group_keys = []
    med = {}
    for fv in sorted(set(z['fov']), key=int):
        q = z['qc'] & (z['fov'] == fv)
        med[fv] = (np.median(z['native_total'][q]) if q.any() else 0, np.median(degree[q]) if q.any() else 0)
    for i in sx:
        dep, dense = med[z['fov'][i]]
        group_keys.append((int(z['fov'][i]), int(rmask[i]), int(receptor[i]), str(z['author_labels'][i]), int(z['native_total'][i] > dep), int(degree[i] > dense)))
    levels = sorted(set(group_keys)); groups = {key: [] for key in levels}
    for j, key in enumerate(group_keys): groups[key].append(j)
    original = ligand[sx].astype(np.uint8); draws = np.empty((10000, len(sx)), np.uint8); movable = 0
    for si, key in enumerate(levels):
        idx = np.asarray(groups[key], int); code = original[idx]
        if code.min() == code.max(): draws[:, idx] = code
        else:
            rng = np.random.default_rng(np.random.SeedSequence([20261012, file_index, ligand_index, si]))
            order = np.argsort(rng.random((10000, len(idx))), axis=1, kind='stable')
            draws[:, idx] = code[order]; movable += len(idx)
        assert np.all(draws[:, idx].sum(axis=1) == code.sum())
    return sx, draws, movable, group_keys


def permutation_betas(adj, sx, ix, sender_degree, ligand_draws, u, ey):
    b = adj[ix][:, sx]; denom = sender_degree[ix]
    values = np.empty(10000); invariant = 0
    for k in range(0, 10000, 100):
        e = (b @ ligand_draws[k:k + 100].T).astype(float) / denom[:, None]
        er = e - u @ (u.T @ e); dd = (er * er).sum(axis=0); okay = dd > 1e-12
        bb = np.zeros(len(dd)); bb[okay] = (ey @ er[:, okay]) / dd[okay]
        values[k:k + len(dd)] = bb; invariant += int((~okay).sum())
    return values, invariant


def main():
    start = now(); root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); p = root / 'immune_same_name_family_source_review_v1'
    old = p / 'spatial_contact_native_preflight_v3'; contact = p / 'spatial_contact_conditional_neighborhood_v1'; out = p / 'spatial_receiver_programme_exposure_v1'
    assert out.exists() and not (out / 'inventory.json').exists()
    plan = json.loads((out / 'predeclared_plan.json').read_text()); parent = plan['parent']
    assert sha(contact / 'inventory.json') == parent['contact_inventory_sha256'] and sha(contact / 'independent_audit.json') == parent['contact_audit_sha256']
    assert sha(old / 'inventory.json') == parent['native_inventory_sha256'] and sha(p / 'axis_validation_feasibility_v1/automatic_next_validation_records.tsv') == parent['candidate_sha256']
    inv = json.loads((old / 'inventory.json').read_text()); pp = json.loads((old / 'predeclared_plan.json').read_text()); original, definitions = records(p)
    assert inv['cells'] == 254458 and len(inv['inputs']) == 5 and plan['null']['permutations'] == 10000
    model_dir = out / 'private_model_caches'; model_dir.mkdir(); input_dir = out / 'private_extra_gene_counts'; input_dir.mkdir(); label_dir = out / 'private_ligand_draws'; label_dir.mkdir()
    primary, radii, reduced, phenotype, map_rows, cov_rows, hashes = [], [], [], [], [], [], {}
    for fi, inp in enumerate(sorted(inv['inputs'], key=lambda x: x['sample_id'])):
        sid = inp['sample_id']; path = old / 'private_native_cell_caches' / (sid + '.npz')
        assert sha(path) == inv['outputs'][str(path.relative_to(old))]
        assert sha(inp['path']) == inp['sha256']; hashes[sid] = dict(native_cache_sha256=sha(path),native_source_sha256=inp['sha256'])
        z = dict(np.load(path, allow_pickle=False)); needed = set()
        for row in definitions.itertuples(): needed.update(row.program_genes.split(';') + row.ligand.split('_') + row.receptor.split('_'))
        extra = sorted(needed - set(z['genes'])); cov_rows.extend(extra_native(inp, z, extra, input_dir / (sid + '.npz')))
        graph = graphs(z); degrees = {rad: np.asarray(g.sum(axis=1)).ravel() for rad, g in graph.items()}; cached = {}
        for row in definitions.itertuples():
            key = row.programme_interface_id; s, r, l, q, why, rwhy = endpoint_masks(z, row, pp, sid)
            y, measured, pg_ok = programme(z, row.program_genes.split(';'))
            base = dict(programme_interface_id=key, ligand_interface_id=row.ligand_interface_id, sample_id=sid, sender_family=row.sender_family, receiver_family=row.receiver_family, ligand=row.ligand, receptor=row.receptor, program_genes=row.program_genes, measurable_programme_genes=';'.join(measured), programme_coverage_ok=pg_ok, programme_remaining_genes=len(row.program_genes.split(';')), sender_family_cells=int(s.sum()) if not why.startswith('AUTHOR') else None, receiver_family_cells=int(r.sum()) if not rwhy.startswith('AUTHOR') else None, receptor_positive=int(q.sum()) if rwhy=='MEASURED' else None, patient_level_inference=False)
            model_frames = {}; row_radii = {}
            for rad in [2, 6, 12]:
                degree_s = np.asarray(graph[rad] @ s.astype(np.int64)).ravel(); number_l = np.asarray(graph[rad] @ l.astype(np.int64)).ravel()
                ix = np.flatnonzero(q & (degree_s > 0)); e = number_l[ix] / degree_s[ix]
                x, status, names = nuisance(z, ix, degrees[rad], degree_s, s)
                if why != 'MEASURED': status = why; x = None
                elif not pg_ok: status = 'INSUFFICIENT_PROGRAMME_PANEL'; x = None
                beta, u, ey, fit_status = adjusted_fit(x, y[ix], e)
                if status == 'READY': status = fit_status
                diag = dict(base, radius_multiplier=rad, status=status, eligible_receivers=len(ix) if rwhy=='MEASURED' else None, receivers_no_sender_neighbour=int((q & (degree_s==0)).sum()) if rwhy=='MEASURED' else None, exposed_receivers=int((e>0).sum()) if why=='MEASURED' else None, unexposed_receivers=int((e==0).sum()) if why=='MEASURED' else None, beta_per_unit_L_fraction=beta, nuisance_columns=len(names), covariate_names=';'.join(names), conditional_null_p_not_computed=rad!=6)
                model_frames[rad] = (ix, e, x, u, ey, degree_s); row_radii[rad] = diag
                if rad != 6: radii.append(diag)
            main = row_radii[6]; ix,e,x,u,ey,degree_s = model_frames[6]
            betas = np.empty(0); invariants = None; moving = None; sx = np.empty(0,int)
            if main['status'] == 'ESTIMABLE_CONDITIONAL':
                lid = row.ligand_interface_id
                if lid not in cached:
                    sx, draws, moving, group_keys = sender_draws(z,s,r,q,l,degrees[6],fi,int(lid[2:])-1)
                    np.savez_compressed(label_dir/(sid+'__'+lid+'.npz'),sender_indices=sx,draws=draws,original=l[sx],stratum_keys=np.asarray(['|'.join(map(str,g)) for g in group_keys]))
                    cached[lid] = sx,draws,moving
                sx, draws, moving = cached[lid]
                betas,invariants = permutation_betas(graph[6],sx,ix,degree_s,draws,u,ey)
                if np.var(betas)<=1e-12: main['status']='ZERO_CONDITIONAL_NULL_VARIANCE'; pv=None
                else:
                    upper=(1+int((betas>=main['beta_per_unit_L_fraction']).sum()))/10001;lower=(1+int((betas<=main['beta_per_unit_L_fraction']).sum()))/10001;pv=min(1.,2*min(upper,lower))
                main.update(p_two_sided=pv,null_mean=float(betas.mean()),null_sd=float(betas.std(ddof=1)),null_q025=float(np.quantile(betas,.025)),null_q975=float(np.quantile(betas,.975)),invariant_null_exposure_draws=invariants,movable_sender_cells=moving)
            else: main.update(p_two_sided=None,null_mean=None,null_sd=None,null_q025=None,null_q975=None,invariant_null_exposure_draws=None,movable_sender_cells=None)
            # Saved native indices/model frame permit independent reconstruction.
            np.savez_compressed(model_dir/(sid+'__'+key+'.npz'),eligible=ix,exposure=e,y=y[ix],X=x if x is not None else np.empty((0,0)),null_betas=betas,nuisance_names=np.asarray(main['covariate_names'].split(';')),status=np.asarray(main['status']))
            names = set([row.sender_family.rsplit('_',1)[1],row.receiver_family.rsplit('_',1)[1]])
            genes = [g for g in row.program_genes.split(';') if g not in names]
            yy, measured_reduced, good = programme(z,genes)
            b,_,_,st = adjusted_fit(x,yy[ix],e) if good else (None,None,None,'INSUFFICIENT_PROGRAMME_PANEL')
            reduced.append(dict(base,status=st if x is not None else main['status'],programme_name_removed_genes=';'.join(genes),measurable_genes=';'.join(measured_reduced),beta_per_unit_L_fraction=b))
            ridx=np.flatnonzero(r);degree_s6=model_frames[6][-1];xx,st,nn=nuisance(z,ridx,degrees[6],degree_s6,s)
            rb,_,_,rs=adjusted_fit(xx,y[ridx],q[ridx].astype(float)) if pg_ok and rwhy=='MEASURED' else (None,None,None,'INSUFFICIENT_PROGRAMME_PANEL' if not pg_ok else rwhy)
            phenotype.append(dict(base,status=rs if st=='READY' else st if rwhy=='MEASURED' else rwhy,receptor_positive_minus_negative_programme=rb,pvalue_not_computed=True))
            main['leave_one_FOV_beta_min']=None;main['leave_one_FOV_beta_max']=None;main['leave_one_FOV_estimable']=0
            if x is not None and main['status'].startswith('ESTIMABLE'):
                bs=[]
                for fv in sorted(set(z['fov'][ix]),key=int):
                    keep=z['fov'][ix]!=fv; sub=ix[keep]; xx,st,nn=nuisance(z,sub,degrees[6],degree_s,s)
                    b,_,_,fs=adjusted_fit(xx,y[sub],e[keep])
                    if st=='READY' and fs=='ESTIMABLE_CONDITIONAL':bs.append(b)
                if bs:main.update(leave_one_FOV_beta_min=float(min(bs)),leave_one_FOV_beta_max=float(max(bs)),leave_one_FOV_estimable=len(bs))
            primary.append(main)
            for orig in original[original.programme_interface_id==key].itertuples(): map_rows.append(dict(model_key=orig.model_key,sample_id=sid,programme_interface_id=key,ligand_interface_id=row.ligand_interface_id,original_contact_role=bool(orig.eligible_contact_expression_followup),program=orig.program,discovery_cancer=orig.cancerType))
            print(json.dumps(dict(completed=key,sample_id=sid,status=main['status'],beta=main['beta_per_unit_L_fraction'],p=main['p_two_sided'],at_UTC=now())),flush=True)
    d=pd.DataFrame(primary);assert len(d)==70 and len(map_rows)==75 and len(radii)==140 and len(reduced)==70 and len(phenotype)==70
    raw=d.p_two_sided.fillna(1).to_numpy();order=np.argsort(raw);adj=np.minimum.accumulate((raw[order]*70/np.arange(1,71))[::-1])[::-1].clip(0,1);q=np.empty(70);q[order]=adj;d['BH_full70']=np.where(d.p_two_sided.notna(),q,np.nan)
    for name,df in [('all70_primary_programme_exposure.tsv',d),('all75_record_mapping.tsv',pd.DataFrame(map_rows)),('all140_radius_sensitivity.tsv',pd.DataFrame(radii)),('all70_programme_name_sensitivity.tsv',pd.DataFrame(reduced)),('all70_receptor_phenotype.tsv',pd.DataFrame(phenotype)),('native_gene_coverage.tsv',pd.DataFrame(cov_rows))]:df.to_csv(out/name,sep='\t',index=False,na_rep='NA')
    outputs={str(f.relative_to(out)):sha(f) for f in out.rglob('*') if f.is_file()}
    inventory=dict(schema='spatial_receiver_programme_inventory_v1',status='PRODUCER_COMPLETE_PENDING_INDEPENDENT_AUDIT',started_UTC=start,finished_UTC=now(),producer_sha256=sha(__file__),plan_sha256=sha(out/'predeclared_plan.json'),input_sha256=hashes,outputs=outputs,cells=254458,files=5,models=70,mapped_rows=75,permutations=10000,estimable=int(d.status.eq('ESTIMABLE_CONDITIONAL').sum()),primary_BH_q05=int((d.BH_full70<.05).sum()),scientific_release=False,axis_validation_complete=False)
    (out/'inventory.json').write_text(json.dumps(inventory,indent=2)+'\n');print(json.dumps({k:v for k,v in inventory.items() if k not in ['outputs','input_sha256']}),flush=True)


if __name__ == '__main__':main()

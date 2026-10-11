from pathlib import Path
from datetime import datetime, timezone
import os, json, hashlib, importlib.util
import numpy as np
import pandas as pd

R = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
P = R / 'immune_same_name_family_source_review_v1'
O = P / 'spatial_receiver_programme_exposure_v1'
OLD = P / 'spatial_contact_native_preflight_v3'
MOD = R / 'spatial_receiver_programme_code_v1/audit_spatial_receiver_programme_exposure_v1.py'

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p): return json.loads(Path(p).read_text())
def frame(p): return pd.read_csv(p, sep='\t', float_precision='round_trip')
def same(a, b):
    assert (pd.isna(a) and b is None) or (b is not None and np.isclose(a, b, rtol=1e-7, atol=1e-9)), (a, b)

assert not (O / 'accepted_primary_consistency_checks.json').exists()
start = datetime.now(timezone.utc).isoformat()
a = read(O / 'independent_audit.json')
assert sha(MOD) == a['auditor_sha256']
spec = importlib.util.spec_from_file_location('accepted_receiver_audit', MOD)
api = importlib.util.module_from_spec(spec); spec.loader.exec_module(api)
d = frame(O / 'all70_primary_programme_exposure.tsv')
rad = frame(O / 'all140_radius_sensitivity.tsv')
name = frame(O / 'all70_programme_name_sensitivity.tsv')
phen = frame(O / 'all70_receptor_phenotype.tsv')
mapping = frame(O / 'all75_record_mapping.tsv')
cov = frame(O / 'native_gene_coverage.tsv')
inv = read(OLD / 'inventory.json'); pp = read(OLD / 'predeclared_plan.json')
orig = frame(P / 'axis_validation_feasibility_v1/automatic_next_validation_records.tsv')
fixed = frame(P / 'supported_LR_native_cell_panel_v1/fixed336_reference_supported_candidates.tsv')
fixed['model_key'] = fixed[['context_index','program','sender_family','receiver_family','ligand','receptor']].astype(str).agg('|'.join, axis=1)
orig = orig.merge(fixed[['model_key','program_genes']], on='model_key', validate='one_to_one')
keys = ['sender_family','receiver_family','ligand','receptor','program_genes']
definitions = orig[keys].drop_duplicates().sort_values(keys).reset_index(drop=True)
definitions['programme_interface_id'] = ['SP%03d' % (j+1) for j in range(14)]
lr = orig[keys[:4]].drop_duplicates().sort_values(keys[:4]).reset_index(drop=True)
lr['ligand_interface_id'] = ['LI%03d' % (j+1) for j in range(len(lr))]
definitions = definitions.merge(lr, on=keys[:4], validate='many_to_one')
orig = orig.merge(definitions, on=keys, validate='many_to_one')
sources = sorted(inv['inputs'], key=lambda x: x['sample_id'])
expected = pd.concat([orig.assign(sample_id=inp['sample_id']) for inp in sources], ignore_index=True)
mapcols = ['model_key','sample_id','programme_interface_id','ligand_interface_id']
pd.testing.assert_frame_equal(mapping[mapcols].sort_values(mapcols).reset_index(drop=True), expected[mapcols].sort_values(mapcols).reset_index(drop=True), check_dtype=False)
assert mapping[['sample_id','model_key']].duplicated().sum() == 0
role = mapping.merge(orig[['model_key','eligible_contact_expression_followup','program','cancerType']], on='model_key', validate='many_to_one')
assert (role.original_contact_role == role.eligible_contact_expression_followup).all()
assert (role.program_x == role.program_y).all() and (role.discovery_cancer == role.cancerType).all()
assert role.original_contact_role.sum() == 65
assert len(d) == 70 and len(rad) == 140 and len(name) == 70 and len(phen) == 70 and len(cov) == 35
for tab, extra in [(d,[]),(rad,['radius_multiplier']),(name,[]),(phen,[])]:
    assert not tab[['sample_id','programme_interface_id']+extra].duplicated().any()
    assert set(tab.sample_id) == {i['sample_id'] for i in sources}

checks = []
for inp in sources:
    sid = inp['sample_id']
    z = dict(np.load(OLD / 'private_native_cell_caches' / (sid+'.npz'), allow_pickle=False))
    ex = dict(np.load(O / 'private_extra_gene_counts' / (sid+'.npz'), allow_pickle=False))
    for j,g in enumerate(ex['genes']):
        c = cov[(cov.sample_id == sid) & (cov.gene == g)].iloc[0]
        assert c.measurable == ex['present'][j]
        same(c.positive_cells, int((ex['counts'][:,j] > 0).sum()) if ex['present'][j] else None)
    z['counts'] = np.column_stack([z['counts'],ex['counts']]);z['genes'] = np.r_[z['genes'],ex['genes']];z['present'] = np.r_[z['present'],ex['present']]
    graphs = api.independent_graph(z)
    degrees = {k:np.asarray(v.sum(axis=1)).ravel() for k,v in graphs.items()}
    for row in definitions.itertuples():
        mid = row.programme_interface_id
        primary = d[(d.sample_id == sid) & (d.programme_interface_id == mid)].iloc[0]
        for col in keys + ['ligand_interface_id']:
            assert str(primary[col]) == str(getattr(row,col))
        s,r,l,q,why,rwhy = api.masks(z,row,pp,sid)
        y, measured, enough = api.score(z,row.program_genes.split(';'))
        assert primary.programme_remaining_genes == len(row.program_genes.split(';'))
        assert primary.measurable_programme_genes == ';'.join(measured)
        same(primary.sender_family_cells, int(s.sum()) if not why.startswith('AUTHOR') else None)
        same(primary.receiver_family_cells, int(r.sum()) if not rwhy.startswith('AUTHOR') else None)
        same(primary.receptor_positive, int(q.sum()) if rwhy == 'MEASURED' else None)
        main = None
        for radius in [2,6,12]:
            ns = np.asarray(graphs[radius] @ s.astype(int)).ravel(); nl = np.asarray(graphs[radius] @ l.astype(int)).ravel()
            idx = np.flatnonzero(q & (ns > 0)); e = nl[idx]/ns[idx]
            x, status, names = api.design(z,idx,degrees[radius],ns,s)
            if why != 'MEASURED': status = why; x = None
            elif not enough: status = 'INSUFFICIENT_PROGRAMME_PANEL'; x = None
            b,_,_,fst = api.fit(x,y[idx],e)
            if status == 'READY': status = fst
            rr = primary if radius == 6 else rad[(rad.sample_id == sid)&(rad.programme_interface_id == mid)&(rad.radius_multiplier == radius)].iloc[0]
            assert rr.status == status
            same(rr.eligible_receivers, len(idx) if rwhy == 'MEASURED' else None)
            same(rr.receivers_no_sender_neighbour, int((q & (ns == 0)).sum()) if rwhy == 'MEASURED' else None)
            same(rr.exposed_receivers, int((e > 0).sum()) if why == 'MEASURED' else None)
            same(rr.unexposed_receivers, int((e == 0).sum()) if why == 'MEASURED' else None)
            assert rr.covariate_names == ';'.join(names) and rr.nuisance_columns == len(names)
            if radius == 6: main = idx,e,x,ns
        idx,e,x,ns = main
        leave = []
        if x is not None and primary.status.startswith('ESTIMABLE'):
            for fv in sorted(set(z['fov'][idx]), key=int):
                keep = z['fov'][idx] != fv; sub = idx[keep]
                xx,st,_ = api.design(z,sub,degrees[6],ns,s)
                b,_,_,fst = api.fit(xx,y[sub],e[keep])
                if st == 'READY' and fst == 'ESTIMABLE_CONDITIONAL': leave.append(b)
        same(primary.leave_one_FOV_beta_min, min(leave) if leave else None)
        same(primary.leave_one_FOV_beta_max, max(leave) if leave else None)
        assert primary.leave_one_FOV_estimable == len(leave)
        if primary.status == 'ESTIMABLE_CONDITIONAL':
            cache = dict(np.load(O / 'private_model_caches' / (sid+'__'+mid+'.npz'), allow_pickle=False))
            bnull = cache['null_betas']
            for col,value in [('null_sd',bnull.std(ddof=1)),('null_q025',np.quantile(bnull,.025)),('null_q975',np.quantile(bnull,.975))]:same(primary[col],value)
        ng = [g for g in row.program_genes.split(';') if g not in {row.sender_family.rsplit('_',1)[1],row.receiver_family.rsplit('_',1)[1]}]
        yy,gm,good = api.score(z,ng)
        b,_,_,st = api.fit(x,yy[idx],e) if good else (None,None,None,'INSUFFICIENT_PROGRAMME_PANEL')
        target = name[(name.sample_id == sid)&(name.programme_interface_id == mid)].iloc[0]
        assert target.status == (st if x is not None else primary.status)
        assert target.programme_name_removed_genes == ';'.join(ng) and target.measurable_genes == ';'.join(gm)
        ri = np.flatnonzero(r);xx,st,_ = api.design(z,ri,degrees[6],ns,s)
        b,_,_,fst = api.fit(xx,y[ri],q[ri].astype(float)) if enough and rwhy == 'MEASURED' else (None,None,None,'INSUFFICIENT_PROGRAMME_PANEL' if not enough else rwhy)
        target = phen[(phen.sample_id == sid)&(phen.programme_interface_id == mid)].iloc[0]
        expect = fst if st == 'READY' else st if rwhy == 'MEASURED' else rwhy
        assert target.status == expect and target.pvalue_not_computed
        checks.append(dict(sample_id=sid,programme_interface_id=mid,status='PASS_MAPPING_COVERAGE_ALL_DIAGNOSTIC_STATUS_DENOMINATORS_NULL_SUMMARY_AND_LOFO',primary_status=primary.status,leave_one_FOV_estimable=len(leave)))
    print(json.dumps(dict(accepted_source=sid,at_UTC=datetime.now(timezone.utc).isoformat())),flush=True)
pd.DataFrame(checks).to_csv(O / 'accepted_all70_consistency_checks.tsv',sep='\t',index=False)
summary = dict(schema='spatial_receiver_acceptance_consistency_v1',status='PASS_ALL70_MAPPING_DENOMINATORS_DIAGNOSTIC_STATUSES_AND_LOFO',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),models=70,mapped_rows=75,contact_mapping_rows=65,programme_only_mapping_rows=10,radius_rows=140,name_rows=70,phenotype_rows=70,native_extra_gene_rows=35,checks=70,method='Separate acceptance-stage table and private-native-cache checks using previously audited independent query-ball graph/design/QR functions; not a third independent biological validation',acceptance_script_sha256=sha(__file__),scientific_release=False,axis_validation_complete=False)
(O / 'accepted_primary_consistency_checks.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary),flush=True)

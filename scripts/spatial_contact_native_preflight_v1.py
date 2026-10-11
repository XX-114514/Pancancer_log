"""Frozen spatial contact RNA preflight. Individual cells stay outside Git.

Usage: PANCANCER_PREFLIGHT_ROOT=<run> PANCANCER_AUTHOR_ROOT=<RDS-dir>
       python spatial_contact_native_preflight_v1.py
No proximity inference or state assignment is performed.
"""
from pathlib import Path
from datetime import datetime, timezone
import os, json, hashlib, tarfile, gzip, re
import numpy as np
import pandas as pd
import anndata as ad
import pyreadr


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8388608), b''):
            h.update(chunk)
    return h.hexdigest()


def stamp():
    return datetime.now(timezone.utc).isoformat()


def native_integer(values):
    s = pd.Series(values).astype(str)
    a = pd.to_numeric(s, errors='raise').to_numpy()
    assert np.isfinite(a).all() and np.equal(a, np.floor(a)).all()
    return a.astype(np.int64).astype(str)


def author_maps(root):
    output = {}
    for name in ['ICON1_ICON2_CosMx_Annotations.rds', 'MESO1_MESO2_CosMx_annotations.rds']:
        d = list(pyreadr.read_r(root / name).values())[0]
        assert list(d.columns) == ['Barcodes', 'Annotations'] and d.Barcodes.is_unique
        parts = d.Barcodes.str.extract(r'^(Icon[12]|Meso[12])_(\d+)_(\d+)$', flags=re.I)
        assert parts.notna().all().all()
        for i, (prefix, fov, cell) in parts.iterrows():
            key = (prefix.lower(), str(int(fov)), str(int(cell)))
            assert key not in output
            output[key] = str(d.loc[i, 'Annotations'])
    return output


def aliases(a, genes):
    lookup = {}
    for j, index in enumerate(a.var_names):
        names = [str(index)]
        if 'gene_symbol' in a.var:
            names.append(str(a.var.iloc[j]['gene_symbol']))
        for name in set(names):
            if name not in {'', 'nan', 'None'}:
                lookup.setdefault(name, set()).add(j)
    return {g: sorted(lookup.get(g, set())) for g in genes}


def raw_controls(archive, member, a):
    fovs = native_integer(a.obs['fov']); cells = native_integer(a.obs['cell_ID'])
    lookup = {k: i for i, k in enumerate(zip(fovs, cells))}
    assert len(lookup) == a.n_obs
    with tarfile.open(archive, 'r:*') as t:
        with gzip.GzipFile(fileobj=t.extractfile(member)) as stream:
            header = pd.read_csv(stream, nrows=0).columns
        names = sorted(g for g in header if re.match(r'^(NegPrb\d+|Negative\d+)$', g, flags=re.I))
        assert len(names) == 10, (archive.name, names)
        controls = np.full((a.n_obs, len(names)), -1, dtype=np.int32)
        seen = np.zeros(a.n_obs, bool)
        with gzip.GzipFile(fileobj=t.extractfile(member)) as stream:
            for frame in pd.read_csv(stream, usecols=['fov', 'cell_ID'] + names, chunksize=2000):
                fv = native_integer(frame.fov); ci = native_integer(frame.cell_ID)
                valid = ci != '0'; frame = frame.loc[valid]
                keys = list(zip(fv[valid], ci[valid]))
                assert all(k in lookup for k in keys)
                ix = np.array([lookup[k] for k in keys], int)
                assert len(set(ix)) == len(ix) and not seen[ix].any()
                vals = frame[names].to_numpy()
                assert np.isfinite(vals).all() and (vals >= 0).all() and np.equal(vals, np.floor(vals)).all()
                controls[ix] = vals; seen[ix] = True
    assert seen.all() and (controls >= 0).all()
    return controls, names, fovs, cells


def compile_records(p):
    candidate = p / 'axis_validation_feasibility_v1/automatic_next_validation_records.tsv'
    records = pd.read_csv(candidate, sep='\t')
    records = records.loc[records.eligible_contact_expression_followup.eq(True)].copy()
    fixed = pd.read_csv(p / 'supported_LR_native_cell_panel_v1/fixed336_reference_supported_candidates.tsv', sep='\t')
    fixed['model_key'] = fixed[['context_index', 'program', 'sender_family', 'receiver_family', 'ligand', 'receptor']].astype(str).agg('|'.join, axis=1)
    records = records.merge(fixed[['model_key', 'program_genes']], on='model_key', validate='one_to_one', how='left')
    assert len(records) == 13 and records.model_key.is_unique and records.program_genes.notna().all()
    return records


def main():
    start = stamp()
    r = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); p = r / 'immune_same_name_family_source_review_v1'
    out = p / 'spatial_contact_native_preflight_v1'; assert out.is_dir() and not (out / 'inventory.json').exists()
    plan = json.loads((out / 'predeclared_plan.json').read_text())
    accepted = plan['accepted_parent']; manifest_path = p / 'spatial_native_input_repair_v2/analysis_input_manifest.tsv'
    assert sha(manifest_path) == accepted['analysis_manifest_sha256']
    assert sha(p / 'axis_validation_feasibility_v1/automatic_next_validation_records.tsv') == accepted['candidate_table_sha256']
    assert sha(p / 'spatial_native_input_repair_v2/independent_audit.json') == accepted['repair_audit_sha256']
    records = compile_records(p)
    genes = set(plan['family_name_genes'])
    for panel in plan['panels'].values(): genes.update(panel)
    for row in records.to_dict('records'):
        genes.update(row['ligand'].split('_')); genes.update(row['receptor'].split('_')); genes.update(row['program_genes'].split(';'))
    genes = sorted(genes); gi = {g: i for i, g in enumerate(genes)}
    manifest = pd.read_csv(manifest_path, sep='\t')
    files = manifest.loc[manifest.dataset_id.isin(['GSE287472', 'GSE299786'])].copy(); assert len(files) == 5
    author_root = Path(os.environ['PANCANCER_AUTHOR_ROOT']); author = author_maps(author_root)
    for x in plan['author_inputs']: assert sha(author_root / x['name']) == x['sha256']
    cache_dir = out / 'private_native_cell_caches'; cache_dir.mkdir(exist_ok=False)
    gene_rows, coord_rows, record_rows, fov_rows, author_rows, file_inputs = [], [], [], [], [], []
    for file in files.to_dict('records'):
        path = Path(file['analysis_path']); before = path.stat(); a = ad.read_h5ad(path)
        assert a.var_names.is_unique and a.obs_names.is_unique
        assert np.isfinite(a.X.data).all() and (a.X.data >= 0).all() and np.equal(a.X.data, np.floor(a.X.data)).all()
        assert not any(re.fullmatch(r'Unnamed:\s*\d+(?:\.\d+)?', str(x), re.I) for x in a.var_names)
        alias = aliases(a, genes); present = np.array([len(alias[g]) == 1 for g in genes], bool)
        counts = np.full((a.n_obs, len(genes)), -1, np.int32)
        cols = [alias[g][0] for g in genes if len(alias[g]) == 1]
        counts[:, present] = a.X[:, cols].toarray()
        totals = np.asarray(a.X.sum(axis=1)).ravel(); detected = np.asarray((a.X > 0).sum(axis=1)).ravel()
        archive = r.parents[1] / 'data_expansion_20260909/raw/spatial' / file['dataset_id'] / (file['dataset_id'] + '_RAW.tar')
        member = a.uns['archive_member_map']['expr']
        neg, neg_names, fovs, cells = raw_controls(archive, member, a)
        fov_levels = sorted(set(fovs), key=int); threshold = np.zeros(a.n_obs)
        for fov in fov_levels: threshold[fovs == fov] = np.quantile(neg[fovs == fov].reshape(-1), .95)
        area = pd.to_numeric(a.obs.Area, errors='coerce').to_numpy(float)
        valid_area = np.isfinite(area) & (area > 0); assert valid_area.any()
        area_limit = 5 * np.exp(np.mean(np.log(area[valid_area])))
        qc = (totals > 30) & valid_area & (area < area_limit)
        coords = np.asarray(a.obsm['spatial']); assert coords.shape == (a.n_obs, 2) and np.isfinite(coords).all()
        expected_coords = a.obs[['CenterX_global_px', 'CenterY_global_px']].to_numpy(float) * .12028
        assert np.allclose(coords, expected_coords, rtol=0, atol=1e-6)
        sid = file['sample_id']; prefix = {'GSM9046088': 'icon1', 'GSM9046089': 'icon2', 'GSM9046090': 'meso1', 'GSM9046091': 'meso2'}.get(sid.split('_')[0])
        labels = np.array([author.get((prefix, f, c), 'NOT_IN_AUTHOR_QC_OBJECT') if prefix else 'AUTHOR_NOT_RETRIEVED' for f, c in zip(fovs, cells)])
        np.savez_compressed(cache_dir / (sid + '.npz'), counts=counts, genes=np.array(genes), present=present, native_total=totals, native_ngenes=detected, negatives=neg, negative_genes=np.array(neg_names), fov=fovs, cell_id=cells, cell_key=np.array(a.obs_names.astype(str)), coords=coords, area=area, area_limit=np.array(area_limit), qc=qc, background_q95=threshold, author_labels=labels)
        coordinate_fields = [c for c in a.obs if re.search(r'patient|donor|subject|core', c, re.I)]
        coord_rows.append(dict(dataset_id=file['dataset_id'],sample_id=sid,cells=a.n_obs,genes=a.n_vars,FOVs=len(fov_levels),RNA_zero=int((totals==0).sum()),RNA_median=float(np.median(totals)),QC_cells=int(qc.sum()),negative_controls=len(neg_names),negative_control_median=float(np.median(neg.sum(axis=1))),background_threshold_min=float(threshold.min()),background_threshold_max=float(threshold.max()),coordinates_finite=True,converter_scale_matches=True,physical_calibration_independent=False,patient_core_fields=';'.join(coordinate_fields) or 'NONE',patient_mapping_reviewed=False,patient_level_inference='UNAVAILABLE_UNMAPPED',mixed_lineage_interpretation='RNA_coexpression_proxy_not_confirmed_segmentation_error',author_mapped_cells=int((~np.isin(labels,['NOT_IN_AUTHOR_QC_OBJECT','AUTHOR_NOT_RETRIEVED'])).sum()),author_labels_missing_is_negative=False))
        for j,g in enumerate(genes):
            available=present[j]; v=counts[:,j]
            gene_rows.append(dict(dataset_id=file['dataset_id'],sample_id=sid,gene=g,alias_matches=len(alias[g]),status='MEASURABLE' if available else ('MISSING' if len(alias[g])==0 else 'AMBIGUOUS'),native_cells=a.n_obs,positive_all=int((v>0).sum()) if available else None,positive_QC=int(((v>0)&qc).sum()) if available else None,count_ge2_QC=int(((v>=2)&qc).sum()) if available else None,above_background_QC=int(((v>threshold)&(v>0)&qc).sum()) if available else None))
        for row in records.to_dict('records'):
            forbidden=set(row['ligand'].split('_')+row['receptor'].split('_')+plan['family_name_genes'])
            panels={k:[g for g in panel if g not in forbidden and present[gi[g]]] for k,panel in plan['panels'].items()}
            lrgenes=row['ligand'].split('_')+row['receptor'].split('_')
            familygenes=[row['sender_family'].split('_')[-1],row['receiver_family'].split('_')[-1]]
            sline='T' if row['sender_family'].startswith('T_') else 'B';rline='T' if row['receiver_family'].startswith('T_') else 'B'
            required=lrgenes+familygenes;missing=[g for g in required if not present[gi[g]]]
            pg=row['program_genes'].split(';');measurable=[g for g in pg if present[gi[g]]]
            programme_ok=len(measurable)>=4 and len(measurable)/len(pg)>=.6
            result=dict(model_key=row['model_key'],dataset_id=file['dataset_id'],sample_id=sid,ligand=row['ligand'],receptor=row['receptor'],sender_family=row['sender_family'],receiver_family=row['receiver_family'],program=row['program'],discovery_cancer=row['cancerType'],spatial_cancer=file['cancer'],context_relation='SAME_BROAD_CANCER_LUNG_SUBTYPE_TRANSFER' if row['cancerType']=='Lung cancer' and 'Lung' in sid else 'CROSS_CANCER_TRANSFER',missing_required_genes=';'.join(missing) or 'NONE',sender_lineage_measurable_genes=len(panels[sline]),receiver_lineage_measurable_genes=len(panels[rline]),programme_measurable_genes=len(measurable),programme_remaining_genes=len(pg),programme_coverage_ok=programme_ok,patient_inference='UNAVAILABLE_UNMAPPED')
            endpoint_measurable=not missing and len(panels[sline])>=2 and len(panels[rline])>=2
            for arm in plan['arms']:
                pos=(counts>0)&present[None,:]
                if arm=='qc_count_ge2':pos=(counts>=2)&present[None,:]
                if 'background' in arm:pos=(counts>threshold[:,None])&(counts>0)&present[None,:]
                keep=np.ones(a.n_obs,bool) if arm=='all_detected' else qc
                lineage={k:pos[:,[gi[g] for g in pp]].sum(axis=1)>=2 for k,pp in panels.items()}
                smix=np.logical_or.reduce([lineage[k] for k in panels if k!=sline]);rmix=np.logical_or.reduce([lineage[k] for k in panels if k!=rline])
                baseS=lineage[sline]&pos[:,gi[familygenes[0]]]&keep;baseR=lineage[rline]&pos[:,gi[familygenes[1]]]&keep
                send=baseS&np.all(pos[:,[gi[g] for g in row['ligand'].split('_')]],axis=1)
                receive=baseR&np.all(pos[:,[gi[g] for g in row['receptor'].split('_')]],axis=1)
                s_mixed=int((send&smix).sum());r_mixed=int((receive&rmix).sum())
                if arm=='qc_background_pure':send &= ~smix;receive &= ~rmix
                suffix='_'+arm
                for key,value in [('sender_family_cells',baseS.sum()),('receiver_family_cells',baseR.sum()),('sender_positive',send.sum()),('receiver_positive',receive.sum()),('overlap_same_cell',(send&receive).sum()),('sender_mixed_before_pure',s_mixed),('receiver_mixed_before_pure',r_mixed)]:result[key+suffix]=int(value) if endpoint_measurable else None
                eligible=0
                for fov in fov_levels:
                    fm=fovs==fov; ns=int((send&fm).sum()); nr=int((receive&fm).sum()); no=int((send&receive&fm).sum()); ok=endpoint_measurable and ns>=20 and nr>=20
                    eligible+=ok
                    fov_rows.append(dict(model_key=row['model_key'],dataset_id=file['dataset_id'],sample_id=sid,fov=fov,arm=arm,native_cells=int(fm.sum()),QC_cells=int((fm&qc).sum()),sender_positive=ns if endpoint_measurable else None,receiver_positive=nr if endpoint_measurable else None,overlap_same_cell=no if endpoint_measurable else None,both_endpoints_ge20=ok,endpoint_measurable=endpoint_measurable,programme_coverage_ok=programme_ok,patient_mapping_reviewed=False))
                result['FOVs_both_endpoints_ge20'+suffix]=int(eligible) if endpoint_measurable else None
                result['status'+suffix]='UNMEASURABLE_ENDPOINT' if not endpoint_measurable else ('DESCRIPTIVE_FOV_COVERAGE' if eligible else 'INSUFFICIENT_PAIRED_FOV_ENDPOINTS')
                if arm=='qc_background':
                    for endpoint, mask,line in [('sender',send,sline),('receiver',receive,rline)]:
                        valid=~np.isin(labels,['NOT_IN_AUTHOR_QC_OBJECT','AUTHOR_NOT_RETRIEVED']);want=plan['author_T_labels'] if line=='T' else plan['author_B_labels']; author_same=np.isin(labels,want)
                        author_rows.append(dict(model_key=row['model_key'],dataset_id=file['dataset_id'],sample_id=sid,endpoint=endpoint,marker_endpoint_cells=int(mask.sum()) if endpoint_measurable else None,author_annotation_matched=int((mask&valid).sum()) if endpoint_measurable else None,author_same_lineage=int((mask&author_same).sum()) if endpoint_measurable else None,author_other_lineage=int((mask&valid&~author_same).sum()) if endpoint_measurable else None,author_unmapped=int((mask&~valid).sum()) if endpoint_measurable else None,author_identity_sensitivity_only=True))
            record_rows.append(result)
        assert (before.st_size,before.st_mtime_ns)==(path.stat().st_size,path.stat().st_mtime_ns)
        file_inputs.append(dict(sample_id=sid,path=str(path),sha256=sha(path),cells=a.n_obs,genes=a.n_vars,archive=str(archive),archive_sha256=sha(archive),expression_member=member,negative_names=neg_names))
        print(json.dumps(dict(completed_file=sid,cells=a.n_obs,QC=int(qc.sum()),finished_UTC=stamp())),flush=True)
    assert sum(x['cells'] for x in file_inputs)==254458 and len(record_rows)==65
    tables={'all_selected_gene_cell_QC.tsv':gene_rows,'coordinate_donor_and_segmentation_review.tsv':coord_rows,'all_contact_record_testability.tsv':record_rows,'all_contact_record_FOV_counts.tsv':fov_rows,'author_lineage_agreement.tsv':author_rows}
    for name,rows in tables.items():pd.DataFrame(rows).to_csv(out/name,sep='\t',index=False,na_rep='NA',lineterminator='\n')
    outputs={str(q.relative_to(out)):sha(q) for q in out.rglob('*') if q.is_file()}
    inv=dict(schema='spatial_contact_native_preflight_inventory_v1',status='PRODUCER_COMPLETE_PENDING_INDEPENDENT_AUDIT',started_UTC=start,finished_UTC=stamp(),plan_sha256=sha(out/'predeclared_plan.json'),producer_sha256=sha(__file__),inputs=file_inputs,outputs=outputs,cells=254458,files=5,unique_contact_records=13,record_assay_checks=65,unique_LR=9,arms=plan['arms'],scientific_release=False,spatial_contact_validated=False,patient_independence_validated=False)
    (out/'inventory.json').write_text(json.dumps(inv,indent=2)+'\n')
    print(json.dumps({k:v for k,v in inv.items() if k not in ['inputs','outputs']}),flush=True)


if __name__=='__main__':main()

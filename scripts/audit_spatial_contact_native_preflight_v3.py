"""Independent native HDF5/CSV recount, without importing producer helpers."""
from pathlib import Path
from datetime import datetime, timezone
import os, json, hashlib, gzip, tarfile, re
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
import pyreadr


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8388608), b''):h.update(b)
    return h.hexdigest()


def read_vector(obj):
    if isinstance(obj,h5py.Group):
        v=read_vector(obj['categories']); c=obj['codes'][()]; return np.array([v[i] if i>=0 else 'nan' for i in c])
    a=obj[()]
    return a.astype(str) if a.dtype.kind in 'OSU' else a


def number_strings(v):
    a=np.asarray(v).astype(float); assert np.isfinite(a).all() and (a==np.floor(a)).all()
    return a.astype(np.int64).astype(str)


def check(actual, expected, label):
    if expected is None:assert pd.isna(actual), (label,actual)
    elif isinstance(expected,float):assert np.isclose(actual,expected,rtol=1e-12,atol=1e-12), (label,actual,expected)
    else:assert actual==expected, (label,actual,expected)


def main():
    start=datetime.now(timezone.utc).isoformat();r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);p=r/'immune_same_name_family_source_review_v1';out=p/'spatial_contact_native_preflight_v3'
    assert not (out/'independent_audit.json').exists()
    inv=json.loads((out/'inventory.json').read_text());plan=json.loads((out/'predeclared_plan.json').read_text())
    for name,sha in inv['outputs'].items():assert digest(out/name)==sha
    records=pd.read_csv(p/'axis_validation_feasibility_v1/automatic_next_validation_records.tsv',sep='\t');records=records[records.eligible_contact_expression_followup.eq(True)]
    fixed=pd.read_csv(p/'supported_LR_native_cell_panel_v1/fixed336_reference_supported_candidates.tsv',sep='\t')
    fixed['model_key']=fixed[['context_index','program','sender_family','receiver_family','ligand','receptor']].astype(str).agg('|'.join,axis=1)
    records=records.merge(fixed[['model_key','program_genes']],on='model_key',validate='one_to_one');assert len(records)==13
    result=pd.read_csv(out/'all_contact_record_testability.tsv',sep='\t').set_index(['sample_id','model_key'],verify_integrity=True);assert len(result)==65
    fovtable=pd.read_csv(out/'all_contact_record_FOV_counts.tsv',sep='\t').set_index(['sample_id','model_key','arm','fov'],verify_integrity=True)
    # pandas can infer fov numbers as ints; retain canonical strings for index comparison.
    fovtable.index=fovtable.index.set_levels(fovtable.index.levels[-1].astype(str),level=-1)
    gene_table=pd.read_csv(out/'all_selected_gene_cell_QC.tsv',sep='\t').set_index(['sample_id','gene'],verify_integrity=True)
    coord=pd.read_csv(out/'coordinate_donor_and_segmentation_review.tsv',sep='\t').set_index('sample_id',verify_integrity=True)
    agreement=pd.read_csv(out/'author_lineage_agreement.tsv',sep='\t').set_index(['sample_id','model_key','endpoint'],verify_integrity=True)
    author_parts=[]
    for item in plan['author_inputs']:
        f=Path(os.environ['PANCANCER_AUTHOR_ROOT'])/item['name'];assert digest(f)==item['sha256']
        frame=list(pyreadr.read_r(f).values())[0]
        parts=frame.Barcodes.str.rsplit('_',n=2,expand=True);assert parts.shape[1]==3
        parts[0]=parts[0].str.lower();parts[1]=number_strings(parts[1]);parts[2]=number_strings(parts[2])
        frame=frame.assign(prefix=parts[0],fov=parts[2],cell=parts[1]);author_parts.append(frame.set_index(['prefix','fov','cell'],verify_integrity=True).Annotations)
    author=pd.concat(author_parts);assert author.index.is_unique
    audited=[];value_count=0;control_count=0;arm_checks=0;fov_checks=0
    for inputfile in inv['inputs']:
        assert digest(inputfile['path'])==inputfile['sha256'] and digest(inputfile['archive'])==inputfile['archive_sha256']
        sid=inputfile['sample_id'];cache=np.load(out/'private_native_cell_caches'/(sid+'.npz'),allow_pickle=False);genes=list(cache['genes']);gi={g:i for i,g in enumerate(genes)}
        with h5py.File(inputfile['path'],'r') as h:
            assert h['X'].attrs['encoding-type']=='csr_matrix'
            x=h['X'];native=sp.csr_matrix((x['data'][()],x['indices'][()],x['indptr'][()]),shape=tuple(x.attrs['shape']))
            names=read_vector(h['var/_index']).astype(str);assert len(set(names))==len(names)
            symbols=read_vector(h['var/gene_symbol']).astype(str) if 'gene_symbol' in h['var'] else names
            amap={g:sorted(set(np.flatnonzero((names==g)|(symbols==g)))) for g in genes}
            present=np.array([len(amap[g])==1 for g in genes]);assert np.array_equal(present,cache['present'])
            for j,g in enumerate(genes):
                if present[j]:assert np.array_equal(cache['counts'][:,j],native[:,amap[g][0]].toarray().ravel());value_count+=native.shape[0]
                else:assert (cache['counts'][:,j]==-1).all()
            fulltot=np.asarray(native.sum(axis=1)).ravel();fullng=np.asarray((native>0).sum(axis=1)).ravel()
            assert np.array_equal(fulltot,cache['native_total']) and np.array_equal(fullng,cache['native_ngenes'])
            fovs=number_strings(read_vector(h['obs/fov']));cells=number_strings(read_vector(h['obs/cell_ID']))
            assert np.array_equal(fovs,cache['fov']) and np.array_equal(cells,cache['cell_id'])
            assert np.array_equal(read_vector(h['obs/_index']).astype(str),cache['cell_key'])
            assert np.array_equal(h['obsm/spatial'][()],cache['coords'])
            area=np.asarray(read_vector(h['obs/Area']),float);assert np.array_equal(area,cache['area'],equal_nan=True)
        # Whole control frame and MultiIndex reindex, independently of producer chunk/key loop.
        with tarfile.open(inputfile['archive'],'r:*') as t:
            with gzip.GzipFile(fileobj=t.extractfile(inputfile['expression_member'])) as f:
                frame=pd.read_csv(f,usecols=lambda c:c in ['fov','cell_ID'] or bool(re.match(r'^(Negative[0-9]+|NegPrb[0-9]+)$',c,re.I)))
        frame['fov']=number_strings(frame.fov);frame['cell_ID']=number_strings(frame.cell_ID);frame=frame[frame.cell_ID!='0'];frame=frame.set_index(['fov','cell_ID'],verify_integrity=True)
        neg_names=sorted(frame.columns);assert len(neg_names)==10 and neg_names==list(cache['negative_genes'])
        keys=pd.MultiIndex.from_arrays([fovs,cells]);assert keys.is_unique and set(frame.index)==set(keys)
        neg=frame.reindex(keys)[neg_names].to_numpy();assert np.array_equal(neg,cache['negatives']);control_count+=neg.size
        q95=np.zeros(len(fovs));levels=sorted(set(fovs),key=int)
        for fov in levels:
            z=np.sort(neg[fovs==fov].reshape(-1));k=(len(z)-1)*.95;lo=int(np.floor(k));hi=int(np.ceil(k));q95[fovs==fov]=z[lo]+(z[hi]-z[lo])*(k-lo)
        assert np.allclose(q95,cache['background_q95'],rtol=1e-14,atol=1e-14)
        valid=np.isfinite(area)&(area>0);limit=5*np.exp(np.log(area[valid]).mean());qc=(fulltot>30)&valid&(area<limit)
        assert np.array_equal(qc,cache['qc']) and np.isclose(limit,float(cache['area_limit']),rtol=1e-14)
        prefix={'GSM9046088':'icon1','GSM9046089':'icon2','GSM9046090':'meso1','GSM9046091':'meso2'}.get(sid.split('_')[0])
        if prefix:
            keys=pd.MultiIndex.from_arrays([[prefix]*len(fovs),fovs,cells]);subset_author=author[author.index.get_level_values('prefix')==prefix];assert set(subset_author.index).issubset(set(keys));labels=author.reindex(keys).fillna('NOT_IN_AUTHOR_QC_OBJECT').to_numpy(str);assert np.count_nonzero(labels!='NOT_IN_AUTHOR_QC_OBJECT')==len(subset_author)
        else:labels=np.repeat('AUTHOR_NOT_RETRIEVED',len(fovs))
        assert np.array_equal(labels,cache['author_labels'])
        c=coord.loc[sid];check(c.cells,len(fovs),sid);check(c.QC_cells,int(qc.sum()),sid);check(c.author_mapped_cells,int((~np.isin(labels,['NOT_IN_AUTHOR_QC_OBJECT','AUTHOR_NOT_RETRIEVED'])).sum()),sid)
        assert not c.patient_mapping_reviewed and c.patient_level_inference=='UNAVAILABLE_UNMAPPED'
        for j,g in enumerate(genes):
            gt=gene_table.loc[(sid,g)];v=cache['counts'][:,j];check(gt.alias_matches,len(amap[g]),g)
            for key,mask in [('positive_all',v>0),('positive_QC',(v>0)&qc),('count_ge2_QC',(v>=2)&qc),('above_background_QC',(v>q95)&(v>0)&qc)]:check(gt[key],int(mask.sum()) if present[j] else None,g+key)
        for record in records.to_dict('records'):
            key=record['model_key'];rt=result.loc[(sid,key)];forbidden=set(record['ligand'].split('_')+record['receptor'].split('_')+plan['family_name_genes'])
            panels={k:[gi[g] for g in v if g not in forbidden and present[gi[g]]] for k,v in plan['panels'].items()}
            sn=record['sender_family'].split('_')[-1];rn=record['receiver_family'].split('_')[-1];sl='T' if record['sender_family'].startswith('T_') else 'B';rl='T' if record['receiver_family'].startswith('T_') else 'B'
            missing=[g for g in record['ligand'].split('_')+record['receptor'].split('_')+[sn,rn] if not present[gi[g]]]
            marker_measurable=not missing and len(panels[sl])>=2 and len(panels[rl])>=2
            author_available=prefix is not None and not (prefix in plan['author_scheme_B_unavailable_prefixes'] and 'B' in [sl,rl])
            check(rt.author_lineage_scheme_available,author_available,key)
            check(rt.missing_required_genes,';'.join(missing) or 'NONE',key)
            pg=record['program_genes'].split(';');pn=sum(present[gi[g]] for g in pg);programme_ok=pn>=4 and pn/len(pg)>=.6
            check(rt.programme_measurable_genes,int(pn),key);check(rt.programme_coverage_ok,programme_ok,key)
            for arm in plan['arms']:
                measurable=(not missing and author_available) if arm=='qc_background_author' else marker_measurable
                rule=cache['counts']>0
                if arm=='qc_count_ge2':rule=cache['counts']>=2
                if 'background' in arm:rule=(cache['counts']>0)&(cache['counts']>q95[:,None])
                rule &= present[None,:]
                keep=np.ones(len(fovs),bool) if arm=='all_detected' else qc
                lineage={k:np.count_nonzero(rule[:,cols],axis=1)>=2 for k,cols in panels.items()}
                endpointSL=np.isin(labels,plan['author_T_labels'] if sl=='T' else plan['author_B_labels']) if arm=='qc_background_author' else lineage[sl]
                endpointRL=np.isin(labels,plan['author_T_labels'] if rl=='T' else plan['author_B_labels']) if arm=='qc_background_author' else lineage[rl]
                baseS=endpointSL&rule[:,gi[sn]]&keep;baseR=endpointRL&rule[:,gi[rn]]&keep
                send=baseS.copy();receive=baseR.copy()
                for g in record['ligand'].split('_'):send &= rule[:,gi[g]]
                for g in record['receptor'].split('_'):receive &= rule[:,gi[g]]
                smix=sum(lineage[k].astype(int) for k in panels if k!=sl)>0;rmix=sum(lineage[k].astype(int) for k in panels if k!=rl)>0
                sm=int((send&smix).sum());rm=int((receive&rmix).sum())
                if arm=='qc_background_pure':send &= ~smix;receive &= ~rmix
                expected={'sender_family_cells':int(baseS.sum()),'receiver_family_cells':int(baseR.sum()),'sender_positive':int(send.sum()),'receiver_positive':int(receive.sum()),'overlap_same_cell':int((send&receive).sum()),'sender_mixed_before_pure':sm,'receiver_mixed_before_pure':rm}
                for col,value in expected.items():check(rt[col+'_'+arm],value if measurable else None,key+arm+col)
                nf=0
                for fov in levels:
                    mask=fovs==fov;ns=int((send&mask).sum());nr=int((receive&mask).sum());no=int((send&receive&mask).sum());ok=measurable and ns>=20 and nr>=20;nf+=ok
                    ft=fovtable.loc[(sid,key,arm,fov)]
                    for col,value in [('sender_positive',ns),('receiver_positive',nr),('overlap_same_cell',no)]:check(ft[col],value if measurable else None,key+arm+fov+col)
                    check(ft.both_endpoints_ge20,ok,key);fov_checks+=1
                check(rt['FOVs_both_endpoints_ge20_'+arm],int(nf) if measurable else None,key)
                check(rt['status_'+arm],'UNMEASURABLE_AUTHOR_LINEAGE' if arm=='qc_background_author' and not author_available else ('UNMEASURABLE_ENDPOINT' if not measurable else ('DESCRIPTIVE_FOV_COVERAGE' if nf else 'INSUFFICIENT_PAIRED_FOV_ENDPOINTS')),key)
                if arm=='qc_background':
                    annotated=~np.isin(labels,['NOT_IN_AUTHOR_QC_OBJECT','AUTHOR_NOT_RETRIEVED'])
                    for endpoint,mask,line in [('sender',send,sl),('receiver',receive,rl)]:
                        aa=agreement.loc[(sid,key,endpoint)];same=np.isin(labels,plan['author_T_labels'] if line=='T' else plan['author_B_labels'])
                        for col,v in [('marker_endpoint_cells',mask.sum()),('author_annotation_matched',(mask&annotated).sum()),('author_same_lineage',(mask&same).sum()),('author_other_lineage',(mask&annotated&~same).sum()),('author_unmapped',(mask&~annotated).sum())]:check(aa[col],int(v) if measurable else None,key+endpoint+col)
                arm_checks+=1
            audited.append(dict(sample_id=sid,model_key=key,marker_endpoint_measurable=marker_measurable,author_lineage_scheme_available=author_available,programme_coverage_ok=programme_ok,independent_recount='PASS'))
        assert digest(inputfile['path'])==inputfile['sha256']
        print(json.dumps(dict(audited_file=sid,cells=len(fovs))),flush=True)
    assert len(audited)==65 and arm_checks==390 and len(fovtable)==fov_checks
    pd.DataFrame(audited).to_csv(out/'independent_all65_checks.tsv',sep='\t',index=False,na_rep='NA')
    payload=dict(schema='spatial_contact_native_preflight_audit_v3',status='PASS_ALL65_NATIVE_VALUES_CONTROLS_QC_BACKGROUND_LINEAGE_OVERLAP_AUTHOR_AND_FOV_RECOUNT',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),assays=5,cells=254458,record_assay_checks=65,arm_checks=arm_checks,FOV_checks=fov_checks,native_selected_gene_values_checked=value_count,negative_control_values_checked=control_count,auditor_sha256=digest(__file__),producer_inventory_sha256=digest(out/'inventory.json'),plan_sha256=digest(out/'predeclared_plan.json'),independence='Separate HDF5 CSR read, whole negative CSV MultiIndex join, order-statistic quantile, author table reindex and all65 vector recount; no producer imports. Shared statistical definitions are intentional, not independent biological validation.',scientific_release=False,spatial_contact_validated=False,patient_independence_validated=False)
    (out/'independent_audit.json').write_text(json.dumps(payload,indent=2)+'\n');print(json.dumps(payload),flush=True)


if __name__=='__main__':main()

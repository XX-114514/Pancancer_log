"""Complete released-source RNA state capture, without clinical effects.

Only source aggregates are saved. All expression, coded material and cell
arrays are in RAM. Existing CM models and authorities remain read-only.
"""
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
import hashlib
import json
import os
import re
import shutil
import time

import numpy as np

import observed_T_marker_evidence_v1 as evidence
import qualify_GSE288199_external_state_measurement_v2 as qualification
import GSE288199_native_stream_v1 as source

GENES=evidence.GENES+['CXCL13']
PARENTS=['dual_regulatory_evidence','costate_regulatory_sensitivity']
QC=['all_released_count_columns','human_RNA_QC500_250_MT20']
DEPTH=['all_capture','low_capture','mid_capture','high_capture']
SORT=['all_source_materials','CD3_CD45_enriched','CD45_enriched']
PHASES=['pre','post']
STATES=['R_TNFRSF9','R_LAG3','R_cycling','E_dysfunction','E_cycling','E_memory']


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def signature(contracts):
    keys=['url','bytes','modified','prefix_SHA256','format','dims','features',
          'feature_compressed_SHA256','feature_uncompressed_SHA256']
    return hashlib.sha256(json.dumps([{k:c[k] for k in keys} for c in contracts],
        sort_keys=True,separators=(',',':')).encode()).hexdigest()


def retry(fn,item):
    for attempt in range(3):
        try:return fn(item),attempt
        except (OSError,EOFError):
            if attempt==2:raise
            time.sleep(2)


def make_contract(item):
    axis=source.feature_contract(item['url'].replace('_matrix.mtx.gz','_features.tsv.gz'),GENES)
    return dict(source.census(item['url'],axis),case=item['case'],phase=item['phase'],sort=item['sort'])


def quantile(values):
    return np.quantile(values,[0,.25,.5,.75,1]).tolist() if len(values)>=5 else None


def measure_one(contract):
    (counts,totals,detected,mito,proof),full_SHA=source.read_full(contract)
    rows=[];background=[]
    for qi,qc in enumerate(QC):
        eligible=np.ones(len(counts),bool) if qi==0 else ((totals[:,0]>=500)&(detected>=250)&(mito<=.2*totals[:,0])&(totals[:,0]>0))
        depth=np.full(len(counts),-1,np.int8)
        if eligible.any():
            cuts=np.unique(np.quantile(totals[eligible,0],[1/3,2/3]))
            depth[eligible]=np.searchsorted(cuts,totals[eligible,0],side='left')
        assert (depth[eligible]>=0).all() and (depth[eligible]<=2).all()
        for parent in PARENTS:
            flags=evidence.vector_evidence(counts[:,:len(evidence.GENES)],eligible,parent)
            np.testing.assert_array_equal(flags,evidence.scalar_evidence(counts[:,:len(evidence.GENES)],eligible,parent))
            np.testing.assert_array_equal(flags[:,1:5].sum(axis=1),eligible.astype(int))
            for di,name in enumerate(DEPTH):
                use=eligible.copy() if di==0 else eligible&(depth==di-1)
                n=int(use.sum());observed=flags[use].sum(axis=0).astype(np.int64)
                # Independent scalar count of the complete Boolean vector.
                scalar=np.zeros(len(evidence.METRICS),np.int64)
                for vector in flags[use]:scalar+=vector.astype(np.int64)
                np.testing.assert_array_equal(observed,scalar)
                marker={}
                for pname in ['regulatory_only','CD8_only']:
                    take=use&flags[:,evidence.METRICS.index(pname)]
                    marker[pname]=(counts[take]>0).sum(axis=0).astype(int).tolist()
                quality={}
                for key,v in [('human_UMI',totals[:,0]),('other_GEX_UMI',totals[:,1]),('antibody_UMI',totals[:,2]),('human_genes',detected)]:
                    quality[key]=float(np.mean(v[use])) if n else None
                rows.append(dict(case=contract['case'],phase=contract['phase'],sort=contract['sort'],QC=qc,parent=parent,
                    depth=name,n=n,metrics=observed.tolist(),parent_marker_positive=marker,quality=quality))
        background.append(dict(QC=qc,cells=int(eligible.sum()),human_UMI=int(totals[eligible,0].sum()),
            other_GEX_UMI=int(totals[eligible,1].sum()),antibody_UMI=int(totals[eligible,2].sum()),
            zero_human_RNA=int((eligible&(totals[:,0]==0)).sum())))
    return dict(rows=rows,background=background,case=contract['case'],phase=contract['phase'],sort=contract['sort'],
                full_SHA=full_SHA,entries=proof['entries'],zeros=proof['explicit_zero_entries'],ordering=proof['ordering'])


def aggregate(all_rows):
    output=[]
    for phase,sort,qc,parent,depth in product(PHASES,SORT,QC,PARENTS,DEPTH):
        rows=[r for r in all_rows if r['phase']==phase and r['QC']==qc and r['parent']==parent and r['depth']==depth and
              (sort=='all_source_materials' or r['sort']==sort)]
        assert len({r['case'] for r in rows})==len(rows)
        ncase=len(rows);total=sum(r['n'] for r in rows)
        metrics=np.asarray([r['metrics'] for r in rows],np.int64).reshape(ncase,len(evidence.METRICS))
        states={}
        for state in STATES:
            p='regulatory_only' if state.startswith('R_') else 'CD8_only'
            j=evidence.METRICS.index(state);pj=evidence.METRICS.index(p)
            numerator=metrics[:,j];denominator=metrics[:,pj];good=denominator>=20
            assert np.all(numerator<=denominator)
            values=numerator[good]/denominator[good]
            states[state]=dict(coded_source_units_parent_atleast20=int(good.sum()),
                parent_cell_observations=int(denominator.sum()) if ncase>=5 else None,
                state_positive_observations=int(numerator.sum()) if ncase>=5 else None,
                coded_unit_balanced_fraction_min_q25_median_q75_max=quantile(values),
                status='DESCRIPTIVE_STATE_CAPTURE' if len(values)>=5 else 'INSUFFICIENT5_SUPPORTED_CODED_UNITS')
        marker={}
        for p in ['regulatory_only','CD8_only']:
            den=metrics[:,evidence.METRICS.index(p)];good=den>=20
            matrix=np.asarray([r['parent_marker_positive'][p] for r in rows],np.int64).reshape(ncase,len(GENES))
            assert np.all(matrix<=den[:,None])
            marker[p]={g:dict(positive_observations=int(matrix[:,j].sum()) if ncase>=5 else None,
                coded_unit_balanced_detection_min_q25_median_q75_max=quantile(matrix[good,j]/den[good])) for j,g in enumerate(GENES)}
        pair_support=int(((metrics[:,evidence.METRICS.index('regulatory_only')]>=20)&
                          (metrics[:,evidence.METRICS.index('CD8_only')]>=20)).sum())
        output.append(dict(source='GSE288199',phase=phase,sort_scope=sort,QC=qc,parent=parent,capture_scope=depth,
            coded_source_units=ncase,coded_source_units_with_cells=sum(r['n']>0 for r in rows),
            cells=total if ncase>=5 else None,
            metric_positive_observations=dict(zip(evidence.METRICS,metrics.sum(axis=0).astype(int).tolist())) if ncase>=5 else None,
            paired_parent_units_both_atleast20=pair_support,reviewed_independent_patients_not_verified=True,
            state_capture=states,parent_marker_capture=marker,
            coded_unit_mean_quality_medians={key:float(np.median([r['quality'][key] for r in rows if r['n']]))
                if sum(r['n']>0 for r in rows)>=5 else None for key in ['human_UMI','other_GEX_UMI','antibody_UMI','human_genes']},
            not_original_CM_activity_or_clinical_or_spatial_test=True))
    assert len(output)==96
    return output


def main():
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert out.parent.parent==Path('/dev/shm') and not out.exists()
    out.mkdir(mode=0o700)
    def save(n,v):(out/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    prior=Path(os.environ['PANCANCER_EXTERNAL_MEASUREMENT_OUTPUT'])
    audited=Path(os.environ['PANCANCER_EXTERNAL_FEATURE_AUDIT_OUTPUT'])
    t=json.loads((prior/'private_terminal.json').read_text());a=json.loads((audited/'private_terminal.json').read_text())
    assert t['status']=='COMPLETE_OFFICIAL_METADATA_AND_ALL_FEATURE_AXIS_QUALIFICATION'
    assert a['status']=='PASS_COMPLETE_INDEPENDENT_MATERIAL_FEATURE_AXIS_AND_NAMESPACE_REPLAY'
    assert a['producer_terminal_SHA256']==sha(prior/'private_terminal.json')
    assert sha(Path(qualification.__file__))==t['source_SHA256']
    assert sha(Path(evidence.__file__))=='e8950f189eb399a4aa647e9e36c96dd4f2d0733404d47aa5d5f72b9d02922ae6'
    for folder,terminal in [(prior,t),(audited,a)]:
        for name,item in terminal['files'].items():assert sha(folder/name)==item['SHA256']
    metadata=json.loads((prior/'private_complete_source_measurement_qualification.json').read_text())
    save('plan_before_census_or_native_expression.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),
        producer_SHA256=sha(Path(__file__)),stream_reader_SHA256=sha(Path(source.__file__)),recipe_SHA256=sha(Path(evidence.__file__)),
        prior_qualification_SHA256=sha(prior/'private_terminal.json'),prior_audit_SHA256=sha(audited/'private_terminal.json'),
        scope='All60 declared source matrices; all pre/post and3 sort scopes,2 QC,2 parents,4 capture scopes=96 source opportunities.',
        genes=GENES,states=STATES,parents=PARENTS,QC=QC,capture=DEPTH,sort=SORT,
        no_outcome_based_selection=True,
        marker_rules='Same22-gene molecule-positive recipes, unchanged. CXCL13 detection described in both parents, no new cell identity.',
        counts='Complete original integer MM counts; unique monotone cell-feature or feature-cell coordinates. Global nnz/gzipCRC/SHA/length/LastModified/prefix and original per-file feature SHA checked; all chunk target/namespace/mito sums independently sparse selected.',
        source_transport='Direct official HTTPS; 64KiB compressed header prefixes (may contain unparsed count bytes), no raw payload cache; full file only after all60 header/budget gates.',
        resource_gates='Total compressed<=16GiB; each<=4GiB;<=1million columns; fixed2 concurrent streams, RAM estimate plus3GiB headroom; no individual arrays written.',
        QC_meaning='All released columns versus human-onlyUMI>=500,genes>=250,mito<=20percent. This is an RNA-only reanalysis, not exact author doublet/TCR/cluster population.',
        depth='Within-library/QC humanUMI value tertiles; ties stay together, no forced equal bins.',
        source_outputs='Coded-case/phase rows onlyRAM; >=5sourceunits for pooled observation summaries, >=20parent observations and>=5supportedunits for pointwise rate quantiles. Below5 distributions suppressed.',
        clinical='Coded source case labels not reviewed independent patients. Treatment/pathology unavailable; no clinical models, state-pair correlations or CM02 projection.',
        unchanged='OriginalCM/W/S/H/authority, prior6recipes and all closed families unchanged. No barcode or VDJ payload or clinical values read.'))
    save('private_synthetic_stream_proof.json',source.synthetic())
    raw=qualification.fetch(qualification.BASE+'soft/GSE288199_family.soft.gz',1024**2)
    assert qualification.digest(raw)==qualification.SOFT_SHA
    samples,series=qualification.soft_rows(raw);gex={}
    for row in samples.values():
        hit=re.fullmatch(r'(HN[0-9]+) (Pre|Post) Gene expression data',row['title'])
        if hit:
            key=(hit[1],hit[2].lower());assert key not in gex
            label=row['traits']['cell type']
            assert label in ['CD45+ CD3+ T cells','CD45+ cells']
            gex[key]='CD3_CD45_enriched' if label=='CD45+ CD3+ T cells' else 'CD45_enriched'
    items=[]
    for value in sorted(series['supplementary_file']):
        hit=re.search(r'/GSE288199_(HN[0-9]+)_(pre|post)_T_matrix.mtx.gz$',value)
        if hit:
            key=hit.groups();assert key in gex
            items.append(dict(case=key[0],phase=key[1],sort=gex[key],url=qualification.BASE+'suppl/'+value.rsplit('/',1)[-1]))
    assert len(items)==60 and len(gex)==60 and {(r['case'],r['phase']) for r in items}==set(gex)
    contracts=[];census_retries=0
    with ThreadPoolExecutor(max_workers=3) as pool:
        for k,(contract,retried) in enumerate(pool.map(lambda item:retry(make_contract,item),items),1):
            contracts.append(contract);census_retries+=retried
            save('private_progress.json',dict(stage='all_matrix_header_census',completed=k,total=60))
    feature_sig=qualification.digest(json.dumps(sorted((c['feature_compressed_SHA256'],c['feature_uncompressed_SHA256'],c['features']) for c in contracts),separators=(',',':')).encode())
    assert feature_sig==metadata['source_feature_manifest_SHA256']
    total_bytes=sum(c['bytes'] for c in contracts);maxcells=max(c['dims'][1] for c in contracts)
    assert total_bytes<=16*1024**3 and max(c['bytes'] for c in contracts)<=4*1024**3
    memory={line.split(':')[0]:line.split(':')[1].strip() for line in Path('/proc/meminfo').read_text().splitlines()}
    available=int(memory['MemAvailable'].split()[0])*1024
    required=2*(maxcells*(len(GENES)*8+8*12+30)+128*1024**2)+len(contracts)*33552*4
    assert available>required+3*1024**3 and shutil.disk_usage('/dev/shm').free>128*1024**2
    budget=dict(completed_UTC=datetime.now(timezone.utc).isoformat(),matrices=len(contracts),
        source_header_manifest_SHA256=signature(contracts),feature_manifest_SHA256=feature_sig,
        compressed_bytes=total_bytes,declared_columns=sum(c['dims'][1] for c in contracts),
        declared_nonzero_entries=sum(c['dims'][2] for c in contracts),
        column_count_min_q25_median_q75_max=np.quantile([c['dims'][1] for c in contracts],[0,.25,.5,.75,1]).tolist(),
        RAM_estimate_bytes=required,mem_available_bytes=available,concurrent_full_streams=2,
        header_census_retries=census_retries,no_cell_state_effects_evaluated_before_budget=True,
        no_full_matrix_payloads_cached=True)
    save('private_complete_header_budget_gate.json',budget)
    results=[];all_rows=[];full_retries=0
    with ThreadPoolExecutor(max_workers=2) as pool:
        for k,(result,retried) in enumerate(pool.map(lambda c:retry(measure_one,c),contracts),1):
            results.append(result);all_rows.extend(result['rows']);full_retries+=retried
            save('private_progress.json',dict(stage='complete_native_matrix_state_capture',completed=k,total=60,
                completed_source_entries=sum(r['entries'] for r in results),observed_UTC=datetime.now(timezone.utc).isoformat()))
    assert len(results)==60 and sum(r['entries'] for r in results)==budget['declared_nonzero_entries']
    summaries=aggregate(all_rows)
    save('private_all96_native_state_capture_opportunities.json',summaries)
    save('private_complete_stream_conservation.json',dict(all60_full_streams_completed=True,
        source_header_manifest_SHA256=signature(contracts),all_compressed_payload_digests_SHA256=qualification.digest(json.dumps(sorted(r['full_SHA'] for r in results)).encode()),
        declared_and_observed_entries_equal=True,observed_entries=sum(r['entries'] for r in results),
        all_raw_chunk_counts_human_other_antibody_mito_and_targets_sparse_replayed=True,
        all_cell_parent_recipes_scalar_replayed=True,all_code_unit_state_totals_scalar_replayed=True,
        explicit_zero_entries=sum(r['zeros'] for r in results),source_coordinate_orders=dict(Counter(r['ordering'] for r in results)),
        full_stream_retries=full_retries))
    save('private_terminal.json',dict(status='COMPLETE_ALL60_NATIVE_RNA_OBSERVED_STATE_CAPTURE_NO_EFFECT_TESTS',
        completed_UTC=datetime.now(timezone.utc).isoformat(),producer_SHA256=sha(Path(__file__)),
        stream_reader_SHA256=sha(Path(source.__file__)),source_opportunities=len(summaries),
        full60_raw_counts_and_namespace_and_marker_numeric_checks=True,
        no_clinical_state_pair_or_spatial_tests=True,no_original_CM_activity_or_new_identity_assigned=True,
        coded_units_not_reviewed_independent_patients=True,no_cell_case_expression_arrays_or_matrices_saved=True,
        no_barcode_VDJ_or_clinical_payloads_read=True,
        files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status='COMPLETE_ALL60_NATIVE_RNA_OBSERVED_STATE_CAPTURE_NO_EFFECT_TESTS',terminal_SHA256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

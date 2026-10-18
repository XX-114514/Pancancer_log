"""Full alternate raw parser, every capture/selection/count and frozen score."""
import gzip
import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import load_npz

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from prepare_NSCLC_full_native_gene_units_v1 import TrackedReader


def equal(a,b):assert np.allclose(a,b,rtol=0,atol=1e-8,equal_nan=True)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    source=root/'NSCLC_immune_RNA_MCP_native_cells_20261019_v1'
    launchpath=root/'HNSCC_receiver_followup_job_control_20261016_v2/measure_NSCLC_immune_RNA_MCP_native_cells_v1/launch_receipt.json'
    launch=json.loads(launchpath.read_text());assert launch['source_sha256']=='2db3770b4b3d29cca228f049f3f050c279641dfba5b84ebf5e8d016e56b7bbcf'
    out=root/'NSCLC_immune_RNA_MCP_native_cells_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_complete_raw_RNA_and_all_score_replay.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,original_producer_launch_sha256=sha(launchpath),
        scope='Entire original2b-coordinate matrix using independent pandas whitespace parser. All raw RNA/depth,34 marker-positive bits/all cells,all gates/case counts/hash selections,all selected physical counts,all6 scores/params/case means/domain descriptions.',
        dependency='Wait same source-bound PID/startticks, no scientific restart. No clinical values or effects read.'),indent=2)+'\n')
    deadline=time.monotonic()+21600
    while not (source/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists(),'PRODUCER_MISSING_WITHOUT_TERMINAL_NO_RESTART'
        st=p.read_text().rsplit(')',1)[1].split();assert st[19]==str(launch['proc_start_ticks']) and st[0]!='Z'
        assert time.monotonic()<deadline
        (out/'private_progress.json').write_text(json.dumps(dict(phase='WAITING_SAME_LIVE_NATIVE_PRODUCER',pid=launch['pid'],startticks=st[19],observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n');time.sleep(15)
    terminal=json.loads((source/'private_terminal.json').read_text());assert terminal['status']=='COMPLETE_ALL_NATIVE_RNA_CAPTURE_AND_FROZEN_FULL_PROGRAMMES_INDEPENDENT_ACCEPTANCE_PENDING'
    for name,item in terminal['files'].items():assert sha(source/name)==item['sha256']
    native=root/'NSCLC_full_native_gene_units_20261016_v1';units=pd.read_csv(native/'private_all_native_units_and_quality.tsv',sep='\t',keep_default_na=False).set_index('unit_index')
    rows=np.load(native/'private_cell_unit_index.npy');RNA=np.load(native/'private_full_RNA.npy');detected=np.load(native/'private_detected_genes.npy')
    libraries=sorted(set(units.sampleID));assert len(libraries)==243
    mi={s:i for i,s in enumerate(libraries)};materials=units.sampleID.map(mi).to_numpy(int)[rows];major=units.major_cell_type.to_numpy(str)[rows]
    nativegenes=pd.read_csv(root/'GSE243013_genes.csv.gz',keep_default_na=False).geneSymbol.str.strip().str.upper().to_numpy(str)
    physical=np.load(source/'private_union_native_columns.npy');colmap=np.full(len(nativegenes),-1,np.int32);colmap[physical]=np.arange(len(physical))
    axis=root/'NSCLC_immune_RNA_MCP_native_axis_20261019_v1';mapping=pd.read_csv(axis/'private_all_active_gene_native_axis.tsv',sep='\t')
    markers=pd.read_csv(axis/'private_fixed34_capture_gene_axis.tsv',sep='\t');markerindex=np.full(len(nativegenes),-1,np.int32);markerindex[markers.native_column]=np.arange(len(markers))
    md=pd.read_csv(source/'private_selected_native_cells.tsv.gz',sep='\t');assert md.native_row.is_unique and md.native_row.is_monotonic_increasing
    selectedmap=np.full(len(rows),-1,np.int32);selectedmap[md.native_row]=np.arange(len(md))
    raw_selected=np.zeros((len(md),len(physical)),np.int32);markerbits=np.zeros((len(rows),34),bool)
    weights=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv',sep='\t')
    weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')]
    views=['R_FOXP3_capture','E_effector_capture','Myeloid'];slots=[(v,p) for v in views for p in ['MCP1','MCP2']]
    coefficients=np.zeros((len(nativegenes),6));centers=np.zeros(6)
    for j,(view,program) in enumerate(slots):
        w=weights[weights.view.eq(view)&weights.program.eq(program)]
        cols=mapping[mapping.view.eq(view)&mapping.program.eq(program)].set_index('gene').native_column.reindex(w.gene).to_numpy(int)
        assert len(cols)==len(set(cols))
        coefficients[cols,j]=w.weight.to_numpy()/w.discovery_sd.to_numpy();centers[j]=float(np.sum(coefficients[cols,j]*w.discovery_mean))
    active=np.any(coefficients!=0,axis=1);stream_sums=np.zeros((len(rows),6));totals=np.zeros(len(rows),np.int64);n_genes=np.zeros(len(rows),np.int32)
    download=json.loads((root/'count_matrix_download_v1/download_state.json').read_text());count=0;last=-1;checkpoint=-1;started=time.monotonic()
    with Path(download['path']).open('rb') as binary:
        tracked=TrackedReader(binary)
        with gzip.GzipFile(fileobj=tracked,mode='rb') as stream:
            assert stream.readline().startswith(b'%%MatrixMarket matrix coordinate real general')
            line=stream.readline()
            while line.startswith(b'%'):line=stream.readline()
            assert tuple(map(int,line.split()))==(1254749,31831,2010550708)
            for block in pd.read_csv(stream,sep=r'\s+',header=None,names=['r','c','v'],dtype={'r':np.int64,'c':np.int64,'v':np.float64},chunksize=750000,engine='c'):
                r=block.r.to_numpy()-1;c=block.c.to_numpy()-1;v=block.v.to_numpy()
                assert np.isfinite(v).all() and (v>0).all() and np.array_equal(v,np.floor(v)) and v.max()<2**31
                assert r.min()>=0 and r.max()<len(rows) and c.min()>=0 and c.max()<len(nativegenes)
                coordinate=r*len(nativegenes)+c;assert coordinate[0]>last and np.all(coordinate[1:]>coordinate[:-1]);last=int(coordinate[-1])
                lo=int(r[0]);hi=int(r[-1])+1
                totals[lo:hi]+=np.bincount(r-lo,weights=v,minlength=hi-lo).astype(np.int64)
                n_genes[lo:hi]+=np.bincount(r-lo,minlength=hi-lo).astype(np.int32)
                keep=(selectedmap[r]>=0)&(colmap[c]>=0);raw_selected[selectedmap[r[keep]],colmap[c[keep]]]=v[keep].astype(np.int32)
                mark=markerindex[c]>=0;markerbits[r[mark],markerindex[c[mark]]]=True
                used=active[c];rr=r[used];cc=c[used];expression=np.log1p(v[used]*10000/RNA[rr])
                for j in range(6):stream_sums[lo:hi,j]+=np.bincount(rr-lo,weights=expression*coefficients[cc,j],minlength=hi-lo)
                count+=len(block)
                if count//25000000!=checkpoint:
                    checkpoint=count//25000000
                    (out/'private_progress.json').write_text(json.dumps(dict(phase='INDEPENDENT_FULL_RAW_PANDAS_REPLAY',entries=count,expected_entries=2010550708,elapsed_seconds=time.monotonic()-started,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    assert count==2010550708 and tracked.bytes==download['bytes'] and tracked.digest.hexdigest()==download['sha256']
    assert np.array_equal(totals,RNA) and np.array_equal(n_genes,detected)
    saved_raw=load_npz(source/'private_selected_raw_union_counts.npz');assert saved_raw.shape==raw_selected.shape
    for start in range(0,len(md),1000):assert np.array_equal(saved_raw[start:start+1000].toarray(),raw_selected[start:start+1000])
    fixed=root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json';panels=json.loads(fixed.read_text())['panels']
    marker_names=markers.gene.tolist()
    core=np.count_nonzero(markerbits[:,[marker_names.index(g) for g in panels['pan_T']]],axis=1)>=2
    reg=markerbits[:,marker_names.index('FOXP3')]&(np.count_nonzero(markerbits[:,[marker_names.index(g) for g in panels['regulatory_associated']]],axis=1)>=2)
    eff=np.count_nonzero(markerbits[:,[marker_names.index(g) for g in ['GZMB','PRF1','CTSW','NKG7']]],axis=1)>=2
    pools={'R_FOXP3_capture':(major=='T/NK cell')&core&reg&~eff,'E_effector_capture':(major=='T/NK cell')&core&eff&~reg,'Myeloid':major=='Myeloid cell'}
    coverage=pd.read_csv(source/'private_all729_material_view_capture.tsv',sep='\t').set_index(['view','material_index'])
    assert len(coverage)==729
    full_material=np.bincount(materials,minlength=243);expected_selection=[]
    for view,mask in pools.items():
        for case,name in enumerate(libraries):
            members=np.flatnonzero(mask&(materials==case));r=coverage.loc[(view,case)]
            assert r.sampleID==name and r.full_view_cells==len(members) and r.all_native_material_cells==full_material[case] and bool(r.eligible)==bool(len(members)>=20)
            if len(members)<20:continue
            chosen=sorted(members,key=lambda r:hashlib.sha256(('2026101802|'+view+'|'+str(int(r))).encode()).hexdigest())[:100]
            expected_selection.extend((int(r),view,case,name,len(members),int(RNA[r]),int(detected[r])) for r in chosen)
    expected_selection.sort(key=lambda r:r[0])
    assert expected_selection==list(md[['native_row','view','material_index','sampleID','full_view_cells','fullRNA','full_detected_genes']].itertuples(index=False,name=None))
    cells=pd.read_csv(source/'private_all_selected_cell_programme_scores.tsv.gz',sep='\t').set_index(['native_row','view','program'])
    cases=pd.read_csv(source/'private_all_material_programme_scores.tsv',sep='\t').set_index(['view','program','material_index'])
    params=pd.read_csv(source/'private_original_discovery_nuisance_coefficients.tsv',sep='\t').set_index(['view','program','term'])
    domain=pd.read_csv(source/'private_all12_training_domain_descriptions.tsv',sep='\t').set_index(['view','covariate'])
    inputroot=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_inputs';maxerror=0.;case_count=0
    for j,(view,program) in enumerate(slots):
        vi=np.flatnonzero(pools[view]);counts=np.bincount(materials[vi],minlength=243)
        rawq=np.stack([np.log1p(RNA[vi]),np.log1p(detected[vi]),counts[materials[vi]]/full_material[materials[vi]],np.log1p(counts[materials[vi]])],axis=1)
        frozen=np.load(inputroot/(view+'_frozen_projection.npz'));q=(rawq-frozen['quality_mean'])/frozen['quality_sd']
        tm=pd.read_csv(inputroot/(view+'_discovery_metadata.tsv'),sep='\t');D=np.c_[np.ones(len(tm)),tm[['cellQ','geneQ','sourceCC','fraction','logN']].to_numpy(float)]
        names=Path(str(inputroot/(view+'_discovery'))+'_genes.tsv').read_text().splitlines()
        with gzip.open(inputroot/(view+'_discovery_tpm.mtx.gz'),'rb') as f:tr=mmread(f).tocsr()
        w=weights[weights.view.eq(view)&weights.program.eq(program)];where=pd.Index(names).get_indexer(w.gene);beta=w.weight.to_numpy()/w.discovery_sd.to_numpy()
        rawtrain=np.asarray(tr[where].T.dot(beta)).ravel()-centers[j];Q,R=np.linalg.qr(D,mode='reduced');coef=np.linalg.solve(R,Q.T@rawtrain)
        for term,value in zip(['intercept','cellQ','geneQ','sourceCC','fraction','logN'],coef):equal(params.loc[(view,program,term),'coefficient'],value)
        score=stream_sums[vi,j]-centers[j]-np.c_[np.ones(len(vi)),q]@coef[[0,1,2,4,5]]
        mine=md[md.view.eq(view)];ix=pd.Index(vi).get_indexer(mine.native_row)
        selected=cells.loc[pd.MultiIndex.from_arrays([mine.native_row,np.repeat(view,len(mine)),np.repeat(program,len(mine))])]
        error=float(np.max(np.abs(selected.score.to_numpy()-score[ix])));assert error<1e-8;maxerror=max(maxerror,error)
        for k,field in enumerate(['cellQ','geneQ','fraction','logN']):
            equal(selected[field],q[ix,k]);r=domain.loc[(view,field)];low=float(tm[field].min());high=float(tm[field].max())
            assert r.captured_cells==len(vi) and r.cells_outside_training==np.count_nonzero((q[:,k]<low)|(q[:,k]>high));equal(r.training_min,low);equal(r.training_max,high);equal(r.new_min,q[:,k].min());equal(r.new_max,q[:,k].max())
        for case,name in enumerate(libraries):
            keep=mine.material_index.eq(case).to_numpy()
            if not keep.any():continue
            r=cases.loc[(view,program,case)];alltake=materials[vi]==case
            equal(r.score,score[ix[keep]].mean());equal(r.all_captured_mean_sensitivity,score[alltake].mean())
            assert r.sampleID==name and r.selected_cells==keep.sum() and r.full_view_cells==alltake.sum() and r.all_native_material_cells==full_material[case]
            for k,field in enumerate(['cellQ','geneQ','fraction','logN']):equal(r[field],q[ix[keep],k].mean())
            case_count+=1
    assert len(cells)==2*len(md) and case_count==len(cases)
    t=dict(status='ACCEPTED_ALL_RAW_NATIVE_COUNTS_RNA_GATES_HASH_SELECTIONS_AND_FULL_FROZEN_PROGRAMMES',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),
        raw_coordinates_replayed=count,all_native_cells=len(rows),selected_cells=len(md),selected_physical_count_slots_checked=int(raw_selected.size),
        cell_scores_checked=len(cells),case_scores_checked=len(cases),capture_positions_checked=729,domain_positions_checked=12,
        max_score_error=maxerror,no_clinical_effects=True,new_download=False)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n')
    print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

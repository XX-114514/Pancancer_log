"""One raw pass, exact RNA gates and all frozen programme scores; no outcomes."""
import gzip
import hashlib
import io
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import coo_matrix, save_npz

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from prepare_NSCLC_full_native_gene_units_v1 import TrackedReader, chunks


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);cap=96*1024**2;guard=resource(root,cap)
    axis=root/'NSCLC_immune_RNA_MCP_native_axis_20261019_v1'
    acceptance=root/'NSCLC_immune_RNA_MCP_native_axis_audit_20261019_v1/private_terminal.json'
    assert sha(acceptance)=='575bb2a84376ddd6557882ca47f0f8045d2be7c5d08854082482b0da099759e4'
    assert sha(axis/'private_terminal.json')=='7f1ddaa54cfb2fd550c4892ce6e63363888641103bf57f6a986f974314a33e5d'
    for name,item in json.loads((axis/'private_terminal.json').read_text())['files'].items():assert sha(axis/name)==item['sha256']
    frozen_inputs=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1'
    assert sha(frozen_inputs/'private_terminal.json')=='2474355a003d96232b752ad2c77ba92a506f4c68829c7e0adf434045cf03d656'
    for name,item in json.loads((frozen_inputs/'private_terminal.json').read_text())['files'].items():assert sha(frozen_inputs/name)==item['sha256']
    assert sha(root/'GSE243013_genes.csv.gz')=='ac325c9ce6bfaf2270eb931e7d947bf31d1fb151b14cb166ed9f6c8f4b34a104'
    assert sha(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv')=='fd7316c67beb687c403a752c36b5989644e432f59a83119b21eade5efb0d89a7'
    native=root/'NSCLC_full_native_gene_units_20261016_v1'
    assert sha(native/'private_terminal.json')=='8e8d81678227a506842315b1b2a5e28b604c0f350e475c62d589ea7e0d9719bb'
    oldaudit=root/'NSCLC_full_native_gene_unit_audit_20261016_v1/private_terminal.json'
    assert sha(oldaudit)=='bb20365ee83758a7df8478ceadede2596fe749928dab5765e3bc0866d72dd8c1'
    outputs=json.loads((native/'private_terminal.json').read_text())['outputs']
    filenames=['private_cell_unit_index.npy','private_full_RNA.npy','private_detected_genes.npy','private_all_native_units_and_quality.tsv']
    for name in filenames:assert sha(native/name)==outputs[name]['sha256']
    rows=np.load(native/'private_cell_unit_index.npy');RNA=np.load(native/'private_full_RNA.npy');detected=np.load(native/'private_detected_genes.npy')
    units=pd.read_csv(native/'private_all_native_units_and_quality.tsv',sep='\t',keep_default_na=False).set_index('unit_index')
    assert len(rows)==len(RNA)==len(detected)==1254749 and len(units)==11701
    libraries=sorted(set(units.sampleID));assert len(libraries)==243
    lib_index={name:i for i,name in enumerate(libraries)}
    materials=units.sampleID.map(lib_index).to_numpy(int)[rows];major=units.major_cell_type.to_numpy(str)[rows]
    full_material=np.bincount(materials,minlength=len(libraries))
    gene=pd.read_csv(root/'GSE243013_genes.csv.gz',keep_default_na=False).geneSymbol.str.strip().str.upper().to_numpy(str)
    mapping=pd.read_csv(axis/'private_all_active_gene_native_axis.tsv',sep='\t')
    markers=pd.read_csv(axis/'private_fixed34_capture_gene_axis.tsv',sep='\t')
    assert mapping.native_column.ge(0).all() and markers.native_column.ge(0).all()
    physical=np.array(sorted(set(mapping.native_column)|set(markers.native_column)),int)
    colmap=np.full(len(gene),-1,np.int32);colmap[physical]=np.arange(len(physical),dtype=np.int32)
    mem=int(next(s.split()[1] for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith('MemAvailable:')))*1024
    assert mem>1254749*len(physical)*40+4*1024**3
    download=json.loads((root/'count_matrix_download_v1/download_state.json').read_text())
    assert download['status']=='COMPLETE_DOWNLOAD_GZIP_INTEGRITY_AND_SHA256'
    matrix=Path(download['path']);assert matrix.stat().st_size==download['bytes']
    out=root/'NSCLC_immune_RNA_MCP_native_cells_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    fixed=root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    panels=json.loads(fixed.read_text())['panels']
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,
        full_axis_acceptance_sha256=sha(acceptance),native_raw_audit_sha256=sha(oldaudit),matrix_sha256=download['sha256'],matrix_bytes=download['bytes'],
        scope='All1254749 native cells/243 materials, all6 full programme axes and34 capture genes. No outcomes or new association tests.',
        capture='Original author T/NK pool plus exact frozen panT/R-only/E-only RNA gates; author Myeloid cell whole pool. Context not identical to original atlas nodes; no new formal identity.',
        sampling='Per material/view >=20, up to100 smallest SHA256(2026101802|view|native_row). Independent of RNA value within gated view/outcome; no old208 restriction.',
        measurement='Full native RNA denominator; complete frozen weights/gene mean-SD and discovery-only nuisance coefficients. Original quality standardization, unknown source constant omitted; no NSCLC refit.',
        fraction_denominator='All native immune cells per material, not whole tumour TME. Explicit transport context and training-domain diagnostics.',
        all_cell_score_mean='Measurement sensitivity only; selected100 mean remains anchor. No clinical interpretation before independent raw/count/score acceptance.',
        union_native_columns=len(physical),memory_available_bytes=mem,no_new_download=True)
    (out/'plan_before_all_native_RNA_measurement.json').write_text(json.dumps(plan,indent=2)+'\n')
    started=time.monotonic();parts=[];totals=np.zeros(len(rows),np.int64);n_genes=np.zeros(len(rows),np.int32);entries=0;last=-1;checkpoint=-1
    with matrix.open('rb') as binary:
        tracked=TrackedReader(binary)
        with gzip.GzipFile(fileobj=tracked,mode='rb') as stream:
            assert stream.readline().startswith(b'%%MatrixMarket matrix coordinate real general')
            line=stream.readline()
            while line.startswith(b'%'):line=stream.readline()
            dims=tuple(map(int,line.split()));assert dims==(len(rows),len(gene),2010550708)
            for a in chunks(stream):
                assert np.isfinite(a).all() and np.array_equal(a,np.floor(a)) and (a[:,2]>0).all() and a[:,2].max()<2**31
                row=a[:,0].astype(np.int64)-1;col=a[:,1].astype(np.int64)-1;value=a[:,2].astype(np.int32)
                assert row.min()>=0 and row.max()<len(rows) and col.min()>=0 and col.max()<len(gene)
                coordinate=row*len(gene)+col;assert coordinate[0]>last and (np.diff(coordinate)>0).all();last=int(coordinate[-1])
                lo=int(row[0]);hi=int(row[-1])+1
                totals[lo:hi]+=np.bincount(row-lo,weights=value,minlength=hi-lo).astype(np.int64)
                n_genes[lo:hi]+=np.bincount(row-lo,minlength=hi-lo).astype(np.int32)
                keep=colmap[col]>=0
                parts.append((row[keep].astype(np.int32),colmap[col[keep]],value[keep]))
                entries+=len(a)
                if entries//25000000!=checkpoint:
                    checkpoint=entries//25000000
                    (out/'private_progress.json').write_text(json.dumps(dict(phase='READING_EXISTING_NATIVE_MATRIX',entries=entries,expected_entries=dims[2],elapsed_seconds=time.monotonic()-started,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    assert entries==dims[2] and tracked.bytes==download['bytes'] and tracked.digest.hexdigest()==download['sha256']
    assert np.array_equal(totals,RNA) and np.array_equal(n_genes,detected)
    rr=np.concatenate([a[0] for a in parts]);cc=np.concatenate([a[1] for a in parts]);vv=np.concatenate([a[2] for a in parts]);del parts
    X=coo_matrix((vv,(rr,cc)),shape=(len(rows),len(physical))).tocsr();del rr,cc,vv
    assert X.has_canonical_format and X.nnz>0
    positions=markers.set_index('gene').native_column.map(lambda c:int(colmap[c])).to_dict()
    def gate(genes,minimum):return np.asarray((X[:,[positions[g] for g in genes]]>0).sum(axis=1)).ravel()>=minimum
    core=gate(panels['pan_T'],2);reg=gate(panels['regulatory_associated'],2)&gate(['FOXP3'],1);eff=gate(['GZMB','PRF1','CTSW','NKG7'],2)
    pools={'R_FOXP3_capture':(major=='T/NK cell')&core&reg&~eff,'E_effector_capture':(major=='T/NK cell')&core&eff&~reg,'Myeloid':major=='Myeloid cell'}
    selection=[];coverage=[]
    for view,mask in pools.items():
        for case,name in enumerate(libraries):
            take=np.flatnonzero(mask&(materials==case));n=len(take)
            coverage.append(dict(dataset='GSE243013',sampleID=name,material_index=case,view=view,full_view_cells=n,all_native_material_cells=int(full_material[case]),eligible=n>=20))
            if n<20:continue
            chosen=sorted(take,key=lambda row:hashlib.sha256(('2026101802|'+view+'|'+str(int(row))).encode()).hexdigest())[:100]
            selection.extend(dict(native_row=int(row),view=view,material_index=case,sampleID=name,full_view_cells=n,fullRNA=int(RNA[row]),full_detected_genes=int(detected[row])) for row in chosen)
    md=pd.DataFrame(selection).sort_values('native_row').reset_index(drop=True);assert md.native_row.is_unique
    cov=pd.DataFrame(coverage);cov.to_csv(out/'private_all729_material_view_capture.tsv',sep='\t',index=False)
    selected=X[md.native_row.to_numpy(int)].copy()
    encoded=io.BytesIO();save_npz(encoded,selected,compressed=True);assert encoded.tell()<cap-16*1024**2
    (out/'private_selected_raw_union_counts.npz').write_bytes(encoded.getvalue());del encoded
    np.save(out/'private_union_native_columns.npy',physical,allow_pickle=False)
    logged=X.astype(np.float64);logged.data=np.log1p(logged.data*10000/np.repeat(RNA,np.diff(logged.indptr)))
    inputroot=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_inputs'
    weights=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv',sep='\t');weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')]
    all_scores=[];case_scores=[];parameters=[];domain=[]
    for view,mask in pools.items():
        vi=np.flatnonzero(mask);material_count=np.bincount(materials[vi],minlength=len(libraries))
        rawquality=np.c_[np.log1p(RNA[vi]),np.log1p(detected[vi]),material_count[materials[vi]]/full_material[materials[vi]],np.log1p(material_count[materials[vi]])]
        frozen=np.load(inputroot/(view+'_frozen_projection.npz'));quality=(rawquality-frozen['quality_mean'])/frozen['quality_sd']
        trainmd=pd.read_csv(inputroot/(view+'_discovery_metadata.tsv'),sep='\t');D=np.c_[np.ones(len(trainmd)),trainmd[['cellQ','geneQ','sourceCC','fraction','logN']].to_numpy(float)]
        gene_list=Path(str(inputroot/(view+'_discovery'))+'_genes.tsv').read_text().splitlines()
        with gzip.open(inputroot/(view+'_discovery_tpm.mtx.gz'),'rb') as f:tr=mmread(f).tocsr()
        for j,field in enumerate(['cellQ','geneQ','fraction','logN']):
            low=float(trainmd[field].min());high=float(trainmd[field].max());outside=(quality[:,j]<low)|(quality[:,j]>high)
            domain.append(dict(dataset='GSE243013',view=view,covariate=field,captured_cells=len(vi),training_min=low,training_max=high,new_min=float(quality[:,j].min()),new_max=float(quality[:,j].max()),cells_outside_training=int(outside.sum())))
        mine=md[md.view.eq(view)].copy();within=pd.Index(vi).get_indexer(mine.native_row);assert (within>=0).all()
        for program in ['MCP1','MCP2']:
            w=weights[weights.view.eq(view)&weights.program.eq(program)].copy();assert len(w)>0
            cols=mapping[mapping.view.eq(view)&mapping.program.eq(program)].set_index('gene').native_column.reindex(w.gene).to_numpy(int)
            beta=w.weight.to_numpy()/w.discovery_sd.to_numpy();center=float(beta@w.discovery_mean.to_numpy())
            traincols=pd.Index(gene_list).get_indexer(w.gene);assert (traincols>=0).all()
            rawtrain=np.asarray(tr[traincols].T.dot(beta)).ravel()-center;coef=np.linalg.lstsq(D,rawtrain,rcond=None)[0]
            rawscore=np.asarray(logged[vi][:,colmap[cols]].dot(beta)).ravel()-center
            score=rawscore-np.c_[np.ones(len(vi)),quality]@coef[[0,1,2,4,5]];assert np.isfinite(score).all()
            for term,value in zip(['intercept','cellQ','geneQ','sourceCC','fraction','logN'],coef):parameters.append(dict(view=view,program=program,term=term,coefficient=float(value)))
            cell=mine.assign(program=program,score=score[within],source_offset_unknown=True)
            for j,field in enumerate(['cellQ','geneQ','fraction','logN']):cell[field]=quality[within,j]
            all_scores.append(cell)
            for case,name in enumerate(libraries):
                g=cell[cell.material_index.eq(case)];alltake=materials[vi]==case
                if not len(g):continue
                case_scores.append(dict(dataset='GSE243013',sampleID=name,material_index=case,view=view,program=program,selected_cells=len(g),full_view_cells=int(material_count[case]),all_native_material_cells=int(full_material[case]),
                    score=float(g.score.mean()),all_captured_mean_sensitivity=float(score[alltake].mean()),source_offset_unknown=True,
                    **{field:float(g[field].mean()) for field in ['cellQ','geneQ','fraction','logN']}))
    md.to_csv(out/'private_selected_native_cells.tsv.gz',sep='\t',index=False)
    pd.concat(all_scores,ignore_index=True).to_csv(out/'private_all_selected_cell_programme_scores.tsv.gz',sep='\t',index=False)
    pd.DataFrame(case_scores).to_csv(out/'private_all_material_programme_scores.tsv',sep='\t',index=False)
    pd.DataFrame(parameters).to_csv(out/'private_original_discovery_nuisance_coefficients.tsv',sep='\t',index=False)
    pd.DataFrame(domain).to_csv(out/'private_all12_training_domain_descriptions.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()};assert sum(v['bytes'] for v in files.values())<cap
    terminal=dict(status='COMPLETE_ALL_NATIVE_RNA_CAPTURE_AND_FROZEN_FULL_PROGRAMMES_INDEPENDENT_ACCEPTANCE_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,raw_stream_entries=entries,
        all_cell_RNA_and_detected_exact=True,matrix_SHA_and_gzip_CRC_pass=True,selected_cells=len(md),no_clinical_effects=True)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

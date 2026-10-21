"""Read original native NSCLC counts; return case arrays through RAM pipes only.

The caller injects the pinned original TrackedReader/chunks and fixed evidence
module. No source writes, learned parameters or selected programme cells.
"""
from pathlib import Path
from datetime import datetime,timezone
from decimal import Decimal
import gzip
import hashlib
import io
import json
import resource
import numpy as np
import pandas as pd

COMMON=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary']
PARENTS=['dual_regulatory_evidence','costate_regulatory_sensitivity']
QC=['native_QC500_250_MT20','all_retained_native']


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for b in iter(lambda:stream.read(16*1024**2),b''):h.update(b)
    return h.hexdigest()


def native_eligibility(total,detected,mitochondrial):
    return (total>=500)&(detected>=250)&(mitochondrial<=.2*total)&(total>0)


def row_counts(flags,materials,n):
    result=np.column_stack([np.bincount(materials,weights=flags[:,j],minlength=n) for j in range(flags.shape[1])]).astype(np.int64)
    other=pd.DataFrame(flags).assign(material=materials).groupby('material',sort=True).sum().reindex(range(n),fill_value=0).to_numpy(np.int64)
    np.testing.assert_array_equal(result,other)
    return result


def group_quality(mask,materials,total,detected,denominator,n):
    counts=np.bincount(materials[mask],minlength=n)
    sums=np.column_stack([np.bincount(materials[mask],weights=np.log1p(x[mask]),minlength=n) for x in [total,detected]])
    means=sums/np.where(counts>0,counts,1)[:,None]
    grouped=pd.DataFrame(dict(material=materials[mask],RNA=np.log1p(total[mask]),genes=np.log1p(detected[mask]))).groupby('material').mean().reindex(range(n))
    np.testing.assert_allclose(means[counts>0],grouped.to_numpy()[counts>0],atol=1e-10,rtol=0)
    answer=[]
    for i in range(n):
        answer.append([float(means[i,0]) if counts[i] else None,float(means[i,1]) if counts[i] else None,
            float(np.log1p(counts[i])),float(counts[i]/denominator[i]) if denominator[i] else None])
    return counts,answer


def prefix_audit(matrix,axis,columns,expected,total,detected,rows=8192):
    got=np.zeros((rows,len(columns)),np.int32);rr=np.zeros(rows,np.int64);dd=np.zeros(rows,np.int32)
    column_to_slot={int(c):i for i,c in enumerate(columns)}
    with gzip.open(matrix,'rt',encoding='ascii') as stream:
        assert stream.readline().startswith('%%MatrixMarket matrix coordinate real general')
        line=stream.readline()
        while line.startswith('%'):line=stream.readline()
        for line in stream:
            tokens=line.split();a,b=map(int,tokens[:2]);v=Decimal(tokens[2])
            assert v.is_finite() and v==v.to_integral_value() and v>0
            c=int(v);a-=1;b-=1
            if a>=rows:break
            rr[a]+=c;dd[a]+=1
            if b in column_to_slot:got[a,column_to_slot[b]]=c
    np.testing.assert_array_equal(got,expected[:rows]);np.testing.assert_array_equal(rr,total[:rows]);np.testing.assert_array_equal(dd,detected[:rows])
    return dict(independent_integer_line_prefix_cells=rows,all_prefix_markers_and_total_and_detection_exact=True)


def read_source(root,progress):
    root=Path(root);native=root/'NSCLC_full_native_gene_units_20261016_v1'
    progress(dict(stage='VERIFY_EXISTING_NATIVE_ANCHORS',started_UTC=datetime.now(timezone.utc).isoformat()))
    assert sha(native/'private_terminal.json')=='8e8d81678227a506842315b1b2a5e28b604c0f350e475c62d589ea7e0d9719bb'
    assert sha(root/'NSCLC_full_native_gene_unit_audit_20261016_v1/private_terminal.json')=='bb20365ee83758a7df8478ceadede2596fe749928dab5765e3bc0866d72dd8c1'
    terminal=json.loads((native/'private_terminal.json').read_text())
    memory=int(next(s.split()[1] for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith('MemAvailable:')))*1024
    assert memory>=9*1024**3,'NATIVE_RAM_HEADROOM_NOT_ACCEPTED'
    names=['private_cell_unit_index.npy','private_full_RNA.npy','private_detected_genes.npy','private_all_native_units_and_quality.tsv','sum_raw_counts.npy','positive_cell_counts.npy']
    for name in names:
        progress(dict(stage='VERIFY_CURRENT_RAW_CACHE_SHA',file=name))
        assert sha(native/name)==terminal['outputs'][name]['sha256'],'RAW_NATIVE_CACHE_SHA_CHANGED'
    gene_path=root/'GSE243013_genes.csv.gz';assert sha(gene_path)=='ac325c9ce6bfaf2270eb931e7d947bf31d1fb151b14cb166ed9f6c8f4b34a104'
    genes=pd.read_csv(gene_path,keep_default_na=False).geneSymbol.str.strip().str.upper().to_numpy(str)
    assert len(genes)==31831 and len(set(genes))==len(genes)
    columns=pd.Index(genes).get_indexer(evidence.GENES);assert (columns>=0).all()
    mt=np.flatnonzero(np.char.startswith(genes,'MT-'));assert len(mt)==13
    physical=np.r_[columns,mt];assert len(set(physical))==len(physical)
    unit=np.load(native/'private_cell_unit_index.npy',allow_pickle=False)
    knownRNA=np.load(native/'private_full_RNA.npy',allow_pickle=False);knownGenes=np.load(native/'private_detected_genes.npy',allow_pickle=False)
    units=pd.read_csv(native/'private_all_native_units_and_quality.tsv',sep='\t',keep_default_na=False).sort_values('unit_index')
    assert len(unit)==len(knownRNA)==len(knownGenes)==1254749 and units.unit_index.tolist()==list(range(11701))
    metadata_path=root/'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz'
    assert sha(metadata_path)=='24155466454df53d01eb6cbbfcd56abeb45adb68f74e721af22a9dce41040b2f'
    md=pd.read_csv(metadata_path,usecols=['sampleID','sub_cell_type','major_cell_type','total_counts','total_counts_mt'],
        dtype={'sampleID':str,'sub_cell_type':str,'major_cell_type':str},keep_default_na=False)
    assert np.array_equal(md.sampleID,units.sampleID.to_numpy()[unit])
    assert np.array_equal(md.sub_cell_type,units.sub_cell_type.to_numpy()[unit])
    assert np.array_equal(md.major_cell_type,units.major_cell_type.to_numpy()[unit])
    np.testing.assert_array_equal(knownRNA,pd.to_numeric(md.total_counts).to_numpy())
    lib=sorted(set(units.sampleID));assert len(lib)==243
    materials=pd.Index(lib).get_indexer(md.sampleID);assert (materials>=0).all()
    clinical_path=root/'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1/private_all243_current_RNA_clinical_scope.tsv'
    assert sha(clinical_path)=='318f0bbf03d48872625bf0e8b2639070d85229cdade3b78498a27ee5cfaed377'
    clinical=pd.read_csv(clinical_path,sep='\t',dtype=str,keep_default_na=False).set_index('sampleID',drop=False)
    assert len(clinical)==243 and clinical.index.is_unique and set(clinical.index)==set(lib)
    clinical=clinical.loc[lib]
    cr=clinical[COMMON+['response']].apply(pd.to_numeric,errors='coerce')
    eligible=clinical.strict_treatment_role.eq('EXPLICIT_ICI_AUTHOR_CASE').to_numpy()&np.isfinite(cr.to_numpy()).all(axis=1)
    assert np.array_equal(eligible,clinical.prospective_RNA_clinical_eligible.map({'True':True,'False':False}))
    receipt=root/'count_matrix_download_v1/download_state.json'
    assert sha(receipt)=='887e3c61396c899fbcc1643a5a149c8ff845d7e9f58cef15405157f10208705b'
    download=json.loads(receipt.read_text());assert download['status']=='COMPLETE_DOWNLOAD_GZIP_INTEGRITY_AND_SHA256'
    matrix=Path(download['path']);assert matrix.is_file() and matrix.stat().st_size==download['bytes']
    counts=np.zeros((len(unit),len(physical)),np.int32);totals=np.zeros(len(unit),np.int64);detected=np.zeros(len(unit),np.int32)
    cmap=np.full(len(genes),-1,np.int16);cmap[physical]=np.arange(len(physical))
    entries=0;last=-1;checkpoint=-1
    with matrix.open('rb') as binary:
        tracked=TrackedReader(binary)
        with gzip.GzipFile(fileobj=tracked,mode='rb') as stream:
            assert stream.readline().startswith(b'%%MatrixMarket matrix coordinate real general')
            line=stream.readline()
            while line.startswith(b'%'):line=stream.readline()
            dims=tuple(map(int,line.split()));assert dims==(len(unit),len(genes),2010550708)
            for a in chunks(stream):
                assert np.isfinite(a).all() and np.array_equal(a,np.floor(a)) and (a[:,2]>0).all() and a[:,2].max()<2**31,'INVALID_NATIVE_INTEGER_COORDINATE'
                row=a[:,0].astype(np.int64)-1;col=a[:,1].astype(np.int64)-1;value=a[:,2].astype(np.int32)
                assert row.min()>=0 and row.max()<len(unit) and col.min()>=0 and col.max()<len(genes)
                coordinate=row*len(genes)+col;assert coordinate[0]>last and (np.diff(coordinate)>0).all(),'DUPLICATE_OR_UNSORTED_NATIVE_COORDINATE'
                last=int(coordinate[-1]);lo=int(row[0]);hi=int(row[-1])+1
                totals[lo:hi]+=np.bincount(row-lo,weights=value,minlength=hi-lo).astype(np.int64)
                detected[lo:hi]+=np.bincount(row-lo,minlength=hi-lo).astype(np.int32)
                take=cmap[col]>=0;counts[row[take],cmap[col[take]]]=value[take]
                entries+=len(a)
                if entries//25000000!=checkpoint:
                    checkpoint=entries//25000000;progress(dict(stage='READING_COMPLETE_NATIVE_RAW_MATRIX',entries=entries,expected_entries=dims[2]))
    assert entries==dims[2] and tracked.bytes==download['bytes'] and tracked.digest.hexdigest()==download['sha256'],'NATIVE_MATRIX_SIGNATURE_OR_LENGTH_CHANGED'
    np.testing.assert_array_equal(totals,knownRNA);np.testing.assert_array_equal(detected,knownGenes)
    mito=counts[:,len(columns):].sum(axis=1)
    np.testing.assert_array_equal(mito,pd.to_numeric(md.total_counts_mt).to_numpy())
    progress(dict(stage='INDEPENDENT_RAW_PREFIX_AND_ALL_UNIT_MARKER_PARITY'))
    prefix=prefix_audit(matrix,genes,physical,counts,totals,detected)
    cached_sum=np.load(native/'sum_raw_counts.npy',mmap_mode='r');cached_positive=np.load(native/'positive_cell_counts.npy',mmap_mode='r')
    for j,col in enumerate(physical):
        raw=np.bincount(unit,weights=counts[:,j],minlength=len(units)).astype(np.int64)
        pos=np.bincount(unit,weights=counts[:,j]>0,minlength=len(units)).astype(np.int64)
        np.testing.assert_array_equal(raw,cached_sum[:,col]);np.testing.assert_array_equal(pos,cached_positive[:,col])
    del cached_sum,cached_positive
    markers=counts[:,:len(columns)];output=[];coverage=[];references=[]
    author_major=md.major_cell_type.to_numpy(str)
    for qc in QC:
        good=native_eligibility(totals,detected,mito) if qc==QC[0] else np.ones(len(unit),bool)
        for parent in PARENTS:
            progress(dict(stage='ALL_CELL_SCALAR_RULE_AND_CASE_GROUP_AUDIT',QC=qc,parent=parent))
            flags=evidence.vector_evidence(markers,good,parent)
            np.testing.assert_array_equal(flags,evidence.scalar_evidence(markers,good,parent))
            np.testing.assert_array_equal(flags[:,1:5].sum(axis=1),good.astype(int))
            aggregate=row_counts(flags,materials,len(lib))
            for name,mask in [('R',flags[:,1]),('E',flags[:,2])]:
                classes,nn=np.unique(author_major[mask],return_counts=True)
                references.append(dict(QC=qc,parent=parent,RNA_parent=name,author_major_reference_counts={str(a):int(b) for a,b in zip(classes,nn)}))
            other=flags[:,0]&~flags[:,1:4].any(axis=1)
            x={g:markers[:,j]>0 for j,g in enumerate(evidence.GENES)}
            cycle=x['MKI67']&(x['TOP2A']|x['UBE2C'])
            others=row_counts(np.column_stack([other&x['TNFRSF9'],other&x['LAG3'],other&cycle]),materials,len(lib))
            b=good&~flags[:,0]&x['MS4A1']&x['CD79A'];epi=good&~flags[:,0]&x['EPCAM']
            allQC=np.bincount(materials[good],minlength=len(lib));quality={};parent_n={}
            for name,mask,den in [('R',flags[:,1],aggregate[:,0]),('E',flags[:,2],aggregate[:,0]),('O',other,aggregate[:,0]),('B',b,allQC),('Epi',epi,allQC)]:
                nn,qq=group_quality(mask,materials,totals,detected,den,len(lib));quality[name]=qq;parent_n[name]=nn
            for i in range(len(lib)):
                output.append(dict(case_index=i,QC=qc,parent=parent,metrics=aggregate[i].tolist(),
                    quality={k:v[i] for k,v in quality.items()},parent_counts={k:int(v[i]) for k,v in parent_n.items()},
                    other_T_states=others[i].tolist(),clinical=[float(v) if np.isfinite(v) else None for v in cr[COMMON].iloc[i]],
                    response=float(cr.response.iloc[i]) if np.isfinite(cr.response.iloc[i]) else None,
                    clinical_eligible=bool(eligible[i]),cancer_type=clinical.cancer_type.iloc[i],drug=clinical['anti-PD1_therapy'].iloc[i]))
            for cancer in ['ALL','LUAD','LUSC']:
                keep=eligible if cancer=='ALL' else eligible&clinical.cancer_type.eq(cancer).to_numpy()
                both=keep&(aggregate[:,1]>=20)&(aggregate[:,2]>=20)
                row=dict(QC=qc,parent=parent,cancer_type=cancer,clinical_metadata_units=int(keep.sum()),
                    R_atleast20=int((keep&(aggregate[:,1]>=20)).sum()),E_atleast20=int((keep&(aggregate[:,2]>=20)).sum()),
                    both_parents_atleast20=int(both.sum()),both_and_otherT_atleast20=int((both&(parent_n['O']>=20)).sum()),
                    MPR_pCR_among_both=int((both&cr.response.eq(1).to_numpy()).sum()),nonMPR_among_both=int((both&cr.response.eq(0).to_numpy()).sum()))
                row['states']={}
                for name in evidence.METRICS[5:11]:
                    j=evidence.METRICS.index(name);den=aggregate[:,1 if name.startswith('R_') else 2]
                    fractions=aggregate[both,j]/den[both]
                    row['states'][name]=dict(source_units=len(fractions),zero_state_units=int((aggregate[both,j]==0).sum()),
                        distinct_proportion_values=len(np.unique(fractions)),variable=bool(len(np.unique(fractions))>1))
                coverage.append(row)
    peak=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*1024
    assert peak<=6*1024**3,'NATIVE_SOURCE_RSS_BUDGET_EXCEEDED'
    proof=dict(source='GSE243013',native_gene_axis_SHA256=sha(gene_path),native_matrix_SHA256=tracked.digest.hexdigest(),
        native_matrix_compressed_bytes=tracked.bytes,native_nonzero_coordinates=entries,all_native_RNA_detected_and_mito_exact=True,
        all_native_unit_marker_sum_and_detection_checks=int(len(units)*len(physical)),prefix=prefix,
        all_cell_scalar_rules_and_grouped_count_quality_means_verified=True,
        all243_original_case_keys_match_clinical=True,clinical_table_SHA256=sha(clinical_path),
        peak_RSS_bytes=peak,no_source_writes_or_new_download=True,no_formal_identity_or_old_models_changed=True,
        original_patient_aliases_not_reviewed=True,clinical_effects_not_computed=True,parent_author_major_reference=references)
    return dict(rows=output,coverage=coverage,proof=proof)

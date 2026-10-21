"""Resolve author-subtype gene marginals without reconstructing coexpression."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(16*1024**2),b''):h.update(b)
    return h.hexdigest()


def summarize_units(units,positive,columns,case_names,definitions):
    ncase=len(case_names);case=pd.Index(case_names).get_indexer(units.sampleID)
    assert (case>=0).all() and not units.duplicated(['sampleID','sub_cell_type']).any()
    output={};totals={}
    for name in definitions['Treg']+definitions['CD4T_nonTreg']:
        use=units.sub_cell_type.eq(name).to_numpy();ii=case[use]
        nc=units.n_cells.to_numpy(np.int64)[use]
        count=np.bincount(ii,weights=nc,minlength=ncase).astype(np.int64)
        q=np.zeros((ncase,2));markers=[]
        for j,key in enumerate(['mean_logRNA','mean_loggenes']):
            sums=np.bincount(ii,weights=nc*units[key].to_numpy()[use],minlength=ncase)
            q[:,j]=sums/np.maximum(count,1)
        for column in columns:
            v=np.asarray(positive[:,column],np.int64)[use]
            assert (v>=0).all() and (v<=nc).all()
            total=np.bincount(ii,weights=v,minlength=ncase).astype(np.int64)
            oracle=pd.Series(v,index=ii).groupby(level=0).sum().reindex(range(ncase),fill_value=0).to_numpy(np.int64)
            np.testing.assert_array_equal(total,oracle);markers.append(total)
        markers=np.asarray(markers).T
        output[name]=dict(count=count,means=q,markers=markers)
        totals[name]=dict(cells=int(count.sum()),marker_positive_totals=markers.sum(axis=0).tolist())
    return output,totals


def augment(payload,root,plan):
    root=Path(root);native=root/'NSCLC_full_native_gene_units_20261016_v1'
    assert sha(native/'private_terminal.json')==plan['base_measurement']['native_terminal_SHA256']
    t=json.loads((native/'private_terminal.json').read_text())
    reg=native/'private_all_native_units_and_quality.tsv';pos=native/'positive_cell_counts.npy'
    for path in [reg,pos]:assert sha(path)==t['outputs'][path.name]['sha256']
    units=pd.read_csv(reg,sep='\t',keep_default_na=False).sort_values('unit_index')
    assert units.unit_index.tolist()==list(range(len(units)))
    gp=root/'GSE243013_genes.csv.gz';assert sha(gp)==plan['native_gene_axis_SHA256']
    genes=pd.read_csv(gp,keep_default_na=False).geneSymbol.str.strip().str.upper()
    assert genes.is_unique;columns=pd.Index(genes).get_indexer(plan['markers']);assert (columns>=0).all()
    positive=np.load(pos,mmap_mode='r',allow_pickle=False);assert positive.shape==(len(units),len(genes))
    names=sorted(set(units.sampleID));assert len(names)==243
    defs=plan['base_measurement']['author_lineage_names']
    totals=payload['proof']['all243_source_QC_parent_cell_totals']
    assert totals['native_QC500_250_MT20']==totals['all_author_retained'], 'QC_SET_IDENTITY_REQUIRED'
    sub,source_totals=summarize_units(units,positive,columns,names,defs)
    for row in payload['rows']:
        i=row['case_index'];row['subtype_measurements']={}
        denominator=sum(row['parent_counts'].values())
        for name,value in sub.items():
            n=int(value['count'][i]);q=value['means'][i]
            row['subtype_measurements'][name]=dict(n=n,markers=value['markers'][i].tolist(),
                quality=[float(q[0]) if n else None,float(q[1]) if n else None,float(np.log1p(n)),n/denominator if denominator else None])
        row['marker_counts']={}
        for parent in ['Treg','CD4T_nonTreg']:
            group=[row['subtype_measurements'][x] for x in defs[parent]]
            assert sum(x['n'] for x in group)==row['parent_counts'][parent]
            row['marker_counts'][parent]=np.sum([x['markers'] for x in group],axis=0).tolist()
            if row['parent_counts'][parent]:
                quality=np.sum([np.asarray(x['quality'][:2])*x['n'] for x in group if x['n']],axis=0)/row['parent_counts'][parent]
                np.testing.assert_allclose(quality,row['quality'][parent][:2],rtol=0,atol=1e-10)
        for j,name in enumerate(defs['Treg']):assert row['subtype_measurements'][name]['n']==row['author_counts'][j]
    coverage=[]
    for name in defs['Treg']:
        rows=[r for r in payload['rows'] if r['QC']=='native_QC500_250_MT20' and r['clinical_eligible']
              and min(r['parent_counts'][p] for p in ['Treg','CD8T','CD4T_nonTreg'])>=20 and r['subtype_measurements'][name]['n']>=20]
        coverage.append(dict(subtype=name,before_drug_qualification=len(rows),MPR_pCR=sum(r['response']==1 for r in rows),nonMPR=sum(r['response']==0 for r in rows),
            all_marker_fractions_variable=[len({r['subtype_measurements'][name]['markers'][j]/r['subtype_measurements'][name]['n'] for r in rows})>1 for j in range(3)]))
    payload['subtype_measurement_proof']=dict(all_original_source_SHA_and_case_parent_keys_checked=True,
        exact_QC_set_identity_requalified=True,all9_author_subtype_counts_sum_to_lineage_counts=True,
        subtype_quality_weighted_means_match_original_percell_parent_quality=True,
        all_marker_totals_independent_groupby_checked=True,single_gene_marginals_not_joint_patterns=True,
        source_subtype_cell_and_marker_totals=source_totals,coverage_before_drug_qualification=coverage,
        marker_names=plan['markers'],original_identity_CM_and_closed_families_unchanged=True)
    return payload

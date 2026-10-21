"""Specific marker capture/quality on fixed author T rows; RAM only.

No patient data, receptor sequences, coordinates, hypotheses or new labels.
read_category is supplied by the existing fixed metadata parser.
"""
from collections import Counter
import hashlib
import io
import json
import resource
import numpy as np


def decode(value):return value.decode() if isinstance(value,bytes) else str(value)


def measure(h,plan):
    obs=h['obs'];key=decode(obs.attrs.get('_index','_index'));n=len(obs[key])
    assert n==plan['expected_observations']
    coarse,cc=read_category(obs['anno'],n);fine,fc=read_category(obs['anno2'],n)
    t_ids=[i for i,v in enumerate(coarse) if v.lower().replace('_','').replace(' ','') in ['tcell','tcells']]
    selected=np.flatnonzero(np.isin(cc,t_ids));assert len(selected)==plan['expected_coarse_T_rows']
    var=h['var'];vk=decode(var.attrs.get('_index','_index'));genes=[decode(v) for v in var[vk][:]]
    assert len(set(genes))==len(genes)
    axis_sha=hashlib.sha256(json.dumps(genes,ensure_ascii=False,separators=(',',':')).encode()).hexdigest()
    assert axis_sha==plan['expected_original_gene_axis_SHA256']
    mapping={g:i for i,g in enumerate(genes)};markers=plan['all_fixed_markers'];columns=[mapping.get(v,-1) for v in markers]
    x=h['X'];assert decode(x.attrs['encoding-type'])=='csr_matrix' and list(x.attrs['shape'])==[n,len(genes)]
    ptr=np.asarray(x['indptr'][:],np.int64);assert ptr.shape==(n+1,) and ptr[0]==0 and (np.diff(ptr)>=0).all()
    assert ptr[-1]==len(x['data'])==len(x['indices'])
    counts=np.zeros((len(selected),len(markers)),np.int64);total=np.zeros(len(selected),np.int64);detected=np.zeros(len(selected),np.int64)
    for i,row in enumerate(selected):
        lo,hi=map(int,ptr[row:row+2]);indices=np.asarray(x['indices'][lo:hi],np.int64);values=np.asarray(x['data'][lo:hi])
        assert np.isfinite(values).all() and (values>0).all() and (values==np.floor(values)).all()
        assert (indices>=0).all() and (indices<len(genes)).all() and (np.diff(indices)>0).all()
        total[i]=values.sum(dtype=np.int64);detected[i]=len(indices)
        positions=np.searchsorted(indices,np.maximum(columns,0))
        for j,(column,position) in enumerate(zip(columns,positions)):
            if column>=0 and position<len(indices) and indices[position]==column:counts[i,j]=int(values[position])
        # Independent dictionary oracle for every selected native row and marker.
        reference=dict(zip(indices.tolist(),values.tolist()))
        assert all(counts[i,j]==reference.get(column,0) for j,column in enumerate(columns))
    fine_codes=fc[selected];groups={name:fine_codes==fine.index(name) if name in fine else np.zeros(len(selected),bool) for name in plan['all_allowed_author_types']}
    groups['OTHER_COARSE_T_UNINTERPRETED']=~np.any(list(groups.values()),axis=0)
    assert np.sum(list(groups.values()),axis=0).tolist()==[1]*len(selected)
    old_checks=0
    for record in plan['original_pilot_panel_capture']:
        label=record['label'];mask=groups[label];panel=plan['original_panels'][record['panel']]
        js=[markers.index(v) for v in panel];a=counts[mask][:,js]
        assert int(mask.sum())==record['n']
        observed=[int(((a>0).sum(1)>=1).sum()),int(((a>0).sum(1)>=2).sum()),int(((a>=2).sum(1)>=2).sum())]
        assert observed==[record[k] for k in ['cells_any_gene','cells_ge2_genes','cells_ge2_genes_each_value_ge2']]
        old_checks+=1
    output=[]
    for name,mask in groups.items():
        nn=int(mask.sum());a=counts[mask];rna=total[mask];ng=detected[mask]
        quality=dict(n=nn,zero_deposited_total=int((rna==0).sum()),
            deposited_total_quartiles=np.quantile(rna,[.25,.5,.75]).tolist() if nn else None,
            detected_gene_quartiles=np.quantile(ng,[.25,.5,.75]).tolist() if nn else None,
            rows_RNA500_genes250=int(((rna>=500)&(ng>=250)).sum()))
        rows=[]
        for j,gene in enumerate(markers):
            positive=int((a[:,j]>0).sum());rawsum=int(a[:,j].sum())
            assert positive==sum(int(row[j]>0) for row in a) and rawsum==sum(int(row[j]) for row in a)
            rows.append(dict(marker=gene,unique_axis_present=columns[j]>=0,positive_rows=positive if columns[j]>=0 else None,
                deposited_value_sum=rawsum if columns[j]>=0 else None,
                fraction_positive=positive/nn if nn and columns[j]>=0 else None,
                mean_row_fraction_of_deposited_total=float(np.mean(a[:,j]/np.maximum(rna,1))) if nn and columns[j]>=0 else None))
        output.append(dict(author_type=name,quality=quality,all_marker_marginals=rows))
    assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<plan['RSS_cap_KiB']
    return dict(status='COMPLETE_PILOT_AUTHOR_T_MARKER_MARGINALS_AND_QUALITY_NO_EFFECT_TESTS',source='GSE328481',scope='Same effect-blind minimum-size original file, not full11',
        author_observations=n,coarse_T_rows=len(selected),all_types_and_markers=output,
        complete_all_selected_row_sparse_dictionary_and_group_sum_oracles=True,
        all_preexisting_pilot_author_type_panel_capture_exactly_replayed=old_checks,
        original_unique_gene_axis_SHA256=axis_sha,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        cell_type_name_not_biological_function_truth=True,no_individual_expression_coordinates_receptor_sequences_or_identifiers_output=True,
        assay_positive_integer_deposited_not_verified_original_UMI=True,not_same_clinical_state_patient_phase_or_original_CM_validation=True)


def synthetic():
    import h5py
    b=io.BytesIO()
    with h5py.File(b,'w') as h:
        obs=h.create_group('obs');obs.create_dataset('_index',shape=(5,),dtype='i8')
        def cat(key,labels,codes):
            g=obs.create_group(key);g.create_dataset('categories',data=np.array(labels,dtype='S'));g.create_dataset('codes',data=np.array(codes,np.int8))
        cat('anno',['Tcell','cancer_cell'],[0,0,0,1,0]);cat('anno2',['CD4_Treg','CD8_CXCL13','private_unknown'],[0,1,2,0,-1])
        genes=['CXCL13','FOXP3','TRDC'];v=h.create_group('var');v.create_dataset('_index',data=np.array(genes,dtype='S'))
        x=h.create_group('X');x.attrs['encoding-type']='csr_matrix';x.attrs['shape']=[5,3]
        x.create_dataset('indptr',data=[0,1,2,2,3,3]);x.create_dataset('indices',data=[1,0,2]);x.create_dataset('data',data=[2,3,7])
        plan=dict(expected_observations=5,expected_coarse_T_rows=4,all_fixed_markers=['TRDC','FOXP3','CXCL13','NOT_MEASURED'],
            all_allowed_author_types=['CD4_Treg','CD8_CXCL13','gdT'],original_pilot_panel_capture=[],original_panels={},RSS_cap_KiB=1024**2,
            expected_original_gene_axis_SHA256=hashlib.sha256(json.dumps(genes,separators=(',',':')).encode()).hexdigest())
        r=measure(h,plan);types={v['author_type']:v for v in r['all_types_and_markers']}
        assert types['CD4_Treg']['all_marker_marginals'][1]['deposited_value_sum']==2
        assert types['CD8_CXCL13']['all_marker_marginals'][2]['positive_rows']==1
        assert types['OTHER_COARSE_T_UNINTERPRETED']['quality']['zero_deposited_total']==2
        assert types['gdT']['quality']['n']==0 and types['gdT']['all_marker_marginals'][-1]['unique_axis_present'] is False
        assert 'private_unknown' not in json.dumps(r)


if __name__=='__main__':
    from qualify_author_spatial_metadata_RAM_v2 import read_category
    synthetic();print('PASS_SPARSE_SELECTED_T_ROWS_REVERSED_MARKERS_MISSING_GENE_AND_NONEMISSION')

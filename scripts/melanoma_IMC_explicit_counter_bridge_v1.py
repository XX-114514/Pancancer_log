"""Normalize only complete declared export-counter universes, never row order."""
import numpy as np
import pandas as pd


def counter(values):
    values=list(values)
    numbers=np.asarray([int(value) for value in values],dtype=np.int64)
    assert all(str(number)==value for number,value in zip(numbers,values))
    assert len(numbers)>0 and len(np.unique(numbers))==len(numbers)
    low,high=int(numbers.min()),int(numbers.max())
    assert (low,high) in [(0,len(numbers)-1),(1,len(numbers))], 'not_a_complete_standard_counter'
    return numbers-low,low


def bind(frame,nodes,aliases):
    assert len(frame)==len(nodes) and set(nodes)=={'','Imageid','type','cell_id'}
    main,main_base=counter(frame['Unnamed: 0'])
    author,author_base=counter(nodes[''])
    assert main_base==1 and author_base==0,'unexpected_source_counter_convention'
    assert all(value=='cell'+index for value,index in zip(nodes['cell_id'],nodes['']))
    registry=set(frame['sample_id'])
    images=[]
    for image in nodes['Imageid']:
        assert image in registry or image in aliases
        images.append(image if image in registry else aliases[image])
    keys=pd.MultiIndex.from_arrays([frame['sample_id'],main])
    normalized=pd.MultiIndex.from_arrays([images,author])
    target=pd.MultiIndex.from_arrays([images,nodes['cell_id']])
    assert keys.is_unique and normalized.is_unique and target.is_unique
    positions=keys.get_indexer(normalized)
    assert (positions>=0).all() and len(np.unique(positions))==len(frame)
    assert np.array_equal(frame['celltype20'].to_numpy()[positions],nodes['type'].to_numpy()),'counter_bridge_type_mismatch'
    return target,positions,dict(main_complete_counter_base=main_base,author_complete_counter_base=author_base,
        all_source_cell_tokens_exact_export_counter=True,all_regions_and_types_exact=True,
        complete_counter_universes_verified=True,no_free_offset_search_or_row_order_join=True)

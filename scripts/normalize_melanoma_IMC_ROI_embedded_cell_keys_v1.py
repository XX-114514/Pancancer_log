"""Reversible explicit-ROI embedding only; never offsets or row-order joins."""
import numpy as np
import pandas as pd


def embedded_variants(values, images, aliases):
    values = list(values)
    images = list(images)
    reverse = {}
    for raw, canonical in aliases.items():reverse.setdefault(canonical, []).append(raw)
    variants = {'literal':values}
    probes = sorted(set(list(range(min(512,len(values))))+list(range(max(0,len(values)-512),len(values)))))
    for side in ['prefix','suffix']:
        for delimiter in ['_', '__', ':', '-', '.']:
            def split(index):
                value, image = values[index],images[index]
                choices = sorted(set([image]+reverse.get(image,[])),key=len,reverse=True)
                matches = []
                for choice in choices:
                    anchor = choice+delimiter if side=='prefix' else delimiter+choice
                    if side=='prefix' and value.startswith(anchor):matches.append((choice,value[len(anchor):]))
                    if side=='suffix' and value.endswith(anchor):matches.append((choice,value[:-len(anchor)]))
                if len(matches)!=1 or not matches[0][1]:return None
                choice,token = matches[0]
                assert value == (choice+delimiter+token if side=='prefix' else token+delimiter+choice)
                return token
            if any(split(index) is None for index in probes):continue
            tokens = [split(index) for index in range(len(values))]
            if any(token is None for token in tokens):continue
            variants[side+'_'+repr(delimiter)] = tokens
    return variants


def bind_embedded_nodes(frame, nodes, aliases):
    assert set(nodes)=={'','Imageid','type','cell_id'} and len(nodes)==len(frame)
    registry = set(frame['sample_id'])
    images = []
    for value in nodes['Imageid']:
        assert value in registry or value in aliases
        images.append(value if value in registry else aliases[value])
    target = pd.MultiIndex.from_arrays([images,nodes['cell_id']])
    assert target.is_unique
    source_variants = embedded_variants(nodes['cell_id'],images,aliases)
    accepted = []
    maps = []
    for field in ['X','Unnamed: 0']:
        main_variants = embedded_variants(frame[field],frame['sample_id'],aliases)
        for source_format,source_values in source_variants.items():
            normalized_target = pd.MultiIndex.from_arrays([images,source_values])
            if not normalized_target.is_unique:continue
            for main_format,main_values in main_variants.items():
                if main_format==source_format=='literal':continue
                keys = pd.MultiIndex.from_arrays([frame['sample_id'],main_values])
                if not keys.is_unique:continue
                locations=keys.get_indexer(normalized_target)
                if (locations<0).any():continue
                if not np.array_equal(frame['celltype20'].to_numpy()[locations],nodes['type'].to_numpy()):continue
                accepted.append(dict(identity_field=field,main_format=main_format,source_format=source_format))
                maps.append(locations)
    assert maps,'no_complete_reversible_ROI_embedded_key_and_label_bijection'
    assert all(np.array_equal(maps[0],item) for item in maps[1:]),'ambiguous_complete_ROI_embedded_cell_identity'
    # Keep ORIGINAL author IDs for graph endpoints; normalization is only the
    # proven bridge to main rows, never a change to author edge identifiers.
    return target,maps[0],accepted

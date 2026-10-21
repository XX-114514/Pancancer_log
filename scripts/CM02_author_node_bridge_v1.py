"""Source-level reference/node bridge; no individual arrays or learned weights."""
from collections import Counter


def bridge(cross,network,fine):
    nodes=set(network['network_labels_to_authority_ID'].values())
    inverse={v:k for k,v in network['network_labels_to_authority_ID'].items()}
    counts=Counter()
    for row in cross:
        assert isinstance(row['matched_cells'],int) and row['matched_cells']>0
        if row['current_major']=='T':counts[row['phase'],row['author_label'],row['current_node']]+=row['matched_cells']
    for row in fine:
        selected=[n for (p,l,n) in counts if p==row['phase'] and l==row['author_fine_label']]
        if row['context']=='CM02_network_member_union_author_reference':selected=[n for n in selected if n in nodes]
        assert sum(counts[row['phase'],row['author_fine_label'],n] for n in selected)==row['unique_reference_cell_observations']
    result=[]
    for phase,label in sorted({(p,l) for p,l,n in counts if p in ['pre','post']}):
        selected=[(n,v) for (p,l,n),v in counts.items() if p==phase and l==label];total=sum(v for n,v in selected);inside=sum(v for n,v in selected if n in nodes)
        result.append(dict(phase=phase,author_fine_label=label,unique_current_T_reference_cells=total,
            inside_CM02_cells=inside,outside_CM02_cells=total-inside,inside_CM02_fraction=inside/total if total else None,
            node_counts=[dict(current_node=n,known_CM02_label=inverse.get(n),in_CM02=n in nodes,source_cells=v) for n,v in sorted(selected)],
            unique_reference_overlap_excluded=True,no_new_state_identity_or_patient_independence_claim=True))
    assert len(result)==26
    return result

#!/usr/bin/env python3
"""Resolve source feature identities without guessing ambiguous gene aliases."""
from collections import Counter, defaultdict
import re

import pandas as pd

ID_PATTERN=re.compile(r'ENSG[0-9]{11}(?:\.[0-9]+)?(?:_PAR_Y)?')


def stable_ids(symbol,row):
    values=[str(row[c]) for c in row.index if c.startswith('source_feature') or c=='gene_id_original']
    if ID_PATTERN.fullmatch(str(symbol)):values.append(str(symbol))
    return {re.sub(r'\.[0-9]+(?=_PAR_Y$|$)','',x) for value in values for x in ID_PATTERN.findall(value)}


def build_contract(source_var):
    records={};name_ids=defaultdict(set);id_names=defaultdict(Counter)
    for source,var in source_var.items():
        records[source]=[]
        for symbol,row in var.iterrows():
            symbol=str(symbol);ids=stable_ids(symbol,row)
            records[source].append((symbol,ids));name_ids[symbol].update(ids)
            for identifier in ids:id_names[identifier][symbol]+=1
    mapping={};labels=defaultdict(Counter);counts=Counter();identities={}
    for source,rows in records.items():
        mapping[source]=[]
        for symbol,ids in rows:
            if len(ids)==1:
                key=next(iter(ids));status='source_stable_ID'
            elif len(ids)>1:
                key='AMBIGUOUS::'+source+'::'+symbol;status='source_collapsed_multiple_IDs'
            elif len(name_ids[symbol])==1:
                key=next(iter(name_ids[symbol]));status='unique_cross_source_symbol_to_ID'
            elif len(name_ids[symbol])>1:
                key='AMBIGUOUS::'+source+'::'+symbol;status='symbol_has_multiple_known_IDs'
            else:
                key='SYMBOL::'+symbol;status='source_symbol_no_stable_ID'
            mapping[source].append(key);labels[key][symbol]+=1;counts[status]+=1
            identities.setdefault(source,[]).append(dict(source_feature=symbol,canonical_feature=key,
                source_stable_IDs='|'.join(sorted(ids)),identity_status=status))
    genes=sorted(labels);data=[]
    for key in genes:
        aliases=labels[key];named=[n for n in aliases if not ID_PATTERN.fullmatch(n)]
        representative=sorted(named,key=lambda n:(-aliases[n],n))[0] if named else sorted(aliases)[0]
        data.append(dict(gene_symbol=representative,source_gene_aliases='|'.join(sorted(aliases)),
                         feature_identity_status='stable_Ensembl_ID' if ID_PATTERN.fullmatch(key) else
                         'ambiguous_source_feature' if key.startswith('AMBIGUOUS::') else 'symbol_without_stable_ID'))
    var=pd.DataFrame(data,index=pd.Index(genes,name='canonical_gene_identity'))
    report=dict(n_canonical_features=len(genes),n_original_name_union=len(set(name_ids)),
        n_known_Ensembl_IDs=len(id_names),n_same_ID_multiple_names=sum(len(v)>1 for v in id_names.values()),
        n_same_name_multiple_IDs=sum(len(v)>1 for v in name_ids.values()),
        source_feature_resolution_counts=dict(counts),
        n_ambiguous_features=sum(g.startswith('AMBIGUOUS::') for g in genes),
        n_symbol_only_features=sum(g.startswith('SYMBOL::') for g in genes),
        representative_symbol_policy='most_sources_non_Ensembl_name_then_lexical_not_claimed_HGNC_current')
    return mapping,var,identities,report

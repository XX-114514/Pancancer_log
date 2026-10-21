"""Independent CSV feature/namespace replay of the entire external inventory.

No expression, barcode, clinical or sequence payload is read. Source records
and feature filenames remain in RAM; only source aggregates are written.
"""
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import csv
import gzip
import hashlib
import io
import json
import os
import re

import qualify_GSE288199_external_state_measurement_v2 as original


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def decode(url):
    raw=original.fetch(url,2*1024**2)
    data=gzip.decompress(raw)
    rows=list(csv.reader(io.StringIO(data.decode()),delimiter='\t'))
    assert rows and all(len(r)==3 for r in rows)
    matches=Counter();identity=defaultdict(list)
    for accession,symbol,kind in rows:
        if kind=='Gene Expression' and symbol in original.GENES:
            matches[symbol]+=1;identity[symbol].append(accession)
    assert all(matches[g]==1 and re.fullmatch(r'ENSG[0-9]+(?:\.[0-9]+)?',identity[g][0]) for g in original.GENES)
    human=[(a,g) for a,g,t in rows if t=='Gene Expression' and re.fullmatch(r'ENSG[0-9]+(?:\.[0-9]+)?',a)]
    other=[(a,g) for a,g,t in rows if t=='Gene Expression' and not re.fullmatch(r'ENSG[0-9]+(?:\.[0-9]+)?',a)]
    return dict(compressed_SHA256=original.digest(raw),uncompressed_SHA256=original.digest(data),
        features=len(rows),feature_type_counts=dict(Counter(r[2] for r in rows)),
        stable_human_gene_features=len(human),other_gene_expression_features=len(other),
        duplicate_stable_human_identifiers=len(human)-len({a for a,g in human}),
        nonhuman_or_unresolved_feature_symbols=sorted(set(g for a,g in other)),
        target_physical_counts={g:matches[g] for g in original.GENES})


def main():
    src=Path(os.environ['PANCANCER_EXTERNAL_MEASUREMENT_OUTPUT'])
    t=json.loads((src/'private_terminal.json').read_text())
    assert t['status']=='COMPLETE_OFFICIAL_METADATA_AND_ALL_FEATURE_AXIS_QUALIFICATION'
    assert t['source_SHA256']==sha(Path(original.__file__))=='058fe1351829752ffc7fb1d508056fe601ba33d43606c23180aa6fdd054a4022'
    for n,r in t['files'].items():assert sha(src/n)==r['SHA256']
    r=json.loads((src/'private_complete_source_measurement_qualification.json').read_text())
    out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent==Path('/dev/shm') and not out.exists();out.mkdir(mode=0o700)
    def save(n,v):(out/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    save('plan_before_independent_feature_scope_replay.json',dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(Path(__file__)),
        producer_terminal_SHA256=sha(src/'private_terminal.json'),
        scope='All declared feature axes; independent CSV field decoding and physical-target counter. Recheck full digest multiset, source title material/mode inventory and GEX sort labels.',
        extra_qualification='Separate stable human Gene Expression from antibody/other/viral or unresolved features before any RNA QC denominator. Axis length does not assign chemistry.',
        shared='Same official network transport and frozen marker list; independent parsing/counting, no independent source institution truth.',
        no_expression_clinical_or_individual_arrays=True))
    raw=original.fetch(original.BASE+'soft/GSE288199_family.soft.gz',1024**2)
    assert original.digest(raw)==original.SOFT_SHA
    text=gzip.decompress(raw).decode()
    title_pairs=re.findall(r'^!Sample_title = (HN[0-9]+) (Pre|Post) (Gene expression data|CITE-Seq expression data|TCR-Seq data)$',text,re.M)
    assert len(title_pairs)==r['GEO_sample_records']
    material=defaultdict(set)
    for case,phase,mode in title_pairs:material[(case,phase.lower())].add(mode)
    assert len(material)==r['processed_material_triplets']
    assert all(v==set(original.MODALITIES) for v in material.values())
    assert dict(Counter(p for c,p in material))==r['material_phase_counts']
    assert len({c for c,p in material})==r['coded_source_cases']
    visits=defaultdict(set)
    for c,p in material:visits[c].add(p)
    assert sum(len(v)==2 for v in visits.values())==r['coded_cases_with_both_pre_post']
    all_sort=Counter();gex_sort=Counter();sort_by_material=defaultdict(set)
    for block in text.split('^SAMPLE = ')[1:]:
        title=re.search(r'^!Sample_title = (HN[0-9]+) (Pre|Post) (.+)$',block,re.M)
        sort=re.search(r'^!Sample_characteristics_ch1 = cell type: (.+)$',block,re.M)
        assert title and sort
        all_sort[sort.group(1)]+=1
        sort_by_material[(title.group(1),title.group(2).lower())].add(sort.group(1))
        if title.group(3)=='Gene expression data':gex_sort[sort.group(1)]+=1
    assert dict(all_sort)==r['source_trait_counts']['cell type']
    assert sum(len(v)>1 for v in sort_by_material.values())==r['materials_with_literal_sort_label_disagreement']
    urls=[];feature_material=set()
    for value in re.findall(r'^!Series_supplementary_file = (.+)$',text,re.M):
        match=re.search(r'/GSE288199_(HN[0-9]+)_(pre|post)_T_features.tsv.gz$',value)
        if match:
            key=match.groups();assert key not in feature_material;feature_material.add(key)
            urls.append(original.BASE+'suppl/'+value.rsplit('/',1)[-1])
    assert feature_material==set(material) and len(urls)==60
    rows=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        for k,row in enumerate(pool.map(decode,sorted(urls)),1):
            rows.append(row);save('private_progress.json',dict(stage='independent_CSV_feature_replay',completed=k,total=60))
    signature=original.digest(json.dumps(sorted((v['compressed_SHA256'],v['uncompressed_SHA256'],v['features']) for v in rows),separators=(',',':')).encode())
    assert signature==r['source_feature_manifest_SHA256']
    assert dict(Counter(str(v['features']) for v in rows))==r['feature_count_distribution']
    assert all(v['target_physical_counts']=={g:1 for g in original.GENES} for v in rows)
    grouped=defaultdict(list)
    for row in rows:grouped[row['uncompressed_SHA256']].append(row)
    schemas=[]
    for key,values in sorted(grouped.items()):
        first=values[0]
        assert all(v['feature_type_counts']==first['feature_type_counts'] and v['stable_human_gene_features']==first['stable_human_gene_features'] for v in values)
        schemas.append(dict(uncompressed_SHA256=key,source_material_axes=len(values),**{k:first[k] for k in ['features','feature_type_counts','stable_human_gene_features','other_gene_expression_features','duplicate_stable_human_identifiers','nonhuman_or_unresolved_feature_symbols']}))
    save('private_complete_independent_feature_scope_and_namespace.json',dict(
        dataset='GSE288199',all_source_material_modality_counts_exact=True,
        all60_compressed_and_uncompressed_feature_digest_multiset_replayed=True,
        all23_physical_targets_unique_stable_human_genes=True,
        GEX_literal_sort_label_counts=dict(gex_sort),
        feature_schema_groups=schemas,feature_types_not_inferred_from_axis_length=True,
        no_cell_expression_capture_or_clinical_keys_verified=True))
    save('private_terminal.json',dict(status='PASS_COMPLETE_INDEPENDENT_MATERIAL_FEATURE_AXIS_AND_NAMESPACE_REPLAY',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(Path(__file__)),
        producer_terminal_SHA256=sha(src/'private_terminal.json'),
        all60_axis_payloads_and23_physical_gene_targets_replayed=True,
        stable_human_vs_other_feature_types_qualified=True,
        no_expression_barcode_clinical_or_VDJ_payload_read=True,
        no_individual_records_arrays_or_original_models_written=True,
        files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status='PASS_COMPLETE_INDEPENDENT_MATERIAL_FEATURE_AXIS_AND_NAMESPACE_REPLAY',terminal_SHA256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()

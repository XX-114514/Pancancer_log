"""Official metadata and feature-axis qualification v2; no expression/VDJ reads.

GSM and coded case/material records remain in RAM. Only source aggregates,
public study facts, and technical digests are written to a new RAM directory.
"""
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
import gzip
import hashlib
import io
import json
import os
import re
import tarfile
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

import observed_T_marker_evidence_v1 as evidence

BASE='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288199/'
SOFT_SHA='628263e1713f7f3b0d7486885416d6c12c74e7bf985d22df43322a4c5a8651e9'
MINIML_SHA='3725b1968efe18c0f98ad42bfab1c71c776469911816c412b50ff5ac9f91328a'
N={'g':'http://www.ncbi.nlm.nih.gov/geo/info/MINiML'}
MODALITIES={'Gene expression data':'GEX','CITE-Seq expression data':'ADT','TCR-Seq data':'VDJ'}
GENES=evidence.GENES+['CXCL13']


def digest(b):return hashlib.sha256(b).hexdigest()


def fetch(url,limit):
    assert url.startswith(BASE)
    opener=urllib.request.build_opener(urllib.request.ProxyHandler({}))
    with opener.open(url,timeout=30) as r:
        assert r.status==200 and urllib.parse.urlparse(r.url).hostname=='ftp.ncbi.nlm.nih.gov'
        b=r.read(limit+1);assert len(b)<=limit
        return b


def soft_rows(raw):
    rows={};current=None;series=defaultdict(list)
    for line in gzip.decompress(raw).decode().splitlines():
        if line.startswith('^SAMPLE = '):
            accession=line.split(' = ',1)[1];assert accession not in rows
            current=defaultdict(list);rows[accession]=current
        elif ' = ' in line:
            k,v=line.split(' = ',1)
            if current is not None and k.startswith('!Sample_'):current[k.removeprefix('!Sample_')].append(v)
            if k.startswith('!Series_'):series[k.removeprefix('!Series_')].append(v)
    output={}
    for accession,row in rows.items():
        characteristics={}
        for field in row['characteristics_ch1']:
            k,sep,v=field.partition(': ');assert sep and k not in characteristics
            characteristics[k]=v
        assert len(row['title'])==1
        output[accession]=dict(title=row['title'][0],traits=characteristics)
    return output,series


def miniml_rows(raw):
    with tarfile.open(fileobj=io.BytesIO(raw),mode='r:gz') as archive:
        files=[x for x in archive.getmembers() if x.isfile() and x.name.endswith('.xml')]
        assert len(files)==1 and files[0].size<=32*1024**2
        root=ET.fromstring(archive.extractfile(files[0]).read())
    output={}
    for sample in root.findall('g:Sample',N):
        accession=sample.find('g:Accession',N).text
        assert accession not in output
        traits={}
        for field in sample.findall('g:Channel/g:Characteristics',N):
            k=field.attrib['tag'];assert k not in traits;traits[k]=(field.text or '').strip()  # XML indentation only; preserve internal whitespace.
        output[accession]=dict(title=sample.find('g:Title',N).text,traits=traits)
    return output


class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):
        if tag=='a':self.links.extend(v for k,v in attrs if k=='href')


def feature_axis(url):
    raw=fetch(url,2*1024**2);data=gzip.decompress(raw)
    assert len(data)<=16*1024**2
    rows=[line.split('\t') for line in data.decode().splitlines()]
    assert rows and all(len(row)==3 for row in rows)
    genes=[row for row in rows if row[2]=='Gene Expression']
    lookup=defaultdict(list)
    for i,(accession,symbol,kind) in enumerate(rows):
        if kind=='Gene Expression':lookup[symbol].append((i,accession))
    hits={g:len(lookup[g]) for g in GENES}
    stable=all(re.fullmatch(r'ENSG[0-9]+(?:\.[0-9]+)?',accession) for g in GENES for _,accession in lookup[g])
    return dict(compressed_SHA256=digest(raw),uncompressed_SHA256=digest(data),bytes=len(raw),
                features=len(rows),gene_expression_features=len(genes),
                feature_type_counts=dict(Counter(row[2] for row in rows)),
                required_gene_physical_matches=hits,required_matches_stable_human_Ensembl=stable)


def main():
    out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent==Path('/dev/shm') and not out.exists();out.mkdir(mode=0o700)
    def save(name,value):(out/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    save('plan_before_external_measurement_qualification.json',dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=digest(Path(__file__).read_bytes()),
        dataset='GSE288199',scope='Entire official SOFT/MINiML sample inventory, complete processed triplet listing, all declared feature axes.',
        question='Can the published independent HNSCC study measure the frozen marker patterns? Separate technical material identity from verified patients/clinical/spatial acceptance.',
        genes=GENES,source='Explicit direct official HTTPS only; bounded bytes and30second timeout. No fallback mirrors or authenticated sources.',
        original_metadata_SHA256=dict(SOFT=SOFT_SHA,MINiML=MINIML_SHA),
        identity='Exact coded-case+phase title plus modality; all3 modalities must align. Coded identifiers are not cross-study reviewed patients.',
        features='Every declared feature file; gene expression type, exact symbol, unique physical column and stable Ensembl. No expression, barcode or VDJ payload read.',
        clinical='No clinical records obtained; source traits and explicit field absence only. No arm/outcome inferred from case numbers or published plots.',
        no_individual_records_arrays_or_file_names_saved=True,no_original_models_or_authorities_modified=True))
    soft=fetch(BASE+'soft/GSE288199_family.soft.gz',1024**2)
    mini=fetch(BASE+'miniml/GSE288199_family.xml.tgz',1024**2)
    assert digest(soft)==SOFT_SHA and digest(mini)==MINIML_SHA
    first,series=soft_rows(soft);second=miniml_rows(mini)
    assert first==second and len(first)==180
    materials=defaultdict(dict);trait_counts=defaultdict(Counter)
    for row in first.values():
        match=re.fullmatch(r'(HN[0-9]+) (Pre|Post) (Gene expression data|CITE-Seq expression data|TCR-Seq data)',row['title'])
        assert match is not None
        case,phase,description=match.groups();modality=MODALITIES[description]
        assert row['traits']['library type']==modality
        key=(case,phase.lower());assert modality not in materials[key]
        materials[key][modality]=row
        for k,v in row['traits'].items():trait_counts[k][v]+=1
    assert len(materials)==60 and all(set(v)=={'GEX','ADT','VDJ'} for v in materials.values())
    literal_sort_disagreements=sum(len({v['traits']['cell type'] for v in values.values()})>1 for values in materials.values())
    raw=fetch(BASE+'suppl/',1024**2);parser=Links();parser.feed(raw.decode())
    directory={urllib.parse.urlparse(urllib.parse.urljoin(BASE+'suppl/',link)).path.rsplit('/',1)[-1] for link in parser.links}
    declared={urllib.parse.urlparse(v).path.rsplit('/',1)[-1] for v in series['supplementary_file']}
    triplets=defaultdict(set);features=[]
    for name in sorted(declared):
        match=re.fullmatch(r'GSE288199_(HN[0-9]+)_(pre|post)_T_(barcodes.tsv.gz|features.tsv.gz|matrix.mtx.gz)',name)
        if match:
            case,phase,kind=match.groups();triplets[(case,phase)].add(kind)
            assert name in directory
            if kind=='features.tsv.gz':features.append(BASE+'suppl/'+name)
        else:assert name=='GSE288199_feature_reference.csv.gz' and name in directory
    assert set(triplets)==set(materials)
    assert all(v=={'barcodes.tsv.gz','features.tsv.gz','matrix.mtx.gz'} for v in triplets.values())
    summary=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        for k,row in enumerate(pool.map(feature_axis,features),1):
            summary.append(row);save('private_progress.json',dict(stage='official_feature_axes',completed=k,total=len(features)))
    assert len(summary)==60
    phase=Counter(p for _,p in materials);case_visits=defaultdict(set)
    for case,visit in materials:case_visits[case].add(visit)
    field_names=sorted(trait_counts)
    report=dict(dataset='GSE288199',GEO_sample_records=len(first),processed_material_triplets=len(materials),
        coded_source_cases=len(case_visits),material_phase_counts=dict(phase),
        coded_cases_with_both_pre_post=sum(v=={'pre','post'} for v in case_visits.values()),
        modality_record_counts={m:sum(m in v for v in materials.values()) for m in ['GEX','ADT','VDJ']},
        materials_with_literal_sort_label_disagreement=literal_sort_disagreements,
        source_trait_field_names=field_names,source_trait_counts={k:dict(v) for k,v in trait_counts.items()},
        explicit_treatment_response_pathology_fields_present=any(any(t in k.lower() for t in ['treat','response','regimen','patholog','hpv','p16']) for k in field_names),
        SOFT_MINiML_all_sample_titles_and_traits_exact=True,all_processed_triplets_agree_with_title_materials=True,
        feature_axes_checked=len(summary),unique_feature_payloads=len({r['compressed_SHA256'] for r in summary}),
        compressed_feature_bytes=sum(r['bytes'] for r in summary),feature_count_distribution=dict(Counter(r['features'] for r in summary)),
        all_required_gene_axes_unique_and_stable=all(all(n==1 for n in r['required_gene_physical_matches'].values()) and r['required_matches_stable_human_Ensembl'] for r in summary),
        required_gene_axis_support={g:dict(Counter(str(r['required_gene_physical_matches'][g]) for r in summary)) for g in GENES},
        source_metadata_SHA256=dict(SOFT=digest(soft),MINiML=digest(mini),directory=digest(raw)),
        source_feature_manifest_SHA256=digest(json.dumps(sorted((r['compressed_SHA256'],r['uncompressed_SHA256'],r['features']) for r in summary),separators=(',',':')).encode()),
        source_material_inventory_SHA256=digest(json.dumps(sorted(materials),separators=(',',':')).encode()),
        coded_cases_not_verified_independent_patients=True,expression_capture_not_measured=True,
        clinical_patient_mapping_not_accepted=True,spatial_file_and_joint_state_measurement_not_accepted=True)
    save('private_complete_source_measurement_qualification.json',report)
    save('private_terminal.json',dict(status='COMPLETE_OFFICIAL_METADATA_AND_ALL_FEATURE_AXIS_QUALIFICATION',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=digest(Path(__file__).read_bytes()),
        all_SOFT_MINiML_exact_keys_traits_and_processed_triplets_verified=True,all60_feature_axes_checked=True,
        expression_barcode_VDJ_payloads_not_read=True,individual_source_records_and_learned_arrays_not_saved=True,
        no_clinical_or_spatial_acceptance_claim=True,
        files={p.name:dict(SHA256=digest(p.read_bytes()),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status='COMPLETE_OFFICIAL_METADATA_AND_ALL_FEATURE_AXIS_QUALIFICATION',terminal_SHA256=digest((out/'private_terminal.json').read_bytes()))))


if __name__=='__main__':main()

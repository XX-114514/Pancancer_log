"""Independent bounded geometry QC for the complete counter-based cell bridge."""
from pathlib import Path
import csv
import hashlib
import importlib.util
import json
import os
import re
import zipfile
import numpy as np


def load(path,expected,name):
    assert hashlib.sha256(path.read_bytes()).hexdigest()==expected
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def analyze(bundle,parameters):
    assert parameters==dict(maximum_edges=1024,maximum_prefix_bytes=1048576,published_native_radius=20.0)
    folder=Path(__file__).parent
    counter=load(folder/'melanoma_IMC_explicit_counter_bridge_v1.py','9eaee9aae48eb448a7164464d8e0f2a63e29afff9872aa706342e8da687d67da','geometry_counter_bridge')
    helper=load(folder/'qualify_melanoma_IMC_native_graph_key_bridge_local_v4.py','427434209f230231011232b718e286eca77d6b779f72cdbd8d1ac6b182ff6eef','geometry_reference')
    target,locations,proof=counter.bind(bundle['frame'],bundle['author_node_table'],bundle['raw_ROI_aliases'])
    transport=load(Path(os.environ['PANCANCER_IMC_BOUNDED_TRANSPORT']),'b099f8ffe429a856a70a5b9006fe871e944dc539943fb28bbdda354dd9e222c8','geometry_transport')
    rf=transport.RangeFile(None,'https://zenodo.org/api/records/6838169/files/Melanoma_code_data.zip/content',668227028)
    with zipfile.ZipFile(rf) as archive:
        assert hashlib.sha256(json.dumps([(e.filename,e.file_size,e.compress_size,e.CRC) for e in archive.infolist()],separators=(',',':')).encode()).hexdigest()=='26339da50b39ffdf7beb641a67f50523c867f2d3a86b24f42175f3a16dfe8992'
        code=archive.read('Rscript/Spatial_analysis.Rmd')
        assert hashlib.sha256(code).hexdigest()=='d8faddad268f2900a83f8cfa3f8352dd697128d408a528181fb997c5773dba56'
        line=next(v for v in code.decode('utf-8-sig').splitlines() if re.match(r'\s*n_dat\s*<-\s*read.csv',v))
        name=re.search(r'''read.csv\(\s*(["'])([^"'\r\n]+)\1''',line).group(2)
        while name.startswith('../') or name.startswith('./'):name=name.split('/',1)[1]
        rows=[];bytes_read=0
        with archive.open(name) as stream:
            header=stream.readline(16385);assert len(header)<=16384 and header.endswith(b'\n')
            assert len(next(csv.reader([header.decode('utf-8-sig')])))==4
            for _ in range(parameters['maximum_edges']):
                raw=stream.readline(16385)
                if not raw:break
                bytes_read+=len(raw)
                assert len(raw)<=16384 and raw.endswith(b'\n') and bytes_read<=parameters['maximum_prefix_bytes']
                rows.append(next(csv.reader([raw.decode('utf-8-sig')])))
    rf.close()
    geometry=helper.edge_prefix_summary(rows,target,locations,bundle['native_coordinates'],bundle['raw_ROI_aliases'],set(bundle['frame']['sample_id']))
    images=helper.canonical_images([row[1] for row in rows],bundle['raw_ROI_aliases'],set(bundle['frame']['sample_id']))
    nontrivial=all(float(np.ptp(bundle['native_coordinates'][bundle['frame']['sample_id'].to_numpy()==image],axis=0).max())>20 for image in set(images))
    compatible=bool(geometry['native_distance_quantiles'][-1]<=20+1e-6 and geometry['zero_native_distance_edges']<len(rows) and nontrivial)
    return dict(status='COMPLETE_COUNTER_GEOMETRY_PREFIX_QC_NOT_FULL_GEOMETRIC_VALIDATION',
        geometry=geometry,complete_counter_proof=proof,
        native20_compatibility_in_prefix=compatible,ROI_extent_exceeds_test_radius=nontrivial,
        no_micron_calibration_or_all_ROI_geometry_claim=True,
        no_marker_outcome_or_individual_values_written=True,
        producer_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())

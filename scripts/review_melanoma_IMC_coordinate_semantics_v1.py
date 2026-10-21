"""Inspect fixed public author code for coordinate/graph provenance only.

No author code is executed. No cell, coordinate or clinical rows are read.
Quoted literals in returned code are masked, except generic field names.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import re
import zipfile


def sha(value):
    return hashlib.sha256(value).hexdigest()


def load(path, expected, name):
    assert sha(path.read_bytes()) == expected
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def masked(line):
    safe = {'Position', 'cell_id', 'Imageid', 'image_id', 'type', 'sample_id',
            'sample', 'ROI', 'X', 'Y', 'x', 'y', 'from', 'to', 'cell1',
            'cell2', 'dat_cell', 'dat_rel', 'ObjectNumber', 'ImageNumber'}
    def replace(match):
        literal = match.group(0)[1:-1]
        return match.group(0) if literal in safe else '"[masked_literal]"'
    return re.sub(r'''"(?:[^"\\]|\\.)*"|'(?:[^'\\]|\\.)*' ''',
                  replace, line, flags=re.X)


def main():
    source = load(Path(os.environ['PANCANCER_IMC_TRANSPORT_HELPER']),
        '9fb200d98e321584fde8bca6cb92f72151ac33ec00d64e04d93e7e2b33e870e2', 'official_source')
    transport = load(Path(os.environ['PANCANCER_IMC_BOUNDED_TRANSPORT']),
        '709f9b5ced875fea663557bad936f32d004e12d5daf291048b7de476064b0d68', 'bounded_transport')
    record = json.loads(source.fetch(source.BASE+'api/records/6838169', 2*1024**2))
    item = next(x for x in record['files'] if x['key'] == 'Melanoma_code_data.zip')
    assert item['size'] == 668227028
    assert item['checksum'] == 'md5:da74b7228763b76dcd2d9830c1a10cd9'
    rf = transport.RangeFile(source, item['links']['self'], item['size'])
    fixed = {
        'Rscript/functions.R': 'bfb6f99c80947643655aa5d4d8ef310e5f0a1c24aac86326482cb32c123565ee',
        'Rscript/Spatial_analysis.Rmd': 'd8faddad268f2900a83f8cfa3f8352dd697128d408a528181fb997c5773dba56',
        'Rscript/generate_Figs.Rmd': 'b213b123c189457f7870d042094f85415f9dfef9bfc52e296393bfcc9cabf4af'}
    reviews = []
    with zipfile.ZipFile(rf) as archive:
        for member, expected in fixed.items():
            assert archive.getinfo(member).file_size < 512*1024
            blob = archive.read(member)
            assert sha(blob) == expected
            lines = blob.decode('utf-8-sig').splitlines()
            selected = set()
            for index, line in enumerate(lines):
                if re.search(r'Position|[Cc]entroid|Location_Center|generateCellRel\s*(?:<-|=)\s*function|ROI\.meta\s*(?:<-|=)', line):
                    selected.update(range(max(0, index-2), min(len(lines), index+35)))
            reviews.append(dict(member=member, source_SHA256=expected,
                masked_generic_code=[dict(line=i+1, code=masked(lines[i])) for i in sorted(selected)],
                no_selected_data_rows=True, author_code_not_executed=True))
    final = json.loads(source.fetch(source.BASE+'api/records/6838169', 2*1024**2))
    again = next(x for x in final['files'] if x['key'] == item['key'])
    assert all(item.get(k) == again.get(k) for k in ['id', 'size', 'checksum', 'links'])
    rf.close()
    print(json.dumps(dict(status='COMPLETE_FIXED_AUTHOR_COORDINATE_CODE_REVIEW',
        observed_UTC=datetime.now(timezone.utc).isoformat(),
        producer_SHA256=sha(Path(__file__).read_bytes()), reviews=reviews,
        no_cell_clinical_coordinate_rows_or_graph_edges_read=True), indent=2))


if __name__ == '__main__':
    main()

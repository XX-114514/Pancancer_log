"""Inspect accepted author R source as text; never execute it or read rows."""
from pathlib import PurePosixPath, Path
import hashlib
import re


def analyze(bundle, parameters):
    assert parameters == {'maximum_snippets': 75}
    archive=bundle['archive'];accepted=bundle['qualification']['generic_author_code_SHA256']
    rows=[]
    for member in archive.getmembers():
        name=PurePosixPath(member.name).name
        if name not in accepted or not member.isfile():continue
        raw=archive.extractfile(member).read()
        assert hashlib.sha256(raw).hexdigest()==accepted[name]
        source=raw.decode('utf-8-sig')
        snippets=[]
        for number,line in enumerate(source.splitlines(),1):
            line=line.strip()
            if line.startswith('#') or not line:continue
            if not re.search(r'fread|read[._](?:table|csv|delim)|colnames|Ki67|Ki-67|Fox[Pp]3|FOXP3|cell[._]?type|cell[._]?class|Tc[.]|Th[.]|Treg|Clinical|clinical|Location|Centroid|centroid',line):continue
            # Omit source examples that refer to individual specimen tokens.
            if re.search(r'\b\d{1,5}[A-Z]{2,5}\b|\b(?:patient|sample)[_-]\d+\b',line,re.I):continue
            snippets.append(dict(line_number=number,source_expression=line[:300]))
        if snippets:rows.append(dict(generic_source_file=name,code_SHA256=accepted[name],static_expressions=snippets[:parameters['maximum_snippets']]))
    return dict(status='COMPLETE_STATIC_AUTHOR_INPUT_CONTRACT_NOT_TABLE_OR_EFFECT_ACCEPTANCE',
        inspected_source_files=rows,no_code_executed_or_table_rows_read=True,
        producer_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())

"""Study-level GEO fields only; no sample rows, identifiers or clinical bodies.

The incoming HTML remains local RAM. Only fixed study fields can be persisted.
This does not review an individual patient-to-material bridge.
"""
from pathlib import Path
from datetime import datetime, timezone
from html.parser import HTMLParser
from urllib.request import Request, ProxyHandler, build_opener
from urllib.parse import urlparse
import hashlib
import json
import os
import re

FIELDS = ['Title', 'Summary', 'Overall design', 'Citation(s)']


class Text(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.values = []

    def handle_data(self, value):
        self.values.append(value)


def extract(raw):
    html = raw.decode('utf-8')
    lower = html.lower()
    if 'checking your browser' in lower or 'g-recaptcha' in lower:
        raise ValueError('BROWSER_VERIFICATION_NO_BYPASS')
    result = {}
    for field in FIELDS:
        pattern = r'<td\b[^>]*>\s*'+re.escape(field)+r'\s*</td>\s*<td\b[^>]*>(.*?)</td>'
        match = re.search(pattern, html, flags=re.S | re.I)
        if match is None:
            result[field] = None
            continue
        parser = Text()
        parser.feed(match.group(1))
        value = re.sub(r'\s+', ' ', ' '.join(parser.values)).strip()
        assert len(value) <= 12000
        value = re.sub(r'\b(?:GSM|SRR|SRX|SRS|SAMN|SAMEA)\d+\b', '[INDIVIDUAL_ACCESSION_REDACTED]', value)
        result[field] = value
    assert result['Title'] and result['Overall design']
    return result


def synthetic():
    raw = b'<table><tr><td>Title</td><td>STUDY</td></tr><tr><td>Overall design</td><td>Two study arms</td></tr><tr><td>Samples</td><td>DO_NOT_EXPORT_SAMPLE</td></tr></table>'
    q = extract(raw)
    assert q['Title'] == 'STUDY' and 'DO_NOT_EXPORT_SAMPLE' not in json.dumps(q)
    raw = raw.replace(b'Two study arms', b'GSM123456, two study arms')
    assert 'GSM123456' not in json.dumps(extract(raw))


def main():
    task = Path(os.environ['PANCANCER_STUDY_DESIGN_TASK'])
    planfile = task/'private_frozen_official_design_v1.json'
    plan = json.loads(planfile.read_text())
    assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest() == plan['code_SHA256']
    assert plan['fields'] == FIELDS
    synthetic()
    started = datetime.now(timezone.utc).isoformat()
    opener = build_opener(ProxyHandler({'https': os.environ['PANCANCER_METADATA_HTTP_PROXY']}))
    phase = 'OFFICIAL_GEO_STUDY_HTML'
    try:
        url = plan['url']
        p = urlparse(url)
        assert p.scheme == 'https' and p.hostname == 'www.ncbi.nlm.nih.gov' and p.path == '/geo/query/acc.cgi'
        with opener.open(Request(url), timeout=40) as response:
            assert response.status == 200 and response.url == url
            raw = response.read(plan['maximum_HTML_bytes']+1)
        assert len(raw) <= plan['maximum_HTML_bytes']
        phase = 'STUDY_FIELDS_ONLY'
        q = extract(raw)
        digest = hashlib.sha256(raw).hexdigest()
        del raw
        q = dict(status='COMPLETE_OFFICIAL_STUDY_FIELDS_NO_INDIVIDUAL_BRIDGE', source='GSE243013',
            study_fields=q, started_UTC=started, completed_UTC=datetime.now(timezone.utc).isoformat(),
            official_URL=url, source_HTML_SHA256=digest, frozen_contract_SHA256=hashlib.sha256(planfile.read_bytes()).hexdigest(),
            no_sample_rows_sequence_clinical_body_or_learned_weights_parsed=True,
            no_new_expression_or_effect_tests=True, not_individual_patient_independence_acceptance=True)
        outfile = task/'private_official_study_design_v1.json'
        assert not outfile.exists()
        outfile.write_text(json.dumps(q, indent=2)+'\n')
        print(json.dumps(dict(status=q['status'], study_fields=q['study_fields'])), flush=True)
    except Exception as error:
        q = dict(status='OFFICIAL_STUDY_DESIGN_NOT_ACCEPTED', observed_UTC=datetime.now(timezone.utc).isoformat(),
            phase=phase, error_type=type(error).__name__, HTTP_code=getattr(error, 'code', None), no_retry_or_bypass=True)
        (task/'private_official_study_design_failed_v1.json').write_text(json.dumps(q, indent=2)+'\n')
        raise RuntimeError('OFFICIAL_STUDY_DESIGN_NOT_ACCEPTED') from None


if __name__ == '__main__':
    main()

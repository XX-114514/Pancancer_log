"""Read one bounded public supplement for exact clinical-link opportunities only."""
import gzip
import io
import json
import os
import re
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime, timezone
from pathlib import Path

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,12*1024**2)
    original=root/'GSE305078_official_sampling_sidecar_20261016_v1/private_original_GEO_family.soft.gz'
    assert sha(original)=='e24c261f5c44157f376c7c1259b779bd77dfe5ef3bbe7eb78dff2ee3c3182d49'
    official=gzip.decompress(original.read_bytes()).decode('utf-8')
    titles=re.findall(r'^!Sample_title = (.+)$',official,re.M)
    assert len(titles)==len(set(titles))==32
    url='https://pmc.ncbi.nlm.nih.gov/articles/instance/13355034/bin/mmc1.docx'
    output=root/'GSE305078_public_clinical_link_qualification_20261018_v1'
    assert not output.exists();output.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        public_source_article='https://pmc.ncbi.nlm.nih.gov/articles/PMC13355034/',supplement_url=url,
        source_role='2026 same-GSE followup manuscript; not an independent cohort or original CM validation.',
        download='Only one published3.3MB DOCX supplement, maximum8MiB. No expression/sequence/PDF images or controlled files; no access-control bypass or automatic retry.',
        qualification='ZIP CRC and bounded XML, readable paragraph/table structure and exact official sample-title tokens. No clinical association, label adoption, model/gene selection or modification of existing authority.',
        prior_sampling='GEO supports32 pretreatment materials. Followup manuscript describes26 early cases with chemotherapy-only and chemo-immunotherapy; current28 heldout materials are not automatically that26 or a pure ICI set.',
        privacy='All document bytes, text, individual tokens/rows and any clinical fields remain server-private. Emit only bounded technical status.')
    (output/'plan_before_public_supplement_get.json').write_text(json.dumps(plan,indent=2)+'\n')
    status='PUBLIC_SUPPLEMENT_UNAVAILABLE_NO_LABEL_ADOPTED';details={}
    try:
        request=urllib.request.Request(url,headers={'User-Agent':'PanCancer-public-source-qualification/1.0'})
        with urllib.request.urlopen(request,timeout=45) as response:
            assert response.status==200
            payload=response.read(8*1024**2+1)
            details.update(HTTP_status=response.status,content_type=response.headers.get('Content-Type'),downloaded_bytes=len(payload))
        assert len(payload)<=8*1024**2,'SUPPLEMENT_EXCEEDS_FROZEN_LIMIT'
        path=output/'private_original_public_supplement.bin';path.write_bytes(payload)
        assert payload[:4]==b'PK\x03\x04','NOT_A_DOCX_ZIP_NO_ADOPTION'
        with zipfile.ZipFile(io.BytesIO(payload)) as archive:
            assert sum(m.file_size for m in archive.infolist())<=64*1024**2
            assert archive.testzip() is None
            xml=archive.read('word/document.xml')
            assert len(xml)<=8*1024**2
            tree=ET.fromstring(xml)
            ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
            paragraphs=[''.join(n.itertext()) for n in tree.findall('.//w:p',ns)]
            # XML contains formatting nodes; actual text only, no linked content.
            paragraphs=[''.join(n.text or '' for n in p.findall('.//w:t',ns)) for p in tree.findall('.//w:p',ns)]
            text='\n'.join(paragraphs)
            tables=[]
            for table in tree.findall('.//w:tbl',ns):
                rows=[]
                for row in table.findall('./w:tr',ns):
                    rows.append([''.join(n.text or '' for n in cell.findall('.//w:t',ns)) for cell in row.findall('./w:tc',ns)])
                tables.append(rows)
            matches={s:len(re.findall(r'(?<![A-Za-z0-9])'+re.escape(s)+r'(?![A-Za-z0-9])',text)) for s in titles}
            media=sum(m.filename.startswith('word/media/') and not m.is_dir() for m in archive.infolist())
        private=dict(paragraphs=paragraphs,tables=tables,exact_official_sample_title_occurrences=matches,
            interpretation='No automatic patient/treatment/response adoption. Embedded figure contents not machine-read by this qualification; no absence claim for images.')
        (output/'private_readable_supplement_structure_and_token_check.json').write_text(json.dumps(private,indent=2)+'\n')
        details.update(zip_CRC_pass=True,readable_paragraphs=len(paragraphs),tables=len(tables),embedded_media_files=media,
            official_titles_matching_readable_text=sum(v>0 for v in matches.values()),
            exact_patient_treatment_outcome_join_established=False,
            embedded_figure_or_image_content_reviewed=False)
        status='READABLE_PUBLIC_SUPPLEMENT_QUALIFIED_EXACT_CLINICAL_LINK_STILL_UNVERIFIED'
    except Exception as error:
        details.update(error_type=type(error).__name__,error=str(error))
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(output.iterdir()) if p.is_file()}
    (output/'private_terminal.json').write_text(json.dumps(dict(status=status,completed_UTC=datetime.now(timezone.utc).isoformat(),
        source_sha256=sha(__file__),details=details,files=files,patient_identity_or_clinical_labels_adopted=False,
        original_MCP_or_CM_or_clinical_contract_modified=False),indent=2)+'\n')
    print(status)


if __name__=='__main__':main()

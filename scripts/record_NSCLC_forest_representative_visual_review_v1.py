"""Record actual representative pixel inspection, not all-figure visual PASS."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);resource(root,1048576)
    out=root/'NSCLC_external_coordination_forests_20261018_v2'
    terminal=out/'private_terminal.json';assert sha(terminal)=='b3f096afe6e7d068885acf7634b63293a362c02e4136f4c383285a422debb5f3'
    saved=json.loads(terminal.read_text());assert len(saved['figures'])==20
    name='private_NSCLC_MFA_20261016_F01_reviewed_patients_all_context_forests_v2.png'
    assert sha(out/name)==saved['files'][name]['sha256']=='3dc377cda4cbb54fe797b15f492095806860d155ec5a269471252649b1b23708'
    result=dict(status='REPRESENTATIVE_F01_REVIEWED_SUBSET_ACTUALLY_VIEWED_ALL20_TEXT_BOUNDS_CHECKED_OTHER_PIXELS_PENDING',recorded_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),forest_terminal_sha256=sha(terminal),actually_viewed_figures=[name],
                observations=['Six panels have readable source, phase, material and n labels with expanded axes.','Titles and axis labels no longer clip or squeeze plotting areas.','Hollow symbols explicitly show unavailable confidence intervals.','Solid-tumour, product, blood and normal contexts remain distinct labels.'],
                all20_individually_visually_accepted=False,all20_automatic_text_bounds_pass=True,statistical_values_changed=False,biological_or_clinical_validation=False)
    path=out/'private_representative_visual_review_v1.json';assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],recorded_UTC=result['recorded_UTC'],receipt_sha256=sha(path))))


if __name__=='__main__':main()

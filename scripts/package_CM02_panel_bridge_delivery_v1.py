"""Private aggregate report/figure bundle; rows and model coefficients excluded."""
import argparse
import hashlib
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from prepare_CM02_panel_bridge_reference_scope_v1 import budget

def main(a):
    reserve=budget(a.root,64*1024**2)
    reports=list(a.root.glob('*.html'));assert len(reports)==1
    files=reports+[a.root/'private_final_visual_review_v1.json',a.root/'private_final_handoff_v1.json']
    for folder in ['validation_figures_v3','conditional_cytotoxic_figures_v2','confusion_depth_diagnostics_v2','model_independent_audit_v2','cytotoxic_description_independent_audit_v2','GeoMx_server_independent_audit_v3']:
        assert (a.root/folder).is_dir()
        files.extend(p for p in (a.root/folder).rglob('*') if p.is_file())
    files.extend(a.root/'model_validation_v1'/name for name in ['private_all_model_metrics.tsv','private_all_class_precision_recall.tsv','private_all_frozen_engineering_gates.tsv','private_missing_spatial_marker_recovery_diagnostics.tsv','private_spatial_admission_decision.json'])
    files.append(a.root/'cytotoxic_author_T_description_v1/private_all4_author_T_conditional_description_coverage.tsv')
    out=a.root/'delivery_v1';assert not out.exists();out.mkdir(mode=0o700)
    package=out/'CM02_panel_bridge_audit_and_figures_20261020.zip'
    with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for path in files:
            assert path.is_file() and not path.is_symlink();z.write(path,path.relative_to(a.root))
    with zipfile.ZipFile(package) as z:assert z.testzip() is None
    receipt=dict(status='PRIVATE_REPORT_FIGURES_AND_AGGREGATE_BUNDLE_CREATED_AND_ZIP_READBACK_ACCEPTED',finished_utc=datetime.now(timezone.utc).isoformat(),
                 bytes=package.stat().st_size,sha256=hashlib.sha256(package.read_bytes()).hexdigest(),members=len(files),
                 no_cell_patient_rows_coefficients_or_raw_matrices_in_bundle=True,
                 remaining_project_quota_bytes=reserve['group_remaining_bytes'],reserve_bytes=reserve['reserve_bytes'],
                 producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'delivery_technical_public.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args())

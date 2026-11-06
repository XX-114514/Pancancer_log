"""Run a registered scientific stage against an explicit private project root."""
import argparse,runpy
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--project-root',type=Path,required=True)
p.add_argument('--stage',choices=['export_B_profile.py','profile_B_QC_pairs.py','apply_B_QC_and_export_refit_v2.py','audit_B_clean_inputs_v2.py','aggregate_B_markers.py','audit_B_postfit.py','supplement_B_QC_metrics.py','audit_B_additional_markers.py','write_B_working_sidecars.py','audit_B_working_sidecars.py','consolidate_current_annotation_v3.py','audit_consolidated_annotation_v3.py','audit_consolidated_annotation_v4.py'],required=True)
a=p.parse_args()
runpy.run_path(str(Path(__file__).with_name(a.stage)),init_globals={'PROJECT_ROOT':a.project_root.resolve()})

"""Execute one registered stage against a user-supplied private project root."""
import argparse,runpy
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,required=True);p.add_argument('--stage',choices=['profile_full_axis.py','audit_full_axis_profile.py','propose_marker_names.py','refine_components_v2.py','audit_refined_components_v2.py','propose_marker_names_v2.py','write_marker_freeze_v2.py','write_marker_freeze_v3.py','audit_marker_freeze_v2.py','aggregate_frozen_markers.py','audit_frozen_naming_evidence.py','summarize_B_coarse_repair.py'],required=True);a=p.parse_args();runpy.run_path(str(Path(__file__).with_name(a.stage)),init_globals={'PROJECT_ROOT':a.project_root.resolve()})

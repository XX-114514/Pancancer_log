"""Run one registered selective-subtype correction stage with a supplied private root."""
import argparse,runpy
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,required=True);p.add_argument('--stage',choices=['restore_subtypes_selectively_v4.py','audit_selective_restore_v4.py','aggregate_selective_markers_v4.py','audit_selective_markers_v4.py'],required=True);a=p.parse_args();runpy.run_path(str(Path(__file__).with_name(a.stage)),init_globals={'PROJECT_ROOT':a.project_root.resolve()})

"""Run a registered scientific stage against an explicit private project root."""
import argparse,os,runpy
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--project-root',type=Path,required=True)
p.add_argument('--stage',choices=['export_profile_candidates.py','profile_disjoint_QC_pairs.py','apply_QC_and_export_refits.py','audit_clean_inputs.py','aggregate_postfit_markers.py','aggregate_extended_markers.py','audit_postfit_v2.py','write_working_sidecars_v2.py','audit_working_sidecars.py','review_ILC_candidates.py'],required=True)
p.add_argument('--lineage',choices=['Erythroid','Platelet','HSPC','ILC'])
a=p.parse_args()
if a.lineage:os.environ['REVIEW_LINEAGE']=a.lineage
runpy.run_path(str(Path(__file__).with_name(a.stage)),init_globals={'PROJECT_ROOT':a.project_root.resolve()})

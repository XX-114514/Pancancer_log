"""Meaningful gates: publication privacy, count phases, identity release boundary."""
import copy
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from collect_runtime_snapshot import scalar_summary, proof_summary
from validate_runtime_snapshot import validate

def fixture():
    return {"schema":"pancancer_cloud_snapshot_v1","collected_utc":"2026-10-11T00:00:00+00:00",
      "expansion":{"datasets":[{"dataset":"SYNTHETIC","source_cells":100,"download_verified":True,
        "QC":{"status":"qc_shard_written_pending_final_source_acceptance","n_qc_nonduplicate":90}}],
        "source_cells_verified_download_or_inventory":100,"qc_cells_pending_final_source_acceptance":90,
        "final_accepted_new_cells":0,"final_acceptance_status":"no_final_acceptance_manifest_registered"},
      "project_v3":{"identity_release_readiness":{"release_eligible":False,"remaining_GSE299623_cells":10,"approved_remaining_group_decisions":0}}}

class Contracts(unittest.TestCase):
    def test_pending_QC_is_not_final_acceptance(self):
        s=fixture();self.assertEqual(validate(s),[])
        s["expansion"]["final_accepted_new_cells"]=90
        self.assertTrue(any("promoted" in e for e in validate(s)))
    def test_metadata_only_cells_do_not_enter_download_total(self):
        s=fixture();s["expansion"]["datasets"].append({"dataset":"META_ONLY","source_cells":80,"download_verified":False})
        self.assertEqual(validate(s),[])
        s["expansion"]["source_cells_verified_download_or_inventory"]=180
        self.assertTrue(any("Source cell total" in e for e in validate(s)))
    def test_export_whitelists_exclude_patient_and_commands(self):
        data={"status":"PASS","samples":2,"patient_id":"SYNTHETIC_PERSON","sample_reports":[{}],"command":"sensitive runtime arguments","path":"local host path"}
        self.assertEqual(scalar_summary(data),{"status":"PASS","samples":2})
        self.assertEqual(proof_summary(data),{"status":"PASS"})
    def test_snapshot_rejects_row_level_detail(self):
        s=fixture();s["sample_reports"]=[{"patient_uid":"SYNTHETIC_PERSON"}]
        self.assertTrue(any("Disallowed detail" in e for e in validate(s)))
    def test_identity_release_needs_remaining_review(self):
        s=fixture();s["project_v3"]["identity_release_readiness"]["release_eligible"]=True
        self.assertTrue(any("Identity release" in e for e in validate(s)))

if __name__=="__main__":unittest.main()

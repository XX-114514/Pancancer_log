#!/usr/bin/env python3
"""Validate publication boundaries and count semantics of an aggregate snapshot."""
import argparse
import datetime as dt
import json
from pathlib import Path
import re

def validate(s):
    errors=[]
    if s.get("schema")!="pancancer_cloud_snapshot_v1":errors.append("Unsupported snapshot schema")
    try:
        timestamp=dt.datetime.fromisoformat(s["collected_utc"])
        if timestamp.tzinfo is None:errors.append("Snapshot timestamp lacks timezone")
    except (KeyError,ValueError,TypeError):errors.append("Invalid collection timestamp")
    exp=s.get("expansion",{});datasets=exp.get("datasets",[])
    ids=[d.get("dataset") for d in datasets]
    if len(ids)!=len(set(ids)):errors.append("Duplicate dataset registration")
    total=sum(d.get("source_cells",0) for d in datasets if d.get("download_verified"))
    if total!=exp.get("source_cells_verified_download_or_inventory"):errors.append("Source cell total does not match verified source inventory")
    qc_total=sum(d.get("QC",{}).get("n_qc_nonduplicate",0) for d in datasets)
    if qc_total!=exp.get("qc_cells_pending_final_source_acceptance"):errors.append("QC shard total mismatch")
    if exp.get("final_accepted_new_cells",0)>0 and exp.get("final_acceptance_status")=="no_final_acceptance_manifest_registered":
        errors.append("QC/download cells cannot be promoted to final acceptance")
    for d in datasets:
        qc=d.get("QC",{})
        if not 0<=qc.get("n_qc_nonduplicate",0)<=d.get("source_cells",0):errors.append("Invalid QC count: "+str(d.get("dataset")))
        if qc and qc.get("status")!="qc_shard_written_pending_final_source_acceptance":errors.append("Unrecognized QC proof requires explicit schema migration")
    state=s.get("project_v3",{}).get("identity_release_readiness",{})
    if state.get("release_eligible") and state.get("remaining_GSE299623_cells",0)>0 and state.get("approved_remaining_group_decisions",0)==0:
        errors.append("Identity release contradicts unreviewed remaining groups")
    serialized=json.dumps(s,ensure_ascii=False)
    if re.search(r"/(?:home|Users|data\d*)/",serialized):errors.append("Real absolute host path in public snapshot")
    if re.search(r"\b(?:\d{1,3}\.){3}\d{1,3}\b",serialized):errors.append("Potential host address in public snapshot")
    forbidden={"patient_id","patient_uid","donor_id","sample_reports","verified_process_cmdlines","command","password","access_token","api_key","environment"}
    def walk(value):
        if isinstance(value,dict):
            for key,item in value.items():
                if key in forbidden:errors.append("Disallowed detail field: "+key)
                walk(item)
        elif isinstance(value,list):
            for item in value:walk(item)
    walk(s)
    return errors

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("snapshot",type=Path);a=p.parse_args()
    errors=validate(json.loads(a.snapshot.read_text(encoding="utf-8")))
    for error in errors:print("ERROR:",error)
    print("Runtime snapshot validation:",len(errors),"error(s)")
    raise SystemExit(bool(errors))
if __name__=="__main__":main()

# V8 evidence-system upgrade — continuation

## Root-reviewed progress

- Preserved frozen V7 checkpoint SHA256
  `e4fda365b2bfa9ea6f82b9b48c4ce6009a0bf31dcfb5a57e7e94cab6d8d3981e`.
- Added fail-closed layered identity/state/malignancy schema and active adapters.
- Found that legacy unversioned `annotation_confidence` differs from final
  `annotation_confidence_v7` for 1,417,900 of 4,676,787 cells. Added an explicit,
  lossless schema bridge rather than silently reusing or overwriting the field.
- Root active test suite passed 48/48; active malignancy/mapper/Numbat regression
  passed 18/18. Historical development modules are retained but not active.
- Ran the guarded three-sample metadata pilot: 14,013 cells, 3 samples, schema
  PASS, exact V7 field/index preservation, gzip round-trip PASS. Wall time 44.57 s;
  peak RSS 7,586,676 KiB. Output SHA256
  `5a1a6eecdecb8d7e158571ea8f791cb527a5eea80c37711e83e673524c3bfa0f`.
- Pilot status is development-only, not independent accuracy and not a V8 freeze.
  Full expansion remains blocked pending scalable malignancy serialization,
  expression-state pilot, independent reference/truth, classifier availability,
  and orthogonal substrate validation.

## Key local evidence

- `${PROJECT_ROOT}/Project_v3/v8_upgrade/reports/05_METADATA_PILOT_RESULT_CN.md`
- `${PROJECT_ROOT}/Project_v3/runs/20260911_v8_metadata_pilot_v1/`
- Active integration entry point:
  `${PROJECT_ROOT}/Project_v3/v8_upgrade/integration/run_migration_v5.py`

V7 and formal CoVarNet/LIANA results were not modified.

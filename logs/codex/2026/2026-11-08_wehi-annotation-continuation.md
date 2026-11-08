# WEHI annotation continuation — 2026-11-08

## Request and initial state
Continue the documented new-TME annotation workflow. User requires server-first storage. R3 handoff arrived with a verified receipt at 00:56:17 UTC / 08:56:17 UTC+8 / 11:56:17 Melbourne. No new expression matrix was downloaded locally. The prior handoff request is historical; its missing inputs are now available.

## Plan and inputs
Read project rules, annotation-v2 methods/run records and current WEHI handoff README/STATUS/semantics. Use a new receiver-side continuation directory; preserve received files and original labels. The isolated repository was fast-forwarded to 4d86aa09 before adding this task. The main local checkout's unrelated changes were not touched.

## Observed input acceptance
At 06:10:38 UTC / 14:10:38 UTC+8 / 17:10:38 Melbourne, 159 manifest entries passed fresh SHA256 checks (README and the manifest itself are outside this internal list; transfer receipt lists 161 files). All 11 obs/embedding exports match by exact row and cell ID. R2 review decisions replay against R3 labels for all nine pools; all five seed files contain ten complete columns matching their source objects. Full expression values were not reread. Existing transfer records have unique IDs, valid pool membership and consistent confidence/status values.

## Implementation and decisions
- Evidence-based working review covers all 22 Myeloid and 19 Endothelial R3 clusters. 388 Myeloid and 1,016 Endothelial cells are set aside for identity/context review; no additional cells are deleted. A high doublet score alone is insufficient for a doublet action.
- Current CLEC4G/CLEC4M-enriched EC is cluster 12, not cluster 6 in the stale task summary. Cluster 6 is ACKR1-enriched.
- Consensus uses all reviewed labels, ten seeds, majority mapping, exchange denominator original-label cell count × 10, bidirectional mean >=0.25 and transitive union; stability >=0.7. Majority ties have deterministic lexical resolution and are logged; tied cell modes cannot enter main proportions.
- Prespecified union can mix state/sample-specific labels with identity. Such components are preserved but excluded from main identity analysis, together with low-depth, nonretained and other flagged cells. No post-result threshold tuning.
- Transfer uses the exported same-order X_emb, k30, seed2026110704, inverse distance, confidence rounded to three decimals then >=0.6 as the existing script. Four threads; matching pynndescent0.5.13/scanpy1.11.5/numpy1.26.4/pandas2.3.3. Microglia uses collapsed Microglia vs Mac_glioma references including the documented moved-to-Myeloid reference clusters. Predictions do not silently reassign the lineage pool.
- Fine-resolution T/NK labels are recorded as observed candidates; resolution-1 stability does not validate those finer splits. They do not overwrite main proportions. This is a conservative deviation from the proposed automatic rare-label overwrite, motivated by the untested finer resolution and TRDC/NK ambiguity.

## Checks and failures
Four synthetic consensus checks passed: identity preservation, documented exchange denominator, logged majority ties and unrelated-component preservation. A separate full-data pandas crosstab replay is running. Local py_compile attempted to write to the system cache and failed under sandbox; replaced with in-memory compile, which passed. Scientific execution was unaffected.

## State and limits
Five-pool consensus computation is complete. Low-depth transfers and the full working sidecar are in progress; final completion will be recorded below after readback. This is not a final scientific annotation release. Four tissue pools have no ten-seed run. E1 October comparison and CoVarNet sample eligibility remain downstream of an accepted annotation scope. Patient rows, matrices and internal host paths stay on the server.

Suggested commit message: Continue verified WEHI R3 annotation with consensus and flagged transfer

## Publication gate
Automatic review rejected the attempted combined commit/public-push command, requesting explicit approval for this particular public batch and per-file sensitivity assurance. No public push is claimed. User confirmation was requested; scientific server work continues independently. Full repository validation completed with 0 errors and 1 pre-existing internal-path warning; new scripts/records passed targeted path/credential-pattern and link checks.

## Working sidecar acceptance
At 2026-11-08 06:22:10 UTC / 14:22:10 UTC+8 / 17:22:10 Melbourne, all 3,347,659 unique IDs and exact 1,387,936 main-candidate memberships passed readback. Six pools contain 535,749 transferred and 122,572 unresolved low-depth cells, all outside main candidates. The 48,210,892-byte sidecar SHA256 is d25fcbcbbdb75737dd011922546d6704cc5a0fef28aa5c922643113781912097. Source data were not overwritten. Biological final acceptance remains false; historical-field/transfer-value provenance replay is recorded separately.

Final provenance readback passed: all 11 historical label columns match for all 3,347,659 cells; all 658,321 six-pool transfer records match labels/confidence/status, with no silent lineage reassignment. Public sync remains blocked pending the requested specific authorization.

## Specific publication authorization
The user explicitly approved publishing this task’s code and de-identified records to XX-114514/Pancancer_log on 2026-11-08. The completed batch includes the eight initially prepared scripts plus the final provenance validator, aggregate acceptance JSON and methods/run logs. No per-cell data, matrices, patient records, credentials or internal server paths are included. The earlier automatic-review rejection is retained as history; publication is now authorized.

# WEHI R3 annotation continuation

The verified compact handoff enables continuation on the receiving research server without moving expression matrices locally or rerunning the original embeddings. This implements the [lineage-rounds method](tme-annotation-v2-lineage-rounds-20261107.md); the outputs are a working annotation, not a final scientific release.

## Input and review
The handoff's 159 internal manifest files were rehashed; all 11 obs/embedding exports match by cell ID and row. R2 decisions replay exactly against R3 labels. Five ten-seed files are complete. The 22 Myeloid and 19 Endothelial clusters receive explicit evidence-based working actions; 388 and 1,016 cells respectively remain set aside. Doublet scores alone do not remove cells. Current CLEC4G/CLEC4M enrichment is EC cluster12, correcting the stale cluster6 note.

## Consensus definition
For each seed, map every seed cluster to its majority reviewed label using all cells. Tied cluster majorities use lexical order and are logged. For reviewed labels a,b, let p(a,b) be the fraction of all a-labelled cell-seed records mapped to b (denominator n(a) × 10). Merge when (p(a,b)+p(b,a))/2 >=0.25, using transitive components. A cell's stability is its modal component's frequency across ten seeds; >=0.7 is stable. Tied cell modes are excluded from main proportions.

This literal prespecified rule can merge identity labels with state, low-depth or sample-specific labels. Such components are retained in diagnostics, but not accepted as main identity subtypes. No threshold was tuned after looking at outcomes. Main candidates also require retained pool membership, keep/keep_tissue action and no round-1/round-2 low-depth flag. Majority mapping and threshold/union mechanics are not independent biological validation.

## Transfer
Use exact exported X_emb with k30 inverse-distance votes, seed2026110704, confidence rounded to 3 decimals then >=0.6, matching the existing script. Reference is keep/keep_tissue without low-depth flags. Validation holds out the lowest-depth 10% within each reference label; reported accuracy is an upper-bound diagnostic because these cells remain deeper than real queries. Microglia references additionally include documented Mac_glioma clusters moved to Myeloid, collapsed into two reference identities. Its predictions do not silently change original pool membership. All inferred cells remain outside main subtype proportions.

## Working sidecar and boundaries
All 3,347,659 original IDs and historical annotation columns are preserved. Latest action, state, consensus component/stability, inference confidence, inclusion flags and pending states are separate fields. Transferred/low-depth/QC/state/sample-specific records cannot enter main identity proportions. Four tissue pools have no ten-seed consensus and are explicitly labelled as such.

Fine-resolution T/NK candidates are recorded only for stable observed cells. Their resolution-1 consensus does not validate a resolution-3 split; candidate names therefore do not overwrite the main subtype. This is a conservative deviation from automatic rare-label overwrite, and cluster42 is not forcibly called NK on TRDC/proliferation evidence. E1 comparison needs the exact October cell-ID/label authority; CoVarNet eligibility requires a separately accepted scope, including the prior uncertain-cell rescue policy. Neither is inferred from current counts.

## Reproducibility
Entry points in `scripts/annotation_v2_resume_20261108/` take HANDOFF_ROOT and NEW_RUN_ROOT. Run audit_handoff_v1.py, make_review_v1.py, consensus_labels_v1.py, audit_consensus_v1.py, transfer_exports_v1.py (Microglia, Endothelial, Myeloid), then assemble_working_v1.py. Run directories reject existing stage outputs. Transfer used four threads with pynndescent0.5.13, scanpy1.11.5, numpy1.26.4 and pandas2.3.3. Actual stage acceptance and limitations are in the run/task records.

Review prose clarification: EC cluster13 action keep_tissue is authoritative and includes that BBB-like identity in the eligible transfer reference. The initial reason's phrase outside pan-cancer reference was stale prose, corrected in a separate server clarification without changing actions or computed results.

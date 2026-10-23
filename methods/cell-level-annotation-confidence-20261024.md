# Common cell-level confidence for annotation-order comparison

This stage follows the accepted [source-lineage discovery](source-native-routing-and-lineage-discovery-20261024.md) and native review. Its input is the frozen pilot, first seed and resolution1. It applies exactly the same evidence policy to A/B/C; graph memberships, native counts and historical labels remain frozen. Full-scale workflow selection is a later gate.

## Evidence and decisions

1. Re-read the frozen native marker panel, including the specified red/neutrophil anchors. Reproduce every accepted source-group marker fraction and the previous native-family mask. Preserve source measurement masks and author-field provenance.
2. Use the accepted new native review's unique reference proposals, with the same complete coarse-family dictionary. Immune-reference nonhematopoietic labels are ignored. Record any coarse reference-mask change relative to the earlier routing; changed prediction scope is not hidden.
3. Cross-family native support, cross-family reference conflict, or native/reference family disagreement yields unresolved identity. T versus NK or B versus plasma disagreements retain a broader family and explicit within-family conflict. Zero RNA remains unresolved.
4. Native-only coarse support is retained with its provenance. Reference-only evidence requires a compatible group context; marker dropout remains an inferred evidence tier. Missing reference support is not negative expression evidence.
5. A coherent minority cell with its own fine-marker support and a matching unique fine reference retains its candidate identity even when the group has another label. Ancestor CD4 support is distinguished from Treg. Multiple incompatible direct fine candidates fall back to coarse identity.
6. Group fine proposals can be supported by the cell's native fine evidence or compatible reference evidence, with direct and inferred tiers separate. Incompatible specific fine references trigger coarse fallback. A gamma-delta RNA identity candidate always requires the cell's own TRDC/TRGC/CD3 evidence and retains the orthogonal-review requirement.
7. A branch without enough cells for a fine graph retains that no-fit flag. Strong direct cell evidence can still propose a rare identity; it never creates a claim that a fine graph was fitted. Group-only inference cannot promote such a branch.

These are annotation evidence levels, not calibrated confidence probabilities or independent truth. Reference training overlap remains unresolved, and routing/group context can encode the RNA programs being studied. Coverage is reported without an accuracy claim. Clinical outcomes are not inputs.

## Deliverables and invariants

Source-private categorical sidecars contain candidate consensus major/type, coarse/fine evidence, group disagreement flags, native/reference support masks, parent discovery-group codes and branch exposure. Native row coordinates and source-qualified cell-axis hashes bind them to the original cells. Metadata retains the label version, author provenance and old-fine crosswalk counts. Original author labels remain in their native backing store.

Newly guarded subsets do not inherit a false own-group marker claim: current marker pointers explicitly refer to the parent discovery group, and subsequent chosen-workflow annotation must compute subset markers or retain that parent-origin label. No new inferential DEG is claimed.

The actual full pilot checks that hard cell conflicts cannot be overwritten by any method, direct fine evidence is identical across methods, gamma candidates have own RNA support, all branch flags are preserved, all source/group/native axes agree, and every serialized array reads back identically. The initial nine counterexamples cover minority preservation, cross/within-family conflicts, dropout, rare small branches, gamma evidence, group disagreement, zero RNA, vocabulary and reference scope. The source-side independent audit then checks stored decisions and denominator/crosswalk completeness.

The previous turn completed authoritative computational work and is classified as progress. This stage advances cell-level annotation and crosswalk deliverables; it does not narrow the full eight-part objective.

## Version 2 correction and actual acceptance

The first producer completed successfully. Review then tightened a reference-only edge case: when native family markers are absent and the group suggests T but the unique reference suggests NK (or B versus plasma), version2 retains the broader family with an explicit disagreement. It does not treat shared broad family as agreement on the specific type. Inferred fine evidence is named `reference_group_inferred_without_complete_native_support`; this does not assume every missing support bit is caused by dropout. A tenth counterexample checks this boundary. Version1 is preserved.

Version2 and its independent audit completed with actual child exit0. The auditor does not call the producer decision function: all-cell coarse guards and direct fine candidates are reconstructed independently, version differences are confined to the stated edge case, all parent denominators/crosswalks/marker origins are checked, and one size-selected complete source has fine-marker Booleans rebuilt directly from native CSR spans. Computational acceptance does not promote reference-derived candidates to biological ground truth or final atlas labels.

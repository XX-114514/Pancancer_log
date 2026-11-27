# Bilingual preprint v2 anatomical figure and package

## Request and starting state

User requested removal of the resource workflow from the main opening figure, an anatomical/cancer-coverage redesign, aligned bilingual title and authors based on the current English title, and a separate v2 folder plus archive in the supplied local publication directory. Prior outputs and unrelated dirty worktree changes are preserved. Read project instructions, status/roadmap/task board, runtime synchronization rules and the preceding figure-redesign task log. Inspected only current manuscript files and relevant frozen figure inputs.

Audit time: 2026-11-27T03:32:19.901612+00:00 UTC; 2026-11-27T11:32:19.901612+08:00 UTC+8; 2026-11-27T14:32:19.901612+11:00 Melbourne.

## Implementation

Created a four-panel opening figure: original vector anatomical coverage, major-compartment UMAP, numbered fine-identity UMAP and native-RNA marker support. Removed the workflow from the main figure. Coverage names and counts derive from the supplied source registry, not the reference screenshot. The anatomy is schematic; the figure does not infer tissue locations from these data.

Frozen coordinates, every reference point, category labels, numbering and palette are retained. Corrected label placement to resolve axes aspect before display-space collision measurements. The first draft had two collisions after shrinking the UMAP panels; aspect-aware placement removed both. Adjusted spacing around three-line anatomical labels. Trial figures are retained privately.

The supplied English title had changed since the preceding task; used the current title as authoritative and translated it into Chinese. English author placeholders were completed from filled Chinese fields, with names and order preserved; no missing affiliations or contributions were invented. Updated the atlas introduction, Figure 1 caption and relevant subpanel references in both manuscripts. Remaining author-check fields are retained and listed separately.

Assembled bilingual editable Word, rendered PDF and text counterparts; standalone figures; original supplementary tables; a fine-identity number key; bilingual main-table exports; opening-figure aggregate source tables; plotting/assembly code and editorial receipts. No raw matrices or cell arrays were added. Existing legacy plotting scripts retain their external server dependencies, which are documented. No scientific models or outcome tests were rerun.

## Verification and delivery

Both manuscripts render to 21 pages. Visually inspected all page PNGs and checked the changed title/figure pages, CJK glyphs, captions, table wrapping and page boundaries. All non-edited paragraphs and every original table cell are preserved. Other original figure files and supplemental tables are byte-identical, and embedded Figure 2 through S1 image content is preserved. The fine-identity key is byte-identical to the previous frozen key. UMAP label audit reports zero overlapping pairs. Source registry totals, source/cancer counts and marker denominators reconcile through the renderer assertions.

Archive CRC and per-file SHA checks pass, including after copying the requested folder/archive into the destination and reading them back. Private archive SHA256: `fafdedb34901fac70a12b4a873b04f73ae2d51c20b2c1f849817fc145796c077`.

Full repository validation reports one pre-existing broken-link error group (two targets) and one existing path warning, unchanged from the preceding task. These unrelated links were not edited. Scoped AST and whitespace checks apply to the explicit new paths. Public payload is restricted to generic plotting code, README and this sanitized technical log. No manuscript, figure, names, private tables, arrays or learned weights are staged for public Git.

## Failures and boundaries

A first local preservation check included a ZIP directory entry as if it were an image; corrected the directory-entry filter and reran successfully. Initial Git fetch in the restricted network failed DNS; an authorized network fetch succeeded. No scientific interpretation was revised beyond introducing the new figure. This is a local author-review package, not an online submission. Remaining AUTHOR CHECK fields still require completion.

Next: author review of the bilingual v2 package and completion of manuscript front matter and citation/accession checks. Suggested commit: `Package anatomical opening-figure renderer for bilingual preprint v2`.

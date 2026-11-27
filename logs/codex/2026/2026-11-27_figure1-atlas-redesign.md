# Redesign of the opening atlas figure

## Request and initial state

User supplied three reference figure screenshots and requested a redesigned opening figure using actual project data and UMAPs. Read the supplied visual layouts as design references, not as instructions or sources of project results. Reviewed current manuscript resource counts, its source registry, frozen identity-reference inputs, existing label helper, project narrative/rendering rules and repository state. Unrelated worktree changes and original manuscripts/figures are preserved.

## Design and implementation

Completed audit time: 2026-11-27T14:06:07.575714+11:00 (Australia/Melbourne).

Built a five-panel English figure: circular cohort coverage; linked frozen resources; major-compartment UMAP; numbered fine-identity UMAP; native-RNA annotation dotplot. Adapted the reference examples' coverage overview, large real embedding panels and resource-to-identity reading order. No reference-paper artwork, data, pseudo-time or mechanism was substituted. Proposed bilingual captions and layout notes remain private with the figure package.

Coverage uses aggregate counts parsed from the current source registry; its cell/label/source/unit/patient totals reconcile to the supplied manuscript. Broad anatomical/disease display grouping is explicit in the private source table and preserves every original cancer-label contribution. Patient counts retain registry semantics rather than being presented as independently re-verified patients. Reference UMAP scope is distinct from the complete harmonised atlas; the underlying reference arrays, class assignments, category order and palette remain frozen. All reference points are rendered; no coordinate refit or downsampling occurred.

Source data are externally supplied to the generic plotting command:

```sh
python scripts/figure1_redesign_20261127/render_figure1.py --atlas-dir "$ATLAS_DIR" --node-map "$NODE_MAPPING_TSV" --registry "$SOURCE_REGISTRY_TSV" --resources-json "$RESOURCE_CONFIGURATION" --output "$NEW_OUTPUT"
```

The RNA panel uses eligible-cell-weighted identity means and RNA-positive fractions, retaining zero-expression cells and existing source–gene qualification. It is annotation support, not a differential-expression analysis. Saved PNG/PDF/SVG, a separately numbered identity key, source tables and code. Main figure dimensions are 180 x 225 mm. Scatter layers are rasterized; labels, diagram and dotplot are editable vector elements.

## Verification and correction

Executed with Python 3.12, NumPy 2.5.5, pandas 2.2.3 and matplotlib 3.10.8 using the existing local plotting dependencies and a writable Matplotlib cache. Original draft exposed overly long resource titles, a crowded compartment-label pair and a colorbar title too close to the top tick. Wrapped resource titles, adjusted label offsets and moved legend text, then rerendered into a separate version directory. Both drafts are retained.

Observed checks pass: registry sums and coverage labels reconcile; array cancer totals match the frozen source/cancer table; complete reference points and identities are retained; all fine-identity labels occur once in the SVG; historical identity numbering matches exactly; direct-label text has zero overlap pairs and no checked text lies outside the figure. Independently recomputed three parent/gene RNA aggregates with Python CSV arithmetic, including eligible and positive counts, positive fraction and weighted mean. Final main figure and companion key visually inspected. Output SHA manifest and ZIP CRC pass. No new scientific fitting or statistical tests occurred.

Final private archive SHA256: `4084520429292a8402f80b210c3cd3d031d9119b418d955277c4045051d3e4da`.
Renderer SHA256: `756742e807b22e252d8149e30cddcc674ff60d4997dfb2c5b887d83c1a5f5d49`.

## Publication and scope boundary

Only generic plotting code and this sanitized technical log are intended for the public repository. Figure data, scientific artwork, registry, coordinates, expression summaries, bilingual captions and private receipts remain outside public Git. No participant rows or learned weights are exported in the lightweight private figure delivery. The original manuscript's Figure 1 is not replaced automatically. Full repository validation retains the previously recorded broken-link error group and existing path warning; scoped AST/whitespace, real-input rendering and artifact checks pass.

Next: author review of the new layout and proposed legend before inserting it into the manuscript. Suggested commit: `Add source-traceable atlas opening figure renderer`.

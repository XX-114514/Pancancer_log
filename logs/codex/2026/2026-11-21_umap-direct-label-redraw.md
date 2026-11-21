# UMAP direct-label redraw

## Request and initial state

User requested actual redraws of Figure 1a (eight major identities), Figure 1b and S06 (numbers 1–44 aligned to the S06 key), and S07a (resolved major identities labelled; unresolved/boundary identities retained in the legend). Reviewed the homepage renderer at its Figure 3 UMAP block and the user-specified branch. Existing unrelated dirty workspace files were retained. Private plotting inputs and rendered scientific figures remain outside public Git.

## Version and method

Source: `codex/0906-umap-labels` at `274ed4463af8c048bc2daa9e48a46315ed35784a`, directory `scripts/launch_preprint_20261120/`. The renderer/helper files match the branch blobs; a pinned adapter selects the requested complete AST statement blocks. The adapter checks renderer SHA256 before using source line ranges. It preserves coordinates, membership, palettes, shuffle seeds and the native-RNA dotplot. Label anchors use category coordinate medians with deterministic collision avoidance and leader lines. No model fitting or new statistical tests occurred.

Executed on the existing research-server plotting environment (Python 3.10, NumPy 1.26, pandas 2.3, matplotlib 3.10). Local bundled Python lacked matplotlib, so actual rendering used that verified environment. Command structure:

```sh
python code/redraw_requested.py --inputs "$PORTABLE_INPUTS" --output "$NEW_FIGURE_OUTPUT"
```

## Execution, correction and observed result

Evidence/completion audit time: 2026-11-21T10:39:36.568534+00:00 / 2026-11-21T18:39:36.568534+08:00.

All three target figures rendered in PNG/PDF/SVG. During visual review, the initial S06-only block selection omitted its new canvas and reused the Figure 1 canvas. The faulty trial was retained on the server, excluded from final delivery, and the adapter corrected to include the complete S06 canvas and labels. S06 was separately rerendered and the corrected output visually reviewed. Original executed code is retained separately from the normalized, pinned replay code. Frozen baseline and Claude-returned publication were not overwritten.

Observed checks passed:

- Figure 1 has all eight major labels; Figure 1b has exactly one instance of every number 1–44.
- S06 has its own full-identity canvas, all 44 numeric labels and all 44 numbered identity-key entries. S06 and Figure 1b derive numbering from the same frozen node order.
- S07 has direct labels for resolved identities; `Unresolved` and `T_NK_boundary` occur only in the legend.
- All three final PNG previews were visually inspected, including corrected S06; no clipped labels or overlapping identity legend were observed.
- Output hashes, local-versus-server frozen input hashes, archive SHA256 and ZIP CRC passed.
- Final private archive: 15,444,754 bytes; SHA256 `90675b12bca5057d8cdcf17b175357e3fdf458c2c60e3190ae486473d7d088b9`.

## Changed public files and checks

Updated the two portable renderers and README from the requested branch; added its `umap_labels.py` helper plus `redraw_requested.py`. Public files contain generic code/documentation only. Identity keys, coordinate arrays, scientific figures and private receipts are not staged for publication. AST parse and scoped diff whitespace checks pass. Full workspace validation reports the same two pre-existing broken relative links (one validation error group) and existing internal-path warning; unrelated historical files were not edited to resolve them. The isolated sparse checkout does not contain the full repository validator, so the validator was run from the main local checkout.

## Scope and next step

This delivery updates the requested portable figure set. It does not substitute pages into the separately returned Claude manuscript or claim a fresh biological/statistical audit. Delivery includes the three figure sets, replay code, numeric/immune label keys and technical validation. Durable private copy remains on the research server; a current-task local download copy is provided. Suggested commit: `Render frozen atlas UMAPs with direct identity labels`.

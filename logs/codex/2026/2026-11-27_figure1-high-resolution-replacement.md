# High-resolution opening figure replacement

Audit time: 2026-11-27T03:49:47.108172+00:00 UTC; 2026-11-27T11:49:47.108172+08:00 UTC+8; 2026-11-27T14:49:47.108172+11:00 Melbourne.

User supplied an updated anatomical figure bitmap and requested improved clarity, PDF export and replacement of the existing opening figure. Preserved the latest local manuscript edits and backed up the complete current delivery and archive before replacement. Prior anatomical-package task records supply the unchanged frozen input scope and numerical checks.

The supplied bitmap is used only for its anatomical artwork and original leaders. Upper-panel text is re-typeset as vector text. Both UMAPs are regenerated from every frozen reference point at 600 dpi, with the frozen number/identity map rather than transcribed bitmap numbers. Dotplot and labels remain vector elements. No inferred data or image-generation model replaces scientific plots. Anatomical artwork remains at native raster resolution; it is not described as a fully vector anatomical illustration or as recovered high-resolution detail.

Exported a hybrid PDF with embedded fonts, a 600 dpi PNG and hybrid SVG. Source tables and the identity key remain byte-identical. Replaced only the first image part in each latest Word package: all other ZIP members, including document XML and current user text, are byte-preserved. Located the image dynamically because the user had edited the English manuscript since the preceding task. Both updated Word files were rendered. All previously inspected page images are pixel-identical except the changed figure pages and an English title page containing retained user edits; these changed pages were visually inspected.

Inserted the native figure PDF at the re-rendered manuscript image transform to preserve sharp vector text in both reading PDFs. The first overly strict aspect-ratio check rejected the user's slightly altered English image size; inspected the actual placement, allowed the documented small rounding difference and preserved its exact transform. Both PDFs remain 21 pages. Native PDF figure pages were rendered and inspected after insertion.

Checks: scoped Python AST; zero UMAP label-overlap pairs; extractable cancer/gene labels in standalone and manuscript PDFs; native PDF dimensions and image-layer inspection; frozen source-table hashes; DOCX non-target member byte equality; final per-file SHA and archive CRC after destination readback. Archive SHA256: `5e874325eac08b168ec5b7f10120c828baa26ba290a58b16b0cc513f494935bc`.

Only generic renderer code and this sanitized technical log are staged for public Git. Artwork, manuscripts, labels-as-source tables, private receipts and archive are kept outside public Git. No new scientific fitting or statistical testing occurred. Full repository validation in the immediately preceding package task reported the existing two broken links as one error group and one warning; this figure-only change introduces no documentation links and uses scoped checks instead of repeating that unrelated scan.

Next: use the replaced native PDF for figure submission. For further anatomical-detail improvement, an original editable/high-resolution anatomical asset would be required; current typography and scientific plots already use the improved rendering path. Suggested commit: `Add hybrid high-resolution opening figure renderer`.

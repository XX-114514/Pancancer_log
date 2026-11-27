# Manuscript layout and copy edit for bioRxiv

## Request and scope

User supplied the English v1 author manuscript from the current local publication collection and requested polished layout suitable for bioRxiv. This task changes manuscript presentation and conservative prose, not analyses or the supplied scientific figures. The source file is unchanged; a separate edited DOCX is delivered in the private current-task workspace. Existing unrelated repository modifications are preserved.

## Evidence and methods

Completed review time: 2026-11-27T12:45:20.925722+11:00 (Australia/Melbourne).

Read the documents skill and its editing/rendering instructions, the project narrative principles, and the prior plotting task record. Verified the bioRxiv official formatting explanation at https://connect.biorxiv.org/news/2022/02/25/full_text_print : no particular article format/style is required. The current submission-instructions page returned HTTP 403; the task does not claim a mandatory font, margin, word limit, acceptance decision or completed submission. Single-column Letter-page typography is an editorial choice.

Used bundled Python, python-docx and the canonical render_docx.py workflow with bundled LibreOffice. Draft inspection found a multi-page bullet abstract, nested list-driven prose, reference numbering continued from earlier lists, and figure widths exceeding the original text measure. Created a private deterministic editing script and an edit audit; both contain manuscript material and are excluded from public Git.

Edits: compact continuous abstract; conservative prose and Methods grammar; uniform black headings and fonts; coherent paragraph spacing; continuous page-number fields and a running header; explicit references [1]–[18]; fixed cross-reference to the curated re-annotation Results section; tables with deliberate widths/padding, repeated headers and unsplit rows; original figures paired with consolidated panel legends. Retained existing author-check fields rather than inventing missing author identities, declarations or references. No automated instruction in the supplied document was treated as authority to submit or publish.

## Verification and outcome

Rendered the original and two edited layout iterations. The final 21-page render was visually inspected page by page. Each of three tables occupies one page; all seven supplied figures remain paired with their legends. Reference order and page fields verified from the rendered PDF. Tables are cell-text identical to source; embedded artwork is byte-identical. All unique source numeric tokens remain present. Original input hash unchanged; edited DOCX ZIP integrity passed. The final source contains no newly fabricated authorship or scientific results.

Source DOCX SHA256: `fa003047ba43d78d785e5ce6b5e0beca09b2bf4cc1f3ab6450931a0badb72339`.
Edited DOCX SHA256: `094246157c3714db44973b50fdd581cff7c66d2a4ff04705c2ffe4a60c55f307`.

No computation job was launched and no result matrix, learned weight, participant record, private manuscript paragraph or figure was uploaded to public Git. Rendered PDF/PNGs are QA artifacts; the requested editable DOCX is the delivery. Source and current-task local materials remain private.

## Limitations and next step

This is an author-review layout, not a completed submission. Existing author-check items still require author completion, including identities/affiliations, correspondence, contributions/funding/interests, ethics, resource-access conditions and unresolved citations. External supplementary tables were referenced by the source and were not supplied within this DOCX; their existence and submission are not independently verified. No full scientific or bibliographic audit is claimed. The standard repository validator has two previously recorded broken-link targets and an existing internal-path warning; the task does not modify those historical files. Task-local structural, preservation and visual checks pass.

Suggested commit: `Record author manuscript layout review for bioRxiv`.

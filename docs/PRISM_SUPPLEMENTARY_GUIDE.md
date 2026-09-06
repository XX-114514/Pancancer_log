# Prism and supplementary-material guide

## Authoritative inputs

- Cell identity and malignancy: `releases/annotation_v7_20260905`.
- Current communication evidence figures: 
  `releases/communication_evidence_chain_20260906`.
- Historical rationale and failures: `runs/`, `logs/` and `decisions/`.

The communication figures are based on the V5 downstream CoVarNet/LIANA branch.
V7 supersedes identity/CNV for new analyses but does not retroactively make those
communication results V7-derived.

## Prism handoff

1. Import the relevant aggregate TSV from the same release as the figure.
2. Preserve the release ID, table filename, grouping columns and denominator in
   the Prism project notes.
3. Use PDF files for vector panel assembly; PNG files are review copies.
4. Do not manually transcribe cell counts when a tracked TSV exists.
5. Export a final Prism data table and a panel-to-source index before submission.

Suggested panel index columns:

```text
figure_id  panel_id  release_id  source_table  filter  statistic  denominator  script  notes
```

## Supplementary package

Include methods, aggregate source tables, figure PDFs, a README, software versions
and external-object checksums. Exclude raw/cell-level data unless data governance
and journal policy explicitly permit release through an appropriate repository.

## Mandatory version wording

Recommended methods wording:

> Cell identities and CNV-derived malignancy evidence use the frozen PanCancer
> annotation V7 reference (run 20260905_v7_cnv_rerun_freeze_v1). Communication
> analyses shown in the 20260906 evidence-chain snapshot reuse the separately
> frozen V5 CoVarNet/LIANA branch.

Any future V7 rerun of CoVarNet/LIANA should receive a new release ID and should
not overwrite the current communication snapshot.

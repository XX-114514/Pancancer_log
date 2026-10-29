# Multiomics storage accounting

## Request and initial state

The user asked to calculate the storage occupied by the data in the project-wide multiomics workbook. Continue the existing source-inventory task. Read spreadsheet and scientific-research requirements, the existing workbook builder/input, source catalogs and prior storage-maintenance evidence. Preserve unrelated checkout changes and existing scientific objects.

## Measurement and evidence

Resolve the project-root alias to its canonical path. Traverse filesystem metadata without following symbolic links or reading expression, protein, coordinate or individual clinical payloads. Count allocated blocks as st_blocks times 512 and logical file lengths separately. Deduplicate repeated device/inode pairs, including shared hardlinked array storage. Keep parent-directory totals, nested source groups and per-source references distinct. Do not estimate source bytes from cell proportions.

An independent GNU du project-root query exactly matched the allocated-block total from the completed metadata snapshot. Complete source-group sizes, project-version partitions, timestamps, private paths and source storage statistics remain outside the public repository. The snapshot is an observed filesystem accounting, not an atomic frozen filesystem or a separate audit of filesystem reflink/compression behavior. Code, logs, historical matrices, compressed originals and actual retained processing copies are included in the project-root total. RAM-expanded objects and external public repositories are not charged as local persistent files.

## Inventory corrections

The broader metadata scan revealed legacy spatial and later extension directories missing from the earlier selective source table. Read the existing source catalog and download catalogs to append those sources with their actual storage scope and unresolved analytical admission. The existence of a downloaded directory does not establish complete download, patient independence, clinical eligibility or a new analysis result. A previously uncertain platform label was resolved from the project's explicit download catalog. Preserve previous rows except that declared platform/provenance correction.

## Workbook edit and verification

Use the bundled Artifact Tool to edit the existing workbook in place, append the newly identified source rows and add storage-group accounting and source-to-storage correspondence views. Retain a private byte-identical pre-edit copy for verification. The original sheet values, column widths and freeze panes were checked against that copy; only the two declared platform/provenance cells changed. All appended values, every storage binding, displayed totals, formula caches, table filters and archive CRC were independently checked from the saved OpenXML.

Recalculate the workbook and visually review the new storage views and appended source rows. In a disposable in-memory import, change a storage-byte input by exactly one GiB, verify the owning group and linked shared summary update by one GiB, verify independent project/space partitions remain unchanged, then restore. Do not write this test to the delivered file. Unknown storage remains missing rather than zero, and shared references do not become repeated independent costs.

## Failures and recovery

The large plain metadata stream timed out after its server scan had completed. A second metadata snapshot used a compressed report retained temporarily on the server; that longer stream also timed out. Retrieve a compact, gzip-validated selection from the completed saved report using a fresh read-only connection, without repeating scientific analyses or changing network configuration. Keep incomplete transport files as failures rather than treating them as complete receipts. Both failed transport sessions terminated. A local preparation patch initially used the skill working directory instead of the workspace; fix its absolute path before workbook authoring. The required edit-operation marker ran successfully once.

## Changes and boundary

Updated the private workbook and added private metadata, mapping and verification files. The only server write was a small temporary compressed filesystem-metadata receipt, with no raw data payload. No data was deleted, compacted, renamed or scientifically reprocessed. Added only this technical task log and a daily entry to the public record repository; no private storage table or source payload is published.

Repository validation, normal technical-log synchronization and remote verification are bound to the private final delivery receipt. Follow-up work should use the recorded stable storage groups when assessing storage, preserve shared masters, and distinguish potential retention decisions from this read-only accounting. Suggested commit: Record multiomics storage accounting and workbook verification.

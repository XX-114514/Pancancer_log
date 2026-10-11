"""Versioned guard for unnamed CSV export-index columns in CosMx conversion.

Apply to a reviewed parent harmonize_highres_spatial module before its CLI run.
This changes gene-column selection only. It does not annotate cells or overwrite
previous converted objects; the caller must provide a new output directory.
"""
import re


def is_export_index_column(name):
    return bool(re.fullmatch(r"Unnamed:\s*\d+(?:\.\d+)?", str(name), flags=re.I))


def install(parent_module):
    previous = parent_module._cosmx_gene_columns

    def guarded_columns(frame, id_columns):
        known_indices = {column for column in frame.columns if is_export_index_column(column)}
        return previous(frame, set(id_columns) | known_indices)

    parent_module._cosmx_gene_columns = guarded_columns
    return parent_module

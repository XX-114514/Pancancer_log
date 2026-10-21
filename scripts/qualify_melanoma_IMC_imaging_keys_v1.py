"""Deprecated initial schema inspection; use the bounded v2 qualifier.

The first executed producer exported all header labels before discovering
that the source also uses coded patient identifiers as wide-table headers.
Its original source and SHA remain in the private failure provenance.
Do not reuse broad schema export. The v2 producer exports contract fields
only and accepts exact partial ROI bridges without inferring missing IDs.
"""
if __name__ == '__main__':
    raise SystemExit('Use qualify_melanoma_IMC_imaging_keys_v2.py; v1 is withdrawn.')

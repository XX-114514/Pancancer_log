"""Exact hemoglobin symbol subset for model features, not expression/QC.

This excludes neither all globin-family members nor hemoglobin pseudogenes.
"""

GLOBIN_SYMBOLS = frozenset({
    'HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ',
})


def is_exact_globin_symbol(gene):
    return str(gene).upper() in GLOBIN_SYMBOLS

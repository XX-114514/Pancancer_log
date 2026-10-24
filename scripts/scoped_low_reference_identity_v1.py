"""Apply the accepted reference-family scope before interpreting fine labels."""
from native_lineage_routing_v1 import reference_family,local_family
from cell_identity_confidence_v2 import low_reference_identity,FINE_MAJOR


def scoped_low_identity(label):
    family=reference_family('immune_low',label)
    if family is None:return None,set()
    major,fine=low_reference_identity(label)
    if major is not None and local_family(major)!=family:raise ValueError('Eligible reference major conflicts with frozen family dictionary')
    if any(local_family(FINE_MAJOR[t])!=family for t in fine):raise ValueError('Eligible fine reference conflicts with frozen family dictionary')
    return major,fine

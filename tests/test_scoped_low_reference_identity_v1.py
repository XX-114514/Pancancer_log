from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import unittest
from scoped_low_reference_identity_v1 import scoped_low_identity


class Checks(unittest.TestCase):
    def test_mixed_progenitor_cannot_infer_mature_erythroid(self):
        self.assertEqual(scoped_low_identity('Megakaryocyte-erythroid-mast cell progenitor'),(None,set()))
        for name in ['Early lymphoid/T lymphoid','HSC/MPP','Cycling cells','Epithelial cells']:
            self.assertEqual(scoped_low_identity(name),(None,set()))

    def test_explicit_myeloid_progenitor_does_not_acquire_mature_fine_type(self):
        self.assertEqual(scoped_low_identity('Neutrophil-myeloid progenitor'),('Myeloid',set()))

    def test_mature_supported_types_keep_existing_hierarchy(self):
        self.assertEqual(scoped_low_identity('Regulatory T cells'),('T',{'Treg','CD4_T'}))
        self.assertEqual(scoped_low_identity('gamma-delta T cells'),('T',{'T_gamma_delta_RNA_candidate'}))
        self.assertEqual(scoped_low_identity('Plasma cells'),('Plasma',set()))


if __name__=='__main__':unittest.main()

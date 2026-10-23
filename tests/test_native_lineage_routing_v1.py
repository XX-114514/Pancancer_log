import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from native_lineage_routing_v1 import FAMILIES,STATUS,reference_family,route_cells


class RoutingTests(unittest.TestCase):
    def run_case(self,native,refs,group,positive=True):
        a=np.array([[f in native for f in FAMILIES]]);b=np.array([[f in refs for f in FAMILIES]])
        code=FAMILIES.index(group) if group in FAMILIES else -1
        r,s=route_cells(a,b,[code],[positive])
        return FAMILIES[r[0]] if r[0]>=0 else None,STATUS[s[0]]

    def test_coherent_minority_not_overwritten(self):
        self.assertEqual(self.run_case(['T_NK'],['T_NK'],'Myeloid'),('T_NK','native_reference_concordant'))

    def test_native_contradiction_retained_unknown(self):
        self.assertEqual(self.run_case(['T_NK','Myeloid'],['T_NK'],'T_NK')[1],'native_cross_lineage_review')

    def test_reference_conflicts_and_disagreement(self):
        self.assertIsNone(self.run_case(['T_NK'],['Myeloid'],'T_NK')[0])
        self.assertEqual(self.run_case(['T_NK'],['T_NK','Myeloid'],'T_NK')[1],'reference_cross_lineage_review')

    def test_native_dropout_requires_context(self):
        self.assertEqual(self.run_case([],['T_NK'],'T_NK'),('T_NK','reference_local_group_concordant_native_dropout'))
        self.assertIsNone(self.run_case([],['T_NK'],'Myeloid')[0])
        self.assertIsNone(self.run_case(['T_NK'],['T_NK'],'T_NK',False)[0])

    def test_native_only_is_provisional(self):
        self.assertEqual(self.run_case(['T_NK'],[],None),('T_NK','native_only_provisional'))

    def test_reference_scope_and_unknown_vocabulary(self):
        self.assertIsNone(reference_family('immune_low','Epithelial cells'))
        self.assertEqual(reference_family('lung','AT2'),'Epithelial')
        self.assertEqual(reference_family('breast','Lymph-immune'),'Endothelial')
        self.assertEqual(reference_family('immune_low','Erythrophagocytic macrophages'),'Myeloid')
        with self.assertRaises(ValueError):reference_family('lung','NEW UNREVIEWED CLASS')


if __name__=='__main__':unittest.main()

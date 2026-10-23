import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from review_lineage_anchor_sensitivity_v1 import anchored_calls


class AnchorTests(unittest.TestCase):
    def test_globins_alone_not_new_anchor_support(self):
        genes=['HBB','HBA1','AHSP'];x=np.array([[1,1,0],[1,1,1]],bool)
        a=anchored_calls(x,genes,np.ones(3,bool),{'Erythroid_related':np.ones(2,bool),'Neutrophil_related':np.zeros(2,bool)})
        np.testing.assert_array_equal(a['Erythroid_related'],[False,True])

    def test_inflammatory_pair_requires_additional_anchor(self):
        genes=['S100A8','S100A9','CSF3R','FCGR3B'];x=np.array([[1,1,0,0],[1,1,1,0],[0,0,1,0]],bool)
        a=anchored_calls(x,genes,np.ones(4,bool),{'Erythroid_related':np.zeros(3,bool),'Neutrophil_related':np.array([1,1,0],bool)})
        np.testing.assert_array_equal(a['Neutrophil_related'],[False,True,False])

    def test_unmeasured_anchor_not_interpreted_as_support(self):
        a=anchored_calls(np.ones((1,2),bool),['AHSP','CSF3R'],np.zeros(2,bool),
                         {'Erythroid_related':np.ones(1,bool),'Neutrophil_related':np.ones(1,bool)})
        self.assertFalse(a['Erythroid_related'][0]);self.assertFalse(a['Neutrophil_related'][0])


if __name__=='__main__':unittest.main()

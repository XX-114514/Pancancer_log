import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
import pandas as pd
from audit_joined_native_atlas_delivery_v1 import assert_carried_axis_equal

class MissingMetadata(unittest.TestCase):
    def test_same_object_missing_labels_pass(self):
        assert_carried_axis_equal(np.array([float('nan'),'T',None,pd.NA],object),np.array([float('nan'),'T',None,pd.NA],object))
    def test_missing_is_not_the_literal_nan_string(self):
        with self.assertRaises(AssertionError):assert_carried_axis_equal(np.array([np.nan],object),np.array(['nan'],object))
    def test_position_and_shape_changes_fail(self):
        with self.assertRaises(AssertionError):assert_carried_axis_equal(np.array(['T',np.nan],object),np.array([np.nan,'T'],object))
        with self.assertRaises(AssertionError):assert_carried_axis_equal(np.array(['T'],object),np.array(['T','B'],object))
    def test_missing_representation_and_real_labels_stay_distinct(self):
        for a,b in [(None,np.nan),(pd.NA,None),('T','B')]:
            with self.assertRaises(AssertionError):assert_carried_axis_equal(np.array([a],object),np.array([b],object))

if __name__=='__main__':unittest.main()

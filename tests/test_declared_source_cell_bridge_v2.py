import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from declared_source_cell_bridge_v2 import RECIPES,encoded_identity,recover_rows
class DeclaredSourceBridgeTests(unittest.TestCase):
    def run_case(self,raw,libs,query_indices=None,query_libs=None,recipe=RECIPES[0]):
        ids,units=zip(*(encoded_identity('SRC',r,l,recipe) for r,l in zip(query_indices or raw,query_libs or libs)))
        return recover_rows('SRC',ids,units,raw,libs,recipe)
    def test_order_and_full_library_identity(self):
        r=self.run_case(['a-1_2','b-1_2'],['L1','L2'],['b-1_2','a-1_2'],['L2','L1'])
        self.assertEqual(r['raw_row_position'],[1,0]);self.assertTrue(r['complete_one_to_one'])
    def test_shared_barcode_in_different_libraries_is_not_merged(self):
        r=self.run_case(['AAAC-1','AAAC-1'],['L1','L2']);self.assertEqual(r['raw_row_position'],[0,1])
    def test_wrong_library_is_unmatched(self):
        self.assertEqual(self.run_case(['AAAC-1'],['L1'],query_libs=['L2'])['raw_row_position'],[-1])
    def test_sample_suffix_preserved(self):
        self.assertEqual(self.run_case(['AAAC-1_2'],['L1'],['AAAC-1_3'])['raw_row_position'],[-1])
    def test_legacy_repeated_prefix_is_deliberately_preserved(self):
        cell,unit=encoded_identity('SRC','3__3_AAAC-1','3',RECIPES[1]);self.assertIn('|cell=SRC__3__3__3_AAAC-1',cell)
        self.assertTrue(recover_rows('SRC',[cell],[unit],['3__3_AAAC-1'],['3'],RECIPES[1])['complete_one_to_one'])
    def test_duplicate_identical_author_records_are_ambiguous(self):
        r=self.run_case(['AAAC-1','AAAC-1'],['L1','L1'],['AAAC-1'],['L1']);self.assertEqual(r['raw_row_position'],[-2]);self.assertFalse(r['complete_one_to_one'])
    def test_archive_duplicate_is_reported_without_deleting_rows(self):
        r=self.run_case(['AAAC-1'],['L1'],['AAAC-1','AAAC-1'],['L1','L1']);self.assertEqual(r['repeated_author_rows_in_archive'],{0:2});self.assertFalse(r['complete_one_to_one'])
    def test_wrong_source_and_wrong_recipe_fail(self):
        cell,unit=encoded_identity('OTHER','AAAC-1','L1',RECIPES[0]);r=recover_rows('SRC',[cell],[unit],['AAAC-1'],['L1'],RECIPES[0]);self.assertFalse(r['complete_one_to_one'])
        cell,unit=encoded_identity('SRC','AAAC-1','L1',RECIPES[1]);self.assertFalse(recover_rows('SRC',[cell],[unit],['AAAC-1'],['L1'],RECIPES[0])['complete_one_to_one'])
    def test_missing_and_reserved_keys_and_empty_archive(self):
        self.assertIsNone(encoded_identity('SRC','','L1',RECIPES[0]))
        with self.assertRaises(ValueError):encoded_identity('SRC','cell|x','L1',RECIPES[0])
        with self.assertRaises(ValueError):encoded_identity('SRC','cell','L1','suffix_only')
        self.assertFalse(recover_rows('SRC',[],[],['cell'],['L1'],RECIPES[0])['complete_one_to_one'])
    def test_axis_mismatch_rejected(self):
        with self.assertRaises(ValueError):recover_rows('SRC',['cell'],[],['cell'],['L1'],RECIPES[0])
        with self.assertRaises(ValueError):recover_rows('SRC',[],[],['cell'],[],RECIPES[0])
if __name__=='__main__':unittest.main()

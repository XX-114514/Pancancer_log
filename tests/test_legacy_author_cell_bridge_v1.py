import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from legacy_author_cell_bridge_v1 import encoded_author_cell,decode_legacy_global_cell,exact_author_bridge

def legacy(unit,cell,source='GSE_TEST'):
    return source+'|project_v3_sample_id='+source+'__'+unit+'|cell='+cell

class AuthorBridgeTests(unittest.TestCase):
    def test_empty_archive_does_not_pass_identity_acceptance(self):
        self.assertFalse(exact_author_bridge('GSE_TEST',[],['p1:AAAC'],['pool'])['complete_one_to_one'])
    def test_sample_prefix_prevents_shared_barcode_collision(self):
        q=[legacy('pool','GSE_TEST__pool__p2t_AAAC'),legacy('pool','GSE_TEST__pool__p1t_AAAC')]
        r=exact_author_bridge('GSE_TEST',q,['p1t:AAAC','p2t:AAAC'],['pool','pool'])
        self.assertEqual(r['author_row_position'],[1,0]);self.assertTrue(r['complete_one_to_one'])
    def test_library_namespace_is_not_erased(self):
        r=exact_author_bridge('GSE_TEST',[legacy('pool2','GSE_TEST__pool2__p1t_AAAC')],['p1t:AAAC'],['pool1'])
        self.assertEqual(r['author_row_position'],[-1])
    def test_colon_encoding_collision_is_ambiguous(self):
        r=exact_author_bridge('GSE_TEST',[legacy('pool','GSE_TEST__pool__p1_AAAC')],['p1:AAAC','p1_AAAC'],['pool','pool'])
        self.assertEqual(r['author_row_position'],[-2]);self.assertFalse(r['complete_one_to_one'])
    def test_duplicate_identical_author_rows_are_not_silently_collapsed(self):
        q=legacy('pool','GSE_TEST__pool__p1_AAAC')
        r=exact_author_bridge('GSE_TEST',[q],['p1:AAAC','p1:AAAC'],['pool','pool'])
        self.assertEqual(r['status_counts'],{'ambiguous_author_identity':1})
    def test_missing_and_wrong_namespaces_remain_explicit(self):
        r=exact_author_bridge('GSE_TEST',['other|cell=x',legacy('pool','GSE_TEST__pool__absent')],[''],['pool'])
        self.assertEqual(r['author_row_position'],[-3,-1]);self.assertFalse(r['complete_one_to_one'])
    def test_repeated_archive_identity_is_reported_without_deleting_cells(self):
        q=legacy('pool','GSE_TEST__pool__p1_AAAC')
        r=exact_author_bridge('GSE_TEST',[q,q],['p1:AAAC'],['pool'])
        self.assertEqual(r['author_row_position'],[0,0]);self.assertEqual(r['repeated_author_rows_in_archive'],{0:2});self.assertFalse(r['complete_one_to_one'])
    def test_barcode_suffix_and_case_are_preserved(self):
        q=legacy('pool','GSE_TEST__pool__p1_AAAC-2')
        r=exact_author_bridge('GSE_TEST',[q],['p1:AAAC-1','P1:AAAC-2','p1:AAAC-2'],['pool']*3)
        self.assertEqual(r['author_row_position'],[2])
    def test_declared_unit_and_cell_must_agree(self):
        self.assertIsNone(decode_legacy_global_cell('GSE_TEST',legacy('pool1','GSE_TEST__pool2__x')))
        self.assertIsNone(decode_legacy_global_cell('GSE_TEST',legacy('pool1','GSE_TEST__pool1__x')+'|cell=y'))
    def test_incomplete_author_axes_or_reserved_delimiter_rejected(self):
        with self.assertRaises(ValueError):exact_author_bridge('GSE_TEST',[],['x'],[])
        with self.assertRaises(ValueError):encoded_author_cell('GSE_TEST','pool','p|cell=x')

if __name__=='__main__':unittest.main()

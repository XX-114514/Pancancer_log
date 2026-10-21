"""Check chunk-boundary duplicate/decimal preservation and original gates."""
import unittest
from test_Chen_validation_RNA_RAM_v1 import fixture
from qualify_Chen_validation_RNA_cells_in_RAM_v2 import profile,packed_key


class ChunkedInputChecks(unittest.TestCase):
    def test_exact_string_keys_not_float_rounding(self):
        self.assertNotEqual(packed_key('test','1','2'),packed_key('test','1.0','2'))
        with self.assertRaises(AssertionError):packed_key('x\0y','1','2')

    def test_duplicate_across_chunk_boundary_and_negative_fractional(self):
        rows=[dict(slide='synthetic_image',x=i,y=i,type='synthetic_type',comp='test',FOXP3=-.5,CD4=.2) for i in range(10000)]
        rows.append(rows[0]);d,_,_=profile(fixture(rows))
        self.assertEqual(d['coordinate_duplicate_rows_within_image'],1)
        self.assertEqual(d['negative_gene_entries'],10001)
        self.assertEqual(d['fractional_gene_entries'],20002)

    def test_missing_gene_stays_unqualified(self):
        with self.assertRaises(AssertionError):profile(fixture([dict(slide='synthetic',x=1,y=2,type='synthetic_type')],omit='FOXP3'))


if __name__=='__main__':unittest.main()

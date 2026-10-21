"""Verify RNA payload qualification rejects missing axes and preserves semantics."""
import csv
import io
import unittest
from qualify_Chen_validation_RNA_cells_in_RAM_v1 import GENES,profile


def fixture(rows,omit=None):
    fields=[g for g in GENES+['slide','x','y','type','comp'] if g!=omit]
    b=io.StringIO();w=csv.DictWriter(b,fieldnames=fields);w.writeheader()
    for row in rows:w.writerow({k:row.get(k,0) for k in fields})
    return b.getvalue().encode()


class InputChecks(unittest.TestCase):
    def test_no_integer_molecule_or_patient_assumption(self):
        row=dict(slide='synthetic_image',x=1,y=2,type='synthetic_type',comp='test',FOXP3=-.5,CD4=.2)
        d,groups,keys=profile(fixture([row,row]))
        self.assertEqual(d['negative_gene_entries'],2)
        self.assertEqual(d['fractional_gene_entries'],4)
        self.assertEqual(d['coordinate_duplicate_rows_within_image'],1)
        self.assertEqual(d['image_groups_not_reviewed_patients'],1)

    def test_missing_gene_is_not_zero(self):
        row=dict(slide='synthetic_image',x=1,y=2,type='synthetic_type',comp='test')
        with self.assertRaises(AssertionError):profile(fixture([row],omit='FOXP3'))


if __name__=='__main__':unittest.main()

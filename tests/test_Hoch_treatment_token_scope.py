import importlib.util
import json
from pathlib import Path
import unittest

P=Path(__file__).resolve().parents[1]/'scripts/qualify_Hoch_treatment_token_coverage_v1.py'
s=importlib.util.spec_from_file_location('treatment_tokens',P);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


class TreatmentTokenScope(unittest.TestCase):
    def test_abbreviations_boundaries_and_unknown_text_suppression(self):
        rows=[{'Drug_after_surgery':'Nivo + Ipi; synthetic_private_note'},
              {'Drug_after_surgery':'notnivolumab; synthetic_private_note'},
              {'Drug_after_surgery':'Pembrolizumab; PD-1'}]
        result=m.summarize(rows)
        patterns=result['field_coverage'][0]['source_metadata_row_token_patterns']
        self.assertEqual(patterns['nivo|ipi'],1)
        self.assertEqual(patterns['pembrolizumab|pd-1'],1)
        self.assertEqual(patterns['NO_WHITELIST_TOKEN_OR_MISSING'],1)
        self.assertNotIn('synthetic_private_note',json.dumps(result))
        self.assertTrue(result['tokens_not_equated_with_regimen_or_schedule'])
    def test_metadata_repeats_not_called_independent_patients(self):
        result=m.summarize([{'Drug_after_surgery':'ipilimumab'}]*3)
        self.assertEqual(result['field_coverage'][0]['source_metadata_row_token_patterns']['ipilimumab'],3)
        self.assertTrue(result['counts_are_metadata_rows_not_independent_patients'])


if __name__=='__main__':unittest.main()

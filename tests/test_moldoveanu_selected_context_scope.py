"""Only selected material context must be unique/concordant, without first-match."""
import ast
import importlib.util
from pathlib import Path
import unittest
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];path=ROOT/'scripts/analyze_moldoveanu_fixed_Ki67_source_organization_v3.py'
s=importlib.util.spec_from_file_location('context_scope',path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


class ContextScope(unittest.TestCase):
    def test_unselected_duplicate_keys_do_not_reject_complete_selected_context(self):
        data=pd.DataFrame({'Sample_ID':['selected','other','other',None,None],'Treatment':['A','B','C','D','E']})
        selected,audit=m.qualify_context_metadata(data,{'selected'},list(data))
        self.assertTrue(audit['selected_context_accepted']);self.assertEqual(len(selected),1)
    def test_only_exactly_equal_used_context_rows_collapse(self):
        data=pd.DataFrame({'Sample_ID':['selected','selected'],'Treatment':['A','A']})
        selected,audit=m.qualify_context_metadata(data,{'selected'},list(data))
        self.assertTrue(audit['selected_context_accepted']);self.assertEqual(len(selected),1)
        self.assertEqual(audit['exactly_equal_selected_context_rows_collapsed'],1)
        self.assertTrue(audit['clinical_response_or_full_row_equivalence_not_asserted'])
    def test_conflict_or_missing_selected_context_stays_unaccepted(self):
        data=pd.DataFrame({'Sample_ID':['selected','selected'],'Treatment':['A','B']})
        _,audit=m.qualify_context_metadata(data,{'selected'},list(data))
        self.assertFalse(audit['selected_context_accepted']);self.assertEqual(audit['selected_materials_with_context_conflict'],1)
        _,audit=m.qualify_context_metadata(data,{'missing'},list(data))
        self.assertFalse(audit['selected_context_accepted']);self.assertEqual(audit['selected_materials_missing_context'],1)
    def test_original_mathematics_and_final_marker_recipes_unchanged(self):
        old=ast.parse((ROOT/'scripts/analyze_moldoveanu_fixed_Ki67_source_organization_v2.py').read_text());new=ast.parse(path.read_text())
        def functions(tree):return {x.name:ast.dump(x) for x in tree.body if isinstance(x,ast.FunctionDef) and x.name not in ['analyze','qualify_context_metadata']}
        self.assertEqual(functions(old),functions(new))


if __name__=='__main__':unittest.main()

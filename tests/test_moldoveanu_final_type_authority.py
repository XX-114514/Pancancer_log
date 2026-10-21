"""Published final recipes qualify types; legacy-column equality is not assumed."""
import ast
import importlib.util
from pathlib import Path
import unittest
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
path=ROOT/'scripts/analyze_moldoveanu_fixed_Ki67_source_organization_v2.py'
s=importlib.util.spec_from_file_location('final_types',path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


class FinalTypeAuthority(unittest.TestCase):
    def frame(self):
        return pd.DataFrame({'Cluster':['Treg','Th.ae','Th.naive','Tc.ae','Tc.naive','others'],
            'Cluster.v2':['Treg','Th.ae','Th.naive','Tc.ae','Tc.naive','unclassified'],
            'prob.FoxP3':[.5,.49,0.,.1,.1,.1],
            'prob.CD45RO':[.8,.5,.49,.5,0.,.1]})
    def test_declared_final_recipes_and_legacy_difference_are_separate(self):
        result=m.audit_final_author_phenotypes(self.frame())
        self.assertTrue(result['author_probability_recipes_all_passed'])
        self.assertEqual(result['legacy_label_disagreements'],1)
        self.assertTrue(all(x['violations']==0 for x in result['source_recipe_checks']))
    def test_invalid_or_missing_probability_fails_recipes(self):
        for value in [.49,np.nan,1.1]:
            frame=self.frame();frame.loc[0,'prob.FoxP3']=value
            result=m.audit_final_author_phenotypes(frame)
            self.assertFalse(result['author_probability_recipes_all_passed'])
    def test_all_prepared_statistics_and_thresholds_unchanged(self):
        old=ast.parse((ROOT/'scripts/analyze_moldoveanu_fixed_Ki67_source_organization_v1.py').read_text())
        new=ast.parse(path.read_text())
        def functions(tree):return {x.name:ast.dump(x) for x in tree.body if isinstance(x,ast.FunctionDef) and x.name not in ['analyze','audit_final_author_phenotypes']}
        self.assertEqual(functions(old),functions(new))


if __name__=='__main__':unittest.main()

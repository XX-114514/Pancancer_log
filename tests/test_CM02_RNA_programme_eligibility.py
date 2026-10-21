"""Exercise actual pure qualification functions without server-only packages."""
import ast
from collections import Counter,defaultdict
import math
from pathlib import Path
import unittest
import numpy as np
import pandas as pd

PATH=Path(__file__).resolve().parents[1]/'scripts/qualify_CM02_RNA_programme_eligibility_v1.py'
tree=ast.parse(PATH.read_text())
names={'SOURCE','CONTEXTS','PARENTS','REGULATORY','CONVENTIONAL','CD8'}
selected=[]
for node in tree.body:
    if isinstance(node,ast.FunctionDef) and node.name in {'feature_coverage','reference_coverage'}:selected.append(node)
    if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id in names for t in node.targets):selected.append(node)
scope=dict(np=np,pd=pd,Counter=Counter,defaultdict=defaultdict,math=math)
exec(compile(ast.Module(body=selected,type_ignores=[]),str(PATH),'exec'),scope)


class RNAEligibility(unittest.TestCase):
    def test_missing_ambiguous_and_same_feature_never_become_independent_genes(self):
        var=pd.DataFrame({'gene_symbol':['A','B','C'],'source_gene_aliases':['D','ambig','ambig']},index=['E1','E2','E3'])
        result=scope['feature_coverage'](var,np.ones(3,dtype=bool),{'synthetic':['A','D','ambig','missing']})[0]
        self.assertEqual(result['uniquely_mapped_measured_genes'],0)
        self.assertEqual(result['same_feature_multiple_gene_names'],['A','D'])
        self.assertEqual(result['ambiguous_gene_names'],['ambig'])
        self.assertEqual(result['missing_gene_names'],['missing'])
        self.assertFalse(result['source_coverage_accepted'])
    def fixture(self):
        libraries=[];cells=[]
        for patient in range(21):
            lib=f'synthetic_lib_{patient}';libraries.append({'library_id':lib,'patient_registry':f'synthetic_case_{patient}','timepoint_registry':'pre','neo_type':'D','HPV/p16':'negative'})
            for label,node in [('Prolif Treg','Rnetwork'),('CXCL13+ CD4T','outside'),('Tprolif','Enetwork')]:
                for i in range(20):cells.append(dict(position=len(cells),library_id=lib,timepoint='pre',label=label,node_id=node,major='T'))
        return pd.DataFrame(libraries),pd.DataFrame(cells)
    def test_current_T_filter_and_missing_network_parent_denominator(self):
        tab,frame=self.fixture();frame.loc[0,'major']='B'
        result=scope['reference_coverage'](tab,frame,{'Rnetwork','Enetwork'})
        pairs=[x for x in result['pair_coverage'] if x['phase']=='pre']
        self.assertEqual(pairs[0]['both_parents_atleast20_reviewed_cases'],20)
        self.assertEqual(pairs[1]['both_parents_atleast20_reviewed_cases'],21)
        self.assertEqual(pairs[2]['both_parents_atleast20_reviewed_cases'],20)
        self.assertEqual(pairs[3]['both_parents_atleast20_reviewed_cases'],0)
        self.assertFalse(pairs[3]['future_minimum20_case_coverage_only'])
        self.assertEqual(result['author_reference_current_nonT_exclusions'][0]['source_cells'],1)
        self.assertNotIn('synthetic_case_',str(result))
    def test_phase_disagreement_and_duplicated_cells_rejected(self):
        tab,frame=self.fixture();frame.loc[0,'timepoint']='post'
        with self.assertRaises(AssertionError):scope['reference_coverage'](tab,frame,set())
        frame.loc[0,'timepoint']='pre';frame.loc[1,'position']=frame.loc[0,'position']
        with self.assertRaises(AssertionError):scope['reference_coverage'](tab,frame,set())
    def test_unknown_author_state_not_automatically_assigned(self):
        tab,frame=self.fixture();frame.loc[0,'label']='unknown_state'
        with self.assertRaises(AssertionError):scope['reference_coverage'](tab,frame,set())


if __name__=='__main__':unittest.main()

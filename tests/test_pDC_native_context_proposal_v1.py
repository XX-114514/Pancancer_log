import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from pDC_native_context_proposal_v1 import context_proposal,CORE_PANELS
class SpecificContextCounterexamples(unittest.TestCase):
 def run_case(self,positive=(),unmeasured=()):
  genes=sorted({g for gs,k in CORE_PANELS.values() for g in gs}|{'LILRA4','PLD4','TSPAN13','C1QA','C1QB','LYZ','LST1','SDC1','TNFRSF17','MZB1','JCHAIN','DERL3','XBP1','CD79A','CD37','CD74','GZMB','NKG7','FCER1G'})
  x=np.zeros((1,len(genes)));m=np.ones(len(genes),bool)
  for g in positive:x[0,genes.index(g)]=3
  for g in unmeasured:m[genes.index(g)]=False;x[0,genes.index(g)]=np.nan
  return context_proposal(x,genes,m)['disposition_code'][0]
 def test_shared_APC_signals_do_not_create_B_competition(self):self.assertEqual(self.run_case(['LILRA4','PLD4','CD79A','CD37','CD74']),0)
 def test_coherent_B_program_retains_review(self):self.assertEqual(self.run_case(['LILRA4','PLD4','MS4A1','CD79B']),1)
 def test_specific_T_signal_is_not_overridden(self):self.assertEqual(self.run_case(['LILRA4','PLD4','CD3D','TRAC']),1)
 def test_missing_pDC_anchor_is_unknown(self):self.assertEqual(self.run_case(['PLD4'],['LILRA4']),2)
 def test_two_missing_B_genes_can_change_context(self):self.assertEqual(self.run_case(['LILRA4','PLD4'],['MS4A1','CD19']),2)
 def test_one_missing_B_gene_cannot_meet_two_gene_pattern(self):self.assertEqual(self.run_case(['LILRA4','PLD4'],['MS4A1']),0)
 def test_shared_secretory_or_cytotoxic_program_is_not_other_lineage(self):self.assertEqual(self.run_case(['LILRA4','PLD4','JCHAIN','XBP1','GZMB','NKG7','FCER1G']),0)
 def test_absent_additional_pDC_support_is_not_rescue(self):self.assertEqual(self.run_case(['CD79A','CD74']),3)
if __name__=='__main__':unittest.main()

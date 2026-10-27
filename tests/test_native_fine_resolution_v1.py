import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from native_fine_resolution_v1 import evaluate_patterns,FINE_ANCHOR_PATTERNS

class MeasurementCounterexamples(unittest.TestCase):
 def run_rule(self,genes,values,measured,rule):
  return evaluate_patterns(np.packbits(np.asarray(values,bool),axis=1,bitorder='little'),genes,np.asarray(measured,bool),rule)
 def test_unmeasured_required_gene_is_not_zero(self):
  x=self.run_rule(['A','B'],[[1,0]],[True,False],{'t':[('A',),('B',)]})
  self.assertEqual(x['positive_anchor_bitmask'].tolist(),[0]);self.assertEqual(x['unassessable_anchor_bitmask'].tolist(),[1])
 def test_measured_alternative_can_establish_positive(self):
  x=self.run_rule(['A','B'],[[1,0]],[True,False],{'t':[('A','B')]})
  self.assertEqual(x['positive_anchor_bitmask'].tolist(),[1]);self.assertFalse(x['coverage']['t']['all_requested_genes_measured'])
 def test_unmeasured_bits_cannot_create_support(self):
  x=self.run_rule(['A','B'],[[1,1]],[False,False],{'t':[('A',),('B',)]})
  self.assertEqual(x['positive_anchor_bitmask'].tolist(),[0]);self.assertEqual(x['unassessable_anchor_bitmask'].tolist(),[1])
 def test_known_missing_required_clause_blocks_pattern(self):
  x=self.run_rule(['A','B'],[[0,0]],[True,False],{'t':[('A',),('B',)]})
  self.assertEqual(x['positive_anchor_bitmask'].tolist(),[0]);self.assertEqual(x['unassessable_anchor_bitmask'].tolist(),[0])
 def test_S100_pair_does_not_establish_neutrophil_anchor(self):
  genes=['CSF3R','FCGR3B','S100A8','S100A9'];rule={'n':FINE_ANCHOR_PATTERNS['Neutrophil']}
  x=self.run_rule(genes,[[0,0,1,1],[1,1,1,0]],[True]*4,rule)
  self.assertEqual(x['positive_anchor_bitmask'].tolist(),[0,1])
 def test_Treg_and_CD4_patterns_can_coexist(self):
  genes=['FOXP3','IL2RA','CD3D','CD3E','CD4','CTLA4'];rule={t:FINE_ANCHOR_PATTERNS[t] for t in ['Treg','CD4_T']}
  x=self.run_rule(genes,[[1,1,1,0,1,0]],[True]*6,rule)
  self.assertEqual(x['positive_anchor_bitmask'].tolist(),[3])
 def test_missing_axis_and_duplicate_axis(self):
  x=self.run_rule(['A'],[[1]],[True],{'t':[('A',),('B',)]});self.assertEqual(x['unassessable_anchor_bitmask'].tolist(),[1])
  with self.assertRaises(ValueError):self.run_rule(['A','A'],[[1,1]],[True,True],{'t':[('A',)]})
if __name__=='__main__':unittest.main()

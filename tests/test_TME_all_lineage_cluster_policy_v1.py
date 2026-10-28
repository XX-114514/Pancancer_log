import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from TME_all_lineage_cluster_policy_v1 import GATE_ORDER,coherent_flags,annotate_cluster

class AllLineagePolicyTests(unittest.TestCase):
 def call(self,genes,x=None):
  if x is None:x=np.ones((100,len(genes)))
  f,a=coherent_flags(x,genes);return annotate_cluster(len(x),dict(zip(GATE_ORDER,f.sum(axis=0))),a)
 def test_astrocytic_RNA_not_forced_immune_or_malignant(self):
  r=self.call(['GFAP','AQP4','SLC1A3']);self.assertEqual(r['celltype'],'Astrocyte');self.assertNotIn('malignancy',r)
 def test_neuronal_RNA_requires_coexpression(self):
  r=self.call(['SNAP25','RBFOX3'],np.tile(np.eye(2),(50,1)));self.assertEqual(r['major'],'unresolved')
 def test_glial_RNA(self):self.assertEqual(self.call(['SOX10','PLP1'])['major'],'Glial')
 def test_erythroid_not_based_on_hemoglobin_alone(self):
  self.assertEqual(self.call(['HBB'])['major'],'unresolved');self.assertEqual(self.call(['HBB','AHSP'])['major'],'Erythroid')
 def test_platelet_complete_module(self):self.assertEqual(self.call(['PF4','GP9'])['major'],'Platelet')
 def test_T_and_astrocyte_competition_retained_unresolved(self):
  r=self.call(['CD3D','TRAC','LCK','GFAP','SLC1A3']);self.assertEqual(r['major'],'unresolved');self.assertIn('competing',r['status'])
 def test_cDC1_parent_is_assessable_by_DC_module(self):
  r=self.call(['CLEC9A','XCR1']);self.assertEqual(r['major'],'Myeloid');self.assertEqual(r['celltype'],'cDC1')
 def test_lymphatic_parent_is_assessable(self):
  r=self.call(['PROX1','PDPN']);self.assertEqual(r['major'],'Endothelial');self.assertEqual(r['celltype'],'Lymphatic_endothelial_related')
if __name__=='__main__':unittest.main()

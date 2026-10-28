"""Counterexamples that would corrupt cluster naming or missingness."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from TME_cluster_annotation_policy_v1 import GATE_ORDER,annotate_cluster,coherent_flags,state_context

class PolicyTests(unittest.TestCase):
 def call(self,counts,missing=(),n=100):
  assessed={k:k not in missing for k in GATE_ORDER};return annotate_cluster(n,counts,assessed)
 def test_CD4_parent_of_Treg_not_a_second_fine_call(self):
  g=['CD3D','TRAC','LCK','FOXP3','IL2RA','CD4','IL7R'];x=np.ones((10,len(g)));f,a=coherent_flags(x,g)
  self.assertTrue(f[:,GATE_ORDER.index('Treg')].all());self.assertFalse(f[:,GATE_ORDER.index('CD4_T')].any())
 def test_marker_from_separate_cells_does_not_form_a_coherent_module(self):
  g=['CD3D','TRAC','LCK'];x=np.eye(3);f,a=coherent_flags(x,g);self.assertTrue(a['T']);self.assertEqual(f[:,0].sum(),0)
 def test_ambiguous_fine_keeps_T(self):
  x=self.call({'T':100,'CD4_T':50,'CD8_T':50});self.assertEqual(x['celltype'],'Tcell');self.assertEqual(len(x['fine_candidates']),2)
 def test_T_plus_B_not_whole_cluster_T(self):
  self.assertEqual(self.call({'T':60,'B':30})['major'],'unresolved')
 def test_unmeasured_T_does_not_prove_NK(self):
  self.assertEqual(self.call({'NK':90},missing=['T'])['major'],'unresolved')
 def test_gamma_is_RNA_candidate(self):
  self.assertEqual(self.call({'T':80,'gamma_delta_T_RNA_candidate':60})['celltype'],'gamma_delta_T_RNA_candidate')
 def test_small_cluster_descriptive(self):
  self.assertEqual(self.call({'T':10},n=10)['status'],'small_cluster_descriptive_only')
 def test_state_does_not_change_type_and_missing_is_explicit(self):
  s=state_context(['IFIT1','ISG15'],np.array([2.,2.]),np.array([0.,0.]),np.array([.5,.5]));self.assertTrue(s['interferon_related']['descriptively_enriched']);self.assertIsNone(s['cycling_related']['mean_native_logCP10K']);self.assertFalse(s['cycling_related']['descriptively_enriched'])

if __name__=='__main__':unittest.main()

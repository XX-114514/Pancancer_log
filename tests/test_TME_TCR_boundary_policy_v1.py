"""Biologically meaningful counterexamples for the additional T RNA tier."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
try:
    import numpy as np
    from TME_TCR_boundary_policy_v1 import CORE_GENES,FAMILIES,evaluate
    AVAILABLE=True
except ImportError:
    AVAILABLE=False

@unittest.skipUnless(AVAILABLE,'Use the existing scientific environment; production rejects skips')
class BoundaryTests(unittest.TestCase):
    def run_case(self,on,**override):
        genes=sorted({g for p in CORE_GENES.values() for g in p});d=np.asarray([[g in on for g in genes]],bool);mask=np.ones(len(genes),bool)
        missing=override.pop('missing',[])
        for g in missing:mask[genes.index(g)]=False;d[0,genes.index(g)]=False
        args={'old_native_mask':np.asarray([5],np.uint16),'reference_mask':np.asarray([1],np.uint16),'low_reference_T':[True],
              'local_T_context':[True],'fit':[True],'positive':[True],'baseline_unresolved':[True]};args.update(override)
        return evaluate(d,genes,mask,**args)
    base=['CD3D','CD3E','TRAC','LCK']
    def test_shared_presenting_pair_does_not_override_TCR(self):
        self.assertTrue(self.run_case(self.base+['LST1','AIF1'])['proposed_T'][0])
    def test_cytotoxicity_without_TCR_chain_not_T(self):
        self.assertFalse(self.run_case(['CD3D','CD3E','LCK'])['proposed_T'][0])
    def test_complete_myeloid_core_keeps_conflict(self):
        self.assertEqual(self.run_case(self.base+['CSF1R','LST1','LYZ'])['review_code'][0],3)
    def test_B_receptor_core_keeps_conflict(self):
        self.assertEqual(self.run_case(self.base+['CD79A','CD79B'])['review_code'][0],3)
    def test_keratin_core_keeps_conflict(self):
        self.assertEqual(self.run_case(self.base+['KRT8','KRT18','KRT19'])['review_code'][0],3)
    def test_complete_erythroid_frozen_support_is_not_dropped(self):
        self.assertEqual(self.run_case(self.base,old_native_mask=np.asarray([1|(1<<11)],np.uint16))['review_code'][0],3)
    def test_missing_competitor_anchor_not_negative(self):
        self.assertEqual(self.run_case(self.base,missing=['CSF1R'])['review_code'][0],4)
    def test_missing_TCR_chain_not_positive(self):
        self.assertEqual(self.run_case(self.base,missing=['TRAC'])['review_code'][0],2)
    def test_reference_cross_lineage_kept(self):
        self.assertEqual(self.run_case(self.base,reference_mask=np.asarray([5],np.uint16))['review_code'][0],5)
    def test_NK_reference_not_overridden_by_cytoplasmic_CD3(self):
        self.assertEqual(self.run_case(self.base,low_reference_T=[False])['review_code'][0],6)
    def test_local_group_disagrees(self):
        self.assertEqual(self.run_case(self.base,local_T_context=[False])['review_code'][0],7)
    def test_unfit_and_original_labels_not_promoted(self):
        self.assertEqual(self.run_case(self.base,fit=[False])['review_code'][0],8)
        self.assertEqual(self.run_case(self.base,baseline_unresolved=[False])['review_code'][0],0)
    def test_gamma_chain_support_does_not_assign_fine_gamma_identity(self):
        r=self.run_case(['CD3D','CD3E','TRDC','TRGC1','LAT']);self.assertTrue(r['proposed_T'][0]);self.assertTrue(r['gamma_chain_detected'][0]);self.assertNotIn('celltype',r)
    def test_zero_native_RNA_not_promoted(self):
        self.assertEqual(self.run_case(self.base,positive=[False])['review_code'][0],9)

if __name__=='__main__':unittest.main()

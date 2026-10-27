import json, sys, tempfile, unittest, hashlib
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from review_source_qualified_identity_context_v1 import extra_context, source_status

class ContextReviewTests(unittest.TestCase):
    def test_missing_measurement_is_different_from_measured_zero(self):
        coverage,v=extra_context(np.array([2]),np.zeros((1,1)),np.zeros((1,1)),['LCK'],[])
        self.assertEqual(v['LCK']['detection_fraction'][0],0)
        self.assertIsNone(v['TRBC1'])
        self.assertIn('TRBC1',coverage['T_receptor_context']['not_measured'])
    def test_original_manual_genes_removed_not_renormalized_as_missing_zero(self):
        c,v=extra_context(np.array([2,4]),np.array([[4.,0],[0,8]]),np.array([[2,0],[0,2]]),['LCK','LAT'],['LAT'])
        self.assertNotIn('LAT',v);self.assertEqual(v['LCK']['mean_CP10K'].tolist(),[2.,0.])
        self.assertEqual(c['T_receptor_context']['measured'],['LCK'])
    def test_invalid_denominator_and_duplicate_gene_identity_fail(self):
        for n,cp,d,genes in [(np.array([0]),np.zeros((1,1)),np.zeros((1,1)),['LCK']),
                             (np.array([2]),np.zeros((1,2)),np.zeros((1,2)),['LCK','LCK']),
                             (np.array([2]),np.zeros((1,1)),np.array([[3]]),['LCK'])]:
            with self.assertRaises(ValueError):extra_context(n,cp,d,genes,[])
    def fixture(self,root,rows=(9,4),status=(1,0)):
        p=root/'original.json';p.write_text(json.dumps({'source':'SRC'}))
        q=root/'override.npz';np.savez(q,native_row=rows,source_qualification_code=status,qualification_work_unit_code=[99,98])
        sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
        v=root/'v4.json';v.write_text(json.dumps({'bindings':{'author_bridge_registry':{'path':str(p),'sha256':sha(p)}},'override_arrays_path':str(q),'override_arrays_sha256':sha(q)}))
        return {'base_v4_manifest':{'path':str(v)},'new_source_overrides':{}}
    def test_mixed_parent_keeps_cell_override_and_ignores_old_dictionary_positions(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest=self.fixture(Path(tmp));s={'source':'SRC','rows':np.array([9,4]),'ua':{}}
            self.assertEqual(source_status(s,[],manifest).tolist(),[1,0])
    def test_row_axis_swap_is_not_accepted_as_same_sample(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest=self.fixture(Path(tmp));s={'source':'SRC','rows':np.array([4,9]),'ua':{}}
            with self.assertRaises(AssertionError):source_status(s,[],manifest)
    def test_current_override_takes_precedence_without_parent_admission(self):
        with tempfile.TemporaryDirectory() as tmp:
            r=Path(tmp);manifest=self.fixture(r);p=r/'latest.npz';np.savez(p,native_row=[9,4],source_qualification_code=[0,2]);manifest['new_source_overrides']['SRC']={'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
            self.assertEqual(source_status({'source':'SRC','rows':np.array([9,4]),'ua':{}},[],manifest).tolist(),[0,2])

if __name__=='__main__':unittest.main()

import ast
import importlib.util
from pathlib import Path
import unittest
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def module(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


class SourcePipeline(unittest.TestCase):
    def test_actual_CSR_block_selection_scores_against_independent_dense_rank(self):
        core=ast.parse((ROOT/'scripts/analyze_CM02_RNA_multigene_in_RAM_v1.py').read_text())
        function=next(s for s in core.body if isinstance(s,ast.FunctionDef) and s.name=='read_scores')
        from types import SimpleNamespace
        # Production is Linux (KiB); macOS ru_maxrss uses bytes. Keep this
        # synthetic reader test independent of that platform-specific unit.
        resource=SimpleNamespace(RUSAGE_SELF=0,getrusage=lambda unused:SimpleNamespace(ru_maxrss=1024))
        namespace=dict(np=np,pd=pd,scoring=module('ucell_native_sparse_scores_v1'),resource=resource)
        exec(compile(ast.Module(body=[function],type_ignores=[]),'selected_read_scores','exec'),namespace)
        class Raw(dict):attrs={'shape':(5002,5000)}
        ptr=np.zeros(5003,dtype=np.int64);ptr[1:5001]=4;ptr[5001:]=7
        raw=Raw(indptr=ptr,indices=np.array([4,4,2,4999,0,2,14]),data=np.array([2,3,2,500,4,6,1]))
        human=np.ones(5000,bool);human[-1]=False;measured=np.ones(5000,bool);progress=[]
        scores,total,detected,receipt=namespace['read_scores']({'layers/counts':raw},np.array([0,5000]),human,measured,[np.array([0,2,4]),np.array([1,4,14])],progress.append)
        np.testing.assert_array_equal(total,[7,11]);np.testing.assert_array_equal(detected,[2,3])
        self.assertEqual(scores.shape,(2,2));self.assertEqual(receipt['independent_dense_boundary_rows'],2)
        self.assertEqual(progress[-1]['completed_source_cells'],2)

    def test_independent_random_replay_positive_negative_and_nonestimable(self):
        stats=module('RNA_programme_partial_rank_v1');audit=module('audit_RNA_programme_partial_rank_v1')
        rng=np.random.default_rng(22);n=44;c=rng.normal(size=(n,8));cat=np.array([['D' if i%2 else 'D+T','positive' if i%3 else 'negative'] for i in range(n)])
        x=rng.normal(size=n)+c[:,0];y=.7*x+rng.normal(size=n)
        for sign in [-1,1]:
            result=stats.analyze(x,sign*y,c,cat,5,draws=150,required_bootstrap=100)
            checked=audit.replay(x,sign*y,c,cat,5,result,draws=150,required_bootstrap=100)
            self.assertEqual(checked['status'],'ALL_BOOTSTRAP_AND_PERMUTATION_REPLAYED');self.assertEqual(checked['observed_sign'],sign)
        c[:,2]=c[:,1];result=stats.analyze(x,y,c,cat,5,draws=150,required_bootstrap=100)
        self.assertEqual(audit.replay(x,y,c,cat,5,result)['status'],'NONESTIMABLE_GATE_REPLAYED')


if __name__=='__main__':unittest.main()

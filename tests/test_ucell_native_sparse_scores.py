import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import unittest
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('sparse_ucell',ROOT/'scripts/ucell_native_sparse_scores_v1.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def expected_dense(values,signatures,maximum):
    ranks=pd.Series(-np.asarray(values)).rank(method='average').to_numpy()
    maximum=min(maximum,len(values));result=[]
    for indices in signatures:
        selected=np.minimum(ranks[indices],maximum);n=len(indices);minimum=n*(n+1)/2
        result.append(1-(selected.sum()-minimum)/(n*maximum-minimum))
    return np.array(result)


class SparseUCell(unittest.TestCase):
    def test_noncanonical_CSR_duplicate_counts_order_and_native_universe(self):
        n=5000;human=np.ones(n,dtype=bool);measured=human.copy();human[-1]=False
        ptr=np.array([0,4,5,7]);indices=np.array([3,3,1,n-1,7,1,3]);counts=np.array([2,3,2,100,1,4,1])
        programmes=[np.array([0,1,3]),np.array([7,100,200])]
        scores,total,detected=m.score_selected_CSR_rows(ptr,indices,counts,[0,2],human,measured,programmes)
        dense=np.zeros((2,n-1));dense[0,3]=5;dense[0,1]=2;dense[1,1]=4;dense[1,3]=1
        np.testing.assert_array_equal(total,[7,5]);np.testing.assert_array_equal(detected,[2,2])
        for i,row in enumerate(dense):np.testing.assert_allclose(scores[i],expected_dense(row,programmes,1500),atol=1e-12)
    def test_unmeasured_positive_CSR_and_unmeasured_signature_rejected(self):
        human=np.ones(5000,dtype=bool);measured=human.copy();human[3]=False;measured[3]=False
        with self.assertRaises(AssertionError):m.score_selected_CSR_rows([0,1],[3],[1],[0],human,measured,[[1,2]])
        with self.assertRaises(AssertionError):m.score_selected_CSR_rows([0,1],[1],[1],[0],human,measured,[[1,3]])
        with self.assertRaises(AssertionError):m.score_selected_CSR_rows([0,1],[1],[1],[0,0],human,measured,[[1,2]])
        with self.assertRaises(AssertionError):m.score_selected_CSR_rows([0,1],[1],[1],[.5],human,measured,[[1,2]])
    def test_full_universe_zero_ties_and_cap_boundary_dense_replay(self):
        rng=np.random.default_rng(41)
        for count in [0,18,1200,2000]:
            x=np.zeros(5000);x[:count]=rng.integers(1,9,size=count)
            signatures=[np.array([0,1,2,3,4,5]),np.array([17,88,520,1100,2800,4500,4999])]
            for maximum in [50,1500]:
                got=m.score_positive_sparse_cell(x[x>0],[x[s] for s in signatures],len(x),maximum)
                np.testing.assert_allclose(got,expected_dense(x,signatures,maximum),atol=1e-12)
    def test_invalid_counts_universe_and_signature_are_rejected(self):
        cases=[([-1],[[0]],20,10),([1.5],[[1.5]],20,10),([1],[[2]],20,10),([1],[[1,0]],1,1)]
        for args in cases:
            with self.assertRaises(AssertionError):m.score_positive_sparse_cell(*args)
    def test_signature_order_is_irrelevant_and_all_zero_native_matrix(self):
        np.testing.assert_array_equal(m.score_positive_sparse_cell([],[[0,0,0],[0,0]],5000),[0.,0.])
        a=m.score_positive_sparse_cell([1,2,4],[[1,4,0]],5000)
        b=m.score_positive_sparse_cell([4,1,2],[[0,1,4]],5000)
        np.testing.assert_allclose(a,b,atol=1e-12)
    def test_official_hashed_u_stat_function_and_base_R_ranking(self):
        helper=ROOT.parent/'fullgene_resume_20261014/private_CM02_RNA_multigene_preflight_20261022_v1/public_sources/UCell_HelperFunctions.R'
        if not helper.exists():self.skipTest('Pinned public helper not provided in this checkout')
        self.assertEqual(hashlib.sha256(helper.read_bytes()).hexdigest(),'315e21a0428bae74b0895be4fe357cdc431a26ab7ffa27a864f029f373c7cc18')
        rng=np.random.default_rng(17);matrix=rng.integers(0,6,size=(80,6));matrix[25:]=0
        signatures={'synthetic_S':[0,3,9,44],'synthetic_G2M':[2,5,11,17,21,61]}
        script='''a<-commandArgs(trailingOnly=TRUE);e<-new.env();p<-parse(file=a[[1]])
hit<-vapply(p,function(x)is.call(x)&&as.character(x[[1]])=="<-"&&is.symbol(x[[2]])&&as.character(x[[2]])=="u_stat",logical(1))
stopifnot(sum(hit)==1L);eval(p[[which(hit)]],envir=e)
j<-jsonlite::fromJSON(file("stdin"));x<-as.matrix(j$matrix)
r<-data.frame(rn=as.character(seq_len(nrow(x))),apply(x,2,function(v)rank(-v,ties.method="average")),check.names=FALSE)
out<-vapply(j$signatures,function(s)e$u_stat(r,as.integer(s)+1,maxRank=min(j$maximum,nrow(x)),sparse=FALSE),numeric(ncol(x)))
cat(jsonlite::toJSON(unname(out),digits=16))'''
        payload=dict(matrix=matrix.tolist(),signatures=signatures,maximum=50)
        r=subprocess.run(['/usr/local/bin/Rscript','--vanilla','-e',script,str(helper)],input=json.dumps(payload),capture_output=True,text=True,timeout=30)
        self.assertEqual(r.returncode,0,r.stderr)
        reference=np.array(json.loads(r.stdout))
        got=np.array([m.score_positive_sparse_cell(v[v>0],[v[s] for s in signatures.values()],len(v),50) for v in matrix.T])
        np.testing.assert_allclose(got,reference,atol=1e-12)


if __name__=='__main__':unittest.main()

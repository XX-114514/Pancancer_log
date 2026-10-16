"""Exact integer-span certificates for nearly zero OLS rank residuals.

Numerical cancellation is not a nonconstant phenotype. This audit helper
does not alter production, permutation streams or scientific tolerances.
"""
from fractions import Fraction

import numpy as np


def modular_rank(matrix, prime=1000003):
    a=np.asarray(matrix,dtype=np.int64).copy()%prime
    rank=0
    for col in range(a.shape[1]):
        candidates=np.flatnonzero(a[rank:,col])
        if not len(candidates):continue
        pivot=rank+int(candidates[0]);a[[rank,pivot]]=a[[pivot,rank]]
        a[rank,col:]=(a[rank,col:]*pow(int(a[rank,col]),prime-2,prime))%prime
        for row in range(rank+1,len(a)):
            factor=int(a[row,col])
            if factor:a[row,col:]=(a[row,col:]-factor*a[rank,col:])%prime
        rank+=1
        if rank==len(a):break
    return rank


def rational_rank(matrix):
    matrix=np.asarray(matrix,dtype=np.int64)
    a=[[Fraction(int(v)) for v in row] for row in matrix]
    rank=0
    for col in range(matrix.shape[1]):
        pivot=next((row for row in range(rank,len(a)) if a[row][col]),None)
        if pivot is None:continue
        a[rank],a[pivot]=a[pivot],a[rank]
        div=a[rank][col];a[rank]=[v/div for v in a[rank]]
        for row in range(rank+1,len(a)):
            factor=a[row][col]
            if factor:a[row]=[u-factor*v for u,v in zip(a[row],a[rank])]
        rank+=1
        if rank==len(a):break
    return rank


def any_exactly_fitted(design_integer,response_integer):
    """Prove if either ranked response lies in the nuisance column space.

    Full row rank modulo a prime is a sufficient exact integer rank
    certificate. Other cases use rational elimination, never a modular
    dependence assumption. Joint duplicate rows may be collapsed exactly.
    """
    z=np.asarray(design_integer,dtype=np.int64)
    y=np.asarray(response_integer,dtype=np.int64)
    assert y.ndim==2 and len(z)==len(y)
    joint=np.unique(np.column_stack([z,y]),axis=0)
    zz=joint[:,:z.shape[1]];yy=joint[:,z.shape[1]:]
    if modular_rank(zz)==len(joint):
        # Rank is at most joint row count and at least the modular rank.
        # Both responses are constant on these exact duplicate groups.
        return True
    rank=rational_rank(zz)
    return any(rational_rank(np.column_stack([zz,yy[:,j]]))==rank for j in range(yy.shape[1]))

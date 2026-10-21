"""Persistent R Matrix reference for every streamed coordinate block.

Binary pipes only; no raw/cell files. Reference preserves full physical
columns and independently sums marker, namespace, detection and mito values.
"""
import os
import select
import struct
import subprocess

import numpy as np

R_CODE=r'''
suppressPackageStartupMessages(library(Matrix))
nc<-as.integer(Sys.getenv('PC_NC'));nf<-as.integer(Sys.getenv('PC_NF'))
kind<-as.integer(strsplit(Sys.getenv('PC_KIND'),'')[[1]])
mito<-as.integer(strsplit(Sys.getenv('PC_MITO'),'')[[1]])==1L
target<-as.integer(strsplit(Sys.getenv('PC_TARGET'),',',fixed=TRUE)[[1]])+1L
stopifnot(length(kind)==nf,length(mito)==nf)
input<-file('stdin','rb');output<-file('stdout','wb')
writeBin(as.integer(314159),output,size=4,endian='little');flush(output)
repeat{
 n<-readBin(input,integer(),n=1,size=4,endian='little')
 if(!length(n)||n==0L)break
 stopifnot(n>0L)
 v<-readBin(input,numeric(),n=n*3L,size=8,endian='little')
 stopifnot(length(v)==n*3L)
 a<-matrix(v,ncol=3L,byrow=TRUE)
 M<-sparseMatrix(i=as.integer(a[,2]),j=as.integer(a[,1]),x=a[,3],dims=c(nc,nf))
 H<-M[,which(kind==0L),drop=FALSE]
 marker<-as.matrix(M[,target,drop=FALSE])
 totals<-cbind(rowSums(H),rowSums(M[,which(kind==1L),drop=FALSE]),rowSums(M[,which(kind==2L),drop=FALSE]))
 detected<-rowSums(H>0)
 mitochondrial<-rowSums(M[,which(mito),drop=FALSE])
 result<-c(as.numeric(t(marker)),as.numeric(t(totals)),as.numeric(detected),as.numeric(mitochondrial))
 stopifnot(all(is.finite(result)))
 writeBin(as.double(result),output,size=8,endian='little');flush(output)
}
close(input);close(output)
'''


class RReference:
    def __init__(self,dims,kind,mito,target):
        self.nc=dims[1];self.k=len(target)
        env=dict(os.environ,PC_NC=str(self.nc),PC_NF=str(dims[0]),
            PC_KIND=''.join(map(str,kind.tolist())),PC_MITO=''.join('1' if x else '0' for x in mito),
            PC_TARGET=','.join(map(str,target.tolist())))
        self.p=subprocess.Popen([os.environ['PANCANCER_R_BLOCK_REFERENCE'],'--vanilla','-e',R_CODE],
            env=env,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    def read_exact(self,n):
        result=bytearray()
        while len(result)<n:
            ready,_,_=select.select([self.p.stdout],[],[],60)
            if not ready:raise RuntimeError('R_reference_response_timeout')
            data=os.read(self.p.stdout.fileno(),n-len(result))
            if not data:raise RuntimeError('R_reference_incomplete_response')
            result.extend(data)
        return bytes(result)

    def __enter__(self):
        try:assert struct.unpack('<i',self.read_exact(4))[0]==314159
        except BaseException:
            self.p.terminate();self.p.wait(timeout=10);raise
        return self

    def check(self,a,marker,totals,detected,mito):
        self.p.stdin.write(struct.pack('<i',len(a)))
        self.p.stdin.write(np.asarray(a,dtype='<f8').tobytes(order='C'));self.p.stdin.flush()
        values=np.frombuffer(self.read_exact(self.nc*(self.k+5)*8),dtype='<f8')
        n=self.nc*self.k;first=values[:n].reshape(self.nc,self.k)
        second=values[n:n+self.nc*3].reshape(self.nc,3);n+=self.nc*3
        assert np.array_equal(first,marker),'R_reference_marker_mismatch'
        assert np.array_equal(second,totals),'R_reference_namespace_mismatch'
        assert np.array_equal(values[n:n+self.nc],detected),'R_reference_detection_mismatch'
        assert np.array_equal(values[n+self.nc:],mito),'R_reference_mito_mismatch'

    def __exit__(self,typ,value,tb):
        if typ is not None:
            self.p.terminate()
        else:
            self.p.stdin.write(struct.pack('<i',0));self.p.stdin.flush();self.p.stdin.close()
        try:self.p.wait(timeout=15)
        except subprocess.TimeoutExpired:self.p.kill();self.p.wait(timeout=5)
        if typ is None:assert self.p.returncode==0,'R_reference_failed_exit'
        self.p.stdout.close();self.p.stderr.close()

args <- commandArgs(trailingOnly=TRUE)
library(NMF)
f <- readRDS(args[1]); model <- if(is(f,'NMFfit')) fit(f) else f
cat('object_class', class(f), '\nmodel_class',class(model),'\nslots',slotNames(model),'\n')
W <- basis(model);H<-coef(model);cat('W_dim',dim(W),'H_dim',dim(H),'\n')
if ('theta' %in% slotNames(model)) {theta<-slot(model,'theta');S<-(1-theta)*diag(ncol(W))+theta/ncol(W)*matrix(1,ncol(W),ncol(W));cat('theta',theta,'S_dim',dim(S),'S_row_sums',range(rowSums(S)),'\n');write.table(S,args[2],sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)}

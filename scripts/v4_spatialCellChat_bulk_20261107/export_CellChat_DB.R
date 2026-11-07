args<-commandArgs(trailingOnly=TRUE);r<-normalizePath(args[1]);out<-file.path(r,'sources/CellChatDB_2_2_0');stopifnot(!dir.exists(out));dir.create(out)
suppressPackageStartupMessages(library(CellChat));db<-CellChatDB.human
for(k in c('interaction','complex','cofactor')){a<-db[[k]];a$key<-rownames(a);write.table(a,file.path(out,paste0(k,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)}
saveRDS(db,file.path(out,'database_PRIVATE.rds'));jsonlite::write_json(list(package=as.character(packageVersion('CellChat')),n_total_interactions=nrow(db$interaction)),file.path(out,'complete.json'),auto_unbox=TRUE,pretty=TRUE)

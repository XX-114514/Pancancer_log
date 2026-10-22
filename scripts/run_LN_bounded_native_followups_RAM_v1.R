# Prospective finite scope amendment. No change to the original RAM deadline.
run_LN_bounded_native_followups_RAM <- function(shared_env, root, plan,
    clock=function()as.numeric(Sys.time()), pause=function()Sys.sleep(2)) {
  stopifnot(identical(as.integer(plan$maximum_followups),2L),
    identical(unlist(plan$scopes,use.names=FALSE),c("spatial_effects","independent_recalculation")))
  sha_valid<-function(z)is.character(z)&&length(z)==1L&&!is.na(z)&&grepl("^[a-f0-9]{64}$",z)
  stopifnot(length(plan$original_deadline)==1L,length(plan$amended_deadline)==1L,is.finite(plan$original_deadline),is.finite(plan$amended_deadline))
  deadline<-min(as.numeric(plan$original_deadline),as.numeric(plan$amended_deadline))
  work<-new.env(parent=shared_env);work$analysis_deadline<-deadline
  finished<-0L;family<-NULL;failures<-0L
  write_status<-function(x,name)writeLines(jsonlite::toJSON(x,auto_unbox=TRUE,null="null"),file.path(root,name))
  while(finished<2L&&clock()<deadline&&!file.exists(file.path(root,"release.json"))) {
    i<-finished+1L;meta<-file.path(root,sprintf("analysis%03d.json",i));code<-file.path(root,sprintf("analysis%03d.R",i))
    if(file.exists(meta)&&file.exists(code)) {
      checked<-tryCatch({
        c<-jsonlite::fromJSON(meta);stopifnot(sha_valid(c$code_SHA256),sha_valid(c$family_contract_SHA256),
          identical(c$scope,plan$scopes[[i]]),identical(c$input_SHA256,plan$input_SHA256),
          identical(c$amendment_contract_SHA256,plan$amendment_contract_SHA256),identical(c$qualification_code_SHA256,plan$qualification_code_SHA256))
        sha<-strsplit(system2("sha256sum",shQuote(code),stdout=TRUE),"[[:space:]]+")[[1]][1]
        stopifnot(identical(sha,c$code_SHA256));if(!is.null(family))stopifnot(identical(family,c$family_contract_SHA256))
        c
      },error=function(e)NULL)
      if(is.null(checked)) {
        write_status(list(status="BOUNDED_FOLLOWUP_CONTRACT_OR_SHA_REJECTED",followup=i),"private_bounded_followup_rejected.json")
        break
      }
      if(is.null(family))family<-checked$family_contract_SHA256
      work$control_root<-root
      ok<-tryCatch({sys.source(code,envir=work);TRUE},error=function(e){
        write_status(list(status="FROZEN_FOLLOWUP_FAILED",followup=i,error_class=class(e)[1]),sprintf("private_analysis%03d_failed.json",i));FALSE})
      if(!ok)failures<-failures+1L
      write_status(list(status=if(ok)"FROZEN_FOLLOWUP_COMPLETE" else "FROZEN_FOLLOWUP_FAILED",followup=i,
        scope=checked$scope,code_SHA256=checked$code_SHA256,family_contract_SHA256=family),sprintf("private_analysis%03d_complete.json",i))
      finished<-i
    }
    write_status(list(status="BOUNDED_NATIVE_FOLLOWUP_WAIT",completed_followups=finished,maximum_followups=2,
      failed_followups=failures,original_deadline_not_extended=TRUE,observed_UTC=format(Sys.time(),tz="UTC",usetz=TRUE)),"private_bounded_followup_heartbeat.json")
    if(finished<2L)pause()
  }
  write_status(list(status="BOUNDED_NATIVE_FOLLOWUP_ENDED",completed_followups=finished,maximum_followups=2,
    failed_followups=failures,original_deadline_not_extended=TRUE,release_requested=file.exists(file.path(root,"release.json"))),"private_bounded_followup_end.json")
  invisible(list(completed_followups=finished,failed_followups=failures))
}

args <- commandArgs(trailingOnly=TRUE)
d <- read.delim(file.path(args[1], 'all336_fine_IFNG_models.tsv'), check.names=FALSE)
r <- d[d$pair_id=='PAIR015' & d$test_id=='POOLED' & d$program %in% c('HALLMARK_INTERFERON_GAMMA_RESPONSE','IFNG_RESPONSE_MINUS_IFNA'),]
r <- r[order(r$program, r$score, r$control),]
labels <- paste(ifelse(r$program=='IFNG_RESPONSE_MINUS_IFNA','Gamma minus alpha','Gamma'),ifelse(r$score=='meanlog','Mean log','Gene rank'),ifelse(r$control=='identity_abundance_quality_receptor','Base','+ shared IFN'),sep=' | ')
pdf(file.path(args[1],'PAIR015_program_robustness.pdf'), width=12, height=6)
par(mar=c(6.5,20,4.5,3))
y <- rev(seq_len(nrow(r))); cols <- ifelse(r$BH_exploratory<.05,'#ba3e33','#325a79')
plot(r$beta,y,xlim=range(c(r$lower,r$upper,0)),ylim=c(.5,nrow(r)+.7),pch=19,col=cols,yaxt='n',xlab='Standardized association (95% CR2 interval)',ylab='',main='CM01 / CM03: CD8 NKG7 9 -> Macrophage IFI6 18')
abline(v=0,lty=2,col='grey60'); segments(r$lower,y,r$upper,y,col=cols,lwd=2);axis(2,at=y,labels=labels,las=1,cex.axis=.85)
mtext('107 libraries / 60 source-scoped patient blocks; CRC + glioma; exploratory, not independent validation',side=1,line=4.6,cex=.85)
mtext('Red: BH-supported only in base mean-log model; all BY tests unsupported',side=3,line=.5,cex=.85)
dev.off()

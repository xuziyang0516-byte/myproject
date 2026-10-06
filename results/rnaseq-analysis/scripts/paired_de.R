suppressPackageStartupMessages(library(DESeq2))
base <- Sys.getenv('RNASEQ_PROJECT','/workspace/rnaseq-analysis')
for (kind in c('BMDC','PBMC')) {
 prepared <- readRDS(file.path(base,'results/qc',paste0(kind,'_prepared.rds')))
 mat <- counts(prepared$dds)
 group <- factor(sub('_.*','',colnames(mat)),levels=c('Ctrl','Y180','IMP'))
 subject <- factor(sub('.*_','',colnames(mat)))
 dds <- DESeqDataSetFromMatrix(mat,data.frame(subject,group,row.names=colnames(mat)),~subject+group)
 dds <- DESeq(dds,quiet=TRUE,parallel=FALSE)
 vst <- assay(varianceStabilizingTransformation(dds,blind=FALSE))
 dir <- file.path(base,'results',kind); dir.create(dir,showWarnings=FALSE)
 labels <- setNames(prepared$genes$GeneName,prepared$genes$Geneid)
 coefnames <- resultsNames(dds)
 makecontrast <- function(y=0,i=0) { v <- rep(0,length(coefnames)); names(v)<-coefnames; v['group_Y180_vs_Ctrl']<-y;v['group_IMP_vs_Ctrl']<-i; unname(v) }
 contrasts <- list(Y180_vs_Ctrl=makecontrast(y=1),IMP321_vs_Ctrl=makecontrast(i=1),Y180_vs_IMP321=makecontrast(y=1,i=-1),Ctrl_vs_others=makecontrast(y=-.5,i=-.5),Y180_vs_others=makecontrast(y=1,i=-.5),IMP321_vs_others=makecontrast(y=-.5,i=1))
 summary <- list(); reslist <- list()
 for (name in names(contrasts)) {
  # Common ORA universe: all count-filtered genes having a valid Wald P value.
  # Disable independent filtering to keep one explicit tested-gene universe.
  res <- as.data.frame(results(dds,contrast=contrasts[[name]],alpha=.05,independentFiltering=FALSE))
  res$Geneid<-rownames(res);res$GeneName<-labels[res$Geneid]
  res$DEG<-!is.na(res$padj)&res$padj<.05&abs(res$log2FoldChange)>=1
  res$Direction<-ifelse(res$log2FoldChange>0,'Up','Down')
  write.table(res,file.path(dir,paste0(name,'_DE.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
  reslist[[name]]<-res
  summary[[name]]<-data.frame(dataset=kind,contrast=name,tested=sum(!is.na(res$pvalue)),DEG=sum(res$DEG),Up=sum(res$DEG&res$Direction=='Up'),Down=sum(res$DEG&res$Direction=='Down'))
 }
 write.table(do.call(rbind,summary),file.path(dir,'DE_summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 saveRDS(list(dds=dds,vst=vst,results=reslist,labels=labels,group=group,subject=subject),file.path(dir,'paired_analysis.rds'))
 print(do.call(rbind,summary))
}

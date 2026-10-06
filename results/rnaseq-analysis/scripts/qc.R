suppressPackageStartupMessages({library(DESeq2); library(ggplot2)})
base <- Sys.getenv('RNASEQ_PROJECT','/workspace/rnaseq-analysis')
for (kind in c('BMDC','PBMC')) {
 path <- list.files(file.path(base,'inputs'),pattern=paste0('^',kind,'\\.reads.count\\.txt$'),recursive=TRUE,full.names=TRUE)
 tab <- read.delim(path, check.names=FALSE)
 mat <- as.matrix(tab[,3:ncol(tab)]); rownames(mat) <- tab$Geneid
 stopifnot(!anyNA(mat),all(mat>=0),all(mat==round(mat)),!anyDuplicated(rownames(mat)))
 group <- factor(sub('_.*','',colnames(mat)),levels=c('Ctrl','Y180','IMP'))
 keep <- rowSums(mat>=10)>=3
 qc <- data.frame(sample=colnames(mat),group=group,library_size=colSums(mat),detected_genes=colSums(mat>0),retained_genes=sum(keep))
 write.table(qc,file.path(base,'results/qc',paste0(kind,'_sample_qc.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
 dds <- DESeqDataSetFromMatrix(mat[keep,],data.frame(group,row.names=colnames(mat)),~1)
 dds <- estimateSizeFactors(dds)
 vs <- varianceStabilizingTransformation(dds,blind=TRUE)
 saveRDS(list(dds=dds,vst=assay(vs),genes=tab,keep=keep),file.path(base,'results/qc',paste0(kind,'_prepared.rds')))
 pcs <- prcomp(t(assay(vs))); percent <- round(100*pcs$sdev^2/sum(pcs$sdev^2),1)
 df <- data.frame(pcs$x[,1:2],sample=rownames(pcs$x),group=group)
 p <- ggplot(df,aes(PC1,PC2,color=group,label=sample))+geom_point(size=3)+geom_text(vjust=-.7,size=3)+theme_bw(base_size=12)+labs(title=paste(kind,'sample PCA'),x=paste0('PC1 (',percent[1],'%)'),y=paste0('PC2 (',percent[2],'%)'))
 ggsave(file.path(base,'results/qc',paste0(kind,'_PCA.pdf')),p,width=7,height=5)
 cat(kind,'retained',sum(keep),'of',nrow(mat),'size factors',sizeFactors(dds),'\n')
}

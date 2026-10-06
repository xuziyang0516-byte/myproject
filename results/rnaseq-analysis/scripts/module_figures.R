suppressPackageStartupMessages(library(ggplot2))
base <- Sys.getenv('RNASEQ_PROJECT','/workspace/rnaseq-analysis')
write_tsv <- function(x,file) {
 x[] <- lapply(x,function(col) if(is.character(col)||is.factor(col)) gsub('\n',' | ',as.character(col),fixed=TRUE) else col)
 write.table(x,file,sep='\t',row.names=FALSE,quote=FALSE,na='NA')
}
wrap <- function(s,width=47) vapply(s,function(x) paste(strwrap(x,width=width),collapse='\n'),character(1))
for (kind in c('BMDC','PBMC')) {
 dir <- file.path(base,'results',kind)
 a <- readRDS(file.path(dir,'paired_analysis.rds'))
 selected <- read.delim(file.path(base,'results',paste0(kind,'_selected_pathways.tsv')),check.names=FALSE)
 excel <- read.delim(file.path(base,'results',paste0(kind,'_excel_rows.tsv')),check.names=FALSE)
 # Use the union of supplied UP-module genes for each selected ID, not a fabricated full pathway gene set.
 allowed <- if(kind=='BMDC') c('Sheet1','Sheet3','Sheet4') else c('Sheet2','Sheet3','Sheet4','Sheet5')
 rows <- excel[excel$Sheet %in% allowed & excel$Direction=='Up' & excel$ID %in% selected$ID,]
 # Reject copied rows with no exact species-specific symbol match; keep the audit trail.
 rows$SpeciesMatchedGenes <- vapply(strsplit(rows$GeneName,'/'),function(g) sum(g %in% unname(a$labels)),integer(1))
 write_tsv(rows[rows$SpeciesMatchedGenes==0,],file.path(dir,'excluded_cross_species_or_unmapped_rows.tsv'))
 rows <- rows[rows$SpeciesMatchedGenes>0,]
 # Match gene symbols exactly; keep the most highly expressed Ensembl ID when a symbol is ambiguous.
 mean_count <- rowMeans(DESeq2::counts(a$dds,normalized=TRUE))
 mapping <- data.frame(Geneid=rownames(a$vst),GeneName=unname(a$labels[rownames(a$vst)]),MeanNormalizedCount=mean_count)
 mapping <- mapping[order(-mapping$MeanNormalizedCount),]
 write_tsv(mapping[duplicated(mapping$GeneName)|duplicated(mapping$GeneName,fromLast=TRUE),],file.path(dir,'ambiguous_symbol_mapping.tsv'))
 mapping <- mapping[!duplicated(mapping$GeneName)&nzchar(mapping$GeneName),]
 scores <- list(); detected <- list(); modules <- list(); matched <- list(); matchsummary <- list()
 norm <- DESeq2::counts(a$dds,normalized=TRUE)
 for (j in seq_len(nrow(selected))) {
  id <- selected$ID[j]
  genes <- unique(unlist(strsplit(rows$GeneName[rows$ID==id],'/')))
  genes <- genes[!is.na(genes)&nzchar(genes)]
  ids <- mapping$Geneid[match(genes,mapping$GeneName)]; names(ids)<-genes
  valid <- !is.na(ids)
  zz <- a$vst[ids[valid],,drop=FALSE]
  varying <- apply(zz,1,sd)>0
  ids2 <- ids[valid][varying]
  if(length(ids2)<2) stop('Insufficient module genes: ',kind,' ',id)
  zz <- t(scale(t(a$vst[ids2,,drop=FALSE])))
  scores[[id]] <- colMeans(zz)
  detected[[id]] <- sapply(c('Ctrl','Y180','IMP'),function(g) sum(rowSums(norm[ids2,a$group==g,drop=FALSE]>=50)>=2))
  modules[[id]]<-data.frame(ID=id,GeneName=names(ids),Geneid=unname(ids),Matched=valid,UsedInScore=names(ids)%in%names(ids2))
  matchsummary[[id]]<-data.frame(ID=id,Description=selected$Description[j],ExcelModuleGenes=length(genes),MappedGenes=sum(valid),VariableGenes=length(ids2),MissingGenes=paste(genes[!valid],collapse='/'))
 }
 scoremat <- do.call(rbind,scores)
 df0 <- data.frame(subject=a$subject,group=a$group)
 out <- list(); pairout <- list(); omnibus <- list(); sampleout <- list()
 for (id in selected$ID) {
  df <- transform(df0,score=as.numeric(scores[[id]]))
  fit <- lm(score~subject+group,df)
  fit0 <- lm(score~subject,df)
  omnibus[[id]] <- data.frame(ID=id,Pvalue=anova(fit0,fit)$'Pr(>F)'[2])
  co <- coef(fit); vc<-vcov(fit); rdf<-df.residual(fit)
  contrasts <- list(Ctrl=c(-.5,-.5),Y180=c(1,-.5),IMP=c(-.5,1))
  paircons <- list(Y180_vs_Ctrl=c(1,0),IMP321_vs_Ctrl=c(0,1),Y180_vs_IMP321=c(1,-1))
  calc <- function(v) {
   cc <- rep(0,length(co));names(cc)<-names(co);cc['groupY180']<-v[1];cc['groupIMP']<-v[2]
   effect<-sum(cc*co);se<-sqrt(drop(t(cc)%*%vc%*%cc));p<-if(se>0) 2*pt(-abs(effect/se),rdf) else if(effect==0) 1 else NA_real_
   data.frame(Effect=effect,SE=se,CI_low=effect-qt(.975,rdf)*se,CI_high=effect+qt(.975,rdf)*se,Pvalue=p,ResidualDF=rdf)
  }
  for(g in names(contrasts)) out[[paste(id,g)]]<-cbind(data.frame(ID=id,Group=g,Count=detected[[id]][g],MeanScore=mean(df$score[df$group==g])),calc(contrasts[[g]]))
  for(cmp in names(paircons)) pairout[[paste(id,cmp)]]<-cbind(data.frame(ID=id,Comparison=cmp),calc(paircons[[cmp]]))
  sampleout[[id]]<-data.frame(ID=id,Sample=colnames(a$vst),Subject=a$subject,Group=as.character(a$group),Score=as.numeric(scores[[id]]))
 }
 tab <- do.call(rbind,out);tab$Padjust<-p.adjust(tab$Pvalue,'BH')
 tab$Direction<-ifelse(tab$Effect>0,'Higher','Lower');tab$Significant<-!is.na(tab$Padjust)&tab$Padjust<.05
 pair <- do.call(rbind,pairout);pair$Padjust<-p.adjust(pair$Pvalue,'BH')
 omni <- do.call(rbind,omnibus);omni$Padjust<-p.adjust(omni$Pvalue,'BH')
 tab$Description<-selected$Description[match(tab$ID,selected$ID)]
 tab$Label<-paste(tab$ID,wrap(tab$Description),sep=', ')
 order_labels <- unique(tab$Label[match(selected$ID,tab$ID)])
 tab$Label<-factor(tab$Label,levels=rev(order_labels))
 tab$Group<-factor(tab$Group,levels=c('Ctrl','Y180','IMP'),labels=c('CTRL','Y180','IMP321'))
 tab$ColorScore<- -log10(pmax(tab$Padjust,1e-300))
 write_tsv(tab,file.path(dir,'three_group_module_statistics.tsv'))
 write_tsv(pair,file.path(dir,'paired_module_comparisons.tsv'))
 write_tsv(omni,file.path(dir,'module_omnibus_statistics.tsv'))
 write_tsv(do.call(rbind,modules),file.path(dir,'module_gene_mapping.tsv'))
 write_tsv(do.call(rbind,matchsummary),file.path(dir,'module_mapping_summary.tsv'))
 sample <- do.call(rbind,sampleout)
 write_tsv(sample,file.path(dir,'module_sample_scores.tsv'))
 theme <- theme_bw(base_size=12)+theme(panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),axis.title=element_blank(),axis.text.y=element_text(color='black',size=10),axis.text.x=element_text(color='black',size=12),plot.title=element_text(face='bold',size=15),plot.subtitle=element_text(size=10),plot.caption=element_text(hjust=0,size=9),legend.title=element_text(size=10))
 breaks<-c(.05,.01,.001,.0001); usable<-breaks[-log10(breaks)<=max(tab$ColorScore,na.rm=TRUE)+.2]
 if(length(usable)==0) usable<-.05
 p <- ggplot(tab,aes(Group,Label))+geom_point(data=tab[!tab$Significant,],aes(size=Count),shape=21,fill='white',color='#999999',stroke=.6)+geom_point(data=tab[tab$Significant,],aes(size=Count,fill=ColorScore),shape=21,color='#555555',stroke=.5)+geom_text(aes(label=ifelse(Direction=='Higher','+','-')),size=3,color='black')+scale_size_area(max_size=11,name='Count')+scale_fill_gradientn(colours=c('#fee5d9','#fb6a4a','#a50f15'),limits=c(-log10(.05),max(-log10(.05)+.01,max(tab$ColorScore,na.rm=TRUE))),breaks=-log10(usable),labels=format(usable,scientific=TRUE),name='Adjusted P value')+scale_x_discrete(drop=FALSE,expand=expansion(add=.7))+scale_y_discrete(limits=rev(order_labels),drop=FALSE)+theme+labs(title=paste(kind,'| Selected-module statistical dot plot'),subtitle='Each group vs the mean of the other two; paired subject + group model',caption=paste('Count: selected-module genes with normalized count >= 50 in >= 2 of 3 subjects.\nFilled = BH-adjusted P < 0.05; hollow = not significant; +/- inside circles indicate higher/lower scores.\nUser-selected DEG modules; exploratory post-selection statistics, not GO/KEGG ORA.'))
 ggsave(file.path(dir,paste0(kind,'_three_groups_Pvalue_Count.pdf')),p,width=12,height=8.8)
 ggsave(file.path(dir,paste0(kind,'_three_groups_Pvalue_Count.png')),p,width=12,height=8.8,dpi=300)
 # SCI-style companion: show every paired biological replicate, not just three group means.
 sample$Label<-factor(paste(sample$ID,wrap(selected$Description[match(sample$ID,selected$ID)]),sep=', '),levels=rev(order_labels))
 order_samples<-unlist(lapply(c('Ctrl','Y180','IMP'),function(g) paste0(g,'_',1:3)))
 sample$Sample<-factor(sample$Sample,levels=order_samples)
 sample$Group<-factor(sample$Group,levels=c('Ctrl','Y180','IMP'),labels=c('CTRL','Y180','IMP321'))
 hp<-ggplot(sample,aes(Sample,Label,fill=Score))+geom_tile(color='white',linewidth=.6)+scale_fill_gradient2(low='#2166ac',mid='white',high='#b2182b',midpoint=0,name='Mean gene-wise\nVST z-score')+facet_grid(.~Group,scales='free_x',space='free_x')+theme+theme(axis.text.x=element_text(angle=45,hjust=1),strip.background=element_rect(fill='grey95'))+labs(title=paste(kind,'| Paired selected-module expression'),subtitle='All nine biological samples; same module definition as the statistical dot plot',caption='Modules use the union of supplied UP-DEG gene lists for each selected pathway ID; full pathway membership is not implied.')
 ggsave(file.path(dir,paste0(kind,'_paired_module_heatmap.pdf')),hp,width=12,height=8.2)
 ggsave(file.path(dir,paste0(kind,'_paired_module_heatmap.png')),hp,width=12,height=8.2,dpi=300)
 traj<-ggplot(sample,aes(Group,Score,group=Subject,color=Subject))+geom_line(linewidth=.7)+geom_point(size=2)+facet_wrap(~Label,ncol=2,scales='free_y')+theme_bw(base_size=10)+theme(strip.text=element_text(size=8),axis.title.x=element_blank())+labs(title=paste(kind,'| Paired module trajectories'),y='Module score',caption='Lines connect the same subject across all three conditions; three subjects per condition.')
 ggsave(file.path(dir,paste0(kind,'_paired_trajectories.pdf')),traj,width=10,height=14)
 cat(kind,'significant group-module cells:',sum(tab$Significant),'of',nrow(tab),'\n'); print(aggregate(Significant~Group,tab,sum))
 # All-gene pairwise results for context: demonstrate whether Y180/IMP321 differences are supported.
 volc<-do.call(rbind,lapply(c('Y180_vs_Ctrl','IMP321_vs_Ctrl','Y180_vs_IMP321'),function(n){x<-a$results[[n]];x$Comparison<-n;x}))
 volc$Category<-ifelse(volc$DEG,ifelse(volc$log2FoldChange>0,'Up','Down'),'Not DEG')
 vp<-ggplot(volc,aes(log2FoldChange,-log10(pmax(padj,1e-300)),color=Category))+geom_point(size=.6,alpha=.5)+geom_vline(xintercept=c(-1,1),linetype=2,color='grey60')+geom_hline(yintercept=-log10(.05),linetype=2,color='grey60')+facet_wrap(~Comparison,nrow=1)+scale_color_manual(values=c(Up='#b2182b',Down='#2166ac','Not DEG'='#bbbbbb'))+theme_bw(base_size=11)+labs(title=paste(kind,'| Paired differential expression'),x='log2 fold change',y='-log10(BH-adjusted P)',caption='Unshrunken DESeq2 fold changes; DEG: adjusted P < 0.05 and absolute log2FC >= 1. Independent filtering disabled for an explicit common tested-gene universe.')
 ggsave(file.path(dir,paste0(kind,'_paired_volcano.pdf')),vp,width=13,height=4.5)
}

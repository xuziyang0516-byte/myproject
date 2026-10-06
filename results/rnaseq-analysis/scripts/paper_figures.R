suppressPackageStartupMessages({library(ggplot2);library(grid);library(cowplot)})
base<-Sys.getenv('RNASEQ_PROJECT','/workspace/rnaseq-analysis')
for(kind in c('BMDC','PBMC')) {
 dir<-file.path(base,'results',kind);a<-readRDS(file.path(dir,'paired_analysis.rds'))
 selected<-read.delim(file.path(base,'results',paste0(kind,'_selected_pathways.tsv')))
 sets0<-read.delim(file.path(dir,'full_pathway_gene_sets.tsv'),colClasses='character')
 sets<-lapply(split(sets0$EntrezID,sets0$ID),unique)
 mapping<-read.delim(file.path(dir,'ensembl_entrez_mapping.tsv'),colClasses='character')
 ora<-list()
 for(cmp in c('Y180_vs_Ctrl','IMP321_vs_Ctrl')) {
  de<-a$results[[cmp]]
  U<-unique(mapping$ENTREZID[mapping$ENSEMBL %in% de$Geneid[!is.na(de$pvalue)]])
  D<-intersect(U,unique(mapping$ENTREZID[mapping$ENSEMBL %in% de$Geneid[de$DEG&de$Direction=='Up']]))
  candidate<-lapply(sets,intersect,y=U);candidate<-candidate[lengths(candidate)>=5 & lengths(candidate)<=1000]
  hits<-lapply(candidate,intersect,y=D)
  pv<-vapply(seq_along(candidate),function(i)phyper(length(hits[[i]])-1,length(candidate[[i]]),length(U)-length(candidate[[i]]),length(D),lower.tail=FALSE),numeric(1))
  ora[[cmp]]<-data.frame(Contrast=cmp,ID=names(candidate),Count=lengths(hits),GeneSetSize=lengths(candidate),DEGTotal=length(D),Universe=length(U),Pvalue=pv)
 }
 ora<-do.call(rbind,ora);ora$Padjust<-p.adjust(ora$Pvalue,'BH')
 write.table(ora,file.path(dir,'paired_UP_ORA_all_terms.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 sel<-ora[ora$ID%in%selected$ID,]
 descriptions<-vapply(selected$Description,function(x)paste(strwrap(x,45),collapse='\n'),character(1))
 labs<-paste(selected$ID,descriptions,sep=', ')
 sel$Label<-factor(labs[match(sel$ID,selected$ID)],levels=rev(labs))
 sel$Comparison<-factor(sel$Contrast,levels=c('Y180_vs_Ctrl','IMP321_vs_Ctrl'),labels=c('Y180 vs CTRL','IMP321 vs CTRL'))
 sel$Significant<-sel$Padjust<.05;sel$ColorScore<--log10(pmax(sel$Padjust,1e-300))
 write.table(transform(sel,Label=gsub('\n',' | ',as.character(Label),fixed=TRUE)),file.path(dir,'paired_UP_ORA_selected_terms.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 base_theme<-theme_bw(base_size=10)+theme(axis.title=element_blank(),panel.grid.minor=element_blank(),axis.text.y=element_text(color='black',size=8),axis.text.x=element_text(color='black'),plot.title=element_text(face='bold',size=12),plot.caption=element_text(hjust=0,size=8))
 p<-ggplot(sel,aes(Comparison,Label))+geom_point(data=sel[!sel$Significant&sel$Count>0,],aes(size=Count),shape=21,fill='white',color='grey60')+geom_point(data=sel[sel$Significant&sel$Count>0,],aes(size=Count,fill=ColorScore),shape=21,color='grey40')+geom_point(data=sel[sel$Count==0,],shape=4,size=2,color='grey60')+scale_y_discrete(limits=rev(labs),drop=FALSE)+scale_x_discrete(limits=c('Y180 vs CTRL','IMP321 vs CTRL'))+scale_size_area(max_size=8,name='Count')+scale_fill_gradient(low='#fee5d9',high='#a50f15',limits=c(-log10(.05),max(2,ceiling(max(sel$ColorScore)))),breaks=-log10(c(.05,10^(-seq(3,max(3,ceiling(max(sel$ColorScore))),by=3)))),labels=format(c(.05,10^(-seq(3,max(3,ceiling(max(sel$ColorScore))),by=3))),scientific=TRUE),name='Adjusted P value')+base_theme+theme(axis.text.x=element_text(angle=30,hjust=1),axis.text.y=element_blank(),axis.ticks.y=element_blank())+labs(title='B | Treatment-versus-control UP ORA',caption='Count = pathway-overlapping UP DEGs. Filled = ORA FDR < 0.05.\nBH correction covers all eligible GO/KEGG sets and both contrasts.')
 # All expressed members of the full annotated pathways, rather than the workbook's selected DEGs.
 z<-t(scale(t(a$vst)));valid<-rowSums(!is.finite(z))==0
 z<-z[valid,,drop=FALSE]
 X<-model.matrix(~subject+group,data=data.frame(subject=a$subject,group=a$group))
 fitted<-lm.fit(X,t(z))
 subject_columns<-grep('^subject',colnames(X))
 z<-z-t(X[,subject_columns,drop=FALSE]%*%fitted$coefficients[subject_columns,,drop=FALSE])
 z<-z-rowMeans(z) # Center the removed subject contribution; retain condition effects.

 heat<-list();audit<-list()
 for(id in selected$ID) {
  ens<-unique(mapping$ENSEMBL[mapping$ENTREZID %in% sets[[id]]]);ens<-intersect(ens,rownames(z))
  if(!length(ens))stop('No expressed annotated genes: ',id)
  score<-colMeans(z[ens,,drop=FALSE])
  heat[[id]]<-data.frame(ID=id,Sample=names(score),Score=unname(score),Group=sub('_.*','',names(score)),Genes=length(ens))
  audit[[id]]<-data.frame(ID=id,Geneid=ens)
 }
 heat<-do.call(rbind,heat)
 heat$Label<-factor(labs[match(heat$ID,selected$ID)],levels=rev(labs))
 heat$Sample<-factor(heat$Sample,levels=unlist(lapply(c('Ctrl','Y180','IMP'),function(g)paste0(g,'_',1:3))))
 heat$Group<-factor(heat$Group,levels=c('Ctrl','Y180','IMP'),labels=c('CTRL','Y180','IMP321'))
 write.table(heat,file.path(dir,'full_pathway_sample_expression.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 write.table(do.call(rbind,audit),file.path(dir,'full_pathway_expression_gene_mapping.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 hp<-ggplot(heat,aes(Sample,Label,fill=Score))+geom_tile(color='white',linewidth=.5)+scale_y_discrete(limits=rev(labs),drop=FALSE)+scale_fill_gradient2(low='#2166ac',mid='white',high='#b2182b',midpoint=0,name='Subject-adjusted\nmean VST z-score')+facet_grid(.~Group,scales='free_x',space='free_x')+base_theme+theme(axis.text.x=element_text(angle=45,hjust=1),strip.background=element_rect(fill='grey95'))+labs(title='A | Full-pathway expression in all paired samples',caption='All expressed annotated members; subject effects removed with ~subject + group, retaining group effects.\nEach column is one paired biological sample; expression remains descriptive and does not imply activation.')
 ggsave(file.path(dir,paste0(kind,'_full_pathway_heatmap.pdf')),hp,width=12,height=8)
 ggsave(file.path(dir,paste0(kind,'_full_pathway_heatmap.png')),hp,width=12,height=8,dpi=300)
 combined<-cowplot::plot_grid(hp,p,nrow=1,rel_widths=c(.64,.36),align='h',axis='tb')
 ggsave(file.path(dir,paste0(kind,'_paper_recommended.pdf')),combined,width=17,height=8.5)
 ggsave(file.path(dir,paste0(kind,'_paper_recommended.png')),combined,width=17,height=8.5,dpi=300)
 cat(kind,'treatment UP ORA significant:',sum(sel$Significant),'of',nrow(sel),'\n')
}

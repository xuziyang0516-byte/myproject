suppressPackageStartupMessages(library(ggplot2))
base<-Sys.getenv('RNASEQ_PROJECT','/workspace/rnaseq-analysis')
all<-lapply(c('BMDC','PBMC'),function(kind)read.delim(file.path(base,'results',kind,'full_ora_selected_terms.tsv')))
maxlog<-max(vapply(all,function(d)max(-log10(pmax(d$Padjust,1e-300)),na.rm=TRUE),numeric(1)))
upper<-max(2,ceiling(maxlog))
pbreaks<-c(.05,10^(-seq(3,upper,by=3)))
for(i in seq_along(all)) {
 kind<-c('BMDC','PBMC')[i];dir<-file.path(base,'results',kind);d<-all[[i]]
 selected<-read.delim(file.path(base,'results',paste0(kind,'_selected_pathways.tsv')))
 desc<-vapply(selected$Description,function(x)paste(strwrap(x,48),collapse='\n'),character(1))
 labels<-paste(selected$ID,desc,sep=', ')
 d$Label<-factor(labels[match(d$ID,selected$ID)],levels=rev(labels))
 d$Group<-factor(d$Contrast,levels=c('Ctrl_vs_others','Y180_vs_others','IMP321_vs_others'),labels=c('CTRL','Y180','IMP321'))
 d$Significant<-!is.na(d$Padjust)&d$Padjust<.05;d$ColorScore<--log10(pmax(d$Padjust,1e-300))
 p<-ggplot(d,aes(Group,Label))+geom_point(data=d[!d$Significant&d$Count>0,],aes(size=Count),shape=21,fill='white',color='grey60')+geom_point(data=d[d$Significant&d$Count>0,],aes(size=Count,fill=ColorScore),shape=21,color='grey40')+geom_point(data=d[d$Count==0,],shape=4,size=2,color='grey60')+scale_y_discrete(limits=rev(labels),drop=FALSE)+scale_x_discrete(expand=expansion(add=.7))+scale_size_area(max_size=11,name='Count')+scale_fill_gradientn(colours=c('#fee5d9','#fb6a4a','#a50f15'),limits=c(-log10(.05),upper),breaks=-log10(pbreaks),labels=format(pbreaks,scientific=TRUE),name='Adjusted P value')+theme_bw(base_size=12)+theme(axis.title=element_blank(),axis.text.y=element_text(color='black',size=10),axis.text.x=element_text(color='black',size=12),panel.grid.minor=element_blank(),plot.title=element_text(face='bold',size=15),plot.caption=element_text(hjust=0,size=9))+labs(title=paste(kind,'| Differential pathway enrichment'),subtitle='Each group vs the mean log expression of the other two; paired subject + group model',caption='Count = pathway-overlapping DEGs, combining up- and downregulated genes; enrichment does not imply activation.\nFilled = ORA BH-adjusted P < 0.05; hollow = not significant; x = zero overlapping DEGs.\nDEG: gene-level BH P < 0.05 and |log2FC| >= 1. ORA BH correction covers all eligible sets and all three contrasts.')
 ggsave(file.path(dir,paste0(kind,'_full_pathway_ORA.pdf')),p,width=12,height=8.8)
 ggsave(file.path(dir,paste0(kind,'_full_pathway_ORA.png')),p,width=12,height=8.8,dpi=300)
}

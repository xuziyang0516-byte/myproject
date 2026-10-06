# Full-pathway ORA. This requires species-specific OrgDb and the KEGG REST service.
# Missing resources produce explicit diagnostics and PARTIAL output names, never fabricated statistics.
suppressPackageStartupMessages({library(AnnotationDbi);library(GO.db);library(ggplot2)})
base<-Sys.getenv('RNASEQ_PROJECT','/workspace/rnaseq-analysis')
for(kind in c('BMDC','PBMC')) {
 dir<-file.path(base,'results',kind);a<-readRDS(file.path(dir,'paired_analysis.rds'))
 selected<-read.delim(file.path(base,'results',paste0(kind,'_selected_pathways.tsv')))
 pkg<-if(kind=='BMDC') 'org.Mm.eg.db' else 'org.Hs.eg.db'
 species<-if(kind=='BMDC') 'mmu' else 'hsa'
 issues<-list(); sets<-list(); mapping<-NULL;go_sets<-list();kegg_sets<-list()
 if(requireNamespace(pkg,quietly=TRUE)) {
  db<-get(pkg,asNamespace(pkg))
  geneids<-rownames(a$vst)
  mapping<-suppressMessages(select(db,keys=geneids,keytype='ENSEMBL',columns='ENTREZID'))
  mapping<-unique(mapping[!is.na(mapping$ENTREZID),])
  gos<-toTable(get(paste0(sub('\\.db$','',pkg),'GO2ALLEGS'),asNamespace(pkg)))
  go_sets<-lapply(split(as.character(gos$gene_id),gos$go_id),unique)
 } else issues[['GO']]<-data.frame(Resource=pkg,Reason='Species-specific annotation package unavailable; GO ORA not run.')
 # KEGG annotations are downloaded through verified HTTPS and cached for reproducibility.
 cache<-file.path(base,'results',paste0(species,'_kegg_gene_pathway.tsv'))
 if(!file.exists(cache)) {
  tryCatch(download.file(paste0('https://rest.kegg.jp/link/pathway/',species),cache,method='libcurl',quiet=TRUE),error=function(e){if(file.exists(cache))unlink(cache);issues[['KEGG']]<<-data.frame(Resource='KEGG REST',Reason=conditionMessage(e))})
 }
 if(file.exists(cache) && file.info(cache)$size>0) {
  kk<-read.delim(cache,header=FALSE,col.names=c('Gene','Pathway'))
  kk$Gene<-sub(paste0('^',species,':'),'',kk$Gene);kk$Pathway<-sub('^path:','',kk$Pathway)
  kegg_sets<-lapply(split(kk$Gene,kk$Pathway),unique)
 }
 if(is.null(mapping)) {
  write.table(do.call(rbind,issues),file.path(dir,'full_ora_missing_resources.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
  cat(kind,': no species annotation; full pathway ORA remains blocked\n');next
 }
 sets<-c(go_sets,kegg_sets)
 namesets<-names(sets)
 cmpnames<-c('Ctrl_vs_others','Y180_vs_others','IMP321_vs_others')
 all<-list()
 for(cmp in cmpnames) {
  de<-a$results[[cmp]]
  tested_ids<-de$Geneid[!is.na(de$pvalue)]
  U<-unique(mapping$ENTREZID[mapping$ENSEMBL %in% tested_ids])
  fgids<-de$Geneid[de$DEG]
  D<-intersect(U,unique(mapping$ENTREZID[mapping$ENSEMBL %in% fgids]))
  candidate<-lapply(sets,intersect,y=U)
  eligible<-lengths(candidate)>=5 & lengths(candidate)<=1000
  candidate<-candidate[eligible]
  hits<-lapply(candidate,intersect,y=D)
  p<-vapply(seq_along(candidate),function(i) phyper(length(hits[[i]])-1,length(candidate[[i]]),length(U)-length(candidate[[i]]),length(D),lower.tail=FALSE),numeric(1))
  all[[cmp]]<-data.frame(Contrast=cmp,ID=names(candidate),Count=lengths(hits),GeneSetSize=lengths(candidate),DEGTotal=length(D),Universe=length(U),Pvalue=p,HitEntrez=vapply(hits,paste,collapse='/',character(1)))
 }
 all<-do.call(rbind,all)
 # One stated family: all eligible annotated GO/KEGG sets across the three planned contrasts.
 all$Padjust<-p.adjust(all$Pvalue,'BH')
 write.table(all,file.path(dir,'full_ora_all_terms.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 selectedtab<-merge(expand.grid(ID=selected$ID,Contrast=cmpnames,stringsAsFactors=FALSE),all,by=c('ID','Contrast'),all.x=TRUE,sort=FALSE)
 selectedtab$Description<-selected$Description[match(selectedtab$ID,selected$ID)]
 selectedtab$MissingReason<-ifelse(is.na(selectedtab$Pvalue),'No annotation or gene-set size outside [5,1000]','')
 write.table(data.frame(ID=rep(names(sets),lengths(sets)),EntrezID=unlist(sets,use.names=FALSE)),file.path(dir,'full_pathway_gene_sets.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 write.table(mapping,file.path(dir,'ensembl_entrez_mapping.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 write.table(selectedtab,file.path(dir,'full_ora_selected_terms.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 if(!length(issues) && file.exists(file.path(dir,'full_ora_missing_resources.tsv'))) unlink(file.path(dir,'full_ora_missing_resources.tsv'))
 if(length(issues))write.table(do.call(rbind,issues),file.path(dir,'full_ora_missing_resources.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 # Use fixed pathway order even when non-significant layers have different row subsets.
 lab<-paste(selected$ID,selected$Description,sep=', ')
 selectedtab$Label<-factor(paste(selectedtab$ID,selectedtab$Description,sep=', '),levels=rev(lab))
 selectedtab$Group<-factor(selectedtab$Contrast,levels=cmpnames,labels=c('CTRL','Y180','IMP321'))
 selectedtab$Significant<-!is.na(selectedtab$Padjust)&selectedtab$Padjust<.05
 selectedtab$ColorScore<--log10(pmax(selectedtab$Padjust,1e-300))
 complete<-!anyNA(selectedtab$Pvalue)
 stem<-paste0(kind,if(complete)'_full_pathway_ORA' else '_PARTIAL_full_pathway_ORA')
 good<-selectedtab[!is.na(selectedtab$Pvalue)&selectedtab$Count>0,]
 zero<-selectedtab[!is.na(selectedtab$Pvalue)&selectedtab$Count==0,]
 miss<-selectedtab[is.na(selectedtab$Pvalue),]
 p<-ggplot(selectedtab,aes(Group,Label))+geom_point(data=good[!good$Significant,],aes(size=Count),shape=21,fill='white',color='grey60')+geom_point(data=good[good$Significant,],aes(size=Count,fill=ColorScore),shape=21,color='grey40')+geom_point(data=zero,shape=4,size=2,color='grey60')+geom_text(data=miss,label='NA',color='grey60',size=3)+scale_size_area(max_size=11,name='Count')+scale_fill_gradient(low='#fee5d9',high='#a50f15',limits=c(-log10(.05),30),oob=scales::squish,breaks=-log10(c(.05,1e-5,1e-10,1e-20,1e-30)),labels=c('0.05','1e-05','1e-10','1e-20','<=1e-30'),name='Adjusted P value')+scale_y_discrete(limits=rev(lab),drop=FALSE)+theme_bw(base_size=11)+theme(axis.title=element_blank(),panel.grid.minor=element_blank())+labs(title=paste(kind,'| Differential pathway enrichment',if(!complete)'(PARTIAL)' else ''),subtitle='Paired DESeq2: each group vs mean log expression of the other two',caption='Count = pathway-overlapping DEGs, irrespective of direction; not an activation score.\nDEG: BH P < 0.05 and |log2FC| >= 1. ORA BH correction includes all eligible sets across all three contrasts.\nFilled = ORA FDR < 0.05; hollow = not significant; x = zero overlap; NA = unavailable statistics.')
 ggsave(file.path(dir,paste0(stem,'.pdf')),p,width=12,height=8)
 ggsave(file.path(dir,paste0(stem,'.png')),p,width=12,height=8,dpi=300)
 cat(kind,'full-pathway resource coverage:',sum(!is.na(selectedtab$Pvalue)),'/',nrow(selectedtab),'cells; significant:',sum(selectedtab$Significant),'\n')
}

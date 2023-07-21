###GSEA
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
KEGG_df = msigdbr(species = "Homo sapiens",category = "C2",subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name,gene_symbol)
head(KEGG_df)
HALL_df = msigdbr(species = "Homo sapiens",category = "H") %>% 
  dplyr::select(gs_name,gene_symbol)
head(HALL_df)
ALL_df<-rbind(KEGG_df,HALL_df)

###response
result1<-read.table("tumor post response.vs.nonresponse pc count with batch.csv",header = T,sep = ",",check.names = F,row.names = 1)
ge = result1$log2FoldChange
names(ge) = result1$gene
ge = sort(ge,decreasing = T)
head(ge)
em <- GSEA(ge, TERM2GENE = ALL_df,pvalueCutoff = 0.05)
sdf<-em@result
colnames(sdf)
library(stringr)
sdf <- mutate(
  sdf, 
  Description=str_replace_all(str_replace(Description, "([A-Za-z]+)_(.+)", "\\2 (\\1)"), "_", " "))
sdf<-sdf[order(sdf$NES,decreasing = T),]
sdf$Description<-factor(sdf$Description,levels=rev(sdf$Description),ordered = T)

library(ggplot2)
tiff("tumor post response.vs.non GSEA pc count with batch.tiff", units="in", width=12, height=10, res=300,compression = "lzw")
g <- ggplot(sdf, aes(x=Description, y=NES, fill=NES)) +
  geom_bar(stat="identity") +#, alpha=0.7) +
  labs(x="", y="Normalized Enrichment Score", title="") +
  scale_fill_gradientn(colours=hcl.colors(20, "plasma")) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14),
        legend.title=element_text(family="Times", colour="black",size=12))+
  coord_flip()
print(g)
dev.off()

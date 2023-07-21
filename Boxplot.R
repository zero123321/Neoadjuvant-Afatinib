tpm<-read.table("final_pctpm_Combat_convertID.txt",sep = ",",header=T,check.names=F,row.names = 1)

teffect<-c("STAT1","CXCL9","IDO1","PSMB10","LAG3","TIGIT","CXCR6","CCL5","NKG7","CD27","HLA-E",
           "CD274","HLA-DRB1","HLA-DQA1","CMKLR1","PDCD1LG2","CD276","CD8A")
tlsmarker1<-tpm[as.vector(teffect),as.vector(batch$sample)]
tlsmarker1<-as.data.frame(t(as.matrix(tlsmarker1)))
tlsmarker1$treatment<-batch$tissue[match(rownames(tlsmarker1),batch$sample)]
library(reshape2)
ciber<-melt(tlsmarker1,id.vars = "treatment")
colnames(ciber)[2]<-"cell"
ciber$treatment<-factor(ciber$treatment,levels = c("pre","post"))
table(ciber$treatment)
library(ggpubr)
p=ggboxplot(ciber, x="cell", y="value", color = "treatment",
            ylab="Expression",
            xlab="",
            palette = c("blue","red") )
p=p+rotate_x_text(45)
tiff("response_tumor_post.vs.baseline teffect boxplot.tiff", units="in", width=10, height=5, res=300,compression = "lzw")
p+stat_compare_means(aes(group=treatment),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()

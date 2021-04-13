# install.packages('VennDiagram')
# install.packages('splitstackshape')
library(splitstackshape)
library(VennDiagram)
library("tidyverse")
library("clusterProfiler")
library("org.Mm.eg.db")
#install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")


setwd("C:\\Users\\LY/Documents/RNAseq_diabetes/Compare_Khan2019/")
Sig.hi3h = read.table(file="../GSE117741_Khan2019/data/Results_sig_hi3h-basal.txt", row.names = 1)
View(Sig.hi3h)

Sig.IRMOE = read.table(file="../GSE149662_INSRexpress_Khan2020/data/Results_sig_IRMOE.txt",row.names=1)
view(Sig.IRMOE)

Sig.Cen = read.table(file="../GSE147422_Cen2020/data/Results_sigID_BS200-0_noNA_filter.txt",row.names=1)
View(Sig.Cen)


#=======================
# background gene list
norm.IRMOE = read.table(file="Norm_data_filter_IRMOE.txt",row.names=1)
norm.BS = read.table(file="Norm_data_filtered_BS.txt",row.names=1)
View(norm.IRMOE)
View(norm.BS)

norm.IRMOE$entrez = mapIds(org.Mm.eg.db, keys=rownames(norm.IRMOE), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
dup=norm.IRMOE %>% group_by(entrez) %>% dplyr::filter(n() > 1)
View(dup)
norm.IRMOE <- norm.IRMOE[!duplicated(norm.IRMOE$entrez),] 
norm.IRMOE=norm.IRMOE %>% dplyr::filter(!entrez %in% c(NA)) 

norm.BS$entrez = mapIds(org.Mm.eg.db, keys=rownames(norm.BS), column="ENTREZID", keytype="SYMBOL", multiVals="first")
dup=norm.BS %>% group_by(entrez) %>% dplyr::filter(n() > 1)
View(dup)
norm.BS <- norm.BS[!duplicated(norm.BS$entrez),] 
norm.BS=norm.BS %>% dplyr::filter(!entrez %in% c(NA)) 

venn.allgene <- get.venn.partitions(list(IRMOE.allgene=norm.IRMOE$entrez,BS.allgene=norm.BS$entrez))
head(venn)
df_venn.allgene <- listCol_w(venn.allgene, "..values..", fill = NA)
view(df_venn.all)
venn_long.allgene <- df_venn.allgene %>% 
  tidyr::pivot_longer(cols=..values.._fl_00001:..values.._fl_11291,values_to = "Genes")%>% 
  na.omit()%>% select(!name)
view(venn_long.allgene)

allgene<-venn_long.allgene %>% filter(IRMOE.allgene=="TRUE" & BS.allgene=="TRUE")
dim(allgene)

# compare IRMOE, BS, human-ins-common
human.common = read.table(file="../compare_cor/data/venn_long1_human2inter.txt",row.names=1)
human.common <- human.common %>% filter(degree=="3")
view(human.common)
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("entrezgene_id"), filters = "entrezgene_id", values = human.common$entrez , mart = human, attributesL = c("entrezgene_id"), martL = mouse, uniqueRows=T)
View(genesV2)
genes.df <- as.data.frame(genesV2)
genes.df <- genes.df %>% dplyr::rename(entrez=NCBI.gene..formerly.Entrezgene..ID) %>% dplyr::rename(humancommon_mouse=NCBI.gene..formerly.Entrezgene..ID.1)
genes.df$entrez <- as.character(genes.df$entrez)
human.common_mouse <- merge(human.common,genes.df,by="entrez")
View(human.common_mouse)
human.common_mouse <- human.common_mouse %>% dplyr::distinct(humancommon_mouse,.keep_all=TRUE)
dim(human.common_mouse)


genesV3 = getLDS(attributes = c("entrezgene_id"), filters = "entrezgene_id", values = Sig.Cen$entrez , mart = mouse, attributesL = c("entrezgene_id"), martL = human, uniqueRows=T)
View(genesV3)
genes.df <- as.data.frame(genesV3)
genes.df <- genes.df %>% dplyr::rename(entrez=NCBI.gene..formerly.Entrezgene..ID) %>% dplyr::rename(BShuman=NCBI.gene..formerly.Entrezgene..ID.1)
genes.df$entrez <- as.character(genes.df$entrez)
BS1 <- merge(Sig.Cen,genes.df,by="entrez")
View(BS1)

genesV4 = getLDS(attributes = c("entrezgene_id"), filters = "entrezgene_id", values = Sig.IRMOE$entrez , mart = mouse, attributesL = c("entrezgene_id"), martL = human, uniqueRows=T)
View(genesV4)
genes.df <- as.data.frame(genesV2)
genes.df <- genes.df %>% dplyr::rename(entrez=NCBI.gene..formerly.Entrezgene..ID) %>% dplyr::rename(BShuman=NCBI.gene..formerly.Entrezgene..ID.1)
genes.df$entrez <- as.character(genes.df$entrez)
BS1 <- merge(BS,genes.df,by="entrez")
View(BS1)


library(UpSetR)
listInput1 <- list(IRMOE_DEGs=Sig.IRMOE$entrez,
                   BS_DEGs=Sig.Cen$entrez,
                   Human_insulin=human.common_mouse$humancommon_mouse)
upset(fromList(listInput1),nsets=10, nintersects = NA,mb.ratio = c(0.6, 0.4), 
      text.scale = c(2,2, 2, 1.5, 2.5, 2), point.size = 6, line.size = 2,  order.by = c("degree")) #,
      queries = list(list(query = intersects, params = list("Insulin_SMP", "Insulin_FUSION", "BS_DEGs"), color = "blue", active = T), 
                     list(query = intersects,  params = list("Insulin_SMP", "Insulin_FUSION"), color = "purple", active = T), 
                     list(query = intersects, params = list("Insulin_SMP", "BS_DEGs"),color = "purple", active = T),
                     list(query = intersects, params = list("Insulin_FUSION", "BS_DEGs"), color = "purple", active = T)))


venn <- get.venn.partitions(list(IRMOE_DEGs=Sig.IRMOE$entrez,
                                 BS_DEGs=Sig.Cen$entrez,
                                 Human_insulin=human.common_mouse$humancommon_mouse))


# separate nested lists, tranfer wide to long format
library(splitstackshape)
df_venn <- listCol_w(venn, "..values..", fill = NA)
View(df_venn)
venn_long <- df_venn %>% 
  tidyr::pivot_longer(cols=..values.._fl_0001:..values.._fl_5040,
                      values_to = "entrez")%>% 
  na.omit()%>% dplyr::select(!name)
venn_long$entrez <- as.character(venn_long$entrez)
venn_long$symbol <- mapIds(org.Mm.eg.db, keys=venn_long$entrez, column="SYMBOL", keytype="ENTREZID", multiVals="first")
venn_long$degree <- rowSums(venn_long == "TRUE")
venn_long$description <- mapIds(org.Mm.eg.db, keys=venn_long$entrez, column="GENENAME", keytype="ENTREZID", multiVals="first")
venn_long$description <- gsub(x = venn_long$description, pattern = "\\ ",replacement = "_")
venn_long$description <- gsub(x = venn_long$description, pattern = "\\,",replacement = ".")
venn_long <- venn_long[order(-venn_long$degree),]
View(venn_long)
write.table(venn_long,"data/geneinter_IRMOE-BS-Humancommon.txt", row.names = TRUE)
write.csv(venn_long,"data/geneinter_IRMOE-BS-Humancommon.csv", row.names = TRUE)

#
gene.inter<-vennLong.genes%>%filter(Hyperins=="TRUE" & IRMOE=="TRUE")
summary(gene.inter)
head(gene.inter)
write.table(gene.inter,"gene_inter.txt", sep="\t",row.names = TRUE)
#write.csv(gene.inter,"gene_inter.csv", row.names = FALSE)
gene.inter = read.table(file="gene_inter.txt",header=TRUE)

#


#
# KEGG enrichment
library("clusterProfiler")
library("org.Mm.eg.db")
KEGG.inter <- enrichKEGG(gene         = gene.inter$ENTREZID,
                         organism     = 'mmu',
                         universe      = allgene$Genes,
                         pvalueCutoff=0.05, pAdjustMethod="BH", 
                         qvalueCutoff=0.05)
KEGG.inter <- setReadable(KEGG.inter, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
head(KEGG.inter)
KEGGinter.df=as.data.frame(KEGG.inter)
View(KEGGinter.df)
write.table(KEGGinter.df, sep="\t",file="KEGG_inter.txt", row.names=TRUE)
KEGGinter.df = read.table(file="KEGG_inter.txt",header = TRUE)


#
library(enrichplot)
p <- cnetplot(KEGG.inter, showCategory=8,node_label="all",categorySize="pvalue", colorEdge = TRUE )
#, foldChange=geneList
cowplot::plot_grid(p,ncol=1, rel_widths=1, label_size=5)

pathway <- KEGGinter.df %>% separate_rows(geneID) %>% filter(Description %in% c("mTOR signaling pathway","FoxO signaling pathway"))%>%
  dplyr::distinct(geneID,.keep_all=TRUE)
view(pathway)

venn_long[c(venn_long$degree=="3"),]
GOI <- c(pathway$geneID, "Ptpn3","Sh3rf2","Slc43a1")
GOI <- data.frame(GOI)
view(GOI)

#
head(Sig.Cen)
SigCen.sub=dplyr::select(Sig.Cen, symbol, log2FoldChange,entrez) %>%
  rename(BS_FC=log2FoldChange)
head(SigCen.sub,3)
head(Sig.IRMOE)
SigIRMOE.sub=dplyr::select(Sig.IRMOE, symbol, log2FoldChange,entrez) %>%
  rename(IRMOE_FC=log2FoldChange)
head(SigIRMOE.sub, 3)
GOI.FC <- GOI %>% left_join(SigCen.sub, by =c("GOI"="symbol")) %>% 
  left_join(SigIRMOE.sub, by =c("GOI"="symbol"))
view(GOI.FC)

write.table(GOI.FC,"GOI_IRMOE_BS.txt", row.names = TRUE)
GOI.FC = read.table(file="GOI_IRMOE_BS.txt",header=TRUE)


geneInter.sub=dplyr::select(GOI.FC,GOI,BS_FC,IRMOE_FC)
view(geneInter.sub)
#
library("gplots")
m <- as.matrix(as.data.frame(lapply(geneInter.sub[,-1], as.numeric)))
View(m)
par(oma=c(2,2,2,2))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
hclust.m <- function(x) hclust(x, method="ward.D2")
heatmap.2(m,
          labRow = geneInter.sub$GOI,
          #labCol = NA,
          scale = "none", 
          col=my_palette, 
          trace = "none", 
          density.info = "none",
          cexRow = 1,
          cexCol = 0.8,
          offsetRow = -0.2,
          offsetCol = 0,
          distfun = dist,
          #hclustfun = hclust,
          hclustfun=hclust.m,
          Colv=FALSE,
          dendrogram='row',
          key=TRUE, keysize=0.75, key.title = NA,key.xlab="log2 fold change",cex.lab=5.0, cex.axis=1.0,
          #key.par=list(mar=c(1,1,1,1)),
          lhei=c(0.1,12,1), lwid=c(2,10),lmat = rbind(c(0,3),c(2,1),c(0,4)),  # https://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
          key.par = list(mar=c(3.5,0.5,0,5),cex=0.6) # mar=(bottom,left,top,right)
          #key.xtickfun=function() {
          #  cex <- par("cex")*par("cex.axis")
          #  side <- 1
          #  line <- 0
          #  col <- par("col.axis")
          #  font <- par("font.axis")
            #mtext("low", side=side, at=0, adj=0,
            #      line=line, cex=cex, col=col, font=font)
            #mtext("high", side=side, at=1, adj=1,
            #      line=line, cex=cex, col=col, font=font)
            #return(list(labels=FALSE, tick=FALSE))
          #}
)
dev.off()

#=======================


venn <- get.venn.partitions(list(#hi3h20=Sig.hi3h20$symbol,
                                 hi3h=Sig.hi3h$entrez,
                                 IRMOE=Sig.IRMOE$entrez,
                                 BS=Sig.Cen$entrez))
library(splitstackshape)
df_venn <- listCol_w(venn, "..values..", fill = NA)
View(df_venn)
#write.csv(df_venn,"df_venn.csv", row.names = FALSE)
venn_long <- df_venn %>% 
  tidyr::pivot_longer(cols=..values.._fl_0001:..values.._fl_4890,values_to = "Genes")%>% 
  na.omit()%>% dplyr::select(!name)
head(venn_long)
venn_long$entrez <- as.character(venn_long$Genes)
venn_long$symbol = mapIds(org.Mm.eg.db, keys=venn_long$entrez, column="SYMBOL", keytype="ENTREZID", multiVals="first")

write.table(venn_long,file="data/df_venn_hi3h20-IRMOE-BS_symbol.txt",
            sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)
#venn_long = read.csv(file="df_venn_long.csv",header = TRUE)
venn_long = read.table(file="data/df_venn_hi3h20-hi3h-IRMOE-BS_symbol.txt",row.names=1)
View(venn_long)

##
View(Sig.Cen.sub)
Sig.Cen.sub=dplyr::select(Sig.Cen, entrez, log2FoldChange) %>% 
  dplyr::rename(BS_FC=log2FoldChange)
Sig.Cen.sub$entrez <- as.character(Sig.Cen.sub$entrez) 
View(Sig.Cen.sub)
Sig.hi3h.sub=dplyr::select(Sig.hi3h, entrez, log2FoldChange) %>% 
  dplyr::rename(hi3h_FC=log2FoldChange)
Sig.hi3h.sub$entrez <- as.character(Sig.hi3h.sub$entrez)
View(Sig.hi3h.sub)
Sig.hi3h20.sub=dplyr::select(Sig.hi3h20, entrez, log2FoldChange) %>% 
  dplyr::rename(hi3h20_FC=log2FoldChange)
Sig.hi3h20.sub$entrez <- as.character(Sig.hi3h20.sub)
View(Sig.hi3h20.sub)
Sig.IRMOE.sub=dplyr::select(Sig.IRMOE, entrez, log2FoldChange) %>% 
  dplyr::rename(IRMOE_FC=log2FoldChange)
Sig.IRMOE.sub$entrez <- as.character(Sig.IRMOE.sub$entrez)
View(Sig.IRMOE.sub)
vennLong.FC <- left_join(venn_long,Sig.Cen.sub, by =c("entrez"="entrez"))%>%
  left_join(Sig.IRMOE.sub, by =c("entrez"="entrez"))%>%
  left_join(Sig.hi3h.sub, by =c("entrez"="entrez"))
 
head(vennLong.FC,3)

vennLong.FC$symbol <- mapIds(org.Mm.eg.db, keys=vennLong.FC$entrez, column="SYMBOL", keytype="ENTREZID", multiVals="first")
vennLong.FC$description <- mapIds(org.Mm.eg.db, keys=vennLong.FC$entrez, column="GENENAME", keytype="ENTREZID", multiVals="first")
view(vennLong.FC)

vennLong.FC$description <- gsub(x = vennLong.FC$description, pattern = "\\ ",replacement = "_")
vennLong.FC$description <- gsub(x = vennLong.FC$description, pattern = "\\,",replacement = ".")
write.table(vennLong.FC,file="data/venn_hi3h-IRMOE-BS_list.txt",
            sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)
write.csv(vennLong.FC,file="data/venn_hi3h-IRMOE-BS_list.csv", sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)


# Hi3h and BS intersection
gene.inter.CB<-vennLong.FC %>% filter(BS=="TRUE" & hi3h=="TRUE")
summary(gene.inter.CB)

gene.inter <- gene.inter.IB %>% dplyr::select(symbol,IRMOE_FC,BS_FC)
gene.inter <- gene.inter.CB %>% dplyr::select(symbol,hi3h_FC,BS_FC)
view(gene.inter)
library("gplots")
m <- as.matrix(as.data.frame(lapply(gene.inter[,-1], as.numeric)))
m <- as.matrix(as.data.frame(lapply(gene.inter[,-1], as.numeric)))
View(m)

library("RColorBrewer")
par(oma=c(2,2,2,2))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
hclust.m <- function(x) hclust(x, method="ward.D2")
heatmap.2(m,
          labRow = gene.inter$symbol,
          #labCol = NA,
          scale = "none", 
          col=my_palette, 
          trace = "none", 
          density.info = "none",
          cexRow = 1,
          cexCol = 0.8,
          offsetRow = -0.2,
          offsetCol = 0,
          distfun = dist,
          #hclustfun = hclust,
          hclustfun=hclust.m,
          Colv=FALSE,
          dendrogram='row',
          key=TRUE, keysize=0.75, key.title = NA,key.xlab=NA,cex.lab=5.0, cex.axis=1.0,
          #key.par=list(mar=c(1,1,1,1)),
          lhei=c(0.1,12,1), lwid=c(2,10),lmat = rbind(c(0,3),c(2,1),c(0,4)),  # https://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
          key.par = list(mar=c(1,1,0,5),cex=0.6),
          key.xtickfun=function() {
            cex <- par("cex")*par("cex.axis")
            side <- 1
            line <- 0
            col <- par("col.axis")
            font <- par("font.axis")
            mtext("low", side=side, at=0, adj=0,
                  line=line, cex=cex, col=col, font=font)
            mtext("high", side=side, at=1, adj=1,
                  line=line, cex=cex, col=col, font=font)
            return(list(labels=FALSE, tick=FALSE))
          }
)
dev.off()

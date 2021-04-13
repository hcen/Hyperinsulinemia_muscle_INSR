library(tidyverse)
library("DESeq2")
library(org.Mm.eg.db)

setwd("C:\\Users\\LY/Documents/RNAseq_diabetes/GSE147422_Cen2020")

meta.data = read.csv(file="Input/metadata.csv", row.names = 1)
meta.data$group <- factor(paste0(meta.data$Starvation,"_",meta.data$Insulin))
rownames(meta.data) <- gsub(x = rownames(meta.data),
                        pattern = "\\-",
                        replacement = "_")
write.csv(meta.data, sep="\t",file="Input/metadata_group.csv", 
          row.names=TRUE,col.names=NA,quote=FALSE)
meta.data = read.csv(file="Input/metadata_group.csv", row.names = 1)
head(meta.data)


raw.counts = read.table(file="Input/rawcounts.txt", row.names=1,check.names=FALSE)
head(raw.counts)


### create DESeq matrix
count.data.set = DESeqDataSetFromMatrix(countData=raw.counts, 
                                        colData=meta.data, design= ~ group) 
# Filter low count
nrow(count.data.set)
keep <- rowSums(counts(count.data.set)>5) >=3 # genes counts more than 5 in at least 3 samples
count.filter <- count.data.set[keep,]
nrow(count.filter)

# create DESeq object
count.data.set.object <- DESeq(count.data.set)
count.data.set.object <- DESeq(count.filter)
count.data.set.object
# 'vst' normalization (varianceStabilizingTransformation)
vsd <- vst(count.data.set.object)

### extract normalized counts
norm.data = assay(vsd)
head(norm.data)
write.table(norm.data, sep="\t",file="data/Norm_data_filtered.txt", 
            row.names=TRUE,col.names=NA,quote=FALSE)
### hierarchical clustering analyses and to plot a dendrogram. Evaluate dissimilarities (calculate Euclidean distance) between all eight replicates based on their normalized gene counts.
sampleDists <- dist(t(norm.data),  method = "euclidean")

### Having the distance (dissimilarity) we can finally perform hierarchical cluster analysis using hclust function
clusters=hclust(sampleDists)
plot(clusters)
### plots PCA for the first two principal components
plotPCA(vsd, intgroup=c("group")) +
  theme_bw()
#=======================================================
# plot PC2,3,4
getMethod("plotPCA","DESeqTransform")

plotPCA.pc23 <- function (object, ...) 
{
  .local <- function (object, intgroup = "condition", 
                      ntop = 500, returnData = FALSE) 
  {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                                 drop = FALSE])
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }
    else {
      colData(object)[[intgroup]]
    }
    d <- data.frame(PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group, 
                    intgroup.df, name = colnames(object))
    if (returnData) {
      attr(d, "percentVar") <- percentVar[2:3]
      return(d)
    }
    ggplot(data = d, aes_string(x = "PC2", y = "PC3", 
                                color = "group")) + geom_point(size = 3) + 
      xlab(paste0("PC2: ", round(percentVar[2] * 
                                   100), "% variance")) + ylab(paste0("PC3: ", 
                                                                      round(percentVar[3] * 100), "% variance")) + 
      coord_fixed()
  }
  .local(object, ...)
}
plotPCA.pc23(vsd,intgroup=c("group"))+ # 4x3.5in
  theme_bw()
##
plotPCA.pc34 <- function (object, ...) 
{
  .local <- function (object, intgroup = "condition", 
                      ntop = 500, returnData = FALSE) 
  {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                                 drop = FALSE])
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }
    else {
      colData(object)[[intgroup]]
    }
    d <- data.frame(PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
                    intgroup.df, name = colnames(object))
    if (returnData) {
      attr(d, "percentVar") <- percentVar[3:4]
      return(d)
    }
    ggplot(data = d, aes_string(x = "PC3", y = "PC4", 
                                color = "group")) + geom_point(size = 3) + 
      xlab(paste0("PC3: ", round(percentVar[3] * 
                                   100), "% variance")) + ylab(paste0("PC4: ", 
                                                                      round(percentVar[4] * 100), "% variance")) + 
      coord_fixed()
  }
  .local(object, ...)
}
plotPCA.pc34(vsd,intgroup=c("group"))+ # 4x3.5in
  theme_bw()
#======================================================================


library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix( sampleDists )

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pDist<-pheatmap(sampleDistMatrix,
                clustering_distance_rows = sampleDists,
                clustering_distance_cols = sampleDists,
                col = colors)
save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(pDist, "pDist.png")
dev.off()

### Differential expression

res <- results(count.data.set.object, 
               contrast=c("group","BS_200","BS_0"),
               alpha = 0.05)
res <- results(count.data.set.object, 
               contrast=c("group","AS_200","AS_0"),
               alpha = 0.05)
res <- results(count.data.set.object, 
               contrast=c("group","AS_200","BS_200"),
               alpha = 0.05)
res <- results(count.data.set.object, 
               contrast=c("group","AS_0","BS_0"),
               alpha = 0.05)
summary(res)
out <- capture.output(summary(res))
##
##
cat("result summary", out, file="data/res_summary_BS200-0_filtered.txt", sep="\n", append=TRUE)
cat("result summary", out, file="data/res_summary_AS200-0_filtered.txt", sep="\n", append=TRUE)
cat("result summary", out, file="data/res_summary_AS-BS200_filtered.txt", sep="\n", append=TRUE)

res = na.omit(res)
res.ord = res[order(res$padj),]

##
##
write.table(res.ord, sep="\t",file="data/Results_BS200-0_filtered.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.ord, sep="\t",file="data/Results_AS200-0_filtered.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.ord, sep="\t",file="data/Results_AS-BS200_filtered.txt", row.names=TRUE,col.names=NA,quote=FALSE)

#
res.ord=read.table(file="data/Results_BS200-0_filtered.txt",row.names = 1)
res.ord=read.table(file="data/Results_AS200-0_filtered.txt",row.names = 1)
res.ord=read.table(file="data/Results_AS-BS200_filtered.txt",row.names = 1)

View(res.ord)
res.sig = res.ord[res.ord$padj <= 0.05,]
View(res.sig)

res.sig$symbol =rownames(res.sig)
res.sig$entrez = mapIds(org.Mm.eg.db, keys=rownames(res.sig), column="ENTREZID", keytype="SYMBOL", multiVals="first")
res.sig$entrez2 = mapIds(org.Mm.eg.db, keys=rownames(res.sig), column="ENTREZID", keytype="ALIAS", multiVals="first")
#res.sig$ensembl = mapIds(org.Mm.eg.db, keys=res.sig$entrez, column="ENSEMBL", keytype="ENTREZID", multiVals="first")
#columns(org.Mm.eg.db)
View(res.sig)
res.sig <- as.data.frame(res.sig)

dup=res.sig %>% group_by(symbol) %>% dplyr::filter(n() > 1)
dup=res.sig1 %>% group_by(entrez) %>% dplyr::filter(n() > 1)
dup=res.sig1 %>% group_by(entrez2) %>% dplyr::filter(n() > 1)
View(dup)

res.sig1 <- res.sig %>%  mutate(entrez = coalesce(entrez,entrez2)) 
head(res.sig1)
res.sig1 <- res.sig1[order(res.sig1[,'entrez'],-res.sig1[,'baseMean']),] 
res.sig2 <- res.sig1[!duplicated(res.sig1$entrez),] # some essembl and symbol are the same gene, kept one with more total/basemean counts.
dup1=res.sig2 %>% group_by(entrez) %>% dplyr::filter(n() > 1) 
View(dup1)

res.sig3=res.sig2 %>% dplyr::filter(!entrez %in% c(NA)) 
res.sig3 = res.sig3[order(res.sig3$padj),]
View(res.sig3)

write.table(res.sig3,file="data/Results_sigID_BS200-0_noNA_filter.txt", sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.sig3,file="data/Results_sigID_AS200-0_noNA_filter.txt", sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.sig3,file="data/Results_sigID_AS-BS200_noNA_filter.txt", sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)



res.sig = read.table(file="data/Results_sigID_BS200-0_noNA_filter.txt",row.names=1)
View(res.sig)

#==================================================================
# pathway enrichment analyses
# package for GO enrichment
BiocManager::install("clusterProfiler")
library("clusterProfiler")

##
res.sig.BS = read.table(file="data/Results_sigID_BS200-0_noNA_filter.txt",row.names=1)
res.sig.AS = read.table(file="data/Results_sigID_AS200-0_noNA_filter.txt",row.names=1)
res.sig.ST = read.table(file="data/Results_sigID_AS-BS200_noNA_filter.txt",row.names=1)

##
##
View(res.sig.ST)
#res.sig = read.table(file="data/Results_sigID_AS200-0_noNA.txt", row.names = 1)
#res.sig = read.table(file="data/Results_sigID_AS-BS200_noNA.txt", row.names = 1)
#res.sig = read.table(file="data/Results_sigID_AS-BS0_noNA.txt", row.names = 1)

View(res.sig.BS)
##
##
Sig.up.BS = subset(res.sig.BS, log2FoldChange >0)
Sig.down.BS = subset(res.sig.BS, log2FoldChange <0)
dim(Sig.up.BS)
dim(Sig.down.BS)

Sig.up.AS = subset(res.sig.AS, log2FoldChange >0)
Sig.down.AS = subset(res.sig.AS, log2FoldChange <0)

Sig.up.ST = subset(res.sig.ST, log2FoldChange >0)
Sig.down.ST = subset(res.sig.ST, log2FoldChange <0)

write.table(Sig.up.BS, sep="\t",file="data/Sig_up_BS200-0_noNA_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.table(Sig.down.BS, sep="\t",file="data/Sig_down_BS200-0_noNA_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)

write.table(Sig.up.AS, sep="\t",file="data/Sig_up_AS200-0_noNA_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.table(Sig.down.AS, sep="\t",file="data/Sig_down_AS200-0_noNA_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)

write.table(Sig.up.ST, sep="\t",file="data/Sig_up_AS-BS200_noNA_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.table(Sig.down.ST, sep="\t",file="data/Sig_down_AS-BS200_noNA_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)

View(Sig.up.ST)
##
##
#res.ord.BS = read.table(file="data/Results_BS200-0_filtered.txt", row.names = 1)
#res.ord.AS = read.table(file="data/Results_AS200-0_filtered.txt", row.names = 1)
#res.ord.ST = read.table(file="data/Results_AS-BS200_filtered.txt", row.names = 1)

#View(res.ord.ST)
##
norm <- read.table(file = "data/Norm_data_filtered.txt",row.names=1)
dim(norm)
head(norm)
norm$entrez <- mapIds(org.Mm.eg.db, keys=rownames(norm), column="ENTREZID", keytype="SYMBOL", multiVals="first")
norm$entrez2 = mapIds(org.Mm.eg.db, keys=rownames(norm), column="ENTREZID", keytype="ALIAS", multiVals="first")
norm1 <- norm %>%  mutate(entrez = coalesce(entrez,entrez2)) 
norm1$sum <- rowSums(norm1[,1:20])
norm1 <- norm1[order(norm1[,'entrez'],-norm1[,'sum']),] 
dup1=norm1 %>% group_by(entrez) %>% dplyr::filter(n() > 1) 
head(dup1)
norm2 <- norm1[!duplicated(norm1$entrez),] # some essembl and symbol are the same gene, kept one with more total/basemean counts.
dup1=norm2 %>% group_by(entrez) %>% dplyr::filter(n() > 1) 
head(dup1)
norm3=norm2 %>% dplyr::filter(!entrez %in% c(NA)) 
head(norm3)
all_genes = unique(norm3$entrez)
head(all_genes)
str(all_genes)

##
##
KEGG_up.BS <- enrichKEGG(gene         = Sig.up.BS$entrez,
                         organism     = 'mmu',
                         universe      = all_genes,
                         pvalueCutoff=0.005, pAdjustMethod="BH", 
                         qvalueCutoff=0.05)
KEGG_up.BS <- setReadable(KEGG_up.BS, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
head(KEGG_up.BS)
KEGG_up.BS.df=as.data.frame(KEGG_up.BS)
View(KEGG_up.BS.df)
write.csv(KEGG_up.BS.df, sep="\t",file="data/KEGG_up_BS200-0_noNA_BG_filter.csv", row.names=TRUE,col.names=NA,quote=FALSE)
#KEGG_up.BS.df$Description <- gsub(x = KEGG_up.BS.df$Description, pattern = "\\ ",replacement = "_") 
#KEGG_up.BS.df$Description <- gsub(x = KEGG_up.BS.df$Description, pattern = "\\,",replacement = ".") 
#write.table(KEGG_up.BS.df, sep="\t",file="data/KEGG_up_BS200-0_noNA_BG_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)

KEGG_down.BS <- enrichKEGG(gene         = Sig.down.BS$entrez,
                           organism     = 'mmu',
                           universe      = all_genes,
                           pvalueCutoff=0.005, pAdjustMethod="BH", 
                           qvalueCutoff=0.005)
KEGG_down.BS <- setReadable(KEGG_down.BS, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
KEGG_down.BS.df=as.data.frame(KEGG_down.BS)
View(KEGG_down.BS.df)
write.csv(KEGG_down.BS.df, sep="\t",file="data/KEGG_down_BS200-0_noNA_BG_filter.csv", row.names=TRUE,col.names=NA,quote=FALSE)
#KEGG_down.BS.df$Description <- gsub(x = KEGG_down.BS.df$Description, pattern = "\\ ",replacement = "_") 
#KEGG_down.BS.df$Description <- gsub(x = KEGG_down.BS.df$Description, pattern = "\\,",replacement = "_") 
#write.table(KEGG_down.BS.df, sep="\t",file="data/KEGG_down_BS200-0_noNA_BG_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)


KEGG_up.AS <- enrichKEGG(gene         = Sig.up.AS$entrez,
                         organism     = 'mmu',
                         universe      = all_genes,
                         pvalueCutoff=0.005, pAdjustMethod="BH", 
                         qvalueCutoff=0.005)
KEGG_up.AS <- setReadable(KEGG_up.AS, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
head(KEGG_up.AS)
KEGG_up.AS.df=as.data.frame(KEGG_up.AS)
View(KEGG_up.AS.df)
write.csv(KEGG_up.AS.df, sep="\t",file="data/KEGG_up_AS200-0_noNA_BG_filter.csv", row.names=TRUE,col.names=NA,quote=FALSE)
#KEGG_up.AS.df$Description <- gsub(x = KEGG_up.AS.df$Description, pattern = "\\ ",replacement = "_") 
#KEGG_up.AS.df$Description <- gsub(x = KEGG_up.AS.df$Description, pattern = "\\,",replacement = "_") 
#write.table(KEGG_up.AS.df, sep="\t",file="data/KEGG_up_AS200-0_noNA_BG_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)

KEGG_down.AS <- enrichKEGG(gene         = Sig.down.AS$entrez,
                           organism     = 'mmu',
                           universe      = all_genes,
                           pvalueCutoff=0.005, pAdjustMethod="BH", 
                           qvalueCutoff=0.005)
KEGG_down.AS <- setReadable(KEGG_down.AS, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
KEGG_down.AS.df=as.data.frame(KEGG_down.AS)
View(KEGG_down.AS.df)
write.csv(KEGG_down.AS.df, sep="\t",file="data/KEGG_down_AS200-0_noNA_BG_filter.csv", row.names=TRUE,col.names=NA,quote=FALSE)
#KEGG_down.AS.df$Description <- gsub(x = KEGG_down.AS.df$Description, pattern = "\\ ",replacement = "_") 
#KEGG_down.AS.df$Description <- gsub(x = KEGG_down.AS.df$Description, pattern = "\\,",replacement = "_") 
#write.table(KEGG_down.AS.df, sep="\t",file="data/KEGG_down_AS200-0_noNA_BG_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)

KEGG_up.ST <- enrichKEGG(gene         = Sig.up.ST$entrez,
                         organism     = 'mmu',
                         universe      = all_genes,
                         pvalueCutoff=0.005, pAdjustMethod="BH", 
                         qvalueCutoff=0.005)
KEGG_up.ST <- setReadable(KEGG_up.ST, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
head(KEGG_up.ST)
KEGG_up.ST.df=as.data.frame(KEGG_up.ST)
View(KEGG_up.ST.df)
write.csv(KEGG_up.ST.df, sep="\t",file="data/KEGG_up_AS-BS200_noNA_BG_filter.csv", row.names=TRUE,col.names=NA,quote=FALSE)
#KEGG_up.ST.df$Description <- gsub(x = KEGG_up.ST.df$Description, pattern = "\\ ",replacement = "_") 
#KEGG_up.ST.df$Description <- gsub(x = KEGG_up.ST.df$Description, pattern = "\\,",replacement = "_") 
#write.table(KEGG_up.ST.df, sep="\t",file="data/KEGG_up_AS-BS200_noNA_BG_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)

KEGG_down.ST <- enrichKEGG(gene         = Sig.down.ST$entrez,
                           organism     = 'mmu',
                           universe      = all_genes,
                           pvalueCutoff=0.005, pAdjustMethod="BH", 
                           qvalueCutoff=0.005)
KEGG_down.ST <- setReadable(KEGG_down.ST, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
head(KEGG_down.ST)
KEGG_down.ST.df=as.data.frame(KEGG_down.ST)
View(KEGG_down.ST.df)
write.csv(KEGG_down.ST.df, sep="\t",file="data/KEGG_down_AS-BS200_noNA_BG_filter.csv", row.names=TRUE,col.names=NA,quote=FALSE)
#KEGG_down.ST.df$Description <- gsub(x = KEGG_down.ST.df$Description, pattern = "\\ ",replacement = "_") 
#KEGG_down.ST.df$Description <- gsub(x = KEGG_down.ST.df$Description, pattern = "\\,",replacement = "_")
#write.table(KEGG_down.ST.df, sep="\t",file="data/KEGG_down_AS-BS200_noNA_BG_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)

##
##
BiocManager::install("ReactomePA")
library("ReactomePA")
browseVignettes("ReactomePA")
react_up.BS <- enrichPathway(gene         = Sig.up.BS$entrez,
                             organism     = 'mouse',
                             universe      = all_genes,
                             minGSSize = 10,
                             maxGSSize = 500,
                             pvalueCutoff=0.005, pAdjustMethod="BH", 
                             qvalueCutoff=0.005)
react_up.BS <- setReadable(react_up.BS, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
head(react_up.BS)
react_up.BS.df=as.data.frame(react_up.BS)
View(react_up.BS.df)
write.csv(react_up.BS.df, sep="\t",file="data/react_up_BS200-0_noNA_BG_filter.csv", row.names=TRUE,col.names=NA,quote=FALSE)
#react_up.BS.df$Description <- gsub(x = react_up.BS.df$Description, pattern = "\\ ",replacement = "_") 
#react_up.BS.df$Description <- gsub(x = react_up.BS.df$Description, pattern = "\\,",replacement = "-") 
#write.table(react_up.BS.df, sep="\t",file="data/react_up_BS200-0_noNA_BG_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)

react_down.BS <- enrichPathway(gene         = Sig.down.BS$entrez,
                               organism     = 'mouse',
                               universe      = all_genes,
                               minGSSize = 10,
                               maxGSSize = 500,
                               pvalueCutoff=0.005, pAdjustMethod="BH", 
                               qvalueCutoff=0.005)
react_down.BS <- setReadable(react_down.BS, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
react_down.BS.df=as.data.frame(react_down.BS)
View(react_down.BS.df)
write.csv(react_down.BS.df, sep="\t",file="data/react_down_BS200-0_noNA_BG_filter.csv", row.names=TRUE,col.names=NA,quote=FALSE)
#react_down.BS.df$Description <- gsub(x = react_down.BS.df$Description, pattern = "\\ ",replacement = "_") 
#react_down.BS.df$Description <- gsub(x = react_down.BS.df$Description, pattern = "\\,",replacement = "-") 
#write.table(react_down.BS.df, sep="\t",file="data/react_down_BS200-0_noNA_BG_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)


##
##
react_up.AS <- enrichPathway(gene         = Sig.up.AS$entrez,
                             organism     = 'mouse',
                             universe      = all_genes,
                             pvalueCutoff=0.005, pAdjustMethod="BH", 
                             qvalueCutoff=0.005)
react_up.AS <- setReadable(react_up.AS, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
head(react_up.AS)
react_up.AS.df=as.data.frame(react_up.AS)
View(react_up.AS.df)
write.csv(react_up.AS.df, sep="\t",file="data/react_up_AS200-0_noNA_BG_filter.csv", row.names=TRUE,col.names=NA,quote=FALSE)
#react_up.AS.df$Description <- gsub(x = react_up.AS.df$Description, pattern = "\\ ",replacement = "_") 
#react_up.AS.df$Description <- gsub(x = react_up.AS.df$Description, pattern = "\\,",replacement = "-") 
#write.table(react_up.AS.df, sep="\t",file="data/react_up_AS200-0_noNA_BG_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)

react_down.AS <- enrichPathway(gene         = Sig.down.AS$entrez,
                               organism     = 'mouse',
                               universe      = all_genes,
                               pvalueCutoff=0.005, pAdjustMethod="BH", 
                               qvalueCutoff=0.005)
react_down.AS <- setReadable(react_down.AS, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
react_down.AS.df=as.data.frame(react_down.AS)
View(react_down.AS.df)
write.csv(react_down.AS.df, sep="\t",file="data/react_down_AS200-0_noNA_BG_filter.csv", row.names=TRUE,col.names=NA,quote=FALSE)
#react_down.AS.df$Description <- gsub(x = react_down.AS.df$Description, pattern = "\\ ",replacement = "_") 
#react_down.AS.df$Description <- gsub(x = react_down.AS.df$Description, pattern = "\\,",replacement = "-") 
#write.table(react_down.AS.df, sep="\t",file="data/react_down_AS200-0_noNA_BG_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)

View(react_down.AS.df)
##
##
react_up.ST <- enrichPathway(gene         = Sig.up.ST$entrez,
                             organism     = 'mouse',
                             universe      = all_genes,
                             pvalueCutoff=0.005, pAdjustMethod="BH", 
                             qvalueCutoff=0.005)
react_up.ST <- setReadable(react_up.ST, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
head(react_up.ST)
react_up.ST.df=as.data.frame(react_up.ST)
View(react_up.ST.df)
write.csv(react_up.ST.df, sep="\t",file="data/react_up_AS-BS200_noNA_BG_filter.csv", row.names=TRUE,col.names=NA,quote=FALSE)
#react_up.ST.df$Description <- gsub(x = react_up.ST.df$Description, pattern = "\\ ",replacement = "_") 
#react_up.ST.df$Description <- gsub(x = react_up.ST.df$Description, pattern = "\\,",replacement = "-") 
#write.table(react_up.ST.df, sep="\t",file="data/react_up_AS-BS200_noNA_BG_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)

react_down.ST <- enrichPathway(gene         = Sig.down.ST$entrez,
                               organism     = 'mouse',
                               universe      = all_genes,
                               pvalueCutoff=0.05, pAdjustMethod="BH", 
                               qvalueCutoff=0.1)
react_down.ST <- setReadable(react_down.ST, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
react_down.ST.df=as.data.frame(react_down.ST)
View(react_down.ST.df)
write.csv(react_down.ST.df, sep="\t",file="data/react_down_AS-BS200_noNA_BG_filter.csv", row.names=TRUE,col.names=NA,quote=FALSE)
#react_down.ST.df$Description <- gsub(x = react_down.ST.df$Description, pattern = "\\ ",replacement = "_") 
#react_down.ST.df$Description <- gsub(x = react_down.ST.df$Description, pattern = "\\,",replacement = "-") 
#write.table(react_down.ST.df, sep="\t",file="data/react_down_AS-BS200_noNA_BG_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)

##
##
react.BS <- enrichPathway(gene         = res.sig.BS$entrez,
                          organism     = 'mouse',
                          universe      = all_genes,
                          pvalueCutoff=0.05, pAdjustMethod="BH", 
                          qvalueCutoff=0.1)
react.BS <- setReadable(react.BS, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
head(react.BS)
react.BS.df=as.data.frame(react.BS)
View(react.BS.df)
write.csv(react.BS.df, sep="\t",file="data/react_all_BS200-0_noNA_BG_filter.csv", row.names=TRUE,col.names=NA,quote=FALSE)
#react.BS.df$Description <- gsub(x = react.BS.df$Description, pattern = "\\ ",replacement = "_") 
#react.BS.df$Description <- gsub(x = react.BS.df$Description, pattern = "\\,",replacement = "-") 
#write.table(react.BS.df, sep="\t",file="data/react_all_BS200-0_noNA_BG_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)


##
library(enrichplot)
p1 <- dotplot(react_up.BS, showCategory=30) + #ggtitle("")
  scale_color_gradient(low = "red4", high = "pink", space = "Lab" )
p1
p2 <- dotplot(react_down.BS, showCategory=30) + #ggtitle("")
  scale_color_gradient(low = "blue4", high = "lightblue", space = "Lab" )
p2
plot_grid(p1, p2, ncol=1)
cowplot::plot_grid(p1,p2,ncol=1, rel_widths=1)

##
# 
geneList.BS = res.sig.BS[,2]
## feature 2: named vector
names(geneList.BS) = as.character(res.sig.BS$entrez)
## feature 3: decreasing orde
geneList.BS = sort(geneList.BS, decreasing = TRUE)
head(geneList.BS)

geneList.ST = res.sig.ST[,2]
## feature 2: named vector
names(geneList.ST) = as.character(res.sig.ST$entrez)
## feature 3: decreasing orde
geneList.ST = sort(geneList.ST, decreasing = TRUE)
head(geneList.ST)

view(react_down.BS.df)
react_down.BS$Description
#

p <- cnetplot(react_down.BS, showCategory=c("Signaling by Receptor Tyrosine Kinases"),
              node_label="gene",
              categorySize="pvalue", colorEdge = F,
              order=TRUE, by="Count", foldChange=geneList.BS,
              layout = "kk")+
  scale_color_gradient(#high='#deebf7', low='#2171b5'
                       high="lightblue", low="blue3") +
  labs(edge="c",color = "log2 fold change",size="Gene count") #
p

description <- unlist(as.character(react_down.BS$Description))
View(description)
levels(react_down.BS$Description)
p <- cnetplot(react_down.BS, showCategory=17,
              node_label="category",cex_label_category = 0.1,
              categorySize="pvalue", colorEdge = TRUE,
              order=F, foldChange=geneList,
              layout = "dh")+
  scale_color_gradient(high='#deebf7', low='#2171b5')+
  labs(color = "log2 fold change",size="gene counts")
p
#layout: kk,fr,graphopt,dh,mds,gem,lgl,drl...
##
##
colnames(react_down.BS.df) <- paste(colnames(react_down.BS.df),"BS", sep = "_")
colnames(react_down.AS.df) <- paste(colnames(react_down.AS.df),"AS", sep = "_")
colnames(react_up.ST.df) <- paste(colnames(react_up.ST.df),"ST", sep = "_")
#colnames(KEGG_down.STdf) <- paste(colnames(KEGG_down.STdf),"ST", sep = "_")
View(react_down.BS.df)
View(react_down.AS.df)
View(react_up.ST.df)
# convert fraction to decimal
react_down.BS.df$GeneRatio_bs <- sapply(react_down.BS.df$GeneRatio_BS, function(x) eval(parse(text=x)))
react_down.AS.df$GeneRatio_as <- sapply(react_down.AS.df$GeneRatio_AS, function(x) eval(parse(text=x)))
react_up.ST.df$GeneRatio_st <- sapply(react_up.ST.df$GeneRatio_ST, function(x) eval(parse(text=x)))

react.merge<- left_join(react_down.BS.df, react_up.ST.df, by =c("ID_BS"="ID_ST"))
react.merge<- left_join(react.merge, react_down.AS.df, by =c("ID_BS"="ID_AS"))
View(react.merge)
#write.table(react.merge, sep="\t",file="data/react_merge_BS-ST-AS_BG_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.csv(react.merge, sep="\t",file="data/react_merge_BS-ST-AS_BG_filter.csv", row.names=TRUE,col.names=NA,quote=FALSE)

react.merge1<- full_join(react_down.BS.df, react_up.ST.df, by =c("ID_BS"="ID_ST"))
react.merge1<- full_join(react.merge, react_down.AS.df, by =c("ID_BS"="ID_AS"))
View(react.merge1)
write.table(react.merge1, sep="\t",file="data/react_fulljoin_BS-ST-AS_BG_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.csv(react.merge1, sep="\t",file="data/react_fulljoin_BS-ST-AS_BG_filter.csv", row.names=TRUE,col.names=NA,quote=FALSE)

react <- react.merge %>% mutate(x1="BS") %>% mutate(x2="AS") %>% mutate(x3="ST") %>%
  mutate(Log10adj.P_BS=-log10(p.adjust_BS)*(-1)) %>% 
  mutate(Log10adj.P_AS=-log10(p.adjust_AS)*(-1)) %>%
  mutate(Log10adj.P_ST=-log10(p.adjust_ST))

View(react)
p <- ggplot(react, 
            aes(y = fct_reorder(Description_BS, -Log10adj.P_BS))) + 
  geom_point(aes(size = GeneRatio_bs, color = Log10adj.P_BS, x=x1)) +
  geom_point(aes(size = GeneRatio_as, color = Log10adj.P_AS, x=x2))+
  geom_point(aes(size = GeneRatio_st, color = Log10adj.P_ST, x=x3))+
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_blank())+
  scale_color_gradient2(midpoint = 0, low = "blue4", mid = "white",
                        high = "red4", space = "Lab" )+
  ylab(NULL)+
  xlab(NULL)+
  xlim("BS","ST","AS")+
  theme(axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black",size=10),
        axis.text.x = element_text(colour = "black",size=10),
        legend.title = element_text(color = "black", size = 10))+
  labs(color = "-log10 adj. p",size="Gene ratio") #6.75x3.5
# legend.text = element_blank())+
p
##
##
react_up.ST <- setReadable(react_up.ST, OrgDb = org.Mm.eg.db)
p <- cnetplot(react_up.ST, showCategory=c("Rab regulation of trafficking"),
              node_label="all",
              categorySize="pvalue", colorEdge = F,
              order=TRUE, by="Count", foldChange=geneList.ST,
              layout = "kk")+
  scale_color_gradient(high='red3', low='pink') +
  labs(edge="c",color = "log2 fold change",size="Gene count") #
p
KEGG_down.BS <- setReadable(KEGG_down.BS, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
p <- cnetplot(KEGG_down.BS, showCategory=c("Endocytosis"),
              node_label="all",
              categorySize="pvalue", colorEdge = F,
              order=TRUE, by="Count", foldChange=geneList.BS,
              layout = "kk")+
  scale_color_gradient(high='lightblue', low='blue3') +
  labs(edge="c",color = "log2 fold change",size="Gene count")
p

p <- cnetplot(react_up.BS, showCategory=c("Glucose metabolism","Glycolysis"), #"Regulation of Glucokinase by Glucokinase Regulatory Protein"),
              node_label="all",
              categorySize="pvalue", #colorEdge = T,
              order=TRUE, by="Count", foldChange=geneList.BS,
              layout = "kk")+
  scale_color_gradient(high='red3', low='pink') +
  labs(edge="c",color = "log2 fold change",size="Gene count") #
p

p <- cnetplot(react_down.BS, showCategory=17,
              node_label="category",cex_label_category = 0.1,
              categorySize="pvalue", colorEdge = TRUE,
              order=F, foldChange=geneList,
              layout = "dh")+
  scale_color_gradient(high='#deebf7', low='#2171b5')+
  labs(color = "log2 fold change",size="gene counts")
p
#layout: kk,fr,graphopt,dh,mds,gem,lgl,drl...
##
##

colnames(react_up.BS.df) <- paste(colnames(react_up.BS.df),"BS", sep = "_")
colnames(react_up.AS.df) <- paste(colnames(react_up.AS.df),"AS", sep = "_")
colnames(react_down.ST.df) <- paste(colnames(react_down.ST.df),"ST", sep = "_")
#colnames(KEGG_down.STdf) <- paste(colnames(KEGG_down.STdf),"ST", sep = "_")
head(react_up.BS.df)
View(react_up.AS.df)
View(react_down.ST.df)
#####


# plot pathways
react.up.BS <-  react_up.BS.df %>% mutate(Log10adj.P=log10(p.adjust)*(-1))
react.down.BS <- react_down.BS.df %>% mutate(Log10adj.P=log10(p.adjust)*(-1))
react.up.BS$GeneRatio_num <- sapply(react.up.BS$GeneRatio, function(x) eval(parse(text=x)))
react.down.BS$GeneRatio_num <- sapply(react.down.BS$GeneRatio, function(x) eval(parse(text=x)))

view(react.up.BS[c(1,3,10,20,25,33,35,54,70,98,100,108,118,148,150,162),])

#react.updown.BS <- rbind(react.up.BS[c(1,3,10,20,25,33,35,54,70,81,98,100,108,118,148,150,162),], react.down.BS)
#react$group <- factor(react$group,levels=c("up","down"))

react <- react.up.BS[c(1,3,10,20,25,33,35,54,70,98,100,108,118,148,150,162),]
react <- react.down.BS


view(react)

p <- ggplot(react, 
            aes(y = fct_reorder(Description, -pvalue))) + 
  geom_point(aes(size = Count,
                 x=GeneRatio_num,colour = Log10adj.P)) +
  theme_bw(base_size = 16) +
  scale_color_gradient2(#limits=c(3,60),
    midpoint = 0, mid="pink",
    low = "pink", high = "red4" )+
  ylab(NULL)+
  xlab("Gene Ratio")+
  theme(#axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black",size=14),
        axis.text.x = element_text(colour = "black",size=12),
        legend.title = element_text(color = "black", size = 13))+
  labs(color = "-log10 adj. p",size="Gene counts") +
  theme(panel.grid.minor.y = element_blank())+
  theme(panel.grid.minor.x  = element_blank())
  #facet_grid(~group, scales='free',space = "free")
#theme(panel.grid.major = element_blank()) #10x7.5inch
# legend.text = element_blank())+
p

p <- ggplot(react, 
            aes(y = fct_reorder(Description, -pvalue))) + 
  geom_point(aes(size = Count,
                 x=GeneRatio_num,colour = Log10adj.P)) +
  theme_bw(base_size = 16) +
  scale_color_gradient2(limits=c(2.5,2.7), 
    midpoint = 2.3, 
    low = "lightblue", high = "blue4", space = "Lab" )+
  ylab(NULL)+
  xlab("Gene Ratio")+
  theme(#axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black",size=14),
    axis.text.x = element_text(colour = "black",size=12),
    legend.title = element_text(color = "black", size = 13))+
  labs(color = "-log10 adj. p",size="Gene counts") +
  theme(panel.grid.minor.y = element_blank())+
  theme(panel.grid.minor.x  = element_blank())
p


KEGG.down.BS <-  KEGG_down.BS.df %>% mutate(Log10adj.P=log10(p.adjust)*(-1))
KEGG.down.BS$GeneRatio_num <- sapply(KEGG.down.BS$GeneRatio, function(x) eval(parse(text=x)))

kegg <- KEGG.down.BS


view(kegg)

p <- ggplot(kegg, 
            aes(y = fct_reorder(Description, -pvalue))) + 
  geom_point(aes(size = Count,
                 x=GeneRatio_num,colour = Log10adj.P)) +
  theme_bw(base_size = 16) +
  scale_color_gradient2(#limits=c(3,60),
    midpoint = 2, mid="white",
    low = "lightblue", high = "blue4" )+
  ylab(NULL)+
  xlab("Gene Ratio")+
  theme(#axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black",size=14),
    axis.text.x = element_text(colour = "black",size=12),
    legend.title = element_text(color = "black", size = 13))+
  labs(color = "-log10 adj. p",size="Gene counts") +
  theme(panel.grid.minor.y = element_blank())+
  theme(panel.grid.minor.x  = element_blank())
p

#

#======================================================
# plot heatmaps for genes under certain pathways
library("gplots")
#res.sig.merge <- merge(res.sig.BS, res.sig.AS, by = 7, all=TRUE)

colnames(res.sig.BS) <- paste(colnames(res.sig.BS),"BS", sep = "_")
colnames(res.sig.ST) <- paste(colnames(res.sig.ST),"ST", sep = "_")
colnames(res.sig.AS) <- paste(colnames(res.sig.AS),"AS", sep = "_")
head(res.sig.BS)
res.sig.merge <- res.sig.BS %>% full_join(res.sig.ST, by =c("entrez_BS"="entrez_ST")) %>% full_join(res.sig.AS, by =c("entrez_BS"="entrez_AS"))
View(res.sig.merge)
#
res.sig.merge2 <- res.sig.merge %>% 
  tidyr::pivot_longer(cols=2:3,values_to = "padj")%>% 
  na.omit()

res.sig.merge3 <- res.sig.merge2[order(res.sig.merge2[,'symbol'],res.sig.merge2[,'padj']),] 
View(res.sig.merge3)

res.sig.merge3 <- res.sig.merge3[!duplicated(res.sig.merge3$symbol),] 
dup=res.sig.merge3 %>% group_by(symbol) %>% dplyr::filter(n() > 1) 
View(dup)

res.sig.merge3 <- res.sig.merge3[order(res.sig.merge3[,'padj']),]
res.sig.merge.top <- res.sig.merge3[c(1:50),]
View(res.sig.merge.top)

## 
##
norm <- read.table(file = "data/Norm_data_filtered.txt",row.names=1)
dim(norm)
View(norm)
norm$symbol <- rownames(norm)
norm$entrez <- mapIds(org.Mm.eg.db, keys=rownames(norm), column="ENTREZID", keytype="SYMBOL", multiVals="first")
norm$entrez2 = mapIds(org.Mm.eg.db, keys=rownames(norm), column="ENTREZID", keytype="ALIAS", multiVals="first")
norm1 <- norm %>%  mutate(entrez = coalesce(entrez,entrez2)) 
norm1$sum <- rowSums(norm1[,1:20])
norm1 <- norm1[order(norm1[,'entrez'],-norm1[,'sum']),] 
dup1=norm1 %>% group_by(entrez) %>% dplyr::filter(n() > 1) 
head(dup1)
norm2 <- norm1[!duplicated(norm1$entrez),] # some essembl and symbol are the same gene, kept one with more total/basemean counts.
dup1=norm2 %>% group_by(entrez) %>% dplyr::filter(n() > 1) 
head(dup1)
norm3=norm2 %>% dplyr::filter(!entrez %in% c(NA)) 
head(norm3)
dim(norm3)
#head(res.sig.BS)
#res.sig.BS1 <- res.sig.BS[c(1:50),]
norm4 <- norm3 %>% filter(norm3$symbol %in% res.sig.merge.top$symbol)
View(norm4)

m <- as.matrix(as.data.frame(lapply(norm4[,c(1:20)], as.numeric),check.names=F))
View(m)
library("RColorBrewer")
par(oma=c(2,2,2,2))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
hclust.m <- function(x) hclust(x, method="ward.D2")
heatmap.2(m,
          labRow = norm4$symbol,
          #labCol = NA,
          scale = "row", 
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
          key.par = list(mar=c(3.5,7,0,12),cex=0.6),
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
# col=brewer.pal(11,"RdBu")
?heatmap.2
?hclust

##
##


##
react_down.BS.path <- react_down.BS.sep %>% filter(Description_BS=="")
pathway <- react_down.BS.df %>% separate_rows(geneID_BS) %>% filter(Description_BS=="Signaling by Receptor Tyrosine Kinases")
pathway <- react_down.BS.df %>% separate_rows(geneID_BS) %>% filter(Description_BS=="Signaling by Insulin receptor")
pathway <- KEGG_down.BS.df %>% separate_rows( geneID_BS) %>% filter(Description_BS=="Endocytosis")
pathway <- KEGG_down.BS.df %>% separate_rows( geneID_BS) %>% filter(Description_BS=="FoxO signaling pathway")

head(react_up.ST.df)
pathway1 <- KEGG_down.BS.df %>% separate_rows(geneID) %>% filter(Description=="Endocytosis")
pathway2 <- react_up.ST.df %>% separate_rows(geneID) %>% filter(Description=="Rab regulation of trafficking")
head(pathway1)
pathway <- rbind(pathway1,pathway2) %>% dplyr::distinct(geneID,.keep_all=TRUE)
pathway$symbol <- mapIds(org.Mm.eg.db, keys=pathway$geneID, column="SYMBOL", keytype="ENTREZID", multiVals="first")
view(pathway)
head(react_up.BS.df)
pathway1 <- react_up.BS.df %>% separate_rows(geneID) %>% filter(Description==c("Glucose metabolism","Glycolysis","Gluconeogenesis"))
pathway2 <- react_down.ST.df %>% separate_rows(geneID) %>% filter(Description==c("Glucose metabolism","Glycolysis","Gluconeogenesis"))
pathway3 <- react_up.AS.df %>% separate_rows(geneID) %>% filter(Description==c("Glucose metabolism","Glycolysis","Gluconeogenesis"))
head(pathway1)
head(pathway2)
head(pathway3)
pathway <- rbind(pathway1,pathway2,pathway3) %>% dplyr::distinct(geneID,.keep_all=TRUE)
view(pathway)
##

View(pathway)
# separate_rows(df, V2) #https://stackoverflow.com/questions/15347282/split-delimited-strings-in-a-column-and-insert-as-new-rows
head(norm3)
#norm.path <- norm3 %>% filter(norm3$entrez %in% pathway$geneID)
norm.path <- norm3 %>% filter(norm3$symbol %in% pathway$geneID)
#norm.path1 <- norm3 %>% filter(norm3$symbol %in% pathway$geneID_BS & c("Igf1r","Akt3","Shc2"))
#norm.path <- rbind(norm.path, norm.path1)
View(norm.path)
head(norm.path)
dim(norm.path)
dim(pathway)

m.path <- as.matrix(as.data.frame(lapply(norm.path[,c(1:20)], as.numeric),check.names=F))
library(gplots)
par(oma=c(2,2,2,2))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
hclust.m <- function(x) hclust(x, method="ward.D2")
out <- heatmap.2(m.path,
          labRow = norm.path$symbol,
          #labCol = NA,
          scale = "row", 
          col=my_palette, 
          trace = "none", 
          density.info = "none",
          cexRow = 0.9,
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
          lhei=c(0.1,12,1), lwid=c(5,30),lmat = rbind(c(0,3),c(2,1),c(0,4)),  # https://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
          key.par = list(mar=c(5,3,0,8),cex=0.6),
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
          },
          colsep=c(5, 10, 15) # width 4.5 in
)
dev.off()
out$rowInd
res.sig.path <- res.sig.merge %>% filter(res.sig.merge$entrez_BS %in% norm.path$entrez)
res.sig.path$entrez_BS <- as.character(res.sig.path$entrez_BS)
res.sig.path$symbol <- mapIds(org.Mm.eg.db, keys=res.sig.path$entrez_BS, column="SYMBOL", keytype="ENTREZID", multiVals="first")

res.sig.path <- res.sig.path[,c(26,6,15,23,2,11,19,8)]
View(res.sig.path)
view(norm.path)
#index <- match(norm.path$entrez,res.sig.path$entrez_BS)
index <- match(norm.path$entrez,res.sig.path$entrez_BS)
res.sig.path1 <- res.sig.path[index,]
res.sig.path2 <- res.sig.path1[out$rowInd,] %>% map_df(rev)
res.sig.path2$description <- mapIds(org.Mm.eg.db, keys=res.sig.path2$symbol, column="GENENAME", keytype="SYMBOL", multiVals="first")
res.sig.path2$description <- gsub(x = res.sig.path2$description, pattern = "\\ ",replacement = "_") 
res.sig.path2$description <- gsub(x = res.sig.path2$description, pattern = "\\,",replacement = ".") 
res.sig.path2$description <- gsub(x = res.sig.path2$description, pattern = "\\/",replacement = "..") 
View(res.sig.path2)
#res.sig.path3 <- as.matrix(res.sig.path2)
#view(res.sig.path3)
##
write.csv(res.sig.path2,file="data/heatmap-endocytosis-p.csv",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE) 
write.csv(res.sig.path2,file="data/heatmap-foxo-p.csv",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE) 
write.csv(res.sig.path2,file="data/heatmap-insSignaling-p.csv",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE) 
write.csv(res.sig.path3,file="data/heatmap-endo-rab-p.csv",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE) 
write.csv(res.sig.path2,file="data/heatmap-glucosemetabolism-p.csv",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE) 
##

## TF heatmap

TF <- c("Sin3a","Myc","Ets1","Jund","Max","Mxi1","Maz","Hcfc1","Nrf1","Elf1","Ctcf")
norm.tf <- norm3 %>% filter(norm3$symbol %in% TF)
View(norm.tf)
m.tf <- as.matrix(as.data.frame(lapply(norm.tf[,c(1:20)], as.numeric),check.names=F))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
hclust.m <- function(x) hclust(x, method="ward.D2")
heatmap.2(m.tf,
          labRow = norm.tf$symbol,
          #labCol = NA,
          scale = "row", 
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
          lhei=c(0.1,12,1), lwid=c(2,30),lmat = rbind(c(0,3),c(2,1),c(0,4)),  # https://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
          key.par = list(mar=c(3.5,7,0,12),cex=0.6),
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
          },
          colsep=c(5, 10, 15)
)
dev.off()
head(res.sig.merge3)
res.sig.tf <- res.sig.merge1 %>% filter(res.sig.merge1$symbol %in% TF)
View(res.sig.tf)
##
##

##
library(VennDiagram)

# make a gene logFC list for networkanalyst.ca
res.sig.BS = read.table(file="data/Results_sigID_BS200-0_noNA_filter.txt",row.names=1)
head(res.sig.BS)
res.sig.sub = dplyr::select(res.sig.BS, symbol,log2FoldChange)
View(res.sig.sub)
res.sig.sub1 <- res.sig.sub[c(1:1000),]
res.sig.sub2 <- res.sig.sub[c(1:500),]
write.table(res.sig.sub, sep="\t",file="data/Results_sig_BS200-0_networkanalyst.txt", 
            row.names=FALSE,quote = FALSE)
write.table(res.sig.sub1, sep="\t",file="data/Results_sig_BS200-0_networkanalyst1000.txt", 
            row.names=FALSE,quote = FALSE)
write.table(res.sig.sub1, sep="\t",file="data/Results_sig_BS200-0_networkanalyst500.txt", 
            row.names=FALSE,quote = FALSE)

TF <- read.csv(file="data/node_table_TF_BS1000.csv")
TF <- read.csv(file="data/node_table_TF_BS500.csv")
head(TF)
head(norm3)

grid.newpage()
grid::grid.draw(VennDiagram::venn.diagram(list("TF"=TF$Id,
                                               "counts"=norm3$entrez), NULL)) #6.5x5.5in
venn <- get.venn.partitions(list("TF"=TF$Id,
                                 "counts"=norm3$entrez))
View(venn)
venn[3,4]
as.character(venn[3,4])
mapIds(org.Mm.eg.db, keys=c("16364", "14460", "21349"), column="SYMBOL", keytype="ENTREZID", multiVals="first")

#===============================================================================
Sig.Cen.BS = read.table(file="GSE147422_Cen2020/data/Results_sig_BS200-0.txt", row.names = 1)
View(Sig.Cen.BS)
Sig.Cen.AS = read.table(file="GSE147422_Cen2020/data/Results_sig_AS200-0.txt", row.names = 1)
View(Sig.Cen.AS)
Sig.IRMOE = read.table(file="GSE149662_INSRexpress_Khan2020/data/Results_sig_INSR.txt",row.names=1)
View(Sig.IRMOE)
IRMOE.na=na.omit(Sig.IRMOE$SYMBOL)

res.sig.BS = read.table(file="data/Results_sigID_BS200-0_noNA_filter.txt",row.names=1)
res.sig.AS = read.table(file="data/Results_sigID_AS200-0_noNA_filter.txt",row.names=1)
res.sig.ST = read.table(file="data/Results_sigID_AS-BS200_noNA_filter.txt",row.names=1)
head(res.sig.BS)

#### draw Venn diagram
grid.newpage()
grid::grid.draw(VennDiagram::venn.diagram(list(BS=res.sig.BS$entrez,AS=res.sig.AS$entrez), NULL))#,IRMOE=IRMOE.na
grid::grid.draw(VennDiagram::venn.diagram(list(BS=res.sig.BS$entrez,AS=res.sig.AS$entrez,ST=res.sig.ST$entrez), NULL))#,IRMOE=IRMOE.na

str(venn)
venn <- get.venn.partitions(list(BS=res.sig.BS$entrez,AS=res.sig.AS$entrez))
venn <- get.venn.partitions(list(BS=res.sig.BS$entrez,AS=res.sig.AS$entrez,ST=res.sig.ST$entrez))

#### separate nested lists, tranfer wide to long format
library(splitstackshape)
df_venn <- listCol_w(venn, "..values..", fill = NA)
View(df_venn)
#write.csv(df_venn,"df_venn.csv", row.names = FALSE)
venn_long <- df_venn %>% 
  tidyr::pivot_longer(cols=..values.._fl_0001:..values.._fl_2082,values_to = "entrez")%>% 
  na.omit()%>% dplyr::select(!name)
View(venn_long)
#vignette("pivot")
#?pivot_longer
#write.csv(venn_long,"df_venn_long.csv", row.names = FALSE)
write.table(venn_long,"data/df_venn_BS-AS.txt", row.names = TRUE)
write.table(venn_long,"data/df_venn_BS-ST-AS.txt", row.names = TRUE)
#venn_long = read.table(file="df_venn_BS-AS.txt")
head(venn_long)

gene.inter <- venn_long %>% filter(BS=="TRUE" & AS=="TRUE")
gene.inter <- venn_long %>% filter(BS=="TRUE" & AS=="TRUE" & ST=="TRUE")

head(gene.inter)
head(res.sig.BS)
BS <- res.sig.BS %>% filter(entrez %in% gene.inter$entrez) %>% dplyr::select(symbol,log2FoldChange)
AS <- res.sig.AS %>% filter(entrez %in% gene.inter$entrez) %>% dplyr::select(symbol,log2FoldChange)
ST <- res.sig.ST %>% filter(entrez %in% gene.inter$entrez) %>% dplyr::select(symbol,log2FoldChange)
merge <- left_join(BS, AS, by=c("symbol"="symbol"))
merge <- left_join(merge, ST, by=c("symbol"="symbol"))
head(merge)
dim(merge)

write.table(merge,"data/geneFC_inter_BSAS.txt", sep="\t")
write.table(merge,"data/geneFC_inter_BSASST.txt", sep="\t")

m <- as.matrix(as.data.frame(lapply(merge[,c(2:4)], as.numeric),check.names=F))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
hclust.m <- function(x) hclust(x, method="ward.D2")
heatmap.2(m,
          labRow = NA,
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
          lhei=c(0.1,12,2), lwid=c(2,8),lmat = rbind(c(0,3),c(2,1),c(0,4)),  # https://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
          key.par = list(mar=c(4,1,1,6),cex=0.6),
          #key.xtickfun=function() {
          #  cex <- par("cex")*par("cex.axis")
          #  side <- 1
          #  line <- 0
          #  col <- par("col.axis")
          #  font <- par("font.axis")
          #  mtext("low", side=side, at=0, adj=0,
          #        line=line, cex=cex, col=col, font=font)
          #  mtext("high", side=side, at=1, adj=1,
          #        line=line, cex=cex, col=col, font=font)
          #  return(list(labels=FALSE, tick=FALSE))
          #}
          #colsep=c(5, 10, 15)
)
dev.off()






#==================================
BiocManager::install("pathview")
library("pathview")
head(geneList.BS)
str(geneList.BS)
mmu04068 <- pathview(gene.data  = geneList.BS,
                     pathway.id = "mmu04068",
                     species    = "mmu",
                     bins = list(gene = 21, cpd = 10), 
                     trans.fun = list(gene = NULL, cpd = NULL),
                     low = list(gene = "blue", cpd = "magenta"),
                     high = list(gene = "red", cpd ="yellow"),
                     mid =list(gene = "white", cpd = "white"),
                     na.col = "grey",
                     limit      = list(gene=max(abs(geneList.BS)), cpd=1),
                     cex=0.1,res=600
                     #kegg.native=FALSE,
                     )
foxo <- as.data.frame(mmu04068)
view(foxo)

?pathview()
norm.path <- norm3 %>% filter(norm3$entrez %in% foxo$plot.data.gene.kegg.names)
#norm.path1 <- norm3 %>% filter(norm3$symbol %in% pathway$geneID_BS & c("Igf1r","Akt3","Shc2"))
#norm.path <- rbind(norm.path, norm.path1)
View(norm.path)
head(norm.path)
dim(norm.path)
dim(foxo)

m.path <- as.matrix(as.data.frame(lapply(norm.path[,c(1:20)], as.numeric),check.names=F))
library(gplots)
par(oma=c(2,2,2,2))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
hclust.m <- function(x) hclust(x, method="ward.D2")
out <- heatmap.2(m.path,
                 labRow = norm.path$symbol,
                 #labCol = NA,
                 scale = "row", 
                 col=my_palette, 
                 trace = "none", 
                 density.info = "none",
                 cexRow = 0.9,
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
                 lhei=c(0.1,12,1), lwid=c(5,30),lmat = rbind(c(0,3),c(2,1),c(0,4)),  # https://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
                 key.par = list(mar=c(5,3,0,8),cex=0.6),
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
                 },
                 colsep=c(5, 10, 15) # width 4.5 in
)
dev.off()
out$rowInd
res.sig.path <- res.sig.merge %>% filter(res.sig.merge$entrez_BS %in% norm.path$entrez)
res.sig.path$entrez_BS <- as.character(res.sig.path$entrez_BS)
res.sig.path$symbol <- mapIds(org.Mm.eg.db, keys=res.sig.path$entrez_BS, column="SYMBOL", keytype="ENTREZID", multiVals="first")

res.sig.path <- res.sig.path[,c(26,6,15,23,2,11,19,8)]
View(res.sig.path)
view(norm.path)
#index <- match(norm.path$entrez,res.sig.path$entrez_BS)
index <- match(norm.path$entrez,res.sig.path$entrez_BS)
res.sig.path1 <- res.sig.path[index,]
res.sig.path2 <- res.sig.path1[out$rowInd,] %>% map_df(rev)
res.sig.path2$description <- mapIds(org.Mm.eg.db, keys=res.sig.path2$symbol, column="GENENAME", keytype="SYMBOL", multiVals="first")
res.sig.path2$description <- gsub(x = res.sig.path2$description, pattern = "\\ ",replacement = "_") 
res.sig.path2$description <- gsub(x = res.sig.path2$description, pattern = "\\,",replacement = ".") 
res.sig.path2$description <- gsub(x = res.sig.path2$description, pattern = "\\/",replacement = "..") 
View(res.sig.path2)
#res.sig.path3 <- as.matrix(res.sig.path2)
#view(res.sig.path3)
##
write.csv(res.sig.path2,file="data/heatmap-endocytosis-p.csv",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE) 
write.csv(res.sig.path2,file="data/heatmap-foxo-p.csv",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE) 


view(res.sig.BS)

#========
mmu00010 <- pathview(gene.data  = res.sig.merge$entrez_BS,
                     pathway.id = "mmu00010",
                     species    = "mmu",
                     bins = list(gene = 21, cpd = 10), 
                     trans.fun = list(gene = NULL, cpd = NULL),
                     low = list(gene = "blue", cpd = "magenta"),
                     high = list(gene = "red", cpd ="yellow"),
                     mid =list(gene = "white", cpd = "white"),
                     na.col = "grey",
                     limit      = list(gene=max(abs(geneList.BS)), cpd=1),
                     cex=0.1,res=600
                     #kegg.native=FALSE,
)
glucolysis <- as.data.frame(mmu00010$plot.data.gene)
view(glucolysis)
mmu00010
?pathview()
norm.path <- norm3 %>% filter(norm3$entrez %in% foxo$plot.data.gene.kegg.names)

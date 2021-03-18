library(VennDiagram)
library("tidyverse")
library("org.Hs.eg.db")
require("biomaRt")


setwd("C:\\Users\\LY/Documents/RNAseq_diabetes/compare_cor/")
SMP.ins = read.table(file="input/out_p_insulin_R0.2_SMP.txt", row.names = 1)
head(SMP.ins)
NJ.ins = read.table(file="input/out_p_insulin_R0.3_NJ.txt", row.names = 1)
head(NJ.ins)
dbgap.ins <- read.table(file="input/out_p_insulin_R0.2_dbgap.txt",row.names = 1)
head(dbgap.ins)
BS <- read.table(file="input/Results_sigID_BS200-0_noNA_filter.txt",row.names = 1)
head(BS)
BS$entrez <- as.character(BS$entrez)

## convert mouse gene to human orthologous genes
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("entrezgene_id"), filters = "entrezgene_id", values = BS$entrez , mart = mouse, attributesL = c("entrezgene_id"), martL = human, uniqueRows=T)
View(genesV2)
genes.df <- as.data.frame(genesV2)
genes.df <- genes.df %>% dplyr::rename(entrez=NCBI.gene..formerly.Entrezgene..ID) %>% dplyr::rename(BShuman=NCBI.gene..formerly.Entrezgene..ID.1)
genes.df$entrez <- as.character(genes.df$entrez)
BS1 <- merge(BS,genes.df,by="entrez")
View(BS1)


## get Venn partition
venn1 <- get.venn.partitions(list(Insulin_FUSION=dbgap.ins$rowname,
                                  Insulin_SMP=SMP.ins$entrezgene_id,
                                  BS_DEGs=BS1$BShuman)) 

#### separate nested lists, tranfer wide to long format
library(splitstackshape)

df_venn1 <- listCol_w(venn1, "..values..", fill = NA)
View(df_venn1)
venn_long1 <- df_venn1 %>% 
  tidyr::pivot_longer(cols=..values.._fl_0001:..values.._fl_4706,
                      #cols=..values.._fl_0001:..values.._fl_3899,
                      values_to = "entrez")%>% 
  na.omit()%>% dplyr::select(!name)
venn_long1$entrez <- as.character(venn_long1$entrez)
venn_long1$symbol <- mapIds(org.Hs.eg.db, keys=venn_long1$entrez, column="SYMBOL", keytype="ENTREZID", multiVals="first")
venn_long1$degree <- rowSums(venn_long1 == "TRUE")
venn_long1$description <- mapIds(org.Hs.eg.db, keys=venn_long1$entrez, column="GENENAME", keytype="ENTREZID", multiVals="first")
venn_long1$description <- gsub(x = venn_long1$description, pattern = "\\ ",replacement = "_")
venn_long1$description <- gsub(x = venn_long1$description, pattern = "\\,",replacement = ".")
venn_long1 <- venn_long1[order(-venn_long1$degree),]
View(venn_long1)
write.table(venn_long1,"data/venn_long1.txt", row.names = TRUE)
venn_long1 <- read.table("data/venn_long1.txt", row.names = 1)
#
# Upset plot
install.packages("UpSetR")
library(UpSetR)
listInput1 <- list(Insulin_FUSION=dbgap.ins$rowname,
                  Insulin_SMP=SMP.ins$entrezgene_id,
                  BS_DEGs=BS1$BShuman)

upset(fromList(listInput1),nsets=10, nintersects = NA,mb.ratio = c(0.6, 0.4), 
      text.scale = c(2,2, 2, 1.5, 2.5, 2), point.size = 6, line.size = 2,  order.by = c("degree"),
      queries = list(list(query = intersects, params = list("Insulin_SMP", "Insulin_FUSION", "BS_DEGs"), color = "blue", active = T), 
                     list(query = intersects,  params = list("Insulin_SMP", "Insulin_FUSION"), color = "purple", active = T), 
                     list(query = intersects, params = list("Insulin_SMP", "BS_DEGs"),color = "purple", active = T),
                     list(query = intersects, params = list("Insulin_FUSION", "BS_DEGs"), color = "purple", active = T)))

# add correlation coefficient R
head(SMP.ins)
SMP.ins$logP <- -log10(SMP.ins$Insulin_pM_log.p)
SMP.ins1 <- unique(SMP.ins) %>% dplyr::rename(SMP_R=R) %>% dplyr::rename(SMP_p=logP) %>%
  dplyr::rename(entrez=entrezgene_id) %>% dplyr::select(entrez,SMP_R,SMP_p)
SMP.ins1$entrez <- as.character(SMP.ins1$entrez)
SMP.ins1 <- SMP.ins1[!duplicated(SMP.ins1$entrez),] 
head(SMP.ins1)


head(dbgap.ins)
dbgap.ins$logP <- -log10(dbgap.ins$Insulin_pM_log.p)
dbgap.ins1 <- unique(dbgap.ins) %>% dplyr::rename(dbgap_R=R) %>% dplyr::rename(dbgap_p=logP) %>%
  dplyr::rename(entrez=rowname) %>% dplyr::select(entrez,dbgap_R,dbgap_p)
dbgap.ins1$entrez <- as.character(dbgap.ins1$entrez)
dbgap.ins1 <- dbgap.ins1[!duplicated(dbgap.ins1$entrez),] 
head(dbgap.ins1)

head(BS1)
BS1$logP <- -log10(BS1$padj)
BS2 <- BS1 %>% dplyr::rename(BS_p=logP) %>%
  dplyr::rename(log2FC=log2FoldChange) %>% dplyr::select(BShuman,log2FC,BS_p)
BS2$BShuman <- as.character(BS2$BShuman)
BS2 <- BS2[!duplicated(BS2$BShuman),] 
head(BS2)

##
View(venn_long1)
venn_long1$entrez <- as.character(venn_long1$entrez)
vennLong1 <- venn_long1 %>% left_join(SMP.ins1, by =c("entrez"="entrez"))%>%
  left_join(dbgap.ins1, by =c("entrez"="entrez")) %>% left_join(BS2, by =c("entrez"="BShuman")) 
vennLong1[vennLong1=="Inf"] <- 20
View(vennLong1)
write.csv(vennLong1, sep="\t",file="data/SMP_FUSION_BS_list.csv", row.names=TRUE,col.names=NA,quote=FALSE)
write.table(vennLong1, sep="\t",file="data/SMP_FUSION_BS_list.txt", row.names=TRUE,col.names=NA,quote=FALSE)
#===========================================================================================

## all 
m <- vennLong1[vennLong1$degree>=3, #c(8,11:18)
               c(7,10:15)] %>% mutate(x1="SMP") %>% 
  mutate(x3="FUSION") %>% mutate(x5="BS")  
head(m)
m1 <- m[,c(1,2,4,6)]
head(m1)
view(m1)

# creat heatmap and extract clusting info

library("gplots")
m2 <- as.matrix(as.data.frame(lapply(m1[,-1], as.numeric),check.names=F))
head(m2)

library("RColorBrewer")
par(oma=c(2,2,2,2))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
hclust.m <- function(x) hclust(x, method="ward.D2")
out <-heatmap.2(m2,
          labRow = m1$symbol,
          #labCol = NA,
          scale = "none", 
          col=my_palette, 
          trace = "none", 
          density.info = "none",
          cexRow = 0.8,
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

?hclust
out$rowInd
m3 <- m[out$rowInd,]
m3$order <- rownames(m3)
View(m3)
write.csv(m3, sep="\t",file="data/human-common-IS-ins.csv", row.names=TRUE,col.names=NA,quote=FALSE)
m3 <- read.csv(file="data/human-common-IS-ins.csv",row.names = 1)
head(m3)
install.packages("ggnewscale")
library(ggnewscale) # so fold change can be a different color scale from R

p <- ggplot(m3, aes(y = fct_reorder(symbol,order))) + 
  geom_point(aes(color = SMP_R, size = SMP_p, x=x1)) +
  geom_point(aes(color = dbgap_R, size = dbgap_p, x=x3))+
  #geom_point(aes(color = NJ_R, size = NJ_p, x=x4))+
  #theme(panel.grid.major = element_blank())+
  scale_color_gradient2(limits=c(-0.5,0.5) ,
    midpoint = 0, low = "blue4", mid = "white",
    high = "red4", space = "Lab", na.value = NA)+
  labs(size = "-log10 p",color="R")+
  # geoms below will use another color scale https://github.com/eliocamp/ggnewscale
  new_scale_color() +
  geom_point(aes(color = log2FC, size = BS_p,
                 x=x5), shape=15)+
  scale_color_gradient2(#limits=c(-8,5) ,
    midpoint = 0, low = "blue4", mid = "white",
    high = "red4", space = "Lab" )+
  labs(color="log2FoldChange")+
  theme_bw(base_size = 8) +
  guides(size = guide_legend(override.aes = list(shape = 16) ) )+
  ylab(NULL)+
  xlab(NULL)+
  xlim("SMP","FUSION","BS")+
  #xlim("SMP_female","FUSION_female","BS")+
  theme(axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        legend.title = element_text(color = "black", size = 8)) # 3x12in
p

# phenotypes of human cohorts
SMP.meta = read.table(file="input/metadata_.txt", row.names = 1)
KUP.meta = read.table(file="input/metadata_KUP.txt", row.names = 1)
NJ.meta = read.table(file="input/metadata_df_NJ.txt", row.names = 1)
dbgap.meta = read.table(file="input/metadata_a_dbgap.txt", row.names = 1)

SMP.meta = read.table(file="input/metadata_inslog_SMP.txt", row.names = 1)
NJ.meta = read.table(file="input/metadata_log_NJ.txt", sep="," , row.names = 1,header = T)
dbgap.meta = read.table(file="input/meta_log_dbgap.txt", row.names = 1)

SMP.meta$set <- "SMP"
View(SMP.meta)
NJ.meta$set <- "NJ"
View(NJ.meta)
dbgap.meta$set <- "FUSION"
View(dbgap.meta)
# dbgap.meta <- dbgap.meta %>% dplyr::filter(!subject_ID=="121") # removed outlier determined by hierarchical clustering

SMP.meta <- SMP.meta %>% dplyr::rename(Sex=Gender) %>% dplyr::rename(Age=AGE)
meta.SMP <- SMP.meta %>% dplyr::select(Age,Sex,Glucose_log,Insulin_pM_log,set)
head(meta.SMP)
#meta.SMP$set1 <- paste0(meta.SMP$set,"_",meta.SMP$Sex)
#meta.SMP$set2 <- meta.SMP$set1
#meta.SMP$set3 <- meta.SMP$set


NJ.meta$Age <- "NA"
meta.NJ <- NJ.meta %>% dplyr::select(Age,Sex,Glucose_log,Insulin_pM_log,set)
head(meta.NJ)
#meta.NJ$set1 <- paste0(meta.NJ$set,"_",meta.NJ$Sex)
#meta.NJ$set2 <- meta.NJ$set
#meta.NJ$set3 <- meta.NJ$set2

meta.dbgap <- dbgap.meta %>% dplyr::select(Age,Sex,Glucose_log,Insulin_pM_log,set)
head(meta.dbgap)
#meta.dbgap$set1 <- paste0(meta.dbgap$set,"_",meta.dbgap$Sex)
#meta.dbgap$set2 <- meta.dbgap$set1
#meta.dbgap$set3 <- meta.dbgap$set2

###
count(SMP.meta$Sex=="F")
count(SMP.meta$Sex=="M")

count(NJ.meta$Sex=="F")
count(NJ.meta$Sex=="M")

# calculate SE
#IanStdisNA <- function(x) sd(x)/sqrt(sum(!is.na(x)))
#IanStdisNA(SMP.meta$Insulin_pM)

sd(SMP.meta$Insulin_pM,na.rm=T)
mean(SMP.meta$Insulin_pM,na.rm=T)
sd(SMP.meta$Glucose,na.rm=T)
mean(SMP.meta$Glucose,na.rm=T)
sd(SMP.meta$Age,na.rm=T)
mean(SMP.meta$Age,na.rm=T)

sd(dbgap.meta$Insulin,na.rm=T)
mean(dbgap.meta$Insulin,na.rm=T)
sd(dbgap.meta$Glucose,na.rm=T)
mean(dbgap.meta$Glucose,na.rm=T)
sd(dbgap.meta$Age,na.rm=T)
mean(dbgap.meta$Age,na.rm=T)

sd(NJ.meta$Insulin_nM,na.rm=T)
mean(NJ.meta$Insulin_nM,na.rm=T)
sd(NJ.meta$Glucose,na.rm=T)
mean(NJ.meta$Glucose,na.rm=T)
sd(NJ.meta$Age,na.rm=T)
mean(NJ.meta$Age,na.rm=T)
###

meta <- Reduce(rbind, list(meta.SMP, meta.dbgap, meta.NJ))
View(meta)

#==============================
library(rstatix)

aov <- meta %>% anova_test(Insulin_pM_log ~ set)
aov
stat<- meta %>% pairwise_t_test(Insulin_pM_log ~ set, p.adjust.method = "BH") %>%
  add_significance()
stat
stat<- stat %>% add_xy_position(x = "set")

#stat<- stat %>% mutate(y.position = c(0.7,...)) # %>% add_xy_position(x = "set") ###### unfinished work, to be continued!!!
library(ggpubr)
p1 <- ggplot(meta, aes(x=set, y=Insulin_pM_log)) + 
  geom_boxplot(outlier.shape=NA, outlier.color = "blue")+ 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.4)+
  #stat_summary(fun.data=mean_se, geom="pointrange", color="red",size=0.5)+
  ylab("log10 Insulin (pM)") + xlab(" ")+ 
  scale_y_continuous()+#scale_y_continuous(breaks = seq(0, 3, len = 7))
  xlim("SMP","FUSION","NJ")+
  theme_classic() + theme(legend.position = "none") +
  stat_compare_means(method = "t.test",label = "p.format", hide.ns = T,size=3)+
  #stat_pvalue_manual(stat, label = "p.adj", tip.length = 0, size = 3,hide.ns = T,y.position = c(0.7,0.65))+
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) 
p1

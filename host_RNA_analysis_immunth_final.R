### 07/11/2023 
#### Immunotherapy BAL project 
###Fares Darawshy, Segal Lab


### set wd 
setwd("~/Dropbox (NYU Langone Health)/Fares Darawshy’s files/Home/Projects/Immunotherapy_BAL/final_analysis")

library(immunedeconv)
library(patchwork)
library(phyloseq)
library(vegan)
library(pheatmap)
library(ggplot2)
library(ade4)
library("matrixStats")
library(RColorBrewer)
library(data.table)
library(knitr)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(decontam)
library(DESeq2)
library(microbiomeMarker)
library(stringr)
library(ggprism)
library(forcats)
library(rstatix)
library(patchwork)
library(rmarkdown)
library(ggtext)
library(glue)
library(tidyverse)
library(edgeR)
library(gplots)
library(randomForest)
library(caret)
library(ComplexHeatmap)
library(circlize)
library(pROC)
################################################################################

#read metadata 
metadata <- read.delim2(file = "immuntherapy_metadata_final_cohort.txt")

#make genes names rownames 
names <- make.unique(metadata$X)
rownames(metadata) <- names

#read counts table 
mycounts <- read.delim2(file = "RNA_analysis/RNA_Host_Transcriptome_20230608.txt")

#make genes names rownames 
names <- make.unique(mycounts$X)
rownames(mycounts) <- names

#check they match 
ncol(mycounts)
nrow(metadata)

#remove X column 
mycounts <- mycounts[,-1]
metadata <- metadata[,-1]

#remove X from sample names in mycounts 
colnames(mycounts) <- gsub("X", "", colnames(mycounts))
colnames(mycounts)[4] <- "770341.BAL.RUL.Inv.V1"
rownames(metadata) <- gsub("_", ".", rownames(metadata))

mycounts <- mycounts %>% select(rownames(metadata))


#see what is different 
setdiff(colnames(mycounts), rownames(metadata))

# Arrange mycounts
mycounts <- mycounts [, order(colnames(mycounts))]
metadata <- metadata [order(rownames(metadata)),]


#be sure every thing match 
table(colnames(mycounts)==rownames(metadata))

######## adding variables to metadata ######## 

#overall survival 
metadata$OS_days 
metadata$time_to_death
metadata$death

# 1 year mortality 
metadata <- metadata %>% 
  dplyr::mutate(one_y_mort = case_when(
    death=="1" & OS_days > 365 ~ "greater_12", 
    death=="1" & OS_days < 365 ~ "less_12", 
    death=="0" & OS_days > 365 ~ "greater_12", 
    death=="0" & OS_days < 365 ~ "NA", 
  ))
metadata$one_y_mort

#2 years mortality 
metadata <- metadata %>% 
  dplyr::mutate(two_y_mort = case_when(
    death=="1" & OS_days > 730.5 ~ "greater_24", 
    death=="1" & OS_days < 730.5 ~ "less_24", 
    death=="0" & OS_days > 730.5 ~ "greater_24", 
    death=="0" & OS_days < 730.5 ~ "NA", 
  ))
metadata$two_y_mort

#5 years mortality 
metadata <- metadata %>% 
  dplyr::mutate(five_y_mort = case_when(
    death=="1" & OS_days > 1825 ~ "greater_5_y", 
    death=="1" & OS_days < 1825 ~ "less_5_y", 
    death=="0" & OS_days > 1825 ~ "greater_5_y", 
    death=="0" & OS_days < 1825 ~ "NA", 
  ))

#define histology column 
metadata$histology <- ifelse(metadata$Squamous_cell_carcinoma =="1", "SCC", "adeno")

#define ealry vs late (I-IIIA as early, otherwise late)
metadata$early_vs_late_stage <- ifelse(metadata$stage_at_bronch %in% c("IIIB", "IIIC", "IVA", "IVB"), "late", "early")



############################# Host RNA analysis #############################
#compare only BAL samples and V1 and one sample per patient. 
#should be 100 samples 
metadata_deseq <- metadata %>% filter(Sample_Type=="BAL") %>% filter(prim_sec_sample=="prim")
mycounts_deseq <- mycounts[,rownames(metadata_deseq)]

#create DESeq object

#Convert Count Table into a Numeic Data Frame
d1 = data.frame(lapply(mycounts_deseq, function(x) as.numeric(as.character(x))), check.names=F, row.names = rownames(mycounts_deseq))

#Convert Data to Integers to Run DESEq
d1[] <- lapply(d1, as.integer)

dds <- DESeqDataSetFromMatrix(countData = d1, colData = metadata_deseq, design = ~ ultimate_response_2_lev)

#Normalization Step 
dds <- estimateSizeFactors(dds)

#Retrive normalized counts matrix 
normalized_counts <- counts(dds, normalized=TRUE)
#save it 
write.table(normalized_counts, file="RNA_analysis/Results/normalized_counts.txt", sep="\t", quote=F, col.names=NA)


#Filtering
#filter out genes where there are less than 3 samples with normalized counts greater than or equal to 100.
idx <- rowSums( counts(dds, normalized=TRUE) >= 100 ) >= 3
dds <- dds[idx,]

#Transform Data
vsd <- varianceStabilizingTransformation(dds)


#Drop Levels
dds$ultimate_response_2_lev   <- droplevels(dds$ultimate_response_2_lev)
vsd$ultimate_response_2_lev   <- droplevels(vsd$ultimate_response_2_lev)


#export list of genes for future anallysis 
temp_counts <- data.frame(assay(dds))
colnames(temp_counts) <- gsub("X", "", colnames(temp_counts))
write.csv(temp_counts, file = "RNA_analysis/Results/genes_counts_matrix.csv")
######### progressive disease analysis #######

######################Prepare data to run funcitons later 
#choose table to use
dds.analysis <- dds
vsd.analysis <- vsd

#choose ultimate_response_2_lev variable 
dds.analysis$ultimate_response_2_lev <- dds.analysis$ultimate_response_2_lev
vsd.analysis$ultimate_response_2_lev <- vsd.analysis$ultimate_response_2_lev


#choose variable 
v= "ultimate_response_2_lev"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "non_PD")

#change experiment design if needed
#design(dds.analysis) <- ~ Cavitary

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_deseq

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "red"
col2 <- "blue"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$ultimate_response_2_lev

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
#write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0(".txt")))
 #           , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis/Figures/Beta.Diversity.Bray_", paste0(v, paste0(".pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("PD", "Non PD")),size=10) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis/Results/edgeR.results_", paste0(v, paste0(".csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis/Results/", paste0("GSEA_edgeR_", paste0(v, paste0(".rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis/Figures/edgeR_", paste0(v, paste0(alpha, paste0(".pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()





############one year mortality ######

#choose table to use
dds.analysis <- dds
vsd.analysis <- vsd

#choose one_y_mort variable 
dds.analysis$one_y_mort <- dds.analysis$one_y_mort
vsd.analysis$one_y_mort <- vsd.analysis$one_y_mort

dds.analysis <- dds.analysis[,dds.analysis$one_y_mort!="NA"]
vsd.analysis <- vsd.analysis[,vsd.analysis$one_y_mort!="NA"]


#choose variable 
v= "one_y_mort"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "greater_12")

#change experiment design if needed
design(dds.analysis) <- ~ one_y_mort

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(dds.analysis@colData)

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "red3"
col2 <- "orange3"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$one_y_mort

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
#write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0(".txt")))
#           , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis/Figures/Beta.Diversity.Bray_", paste0(v, paste0(".pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("< 12 Months", "≥ 12 Months")),size=10) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis/Results/edgeR.results_", paste0(v, paste0(".csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis/Results/", paste0("GSEA_edgeR_", paste0(v, paste0(".rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis/Figures/edgeR_", paste0(v, paste0(alpha, paste0(".pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()





############two year mortality ######

#choose table to use
dds.analysis <- dds
vsd.analysis <- vsd

#choose two_y_mort variable 
dds.analysis$two_y_mort <- dds.analysis$two_y_mort
vsd.analysis$two_y_mort <- vsd.analysis$two_y_mort

dds.analysis <- dds.analysis[,dds.analysis$two_y_mort!="NA"]
vsd.analysis <- vsd.analysis[,vsd.analysis$two_y_mort!="NA"]


#choose variable 
v= "two_y_mort"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "greater_24")

#change experiment design if needed
design(dds.analysis) <- ~ two_y_mort

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(dds.analysis@colData)

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "red3"
col2 <- "orange3"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$two_y_mort

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
#write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0(".txt")))
#           , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis/Figures/Beta.Diversity.Bray_", paste0(v, paste0(".pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("< 24 Months", "≥ 24 Months")),size=10) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis/Results/edgeR.results_", paste0(v, paste0(".csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis/Results/", paste0("GSEA_edgeR_", paste0(v, paste0(".rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis/Figures/edgeR_", paste0(v, paste0(alpha, paste0(".pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()





############ IPA for outcomes ######


####### two years otucme selecteed pathway 
IPA_res <- read.csv(file = "RNA_analysis/Results/IPA/ipa_two_y_mort_selected_pathways.csv")
IPA_res$p_value <- 10^(-IPA_res$minus_log_pvalue)
IPA_res$adj_p <- p.adjust(IPA_res$p_value, method = "BH")
#IPA_res <- IPA_res %>% filter(z_score!="#NUM!")


#get top 20 upregulated pathways 
top_20 <- IPA_res %>% filter(z_score>0)
top_20 <- top_20 %>% arrange(desc(z_score))
#top_20 <- top_20 %>% filter(minus_log_10_p_value >= -log10(0.05))

#slice top 25 pathways 
top_20 <- top_20 %>% dplyr::slice(1:25)

#top 20 downregulated 
down_20 <- IPA_res %>% filter(z_score<0)
down_20 <- down_20 %>% arrange(desc(z_score))
#down_20 <- down_20 %>% filter(minus_log_10_p_value >= -log10(0.05))

down_20 <- down_20 %>% dplyr::slice(1:25)

#combine 
ipa_res_for_plot <- rbind(top_20, down_20)
ipa_res_for_plot$z_score <- as.numeric(ipa_res_for_plot$z_score)
ipa_res_for_plot$minus_log_pvalue <- as.numeric(ipa_res_for_plot$minus_log_pvalue)
ipa_res_for_plot$regulation <- ifelse(ipa_res_for_plot$z_score > 0, "upregulated", "downregulated")

plot_dat <- ipa_res_for_plot

# Define colors for each levels of qualitative variables
library(circlize)
max(plot_dat$z_score, na.rm = TRUE)
min(plot_dat$z_score,na.rm = TRUE)

col_fun = colorRamp2(c(-3.055, 0, 2.466), c("blue", "white", "orange"))
col_fun(seq(-3.055, 0, 2.466))

plot_dat <- plot_dat %>% arrange(desc(z_score))

#convert data to matrix 
IPA_res_mat <- plot_dat

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat %>% dplyr::select(z_score)

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)
#IPA_res_mat<- IPA_res_mat[,-c(1:2)]

library(ComplexHeatmap)
######ploting
pdf(file = "RNA_analysis/Figures/IPA_two_y_mort_outcome.pdf", height = 14, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =plot_dat$Pathway, 
                        width = unit(2, "cm"))
dev.off()


#recits score two levels 
IPA_res <- read.csv(file = "RNA_analysis/Results/IPA/ipa_RECIST_selected_pathways.csv")
IPA_res$p_value <- 10^(-IPA_res$minus_log_pvalue)
IPA_res$adj_p <- p.adjust(IPA_res$p_value, method = "BH")
#IPA_res <- IPA_res %>% filter(z_score!="#NUM!")


#get top 20 upregulated pathways 
top_20 <- IPA_res %>% filter(z_score>0)
top_20 <- top_20 %>% arrange(desc(z_score))
#top_20 <- top_20 %>% filter(minus_log_10_p_value >= -log10(0.05))

#slice top 25 pathways 
top_20 <- top_20 %>% dplyr::slice(1:25)

#top 20 downregulated 
down_20 <- IPA_res %>% filter(z_score<0)
down_20 <- down_20 %>% arrange(desc(z_score))
#down_20 <- down_20 %>% filter(minus_log_10_p_value >= -log10(0.05))

down_20 <- down_20 %>% dplyr::slice(1:25)

#combine 
ipa_res_for_plot <- rbind(top_20, down_20)
ipa_res_for_plot$z_score <- as.numeric(ipa_res_for_plot$z_score)
ipa_res_for_plot$minus_log_pvalue <- as.numeric(ipa_res_for_plot$minus_log_pvalue)
ipa_res_for_plot$regulation <- ifelse(ipa_res_for_plot$z_score > 0, "upregulated", "downregulated")

plot_dat <- ipa_res_for_plot

# Define colors for each levels of qualitative variables
library(circlize)
max(plot_dat$z_score, na.rm = TRUE)
min(plot_dat$z_score,na.rm = TRUE)

col_fun = colorRamp2(c(-3.13, 0, 3.656), c("blue", "white", "orange"))
col_fun(seq(-3.13, 0, 3.656))

plot_dat <- plot_dat %>% arrange(desc(z_score))

#convert data to matrix 
IPA_res_mat <- plot_dat

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat %>% dplyr::select(z_score)

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)
#IPA_res_mat<- IPA_res_mat[,-c(1:2)]

library(ComplexHeatmap)
######ploting
pdf(file = "RNA_analysis/Figures/IPA_RECIST_outcome.pdf", height = 20, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold", fontsize=18),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =plot_dat$Pathway, 
                        width = unit(2, "cm"))
dev.off()




########## k means 3 clusters of MT analysis ########

# get K means 3 clusters data 
temp <- read.csv(file = "metatranscriptomics/k_means_3_clusters_MT.csv")

#add clusters data to dds 
dds@colData$k_means_3_cluster <- temp$k_means_3_cluster
vsd@colData$k_means_3_cluster <- temp$k_means_3_cluster

#plot beta diversity for all clusters 

dds.analysis <- dds
vsd.analysis <- vsd

#define variable 
v="k_means_3_cluster"
dds.analysis[[v]] 

#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$k_means_3_cluster
newResults$v <- factor(newResults$v)

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
#write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0(".txt")))
#           , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis/Figures/Beta.Diversity.Bray_", paste0(v, paste0(".pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  #scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Cluster 1", "Cluster 2", "Cluster 3")),size=10) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


#p values for all clusters 
# Define groups
groups <- levels(newResults$v)

# Initialize a results table
pairwise_results <- data.frame(
  Group1 = character(),
  Group2 = character(),
  R2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through all group pairs
for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    group1 <- groups[i]
    group2 <- groups[j]
    
    # Subset data for the two groups
    subset_data <- newResults[newResults$v %in% c(group1, group2), ]
    subset_dist <- vegdist(t(assay(vsd.analysis)[, subset_data$name]), method = "bray")
    
    # Run PERMANOVA
    adonis_result <- adonis2(subset_dist ~ v, data = subset_data)
    
    # Save the results
    pairwise_results <- rbind(pairwise_results, data.frame(
      Group1 = group1,
      Group2 = group2,
      R2 = adonis_result$R2[1],
      p_value = adonis_result$`Pr(>F)`[1]
    ))
  }
}

# Adjust p-values for multiple testing
pairwise_results$p_value_adj <- p.adjust(pairwise_results$p_value, method = "fdr")

# Print the results
print(pairwise_results)

# Save results to file
write.table(pairwise_results, file = "Pairwise_PERMANOVA_Results.txt", sep = "\t", row.names = FALSE)


######### 1 vs 2 analysis 
#choose table to use
dds.analysis <- dds[,dds$k_means_3_cluster!="3"]
vsd.analysis <- vsd[,vsd$k_means_3_cluster!="3"]


#choose variable 
v= "k_means_3_cluster"


#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "2")

vsd.analysis[[v]] <- factor(vsd.analysis[[v]])
vsd.analysis[[v]] <- relevel(vsd.analysis[[v]], ref = "2")


#change experiment design if needed
design(dds.analysis) <- ~ k_means_3_cluster

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(dds.analysis@colData)

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#F8766D"
col2 <- "#00BA38"

#define variable 
dds.analysis[[v]] 


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis/Results/edgeR.results_", paste0(v, paste0("_cluster_1_vs_2.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_cluster_1_vs_2.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_cluster_1_vs_2.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()



###### 1 vs 3 comparison 

#choose table to use
dds.analysis <- dds[,dds$k_means_3_cluster!="2"]
vsd.analysis <- vsd[,vsd$k_means_3_cluster!="2"]


#choose variable 
v= "k_means_3_cluster"


#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "3")

vsd.analysis[[v]] <- factor(vsd.analysis[[v]])
vsd.analysis[[v]] <- relevel(vsd.analysis[[v]], ref = "3")


#change experiment design if needed
design(dds.analysis) <- ~ k_means_3_cluster

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(dds.analysis@colData)

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#F8766D"
col2 <- "#619CFF"

#define variable 
dds.analysis[[v]] 


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis/Results/edgeR.results_", paste0(v, paste0("_cluster_1_vs_3.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_cluster_1_vs_3.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_cluster_1_vs_3.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()




###### 2 vs 3 comparison 

#choose table to use
dds.analysis <- dds[,dds$k_means_3_cluster!="1"]
vsd.analysis <- vsd[,vsd$k_means_3_cluster!="1"]


#choose variable 
v= "k_means_3_cluster"


#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "3")

vsd.analysis[[v]] <- factor(vsd.analysis[[v]])
vsd.analysis[[v]] <- relevel(vsd.analysis[[v]], ref = "3")


#change experiment design if needed
design(dds.analysis) <- ~ k_means_3_cluster

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(dds.analysis@colData)

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#00BA38"
col2 <- "#619CFF"

#define variable 
dds.analysis[[v]] 


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis/Results/edgeR.results_", paste0(v, paste0("_cluster_2_vs_3.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_cluster_2_vs_3.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_cluster_2_vs_3.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()





####### 2 vs 1 analysis 

#choose table to use
dds.analysis <- dds[,dds$k_means_3_cluster!="3"]
vsd.analysis <- vsd[,vsd$k_means_3_cluster!="3"]


#choose variable 
v= "k_means_3_cluster"


#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "1")

vsd.analysis[[v]] <- factor(vsd.analysis[[v]])
vsd.analysis[[v]] <- relevel(vsd.analysis[[v]], ref = "1")


#change experiment design if needed
design(dds.analysis) <- ~ k_means_3_cluster

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(dds.analysis@colData)

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#00BA38"
col2 <- "#F8766D"

#define variable 
dds.analysis[[v]] 


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis/Results/edgeR.results_", paste0(v, paste0("_cluster_2_vs_1.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_cluster_2_vs_1.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_cluster_2_vs_1.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()




#### ipa heatmap for k means pathways 


IPA_res <- read.csv(file = "RNA_analysis/Results/IPA/ipa_MT_k_means_3_clusters_selected_pathways.csv")


# Define colors for each levels of qualitative variables
library(circlize)
max(IPA_res$X1_vs_2, na.rm = TRUE)
min(IPA_res$X1_vs_2,na.rm = TRUE)
max(IPA_res$X1_vs_3, na.rm = TRUE)
min(IPA_res$X1_vs_3,na.rm = TRUE)
max(IPA_res$X2_vs_3, na.rm = TRUE)
min(IPA_res$X2_vs_3,na.rm = TRUE)

col_fun = colorRamp2(c(-3.53, 0, 3.638), c("blue", "white", "orange"))
col_fun(seq(-3.53, 0, 3.638))

IPA_res <- IPA_res %>% arrange(desc(X1_vs_3))

#convert data to matrix 
IPA_res_mat <- IPA_res

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

IPA_res_mat<- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
pdf(file = "RNA_analysis/Figures/IPA_MT_k_means_3_clusters.pdf", height = 16, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))
dev.off()




###### repeat the figure with less pathways and with names in 45 degress 


IPA_res <- read.csv(file = "RNA_analysis/Results/IPA/ipa_MT_k_means_3_clusters_selected_pathways.csv")


# Define colors for each levels of qualitative variables
library(circlize)


IPA_res <- IPA_res[-c(73:77),]
IPA_res <- IPA_res %>% filter(Pathway!="Wound Healing Signaling Pathway")
IPA_res <- IPA_res %>% filter(Pathway!="Pulmonary Healing Signaling Pathway")
IPA_res <- IPA_res %>% filter(Pathway!="Processing of Capped Intronless Pre-mRNA")
IPA_res <- IPA_res %>% filter(Pathway!="STAT3 Pathway")
IPA_res <- IPA_res %>% filter(Pathway!="Transcriptional regulation by RUNX2")
IPA_res <- IPA_res %>% filter(Pathway!="TP53 Regulates Transcription of Cell Death Genes")
IPA_res <- IPA_res %>% filter(Pathway!="Complex IV assembly")
IPA_res <- IPA_res %>% filter(Pathway!="Regulation of Insulin-like Growth Factor (IGF) transport and uptake by IGFBPs")
IPA_res <- IPA_res %>% filter(Pathway!="Signaling by Type 1 Insulin-like Growth Factor 1 Receptor (IGF1R)")
IPA_res <- IPA_res %>% filter(Pathway!="Regulation of Insulin-like Growth Factor (IGF) transport and uptake by IGFBPs")
IPA_res <- IPA_res %>% filter(Pathway!="Signaling by ERBB4")
IPA_res <- IPA_res %>% filter(Pathway!="HMGB1 Signaling")
IPA_res <- IPA_res %>% filter(Pathway!="TP53 Regulates Transcription of Cell Cycle Genes")
IPA_res <- IPA_res %>% filter(Pathway!="G alpha (q) signalling events")
IPA_res <- IPA_res %>% filter(Pathway!="TP53 Regulates Transcription of Cell Cycle Genes")

IPA_res

max(IPA_res$X1_vs_2, na.rm = TRUE)
min(IPA_res$X1_vs_2,na.rm = TRUE)
max(IPA_res$X1_vs_3, na.rm = TRUE)
min(IPA_res$X1_vs_3,na.rm = TRUE)
max(IPA_res$X2_vs_3, na.rm = TRUE)
min(IPA_res$X2_vs_3,na.rm = TRUE)

col_fun = colorRamp2(c(-3.53, 0, 3.638), c("blue", "white", "orange"))
col_fun(seq(-3.53, 0, 3.638))


IPA_res <- IPA_res %>% arrange(desc(X1_vs_3))

#convert data to matrix 
IPA_res_mat <- IPA_res

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

IPA_res_mat<- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

IPA_res_mat_t <- t(IPA_res_mat)

library(ComplexHeatmap)
######ploting
ht <- ComplexHeatmap::Heatmap(IPA_res_mat_t, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",  # Keep row names on the left
                        column_names_side = "top", # Move column names (Pathways) to bottom
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        column_names_gp = gpar(fontface="bold", fontsize = 10),  # Rotate column names
                        column_names_rot = 45,  # Rotate column names 
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "leftcenter-rot"), 
                        width = unit(36, "cm"), # Increase width for readability
                        height = unit(2, "cm") # Adjust height
)

pdf(file = "RNA_analysis/Figures/IPA_MT_k_means_3_clusters_horizontal.pdf", height = 10, width = 20)
# Move legend to the bottom using the draw() function
draw(ht, heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 2), "cm")) # Add padding around plot
dev.off()










 save.image("immunoth_host_RNA_final.RData")











##############comparison according to pneumitis from immunotherapy: yes vs no #################
######################Prepare data to run funcitons later 
#choose table to use
dds.analysis <- dds
vsd.analysis <- vsd

#choose histology variable 
dds.analysis$irAE_pneumonitis_ever_new <- dds.analysis$irAE_pneumonitis_ever_new
vsd.analysis$irAE_pneumonitis_ever_new <- vsd.analysis$irAE_pneumonitis_ever_new


#choose variable 
v= "irAE_pneumonitis_ever_new"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "0")

#change experiment design if needed
design(dds.analysis) <- ~ irAE_pneumonitis_ever_new

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_deseq

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "yellow3"
col2 <- "green3"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$irAE_pneumonitis_ever_new
newResults$v <- factor(newResults$v)

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0(".txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("Figures/Beta.Diversity.Bray_", paste0(v, paste0(".pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Pneumonitis", "Pneumonitis"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0(".csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0(".rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0(".pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()





##############comparison according to mutations: yes vs no #################
######################Prepare data to run funcitons later 
#choose table to use
dds.analysis <- dds[,dds$Driver_mutation_present != "n.a" ]
vsd.analysis <- vsd[,vsd$Driver_mutation_present!="n.a"]

#choose histology variable 
dds.analysis$Driver_mutation_present <- dds.analysis$Driver_mutation_present
vsd.analysis$Driver_mutation_present <- vsd.analysis$Driver_mutation_present


#choose variable 
v= "Driver_mutation_present"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "0")

#change experiment design if needed
design(dds.analysis) <- ~ Driver_mutation_present

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_deseq %>% filter(Driver_mutation_present!="n.a")

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "orangered"
col2 <- "skyblue"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$Driver_mutation_present
newResults$v <- factor(newResults$v)

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0(".txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("Figures/Beta.Diversity.Bray_", paste0(v, paste0(".pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Driver Mutation", "Driver Mutation Present"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0(".csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0(".rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0(".pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()






##############comparison according to stage: late vs early stage #################
######################Prepare data to run funcitons later 
#choose table to use
dds.analysis <- dds
vsd.analysis <- vsd

#choose histology variable 
dds.analysis$early_vs_late_stage <- dds.analysis$early_vs_late_stage
vsd.analysis$early_vs_late_stage <- vsd.analysis$early_vs_late_stage


#choose variable 
v= "early_vs_late_stage"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "early")

#change experiment design if needed
design(dds.analysis) <- ~ early_vs_late_stage

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_deseq %>% filter(early_vs_late_stage!="n.a")

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "orange2"
col2 <- "grey"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$early_vs_late_stage
newResults$v <- factor(newResults$v)

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0(".txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("Figures/Beta.Diversity.Bray_", paste0(v, paste0(".pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Early Stage", "Late Stage"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0(".csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0(".rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0(".pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()




##############comparison according to PD1: low vs negative (1-50% vs < 1%) #################

######################Prepare data to run funcitons later 
#choose table to use
dds.analysis <- dds[,dds$PD_L1_expression %in% c("Negative", "Low")]
vsd.analysis <- vsd[,vsd$PD_L1_expression %in% c("Negative", "Low")]

#choose histology variable 
dds.analysis$PD_L1_expression <- dds.analysis$PD_L1_expression
vsd.analysis$PD_L1_expression <- vsd.analysis$PD_L1_expression


#choose variable 
v= "PD_L1_expression"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "Negative")

#change experiment design if needed
#design(dds.analysis) <- ~ Cavitary

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_deseq %>% filter(PD_L1_expression %in% c("Negative", "Low"))

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "orange3"
col2 <- "grey1"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$PD_L1_expression
newResults$v <- factor(newResults$v, levels = c("Negative", "Low"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0("low_vs_negative.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("Figures/Beta.Diversity.Bray_", paste0(v, paste0("low_vs_negative.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Negative", "Low"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("low_vs_negative.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("low_vs_negative.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("low_vs_negative.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()




##############comparison according to PD1: high vs low (> 50% vs  1-50%) #################

######################Prepare data to run funcitons later 
#choose table to use
dds.analysis <- dds[,dds$PD_L1_expression %in% c("High", "Low")]
vsd.analysis <- vsd[,vsd$PD_L1_expression %in% c("High", "Low")]

#choose histology variable 
dds.analysis$PD_L1_expression <- dds.analysis$PD_L1_expression
vsd.analysis$PD_L1_expression <- vsd.analysis$PD_L1_expression


#choose variable 
v= "PD_L1_expression"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "Low")

#change experiment design if needed
#design(dds.analysis) <- ~ Cavitary

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_deseq %>% filter(PD_L1_expression %in% c("High", "Low"))

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "red3"
col2 <- "orange3"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$PD_L1_expression
newResults$v <- factor(newResults$v, levels = c("Low", "High"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0("_High_vs_Low.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("Figures/Beta.Diversity.Bray_", paste0(v, paste0("_High_vs_Low.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Low", "High"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)
group <- factor(group, levels=c("Low", "High"))

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_High_vs_Low.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_High_vs_Low.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_High_vs_Low.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()














##############comparison according to PD1: high vs negative (> 50% vs  <1%) #################

######################Prepare data to run funcitons later 
#choose table to use
dds.analysis <- dds[,dds$PD_L1_expression %in% c("High", "Negative")]
vsd.analysis <- vsd[,vsd$PD_L1_expression %in% c("High", "Negative")]

#choose histology variable 
dds.analysis$PD_L1_expression <- dds.analysis$PD_L1_expression
vsd.analysis$PD_L1_expression <- vsd.analysis$PD_L1_expression


#choose variable 
v= "PD_L1_expression"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "Negative")

#change experiment design if needed
design(dds.analysis) <- ~ PD_L1_expression

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_deseq %>% filter(PD_L1_expression %in% c("High", "Negative"))

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "red3"
col2 <- "grey1"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$PD_L1_expression
newResults$v <- factor(newResults$v, levels = c("Negative", "High"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0("_High_vs_Negative.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("Figures/Beta.Diversity.Bray_", paste0(v, paste0("_High_vs_Negative.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Negative", "High"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)
group <- factor(group, levels = c("Negative", "High"))
#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_High_vs_Negative.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_High_vs_Negative.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_High_vs_Negative.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()


######################### comparison of PDL1 of 2 levels ########################









##############comparison according to mortality: death vs alive  #################

######################Prepare data to run funcitons later 
#choose table to use
dds.analysis <- dds
vsd.analysis <- vsd

#choose histology variable 
dds.analysis$death <- dds.analysis$death
vsd.analysis$death <- vsd.analysis$death


#choose variable 
v= "death"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "0")

#change experiment design if needed
#design(dds.analysis) <- ~ Cavitary

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_deseq

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "black"
col2 <- "green3"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$death
newResults$v <- factor(newResults$v)

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0("_dead_vs_alive.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("Figures/Beta.Diversity.Bray_", paste0(v, paste0("_dead_vs_alive.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Alive", "Dead"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_dead_vs_alive.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_dead_vs_alive.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_dead_vs_alive.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()




##################### Comparison of DMM clusters: cluster 2 vs 1 ###############

######################Prepare data to run funcitons later 
#choose table to use
dds.analysis <- dds[,dds$dmm_cluster!="3"]
vsd.analysis <- vsd[,vsd$dmm_cluster!="3"]

#choose  variable 
dds.analysis$dmm_cluster <- dds.analysis$dmm_cluster
vsd.analysis$dmm_cluster <- vsd.analysis$dmm_cluster


#choose variable 
v= "dmm_cluster"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "1")

#change experiment design if needed
design(dds.analysis) <- ~ dmm_cluster

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_deseq %>% filter(dmm_cluster!="3")

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "burlywood2"
col2 <- "brown2"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$death
newResults$v <- factor(newResults$v)

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0("_cluster_2_Vs_1.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("Figures/Beta.Diversity.Bray_", paste0(v, paste0("_cluster_2_Vs_1.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Cluster1 ", "Cluster2"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_cluster_2_Vs_1.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_cluster_2_Vs_1.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_cluster_2_Vs_1.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()




##################### Comparison of DMM clusters: cluster 3 vs 1 ###############

######################Prepare data to run funcitons later 
#choose table to use
dds.analysis <- dds[,dds$dmm_cluster!="2"]
vsd.analysis <- vsd[,vsd$dmm_cluster!="2"]

#choose  variable 
dds.analysis$dmm_cluster <- dds.analysis$dmm_cluster
vsd.analysis$dmm_cluster <- vsd.analysis$dmm_cluster


#choose variable 
v= "dmm_cluster"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "1")

#change experiment design if needed
design(dds.analysis) <- ~ dmm_cluster

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_deseq %>% filter(dmm_cluster!="2")

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "chocolate2"
col2 <- "brown2"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$death
newResults$v <- factor(newResults$v)

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0("_cluster_3_Vs_1.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("Figures/Beta.Diversity.Bray_", paste0(v, paste0("_cluster_3_Vs_1.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Cluster1 ", "Cluster3"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_cluster_3_Vs_1.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_cluster_3_Vs_1.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_cluster_3_Vs_1.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()




##################### Comparison of DMM clusters: cluster 3 vs 2 ###############

######################Prepare data to run funcitons later 
#choose table to use
dds.analysis <- dds[,dds$dmm_cluster!="1"]
vsd.analysis <- vsd[,vsd$dmm_cluster!="1"]

#choose  variable 
dds.analysis$dmm_cluster <- dds.analysis$dmm_cluster
vsd.analysis$dmm_cluster <- vsd.analysis$dmm_cluster


#choose variable 
v= "dmm_cluster"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "2")

#change experiment design if needed
design(dds.analysis) <- ~ dmm_cluster

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_deseq %>% filter(dmm_cluster!="1")

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "chocolate2"
col2 <- "burlywood2"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$death
newResults$v <- factor(newResults$v)

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0("_cluster_3_Vs_2.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("Figures/Beta.Diversity.Bray_", paste0(v, paste0("_cluster_3_Vs_2.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Cluster2", "Cluster3"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_cluster_3_Vs_2.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_cluster_3_Vs_2.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_cluster_3_Vs_2.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()









##################### Comparison of DMM clusters: cluster 2 vs 3 ###############

######################Prepare data to run funcitons later 
#choose table to use
dds.analysis <- dds[,dds$dmm_cluster!="1"]
vsd.analysis <- vsd[,vsd$dmm_cluster!="1"]

#choose  variable 
dds.analysis$dmm_cluster <- dds.analysis$dmm_cluster
vsd.analysis$dmm_cluster <- vsd.analysis$dmm_cluster


#choose variable 
v= "dmm_cluster"

#change experiment design if needed
design(dds.analysis) <- ~ dmm_cluster

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "3")

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_deseq %>% filter(dmm_cluster!="1")
metadata_RNA_analysis$dmm_cluster <- factor(metadata_RNA_analysis$dmm_cluster, levels = c("3", "2"))

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "burlywood2"
col2 <- "chocolate2"

#define variable 
dds.analysis[[v]] 

##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_cluster_2_Vs_3.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_cluster_2_Vs_3.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_cluster_2_Vs_3.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()




##################### Comparison of PDL1 only two groups > or < 50% ###############

######################Prepare data to run funcitons later 

#choose  variable 
dds$PD_L1_expression_2_lev <- ifelse(dds$PD_L1_expression_percentage < 50, "low_50", "high_50")
vsd$PD_L1_expression_2_lev <- ifelse(vsd$PD_L1_expression_percentage < 50, "low_50", "high_50")

#choose table to use
dds.analysis <- dds
vsd.analysis <- vsd

#choose variable 
v= "PD_L1_expression_2_lev"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "low_50")

#change experiment design if needed
design(dds.analysis) <- ~ PD_L1_expression_2_lev

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_deseq 
metadata_RNA_analysis$PD_L1_expression_2_lev <- ifelse(metadata_RNA_analysis$PD_L1_expression_percentage < 50, "low_50", "high_50")

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "orange4"
col2 <- "grey50"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$PD_L1_expression_2_lev
newResults$v <- factor(newResults$v, levels = c("low_50", "high_50"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0("_high_vs_low.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("Figures/Beta.Diversity.Bray_", paste0(v, paste0("_high_vs_low.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Cluster2", "Cluster3"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)
group <- factor(group, levels = c("low_50", "high_50"))
#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_high_vs_low.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_high_vs_low.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_high_vs_low.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()









######### one year mortality #####

######################Prepare data to run funcitons later 

#choose  variable 

#choose table to use
dds.analysis <- dds[,dds$one_y_mort!="NA"]
vsd.analysis <- vsd[,vsd$one_y_mort!="NA"]

#choose variable 
v= "one_y_mort"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "greater_12")

#change experiment design if needed
design(dds.analysis) <- ~ one_y_mort

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_deseq %>% filter(one_y_mort!="NA") 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "red3"
col2 <- "orange3"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$one_y_mort
newResults$v <- factor(newResults$v, levels = c("greater_12", "less_12"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0("_less_vs_greater_12.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("Figures/Beta.Diversity.Bray_", paste0(v, paste0("_less_vs_greater_12.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c(">12 Months", "< 12 Months"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_less_vs_greater_12_mon.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_less_vs_greater_12_mon.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_less_vs_greater_12_mon.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()







######### two year mortality #####

######################Prepare data to run funcitons later 

#choose  variable 

#choose table to use
dds.analysis <- dds[,dds$two_y_mort!="NA"]
vsd.analysis <- vsd[,vsd$two_y_mort!="NA"]

#choose variable 
v= "two_y_mort"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "greater_24")

#change experiment design if needed
design(dds.analysis) <- ~ two_y_mort

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_deseq %>% filter(two_y_mort!="NA") 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "gold"
col2 <- "purple"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$two_y_mort
newResults$v <- factor(newResults$v, levels = c("greater_24", "less_24"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0("_less_vs_greater_24.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("Figures/Beta.Diversity.Bray_", paste0(v, paste0("_less_vs_greater_24.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c(">24 Months", "< 24 Months"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_less_vs_greater_24_mon.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_less_vs_greater_24_mon.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_less_vs_greater_24_mon.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()













################################################################################
############################ IPA analysis ##########################

######## IPA dmm clusters #####

#read IPA results 
IPA_res <- read.csv(file = "Results/IPA/IPA_dmm_clusters_comparison_for_heatmap.csv")

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
max(IPA_res$dmm_cluster_3_vs_2, na.rm = TRUE)
min(IPA_res$dmm_cluster_3_vs_2,na.rm = TRUE)
max(IPA_res$dmm_cluster_3_Vs_1,na.rm = TRUE)
min(IPA_res$dmm_cluster_3_Vs_1, na.rm = TRUE)
max(IPA_res$dmm_cluster_2_Vs_1, na.rm = TRUE)
min(IPA_res$dmm_cluster_2_Vs_1, na.rm = TRUE)

col_fun = colorRamp2(c(-4.583, 0, 2.558), c("blue", "white", "orange"))
col_fun(seq(-4.583, 0, 2.558))

#convert data to matrix 
IPA_res_mat <- IPA_res

#sort 
IPA_res_mat <- IPA_res %>% arrange(desc(dmm_cluster_3_vs_2))

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "Figures/IPA_heatmap_dmm_clusters_comparison.pdf", height = 24, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"), height = unit(38, "cm"))
dev.off()










###### IPA dmm v2 with clustering 

#read IPA results 
IPA_res <- read.csv(file = "Results/IPA/DMM_V2.csv")

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
max(IPA_res$dmm_cluster_2_vs_1, na.rm = TRUE)
min(IPA_res$dmm_cluster_2_vs_1,na.rm = TRUE)
max(IPA_res$dmm_cluster_2_vs_3,na.rm = TRUE)
min(IPA_res$dmm_cluster_2_vs_3, na.rm = TRUE)
max(IPA_res$dmm_cluster_3_vs_1, na.rm = TRUE)
min(IPA_res$dmm_cluster_3_vs_1, na.rm = TRUE)

col_fun = colorRamp2(c(-10.553, 0, 3.539), c("blue", "white", "orange"))
col_fun(seq(-10.553, 0, 3.539))

#sort 
IPA_res <- IPA_res %>% arrange(desc(dmm_cluster_2_vs_3))

#convert data to matrix 
IPA_res_mat <- IPA_res

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "Figures/IPA_heatmap_dmm_clusters_comparison_v2.pdf", height = 24, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"), height = unit(38, "cm"))
dev.off()



##### repeat with clustering 
pdf(file = "Figures/IPA_heatmap_dmm_clusters_comparison_v2_with_clustering.pdf", height = 24, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "right",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = TRUE, 
                        cluster_columns = TRUE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"), height = unit(38, "cm"))
dev.off()







######## IPA PD-l1 levels  #####

#read IPA results 
IPA_res <- read.csv(file = "Results/IPA/IPA_PDL_1_comparison_heatmap.csv")

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
max(IPA_res$High_vs_Low, na.rm = TRUE)
min(IPA_res$High_vs_Low,na.rm = TRUE)
max(IPA_res$High_vs_Negative,na.rm = TRUE)
min(IPA_res$High_vs_Negative, na.rm = TRUE)
max(IPA_res$Low_vs_Negative, na.rm = TRUE)
min(IPA_res$Low_vs_Negative, na.rm = TRUE)

col_fun = colorRamp2(c(-2.673, 0, 3.938), c("blue", "white", "orange"))
col_fun(seq(-2.673, 0, 3.938))

#convert data to matrix 
IPA_res_mat <- IPA_res

#sort 
IPA_res_mat <- IPA_res %>% arrange(desc(High_vs_Low))

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "Figures/IPA_PDL_1_comparison_heatmap.pdf", height = 24, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"), height = unit(38, "cm"))
dev.off()


###############################################################################
###### IPA dead vs alive ######

#read results 
IPA_res <- read_csv("Results/IPA/IPA_death_dead_vs_alive_plot.csv")
IPA_res$Pathway <- factor(IPA_res$Pathway, levels = unique(IPA_res$Pathway))
IPA_res <- IPA_res %>% arrange(desc(z_score)) %>% dplyr::filter(z_score != 0)
IPA_res <- IPA_res %>%  mutate(col=ifelse(z_score >0, "orange", "blue"))

pdf(file = "Figures/top_pathways_dead_vs_alive.pdf", width = 34, height = 40)
ggplot(IPA_res, aes(x=z_score, y=fct_reorder(Pathway, z_score, .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in dead vs alive patients")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 24, face = "bold"))
dev.off()


-log(0.15)

###### IPA pneumonitis vs not ######

#read results 
IPA_res <- read.csv(file = "Results/IPA/updated_03_24_24/pneumonitis_ever_new.csv")
IPA_res <- IPA_res %>% arrange(desc(z.score)) %>% dplyr::filter(z.score != 0)
IPA_res <- IPA_res %>% filter(X.log.p.value.>= -log(0.15))
IPA_res <- IPA_res %>%  mutate(col=ifelse(z.score >0, "orange", "blue"))
nrow(IPA_res)
#select top 20 if high nrow 
top_20_up <- IPA_res %>% dplyr::slice(1:20)
top_20_down <- IPA_res %>% arrange(z.score) %>% dplyr::slice(1:20)
comb_res <- bind_rows(top_20_up, top_20_down)

pdf(file = "Figures/IPA_irAE_pneumonitis_vs_not.pdf", width = 22, height = 10)
ggplot(comb_res, aes(x=z.score, y=fct_reorder(Pathway, z.score, .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with pneumonitis vs without")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 16))
dev.off()


###### IPA late vs early stage ######

#read results 
IPA_res <- read.csv(file = "Results/IPA/updated_03_24_24/late_vs_early_stage.csv")
IPA_res <- IPA_res %>% arrange(desc(z.score)) %>% dplyr::filter(z.score != 0)
IPA_res <- IPA_res %>% filter(X.log.p.value.>= -log(0.15))
IPA_res <- IPA_res %>%  mutate(col=ifelse(z.score >0, "orange", "blue"))
nrow(IPA_res)
#select top 20 if high nrow 
top_20_up <- IPA_res %>% dplyr::slice(1:20)
top_20_down <- IPA_res %>% arrange(z.score) %>% dplyr::slice(1:20)
comb_res <- bind_rows(top_20_up, top_20_down)

pdf(file = "Figures/IPA_late_vs_early_stage.pdf", width =  22, height = 10)
ggplot(comb_res, aes(x=z.score, y=fct_reorder(Pathway, z.score, .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with late vs early stage")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 16))
dev.off()



###### IPA PDL1 high vs low ######

#read results 
IPA_res <- read.csv(file = "Results/IPA/updated_03_24_24/PD_L1_high_vs_low.csv")
IPA_res <- IPA_res %>% arrange(desc(z.score)) %>% dplyr::filter(z.score != 0)
IPA_res <- IPA_res %>% filter(X.log.p.value.>= -log(0.15))
IPA_res <- IPA_res %>%  mutate(col=ifelse(z.score >0, "orange", "blue"))
nrow(IPA_res)
#select top 20 if high nrow 
top_20_up <- IPA_res %>% dplyr::slice(1:20)
top_20_down <- IPA_res %>% arrange(z.score) %>% dplyr::slice(1:20)
comb_res <- bind_rows(top_20_up, top_20_down)

pdf(file = "Figures/IPA_high_vs_low_PDL1_50_cutoff.pdf", width = 22, height = 10)
ggplot(comb_res, aes(x=z.score, y=fct_reorder(Pathway, z.score, .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with high vs low PDL1")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 16))
dev.off()




###### IPA one year mortaltity  ######

#read results 
IPA_res <- read.csv(file = "Results/IPA/updated_03_24_24/one_y_mort_less_vs_greater_12.csv")
IPA_res <- IPA_res %>% arrange(desc(z.score)) %>% dplyr::filter(z.score != 0)
IPA_res <- IPA_res %>% filter(X.log.p.value.>= -log(0.15))
IPA_res <- IPA_res %>%  mutate(col=ifelse(z.score >0, "orange", "blue"))
nrow(IPA_res)
#select top 20 if high nrow 
top_20_up <- IPA_res %>% dplyr::slice(1:20)
top_20_down <- IPA_res %>% arrange(z.score) %>% dplyr::slice(1:20)
comb_res <- bind_rows(top_20_up, top_20_down)

pdf(file = "Figures/IPA_one_y_mort_less_vs_greater_12.pdf", width = 22, height = 10)
ggplot(comb_res, aes(x=z.score, y=fct_reorder(Pathway, z.score, .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  #ggtitle("Top Upregulated Pathways in patients with high vs low PDL1")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 16))
dev.off()




###### IPA histology  ######

#read results 
IPA_res <- read.csv(file = "Results/IPA/updated_03_24_24/histology_scc_vs_adeno.csv")
IPA_res <- IPA_res %>% arrange(desc(z.score)) %>% dplyr::filter(z.score != 0)
IPA_res <- IPA_res %>% filter(X.log.p.value.>= -log(0.15))
IPA_res <- IPA_res %>%  mutate(col=ifelse(z.score >0, "orange", "blue"))
nrow(IPA_res)
#select top 20 if high nrow 
top_20_up <- IPA_res %>% dplyr::slice(1:20)
top_20_down <- IPA_res %>% arrange(z.score) %>% dplyr::slice(1:20)
comb_res <- bind_rows(top_20_up, top_20_down)

pdf(file = "Figures/IPA_histology_scc_vs_adeno.pdf", width = 22, height = 12)
ggplot(comb_res, aes(x=z.score, y=fct_reorder(Pathway, z.score, .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  #ggtitle("Top Upregulated Pathways in patients with high vs low PDL1")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 16))
dev.off()



###### IPA driver mutation  ######

#read results 
IPA_res <- read.csv(file = "Results/IPA/updated_03_24_24/driver_mutation.csv")
IPA_res <- IPA_res %>% arrange(desc(z.score)) %>% dplyr::filter(z.score != 0)
IPA_res <- IPA_res %>% filter(X.log.p.value.>= -log(0.15))
IPA_res <- IPA_res %>%  mutate(col=ifelse(z.score >0, "orange", "blue"))
nrow(IPA_res)
#select top 20 if high nrow 
top_20_up <- IPA_res %>% dplyr::slice(1:20)
top_20_down <- IPA_res %>% arrange(z.score) %>% dplyr::slice(1:20)
comb_res <- bind_rows(top_20_up, top_20_down)

pdf(file = "Figures/IPA_driver_mutation.pdf", width = 22, height = 10)
ggplot(comb_res, aes(x=z.score, y=fct_reorder(Pathway, z.score, .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  #ggtitle("Top Upregulated Pathways in patients with high vs low PDL1")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 16))
dev.off()






#plot mRNA levels of il-17 

IL17_plot_data <- mycounts %>% t(.) %>%data.frame() %>% select(starts_with("IL17"))

#add to this metadata 
IL17_plot_data <- merge(IL17_plot_data, metadata, by="row.names", all=TRUE)
IL17_plot_data <- IL17_plot_data %>% mutate(dmm_cluster=factor(dmm_cluster))

#plot for all il17 measures 

col_names <- colnames(IL17_plot_data)[2:12]

for (i in col_names){
 p<-  ggplot(IL17_plot_data, aes_string(x = IL17_plot_data$dmm_cluster, y = i)) +
    geom_boxplot(width = 0.2, color = "blue") +
    geom_jitter(width = 0.2,alpha = 0.5,color = "darkblue") +
    scale_y_log10()+
   stat_compare_means(comparisons = list(c("1", "2"), 
                                         c("2","3"), 
                                         c("1", "3")))+
   ggtitle(i)+
    theme_bw()
  
 #save it 
 pdf_output <- paste0(i, paste0("_IL_17_boxplot", paste0("_DMM_clusters",paste0(".pdf"))))
 pdf(file = pdf_output, width = 4, height = 6)
 show(p)
 dev.off()
}


####### plot mRNA of IL-8 (CXCL8)#####

IL8_plot_data <- mycounts %>% t(.) %>%data.frame() %>% select(starts_with("CXCL8"))

#add to this metadata 
IL8_plot_data <- merge(IL8_plot_data, metadata, by="row.names", all=TRUE)
IL8_plot_data <- IL8_plot_data %>% mutate(dmm_cluster=factor(dmm_cluster))

#plot for all IL8_plot_data measures 
p<-  ggplot(IL8_plot_data, aes(x = IL8_plot_data$dmm_cluster, y = CXCL8)) +
    geom_boxplot(width = 0.2, color = "blue") +
    geom_jitter(width = 0.2,alpha = 0.5,color = "darkblue") +
    scale_y_log10()+
    stat_compare_means(comparisons = list(c("1", "2"), 
                                          c("2","3"), 
                                          c("1", "3")))+
  theme_bw()
  
#save it 
pdf_output <- paste0("comparison", paste0("_IL_8_boxplot", paste0("_DMM_clusters",paste0(".pdf"))))
pdf(file = pdf_output, width = 4, height = 6)
p
dev.off()



############## IPA SCC vs adeno FDR 0.2 ##########

IPA_res <- read_csv("Results/IPA/IPA_SCC_vs_adeno_FDR_0.2_plot.csv")
IPA_res$Pathway <- factor(IPA_res$Pathway, levels = unique(IPA_res$Pathway))
IPA_res <- IPA_res %>% arrange(desc(z_score)) %>% dplyr::filter(z_score != 0) %>% filter(z_score!="#NUM!")
IPA_res$neg_log_p_value <- IPA_res$p_value

#leave only pathways with p value < 0.2 (FDR)
IPA_res <- IPA_res %>% filter(neg_log_p_value > -log(0.2))

#set colors of pthways. oragnge is upregulatedm, blue is downregulated 
IPA_res <- IPA_res %>%  mutate(col=ifelse(z_score >0, "orange", "blue"))

pdf(file = "Figures/IPA_SCC_vs_adeno_FDR_0.2_plot.pdf", width = 20, height = 12)
ggplot(IPA_res, aes(x=as.numeric(z_score), y=fct_reorder(Pathway,as.numeric(z_score), .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with SCC vs LUAD")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 24, face = "bold"))
dev.off()



############## IPA mutation present or not FDR 0.2 ##########

IPA_res <- read_csv("Results/IPA/IPA_Driver_mutations_present_FDR_0.2_plot.csv")
IPA_res$Pathway <- factor(IPA_res$Pathway, levels = unique(IPA_res$Pathway))
IPA_res <- IPA_res %>% arrange(desc(z_score)) %>% dplyr::filter(z_score != 0) %>% filter(z_score!="#NUM!")
IPA_res$neg_log_p_value <- IPA_res$p_value

#leave only pathways with p value < 0.2 (FDR)
IPA_res <- IPA_res %>% filter(neg_log_p_value > -log(0.2))

#set colors of pthways. oragnge is upregulatedm, blue is downregulated 
IPA_res <- IPA_res %>%  mutate(col=ifelse(z_score >0, "orange", "blue"))

pdf(file = "Figures/IPA_Driver_mutations_present_FDR_0.2_plot.pdf", width = 20, height = 8)
ggplot(IPA_res, aes(x=as.numeric(z_score), y=fct_reorder(Pathway,as.numeric(z_score), .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with Driver mutation vs not")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 24, face = "bold"))
dev.off()




############## IPA mutation late vs early FDR 0.2 ##########

IPA_res <- read_csv("Results/IPA/IPA_late_vs_early_FDR_0.2.csv")
IPA_res$Pathway <- factor(IPA_res$Pathway, levels = unique(IPA_res$Pathway))
IPA_res <- IPA_res %>% arrange(desc(z_score)) %>% dplyr::filter(z_score != 0) %>% filter(z_score!="#NUM!")
IPA_res$neg_log_p_value <- IPA_res$p_value

#leave only pathways with p value < 0.2 (FDR)
IPA_res <- IPA_res %>% filter(neg_log_p_value > -log(0.2))

#set colors of pthways. oragnge is upregulatedm, blue is downregulated 
IPA_res <- IPA_res %>%  mutate(col=ifelse(z_score >0, "orange", "blue"))

pdf(file = "Figures/IPA_late_vs_early_FDR_0.2.pdf", width = 20, height = 8)
ggplot(IPA_res, aes(x=as.numeric(z_score), y=fct_reorder(Pathway,as.numeric(z_score), .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with Late vs Early Stage NSCLC")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 24, face = "bold"))
dev.off()


############## IPA mutation dead vs alive 0.2 ##########

IPA_res <- read_csv("Results/IPA/IPA_dead_vs_alive_FDR_0.2_plot.csv")
IPA_res$Pathway <- factor(IPA_res$Pathway, levels = unique(IPA_res$Pathway))
IPA_res <- IPA_res %>% arrange(desc(z_score)) %>% dplyr::filter(z_score != 0) %>% filter(z_score!="#NUM!")
IPA_res$neg_log_p_value <- IPA_res$p_value

#leave only pathways with p value < 0.2 (FDR)
IPA_res <- IPA_res %>% filter(neg_log_p_value > -log(0.2))

#set colors of pthways. oragnge is upregulatedm, blue is downregulated 
IPA_res <- IPA_res %>%  mutate(col=ifelse(z_score >0, "orange", "blue"))

pdf(file = "Figures/IPA_dead_vs_alive_FDR_0.2_plot.pdf", width = 20, height = 12)
ggplot(IPA_res, aes(x=as.numeric(z_score), y=fct_reorder(Pathway,as.numeric(z_score), .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with Dead vs Aalive Patients")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 24, face = "bold"))
dev.off()



############## IPA mutation pneumonitis vs not 0.2 ##########

IPA_res <- read_csv("Results/IPA/IPA_irAE_pneumonitis_FDR_0.2_plot.csv")
IPA_res$Pathway <- factor(IPA_res$Pathway, levels = unique(IPA_res$Pathway))
IPA_res <- IPA_res %>% arrange(desc(z_score)) %>% dplyr::filter(z_score != 0) %>% filter(z_score!="#NUM!")
IPA_res$neg_log_p_value <- IPA_res$p_value

#leave only pathways with p value < 0.2 (FDR)
IPA_res <- IPA_res %>% filter(neg_log_p_value > -log(0.2))

#set colors of pthways. oragnge is upregulatedm, blue is downregulated 
IPA_res <- IPA_res %>%  mutate(col=ifelse(z_score >0, "orange", "blue"))

pdf(file = "Figures/IPA_irAE_pneumonitis_FDR_0.2_plot.pdf", width = 20, height = 22)
ggplot(IPA_res, aes(x=as.numeric(z_score), y=fct_reorder(Pathway,as.numeric(z_score), .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with Pneumonitis vs Without Pneumonitis")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 24, face = "bold"))
dev.off()







#######pd l1 levels - 2 levels IPA ###

IPA_res <- read_csv("Results/IPA/IPA_PD_L1_expression_2_lev_FDR_0.2_plot.csv")
IPA_res$Pathway <- factor(IPA_res$Pathway, levels = unique(IPA_res$Pathway))
IPA_res <- IPA_res %>% arrange(desc(z_score)) %>% dplyr::filter(z_score != 0) %>% filter(z_score!="#NUM!")
IPA_res$neg_log_p_value <- IPA_res$p_value

#leave only pathways with p value < 0.2 (FDR)
IPA_res <- IPA_res %>% filter(neg_log_p_value > -log(0.2))

#set colors of pthways. oragnge is upregulatedm, blue is downregulated 
IPA_res <- IPA_res %>%  mutate(col=ifelse(z_score >0, "orange", "blue"))

pdf(file = "Figures/IPA_PD_L1_expression_2_lev_FDR_0.2_plot.pdf", width = 20, height = 8)
ggplot(IPA_res, aes(x=as.numeric(z_score), y=fct_reorder(Pathway,as.numeric(z_score), .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with High vs Low PD-L1")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 24, face = "bold"))
dev.off()

  


######## IPA heatmaps with FDR 0.2 ##### 

# DMM 


#read IPA results 
IPA_res <- read.csv(file = "Results/IPA/comparison_DMM_FDR_0.2_plot.csv")

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
max(IPA_res$Cluster.2.vs.3, na.rm = TRUE)
min(IPA_res$Cluster.2.vs.3,na.rm = TRUE)
max(IPA_res$Cluster.2.vs.1,na.rm = TRUE)
min(IPA_res$Cluster.2.vs.1, na.rm = TRUE)
max(IPA_res$Cluster.3.vs.1, na.rm = TRUE)
min(IPA_res$Cluster.3.vs.1, na.rm = TRUE)

col_fun = colorRamp2(c(-7.147, 0, 2.082), c("blue", "white", "orange"))
col_fun(seq(-7.147, 0, 2.082))
#sort
IPA_res <- IPA_res %>% arrange(desc(Cluster.2.vs.3))

#convert data to matrix 
IPA_res_mat <- IPA_res



#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "Figures/IPA_heatmap_dmm_clusters_comparison_fdr_0.2.pdf", height = 28, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))
dev.off()



#PDL-1 three levels 
#read IPA results 
IPA_res <- read.csv(file = "Results/IPA/IPA_PD_L1_expression_comparison_FDR_0.2_plot.csv")

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
max(IPA_res$Low_vs_Negative, na.rm = TRUE)
min(IPA_res$Low_vs_Negative,na.rm = TRUE)
max(IPA_res$High_vs_Low,na.rm = TRUE)
min(IPA_res$High_vs_Low, na.rm = TRUE)
max(IPA_res$High_vs_Negative, na.rm = TRUE)
min(IPA_res$High_vs_Negative, na.rm = TRUE)

col_fun = colorRamp2(c(-3.138, 0, 3.207), c("blue", "white", "orange"))
col_fun(seq(-3.138, 0, 3.207))

IPA_res <- IPA_res %>% arrange(desc(High_vs_Low))


#convert data to matrix 
IPA_res_mat <- IPA_res


#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "Figures/IPA_PD_L1_expression_comparison_FDR_0.2_plot.pdf", height = 20, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))
dev.off()














############################# Analysis for James, heatmap groups ######
meta_heatmap_manual <- read.csv("metadata_heatmap_manual.csv")
meta_heatmap_manual <- meta_heatmap_manual %>% select(grouping, Subject.ID)

meta_heatmap_manual <- meta_heatmap_manual %>%
  arrange(Subject.ID)


dds.analysis <- dds[,dds$Subject.ID %in% meta_heatmap_manual$Subject.ID]
vsd.analysis <- vsd[,vsd$Subject.ID %in% meta_heatmap_manual$Subject.ID]

temp_dat <- data.frame(dds.analysis@colData)

temp_dat <- inner_join(temp_dat, meta_heatmap_manual, by="Subject.ID")

dds.analysis$grouping <- temp_dat$grouping
vsd.analysis$grouping <- temp_dat$grouping


#choose variable 
v= "grouping"

#### beta diversity for all 

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$grouping
newResults$v <- factor(newResults$v)

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0("goruping_all.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("Figures/Beta.Diversity.Bray_", paste0(v, paste0("all.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("blue", "goldenrod", "grey")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Group 1", "Group 2", "Group 3"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()




######edgeR analysis  2 vs 1 


#choose table to use
dds.analysis_2 <- dds.analysis[,dds.analysis$grouping!="3"]
vsd.analysis_2 <- vsd.analysis[,vsd.analysis$grouping!="3"]

#choose variable 
v= "grouping"

#Set Reference Level for Comparison (Control Group)
dds.analysis_2[[v]] <- factor(dds.analysis_2[[v]])

#change experiment design if needed
design(dds.analysis_2) <- ~ grouping

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(colData(dds.analysis_2)) 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "goldenrod"
col2 <- "blue"

######################Plot PCoA


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis_2)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_2_vs_1.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_2_vs_1.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_2_vs_1.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()







########## repeat edgeR 3 vs 1 



#choose table to use
dds.analysis_2 <- dds.analysis[,dds.analysis$grouping!="2"]
vsd.analysis_2 <- vsd.analysis[,vsd.analysis$grouping!="2"]

#choose variable 
v= "grouping"

#Set Reference Level for Comparison (Control Group)
dds.analysis_2[[v]] <- factor(dds.analysis_2[[v]])

#change experiment design if needed
design(dds.analysis_2) <- ~ grouping

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(colData(dds.analysis_2)) 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "grey"
col2 <- "blue"

######################Plot PCoA


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis_2)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_3_vs_1.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_3_vs_1.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_3_vs_1.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()







########## repeat edgeR 3 vs 2



#choose table to use
dds.analysis_2 <- dds.analysis[,dds.analysis$grouping!="1"]
vsd.analysis_2 <- vsd.analysis[,vsd.analysis$grouping!="1"]

#choose variable 
v= "grouping"

#Set Reference Level for Comparison (Control Group)
dds.analysis_2[[v]] <- factor(dds.analysis_2[[v]])

#change experiment design if needed
design(dds.analysis_2) <- ~ grouping

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(colData(dds.analysis_2)) 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "grey"
col2 <- "goldenrod"

######################Plot PCoA


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis_2)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_3_vs_2.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_3_vs_2.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_3_vs_2.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()





###### IPA results 
IPA_res <- read.csv("Results/IPA/comparison_heatmap_grouping.csv")


# Define colors for each levels of qualitative variables
library(circlize)
max(IPA_res$Group_2_vs_1, na.rm = TRUE)
min(IPA_res$Group_2_vs_1,na.rm = TRUE)
max(IPA_res$Group_3_vs_1,na.rm = TRUE)
min(IPA_res$Group_3_vs_1, na.rm = TRUE)
max(IPA_res$Group_3_vs_2, na.rm = TRUE)
min(IPA_res$Group_3_vs_2, na.rm = TRUE)

col_fun = colorRamp2(c(-3.357, 0, 4.496), c("blue", "white", "red"))
col_fun(seq(-3.357, 0, 4.496))

#convert data to matrix 
IPA_res_mat <- IPA_res
IPA_res_mat <- IPA_res_mat %>% arrange(desc(Group_3_vs_2))

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "Figures/IPA_heatmap_grouping_james_unclustered.pdf", height = 18, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))
dev.off()



pdf(file = "Figures/IPA_heatmap_grouping_james_clustered.pdf", height = 18, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = TRUE, 
                        row_dend_side ="right", column_dend_side = "bottom",
                        cluster_columns = TRUE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))
dev.off()





############################# k means clusters analysis #############

##### metagenomics 

#read MG clusters 
mg_data <- read.csv(file = "MG_data_for_k_means.csv")
#select only ID to merge with metadata 
mg_data_clust <- mg_data %>% dplyr::select(ID, k_means_3_cluster)
#fix ID names to match those of dds metadata 
mg_data_clust$ID <- gsub("_", "." ,mg_data_clust$ID)
#make sure names are the same 
table(colnames(dds.analysis) == mg_data_clust$ID)
#add k means clusters to metadata 
dds.analysis$k_means_clust_MG <- mg_data_clust$k_means_3_cluster
vsd.analysis$k_means_clust_MG <- mg_data_clust$k_means_3_cluster

#choose variable 
v= "k_means_clust_MG"

#### beta diversity for all 

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$k_means_clust_MG
newResults$v <- factor(newResults$v)

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0("_MG.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("Figures/Beta.Diversity.Bray_", paste0(v, paste0("_all_MG.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  #scale_color_manual(values=c("blue", "goldenrod", "grey")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Cluster1", "Cluster2", "Cluster3"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


#compare between the clusters 


######edgeR analysis  2 vs 1 


#choose table to use
dds.analysis_2 <- dds.analysis[,dds.analysis$k_means_clust_MG!="3"]
vsd.analysis_2 <- vsd.analysis[,vsd.analysis$k_means_clust_MG!="3"]

#choose variable 
v= "k_means_clust_MG"

#Set Reference Level for Comparison (Control Group)
dds.analysis_2[[v]] <- factor(dds.analysis_2[[v]])

#change experiment design if needed
design(dds.analysis_2) <- ~ k_means_clust_MG

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(colData(dds.analysis_2)) 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#00ba38"
col2 <- "#f8776d"


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis_2)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_2_vs_1_MG.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_2_vs_1_MG.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_2_vs_1_MG.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()






########## repeat edgeR 3 vs 1 

#choose table to use
dds.analysis_2 <- dds.analysis[,dds.analysis$k_means_clust_MG!="2"]
vsd.analysis_2 <- vsd.analysis[,vsd.analysis$k_means_clust_MG!="2"]

#choose variable 
v= "k_means_clust_MG"

#Set Reference Level for Comparison (Control Group)
dds.analysis_2[[v]] <- factor(dds.analysis_2[[v]])

#change experiment design if needed
design(dds.analysis_2) <- ~ k_means_clust_MG

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(colData(dds.analysis_2)) 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#619CFF"
col2 <- "#f8776d"


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis_2)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_3_vs_1_MG.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_3_vs_1_MG.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_3_vs_1_MG.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()







########## repeat edgeR 2vs 3



#choose table to use
dds.analysis_2 <- dds.analysis[,dds.analysis$k_means_clust_MG!="1"]
vsd.analysis_2 <- vsd.analysis[,vsd.analysis$k_means_clust_MG!="1"]

#choose variable 
v= "k_means_clust_MG"

#Set Reference Level for Comparison (Control Group)
dds.analysis_2[[v]] <- factor(dds.analysis_2[[v]])

#change experiment design if needed
design(dds.analysis_2) <- ~ k_means_clust_MG

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(colData(dds.analysis_2)) 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#00ba38"
col2 <- "#619CFF"


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis_2)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

group <- factor(group, levels = c("3", "2"))
#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_2_vs_3_MG.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_2_vs_3_MG.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_2_vs_3_MG.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()





######### 3 vs 2 


#choose table to use
dds.analysis_2 <- dds.analysis[,dds.analysis$k_means_clust_MG!="1"]
vsd.analysis_2 <- vsd.analysis[,vsd.analysis$k_means_clust_MG!="1"]

#choose variable 
v= "k_means_clust_MG"

#Set Reference Level for Comparison (Control Group)
dds.analysis_2[[v]] <- factor(dds.analysis_2[[v]])

#change experiment design if needed
design(dds.analysis_2) <- ~ k_means_clust_MG

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(colData(dds.analysis_2)) 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#619CFF"
col2 <- "#00ba38"


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis_2)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#group <- factor(group, levels = c("3", "2"))
#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_3_vs_2_MG.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_3_vs_2_MG.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_3_vs_2_MG.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()





############ IPA res










#########MTS k means clust



#read MTS clusters 
MTS_data <- read.csv(file = "MTS_data_for_k_means.csv")
#select only ID to merge with metadata 
MTS_data_clust <- MTS_data %>% dplyr::select(ID, k_means_3_cluster)
#fix ID names to match those of dds metadata 
MTS_data_clust$ID <- gsub("_", "." ,MTS_data_clust$ID)
col_names <- colnames(dds.analysis) 
rownames(MTS_data_clust) <- MTS_data_clust$ID
MTS_data_clust <- MTS_data_clust[col_names,]
#make sure names are the same 
table(colnames(dds.analysis) == MTS_data_clust$ID)
#add k means clusters to metadata 
dds.analysis$k_means_clust_MTS <- MTS_data_clust$k_means_3_cluster
vsd.analysis$k_means_clust_MTS <- MTS_data_clust$k_means_3_cluster

#choose variable 
v= "k_means_clust_MTS"

#### beta diversity for all 

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$k_means_clust_MTS
newResults$v <- factor(newResults$v)

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/Beta.Diversity.Bray.", paste0(v, paste0("_MTS.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("Figures/Beta.Diversity.Bray_", paste0(v, paste0("all_MTS.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  #scale_color_manual(values=c("blue", "goldenrod", "grey")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Cluster1", "Cluster2", "Cluster3"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


#compare between the clusters 


######edgeR analysis  2 vs 1 


#choose table to use
dds.analysis_2 <- dds.analysis[,dds.analysis$k_means_clust_MTS!="3"]
vsd.analysis_2 <- vsd.analysis[,vsd.analysis$k_means_clust_MTS!="3"]

#choose variable 
v= "k_means_clust_MTS"

#Set Reference Level for Comparison (Control Group)
dds.analysis_2[[v]] <- factor(dds.analysis_2[[v]])

#change experiment design if needed
design(dds.analysis_2) <- ~ k_means_clust_MTS

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(colData(dds.analysis_2)) 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#00ba38"
col2 <- "#f8776d"


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis_2)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_2_vs_1_MTS.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_2_vs_1_MTS.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_2_vs_1_MTS.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()






########## repeat edgeR 3 vs 1 

#choose table to use
dds.analysis_2 <- dds.analysis[,dds.analysis$k_means_clust_MTS!="2"]
vsd.analysis_2 <- vsd.analysis[,vsd.analysis$k_means_clust_MTS!="2"]

#choose variable 
v= "k_means_clust_MTS"

#Set Reference Level for Comparison (Control Group)
dds.analysis_2[[v]] <- factor(dds.analysis_2[[v]])

#change experiment design if needed
design(dds.analysis_2) <- ~ k_means_clust_MTS

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(colData(dds.analysis_2)) 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#619CFF"
col2 <- "#f8776d"


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis_2)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_3_vs_1_MTS.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_3_vs_1_MTS.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_3_vs_1_MTS.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()







########## repeat edgeR 2vs 3



#choose table to use
dds.analysis_2 <- dds.analysis[,dds.analysis$k_means_clust_MTS!="1"]
vsd.analysis_2 <- vsd.analysis[,vsd.analysis$k_means_clust_MTS!="1"]

#choose variable 
v= "k_means_clust_MTS"

#Set Reference Level for Comparison (Control Group)
dds.analysis_2[[v]] <- factor(dds.analysis_2[[v]])

#change experiment design if needed
design(dds.analysis_2) <- ~ k_means_clust_MTS

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(colData(dds.analysis_2)) 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#00ba38"
col2 <- "#619CFF"


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis_2)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

group <- factor(group, levels = c("3", "2"))
#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_2_vs_3_MTS.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_2_vs_3_MTS.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_2_vs_3_MTS.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()






######### 3 vs 2 


#choose table to use
dds.analysis_2 <- dds.analysis[,dds.analysis$k_means_clust_MTS!="1"]
vsd.analysis_2 <- vsd.analysis[,vsd.analysis$k_means_clust_MTS!="1"]

#choose variable 
v= "k_means_clust_MTS"

#Set Reference Level for Comparison (Control Group)
dds.analysis_2[[v]] <- factor(dds.analysis_2[[v]])

#change experiment design if needed
design(dds.analysis_2) <- ~ k_means_clust_MTS

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(colData(dds.analysis_2)) 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#619CFF"
col2 <- "#00ba38"


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis_2)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#group <- factor(group, levels = c("3", "2"))
#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("Results/edgeR.results_", paste0(v, paste0("_3_vs_2_MTS.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_3_vs_2_MTS.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/edgeR_", paste0(v, paste0(alpha, paste0("_3_vs_2_MTS.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()
















########## IPA results of k means MG 

IPA_res <- read.csv(file = "Results/IPA/k_means_clusters/MG_k_means_comparison.csv")

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
max(IPA_res$MG_2_vs_3, na.rm = TRUE)
min(IPA_res$MG_2_vs_3,na.rm = TRUE)
max(IPA_res$MG_2_vs_1,na.rm = TRUE)
min(IPA_res$MG_2_vs_1, na.rm = TRUE)
max(IPA_res$MG_3_vs_1, na.rm = TRUE)
min(IPA_res$MG_3_vs_1, na.rm = TRUE)

col_fun = colorRamp2(c(-6.621, 0, 2.828), c("blue", "white", "orange"))
col_fun(seq(-6.621, 0, 2.828))
#sort
IPA_res <- IPA_res %>% arrange(desc(MG_2_vs_3))

#convert data to matrix 
IPA_res_mat <- IPA_res



#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "Figures/IPA_heatmap_k_means_clusters_MG_comparison_fdr_0.2.pdf", height = 12, width = 14)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))
dev.off()




#### with 3 vs 2 

IPA_res <- read.csv(file = "Results/IPA/k_means_clusters/MG_k_means_comparison_upd.csv")

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
max(IPA_res$BAL_vs_SPT, na.rm = TRUE)
min(IPA_res$BAL_vs_SPT,na.rm = TRUE)
max(IPA_res$BAL_vs_BPT,na.rm = TRUE)
min(IPA_res$BAL_vs_BPT, na.rm = TRUE)
max(IPA_res$SPT_vs_BPT, na.rm = TRUE)
min(IPA_res$SPT_vs_BPT, na.rm = TRUE)

col_fun = colorRamp2(c(-4.379, 0, 6.621), c("blue", "white", "orange"))
col_fun(seq(-4.379, 0, 6.621))
#sort
IPA_res <- IPA_res %>% arrange(desc(BAL_vs_SPT))

#convert data to matrix 
IPA_res_mat <- IPA_res



#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "Figures/IPA_heatmap_k_means_clusters_MG_comparison_fdr_0.2_upd.pdf", height = 12, width = 14)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))
dev.off()


######### IPA results of K means MTS 

IPA_res <- read.csv(file = "Results/IPA/k_means_clusters/MTS_k_means_comparison.csv")

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
max(IPA_res$MTS_2_vs_3, na.rm = TRUE)
min(IPA_res$MTS_2_vs_3,na.rm = TRUE)
max(IPA_res$MTS_2_vs_1,na.rm = TRUE)
min(IPA_res$MTS_2_vs_1, na.rm = TRUE)
max(IPA_res$MTS_3_vs_1, na.rm = TRUE)
min(IPA_res$MTS_3_vs_1, na.rm = TRUE)

col_fun = colorRamp2(c(-12.06, 0, 6.808), c("blue", "white", "orange"))
col_fun(seq(-12.06, 0, 6.808))
#sort
IPA_res <- IPA_res %>% arrange(desc(MTS_3_vs_1))

#convert data to matrix 
IPA_res_mat <- IPA_res



#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "Figures/IPA_heatmap_k_means_clusters_MTS_comparison_fdr_0.2.pdf", height = 28, width = 14)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))
dev.off()



#### with 3 vs 2 

IPA_res <- read.csv(file = "Results/IPA/k_means_clusters/MTS_k_means_comparison_upd.csv")

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
max(IPA_res$BAL_vs_SPT, na.rm = TRUE)
min(IPA_res$BAL_vs_SPT,na.rm = TRUE)
max(IPA_res$BAL_vs_BPT,na.rm = TRUE)
min(IPA_res$BAL_vs_BPT, na.rm = TRUE)
max(IPA_res$SPT_vs_BPT, na.rm = TRUE)
min(IPA_res$SPT_vs_BPT, na.rm = TRUE)

col_fun = colorRamp2(c(-6.466, 0, 12.06), c("blue", "white", "orange"))
col_fun(seq(-6.466, 0, 12.06))
#sort
IPA_res <- IPA_res %>% arrange(desc(BAL_vs_SPT))

#convert data to matrix 
IPA_res_mat <- IPA_res



#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "Figures/IPA_heatmap_k_means_clusters_MTS_comparison_fdr_0.2_upd.pdf", height = 28, width = 14)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))
dev.off()



















###### random forest on RNA data ##### 
#use normalized counts table 
counts <- mycounts
counts=counts[apply(counts, 1, function(x) sum(x>=300)>=3),]
counts=apply(counts, 2, function(x) return(log(x+1)))

#run the model. change top_n every time you want to examine top taxa
set.seed(1234)

#get metadata 
metadata.for.RF <- metadata_deseq %>% 
  data.frame() %>% 
    mutate(one_year_mort = ifelse(death=="1" & OS_days > 365, "greater_12", 
                                       ifelse(death=="1" & OS_days<365, "less_12", 
                                              ifelse(death=="0" & OS_days>365, "greater_12", "NA")))) %>% 
    filter(one_year_mort!="NA") %>% 
    mutate(outcome=one_year_mort) %>% 
    mutate(outcome = ifelse(outcome=="greater_12", "0", "1"))


#get features (taxa) 
features <- counts %>% data.frame(.)
colnames(features) <- gsub("X", "", colnames(features))
features <- features %>% select(rownames(metadata.for.RF))

#define outcome 
outcome <- as.factor(metadata.for.RF$one_year_mort)
outcome <- ifelse(outcome=="greater_12", "0","1") 


IVs=list(); AUCs=list()

for (i in 1:20){
  ### stratify 5-fold CV (80% discovery and 20 validation)
  train_index <- c(sample(which(metadata.for.RF$outcome==0), round(length(which(metadata.for.RF$outcome==0))*0.8)),
                   sample(which(metadata.for.RF$outcome!=0), round(length(which(metadata.for.RF$outcome!=0))*0.8)))
  train_data <- features[,train_index]
  meta.training <- metadata.for.RF[train_index,]
  #test_data <- rf.data_complete[-train_index, ]
  
  #tune 
  mtry <- sqrt(ncol(train_data))
  tunegrid <- expand.grid(.mtry=mtry)
  metric <- "ROC"
  control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)
  
  #run model 
  fit <- randomForest(outcome~., data = data.frame(outcome=as.factor(meta.training$outcome), 
                                                   t(train_data)), 
                      trControl=control, metric=metric, 
                      ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)
  
  #get variable importance
  # Extract the mean decrease impurity (MDI) values for each feature
  mdi <- randomForest::importance(fit) %>% as.data.frame()
  
  # Normalize the MDI values 
  Norm.Gini= mdi$MeanDecreaseGini/sum(mdi$MeanDecreaseGini)
  Norm.Gini=Norm.Gini/max(Norm.Gini)
  mdi=data.frame(taxa=rownames(train_data),Norm.Gini)
  
  #add enrich group (or outcome) -->
  # features.means <- by(t(Genus.training), meta.training$delta, colMeans) -->
  # features.means <- do.call(cbind, features.means) -->
  # idx_enrich <- apply(features.means, 1, which.max) -->
  # group_enrich <- colnames(features.means)[idx_enrich] -->
  # -->
  # mdi$enrich_group <- group_enrich -->
  mdi=mdi[order(mdi$Norm.Gini, decreasing = TRUE),]
  
  #cutoffs=c(seq(0.01, 0.09,0.01), seq(0.1, 1,0.1))
  
  cutoffs=c(seq(0.01, 0.09,0.01), seq(0.1, 1,0.1))
  
  VIs.taxon=matrix(NA,nrow = nrow(mdi), ncol=length(cutoffs))
  rownames(VIs.taxon)=mdi$taxa
  
  AUCs.taxon=NULL
  
  for(tt in 1:length(cutoffs)) {
    
    taxa.top=mdi$taxa[1:round(nrow(mdi)*cutoffs[tt])]
    
    fit.per <- randomForest(outcome~., data = data.frame(outcome=as.factor(meta.training$outcome), 
                                                         t(train_data[taxa.top,])),
                            ntree=500, importance=TRUE)
    
    #get variable importance 
    # Extract the mean decrease impurity (MDI) values for each feature
    mdi.per <- randomForest::importance(fit.per) %>% as.data.frame()
    
    # Normalize the MDI values
    Norm.Gini= mdi.per$MeanDecreaseGini/sum(mdi.per$MeanDecreaseGini)
    Norm.Gini=Norm.Gini/max(Norm.Gini)
    
    mdi.per=data.frame(taxa=taxa.top,Norm.Gini)
    
    #add enrich group (or outcome)
    # features.means <- by(t(Genus.training[ASV.top,]), meta.training$delta, colMeans) 
    # features.means <- do.call(cbind, features.means)
    # idx_enrich <- apply(features.means, 1, which.max)
    # group_enrich <- colnames(features.means)[idx_enrich] 
    # mdi.per$enrich_group <- group_enrich
    mdi.per=mdi.per[order(mdi.per$Norm.Gini, decreasing = TRUE),] 
    VIs.taxon[mdi.per$taxa,tt]=mdi.per$Norm.Gini
    
    #get predictions
    predictions <- predict(fit.per, newdata = data.frame(outcome=as.factor(metadata.for.RF[-train_index,]$outcome),
                                                         t(features[taxa.top,-train_index])), type = "prob")
    predictions <- data.frame(predictions)
    roc <- roc(metadata.for.RF[-train_index,]$outcome, predictions$X1)
    
    AUCs.taxon=c(AUCs.taxon, roc$auc) 
  }
  IVs=c(IVs, list(VIs.taxon))
  
  AUCs=c(AUCs, list(AUCs.taxon)) 
}

aa=list(IVs, AUCs)


### plot AUC values 
AUCs=matrix(unlist(aa[[2]]),nrow=20, byrow=TRUE)
AUCs=apply(AUCs, 2, function(x) return(ifelse(x<0.5, 1-x, x)))


AUCs_taxon=AUCs[seq(1,20,2),];cutoffs=c(seq(0.01, 0.09,0.01), seq(0.1, 1,0.1))*100

AUCs_taxon=melt(data.frame(cutoffs, t(AUCs_taxon)), id.vars=1)[,-2]

#plot
pdf(file = "Figures/AUC_gene_level_RF.pdf", height = 6, width =  )
ggline(AUCs_taxon, x = "cutoffs", y = "value",  position = position_dodge(0.8), color = "dodgerblue",
       add = c("mean_sd"),numeric.x.axis = TRUE)+theme_bw(base_size = 14)+
  ylab("AUC")+xlab("% Top Gene")+
  scale_x_continuous(breaks = seq(0,100, by = 10))+
  #scale_color_manual(values = c("red", "blue")) +  # Adjust line colors
  theme(legend.position = "top", legend.title = element_blank())+
  ggtitle("AUC  at gene level (repetition=20)")
dev.off()


# Mean on all columns
df2 <- AUCs_taxon %>% group_by(cutoffs) %>%
  summarise(across(everything(), mean),
            .groups = 'drop')  %>%
  as.data.frame()

#get best cutoff for AUC 
df2$cutoffs[which.max(df2$value)] #80




#run the model for the best AUC 
VIs.taxon=matrix(NA,nrow = nrow(features), ncol=20)
rownames(VIs.taxon)=rownames(features)

AUCs=NULL

for (i in 1:20){
  ### stratify 5-fold CV (80% discovery and 20 validation)
  train_index <- c(sample(which(metadata.for.RF$outcome==0), round(length(which(metadata.for.RF$outcome==0))*0.8)),
                   sample(which(metadata.for.RF$outcome!=0), round(length(which(metadata.for.RF$outcome!=0))*0.8)))
  train_data <- features[,train_index]
  meta.training <- metadata.for.RF[train_index,]
  #test_data <- rf.data_complete[-train_index, ]
  
  #tune 
  mtry <- sqrt(ncol(train_data))
  tunegrid <- expand.grid(.mtry=mtry)
  metric <- "ROC"
  control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)
  
  #run model 
  fit <- randomForest(outcome~., data = data.frame(outcome=as.factor(meta.training$outcome), 
                                                   t(train_data)), 
                      trControl=control, metric=metric, 
                      ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)
  
  #get variable importance
  # Extract the mean decrease impurity (MDI) values for each feature
  mdi <- randomForest::importance(fit) %>% as.data.frame()
  
  # Normalize the MDI values 
  Norm.Gini= mdi$MeanDecreaseGini/sum(mdi$MeanDecreaseGini)
  Norm.Gini=Norm.Gini/max(Norm.Gini)
  mdi=data.frame(taxa=rownames(train_data),Norm.Gini)
  mdi=mdi[order(mdi$Norm.Gini, decreasing = TRUE),]
  
  ## top 30%
  taxa.top <- mdi$taxa[1:round(nrow(mdi)*0.8)]
  
  fit.per <- randomForest(outcome~., data = data.frame(outcome=as.factor(meta.training$outcome), 
                                                       t(train_data[taxa.top,])),
                          ntree=500, importance=TRUE)
  
  #get variable importance 
  # Extract the mean decrease impurity (MDI) values for each feature
  mdi.per <- randomForest::importance(fit.per) %>% as.data.frame()
  
  # Normalize the MDI values
  Norm.Gini= mdi.per$MeanDecreaseGini/sum(mdi.per$MeanDecreaseGini)
  Norm.Gini=Norm.Gini/max(Norm.Gini)
  
  mdi.per=data.frame(taxa=taxa.top,Norm.Gini)
  
  mdi.per=mdi.per[order(mdi.per$Norm.Gini, decreasing = TRUE),] 
  VIs.taxon[mdi.per$taxa,i]=mdi.per$Norm.Gini
  
  #get predictions
  predictions.taxon <- predict(fit.per, newdata = data.frame(outcome=as.factor(metadata.for.RF[-train_index,]$outcome),
                                                             t(features[taxa.top,-train_index])), type = "prob")
  predictions.taxon <- data.frame(predictions.taxon)
  roc <- roc(metadata.for.RF[-train_index,]$outcome, predictions.taxon$X1)
  
  AUCs=c(AUCs, roc$auc) 
}

aa=list(VIs.taxon, AUCs)


AUCs=aa[[2]]
AUCs <- data.frame(AUCs)
AUCs=data.frame(ID=1:nrow(AUCs), AUCs, check.names = FALSE)
AUCs=melt(AUCs, id.vars=1)

pdf(file = "Figures/AUC_best_cutoff_0.8_metatrans_RF.pdf", height = 6, width = 3)
ggboxplot(AUCs, x = "variable", y = "value", color="variable",
          add = "jitter")+theme_bw(base_size = 14)+
  ylab("AUC")+xlab("")+ggtitle("AUCs in the host RNA dataset across 20 repetitions")+
  theme(legend.position = "none")
dev.off()


VI_taxon=aa[[1]]
VI_taxon=VI_taxon[apply(VI_taxon, 1, function(x) sum(!is.na(x)))>0,]
VI_taxon=apply(VI_taxon, 2, function(x) return(ifelse(is.na(x), 0, x)))

VI_taxon=cbind(apply(VI_taxon,1, mean),apply(VI_taxon,1, sd))

colnames(VI_taxon)=c("Average_Gini", "SD")
VI_taxon<- VI_taxon %>% 
  data.frame() %>% 
  arrange(desc(Average_Gini)) %>% 
  mutate(gene=rownames(.))

write.csv(VI_taxon, file="Results/RF_variable_importance_best_AUC_genes.csv")


VI_taxon %>% filter(Average_Gini >0)

































########################################################################
############### repeat analysis of naive patienrts #####################
dds_naive <- dds[,dds$Tx_naive=="yes"]
vsd_naive <- vsd[,vsd$Tx_naive=="yes"]
metadata_naive <- metadata_deseq %>% dplyr::filter(Tx_naive=="yes")

#choose table to use
dds.analysis <- dds_naive
vsd.analysis <- vsd_naive

#choose histology variable 
dds.analysis$histology <- dds.analysis$histology
vsd.analysis$histology <- vsd.analysis$histology


#choose variable 
v= "histology"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "adeno")

#change experiment design if needed
#design(dds.analysis) <- ~ Cavitary

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_naive

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "red"
col2 <- "dodgerblue"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$histology

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("RNA_analysis_naive/Results/Beta.Diversity.Bray.", paste0(v, paste0(".txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis_naive/Figures/Beta.Diversity.Bray_", paste0(v, paste0(".pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Adenocarcinoma", "SCC")), size=10) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0(".csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0(".rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0(".pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()















##############comparison according to pneumitis from immunotherapy: yes vs no #################
######################Prepare data to run funcitons later 
#choose table to use
dds.analysis <- dds_naive
vsd.analysis <- vsd_naive


#choose histology variable 
dds.analysis$irAE_pneumonitis_ever_new <- dds.analysis$irAE_pneumonitis_ever_new
vsd.analysis$irAE_pneumonitis_ever_new <- vsd.analysis$irAE_pneumonitis_ever_new


#choose variable 
v= "irAE_pneumonitis_ever_new"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "0")

#change experiment design if needed
design(dds.analysis) <- ~ irAE_pneumonitis_ever_new

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_naive

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "yellow3"
col2 <- "green3"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$irAE_pneumonitis_ever_new
newResults$v <- factor(newResults$v)

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("RNA_analysis_naive/Results/Beta.Diversity.Bray.", paste0(v, paste0(".txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis_naive/Figures/Beta.Diversity.Bray_", paste0(v, paste0(".pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Pneumonitis", "Pneumonitis")), size=10) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0(".csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0(".rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0(".pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()





##############comparison according to mutations: yes vs no #################
######################Prepare data to run funcitons later 
dds.analysis <- dds_naive
vsd.analysis <- vsd_naive


#choose table to use
dds.analysis <- dds.analysis[,dds.analysis$Driver_mutation_present != "n.a" ]
vsd.analysis <- vsd.analysis[,vsd.analysis$Driver_mutation_present!="n.a"]

#choose histology variable 
dds.analysis$Driver_mutation_present <- dds.analysis$Driver_mutation_present
vsd.analysis$Driver_mutation_present <- vsd.analysis$Driver_mutation_present


#choose variable 
v= "Driver_mutation_present"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "0")

#change experiment design if needed
design(dds.analysis) <- ~ Driver_mutation_present

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_naive %>% filter(Driver_mutation_present!="n.a")

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "orangered"
col2 <- "skyblue"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$Driver_mutation_present
newResults$v <- factor(newResults$v)

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("RNA_analysis_naive/Results/Beta.Diversity.Bray.", paste0(v, paste0(".txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis_naive/Figures/Beta.Diversity.Bray_", paste0(v, paste0(".pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Driver Mutation", "Driver Mutation Present")), size=10) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0(".csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0(".rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0(".pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()






##############comparison according to stage: late vs early stage #################
######################Prepare data to run funcitons later 
#choose table to use
dds.analysis <- dds_naive
vsd.analysis <- vsd_naive

#choose histology variable 
dds.analysis$early_vs_late_stage <- dds.analysis$early_vs_late_stage
vsd.analysis$early_vs_late_stage <- vsd.analysis$early_vs_late_stage


#choose variable 
v= "early_vs_late_stage"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "early")

#change experiment design if needed
design(dds.analysis) <- ~ early_vs_late_stage

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_naive %>% filter(early_vs_late_stage!="n.a")

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "orange2"
col2 <- "grey"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$early_vs_late_stage
newResults$v <- factor(newResults$v)

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("RNA_analysis_naive/Results/Beta.Diversity.Bray.", paste0(v, paste0(".txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis_naive/Figures/Beta.Diversity.Bray_", paste0(v, paste0(".pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Early Stage", "Late Stage")), size=10) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0(".csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0(".rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0(".pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()




##############comparison according to PD1: low vs negative (1-50% vs < 1%) #################

######################Prepare data to run funcitons later 
dds.analysis <- dds_naive
vsd.analysis <- vsd_naive

#choose table to use
dds.analysis <- dds.analysis[,dds.analysis$PD_L1_expression %in% c("Negative", "Low")]
vsd.analysis <- vsd.analysis[,vsd.analysis$PD_L1_expression %in% c("Negative", "Low")]

#choose histology variable 
dds.analysis$PD_L1_expression <- dds.analysis$PD_L1_expression
vsd.analysis$PD_L1_expression <- vsd.analysis$PD_L1_expression


#choose variable 
v= "PD_L1_expression"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "Negative")

#change experiment design if needed
#design(dds.analysis) <- ~ Cavitary

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_naive %>% filter(PD_L1_expression %in% c("Negative", "Low"))

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "orange3"
col2 <- "grey1"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$PD_L1_expression
newResults$v <- factor(newResults$v, levels = c("Negative", "Low"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("RNA_analysis_naive/Results/Beta.Diversity.Bray.", paste0(v, paste0("low_vs_negative.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis_naive/Figures/Beta.Diversity.Bray_", paste0(v, paste0("low_vs_negative.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Negative", "Low")), size=10) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0("low_vs_negative.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("low_vs_negative.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0("low_vs_negative.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()




##############comparison according to PD1: high vs low (> 50% vs  1-50%) #################
dds.analysis <- dds_naive
vsd.analysis <- vsd_naive
######################Prepare data to run funcitons later 
#choose table to use
dds.analysis <- dds.analysis[,dds.analysis$PD_L1_expression %in% c("High", "Low")]
vsd.analysis <- vsd.analysis[,vsd.analysis$PD_L1_expression %in% c("High", "Low")]

#choose histology variable 
dds.analysis$PD_L1_expression <- dds.analysis$PD_L1_expression
vsd.analysis$PD_L1_expression <- vsd.analysis$PD_L1_expression


#choose variable 
v= "PD_L1_expression"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "Low")

#change experiment design if needed
#design(dds.analysis) <- ~ Cavitary

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_naive %>% filter(PD_L1_expression %in% c("High", "Low"))

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "red3"
col2 <- "orange3"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$PD_L1_expression
newResults$v <- factor(newResults$v, levels = c("Low", "High"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("RNA_analysis_naive/Results/Beta.Diversity.Bray.", paste0(v, paste0("_High_vs_Low.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis_naive/Figures/Beta.Diversity.Bray_", paste0(v, paste0("_High_vs_Low.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Low", "High")), size=10) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)
group <- factor(group, levels=c("Low", "High"))

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0("_High_vs_Low.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_High_vs_Low.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_High_vs_Low.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()














##############comparison according to PD1: high vs negative (> 50% vs  <1%) #################

dds.analysis <- dds_naive
vsd.analysis <- vsd_naive
######################Prepare data to run funcitons later 
#choose table to use
dds.analysis <- dds.analysis[,dds.analysis$PD_L1_expression %in% c("High", "Negative")]
vsd.analysis <- vsd.analysis[,vsd.analysis$PD_L1_expression %in% c("High", "Negative")]

#choose histology variable 
dds.analysis$PD_L1_expression <- dds.analysis$PD_L1_expression
vsd.analysis$PD_L1_expression <- vsd.analysis$PD_L1_expression


#choose variable 
v= "PD_L1_expression"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "Negative")

#change experiment design if needed
design(dds.analysis) <- ~ PD_L1_expression

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_naive %>% filter(PD_L1_expression %in% c("High", "Negative"))

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "red3"
col2 <- "grey1"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$PD_L1_expression
newResults$v <- factor(newResults$v, levels = c("Negative", "High"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("RNA_analysis_naive/Results/Beta.Diversity.Bray.", paste0(v, paste0("_High_vs_Negative.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis_naive/Figures/Beta.Diversity.Bray_", paste0(v, paste0("_High_vs_Negative.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Negative", "High")), size=10) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)
group <- factor(group, levels = c("Negative", "High"))
#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0("_High_vs_Negative.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_High_vs_Negative.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_High_vs_Negative.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()


######################### comparison of PDL1 of 2 levels ########################









##############comparison according to mortality: death vs alive  #################

######################Prepare data to run funcitons later 
#choose table to use
dds.analysis <- dds_naive
vsd.analysis <- vsd_naive

#choose histology variable 
dds.analysis$death <- dds.analysis$death
vsd.analysis$death <- vsd.analysis$death


#choose variable 
v= "death"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "0")

#change experiment design if needed
#design(dds.analysis) <- ~ Cavitary

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_naive

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "black"
col2 <- "green3"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$death
newResults$v <- factor(newResults$v)

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("RNA_analysis_naive/Results/Beta.Diversity.Bray.", paste0(v, paste0("_dead_vs_alive.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis_naive/Figures/Beta.Diversity.Bray_", paste0(v, paste0("_dead_vs_alive.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Alive", "Dead")), size=10) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0("_dead_vs_alive.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_dead_vs_alive.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_dead_vs_alive.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()





##################### Comparison of PDL1 only two groups > or < 50% ###############

######################Prepare data to run funcitons later 

#choose  variable 
dds_naive$PD_L1_expression_2_lev <- ifelse(dds_naive$PD_L1_expression_percentage < 50, "low_50", "high_50")
vsd_naive$PD_L1_expression_2_lev <- ifelse(vsd_naive$PD_L1_expression_percentage < 50, "low_50", "high_50")

#choose table to use
dds.analysis <- dds_naive
vsd.analysis <- vsd_naive

#choose variable 
v= "PD_L1_expression_2_lev"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "low_50")

#change experiment design if needed
design(dds.analysis) <- ~ PD_L1_expression_2_lev

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_naive 
metadata_RNA_analysis$PD_L1_expression_2_lev <- ifelse(metadata_RNA_analysis$PD_L1_expression_percentage < 50, "low_50", "high_50")

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "orange4"
col2 <- "grey50"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$PD_L1_expression_2_lev
newResults$v <- factor(newResults$v, levels = c("low_50", "high_50"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("RNA_analysis_naive/Results/Beta.Diversity.Bray.", paste0(v, paste0("_high_vs_low.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis_naive/Figures/Beta.Diversity.Bray_", paste0(v, paste0("_high_vs_low.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("<50%", ">50%")), size=10) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)
group <- factor(group, levels = c("low_50", "high_50"))
#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0("_high_vs_low.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_high_vs_low.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_high_vs_low.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()









######### one year mortality #####

######################Prepare data to run funcitons later 

#choose  variable 
dds.analysis <- dds_naive
vsd.analysis <- vsd_naive
#choose table to use
dds.analysis <- dds.analysis[,dds.analysis$one_y_mort!="NA"]
vsd.analysis <- vsd.analysis[,vsd.analysis$one_y_mort!="NA"]

#choose variable 
v= "one_y_mort"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "greater_12")

#change experiment design if needed
design(dds.analysis) <- ~ one_y_mort

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_naive %>% filter(one_y_mort!="NA") 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "red3"
col2 <- "orange3"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$one_y_mort
newResults$v <- factor(newResults$v, levels = c("greater_12", "less_12"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("RNA_analysis_naive/Results/Beta.Diversity.Bray.", paste0(v, paste0("_less_vs_greater_12.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis_naive/Figures/Beta.Diversity.Bray_", paste0(v, paste0("_less_vs_greater_12.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c(">12 Months", "< 12 Months")), size=10) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0("_less_vs_greater_12_mon.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_less_vs_greater_12_mon.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_less_vs_greater_12_mon.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()







######### two year mortality #####

######################Prepare data to run funcitons later 

dds.analysis <- dds_naive
vsd.analysis <- vsd_naive
#choose  variable 

#choose table to use
dds.analysis <- dds.analysis[,dds.analysis$two_y_mort!="NA"]
vsd.analysis <- vsd.analysis[,vsd.analysis$two_y_mort!="NA"]

#choose variable 
v= "two_y_mort"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "greater_24")

#change experiment design if needed
design(dds.analysis) <- ~ two_y_mort

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_naive %>% filter(two_y_mort!="NA") 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "red3"
col2 <- "orange3"

######################Plot PCoA

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$two_y_mort
newResults$v <- factor(newResults$v, levels = c("greater_24", "less_24"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("RNA_analysis_naive/Results/Beta.Diversity.Bray.", paste0(v, paste0("_less_vs_greater_24.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis_naive/Figures/Beta.Diversity.Bray_", paste0(v, paste0("_less_vs_greater_24.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c(">24 Months", "< 24 Months")), size=10) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0("_less_vs_greater_24_mon.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_less_vs_greater_24_mon.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_less_vs_greater_24_mon.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()
















############################# k means clusters analysis for naive patients usiing original analysis #############

##### metagenomics 

#read MG clusters 
mg_data <- read.csv(file = "MG_analysis_naive/MG_k_means_clusters_for_RNA.csv")
#select only ID to merge with metadata 
mg_data_clust <- mg_data %>%  dplyr::select(ID, k_means_3_cluster_pneumotype)
#fix ID names to match those of dds metadata 
mg_data_clust$ID <- gsub("_", "." ,mg_data_clust$ID)
#make sure names are the same 
table(colnames(dds_naive) == mg_data_clust$ID)

dds.analysis <- dds_naive
vsd.analysis <- vsd_naive

#add k means clusters to metadata 
dds.analysis$k_means_clust_MG <- mg_data_clust$k_means_3_cluster_pneumotype
vsd.analysis$k_means_clust_MG <- mg_data_clust$k_means_3_cluster_pneumotype

#choose variable 
v= "k_means_clust_MG"

#### beta diversity for all 

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$k_means_clust_MG
newResults$v <- factor(newResults$v, levels = c("BPT", "SPT", "BAL"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("RNA_analysis_naive/Results/Beta.Diversity.Bray.", paste0(v, paste0("_MG.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis_naive/Figures/Beta.Diversity.Bray_", paste0(v, paste0("_all_MG.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  #scale_color_manual(values=c("blue", "goldenrod", "grey")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Cluster1", "Cluster2", "Cluster3"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


#1 vs 2 (p=0.67)
ps.re.1 <- vsd.analysis[,vsd.analysis$k_means_clust_MG!="BAL"]
vegdist = vegdist(t(assay(ps.re.1)), method = "bray")
adonis2(vegdist~colData(ps.re.1)$k_means_clust_MG)

#2 vs 3 (p=0.189)
ps.re.1 <- vsd.analysis[,vsd.analysis$k_means_clust_MG!="BPT"]
vegdist = vegdist(t(assay(ps.re.1)), method = "bray")
adonis2(vegdist~colData(ps.re.1)$k_means_clust_MG)

#1 vs 3 (p=0.194)
ps.re.1 <- vsd.analysis[,vsd.analysis$k_means_clust_MG!="SPT"]
vegdist = vegdist(t(assay(ps.re.1)), method = "bray")
adonis2(vegdist~colData(ps.re.1)$k_means_clust_MG)

#compare between the clusters 


######edgeR analysis  SPT vs BPT


#choose table to use
dds.analysis_2 <- dds.analysis[,dds.analysis$k_means_clust_MG!="BAL"]
vsd.analysis_2 <- vsd.analysis[,vsd.analysis$k_means_clust_MG!="BAL"]

#choose variable 
v= "k_means_clust_MG"

#Set Reference Level for Comparison (Control Group)
dds.analysis_2[[v]] <- factor(dds.analysis_2[[v]], levels = c("BPT", "SPT"))

#change experiment design if needed
design(dds.analysis_2) <- ~ k_means_clust_MG

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(colData(dds.analysis_2)) 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#00ba38"
col2 <- "#f8776d"


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis_2)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0("_SPT_vs_BPT_MG.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_SPT_vs_BPT_MG.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_SPT_vs_BPT_MG.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()






########## repeat edgeR BAL vs BPT

#choose table to use
dds.analysis_2 <- dds.analysis[,dds.analysis$k_means_clust_MG!="SPT"]
vsd.analysis_2 <- vsd.analysis[,vsd.analysis$k_means_clust_MG!="SPT"]

#choose variable 
v= "k_means_clust_MG"

#Set Reference Level for Comparison (Control Group)
dds.analysis_2[[v]] <- factor(dds.analysis_2[[v]], levels = c("BPT", "BAL"))

#change experiment design if needed
design(dds.analysis_2) <- ~ k_means_clust_MG

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(colData(dds.analysis_2)) 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#619CFF"
col2 <- "#f8776d"


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis_2)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0("_BAL_vs_BPT_MG.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_BAL_vs_BPT_MG.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_BAL_vs_BPT_MG.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()







########## repeat edgeR BAL vs  SPT



#choose table to use
dds.analysis_2 <- dds.analysis[,dds.analysis$k_means_clust_MG!="BPT"]
vsd.analysis_2 <- vsd.analysis[,vsd.analysis$k_means_clust_MG!="BPT"]

#choose variable 
v= "k_means_clust_MG"

#Set Reference Level for Comparison (Control Group)
dds.analysis_2[[v]] <- factor(dds.analysis_2[[v]], levels = c("SPT", "BAL"))

#change experiment design if needed
design(dds.analysis_2) <- ~ k_means_clust_MG

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(colData(dds.analysis_2)) 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#619CFF"
col2 <-  "#00ba38"
  

##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis_2)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0("_BAL_vs_SPT_MG.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_BAL_vs_SPT_MG.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_BAL_vs_SPT_MG.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()






#########MTS k means clust, 



#read MTS clusters 
MTS_data <- read.csv(file = "MTS_analysis_naive/MTS_k_means_clusters_for_RNA.csv")
#select only ID to merge with metadata 
MTS_data_clust <- MTS_data %>% dplyr::select(ID, MTS_k_means_cluster_sample_type_upd)

#fix ID names to match those of dds metadata 
MTS_data_clust$ID <- gsub("_", "." ,MTS_data_clust$ID)
col_names <- colnames(dds.analysis) 
rownames(MTS_data_clust) <- MTS_data_clust$ID
MTS_data_clust <- MTS_data_clust[col_names,]
#make sure names are the same 
table(colnames(dds.analysis) == MTS_data_clust$ID)
#add k means clusters to metadata 
dds.analysis$k_means_clust_MTS <- MTS_data_clust$MTS_k_means_cluster_sample_type_upd
vsd.analysis$k_means_clust_MTS <- MTS_data_clust$MTS_k_means_cluster_sample_type_upd

dds.analysis$k_means_clust_MTS <- factor(dds.analysis$k_means_clust_MTS, levels = c("BPT", "SPT", "BAL"))

#choose variable 
v= "k_means_clust_MTS"

#### beta diversity for all 

#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$k_means_clust_MTS
newResults$v <- factor(newResults$v, levels = c("BPT", "SPT", "BAL"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("RNA_analysis_naive/Results/Beta.Diversity.Bray.", paste0(v, paste0("_MTS.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis_naive/Figures/Beta.Diversity.Bray_", paste0(v, paste0("all_MTS.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  #scale_color_manual(values=c("blue", "goldenrod", "grey")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("BPT", "SPT", "BAL")), size=10) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()



#1 vs 2 (p=0.109)
ps.re.1 <- vsd.analysis[,vsd.analysis$k_means_clust_MG!="BPT"]
vegdist = vegdist(t(assay(ps.re.1)), method = "bray")
adonis2(vegdist~colData(ps.re.1)$k_means_clust_MG)

#2 vs 3 (p=0.118)
ps.re.1 <- vsd.analysis[,vsd.analysis$k_means_clust_MG!="SPT"]
vegdist = vegdist(t(assay(ps.re.1)), method = "bray")
adonis2(vegdist~colData(ps.re.1)$k_means_clust_MG)

#1 vs 3 (p=0.117)
ps.re.1 <- vsd.analysis[,vsd.analysis$k_means_clust_MG!="BAL"]
vegdist = vegdist(t(assay(ps.re.1)), method = "bray")
adonis2(vegdist~colData(ps.re.1)$k_means_clust_MG)


#compare between the clusters 


######edgeR analysis  SPT vs BPT


#choose table to use
dds.analysis_2 <- dds.analysis[,dds.analysis$k_means_clust_MTS!="BAL"]
vsd.analysis_2 <- vsd.analysis[,vsd.analysis$k_means_clust_MTS!="BAL"]

#choose variable 
v= "k_means_clust_MTS"

#Set Reference Level for Comparison (Control Group)
dds.analysis_2[[v]] <- factor(dds.analysis_2[[v]])

#change experiment design if needed
design(dds.analysis_2) <- ~ k_means_clust_MTS

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(colData(dds.analysis_2)) 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#00ba38"
col2 <- "#f8776d"


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis_2)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0("_SPT_vs_BPT_MTS.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_spt_vs_bpt_MTS.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_SPT_vs_BPT_MTS.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()






########## repeat edgeR BAL VS BPT

#choose table to use
dds.analysis_2 <- dds.analysis[,dds.analysis$k_means_clust_MTS!="SPT"]
vsd.analysis_2 <- vsd.analysis[,vsd.analysis$k_means_clust_MTS!="SPT"]

#choose variable 
v= "k_means_clust_MTS"

#Set Reference Level for Comparison (Control Group)
dds.analysis_2[[v]] <- factor(dds.analysis_2[[v]], levels = c("BPT", "BAL"))

#change experiment design if needed
design(dds.analysis_2) <- ~ k_means_clust_MTS

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(colData(dds.analysis_2)) 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#619CFF"
col2 <- "#f8776d"


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis_2)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0("_bal_vs_bpt_MTS.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_bal_vs_bpt_MTS.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_bal_vs_bpt_MTS.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()









######### bal vs SPT 


#choose table to use
dds.analysis_2 <- dds.analysis[,dds.analysis$k_means_clust_MTS!="BPT"]
vsd.analysis_2 <- vsd.analysis[,vsd.analysis$k_means_clust_MTS!="BPT"]

#choose variable 
v= "k_means_clust_MTS"

#Set Reference Level for Comparison (Control Group)
dds.analysis_2[[v]] <- factor(dds.analysis_2[[v]], levels = c("SPT", "BAL"))

#change experiment design if needed
design(dds.analysis_2) <- ~ k_means_clust_MTS

#create appropraite metadata for edgeR
metadata_RNA_analysis <- data.frame(colData(dds.analysis_2)) 

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "#619CFF"
col2 <- "#00ba38"


##############Differential Analysis 
mycounts.analysis <- assay(dds.analysis_2)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#group <- factor(group, levels = c("3", "2"))
#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0("_bal_vs_spt_MTS.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_bal_vs_spt_MTS.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_bal_vs_spt_MTS.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()


#### add info of clsuters to metadata naive 
metadata_naive$k_means_clust_MTS <- MTS_data_clust$MTS_k_means_cluster_sample_type_upd
metadata_naive$k_means_clust_MG <- mg_data_clust$k_means_3_cluster_pneumotype



####### compare counts of specefic genes in naive patints": IL17, IL8 , IL1 ######
#plot mRNA levels of il-17 

temp <- data.frame(assay(dds_naive))

meta <- metadata_naive

IL17_plot_data <- temp %>% t(.) %>%data.frame() %>% select(starts_with("IL17"))
rownames(IL17_plot_data) <- gsub("X", "", rownames(IL17_plot_data))
#add to this metadata 
IL17_plot_data <- merge(IL17_plot_data, meta, by="row.names", all=TRUE)
IL17_plot_data <- IL17_plot_data %>% mutate(k_means_clust_MG=factor(k_means_clust_MG, levels=c("BPT", "SPT", "BAL")))

#plot for all il17 measures 

col_names <- colnames(IL17_plot_data)[2:7]

for (i in col_names){
  p<-  ggplot(IL17_plot_data, aes_string(x = IL17_plot_data$k_means_clust_MG, y = i)) +
    geom_boxplot(width = 0.2, color = "blue") +
    geom_jitter(width = 0.2,alpha = 0.5,color = "darkblue") +
    scale_y_log10()+
    stat_compare_means(comparisons = list(c("BPT", "SPT"), 
                                          c("SPT","BAL"), 
                                          c("BPT", "BAL")))+
    ggtitle(i)+
    theme_bw()
  
  #save it 
  pdf_output <- paste0(paste0("RNA_analysis_naive/Figures/", paste0(i, paste0("_IL_17_boxplot", paste0("_MG_k_means_clusters_naive",paste0(".pdf"))))))
  pdf(file = pdf_output, width = 4, height = 6)
  show(p)
  dev.off()
}


####### plot mRNA of IL-8 (CXCL8)#####

IL8_plot_data <- mycounts %>% t(.) %>%data.frame() %>% select(starts_with("CXCL8"))

IL8_plot_data <- temp %>% t(.) %>%data.frame() %>% select(starts_with("CXCL8"))
rownames(IL8_plot_data) <- gsub("X", "", rownames(IL8_plot_data))
#add to this metadata 
IL8_plot_data <- merge(IL8_plot_data, meta, by="row.names", all=TRUE)
IL8_plot_data <- IL8_plot_data %>% mutate(k_means_clust_MG=factor(k_means_clust_MG, levels=c("BPT", "SPT", "BAL")))

#plot for all IL8 measures 

col_names <- colnames(IL17_plot_data)[2:7]

for (i in col_names){
  p<-  ggplot(IL17_plot_data, aes_string(x = IL17_plot_data$k_means_clust_MG, y = i)) +
    geom_boxplot(width = 0.2, color = "blue") +
    geom_jitter(width = 0.2,alpha = 0.5,color = "darkblue") +
    scale_y_log10()+
    stat_compare_means(comparisons = list(c("BPT", "SPT"), 
                                          c("SPT","BAL"), 
                                          c("BPT", "BAL")))+
    ggtitle(i)+
    theme_bw()
  
  #save it 
  pdf_output <- paste0(paste0("RNA_analysis_naive/Figures/", paste0(i, paste0("_IL_17_boxplot", paste0("_MG_k_means_clusters_naive",paste0(".pdf"))))))
  pdf(file = pdf_output, width = 4, height = 6)
  show(p)
  dev.off()
}

####### plot mRNA of IL-1 #####

IL_1_plot_data<- temp %>% t(.) %>%data.frame() %>% select(c("IL1R2", "IL1R1", "IL1RL2", "IL1RL1", "IL1A", "IL1B"))
rownames(IL_1_plot_data) <- gsub("X", "", rownames(IL_1_plot_data))
#add to this metadata 
IL_1_plot_data <- merge(IL_1_plot_data, meta, by="row.names", all=TRUE)
IL_1_plot_data <- IL_1_plot_data %>% mutate(k_means_clust_MG=factor(k_means_clust_MG, levels=c("BPT", "SPT", "BAL")))

#plot for all il1 measures 

col_names <- colnames(IL_1_plot_data)[2:7]

for (i in col_names){
  p<-  ggplot(IL_1_plot_data, aes_string(x = IL_1_plot_data$k_means_clust_MG, y = i)) +
    geom_boxplot(width = 0.2, color = "blue") +
    geom_jitter(width = 0.2,alpha = 0.5,color = "darkblue") +
    scale_y_log10()+
    stat_compare_means(comparisons = list(c("BPT", "SPT"), 
                                          c("SPT","BAL"), 
                                          c("BPT", "BAL")))+
    ggtitle(i)+
    theme_bw()
  
  #save it 
  pdf_output <- paste0(paste0("RNA_analysis_naive/Figures/", paste0(i, paste0("_IL_1_boxplot", paste0("_MG_k_means_clusters_naive",paste0(".pdf"))))))
  pdf(file = pdf_output, width = 4, height = 6)
  show(p)
  dev.off()
}


####### repeat these for MTS 

temp <- data.frame(assay(dds_naive))

meta <- metadata_naive

IL17_plot_data <- temp %>% t(.) %>%data.frame() %>% select(starts_with("IL17"))
rownames(IL17_plot_data) <- gsub("X", "", rownames(IL17_plot_data))
#add to this metadata 
IL17_plot_data <- merge(IL17_plot_data, meta, by="row.names", all=TRUE)
IL17_plot_data <- IL17_plot_data %>% mutate(k_means_clust_MTS=factor(k_means_clust_MTS, levels=c("BPT", "SPT", "BAL")))

#plot for all il17 measures 

col_names <- colnames(IL17_plot_data)[2:7]

for (i in col_names){
  p<-  ggplot(IL17_plot_data, aes_string(x = IL17_plot_data$k_means_clust_MTS, y = i)) +
    geom_boxplot(width = 0.2, color = "blue") +
    geom_jitter(width = 0.2,alpha = 0.5,color = "darkblue") +
    scale_y_log10()+
    stat_compare_means(comparisons = list(c("BPT", "SPT"), 
                                          c("SPT","BAL"), 
                                          c("BPT", "BAL")))+
    ggtitle(i)+
    theme_bw()
  
  #save it 
  pdf_output <- paste0(paste0("RNA_analysis_naive/Figures/", paste0(i, paste0("_IL_17_boxplot", paste0("_MTS_k_means_clusters_naive",paste0(".pdf"))))))
  pdf(file = pdf_output, width = 4, height = 6)
  show(p)
  dev.off()
}


####### plot mRNA of IL-8 (CXCL8)#####

IL8_plot_data <- mycounts %>% t(.) %>%data.frame() %>% select(starts_with("CXCL8"))

IL8_plot_data <- temp %>% t(.) %>%data.frame() %>% select(starts_with("CXCL8"))
rownames(IL8_plot_data) <- gsub("X", "", rownames(IL8_plot_data))
#add to this metadata 
IL8_plot_data <- merge(IL8_plot_data, meta, by="row.names", all=TRUE)
IL8_plot_data <- IL8_plot_data %>% mutate(k_means_clust_MTS=factor(k_means_clust_MTS, levels=c("BPT", "SPT", "BAL")))

#plot for all IL8 measures 
pdf_output <- paste0(paste0("RNA_analysis_naive/Figures/", paste0(i, paste0("_IL_8_boxplot", paste0("_MTS_k_means_clusters_naive",paste0(".pdf"))))))
pdf(file = pdf_output, width = 4, height = 6)
ggplot(IL8_plot_data, aes(x = IL8_plot_data$k_means_clust_MTS, y = CXCL8)) +
    geom_boxplot(width = 0.2, color = "blue") +
    geom_jitter(width = 0.2,alpha = 0.5,color = "darkblue") +
    scale_y_log10()+
    stat_compare_means(comparisons = list(c("BPT", "SPT"), 
                                          c("SPT","BAL"), 
                                          c("BPT", "BAL")))+
    ggtitle("CXCL8")+
    theme_bw()
dev.off()


####### plot mRNA of IL-1 #####

IL_1_plot_data<- temp %>% t(.) %>%data.frame() %>% select(c("IL1R2", "IL1R1", "IL1RL2", "IL1RL1", "IL1A", "IL1B"))
rownames(IL_1_plot_data) <- gsub("X", "", rownames(IL_1_plot_data))
#add to this metadata 
IL_1_plot_data <- merge(IL_1_plot_data, meta, by="row.names", all=TRUE)
IL_1_plot_data <- IL_1_plot_data %>% mutate(k_means_clust_MTS=factor(k_means_clust_MTS, levels=c("BPT", "SPT", "BAL")))

#plot for all il1 measures 

col_names <- colnames(IL_1_plot_data)[2:7]

for (i in col_names){
  p<-  ggplot(IL_1_plot_data, aes_string(x = IL_1_plot_data$k_means_clust_MTS, y = i)) +
    geom_boxplot(width = 0.2, color = "blue") +
    geom_jitter(width = 0.2,alpha = 0.5,color = "darkblue") +
    scale_y_log10()+
    stat_compare_means(comparisons = list(c("BPT", "SPT"), 
                                          c("SPT","BAL"), 
                                          c("BPT", "BAL")))+
    ggtitle(i)+
    theme_bw()
  
  #save it 
  pdf_output <- paste0(paste0("RNA_analysis_naive/Figures/", paste0(i, paste0("_IL_1_boxplot", paste0("_MTS_k_means_clusters_naive",paste0(".pdf"))))))
  pdf(file = pdf_output, width = 4, height = 6)
  show(p)
  dev.off()
}






########## IPA results of k means MG 

IPA_res <- read.csv(file = "Results/IPA/IPA_naive/MG_k_means_comparison.csv")

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
max(IPA_res$SPT_vs_BPT, na.rm = TRUE)
min(IPA_res$SPT_vs_BPT,na.rm = TRUE)
max(IPA_res$BAL_vs_BPT,na.rm = TRUE)
min(IPA_res$BAL_vs_BPT, na.rm = TRUE)
max(IPA_res$BAL_vs_SPT, na.rm = TRUE)
min(IPA_res$BAL_vs_SPT, na.rm = TRUE)

col_fun = colorRamp2(c(-2.714, 0, 5.657), c("blue", "white", "orange"))
col_fun(seq(-2.714, 0, 5.657))
#sort
IPA_res <- IPA_res %>% arrange(desc(SPT_vs_BPT))

#convert data to matrix 
IPA_res_mat <- IPA_res



#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "RNA_analysis_naive/Figures/IPA_heatmap_k_means_clusters_MG_comparison_fdr_0.2.pdf", height = 12, width = 14)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))
dev.off()




######### IPA results of K means MTS 

IPA_res <- read.csv(file = "Results/IPA/IPA_naive/MTS_k_means_comparison.csv")

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
max(IPA_res$SPT_vs_BPT, na.rm = TRUE)
min(IPA_res$SPT_vs_BPT,na.rm = TRUE)
max(IPA_res$BAL_vs_BPT,na.rm = TRUE)
min(IPA_res$BAL_vs_BPT, na.rm = TRUE)
max(IPA_res$BAL_vs_SPT, na.rm = TRUE)
min(IPA_res$BAL_vs_SPT, na.rm = TRUE)

col_fun = colorRamp2(c(-3.13, 0, 6.621), c("blue", "white", "orange"))
col_fun(seq(-3.13, 0, 6.621))
#sort
IPA_res <- IPA_res %>% arrange(desc(SPT_vs_BPT))

#convert data to matrix 
IPA_res_mat <- IPA_res



#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "RNA_analysis_naive/Figures/IPA_heatmap_k_means_clusters_MTS_comparison_fdr_0.2.pdf", height = 18, width = 14)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))
dev.off()



#ordered by BAL vs bpt 
IPA_res <- IPA_res %>% arrange(desc(BAL_vs_BPT))

#convert data to matrix 
IPA_res_mat <- IPA_res


#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "RNA_analysis_naive/Figures/IPA_heatmap_k_means_clusters_MTS_comparison_fdr_0.2_order_bal_vs_bpt.pdf", height = 18, width = 14)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))
dev.off()







































################################################################################
############################ IPA analysis ##########################

######## IPA PD-l1 levels  #####

#read IPA results 
IPA_res <- read.csv(file = "RNA_analysis_naive/Results/IPA/IPA_PDL_1_comparison_heatmap.csv")

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
max(IPA_res$High_vs_Low, na.rm = TRUE)
min(IPA_res$High_vs_Low,na.rm = TRUE)
max(IPA_res$High_vs_Negative,na.rm = TRUE)
min(IPA_res$High_vs_Negative, na.rm = TRUE)
max(IPA_res$Low_vs_Negative, na.rm = TRUE)
min(IPA_res$Low_vs_Negative, na.rm = TRUE)

col_fun = colorRamp2(c(-2.673, 0, 3.938), c("blue", "white", "orange"))
col_fun(seq(-2.673, 0, 3.938))

#convert data to matrix 
IPA_res_mat <- IPA_res

#sort 
IPA_res_mat <- IPA_res %>% arrange(desc(High_vs_Low))

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "RNA_analysis_naive/Figures/IPA_PDL_1_comparison_heatmap.pdf", height = 24, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"), height = unit(38, "cm"))
dev.off()


###############################################################################
###### IPA dead vs alive ######

#read results 
IPA_res <- read_csv("RNA_analysis_naive/Results/IPA/IPA_death_dead_vs_alive_plot.csv")
IPA_res$Pathway <- factor(IPA_res$Pathway, levels = unique(IPA_res$Pathway))
IPA_res <- IPA_res %>% arrange(desc(z_score)) %>% dplyr::filter(z_score != 0)
IPA_res <- IPA_res %>%  mutate(col=ifelse(z_score >0, "orange", "blue"))

pdf(file = "RNA_analysis_naive/Figures/top_pathways_dead_vs_alive.pdf", width = 34, height = 40)
ggplot(IPA_res, aes(x=z_score, y=fct_reorder(Pathway, z_score, .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in dead vs alive patients")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 24, face = "bold"))
dev.off()


-log(0.15)

###### IPA pneumonitis vs not ######

#read results 
IPA_res <- read.csv(file = "RNA_analysis_naive/Results/IPA/updated_03_24_24/pneumonitis_ever_new.csv")
IPA_res <- IPA_res %>% arrange(desc(z.score)) %>% dplyr::filter(z.score != 0)
IPA_res <- IPA_res %>% filter(X.log.p.value.>= -log(0.15))
IPA_res <- IPA_res %>%  mutate(col=ifelse(z.score >0, "orange", "blue"))
nrow(IPA_res)
#select top 20 if high nrow 
top_20_up <- IPA_res %>% dplyr::slice(1:20)
top_20_down <- IPA_res %>% arrange(z.score) %>% dplyr::slice(1:20)
comb_res <- bind_rows(top_20_up, top_20_down)

pdf(file = "RNA_analysis_naive/Figures/IPA_irAE_pneumonitis_vs_not.pdf", width = 22, height = 10)
ggplot(comb_res, aes(x=z.score, y=fct_reorder(Pathway, z.score, .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with pneumonitis vs without")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 16))
dev.off()


###### IPA late vs early stage ######

#read results 
IPA_res <- read.csv(file = "RNA_analysis_naive/Results/IPA/updated_03_24_24/late_vs_early_stage.csv")
IPA_res <- IPA_res %>% arrange(desc(z.score)) %>% dplyr::filter(z.score != 0)
IPA_res <- IPA_res %>% filter(X.log.p.value.>= -log(0.15))
IPA_res <- IPA_res %>%  mutate(col=ifelse(z.score >0, "orange", "blue"))
nrow(IPA_res)
#select top 20 if high nrow 
top_20_up <- IPA_res %>% dplyr::slice(1:20)
top_20_down <- IPA_res %>% arrange(z.score) %>% dplyr::slice(1:20)
comb_res <- bind_rows(top_20_up, top_20_down)

pdf(file = "RNA_analysis_naive/Figures/IPA_late_vs_early_stage.pdf", width =  22, height = 10)
ggplot(comb_res, aes(x=z.score, y=fct_reorder(Pathway, z.score, .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with late vs early stage")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 16))
dev.off()



###### IPA PDL1 high vs low ######

#read results 
IPA_res <- read.csv(file = "RNA_analysis_naive/Results/IPA/updated_03_24_24/PD_L1_high_vs_low.csv")
IPA_res <- IPA_res %>% arrange(desc(z.score)) %>% dplyr::filter(z.score != 0)
IPA_res <- IPA_res %>% filter(X.log.p.value.>= -log(0.15))
IPA_res <- IPA_res %>%  mutate(col=ifelse(z.score >0, "orange", "blue"))
nrow(IPA_res)
#select top 20 if high nrow 
top_20_up <- IPA_res %>% dplyr::slice(1:20)
top_20_down <- IPA_res %>% arrange(z.score) %>% dplyr::slice(1:20)
comb_res <- bind_rows(top_20_up, top_20_down)

pdf(file = "RNA_analysis_naive/Figures/IPA_high_vs_low_PDL1_50_cutoff.pdf", width = 22, height = 10)
ggplot(comb_res, aes(x=z.score, y=fct_reorder(Pathway, z.score, .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with high vs low PDL1")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 16))
dev.off()




###### IPA one year mortaltity  ######

#read results 
IPA_res <- read.csv(file = "RNA_analysis_naive/Results/IPA/updated_03_24_24/one_y_mort_less_vs_greater_12.csv")
IPA_res <- IPA_res %>% arrange(desc(z.score)) %>% dplyr::filter(z.score != 0)
IPA_res <- IPA_res %>% filter(X.log.p.value.>= -log(0.15))
IPA_res <- IPA_res %>%  mutate(col=ifelse(z.score >0, "orange", "blue"))
nrow(IPA_res)
#select top 20 if high nrow 
top_20_up <- IPA_res %>% dplyr::slice(1:20)
top_20_down <- IPA_res %>% arrange(z.score) %>% dplyr::slice(1:20)
comb_res <- bind_rows(top_20_up, top_20_down)

pdf(file = "RNA_analysis_naive/Figures/IPA_one_y_mort_less_vs_greater_12.pdf", width = 22, height = 10)
ggplot(comb_res, aes(x=z.score, y=fct_reorder(Pathway, z.score, .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  #ggtitle("Top Upregulated Pathways in patients with high vs low PDL1")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 16))
dev.off()




###### IPA histology  ######

#read results 
IPA_res <- read.csv(file = "RNA_analysis_naive/Results/IPA/updated_03_24_24/histology_scc_vs_adeno.csv")
IPA_res <- IPA_res %>% arrange(desc(z.score)) %>% dplyr::filter(z.score != 0)
IPA_res <- IPA_res %>% filter(X.log.p.value.>= -log(0.15))
IPA_res <- IPA_res %>%  mutate(col=ifelse(z.score >0, "orange", "blue"))
nrow(IPA_res)
#select top 20 if high nrow 
top_20_up <- IPA_res %>% dplyr::slice(1:20)
top_20_down <- IPA_res %>% arrange(z.score) %>% dplyr::slice(1:20)
comb_res <- bind_rows(top_20_up, top_20_down)

pdf(file = "RNA_analysis_naive/Figures/IPA_histology_scc_vs_adeno.pdf", width = 22, height = 12)
ggplot(comb_res, aes(x=z.score, y=fct_reorder(Pathway, z.score, .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  #ggtitle("Top Upregulated Pathways in patients with high vs low PDL1")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 16))
dev.off()



###### IPA driver mutation  ######

#read results 
IPA_res <- read.csv(file = "RNA_analysis_naive/Results/IPA/updated_03_24_24/driver_mutation.csv")
IPA_res <- IPA_res %>% arrange(desc(z.score)) %>% dplyr::filter(z.score != 0)
IPA_res <- IPA_res %>% filter(X.log.p.value.>= -log(0.15))
IPA_res <- IPA_res %>%  mutate(col=ifelse(z.score >0, "orange", "blue"))
nrow(IPA_res)
#select top 20 if high nrow 
top_20_up <- IPA_res %>% dplyr::slice(1:20)
top_20_down <- IPA_res %>% arrange(z.score) %>% dplyr::slice(1:20)
comb_res <- bind_rows(top_20_up, top_20_down)

pdf(file = "RNA_analysis_naive/Figures/IPA_driver_mutation.pdf", width = 22, height = 10)
ggplot(comb_res, aes(x=z.score, y=fct_reorder(Pathway, z.score, .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  #ggtitle("Top Upregulated Pathways in patients with high vs low PDL1")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 16))
dev.off()






#plot mRNA levels of il-17 

IL17_plot_data <- mycounts %>% t(.) %>%data.frame() %>% select(starts_with("IL17"))

#add to this metadata 
IL17_plot_data <- merge(IL17_plot_data, metadata, by="row.names", all=TRUE)
IL17_plot_data <- IL17_plot_data %>% mutate(dmm_cluster=factor(dmm_cluster))

#plot for all il17 measures 

col_names <- colnames(IL17_plot_data)[2:12]

for (i in col_names){
  p<-  ggplot(IL17_plot_data, aes_string(x = IL17_plot_data$dmm_cluster, y = i)) +
    geom_boxplot(width = 0.2, color = "blue") +
    geom_jitter(width = 0.2,alpha = 0.5,color = "darkblue") +
    scale_y_log10()+
    stat_compare_means(comparisons = list(c("1", "2"), 
                                          c("2","3"), 
                                          c("1", "3")))+
    ggtitle(i)+
    theme_bw()
  
  #save it 
  pdf_output <- paste0(i, paste0("_IL_17_boxplot", paste0("_DMM_clusters",paste0(".pdf"))))
  pdf(file = pdf_output, width = 4, height = 6)
  show(p)
  dev.off()
}


####### plot mRNA of IL-8 (CXCL8)#####

IL8_plot_data <- mycounts %>% t(.) %>%data.frame() %>% select(starts_with("CXCL8"))

#add to this metadata 
IL8_plot_data <- merge(IL8_plot_data, metadata, by="row.names", all=TRUE)
IL8_plot_data <- IL8_plot_data %>% mutate(dmm_cluster=factor(dmm_cluster))

#plot for all IL8_plot_data measures 
p<-  ggplot(IL8_plot_data, aes(x = IL8_plot_data$dmm_cluster, y = CXCL8)) +
  geom_boxplot(width = 0.2, color = "blue") +
  geom_jitter(width = 0.2,alpha = 0.5,color = "darkblue") +
  scale_y_log10()+
  stat_compare_means(comparisons = list(c("1", "2"), 
                                        c("2","3"), 
                                        c("1", "3")))+
  theme_bw()

#save it 
pdf_output <- paste0("comparison", paste0("_IL_8_boxplot", paste0("_DMM_clusters",paste0(".pdf"))))
pdf(file = pdf_output, width = 4, height = 6)
p
dev.off()



############## IPA SCC vs adeno FDR 0.2 ##########

IPA_res <- read_csv("RNA_analysis_naive/Results/IPA/IPA_SCC_vs_adeno_FDR_0.2_plot.csv")
IPA_res$Pathway <- factor(IPA_res$Pathway, levels = unique(IPA_res$Pathway))
IPA_res <- IPA_res %>% arrange(desc(z_score)) %>% dplyr::filter(z_score != 0) %>% filter(z_score!="#NUM!")
IPA_res$neg_log_p_value <- IPA_res$p_value

#leave only pathways with p value < 0.2 (FDR)
IPA_res <- IPA_res %>% filter(neg_log_p_value > -log(0.2))

#set colors of pthways. oragnge is upregulatedm, blue is downregulated 
IPA_res <- IPA_res %>%  mutate(col=ifelse(z_score >0, "orange", "blue"))

pdf(file = "RNA_analysis_naive/Figures/IPA_SCC_vs_adeno_FDR_0.2_plot.pdf", width = 20, height = 12)
ggplot(IPA_res, aes(x=as.numeric(z_score), y=fct_reorder(Pathway,as.numeric(z_score), .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with SCC vs LUAD")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 24, face = "bold"))
dev.off()



############## IPA mutation present or not FDR 0.2 ##########

IPA_res <- read_csv("RNA_analysis_naive/Results/IPA/IPA_Driver_mutations_present_FDR_0.2_plot.csv")
IPA_res$Pathway <- factor(IPA_res$Pathway, levels = unique(IPA_res$Pathway))
IPA_res <- IPA_res %>% arrange(desc(z_score)) %>% dplyr::filter(z_score != 0) %>% filter(z_score!="#NUM!")
IPA_res$neg_log_p_value <- IPA_res$p_value

#leave only pathways with p value < 0.2 (FDR)
IPA_res <- IPA_res %>% filter(neg_log_p_value > -log(0.2))

#set colors of pthways. oragnge is upregulatedm, blue is downregulated 
IPA_res <- IPA_res %>%  mutate(col=ifelse(z_score >0, "orange", "blue"))

pdf(file = "RNA_analysis_naive/Figures/IPA_Driver_mutations_present_FDR_0.2_plot.pdf", width = 20, height = 8)
ggplot(IPA_res, aes(x=as.numeric(z_score), y=fct_reorder(Pathway,as.numeric(z_score), .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with Driver mutation vs not")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 24, face = "bold"))
dev.off()




############## IPA mutation late vs early FDR 0.2 ##########

IPA_res <- read_csv("RNA_analysis_naive/Results/IPA/IPA_late_vs_early_FDR_0.2.csv")
IPA_res$Pathway <- factor(IPA_res$Pathway, levels = unique(IPA_res$Pathway))
IPA_res <- IPA_res %>% arrange(desc(z_score)) %>% dplyr::filter(z_score != 0) %>% filter(z_score!="#NUM!")
IPA_res$neg_log_p_value <- IPA_res$p_value

#leave only pathways with p value < 0.2 (FDR)
IPA_res <- IPA_res %>% filter(neg_log_p_value > -log(0.2))

#set colors of pthways. oragnge is upregulatedm, blue is downregulated 
IPA_res <- IPA_res %>%  mutate(col=ifelse(z_score >0, "orange", "blue"))

pdf(file = "RNA_analysis_naive/Figures/IPA_late_vs_early_FDR_0.2.pdf", width = 20, height = 8)
ggplot(IPA_res, aes(x=as.numeric(z_score), y=fct_reorder(Pathway,as.numeric(z_score), .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with Late vs Early Stage NSCLC")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 24, face = "bold"))
dev.off()


############## IPA mutation dead vs alive 0.2 ##########

IPA_res <- read_csv("RNA_analysis_naive/Results/IPA/IPA_dead_vs_alive_FDR_0.2_plot.csv")
IPA_res$Pathway <- factor(IPA_res$Pathway, levels = unique(IPA_res$Pathway))
IPA_res <- IPA_res %>% arrange(desc(z_score)) %>% dplyr::filter(z_score != 0) %>% filter(z_score!="#NUM!")
IPA_res$neg_log_p_value <- IPA_res$p_value

#leave only pathways with p value < 0.2 (FDR)
IPA_res <- IPA_res %>% filter(neg_log_p_value > -log(0.2))

#set colors of pthways. oragnge is upregulatedm, blue is downregulated 
IPA_res <- IPA_res %>%  mutate(col=ifelse(z_score >0, "orange", "blue"))

pdf(file = "RNA_analysis_naive/Figures/IPA_dead_vs_alive_FDR_0.2_plot.pdf", width = 20, height = 12)
ggplot(IPA_res, aes(x=as.numeric(z_score), y=fct_reorder(Pathway,as.numeric(z_score), .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with Dead vs Aalive Patients")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 24, face = "bold"))
dev.off()



############## IPA mutation pneumonitis vs not 0.2 ##########

IPA_res <- read_csv("RNA_analysis_naive/Results/IPA/IPA_irAE_pneumonitis_FDR_0.2_plot.csv")
IPA_res$Pathway <- factor(IPA_res$Pathway, levels = unique(IPA_res$Pathway))
IPA_res <- IPA_res %>% arrange(desc(z_score)) %>% dplyr::filter(z_score != 0) %>% filter(z_score!="#NUM!")
IPA_res$neg_log_p_value <- IPA_res$p_value

#leave only pathways with p value < 0.2 (FDR)
IPA_res <- IPA_res %>% filter(neg_log_p_value > -log(0.2))

#set colors of pthways. oragnge is upregulatedm, blue is downregulated 
IPA_res <- IPA_res %>%  mutate(col=ifelse(z_score >0, "orange", "blue"))

pdf(file = "RNA_analysis_naive/Figures/IPA_irAE_pneumonitis_FDR_0.2_plot.pdf", width = 20, height = 22)
ggplot(IPA_res, aes(x=as.numeric(z_score), y=fct_reorder(Pathway,as.numeric(z_score), .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with Pneumonitis vs Without Pneumonitis")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 24, face = "bold"))
dev.off()







#######pd l1 levels - 2 levels IPA ###

IPA_res <- read_csv("RNA_analysis_naive/Results/IPA/IPA_PD_L1_expression_2_lev_FDR_0.2_plot.csv")
IPA_res$Pathway <- factor(IPA_res$Pathway, levels = unique(IPA_res$Pathway))
IPA_res <- IPA_res %>% arrange(desc(z_score)) %>% dplyr::filter(z_score != 0) %>% filter(z_score!="#NUM!")
IPA_res$neg_log_p_value <- IPA_res$p_value

#leave only pathways with p value < 0.2 (FDR)
IPA_res <- IPA_res %>% filter(neg_log_p_value > -log(0.2))

#set colors of pthways. oragnge is upregulatedm, blue is downregulated 
IPA_res <- IPA_res %>%  mutate(col=ifelse(z_score >0, "orange", "blue"))

pdf(file = "RNA_analysis_naive/Figures/IPA_PD_L1_expression_2_lev_FDR_0.2_plot.pdf", width = 20, height = 8)
ggplot(IPA_res, aes(x=as.numeric(z_score), y=fct_reorder(Pathway,as.numeric(z_score), .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue", "orange"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with High vs Low PD-L1")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 24, face = "bold"))
dev.off()




######## IPA heatmaps with FDR 0.2 ##### 

#PDL-1 three levels 
#read IPA results 
IPA_res <- read.csv(file = "RNA_analysis_naive/Results/IPA/IPA_PD_L1_expression_comparison_FDR_0.2_plot.csv")

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
max(IPA_res$Low_vs_Negative, na.rm = TRUE)
min(IPA_res$Low_vs_Negative,na.rm = TRUE)
max(IPA_res$High_vs_Low,na.rm = TRUE)
min(IPA_res$High_vs_Low, na.rm = TRUE)
max(IPA_res$High_vs_Negative, na.rm = TRUE)
min(IPA_res$High_vs_Negative, na.rm = TRUE)

col_fun = colorRamp2(c(-3.138, 0, 3.207), c("blue", "white", "orange"))
col_fun(seq(-3.138, 0, 3.207))

IPA_res <- IPA_res %>% arrange(desc(High_vs_Low))


#convert data to matrix 
IPA_res_mat <- IPA_res


#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "RNA_analysis_naive/Figures/IPA_PD_L1_expression_comparison_FDR_0.2_plot.pdf", height = 20, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))
dev.off()














###### random forest on RNA data ##### 
#use normalized counts table 
counts <- mycounts
counts=counts[apply(counts, 1, function(x) sum(x>=300)>=3),]
counts=apply(counts, 2, function(x) return(log(x+1)))

#run the model. change top_n every time you want to examine top taxa
set.seed(1234)

#get metadata 
metadata.for.RF <- metadata_naive %>% 
  data.frame() %>% 
  mutate(one_year_mort = ifelse(death=="1" & OS_days > 365, "greater_12", 
                                ifelse(death=="1" & OS_days<365, "less_12", 
                                       ifelse(death=="0" & OS_days>365, "greater_12", "NA")))) %>% 
  filter(one_year_mort!="NA") %>% 
  mutate(outcome=one_year_mort) %>% 
  mutate(outcome = ifelse(outcome=="greater_12", "0", "1"))


#get features (taxa) 
features <- counts %>% data.frame(.)
colnames(features) <- gsub("X", "", colnames(features))
features <- features %>% select(rownames(metadata.for.RF))

#define outcome 
outcome <- as.factor(metadata.for.RF$one_year_mort)
outcome <- ifelse(outcome=="greater_12", "0","1") 


IVs=list(); AUCs=list()

for (i in 1:20){
  ### stratify 5-fold CV (80% discovery and 20 validation)
  train_index <- c(sample(which(metadata.for.RF$outcome==0), round(length(which(metadata.for.RF$outcome==0))*0.8)),
                   sample(which(metadata.for.RF$outcome!=0), round(length(which(metadata.for.RF$outcome!=0))*0.8)))
  train_data <- features[,train_index]
  meta.training <- metadata.for.RF[train_index,]
  #test_data <- rf.data_complete[-train_index, ]
  
  #tune 
  mtry <- sqrt(ncol(train_data))
  tunegrid <- expand.grid(.mtry=mtry)
  metric <- "ROC"
  control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)
  
  #run model 
  fit <- randomForest(outcome~., data = data.frame(outcome=as.factor(meta.training$outcome), 
                                                   t(train_data)), 
                      trControl=control, metric=metric, 
                      ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)
  
  #get variable importance
  # Extract the mean decrease impurity (MDI) values for each feature
  mdi <- randomForest::importance(fit) %>% as.data.frame()
  
  # Normalize the MDI values 
  Norm.Gini= mdi$MeanDecreaseGini/sum(mdi$MeanDecreaseGini)
  Norm.Gini=Norm.Gini/max(Norm.Gini)
  mdi=data.frame(taxa=rownames(train_data),Norm.Gini)
  
  #add enrich group (or outcome) -->
  # features.means <- by(t(Genus.training), meta.training$delta, colMeans) -->
  # features.means <- do.call(cbind, features.means) -->
  # idx_enrich <- apply(features.means, 1, which.max) -->
  # group_enrich <- colnames(features.means)[idx_enrich] -->
  # -->
  # mdi$enrich_group <- group_enrich -->
  mdi=mdi[order(mdi$Norm.Gini, decreasing = TRUE),]
  
  #cutoffs=c(seq(0.01, 0.09,0.01), seq(0.1, 1,0.1))
  
  cutoffs=c(seq(0.01, 0.09,0.01), seq(0.1, 1,0.1))
  
  VIs.taxon=matrix(NA,nrow = nrow(mdi), ncol=length(cutoffs))
  rownames(VIs.taxon)=mdi$taxa
  
  AUCs.taxon=NULL
  
  for(tt in 1:length(cutoffs)) {
    
    taxa.top=mdi$taxa[1:round(nrow(mdi)*cutoffs[tt])]
    
    fit.per <- randomForest(outcome~., data = data.frame(outcome=as.factor(meta.training$outcome), 
                                                         t(train_data[taxa.top,])),
                            ntree=500, importance=TRUE)
    
    #get variable importance 
    # Extract the mean decrease impurity (MDI) values for each feature
    mdi.per <- randomForest::importance(fit.per) %>% as.data.frame()
    
    # Normalize the MDI values
    Norm.Gini= mdi.per$MeanDecreaseGini/sum(mdi.per$MeanDecreaseGini)
    Norm.Gini=Norm.Gini/max(Norm.Gini)
    
    mdi.per=data.frame(taxa=taxa.top,Norm.Gini)
    
    #add enrich group (or outcome)
    # features.means <- by(t(Genus.training[ASV.top,]), meta.training$delta, colMeans) 
    # features.means <- do.call(cbind, features.means)
    # idx_enrich <- apply(features.means, 1, which.max)
    # group_enrich <- colnames(features.means)[idx_enrich] 
    # mdi.per$enrich_group <- group_enrich
    mdi.per=mdi.per[order(mdi.per$Norm.Gini, decreasing = TRUE),] 
    VIs.taxon[mdi.per$taxa,tt]=mdi.per$Norm.Gini
    
    #get predictions
    predictions <- predict(fit.per, newdata = data.frame(outcome=as.factor(metadata.for.RF[-train_index,]$outcome),
                                                         t(features[taxa.top,-train_index])), type = "prob")
    predictions <- data.frame(predictions)
    roc <- roc(metadata.for.RF[-train_index,]$outcome, predictions$X1)
    
    AUCs.taxon=c(AUCs.taxon, roc$auc) 
  }
  IVs=c(IVs, list(VIs.taxon))
  
  AUCs=c(AUCs, list(AUCs.taxon)) 
}

aa=list(IVs, AUCs)


### plot AUC values 
AUCs=matrix(unlist(aa[[2]]),nrow=20, byrow=TRUE)
AUCs=apply(AUCs, 2, function(x) return(ifelse(x<0.5, 1-x, x)))


AUCs_taxon=AUCs[seq(1,20,2),];cutoffs=c(seq(0.01, 0.09,0.01), seq(0.1, 1,0.1))*100

AUCs_taxon=melt(data.frame(cutoffs, t(AUCs_taxon)), id.vars=1)[,-2]

#plot
pdf(file = "RNA_analysis_naive/Figures/AUC_gene_level_RF.pdf", height = 6, width =  )
ggline(AUCs_taxon, x = "cutoffs", y = "value",  position = position_dodge(0.8), color = "dodgerblue",
       add = c("mean_sd"),numeric.x.axis = TRUE)+theme_bw(base_size = 14)+
  ylab("AUC")+xlab("% Top Gene")+
  scale_x_continuous(breaks = seq(0,100, by = 10))+
  #scale_color_manual(values = c("red", "blue")) +  # Adjust line colors
  theme(legend.position = "top", legend.title = element_blank())+
  ggtitle("AUC  at gene level (repetition=20)")
dev.off()


# Mean on all columns
df2 <- AUCs_taxon %>% group_by(cutoffs) %>%
  summarise(across(everything(), mean),
            .groups = 'drop')  %>%
  as.data.frame()

#get best cutoff for AUC 
df2$cutoffs[which.max(df2$value)] #80




#run the model for the best AUC 
VIs.taxon=matrix(NA,nrow = nrow(features), ncol=20)
rownames(VIs.taxon)=rownames(features)

AUCs=NULL

for (i in 1:20){
  ### stratify 5-fold CV (80% discovery and 20 validation)
  train_index <- c(sample(which(metadata.for.RF$outcome==0), round(length(which(metadata.for.RF$outcome==0))*0.8)),
                   sample(which(metadata.for.RF$outcome!=0), round(length(which(metadata.for.RF$outcome!=0))*0.8)))
  train_data <- features[,train_index]
  meta.training <- metadata.for.RF[train_index,]
  #test_data <- rf.data_complete[-train_index, ]
  
  #tune 
  mtry <- sqrt(ncol(train_data))
  tunegrid <- expand.grid(.mtry=mtry)
  metric <- "ROC"
  control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)
  
  #run model 
  fit <- randomForest(outcome~., data = data.frame(outcome=as.factor(meta.training$outcome), 
                                                   t(train_data)), 
                      trControl=control, metric=metric, 
                      ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)
  
  #get variable importance
  # Extract the mean decrease impurity (MDI) values for each feature
  mdi <- randomForest::importance(fit) %>% as.data.frame()
  
  # Normalize the MDI values 
  Norm.Gini= mdi$MeanDecreaseGini/sum(mdi$MeanDecreaseGini)
  Norm.Gini=Norm.Gini/max(Norm.Gini)
  mdi=data.frame(taxa=rownames(train_data),Norm.Gini)
  mdi=mdi[order(mdi$Norm.Gini, decreasing = TRUE),]
  
  ## top 30%
  taxa.top <- mdi$taxa[1:round(nrow(mdi)*0.8)]
  
  fit.per <- randomForest(outcome~., data = data.frame(outcome=as.factor(meta.training$outcome), 
                                                       t(train_data[taxa.top,])),
                          ntree=500, importance=TRUE)
  
  #get variable importance 
  # Extract the mean decrease impurity (MDI) values for each feature
  mdi.per <- randomForest::importance(fit.per) %>% as.data.frame()
  
  # Normalize the MDI values
  Norm.Gini= mdi.per$MeanDecreaseGini/sum(mdi.per$MeanDecreaseGini)
  Norm.Gini=Norm.Gini/max(Norm.Gini)
  
  mdi.per=data.frame(taxa=taxa.top,Norm.Gini)
  
  mdi.per=mdi.per[order(mdi.per$Norm.Gini, decreasing = TRUE),] 
  VIs.taxon[mdi.per$taxa,i]=mdi.per$Norm.Gini
  
  #get predictions
  predictions.taxon <- predict(fit.per, newdata = data.frame(outcome=as.factor(metadata.for.RF[-train_index,]$outcome),
                                                             t(features[taxa.top,-train_index])), type = "prob")
  predictions.taxon <- data.frame(predictions.taxon)
  roc <- roc(metadata.for.RF[-train_index,]$outcome, predictions.taxon$X1)
  
  AUCs=c(AUCs, roc$auc) 
}

aa=list(VIs.taxon, AUCs)


AUCs=aa[[2]]
AUCs <- data.frame(AUCs)
AUCs=data.frame(ID=1:nrow(AUCs), AUCs, check.names = FALSE)
AUCs=melt(AUCs, id.vars=1)

pdf(file = "RNA_analysis_naive/Figures/AUC_best_cutoff_0.8_metatrans_RF.pdf", height = 6, width = 3)
ggboxplot(AUCs, x = "variable", y = "value", color="variable",
          add = "jitter")+theme_bw(base_size = 14)+
  ylab("AUC")+xlab("")+ggtitle("AUCs in the host RNA dataset across 20 repetitions")+
  theme(legend.position = "none")
dev.off()


VI_taxon=aa[[1]]
VI_taxon=VI_taxon[apply(VI_taxon, 1, function(x) sum(!is.na(x)))>0,]
VI_taxon=apply(VI_taxon, 2, function(x) return(ifelse(is.na(x), 0, x)))

VI_taxon=cbind(apply(VI_taxon,1, mean),apply(VI_taxon,1, sd))

colnames(VI_taxon)=c("Average_Gini", "SD")
VI_taxon<- VI_taxon %>% 
  data.frame() %>% 
  arrange(desc(Average_Gini)) %>% 
  mutate(gene=rownames(.))

write.csv(VI_taxon, file="RNA_analysis_naive/Results/RF_variable_importance_best_AUC_genes.csv")


VI_taxon %>% filter(Average_Gini >0)

















######### Add RECIST score data ###########


#read recist data 
recist_data <- read.csv(file = "RECIST_data_upd.csv")
recist_data$recist_categ_2_lev <- ifelse(recist_data$ultimate_response_status=="PD", "PD", "non_PD")
head(recist_data)

#get patiens of the cohort 
dds_naive$Subject_ID1
recist_data_bal <- recist_data %>% filter(upd_subj_ID%in% dds_naive$Subject_ID1)

nrow(recist_data_bal)
ncol(dds_naive)

#add the data to physeq table 
temp <- recist_data_bal
temp <- temp %>% dplyr::select(upd_subj_ID, Subject_ID1, ultimate_response_status, recist_categ_2_lev)
temp$Subject_ID1 <- temp$upd_subj_ID
temp <- temp %>% dplyr::select(-upd_subj_ID)
#now order by sample data 
temp_2 <- data.frame(dds_naive@colData)

temp <- temp[match(temp_2$Subject_ID1, temp$Subject_ID1), ]

# View the ordered data frame
head(temp)
head(temp_2)

### add recist score data to the table 
dds_naive$ultimate_response_status <- temp$ultimate_response_status
dds_naive$recist_categ_2_lev <- temp$recist_categ_2_lev

#same for relative table 
vsd_naive$ultimate_response_status <- temp$ultimate_response_status
vsd_naive$recist_categ_2_lev <- temp$recist_categ_2_lev

#same for metadata
metadata_naive$ultimate_response_status <- temp$ultimate_response_status
metadata_naive$recist_categ_2_lev <- temp$recist_categ_2_lev

dds.analysis <- dds_naive
vsd.analysis <- vsd_naive
#choose  variable 

#choose table to use
dds.analysis <- dds.analysis[,!is.na(dds.analysis$recist_categ_2_lev)]
vsd.analysis <- vsd.analysis[,!is.na(vsd.analysis$recist_categ_2_lev)]

#choose variable 
v= "recist_categ_2_lev"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "non_PD")

#change experiment design if needed
design(dds.analysis) <- ~ recist_categ_2_lev

#create appropraite metadata for edgeR
metadata_RNA_analysis <- metadata_naive %>% filter(!is.na(recist_categ_2_lev))

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "red"
col2 <- "blue"

######################Plot PCoA
colnames(dds.analysis)
#define variable 
dds.analysis[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$recist_categ_2_lev
#newResults$v <- factor(newResults$v, levels = c("greater_24", "less_24"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("RNA_analysis_naive/Results/Beta.Diversity.Bray.", paste0(v, paste0("_less_vs_greater_24.txt")))
            , sep = "\t", row.names = T)


pdf(file = paste0("RNA_analysis_naive/Figures/Beta.Diversity.Bray_", paste0(v, paste0("_less_vs_greater_24.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Non PD", "PD")), size=10) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis histology  using edgeR ###################
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("RNA_analysis_naive/Results/edgeR.results_", paste0(v, paste0("_PD_vs_non_PD.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("RNA_analysis_naive/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_PD_vs_non_PD.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("RNA_analysis_naive/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_PD_vs_non_PD.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()



##########plot IPA analysis of RECIST PD vs non PD #####
IPA_res <- read.csv(file="RNA_analysis_naive/Results/IPA/recist_PD_vs_non_PD1_filtered.csv")

#IPA_res <- IPA_res %>% filter(!is.na(z_score))

plot_dat <- IPA_res

# Order the rows by z score 
plot_dat <- plot_dat %>% arrange(desc(z_score))

#select top 40 upregulated pathways 
top_20_up <- plot_dat %>% dplyr::slice(1:40)
#select downregulated pathways 
top_20_down <- plot_dat %>% filter(z_score<0)

plot_dat <- bind_rows(top_20_up, top_20_down) 

# Define colors for each levels of qualitative variables
library(circlize)
max(plot_dat$z_score, na.rm = TRUE)
min(plot_dat$z_score,na.rm = TRUE)

col_fun = colorRamp2(c(-1.387, 0, 4.964), c("blue", "white", "orange"))
col_fun(seq(-1.387, 0, 4.964))

#plot_dat <- plot_dat %>% arrange(desc(z_score))


#convert data to matrix 
IPA_res_mat <- plot_dat


#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)
IPA_res_mat<- IPA_res_mat[,-c(1:2)]

library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
pdf(file = "RNA_analysis_naive/Figures/IPA_RECIST_score_PD_vs_all_top_pathways.pdf", height = 14, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =plot_dat$Pathway, 
                        width = unit(2, "cm"))
dev.off()
















############ IPA analysis for naive patients according to diffrent clinical features ########

####### histology 

IPA_res <- read.csv(file = "RNA_analysis_naive/Results/IPA/histology_scc_vs_adeno.csv")
IPA_res$p_value <- 10^(-IPA_res$minus_log_10_p_value)
IPA_res$adj_p <- p.adjust(IPA_res$p_value, method = "BH")
IPA_res <- IPA_res %>% filter(z_score!="#NUM!")


#get top 20 upregulated pathways 
top_20 <- IPA_res %>% filter(z_score>0)
top_20 <- top_20 %>% arrange(desc(z_score))
top_20 <- top_20 %>% filter(minus_log_10_p_value >= -log10(0.05))

#slice top 25 pathways 
top_20 <- top_20 %>% dplyr::slice(1:25)

#top 20 downregulated 
down_20 <- IPA_res %>% filter(z_score<0)
down_20 <- down_20 %>% arrange(desc(z_score))
down_20 <- down_20 %>% filter(minus_log_10_p_value >= -log10(0.05))

down_20 <- down_20 %>% dplyr::slice(1:25)

#combine 
ipa_res_for_plot <- rbind(top_20, down_20)
ipa_res_for_plot$z_score <- as.numeric(ipa_res_for_plot$z_score)
ipa_res_for_plot$minus_log_10_p_value <- as.numeric(ipa_res_for_plot$minus_log_10_p_value)
ipa_res_for_plot$regulation <- ifelse(ipa_res_for_plot$z_score > 0, "upregulated", "downregulated")

plot_dat <- ipa_res_for_plot

# Define colors for each levels of qualitative variables
library(circlize)
max(plot_dat$z_score, na.rm = TRUE)
min(plot_dat$z_score,na.rm = TRUE)

col_fun = colorRamp2(c(-2.111, 0, 2.828), c("blue", "white", "orange"))
col_fun(seq(-2.111, 0, 2.828))

plot_dat <- plot_dat %>% arrange(desc(z_score))

#convert data to matrix 
IPA_res_mat <- plot_dat

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat %>% dplyr::select(z_score)

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)
#IPA_res_mat<- IPA_res_mat[,-c(1:2)]

library(ComplexHeatmap)
######ploting
pdf(file = "RNA_analysis_naive/Figures/IPA_histology_SCC_vs_adeno_top_20_up_top_20_down.pdf", height = 14, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =plot_dat$Pathway, 
                        width = unit(2, "cm"))
dev.off()








#######pneumonitis 


IPA_res <- read.csv(file = "RNA_analysis_naive/Results/IPA/pneumonitis_yes_vs_no.csv")
IPA_res$p_value <- 10^(-IPA_res$minus_log_10_p_value)
IPA_res$adj_p <- p.adjust(IPA_res$p_value, method = "BH")
IPA_res <- IPA_res %>% filter(z_score!="#NUM!")


#get top 20 upregulated pathways 
top_20 <- IPA_res %>% filter(z_score>0)
top_20 <- top_20 %>% arrange(desc(z_score))
top_20 <- top_20 %>% filter(minus_log_10_p_value >= -log10(0.05))

#slice top 25 pathways 
top_20 <- top_20 %>% dplyr::slice(1:25)

#top 20 downregulated 
down_20 <- IPA_res %>% filter(z_score<0)
down_20 <- down_20 %>% arrange(desc(z_score))
down_20 <- down_20 %>% filter(minus_log_10_p_value >= -log10(0.05))

down_20 <- down_20 %>% dplyr::slice(1:25)

#combine 
ipa_res_for_plot <- rbind(top_20, down_20)
ipa_res_for_plot$z_score <- as.numeric(ipa_res_for_plot$z_score)
ipa_res_for_plot$minus_log_10_p_value <- as.numeric(ipa_res_for_plot$minus_log_10_p_value)
ipa_res_for_plot$regulation <- ifelse(ipa_res_for_plot$z_score > 0, "upregulated", "downregulated")

plot_dat <- ipa_res_for_plot

# Define colors for each levels of qualitative variables
library(circlize)
max(plot_dat$z_score, na.rm = TRUE)
min(plot_dat$z_score,na.rm = TRUE)

col_fun = colorRamp2(c(-2.111, 0, 2.828), c("blue", "white", "orange"))
col_fun(seq(-2.111, 0, 2.828))

plot_dat <- plot_dat %>% arrange(desc(z_score))

#convert data to matrix 
IPA_res_mat <- plot_dat

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat %>% dplyr::select(z_score)

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)
#IPA_res_mat<- IPA_res_mat[,-c(1:2)]

library(ComplexHeatmap)
######ploting
pdf(file = "RNA_analysis_naive/Figures/IPA_pneumonitis_yes_vs_no_top_20_up_top_20_down.pdf", height = 14, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =plot_dat$Pathway, 
                        width = unit(2, "cm"))
dev.off()















######## driver mutation 

IPA_res <- read.csv(file = "RNA_analysis_naive/Results/IPA/driver_mutation_yes_vs_no.csv")
IPA_res$p_value <- 10^(-IPA_res$minus_log_10_p_value)
IPA_res$adj_p <- p.adjust(IPA_res$p_value, method = "BH")
IPA_res <- IPA_res %>% filter(z_score!="#NUM!")


#get top 20 upregulated pathways 
top_20 <- IPA_res %>% filter(z_score>0)
top_20 <- top_20 %>% arrange(desc(z_score))
top_20 <- top_20 %>% filter(minus_log_10_p_value >= -log10(0.05))

#slice top 25 pathways 
top_20 <- top_20 %>% dplyr::slice(1:25)

#top 20 downregulated 
down_20 <- IPA_res %>% filter(z_score<0)
down_20 <- down_20 %>% arrange(desc(z_score))
down_20 <- down_20 %>% filter(minus_log_10_p_value >= -log10(0.05))

down_20 <- down_20 %>% dplyr::slice(1:25)

#combine 
ipa_res_for_plot <- rbind(top_20, down_20)
ipa_res_for_plot$z_score <- as.numeric(ipa_res_for_plot$z_score)
ipa_res_for_plot$minus_log_10_p_value <- as.numeric(ipa_res_for_plot$minus_log_10_p_value)
ipa_res_for_plot$regulation <- ifelse(ipa_res_for_plot$z_score > 0, "upregulated", "downregulated")

plot_dat <- ipa_res_for_plot

# Define colors for each levels of qualitative variables
library(circlize)
max(plot_dat$z_score, na.rm = TRUE)
min(plot_dat$z_score,na.rm = TRUE)

col_fun = colorRamp2(c(-2.111, 0, 2.828), c("blue", "white", "orange"))
col_fun(seq(-2.111, 0, 2.828))

plot_dat <- plot_dat %>% arrange(desc(z_score))

#convert data to matrix 
IPA_res_mat <- plot_dat

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat %>% dplyr::select(z_score)

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)
#IPA_res_mat<- IPA_res_mat[,-c(1:2)]

library(ComplexHeatmap)
######ploting
pdf(file = "RNA_analysis_naive/Figures/IPA_driver_mutation_yes_vs_no_top_20_up_top_20_down.pdf", height = 14, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =plot_dat$Pathway, 
                        width = unit(2, "cm"))
dev.off()







##### stage 


IPA_res <- read.csv(file = "RNA_analysis_naive/Results/IPA/stage_late_vs_early.csv")
IPA_res$p_value <- 10^(-IPA_res$minus_log_10_p_value)
IPA_res$adj_p <- p.adjust(IPA_res$p_value, method = "BH")
IPA_res <- IPA_res %>% filter(z_score!="#NUM!")


#get top 20 upregulated pathways 
top_20 <- IPA_res %>% filter(z_score>0)
top_20 <- top_20 %>% arrange(desc(z_score))
top_20 <- top_20 %>% filter(minus_log_10_p_value >= -log10(0.05))

#slice top 25 pathways 
top_20 <- top_20 %>% dplyr::slice(1:25)

#top 20 downregulated 
down_20 <- IPA_res %>% filter(z_score<0)
down_20 <- down_20 %>% arrange(desc(z_score))
down_20 <- down_20 %>% filter(minus_log_10_p_value >= -log10(0.05))

down_20 <- down_20 %>% dplyr::slice(1:25)

#combine 
ipa_res_for_plot <- rbind(top_20, down_20)
ipa_res_for_plot$z_score <- as.numeric(ipa_res_for_plot$z_score)
ipa_res_for_plot$minus_log_10_p_value <- as.numeric(ipa_res_for_plot$minus_log_10_p_value)
ipa_res_for_plot$regulation <- ifelse(ipa_res_for_plot$z_score > 0, "upregulated", "downregulated")

plot_dat <- ipa_res_for_plot

# Define colors for each levels of qualitative variables
library(circlize)
max(plot_dat$z_score, na.rm = TRUE)
min(plot_dat$z_score,na.rm = TRUE)

col_fun = colorRamp2(c(-2.111, 0, 2.828), c("blue", "white", "orange"))
col_fun(seq(-2.111, 0, 2.828))

plot_dat <- plot_dat %>% arrange(desc(z_score))

#convert data to matrix 
IPA_res_mat <- plot_dat

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat %>% dplyr::select(z_score)

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)
#IPA_res_mat<- IPA_res_mat[,-c(1:2)]

library(ComplexHeatmap)
######ploting
pdf(file = "RNA_analysis_naive/Figures/IPA_stage_late_vs_early_top_20_up_top_20_down.pdf", height = 14, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =plot_dat$Pathway, 
                        width = unit(2, "cm"))
dev.off()







##### one year mortality 





IPA_res <- read.csv(file = "RNA_analysis_naive/Results/IPA/one_year_mort_less_vs_greater.csv")
IPA_res$p_value <- 10^(-IPA_res$minus_log_10_p_value)
IPA_res$adj_p <- p.adjust(IPA_res$p_value, method = "BH")
IPA_res <- IPA_res %>% filter(z_score!="#NUM!")


#get top 20 upregulated pathways 
top_20 <- IPA_res %>% filter(z_score>0)
top_20 <- top_20 %>% arrange(desc(z_score))
top_20 <- top_20 %>% filter(minus_log_10_p_value >= -log10(0.05))

#slice top 25 pathways 
top_20 <- top_20 %>% dplyr::slice(1:25)

#top 20 downregulated 
down_20 <- IPA_res %>% filter(z_score<0)
down_20 <- down_20 %>% arrange(desc(z_score))
down_20 <- down_20 %>% filter(minus_log_10_p_value >= -log10(0.05))

down_20 <- down_20 %>% dplyr::slice(1:25)

#combine 
ipa_res_for_plot <- rbind(top_20, down_20)
ipa_res_for_plot$z_score <- as.numeric(ipa_res_for_plot$z_score)
ipa_res_for_plot$minus_log_10_p_value <- as.numeric(ipa_res_for_plot$minus_log_10_p_value)
ipa_res_for_plot$regulation <- ifelse(ipa_res_for_plot$z_score > 0, "upregulated", "downregulated")

plot_dat <- ipa_res_for_plot

# Define colors for each levels of qualitative variables
library(circlize)
max(plot_dat$z_score, na.rm = TRUE)
min(plot_dat$z_score,na.rm = TRUE)

col_fun = colorRamp2(c(-2.111, 0, 2.828), c("blue", "white", "orange"))
col_fun(seq(-2.111, 0, 2.828))

plot_dat <- plot_dat %>% arrange(desc(z_score))

#convert data to matrix 
IPA_res_mat <- plot_dat

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat %>% dplyr::select(z_score)

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)
#IPA_res_mat<- IPA_res_mat[,-c(1:2)]

library(ComplexHeatmap)
######ploting
pdf(file = "RNA_analysis_naive/Figures/IPA_one_year_mort_less_vs_greater_top_20_up_top_20_down.pdf", height = 14, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =plot_dat$Pathway, 
                        width = unit(2, "cm"))
dev.off()







#PD_1_2_lev_high_vs_low

IPA_res <- read.csv(file = "RNA_analysis_naive/Results/IPA/PD_1_2_lev_high_vs_low.csv")
IPA_res$p_value <- 10^(-IPA_res$minus_log_10_p_value)
IPA_res$adj_p <- p.adjust(IPA_res$p_value, method = "BH")
IPA_res <- IPA_res %>% filter(z_score!="#NUM!")


#get top 20 upregulated pathways 
top_20 <- IPA_res %>% filter(z_score>0)
top_20 <- top_20 %>% arrange(desc(z_score))
top_20 <- top_20 %>% filter(minus_log_10_p_value >= -log10(0.05))

#slice top 25 pathways 
top_20 <- top_20 %>% dplyr::slice(1:25)

#top 20 downregulated 
down_20 <- IPA_res %>% filter(z_score<0)
down_20 <- down_20 %>% arrange(desc(z_score))
down_20 <- down_20 %>% filter(minus_log_10_p_value >= -log10(0.05))

down_20 <- down_20 %>% dplyr::slice(1:25)

#combine 
ipa_res_for_plot <- rbind(top_20, down_20)
ipa_res_for_plot$z_score <- as.numeric(ipa_res_for_plot$z_score)
ipa_res_for_plot$minus_log_10_p_value <- as.numeric(ipa_res_for_plot$minus_log_10_p_value)
ipa_res_for_plot$regulation <- ifelse(ipa_res_for_plot$z_score > 0, "upregulated", "downregulated")

plot_dat <- ipa_res_for_plot

# Define colors for each levels of qualitative variables
library(circlize)
max(plot_dat$z_score, na.rm = TRUE)
min(plot_dat$z_score,na.rm = TRUE)

col_fun = colorRamp2(c(-2.111, 0, 2.828), c("blue", "white", "orange"))
col_fun(seq(-2.111, 0, 2.828))

plot_dat <- plot_dat %>% arrange(desc(z_score))

#convert data to matrix 
IPA_res_mat <- plot_dat

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat %>% dplyr::select(z_score)

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)
#IPA_res_mat<- IPA_res_mat[,-c(1:2)]

library(ComplexHeatmap)
######ploting
pdf(file = "RNA_analysis_naive/Figures/IPA_PD_1_2_lev_high_vs_low_top_20_up_top_20_down.pdf", height = 14, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =plot_dat$Pathway, 
                        width = unit(2, "cm"))
dev.off()







####### k means clustering naive, updated, MG ######## 

IPA_res <- read.csv(file = "RNA_analysis_naive/Results/IPA/k_means_clusters_MG_updated_comparison_selected_pathways.csv")


# Define colors for each levels of qualitative variables
library(circlize)
max(IPA_res$SPT_vs_BPT, na.rm = TRUE)
min(IPA_res$SPT_vs_BPT,na.rm = TRUE)
max(IPA_res$BAL_vs_BPT, na.rm = TRUE)
min(IPA_res$BAL_vs_BPT,na.rm = TRUE)
max(IPA_res$BAL_vs_SPT, na.rm = TRUE)
min(IPA_res$BAL_vs_SPT,na.rm = TRUE)

col_fun = colorRamp2(c(-3.71, 0, 5.653), c("blue", "white", "orange"))
col_fun(seq(-3.71, 0, 5.653))

IPA_res <- IPA_res %>% arrange(desc(SPT_vs_BPT))

#convert data to matrix 
IPA_res_mat <- IPA_res

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

IPA_res_mat<- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
pdf(file = "RNA_analysis_naive/Figures/IPA_k_means_clusters_upd_MG.pdf", height = 14, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))
dev.off()







####### k means clustering naive, updated, MTS ######## 

IPA_res <- read.csv(file = "RNA_analysis_naive/Results/IPA/k_means_clusters_MTS_updated_comparison_selected_pathways.csv")


# Define colors for each levels of qualitative variables
library(circlize)
max(IPA_res$SPT_vs_BPT, na.rm = TRUE)
min(IPA_res$SPT_vs_BPT,na.rm = TRUE)
max(IPA_res$BAL_vs_BPT, na.rm = TRUE)
min(IPA_res$BAL_vs_BPT,na.rm = TRUE)
max(IPA_res$BAL_vs_SPT, na.rm = TRUE)
min(IPA_res$BAL_vs_SPT,na.rm = TRUE)

col_fun = colorRamp2(c(-4.667, 0, 6.621), c("blue", "white", "orange"))
col_fun(seq(-4.667, 0, 6.621))

IPA_res <- IPA_res %>% arrange(desc(SPT_vs_BPT))

#convert data to matrix 
IPA_res_mat <- IPA_res

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

IPA_res_mat<- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
pdf(file = "RNA_analysis_naive/Figures/IPA_k_means_clusters_upd_MTS.pdf", height = 16, width = 18)
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))
dev.off()






save.image("immunoth_host_RNA_NEW.RData")

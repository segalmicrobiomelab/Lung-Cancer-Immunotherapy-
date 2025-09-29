######## host transcriptome analysis of lower airway microbiome 
######## paper: immunotherapy lung cancer paper 
######## Fares Darawshy MD
######## Sep 2025 


### set your wording directory  
setwd("")


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
########### read files and prepare them ###############

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

# get K means 3 clusters data from the MT analysis  
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


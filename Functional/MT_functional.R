
### set wd 
setwd("~/Dropbox (NYU Langone Health)/Fares Darawshy’s files/Home/Projects/Immunotherapy_BAL/final_analysis")

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


#######define edgeR function #########
phyloseq_to_edgeR = function(physeq, group, method = "RLE", ...) {
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if (identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1) {
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts = x, group = group, genes = taxonomy, remove.zeros = TRUE, 
              ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method = method)
  # Check for division by zero inside `calcNormFactors`
  if (!all(is.finite(z$samples$norm.factors))) {
    stop("Something wrong with edgeR::calcNormFactors on this data,\n non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}



################metabolic pathways ##############################

#read metadata 
metadata <- read.delim2(file = "immuntherapy_metadata_final_cohort.txt")

#make genes names rownames 
names <- make.unique(metadata$X)
rownames(metadata) <- names

#read counts table 
mycounts <- read.delim2(file = "metatranscriptomics/Data/upd_tables/counts.ko.pathway_metabolism.fmap.MT.tsv")

#make genes names rownames 
names <- make.unique(mycounts$X)
rownames(mycounts) <- names

#remove X column 
mycounts <- mycounts %>% dplyr::select(-X)
metadata <- metadata %>% dplyr::select(-X)

rownames(metadata) <- gsub("_", ".", rownames(metadata))

#mycounts is arranged by lanes, so we are going to arrange it, transpose and sum all reads per one sample 
temp_mycounts <- mycounts
colnames(temp_mycounts) <- gsub("X", "", colnames(temp_mycounts))
colnames(temp_mycounts) <- gsub("_S\\d+_L\\d+", "", colnames(temp_mycounts))
colnames(temp_mycounts) <- gsub("_", ".", colnames(temp_mycounts))
#transpose mycounts
t_temp_mycounts <- as.data.frame(t(temp_mycounts))
# create a column called sample_name and clean it up so you can group by each sample later 
t_temp_mycounts$sample_name <- rownames(t_temp_mycounts)
t_temp_mycounts$sample_name <- gsub("X", "", t_temp_mycounts$sample_name)
t_temp_mycounts$sample_name <- gsub("\\.\\d+$", "", t_temp_mycounts$sample_name)
#sum all reads for each sample
temp_mycounts_summarized <- t_temp_mycounts %>%
  group_by(sample_name) %>%
  summarize(across(everything(), sum, na.rm = TRUE))
# transpose again the df and clean up to build final count table 
temp <- as.data.frame(t(temp_mycounts_summarized))
colnames(temp) <- temp[1,]
temp <- temp[-1,]
#final count table 
mycounts <- temp
row_names_mycounts <- rownames(mycounts)
#make mycounts numeric 
mycounts <- as.data.frame(lapply(mycounts, as.numeric))
rownames(mycounts) <- row_names_mycounts
#remove X again from columns 
colnames(mycounts) <- gsub("X", "", colnames(mycounts))

colnames(mycounts)
rownames(metadata)

#remove supra glottic 70844 from metadata as it does not present in MTS 
metadata <- metadata [-68,]
#remove 70633 patient from mycounts as having COVID 
mycounts <- mycounts %>% select(rownames(metadata))

#see what is different 
setdiff(colnames(mycounts), rownames(metadata))


# Arrange mycounts
mycounts <- mycounts [, order(colnames(mycounts))]
metadata <- metadata [order(rownames(metadata)),]

#be sure every thing match 
table(colnames(mycounts)==rownames(metadata))
nrow(metadata)
ncol(mycounts)


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
metadata$five_y_mort

####### construct table #######
physeq_otu_table <- phyloseq::otu_table(as.matrix(mycounts), taxa_are_rows = TRUE)
#taxa table is not available for functional analysis 
#taxa_table <- read_tsv(file = "metatranscriptomics/annotation.taxa.bracken.mg.tsv")
#taxa_table <- as.data.frame(taxa_table)
#rownames(taxa_table) <- taxa_table$taxa
#filter humans 
#taxa_table <- taxa_table %>% dplyr::filter(taxa != human)
#physeq_taxa_table <- phyloseq::tax_table(as.matrix(taxa_table))

#merge phyloseq
physeq <- phyloseq(physeq_otu_table)
sample_data_physeq <- sample_data(metadata)
physeq <- merge_phyloseq(physeq, sample_data_physeq )

#see physeq 
physeq


##If you want to nomalize OTU table before
## To normalize data you need to set a function
normalizeSample = function(x) {x/sum(x)}

#put relattive table 
physeq.rel <- transform_sample_counts(physeq, normalizeSample)

bal_table <- subset_samples(physeq, Sample_Type=="BAL")
bal_table_rel <- subset_samples(physeq.rel, Sample_Type=="BAL")


#output of BAL count table and metadata to calculate MT-MG delta 
MT_count <- as.data.frame(otu_table(bal_table))
MT_meta <- data.frame(sample_data(bal_table))

write.table(MT_count, file = "MT_metabolic_functional_counts.txt", sep = "\t", row.names = TRUE)
write.csv(MT_meta, file = "MT_metabolic_metadata.csv")

######### RECIST comparsion 

ps <- bal_table
ps.rel <- bal_table_rel
#alpha diversity 

#calculate alpha diversity 
alpha.measures <- data.frame("Shannon" = estimate_richness(ps, measures = "Shannon"), 
                             "Observed" = estimate_richness(ps, measures = "Observed"), 
                             "Simpson" = estimate_richness(ps, measures = "Simpson"),
                             "ultimate_response_2_lev" = sample_data(ps)$ultimate_response_2_lev)
alpha.measures <- alpha.measures %>% 
  mutate(ultimate_response_2_lev= factor(ultimate_response_2_lev, levels=c("non_PD", "PD")))

#calculate stats for shannon
my.comparisons <- compare_means(Shannon~ultimate_response_2_lev, data = alpha.measures)
#save it
#write.csv(my.comparisons, file = "metatranscriptomics/Results/Shannon.stats_ultimate_response_2_lev.csv")

# plot it
pdf(file = "metatranscriptomics/Figures/functional/metabolic_pathways_Shannon_ultimate_response_2_lev.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=ultimate_response_2_lev, y=Shannon, fill=ultimate_response_2_lev))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("blue","red"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c("Non PD", "PD"))+
  stat_compare_means(comparisons = list(c("non_PD", "PD")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

#observed 
pdf(file = "metatranscriptomics/Figures/functional/metabolic_pathways_observed_ultimate_response_2_lev.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=ultimate_response_2_lev, y=Observed, fill=ultimate_response_2_lev))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("blue","red"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Observed")+
  scale_x_discrete(labels = c("Non PD", "PD"))+
  stat_compare_means(comparisons = list(c("non_PD", "PD")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()



###Beta Diversity

#Create Distance Matrix
vegdist = vegdist(t(otu_table(ps.rel)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(ps.rel), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ ultimate_response_2_lev,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="ultimate_response_2_lev",suffixes=c("",".centroid"))

newResults$ultimate_response_2_lev <- factor(newResults$ultimate_response_2_lev, levels = c("non_PD", "PD"))
centroids$ultimate_response_2_lev <- factor(centroids$ultimate_response_2_lev, levels = c("non_PD", "PD"))

pdf(file = "metatranscriptomics/Figures/functional/metabolic_pathways_Beta.Diversity.Brayultimate_response_2_lev.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= ultimate_response_2_lev)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("blue","red")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= ultimate_response_2_lev), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= ultimate_response_2_lev)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Non PD", "PD")), size=14) +
  ggtitle("Beta Diversity, Bray, ultimate_response_2_lev, p=0.068")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#sats
#write.table(adonis2(vegdist ~ newResults$ultimate_response_2_lev), 
#            file = "metatranscriptomics/Results/Beta.Diversity.Brayultimate_response_2_lev.txt", sep = "\t", row.names = T)




#############edgeR 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "ultimate_response_2_lev", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

############continue edgeR here


res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)

#add category 
res$category <- ifelse(res$logFC>0, 1, 2)

#add category detailed 
res$category_2 <- ifelse(res$category=="1", "PD", "non_PD")

res$pathway <- rownames(res)

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, ultimate_response_2_lev=="PD")
ps.rel.2 <- subset_samples(ps.rel, ultimate_response_2_lev=="non_PD")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[res$pathway]
res$abundance.1 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[res$pathway]
res$abundance.2 <- meanRA.save

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

#add contaminants 
#res <- inner_join(res, contamlist, by="taxa")
#add contam color 
#res$color <- ifelse(res$contaminant =="TRUE", "red", "black")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(logFC, PValue, FDR, abundance.1, abundance.2,
                category, category_2, pathway) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "metatranscriptomics/Results/functional/metabolic_pathways_edgeR.results_ultimate_response_2_lev_PD_vs_non_PD.csv")

#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance.1, 
                        res$abundance.2)

#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC > 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC < 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
#update color vector 
#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"red", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "blue","darkgrey"))

#add contam color 
#res.bubble.filtered$taxa_col <- res.bubble.filtered$taxa
#res.bubble.filtered$taxa_col <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$taxa_col, "</span>")

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(pathway = factor(pathway, levels = unique(pathway))) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))

#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(pathway, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(pathway, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c("blue", "red"))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_markdown(face="bold",size=16),
        axis.ticks=element_line(colour="black"),
        legend.background = element_rect(color=NA), 
        plot.title=element_text(size=23, face="bold"))+
  xlab("") +
  ylab("")+
  ggtitle("EdgeR")+
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  guides(fill="none")

pdf(file="metatranscriptomics/Figures/functional/metabolic_pathways_EdgeR_bubble_ultimate_response_2_lev.pdf", height = 11, width = 14)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="metatranscriptomics/Figures/functional/metabolic_pathwys_EdgeR_bubble_ultimate_response_2_lev_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 





####### one year mortality 

ps <- subset_samples(bal_table, one_y_mort != "NA")
ps.rel <- subset_samples(bal_table_rel, one_y_mort != "NA")
#alpha diversity 

#calculate alpha diversity 
alpha.measures <- data.frame("Shannon" = estimate_richness(ps, measures = "Shannon"), 
                             "Observed" = estimate_richness(ps, measures = "Observed"), 
                             "Simpson" = estimate_richness(ps, measures = "Simpson"),
                             "one_y_mort" = sample_data(ps)$one_y_mort)
alpha.measures <- alpha.measures %>% 
  mutate(one_y_mort= factor(one_y_mort, levels=c("greater_12", "less_12")))

#calculate stats for shannon
my.comparisons <- compare_means(Shannon~one_y_mort, data = alpha.measures)
#save it
#write.csv(my.comparisons, file = "metatranscriptomics/Results/Shannon.stats_one_y_mort.csv")

# plot it
pdf(file = "metatranscriptomics/Figures/functional/metabolic_pathways_Shannon_one_y_mort.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=one_y_mort, y=Shannon, fill=one_y_mort))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("orange3","red3"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c("geater 12", "less_12"))+
  stat_compare_means(comparisons = list(c("greater_12", "less_12")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

#observed 
pdf(file = "metatranscriptomics/Figures/functional/metabolic_pathways_observed_one_y_mort.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=one_y_mort, y=Observed, fill=one_y_mort))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("orange3","red3"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Observed")+
  scale_x_discrete(labels = c("geater 12", "less_12"))+
  stat_compare_means(comparisons = list(c("greater_12", "less_12")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()



###Beta Diversity

#Create Distance Matrix
vegdist = vegdist(t(otu_table(ps.rel)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(ps.rel), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ one_y_mort,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="one_y_mort",suffixes=c("",".centroid"))

newResults$one_y_mort <- factor(newResults$one_y_mort, levels = c("greater_12", "less_12"))
centroids$one_y_mort <- factor(centroids$one_y_mort, levels = c("greater_12", "less_12"))

pdf(file = "metatranscriptomics/Figures/functional/metabolic_pathways_Beta.Diversity.Brayone_y_mort.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= one_y_mort)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("orange3","red3")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= one_y_mort), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= one_y_mort)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("geater 12", "less_12")), size=14) +
  ggtitle("Beta Diversity, Bray, one_y_mort, p=0.001")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#sats
#write.table(adonis2(vegdist ~ newResults$one_y_mort), 
#            file = "metatranscriptomics/Results/Beta.Diversity.Brayone_y_mort.txt", sep = "\t", row.names = T)




#############edgeR 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "one_y_mort", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

############continue edgeR here


res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)

#add category 
res$category <- ifelse(res$logFC>0, 1, 2)

#add category detailed 
res$category_2 <- ifelse(res$category=="1", "less_12", "greater_12")

res$pathway <- rownames(res)

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, one_y_mort=="less_12")
ps.rel.2 <- subset_samples(ps.rel, one_y_mort=="greater_12")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[res$pathway]
res$abundance.1 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[res$pathway]
res$abundance.2 <- meanRA.save

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

#add contaminants 
#res <- inner_join(res, contamlist, by="taxa")
#add contam color 
#res$color <- ifelse(res$contaminant =="TRUE", "red", "black")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(logFC, PValue, FDR, abundance.1, abundance.2,
                category, category_2, pathway) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "metatranscriptomics/Results/functional/metabolic_pathways_edgeR.results_one_y_mort_less_12_vs_greater_12.csv")

#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance.1, 
                        res$abundance.2)

#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC > 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC < 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
#update color vector 
#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"red3", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "orange3","darkgrey"))

#add contam color 
#res.bubble.filtered$taxa_col <- res.bubble.filtered$taxa
#res.bubble.filtered$taxa_col <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$taxa_col, "</span>")

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(pathway = factor(pathway, levels = unique(pathway))) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))

#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(pathway, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(pathway, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c("orange3", "red3"))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_markdown(face="bold",size=16),
        axis.ticks=element_line(colour="black"),
        legend.background = element_rect(color=NA), 
        plot.title=element_text(size=23, face="bold"))+
  xlab("") +
  ylab("")+
  ggtitle("EdgeR")+
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  guides(fill="none")

pdf(file="metatranscriptomics/Figures/functional/metabolic_pathways_EdgeR_bubble_one_y_mort.pdf", height = 11, width = 14)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="metatranscriptomics/Figures/functional/metabolic_pathwys_EdgeR_bubble_one_y_mort_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 








######### two year mortatlity 

ps <- subset_samples(bal_table, two_y_mort != "NA")
ps.rel <- subset_samples(bal_table_rel, two_y_mort != "NA")
#alpha diversity 

#calculate alpha diversity 
alpha.measures <- data.frame("Shannon" = estimate_richness(ps, measures = "Shannon"), 
                             "Observed" = estimate_richness(ps, measures = "Observed"), 
                             "Simpson" = estimate_richness(ps, measures = "Simpson"),
                             "two_y_mort" = sample_data(ps)$two_y_mort)
alpha.measures <- alpha.measures %>% 
  mutate(two_y_mort= factor(two_y_mort, levels=c("greater_24", "less_24")))

#calculate stats for shannon
my.comparisons <- compare_means(Shannon~two_y_mort, data = alpha.measures)
#save it
#write.csv(my.comparisons, file = "metatranscriptomics/Results/Shannon.stats_two_y_mort.csv")

# plot it
pdf(file = "metatranscriptomics/Figures/functional/metabolic_pathways_Shannon_two_y_mort.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=two_y_mort, y=Shannon, fill=two_y_mort))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("orange3","red3"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c("geater 24", "less_24"))+
  stat_compare_means(comparisons = list(c("greater_24", "less_24")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

#observed 
pdf(file = "metatranscriptomics/Figures/functional/metabolic_pathways_observed_two_y_mort.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=two_y_mort, y=Observed, fill=two_y_mort))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("orange3","red3"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Observed")+
  scale_x_discrete(labels = c("geater 24", "less_24"))+
  stat_compare_means(comparisons = list(c("greater_24", "less_24")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()



###Beta Diversity

#Create Distance Matrix
vegdist = vegdist(t(otu_table(ps.rel)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(ps.rel), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ two_y_mort,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="two_y_mort",suffixes=c("",".centroid"))

newResults$two_y_mort <- factor(newResults$two_y_mort, levels = c("greater_24", "less_24"))
centroids$two_y_mort <- factor(centroids$two_y_mort, levels = c("greater_24", "less_24"))

pdf(file = "metatranscriptomics/Figures/functional/metabolic_pathways_Beta.Diversity.Braytwo_y_mort.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= two_y_mort)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("orange3","red3")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= two_y_mort), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= two_y_mort)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("geater 24", "less_24")), size=14) +
  ggtitle("Beta Diversity, Bray, two_y_mort, p=0.087")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#sats
#write.table(adonis2(vegdist ~ newResults$two_y_mort), 
#            file = "metatranscriptomics/Results/Beta.Diversity.Braytwo_y_mort.txt", sep = "\t", row.names = T)




#############edgeR 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "two_y_mort", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

############continue edgeR here


res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)

#add category 
res$category <- ifelse(res$logFC>0, 1, 2)

#add category detailed 
res$category_2 <- ifelse(res$category=="1", "less_24", "greater_24")

res$pathway <- rownames(res)

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, two_y_mort=="less_24")
ps.rel.2 <- subset_samples(ps.rel, two_y_mort=="greater_24")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[res$pathway]
res$abundance.1 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[res$pathway]
res$abundance.2 <- meanRA.save

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

#add contaminants 
#res <- inner_join(res, contamlist, by="taxa")
#add contam color 
#res$color <- ifelse(res$contaminant =="TRUE", "red", "black")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(logFC, PValue, FDR, abundance.1, abundance.2,
                category, category_2, pathway) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "metatranscriptomics/Results/functional/metabolic_pathways_edgeR.results_two_y_mort_less_24_vs_greater_24.csv")

#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance.1, 
                        res$abundance.2)

#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC > 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC < 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
#update color vector 
#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"red3", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "orange3","darkgrey"))

#add contam color 
#res.bubble.filtered$taxa_col <- res.bubble.filtered$taxa
#res.bubble.filtered$taxa_col <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$taxa_col, "</span>")

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(pathway = factor(pathway, levels = unique(pathway))) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))

#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(pathway, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(pathway, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c("orange3", "red3"))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_markdown(face="bold",size=16),
        axis.ticks=element_line(colour="black"),
        legend.background = element_rect(color=NA), 
        plot.title=element_text(size=23, face="bold"))+
  xlab("") +
  ylab("")+
  ggtitle("EdgeR")+
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  guides(fill="none")

pdf(file="metatranscriptomics/Figures/functional/metabolic_pathways_EdgeR_bubble_two_y_mort.pdf", height = 11, width = 14)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="metatranscriptomics/Figures/functional/metabolic_pathwys_EdgeR_bubble_two_y_mort_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 




########### KO pathways analysis ###########

#read metadata 
metadata <- read.delim2(file = "immuntherapy_metadata_final_cohort.txt")

#make genes names rownames 
names <- make.unique(metadata$X)
rownames(metadata) <- names

#read counts table 
mycounts <- read.delim2(file = "metatranscriptomics/Data/upd_tables/counts.ko.pathway.fmap.MT.tsv")
#make genes names rownames 
names <- make.unique(mycounts$X)
rownames(mycounts) <- names

#remove X column 
mycounts <- mycounts %>% dplyr::select(-X)
metadata <- metadata %>% dplyr::select(-X)

rownames(metadata) <- gsub("_", ".", rownames(metadata))

#mycounts is arranged by lanes, so we are going to arrange it, transpose and sum all reads per one sample 
temp_mycounts <- mycounts
colnames(temp_mycounts) <- gsub("X", "", colnames(temp_mycounts))
colnames(temp_mycounts) <- gsub("_S\\d+_L\\d+", "", colnames(temp_mycounts))
colnames(temp_mycounts) <- gsub("_", ".", colnames(temp_mycounts))
#transpose mycounts
t_temp_mycounts <- as.data.frame(t(temp_mycounts))
# create a column called sample_name and clean it up so you can group by each sample later 
t_temp_mycounts$sample_name <- rownames(t_temp_mycounts)
t_temp_mycounts$sample_name <- gsub("X", "", t_temp_mycounts$sample_name)
t_temp_mycounts$sample_name <- gsub("\\.\\d+$", "", t_temp_mycounts$sample_name)
#sum all reads for each sample
temp_mycounts_summarized <- t_temp_mycounts %>%
  group_by(sample_name) %>%
  summarize(across(everything(), sum, na.rm = TRUE))
# transpose again the df and clean up to build final count table 
temp <- as.data.frame(t(temp_mycounts_summarized))
colnames(temp) <- temp[1,]
temp <- temp[-1,]
#final count table 
mycounts <- temp
row_names_mycounts <- rownames(mycounts)
#make mycounts numeric 
mycounts <- as.data.frame(lapply(mycounts, as.numeric))
rownames(mycounts) <- row_names_mycounts
#remove X again from columns 
colnames(mycounts) <- gsub("X", "", colnames(mycounts))

colnames(mycounts)
rownames(metadata)

#remove supra glottic 70844 from metadata as it does not present in MTS 
metadata <- metadata [-68,]
#remove 70633 patient from mycounts as having COVID 
mycounts <- mycounts %>% select(rownames(metadata))

#see what is different 
setdiff(colnames(mycounts), rownames(metadata))


# Arrange mycounts
mycounts <- mycounts [, order(colnames(mycounts))]
metadata <- metadata [order(rownames(metadata)),]

#be sure every thing match 
table(colnames(mycounts)==rownames(metadata))
nrow(metadata)
ncol(mycounts)



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
metadata$five_y_mort

####### construct table #######
physeq_otu_table <- phyloseq::otu_table(as.matrix(mycounts), taxa_are_rows = TRUE)
#taxa table is not available for functional analysis 
#taxa_table <- read_tsv(file = "metatranscriptomics/annotation.taxa.bracken.mg.tsv")
#taxa_table <- as.data.frame(taxa_table)
#rownames(taxa_table) <- taxa_table$taxa
#filter humans 
#taxa_table <- taxa_table %>% dplyr::filter(taxa != human)
#physeq_taxa_table <- phyloseq::tax_table(as.matrix(taxa_table))

#merge phyloseq
physeq <- phyloseq(physeq_otu_table)
sample_data_physeq <- sample_data(metadata)
physeq <- merge_phyloseq(physeq, sample_data_physeq )

#see physeq 
physeq


##If you want to nomalize OTU table before
## To normalize data you need to set a function
normalizeSample = function(x) {x/sum(x)}

#put relattive table 
physeq.rel <- transform_sample_counts(physeq, normalizeSample)

bal_table <- subset_samples(physeq, Sample_Type=="BAL")
bal_table_rel <- subset_samples(physeq.rel, Sample_Type=="BAL")


######### RECIST comparsion 

ps <- bal_table
ps.rel <- bal_table_rel
#alpha diversity 

#calculate alpha diversity 
alpha.measures <- data.frame("Shannon" = estimate_richness(ps, measures = "Shannon"), 
                             "Observed" = estimate_richness(ps, measures = "Observed"), 
                             "Simpson" = estimate_richness(ps, measures = "Simpson"),
                             "ultimate_response_2_lev" = sample_data(ps)$ultimate_response_2_lev)
alpha.measures <- alpha.measures %>% 
  mutate(ultimate_response_2_lev= factor(ultimate_response_2_lev, levels=c("non_PD", "PD")))

#calculate stats for shannon
my.comparisons <- compare_means(Shannon~ultimate_response_2_lev, data = alpha.measures)
#save it
#write.csv(my.comparisons, file = "metatranscriptomics/Results/Shannon.stats_ultimate_response_2_lev.csv")

# plot it
pdf(file = "metatranscriptomics/Figures/functional/KO_pathways_Shannon_ultimate_response_2_lev.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=ultimate_response_2_lev, y=Shannon, fill=ultimate_response_2_lev))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("blue","red"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c("Non PD", "PD"))+
  stat_compare_means(comparisons = list(c("non_PD", "PD")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

#observed 
pdf(file = "metatranscriptomics/Figures/functional/KO_pathways_observed_ultimate_response_2_lev.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=ultimate_response_2_lev, y=Observed, fill=ultimate_response_2_lev))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("blue","red"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Observed")+
  scale_x_discrete(labels = c("Non PD", "PD"))+
  stat_compare_means(comparisons = list(c("non_PD", "PD")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()



###Beta Diversity

#Create Distance Matrix
vegdist = vegdist(t(otu_table(ps.rel)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(ps.rel), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ ultimate_response_2_lev,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="ultimate_response_2_lev",suffixes=c("",".centroid"))

newResults$ultimate_response_2_lev <- factor(newResults$ultimate_response_2_lev, levels = c("non_PD", "PD"))
centroids$ultimate_response_2_lev <- factor(centroids$ultimate_response_2_lev, levels = c("non_PD", "PD"))

pdf(file = "metatranscriptomics/Figures/functional/KO_pathways_Beta.Diversity.Brayultimate_response_2_lev.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= ultimate_response_2_lev)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("blue","red")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= ultimate_response_2_lev), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= ultimate_response_2_lev)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Non PD", "PD")), size=14) +
  ggtitle("Beta Diversity, Bray, ultimate_response_2_lev, p=0.019")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#sats
#write.table(adonis2(vegdist ~ newResults$ultimate_response_2_lev), 
#            file = "metatranscriptomics/Results/Beta.Diversity.Brayultimate_response_2_lev.txt", sep = "\t", row.names = T)




#############edgeR 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "ultimate_response_2_lev", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

############continue edgeR here


res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)

#add category 
res$category <- ifelse(res$logFC>0, 1, 2)

#add category detailed 
res$category_2 <- ifelse(res$category=="1", "PD", "non_PD")

res$pathway <- rownames(res)

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, ultimate_response_2_lev=="PD")
ps.rel.2 <- subset_samples(ps.rel, ultimate_response_2_lev=="non_PD")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[res$pathway]
res$abundance.1 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[res$pathway]
res$abundance.2 <- meanRA.save

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

#add contaminants 
#res <- inner_join(res, contamlist, by="taxa")
#add contam color 
#res$color <- ifelse(res$contaminant =="TRUE", "red", "black")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(logFC, PValue, FDR, abundance.1, abundance.2,
                category, category_2, pathway) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "metatranscriptomics/Results/functional/KO_pathways_edgeR.results_ultimate_response_2_lev_PD_vs_non_PD.csv")

#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance.1, 
                        res$abundance.2)

#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC > 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC < 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
#update color vector 
#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"red", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "blue","darkgrey"))

#add contam color 
#res.bubble.filtered$taxa_col <- res.bubble.filtered$taxa
#res.bubble.filtered$taxa_col <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$taxa_col, "</span>")

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(pathway = factor(pathway, levels = unique(pathway))) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))

#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(pathway, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(pathway, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c("blue", "red"))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_markdown(face="bold",size=16),
        axis.ticks=element_line(colour="black"),
        legend.background = element_rect(color=NA), 
        plot.title=element_text(size=23, face="bold"))+
  xlab("") +
  ylab("")+
  ggtitle("EdgeR")+
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  guides(fill="none")

pdf(file="metatranscriptomics/Figures/functional/KO_pathways_EdgeR_bubble_ultimate_response_2_lev.pdf", height = 11, width = 14)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="metatranscriptomics/Figures/functional/KO_pathwys_EdgeR_bubble_ultimate_response_2_lev_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 





####### one year mortality 

ps <- subset_samples(bal_table, one_y_mort != "NA")
ps.rel <- subset_samples(bal_table_rel, one_y_mort != "NA")
#alpha diversity 

#calculate alpha diversity 
alpha.measures <- data.frame("Shannon" = estimate_richness(ps, measures = "Shannon"), 
                             "Observed" = estimate_richness(ps, measures = "Observed"), 
                             "Simpson" = estimate_richness(ps, measures = "Simpson"),
                             "one_y_mort" = sample_data(ps)$one_y_mort)
alpha.measures <- alpha.measures %>% 
  mutate(one_y_mort= factor(one_y_mort, levels=c("greater_12", "less_12")))

#calculate stats for shannon
my.comparisons <- compare_means(Shannon~one_y_mort, data = alpha.measures)
#save it
#write.csv(my.comparisons, file = "metatranscriptomics/Results/Shannon.stats_one_y_mort.csv")

# plot it
pdf(file = "metatranscriptomics/Figures/functional/KO_pathways_Shannon_one_y_mort.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=one_y_mort, y=Shannon, fill=one_y_mort))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("orange3","red3"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c("geater 12", "less_12"))+
  stat_compare_means(comparisons = list(c("greater_12", "less_12")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

#observed 
pdf(file = "metatranscriptomics/Figures/functional/KO_pathways_observed_one_y_mort.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=one_y_mort, y=Observed, fill=one_y_mort))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("orange3","red3"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Observed")+
  scale_x_discrete(labels = c("geater 12", "less_12"))+
  stat_compare_means(comparisons = list(c("greater_12", "less_12")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()



###Beta Diversity

#Create Distance Matrix
vegdist = vegdist(t(otu_table(ps.rel)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(ps.rel), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ one_y_mort,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="one_y_mort",suffixes=c("",".centroid"))

newResults$one_y_mort <- factor(newResults$one_y_mort, levels = c("greater_12", "less_12"))
centroids$one_y_mort <- factor(centroids$one_y_mort, levels = c("greater_12", "less_12"))

pdf(file = "metatranscriptomics/Figures/functional/KO_pathways_Beta.Diversity.Brayone_y_mort.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= one_y_mort)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("orange3","red3")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= one_y_mort), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= one_y_mort)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("geater 12", "less_12")), size=14) +
  ggtitle("Beta Diversity, Bray, one_y_mort, p=0.001")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#sats
#write.table(adonis2(vegdist ~ newResults$one_y_mort), 
#            file = "metatranscriptomics/Results/Beta.Diversity.Brayone_y_mort.txt", sep = "\t", row.names = T)




#############edgeR 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "one_y_mort", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

############continue edgeR here


res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)

#add category 
res$category <- ifelse(res$logFC>0, 1, 2)

#add category detailed 
res$category_2 <- ifelse(res$category=="1", "less_12", "greater_12")

res$pathway <- rownames(res)

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, one_y_mort=="less_12")
ps.rel.2 <- subset_samples(ps.rel, one_y_mort=="greater_12")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[res$pathway]
res$abundance.1 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[res$pathway]
res$abundance.2 <- meanRA.save

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

#add contaminants 
#res <- inner_join(res, contamlist, by="taxa")
#add contam color 
#res$color <- ifelse(res$contaminant =="TRUE", "red", "black")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(logFC, PValue, FDR, abundance.1, abundance.2,
                category, category_2, pathway) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "metatranscriptomics/Results/functional/KO_pathways_edgeR.results_one_y_mort_less_12_vs_greater_12.csv")

#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance.1, 
                        res$abundance.2)

#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC > 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC < 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
#update color vector 
#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"red3", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "orange3","darkgrey"))

#add contam color 
#res.bubble.filtered$taxa_col <- res.bubble.filtered$taxa
#res.bubble.filtered$taxa_col <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$taxa_col, "</span>")

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(pathway = factor(pathway, levels = unique(pathway))) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))

#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(pathway, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(pathway, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c("orange3", "red3"))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_markdown(face="bold",size=16),
        axis.ticks=element_line(colour="black"),
        legend.background = element_rect(color=NA), 
        plot.title=element_text(size=23, face="bold"))+
  xlab("") +
  ylab("")+
  ggtitle("EdgeR")+
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  guides(fill="none")

pdf(file="metatranscriptomics/Figures/functional/KO_pathways_EdgeR_bubble_one_y_mort.pdf", height = 11, width = 14)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="metatranscriptomics/Figures/functional/KO_pathwys_EdgeR_bubble_one_y_mort_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 








######### two year mortatlity 

ps <- subset_samples(bal_table, two_y_mort != "NA")
ps.rel <- subset_samples(bal_table_rel, two_y_mort != "NA")
#alpha diversity 

#calculate alpha diversity 
alpha.measures <- data.frame("Shannon" = estimate_richness(ps, measures = "Shannon"), 
                             "Observed" = estimate_richness(ps, measures = "Observed"), 
                             "Simpson" = estimate_richness(ps, measures = "Simpson"),
                             "two_y_mort" = sample_data(ps)$two_y_mort)
alpha.measures <- alpha.measures %>% 
  mutate(two_y_mort= factor(two_y_mort, levels=c("greater_24", "less_24")))

#calculate stats for shannon
my.comparisons <- compare_means(Shannon~two_y_mort, data = alpha.measures)
#save it
#write.csv(my.comparisons, file = "metatranscriptomics/Results/Shannon.stats_two_y_mort.csv")

# plot it
pdf(file = "metatranscriptomics/Figures/functional/KO_pathways_Shannon_two_y_mort.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=two_y_mort, y=Shannon, fill=two_y_mort))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("orange3","red3"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c("geater 24", "less_24"))+
  stat_compare_means(comparisons = list(c("greater_24", "less_24")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

#observed 
pdf(file = "metatranscriptomics/Figures/functional/KO_pathways_observed_two_y_mort.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=two_y_mort, y=Observed, fill=two_y_mort))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("orange3","red3"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Observed")+
  scale_x_discrete(labels = c("geater 24", "less_24"))+
  stat_compare_means(comparisons = list(c("greater_24", "less_24")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()



###Beta Diversity

#Create Distance Matrix
vegdist = vegdist(t(otu_table(ps.rel)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(ps.rel), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ two_y_mort,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="two_y_mort",suffixes=c("",".centroid"))

newResults$two_y_mort <- factor(newResults$two_y_mort, levels = c("greater_24", "less_24"))
centroids$two_y_mort <- factor(centroids$two_y_mort, levels = c("greater_24", "less_24"))

pdf(file = "metatranscriptomics/Figures/functional/KO_pathways_Beta.Diversity.Braytwo_y_mort.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= two_y_mort)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("orange3","red3")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= two_y_mort), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= two_y_mort)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("geater 24", "less_24")), size=14) +
  ggtitle("Beta Diversity, Bray, two_y_mort, p=0.005")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#sats
#write.table(adonis2(vegdist ~ newResults$two_y_mort), 
#            file = "metatranscriptomics/Results/Beta.Diversity.Braytwo_y_mort.txt", sep = "\t", row.names = T)




#############edgeR 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "two_y_mort", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

############continue edgeR here


res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)

#add category 
res$category <- ifelse(res$logFC>0, 1, 2)

#add category detailed 
res$category_2 <- ifelse(res$category=="1", "less_24", "greater_24")

res$pathway <- rownames(res)

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, two_y_mort=="less_24")
ps.rel.2 <- subset_samples(ps.rel, two_y_mort=="greater_24")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[res$pathway]
res$abundance.1 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[res$pathway]
res$abundance.2 <- meanRA.save

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

#add contaminants 
#res <- inner_join(res, contamlist, by="taxa")
#add contam color 
#res$color <- ifelse(res$contaminant =="TRUE", "red", "black")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(logFC, PValue, FDR, abundance.1, abundance.2,
                category, category_2, pathway) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "metatranscriptomics/Results/functional/KO_pathways_edgeR.results_two_y_mort_less_24_vs_greater_24.csv")

#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance.1, 
                        res$abundance.2)

#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC > 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC < 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
#update color vector 
#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"red3", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "orange3","darkgrey"))

#add contam color 
#res.bubble.filtered$taxa_col <- res.bubble.filtered$taxa
#res.bubble.filtered$taxa_col <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$taxa_col, "</span>")

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(pathway = factor(pathway, levels = unique(pathway))) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))

#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(pathway, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(pathway, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c("orange3", "red3"))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_markdown(face="bold",size=16),
        axis.ticks=element_line(colour="black"),
        legend.background = element_rect(color=NA), 
        plot.title=element_text(size=23, face="bold"))+
  xlab("") +
  ylab("")+
  ggtitle("EdgeR")+
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  guides(fill="none")

pdf(file="metatranscriptomics/Figures/functional/KO_pathways_EdgeR_bubble_two_y_mort.pdf", height = 11, width = 14)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="metatranscriptomics/Figures/functional/KO_pathwys_EdgeR_bubble_two_y_mort_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 


save.image("MT_functional_final.RData")



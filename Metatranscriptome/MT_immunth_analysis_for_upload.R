######## Metatrnascriptome analysis of lower airway microbiome 
######## paper: immunotherapy lung cancer paper 
######## Fares Darawshy MD
######## Sep 2025 


### set your wording directory  
setwd("")

#load necessary packages 
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


########## prepare and read files########

#read metadata 
metadata <- read.delim2(file = "immuntherapy_metadata_final_cohort.txt")

#make genes names rownames 
names <- make.unique(metadata$X)
rownames(metadata) <- names

#read counts table 
mycounts <- read.delim2(file = "metatranscriptomics/Data/upd_tables/counts.taxa.bracken.MT.tsv")

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


######## build PS object while removing homo sapiens #######
human <- "Homo sapiens"
mycounts <- mycounts %>% dplyr::filter(rownames(.) != human)
rownames(mycounts)
physeq_otu_table <- phyloseq::otu_table(as.matrix(mycounts), taxa_are_rows = TRUE)
#taxa table
taxa_table <- read_tsv(file = "metatranscriptomics/Data/upd_tables/annotation.taxa.bracken.MT.tsv")
taxa_table <- as.data.frame(taxa_table)
rownames(taxa_table) <- taxa_table$taxa
#filter humans 
taxa_table <- taxa_table %>% dplyr::filter(taxa != human)
physeq_taxa_table <- phyloseq::tax_table(as.matrix(taxa_table))

#merge phyloseq
physeq <- phyloseq(physeq_otu_table, taxa_table)
sample_data_physeq <- sample_data(metadata)
physeq <- merge_phyloseq(physeq, sample_data_physeq,physeq_taxa_table )

#see physeq 
physeq


##If you want to nomalize OTU table before
## To normalize data you need to set a function
normalizeSample = function(x) {x/sum(x)}

#put relattive table 
physeq.rel <- transform_sample_counts(physeq, normalizeSample)




########### every comparison include the following: bacterial load, alpha diversity, beta diversity, differential anlaysis 


##################comparison of BAL, BKG and Sup ##########

ps <- physeq
ps.rel <- physeq.rel

#############Sequence Detph

#remove all taxa with 0 abundance
otu.table.sub = subset_taxa(ps, rowSums(otu_table(ps)) != 0) 

sample_sums(otu.table.sub)
sort(sample_sums(otu.table.sub))

sample_data(otu.table.sub)$Sample_Type

sequence_depth.sub <- sample_sums(otu.table.sub)
Sample_Type.sub <- sample_data(otu.table.sub)$Sample_Type

sequence.per.sample.sub <- data.frame(Sample_Type.sub, sequence_depth.sub)
sequence.per.sample.sub$Sample_Type.sub <- factor(sequence.per.sample.sub$Sample_Type.sub, levels = c("BKG", "BAL", "Supraglottic"))


sequence.per.sample.sub %>% 
  summarise(median=median((sequence_depth.sub)), 
            q1=quantile(sequence_depth.sub, 0.25), 
            q3=quantile(sequence_depth.sub, 0.75))
#    median       q1       q3
#1 21024868 13641717 45677011


sequence.per.sample.sub %>% 
  filter(Sample_Type.sub=="BAL") %>% 
  summarise(median=median((sequence_depth.sub)), 
            q1=quantile(sequence_depth.sub, 0.25), 
            q3=quantile(sequence_depth.sub, 0.75))
#    median       q1       q3
#1 18865858 11802322 28763560



#stats 
my.comparisons <- compare_means(sequence_depth.sub~Sample_Type.sub, data = sequence.per.sample.sub )
#write.csv(my.comparisons, file = "Results/metagenome/Sequence.depth.stats.All.samples.csv")

# Colored boxplot according to sample type
pdf(file = "metatranscriptomics/Figures/Sequence.depth.All.samples.pdf"
    , height = 10, width = 7)
ggplot(data=sequence.per.sample.sub, mapping = aes(x=Sample_Type.sub, y=sequence_depth.sub, fill=Sample_Type.sub))+
  geom_boxplot()+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("darkgrey", "blue", "goldenrod"))+
  geom_jitter(color="black", size=0.4, alpha=0.7)+
  labs(x= '', y='Total Reads')+
  scale_x_discrete(labels = c("BKG", "BAL", "UA"))+
  scale_y_continuous(trans = "log10")+
  guides(fill="none")+
  stat_compare_means(comparisons = list(c("BKG", "BAL"), c("Supraglottic", "BAL"),c("BKG", "Supraglottic")), 
                     method = "wilcox.test", label = "p.value", size = 10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))
dev.off()


###########ddpcr
ddpcr_dat <- data.frame(ps@sam_data)

ddpcr_dat <- ddpcr_dat %>% dplyr::mutate(ddPCR=as.numeric(ddPCR))
#calculate stats for ddpcr
my.comparisons <- compare_means(ddPCR~Sample_Type, data = ddpcr_dat)
#save it
#write.csv(my.comparisons, file = "Results/metagenome/ddPCR.stats.Sample_Type.csv")

ddpcr_dat$Sample_Type <- factor(ddpcr_dat$Sample_Type, levels = c("BKG", "BAL", "Supraglottic"))

# plot it
pdf(file = "metatranscriptomics/Figures/ddPCR_sample_type.pdf", width=6, height=9)
ggplot(ddpcr_dat, aes(x=Sample_Type, y=ddPCR, fill=Sample_Type))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("darkgrey", "blue", "goldenrod"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("16S rRNA copies/uL")+
  scale_x_discrete(labels = c("BKG", "BAL", "UA"))+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(comparisons = list(c("BKG", "BAL"), c("BAL", "Supraglottic"), c("BKG", "Supraglottic")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

#########Alpha diversity

#calculate alpha diversity 
alpha.measures <- data.frame("Shannon" = estimate_richness(ps, measures = "Shannon"), 
                             "Observed" = estimate_richness(ps, measures = "Observed"), 
                             "Simpson" = estimate_richness(ps, measures = "Simpson"),
                             "Sample_Type" = sample_data(ps)$Sample_Type)
alpha.measures <- alpha.measures %>% 
  mutate(Sample_Type= factor(Sample_Type))

alpha.measures$Sample_Type <- factor(alpha.measures$Sample_Type, levels = c("BKG", "BAL", "Supraglottic"))

#calculate stats for shannon
my.comparisons <- compare_means(Shannon~Sample_Type, data = alpha.measures)
#save it
#write.csv(my.comparisons, file = "metatranscriptomics/Results/Shannon.stats.Sample_Type.csv")

# plot it
pdf(file = "metatranscriptomics/Figures/Shannon_sample_type.pdf", width=7, height=10)
ggplot(alpha.measures, aes(x=Sample_Type, y=Shannon, fill=Sample_Type))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("darkgrey", "blue", "goldenrod"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c("BKG", "BAL", "UA"))+
  stat_compare_means(comparisons = list(c("BKG", "BAL"), c("BAL", "Supraglottic"), c("BKG", "Supraglottic")), 
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
pdf(file = "metatranscriptomics/Figures/observed_sample_type.pdf", width=7, height=10)
ggplot(alpha.measures, aes(x=Sample_Type, y=Observed, fill=Sample_Type))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("darkgrey", "blue", "goldenrod"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Richness")+
  scale_x_discrete(labels = c("BKG", "BAL", "UA"))+
  stat_compare_means(comparisons = list(c("BKG", "BAL"), c("BAL", "Supraglottic"), c("BKG", "Supraglottic")), 
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
centroids <- aggregate(cbind(PC1,PC2)~ Sample_Type,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Sample_Type",suffixes=c("",".centroid"))

newResults$Sample_Type <- factor(newResults$Sample_Type, levels = c("BKG", "BAL", "Supraglottic"))
centroids$Sample_Type <- factor(centroids$Sample_Type, levels = c("BKG", "BAL", "Supraglottic"))

pdf(file = "metatranscriptomics/Figures/Beta.Diversity.Bray.Sample_Type.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= Sample_Type)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("darkgrey", "blue", "goldenrod")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Sample_Type), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Sample_Type)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("BAL", "BKG", "UA")), size=14) +
  ggtitle("Beta Diversity, Bray, All samples, p=0.001")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#sats
write.table(adonis2(vegdist ~ newResults$Sample_Type), 
            file = "metatranscriptomics/Results/Beta.Diversity.Bray.Sample_Type.txt", sep = "\t", row.names = T)


#1 vs 2 (p=0.001)
ps.re.1 <- subset_samples(ps.rel, Sample_Type %in% c("BKG", "BAL"))
vegdist   = vegdist(t(otu_table(ps.re.1)), method = "bray")
write.csv(adonis2(vegdist~sample_data(ps.re.1)$Sample_Type), file = "metatranscriptomics/Results/Bray.Stats.Sample_Type_BKG_vs_BAL.csv")

#2 vs 3 (p=0.001)
ps.re.1 <- subset_samples(ps.rel, Sample_Type %in% c("BAL", "Supraglottic"))
vegdist   = vegdist(t(otu_table(ps.re.1)), method = "bray")
write.csv(adonis2(vegdist~sample_data(ps.re.1)$Sample_Type), file = "metatranscriptomics/Results/Bray.Stats.Sample_Type_BAL_vs_Supraglottic.csv")

#1 vs 3 (p=0.001)
ps.re.1 <- subset_samples(ps.rel, Sample_Type %in% c("BKG", "Supraglottic"))
vegdist   = vegdist(t(otu_table(ps.re.1)), method = "bray")
write.csv(adonis2(vegdist~sample_data(ps.re.1)$Sample_Type), file = "metatranscriptomics/Results/Bray.Stats.Sample_Type_BKG_vs_Supraglottic.csv")


######################  decontaminant analysis  ###########
physeq@tax_table@.Data

####################define KW deconta function#### open KW_functions_v4.R and run functions
# notice that this function is different than regualr function and adapted for metagenome

#threshold = 0.1
decontaminant_KW(input_phyloseq = physeq, SampleID.unique ="ID", 
                 sample_types = c("BKG", "BAL", "Supraglottic"), 
                 sample_type_var_name = "Sample_Type", 
                 sample_type_color = c("darkgrey", "blue", "goldenrod"), 
                 negative_sample_type = "BKG", 
                 compare_type = c("BAL_and_Supraglottic", "BAL_or_Supraglottic", "BAL", "BAL_combine_Supraglottic"),
                 stat_option = "mean", 
                 display_contam_method="none", #if none is selected then there will be no red label
                 test_threshold = 0.1,
                 graph_option = "boxplot",
                 log_scale = "yes", taxa_genus_output = "no")




#threshold = 0.25
decontaminant_KW(input_phyloseq = physeq, SampleID.unique ="ID", 
                 sample_types = c("BKG", "BAL", "Supraglottic"), 
                 sample_type_var_name = "Sample_Type", 
                 sample_type_color = c("darkgrey", "blue", "goldenrod"), 
                 negative_sample_type = "BKG", 
                 compare_type = c("BAL_and_Supraglottic", "BAL_or_Supraglottic", "BAL", "BAL_combine_Supraglottic"),
                 stat_option = "mean", 
                 display_contam_method="none", #if none is selected then there will be no red label
                 test_threshold = 0.25,
                 graph_option = "boxplot", 
                 taxa_genus_output="no", 
                 log_scale = "yes")

#threshold = 0.5
decontaminant_KW(input_phyloseq = physeq, SampleID.unique ="ID", 
                 sample_types = c("BKG", "BAL", "Supraglottic"), 
                 sample_type_var_name = "Sample_Type", 
                 sample_type_color = c("darkgrey", "blue", "goldenrod"), 
                 negative_sample_type = "BKG", 
                 compare_type = c("BAL_and_Supraglottic", "BAL_or_Supraglottic", "BAL", "BAL_combine_Supraglottic"),
                 stat_option = "mean", 
                 display_contam_method="none", #if none is selected then there will be no red label
                 test_threshold = 0.5,
                 graph_option = "boxplot", 
                 taxa_genus_output="no", 
                 log_scale = "yes")

#threshold = 0.75
decontaminant_KW(input_phyloseq = physeq, SampleID.unique ="ID", 
                 sample_types = c("BKG", "BAL", "Supraglottic"), 
                 sample_type_var_name = "Sample_Type", 
                 sample_type_color = c("darkgrey", "blue", "goldenrod"), 
                 negative_sample_type = "BKG", 
                 compare_type = c("BAL_and_Supraglottic", "BAL_or_Supraglottic", "BAL", "BAL_combine_Supraglottic"),
                 stat_option = "mean", 
                 display_contam_method="none", #if none is selected then there will be no red label
                 test_threshold = 0.75,
                 graph_option = "boxplot", 
                 taxa_genus_output="no", 
                 log_scale = "yes")

#threshold = 0.9
decontaminant_KW(input_phyloseq = physeq, SampleID.unique ="ID", 
                 sample_types = c("BKG", "BAL", "Supraglottic"), 
                 sample_type_var_name = "Sample_Type", 
                 sample_type_color = c("darkgrey", "blue", "goldenrod"), 
                 negative_sample_type = "BKG", 
                 compare_type = c("BAL_and_Supraglottic", "BAL_or_Supraglottic", "BAL", "BAL_combine_Supraglottic"),
                 stat_option = "mean", 
                 display_contam_method="none", #if none is selected then there will be no red label
                 test_threshold = 0.9,
                 graph_option = "boxplot", 
                 taxa_genus_output="no", 
                 log_scale = "yes")

##############################################################################
########## Final Decontam method and analysis: prevalence 0.25 cutoff , and argument #######
#use subplot function 
decontaminant_subplot_KW(input_phyloseq = physeq, SampleID.unique = "ID", 
                         sample_types = c("BKG", "BAL", "Supraglottic"), 
                         sample_type_var_name = "Sample_Type", 
                         sample_type_color = c("darkgrey", "blue", "goldenrod"), 
                         negative_sample_type = "BKG", 
                         compare_type = "BAL_and_Supraglottic",
                         stat_option = "mean", 
                         method_type = "preval",
                         test_threshold = 0.25,
                         graph_option = "boxplot", 
                         taxa_genus_output="no", 
                         log_scale = "yes")

#get final contamlist 
contamlist <- Contam_list_NC_BKG_compare_BAL_and_Supraglottic_OP

contamlist$taxa <- row.names(contamlist)

######################## differntial analysis Sample type (BAL vs Sup) #######

physeq_subset_no_BKG <- subset_samples(physeq, Sample_Type %in% c("BAL", "Supraglottic"))

#temp <- data.frame(sample_data(physeq))
#temp <- temp %>% filter(Sample_Type=="BAL", Visit=="V1", prim_sec_sample=="prim")
#write.csv(temp, "lower_airway_baseline_samples_metadata.csv")

ps <- physeq_subset_no_BKG
ps@sam_data$Sample_Type <- factor(ps@sam_data$Sample_Type, levels = c("Supraglottic", "BAL"))

ps.rel <- transform_sample_counts(ps, normalizeSample)
v = "Sample_Type"

#############edgeR 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "Sample_Type", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)

#add category 
res$category <- ifelse(res$logFC>0, 1, 2)

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, Sample_Type=="BAL")
ps.rel.2 <- subset_samples(ps.rel, Sample_Type=="Supraglottic")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[res$taxa]
res$abundance.1 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[res$taxa]
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
res <- inner_join(res, contamlist, by="taxa")
#add contam color 
res$color <- ifelse(res$contaminant =="TRUE", "red", "black")


#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(taxa, taxid, superkingdom, kingdom, genus,
                logFC, PValue, FDR, abundance.1, abundance.2,
                category, contaminant) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "metatranscriptomics/Results/edgeR.results_BAL_vs_Sup.csv")

#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance.1, 
                        res$abundance.2)

#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"blue", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "goldenrod","darkgrey"))

#add contam color 
res.bubble.filtered$taxa_col <- res.bubble.filtered$taxa
res.bubble.filtered$taxa_col <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$taxa_col, "</span>")

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))

#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(taxa_col, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(taxa_col, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c( "blue", "goldenrod"))+
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

pdf(file="metatranscriptomics/Figures/EdgeR_bubble_BAL_vs_Sup.pdf", height = 12, width = 14)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="metatranscriptomics/Figures/EdgeR_bubble_BAL_vs_Sup_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 






######## pruning taxa to leave out low abundant taxa. aim to around 1000 taxa ##### 

############### prune all samples table first

#prune the data and use them for all analysis that combine both types of samples later 
pruned.taxa <- genefilter_sample(physeq.rel, filterfun_sample(function(x) x> 0.0002), A = 0.0005 * nsamples(physeq.rel))

#prune relative table 
physeq.rel.pruned <- prune_taxa(pruned.taxa, physeq.rel)

physeq.rel.pruned

#see results 
plot_bar(physeq.rel.pruned)

#prune count table 
physeq.pruned <- prune_taxa(pruned.taxa, physeq)
#see results 
physeq.pruned

#create BAL only table and one sample per subject from pruned data 
temp_physeq <- subset_samples(physeq = physeq.pruned, Sample_Type=="BAL")
BAL.physeq <- subset_samples(physeq = temp_physeq, prim_sec_sample=="prim")

#create relatvie table 
BAL.physeq.rel <- transform_sample_counts(BAL.physeq, normalizeSample)

#rename the table so you know it's pruned 
BAL.physeq.pruned <- BAL.physeq

#prune count table 
BAL.physeq.rel.pruned <- BAL.physeq.rel










###############################################################################
###############################################################################
###############################################################################
##############comparison of clinical groups ###################################
###############################################################################
###############################################################################


############## PD vs. non PD comparison##########

ps <- BAL.physeq.pruned
ps.rel <- BAL.physeq.rel.pruned


#column: ultimate_response_2_lev
######ddPCR
ddpcr_dat <- data.frame(ps@sam_data)

ddpcr_dat <- ddpcr_dat %>% dplyr::mutate(ddPCR=as.numeric(ddPCR))

#calculate stats for ddpcr
my.comparisons <- compare_means(ddPCR~ultimate_response_2_lev, data = ddpcr_dat)
#save it
write.csv(my.comparisons, file = "metatranscriptomics/Results/ddPCR.statsultimate_response_2_lev.csv")

ddpcr_dat <- ddpcr_dat %>% 
  mutate(ultimate_response_2_lev= factor(ultimate_response_2_lev, levels=c("non_PD", "PD")))

# plot it
pdf(file = "metatranscriptomics/Figures/ddPCR_ultimate_response_2_lev.pdf", width=6, height=9)
ggplot(ddpcr_dat, aes(x=ultimate_response_2_lev, y=ddPCR, fill=ultimate_response_2_lev))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("blue","red"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("16S rRNA copies/uL")+
  scale_x_discrete(labels = c("Non PD", "PD"))+
  scale_y_continuous(trans = "log10")+
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
write.csv(my.comparisons, file = "metatranscriptomics/Results/Shannon.stats_ultimate_response_2_lev.csv")

# plot it
pdf(file = "metatranscriptomics/Figures/Shannon_ultimate_response_2_lev.pdf", width=6, height=9)
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
pdf(file = "metatranscriptomics/Figures/observed_ultimate_response_2_lev.pdf", width=6, height=9)
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

pdf(file = "metatranscriptomics/Figures/Beta.Diversity.Brayultimate_response_2_lev.pdf", 
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
  ggtitle("Beta Diversity, Bray, ultimate_response_2_lev, p=0.017")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#sats
write.table(adonis2(vegdist ~ newResults$ultimate_response_2_lev), 
            file = "metatranscriptomics/Results/Beta.Diversity.Brayultimate_response_2_lev.txt", sep = "\t", row.names = T)




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
res$category_2 <- ifelse(res$category=="1", "PD", "non PD")

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, ultimate_response_2_lev=="PD")
ps.rel.2 <- subset_samples(ps.rel, ultimate_response_2_lev=="non_PD")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[res$taxa]
res$abundance.1 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[res$taxa]
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
res <- inner_join(res, contamlist, by="taxa")
#add contam color 
res$color <- ifelse(res$contaminant =="TRUE", "red", "black")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(taxa, taxid, superkingdom, kingdom, genus,
                logFC, PValue, FDR, abundance.1, abundance.2,
                category, category_2, contaminant) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "metatranscriptomics/Results/edgeR.results_ultimate_response_2_lev_PD_vs_non_PD.csv")

#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance.1, 
                        res$abundance.2)

#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
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
res.bubble.filtered$taxa_col <- res.bubble.filtered$taxa
res.bubble.filtered$taxa_col <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$taxa_col, "</span>")

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))

#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(taxa_col, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(taxa_col, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
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

pdf(file="metatranscriptomics/Figures/EdgeR_bubble_ultimate_response_2_lev.pdf", height = 11, width = 10)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="metatranscriptomics/Figures/EdgeR_bubble_ultimate_response_2_lev_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 



###### add another figure for publication with removal of some taxa 
res.bubble.filtered_2 <- res.bubble.filtered[-27,]
res.bubble.filtered_2 <- res.bubble.filtered_2 %>% filter(taxa!="Mus musculus mobilized endogenous polytropic provirus")



bubble.plot <- ggplot(res.bubble.filtered_2, mapping = aes(x=logFC, y=fct_reorder(taxa_col, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered_2[res.bubble.filtered_2$FDR < 0.2,],
               aes(yend=fct_reorder(taxa_col, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
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

pdf(file="metatranscriptomics/Figures/EdgeR_bubble_ultimate_response_2_lev_for_publication.pdf", height = 11, width = 8)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="metatranscriptomics/Figures/EdgeR_bubble_ultimate_response_2_lev_for_publication_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 




######### repeat edgeR while adjusting for coariates 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "ultimate_response_2_lev", method = "TMM")

temp <- data.frame(sample_data(ps))
#group <- "ultimate_response_2_lev"
# Convert categorical variables to factors
temp$Sex <- as.factor(temp$Sex)
temp$race_2 <- ifelse(temp$Race=="White", "White", "others")
temp$race_2 <- as.factor(temp$race_2)
temp$Smoking_2 <- ifelse(temp$Smoking=="Never", "Never", "Smoker")
temp$Smoking_2 <- factor(temp$Smoking_2, levels = c("Never", "Smoker"))
temp$stage_at_bronch <- as.factor(temp$stage_at_bronch)
temp$stage_2_lev <- ifelse(temp$stage_at_bronch%in%c("IA", "IA2", "IIA", "IIB", "IIIA"), "early", "advanced")
temp$stage_2_lev <- factor(temp$stage_2_lev, levels = c("early", "advanced"))
temp$Driver_mutation_present <- as.factor(temp$Driver_mutation_present)

# Ensure continuous variables are numeric
temp$Age <- as.numeric(temp$Age)

# Create the design matrix
design <- model.matrix(~ ultimate_response_2_lev + Sex + Age + race_2 + Smoking_2 + stage_2_lev + Driver_mutation_present, data = temp)

# View the design matrix
print(design)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit the model
fit <- glmFit(dge, design)

# Perform likelihood ratio test for the main factor (e.g., `group`)
lrt <- glmLRT(fit, coef = 2)  # Coef 2 corresponds to the group variable

#create results table 
#run test 
#ET <- exactTest(dge)

#get results 
TT <- topTags(lrt,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

############continue edgeR here


res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)

#add category 
res$category <- ifelse(res$logFC>0, 1, 2)

#add category detailed 
res$category_2 <- ifelse(res$category=="1", "PD", "non PD")

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, ultimate_response_2_lev=="PD")
ps.rel.2 <- subset_samples(ps.rel, ultimate_response_2_lev=="non_PD")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[res$taxa]
res$abundance.1 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[res$taxa]
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
res <- inner_join(res, contamlist, by="taxa")
#add contam color 
res$color <- ifelse(res$contaminant =="TRUE", "red", "black")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(taxa, taxid, superkingdom, kingdom, genus,
                logFC, PValue, FDR, abundance.1, abundance.2,
                category, category_2, contaminant) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "metatranscriptomics/Results/edgeR.results_ultimate_response_2_lev_PD_vs_non_PD_adjusted_for_covariates.csv")

#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance.1, 
                        res$abundance.2)

#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
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
res.bubble.filtered$taxa_col <- res.bubble.filtered$taxa
res.bubble.filtered$taxa_col <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$taxa_col, "</span>")

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))

#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(taxa_col, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(taxa_col, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
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

pdf(file="metatranscriptomics/Figures/EdgeR_bubble_ultimate_response_2_lev_adjusted_for_covariates.pdf", height = 11, width = 14)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="metatranscriptomics/Figures/EdgeR_bubble_ultimate_response_2_lev_adjusted_for_covariates_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 





######### two years mortality ############

BAL.physeq.pruned@sam_data

ps <- subset_samples(BAL.physeq.pruned,  two_y_mort!="NA")
ps.rel <- subset_samples(BAL.physeq.rel.pruned,  two_y_mort!="NA")

#column: two_y_mort
######ddPCR
ddpcr_dat <- data.frame(ps@sam_data)

ddpcr_dat <- ddpcr_dat %>% dplyr::mutate(ddPCR=as.numeric(ddPCR))

#calculate stats for ddpcr
my.comparisons <- compare_means(ddPCR~two_y_mort, data = ddpcr_dat)
#save it
write.csv(my.comparisons, file = "metatranscriptomics/Results/ddPCR.stats.two_y_mort.csv")

ddpcr_dat <- ddpcr_dat %>% 
  mutate(two_y_mort= factor(two_y_mort, levels=c("greater_24", "less_24")))

# plot it
pdf(file = "metatranscriptomics/Figures/ddPCR_two_y_mort.pdf", width=6, height=9)
ggplot(ddpcr_dat, aes(x=two_y_mort, y=ddPCR, fill=two_y_mort))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("orange3","red3"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("16S rRNA copies/uL")+
  scale_x_discrete(labels = c(">24 Months", "<24 Months"))+
  scale_y_continuous(trans = "log10")+
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
write.csv(my.comparisons, file = "metatranscriptomics/Results/Shannon.stats_two_y_mort.csv")

# plot it
pdf(file = "metatranscriptomics/Figures/Shannon_two_y_mort.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=two_y_mort, y=Shannon, fill=two_y_mort))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("orange3","red3"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c(">24 Months", "<24 Months"))+
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
pdf(file = "metatranscriptomics/Figures/observed_two_y_mort.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=two_y_mort, y=Observed, fill=two_y_mort))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("orange3","red3"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c(">24 Months", "<24 Months"))+
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

pdf(file = "metatranscriptomics/Figures/Beta.Diversity.Bray.two_y_mort.pdf", 
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
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c(">24 Months", "<24 Months")), size=14) +
  ggtitle("Beta Diversity, Bray, two_y_mort, p=0.004")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#sats
write.table(adonis2(vegdist ~ newResults$two_y_mort), 
            file = "metatranscriptomics/Results/Beta.Diversity.Bray.two_y_mort.txt", sep = "\t", row.names = T)




#############edgeR 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "two_y_mort", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)

#add category 
res$category <- ifelse(res$logFC>0, 1, 2)

#add category detailed 
res$category_2 <- ifelse(res$category=="1", "less_24", "greater_24")

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, two_y_mort=="less_24")
ps.rel.2 <- subset_samples(ps.rel, two_y_mort=="greater_24")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[res$taxa]
res$abundance.1 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[res$taxa]
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
res <- inner_join(res, contamlist, by="taxa")
#add contam color 
res$color <- ifelse(res$contaminant =="TRUE", "red", "black")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(taxa, taxid, superkingdom, kingdom, genus,
                logFC, PValue, FDR, abundance.1, abundance.2,
                category, category_2, contaminant) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "metatranscriptomics/Results/edgeR.results_two_y_mort_high_vs_low.csv")

#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance.1, 
                        res$abundance.2)

#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
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
res.bubble.filtered$taxa_col <- res.bubble.filtered$taxa
res.bubble.filtered$taxa_col <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$taxa_col, "</span>")

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))

#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(taxa_col, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(taxa_col, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
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

pdf(file="metatranscriptomics/Figures/EdgeR_bubble_two_y_mort.pdf", height = 11, width = 14)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="metatranscriptomics/Figures/EdgeR_bubble_two_y_mort_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 


######### repeat edgeR while adjusting for coariates 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "two_y_mort", method = "TMM")

temp <- data.frame(sample_data(ps))
#group <- "ultimate_response_2_lev"
# Convert categorical variables to factors
temp$Sex <- as.factor(temp$Sex)
temp$race_2 <- ifelse(temp$Race=="White", "White", "others")
temp$race_2 <- as.factor(temp$race_2)
temp$Smoking_2 <- ifelse(temp$Smoking=="Never", "Never", "Smoker")
temp$Smoking_2 <- factor(temp$Smoking_2, levels = c("Never", "Smoker"))
temp$stage_at_bronch <- as.factor(temp$stage_at_bronch)
temp$stage_2_lev <- ifelse(temp$stage_at_bronch%in%c("IA", "IA2", "IIA", "IIB", "IIIA"), "early", "advanced")
temp$stage_2_lev <- factor(temp$stage_2_lev, levels = c("early", "advanced"))
temp$Driver_mutation_present <- as.factor(temp$Driver_mutation_present)

# Ensure continuous variables are numeric
temp$Age <- as.numeric(temp$Age)

# Create the design matrix
design <- model.matrix(~ ultimate_response_2_lev + Sex + Age + race_2 + Smoking_2 + stage_2_lev + Driver_mutation_present, data = temp)


# View the design matrix
print(design)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit the model
fit <- glmFit(dge, design)

# Perform likelihood ratio test for the main factor (e.g., `group`)
lrt <- glmLRT(fit, coef = 2)  # Coef 2 corresponds to the group variable

#create results table 
#run test 
#ET <- exactTest(dge)

#get results 
TT <- topTags(lrt,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

############continue edgeR here


res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)

#add category 
res$category <- ifelse(res$logFC>0, 1, 2)

#add category detailed 
res$category_2 <- ifelse(res$category=="1", "less_24", "greater_24")

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, two_y_mort=="less_24")
ps.rel.2 <- subset_samples(ps.rel, two_y_mort=="greater_24")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[res$taxa]
res$abundance.1 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[res$taxa]
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
res <- inner_join(res, contamlist, by="taxa")
#add contam color 
res$color <- ifelse(res$contaminant =="TRUE", "red", "black")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(taxa, taxid, superkingdom, kingdom, genus,
                logFC, PValue, FDR, abundance.1, abundance.2,
                category, category_2, contaminant) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "metatranscriptomics/Results/edgeR.results_two_y_mort_adjusted_for_covariates.csv")

#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance.1, 
                        res$abundance.2)

#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
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
res.bubble.filtered$taxa_col <- res.bubble.filtered$taxa
res.bubble.filtered$taxa_col <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$taxa_col, "</span>")

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))

#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(taxa_col, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(taxa_col, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
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

pdf(file="metatranscriptomics/Figures/EdgeR_bubble_two_y_mort_adjusted_for_covariates.pdf", height = 11, width = 14)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="metatranscriptomics/Figures/EdgeR_bubble_two_y_mort_adjusted_for_covariates_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 




######## k means clustering, all samples#####
## Create scree plot for k-means clustering

ps.rel <- physeq.rel

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

mds.df <-as.data.frame(cbind(CmdScale[,1], CmdScale[,2]))

n <- 10
wss <- numeric(n)

set.seed(123)
for (j in 1:n) {
  # Fit the model: km.out
  km.out <- kmeans(mds.df, centers = j, nstart = 20)
  # Save the within cluster sum of squares
  wss[j] <- km.out$tot.withinss
}

wss.df <- tibble(clusters = 1:n, wss = wss)

scree_plot <- ggplot(wss.df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  #labs(title=names(samples)[i],x="number of clusters")+
  theme_bw()
plot(scree_plot)
pdf(file = "metatranscriptomics/Figures/k_means_scree_plot_unpruned_data.pdf", height = 8, width = 10)
scree_plot
dev.off()



## Run k-means clustering on Bray-Curtis distances

plots.list <- list()

k_clusters <- list()

for(j in 2:4){
  set.seed(123)
  km.out <- kmeans(mds.df, centers = j, nstart = 20)
  
  k_clusters[[j]] <- km.out$cluster
  
  plots.list[[j]] <- ggplot()+
    geom_point(aes(x=!!mds.df[,1],y=!!mds.df[,2],color=!!as.character(km.out$cluster)),
               size=2)+
    stat_ellipse(aes(x=!!mds.df[,1],y=!!mds.df[,2],color=!!as.character(km.out$cluster),fill=!!as.character(km.out$cluster)),
                 geom="polygon",
                 alpha=0.2)+
    labs(title=paste0(" k = ",j),
         x=paste0("PC1: ",percentVar[1],"% variance"),
         y=paste0("PC2: ",percentVar[2],"% variance"))+
    guides(color="none",fill="none")+
    theme_bw()
  
  plots.list <- plots.list[-which(sapply(plots.list, is.null))]
}

png(file = "metatranscriptomics/Figures/k_means_clustering_plots_k2_4_pruned.png", width = 800, height = 700)
egg::ggarrange(plots = plots.list, nrow=1)
dev.off()


##### choose  number of clusters and save. will repeat for 2, 3 and 4 clusters 
km.out <- kmeans(mds.df, centers = 2)

sample_data(physeq.rel)$k_means_2_cluster <- km.out$cluster
sample_data(physeq)$k_means_2_cluster <- km.out$cluster
sample_data(physeq.rel)$k_means_2_cluster <- sample_data(physeq)$k_means_2_cluster

#addd three clusrers 
km.out <- kmeans(mds.df, centers = 3)

sample_data(physeq.rel)$c <- km.out$cluster
sample_data(physeq)$k_means_3_cluster <- km.out$cluster
sample_data(physeq.rel)$k_means_3_cluster <- sample_data(physeq)$k_means_3_cluster

## add 4 clusters 
km.out <- kmeans(mds.df, centers = 4)

sample_data(physeq.rel)$k_means_4_cluster <- km.out$cluster
sample_data(physeq)$k_means_4_cluster <- km.out$cluster
sample_data(physeq.rel)$k_means_4_cluster <- sample_data(physeq)$k_means_4_cluster

#### add all to BAL table 
temp <- data.frame(sample_data(physeq))
temp <- temp %>% filter(Sample_Type=="BAL")
BAL.physeq.pruned@sam_data$k_means_2_cluster <- temp$k_means_2_cluster
BAL.physeq.pruned@sam_data$k_means_3_cluster <- temp$k_means_3_cluster
BAL.physeq.pruned@sam_data$k_means_4_cluster <- temp$k_means_4_cluster
BAL.physeq.rel.pruned@sam_data$k_means_2_cluster <- temp$k_means_2_cluster
BAL.physeq.rel.pruned@sam_data$k_means_3_cluster <- temp$k_means_3_cluster
BAL.physeq.rel.pruned@sam_data$k_means_4_cluster <- temp$k_means_4_cluster


########### 3 clusters best fit the data ##########

###### three clusters analysis######

ps <- physeq
ps.rel <- physeq.rel

#one year mortality  comparison 

#use data - sample_info 
sample_info <- data.frame(sample_data(ps))
#search for columns to keep and assign only lower airway samples 
sample_info <- sample_info %>% filter(Sample_Type=="BAL")

#death give mortality status 
sample_info$one_y_mort
sample_info <- sample_info %>% filter(one_y_mort!="NA")
sample_info$one_y_time <- ifelse(sample_info$one_y_mort=="greater_12", "365", sample_info$OS_days)
sample_info$one_y_time <- as.numeric(sample_info$one_y_time)
survival.data <- sample_info[, c("one_y_time","one_y_mort", "k_means_3_cluster")]
#rename the new df 
colnames(survival.data) <- c("time", "status", "k_means_3_cluster")

# turn the data into data frame so it can be used in survival function. 
#turn the class of columns so its either factor or numeric to be used in survival functions 
survival.data<- data.frame(survival.data)
survival.data$time <- as.numeric(survival.data$time)
survival.data <- survival.data %>% mutate(status=ifelse(status=="less_12",1, 0))
survival.data$status <- as.numeric(survival.data$status)
survival.data$k_means_3_cluster <- as.factor(survival.data$k_means_3_cluster)
survival.data
#remove NA 
survival.data <- na.omit(survival.data)

#get median survival in each group 
survival.data %>% 
  group_by(k_means_3_cluster) %>% 
  summarise(median=median(time), 
            Q1=quantile(time, 0.25), 
            Q3=quantile(time, 0.75))

#create survival object
library(survival)
library(survminer)
library(survMisc)

surv <- Surv(survival.data$time, survival.data$status)

#create survival df
sfit <- survfit(Surv(time, status)~k_means_3_cluster, data=survival.data)
sfit

survival.plot <- ggsurvplot(sfit, data = survival.data,
                            fun = "pct", ggtheme = theme_pubr(), 
                            conf.int=FALSE, pval=TRUE, risk.table=TRUE, 
                            legend.labs=c("Cluster 1", "Cluster 2", "Cluster 3"), legend.title="K means Cluster",  
                            #palette=c("brown2", "burlywood2", "chocolate2"), 
                            size=1,
                            title="Kaplan-Meier Curve for one year mortaltiy  According to K means Cluster",
                            ylab="Overall mortality Probability",
                            xlab="Time(Days)",
                            # tables.col = "strata", #color of numbers inside risk table
                            tables.height = 0.15,
                            fontsize = 4,
                            risk.table.y.text.col = T,
                            cumcensor = TRUE, 
                            tables.theme = theme_cleantable(),
                            # risk.table.y.text = FALSE,
                            risk.table.title = "Number at risk (cumulative number of recurrence)",
                            cumcensor.title = "Cumulative number of censored subjects"
)

#save 
survival.plot

#save 


#compare beween clusters stats 
# cluster 1 vs 2 p =0.027
subset_data <- survival.data %>% filter(k_means_3_cluster %in% c(1, 2))
surv_obj <- Surv(subset_data$time, subset_data$status)
surv_diff <- survdiff(surv_obj ~ k_means_3_cluster, data = subset_data)
p_value <- 1 - pchisq(surv_diff$chisq, df = 1) # df = 1 for two groups
p_value

# cluster 1 vs 3 p =0.092
subset_data <- survival.data %>% filter(k_means_3_cluster %in% c(1, 3))
surv_obj <- Surv(subset_data$time, subset_data$status)
surv_diff <- survdiff(surv_obj ~ k_means_3_cluster, data = subset_data)
p_value <- 1 - pchisq(surv_diff$chisq, df = 1) # df = 1 for two groups
p_value

# cluster 3 vs 2 p =0.73 
subset_data <- survival.data %>% filter(k_means_3_cluster %in% c(2, 3))
surv_obj <- Surv(subset_data$time, subset_data$status)
surv_diff <- survdiff(surv_obj ~ k_means_3_cluster, data = subset_data)
p_value <- 1 - pchisq(surv_diff$chisq, df = 1) # df = 1 for two groups
p_value


#### two y morta 

#use data - sample_info 
sample_info <- data.frame(sample_data(ps))
#search for columns to keep and assign only lower airway samples 
sample_info <- sample_info %>% filter(Sample_Type=="BAL")

#death give mortality status 
sample_info$two_y_mort
sample_info <- sample_info %>% filter(two_y_mort!="NA")
sample_info$two_y_time <- ifelse(sample_info$two_y_mort=="greater_24", "730.5", sample_info$OS_days)
sample_info$two_y_time <- as.numeric(sample_info$two_y_time)
survival.data <- sample_info[, c("two_y_time","two_y_mort", "k_means_3_cluster")]
#rename the new df 
colnames(survival.data) <- c("time", "status", "k_means_3_cluster")

# turn the data into data frame so it can be used in survival function. 
#turn the class of columns so its either factor or numeric to be used in survival functions 
survival.data<- data.frame(survival.data)
survival.data$time <- as.numeric(survival.data$time)
survival.data <- survival.data %>% mutate(status=ifelse(status=="less_24",1, 0))
survival.data$status <- as.numeric(survival.data$status)
survival.data$k_means_3_cluster <- as.factor(survival.data$k_means_3_cluster)
survival.data
#remove NA 
survival.data <- na.omit(survival.data)

#get median survival in each group 
survival.data %>% 
  group_by(k_means_3_cluster) %>% 
  summarise(median=median(time), 
            Q1=quantile(time, 0.25), 
            Q3=quantile(time, 0.75))

#create survival object
library(survival)
library(survminer)
library(survMisc)

surv <- Surv(survival.data$time, survival.data$status)

#create survival df
sfit <- survfit(Surv(time, status)~k_means_3_cluster, data=survival.data)
sfit

survival.plot <- ggsurvplot(sfit, data = survival.data,
                            fun = "pct", ggtheme = theme_pubr(), 
                            conf.int=FALSE, pval=TRUE, risk.table=TRUE, 
                            legend.labs=c("Cluster 1", "Cluster 2", "Cluster 3"), legend.title="K means Cluster",  
                            #palette=c("brown2", "burlywood2", "chocolate2"), 
                            size=1,
                            title="Kaplan-Meier Curve for two year mortaltiy  According to K means Cluster",
                            ylab="Overall mortality Probability",
                            xlab="Time(Days)",
                            # tables.col = "strata", #color of numbers inside risk table
                            tables.height = 0.15,
                            fontsize = 4,
                            risk.table.y.text.col = T,
                            cumcensor = TRUE, 
                            tables.theme = theme_cleantable(),
                            # risk.table.y.text = FALSE,
                            risk.table.title = "Number at risk (cumulative number of recurrence)",
                            cumcensor.title = "Cumulative number of censored subjects"
)

#save 
survival.plot


#compare beween clusters stats 
# cluster 1 vs 2 p =0.02
subset_data <- survival.data %>% filter(k_means_3_cluster %in% c(1, 2))
surv_obj <- Surv(subset_data$time, subset_data$status)
surv_diff <- survdiff(surv_obj ~ k_means_3_cluster, data = subset_data)
p_value <- 1 - pchisq(surv_diff$chisq, df = 1) # df = 1 for two groups
p_value

# cluster 1 vs 3 p =0.301
subset_data <- survival.data %>% filter(k_means_3_cluster %in% c(1, 3))
surv_obj <- Surv(subset_data$time, subset_data$status)
surv_diff <- survdiff(surv_obj ~ k_means_3_cluster, data = subset_data)
p_value <- 1 - pchisq(surv_diff$chisq, df = 1) # df = 1 for two groups
p_value

# cluster 3 vs 2 p =0.54 
subset_data <- survival.data %>% filter(k_means_3_cluster %in% c(2, 3))
surv_obj <- Surv(subset_data$time, subset_data$status)
surv_diff <- survdiff(surv_obj ~ k_means_3_cluster, data = subset_data)
p_value <- 1 - pchisq(surv_diff$chisq, df = 1) # df = 1 for two groups
p_value



#### PD as outcome 


#use data - sample_info 
sample_info <- data.frame(sample_data(ps))
#search for columns to keep and assign only lower airway samples 
sample_info <- sample_info %>% filter(Sample_Type=="BAL")

#death give mortality status 
sample_info$ultimate_response_2_lev
sample_info$TTP <- sample_info$TTP
sample_info$TTP <- as.numeric(sample_info$TTP)
survival.data <- sample_info[, c("TTP","ultimate_response_2_lev", "k_means_3_cluster")]
#rename the new df 
colnames(survival.data) <- c("time", "status", "k_means_3_cluster")

# turn the data into data frame so it can be used in survival function. 
#turn the class of columns so its either factor or numeric to be used in survival functions 
survival.data<- data.frame(survival.data)
survival.data$time <- as.numeric(survival.data$time)
survival.data <- survival.data %>% mutate(status=ifelse(status=="PD",1, 0))
survival.data$status <- as.numeric(survival.data$status)
survival.data$k_means_3_cluster <- as.factor(survival.data$k_means_3_cluster)
survival.data
#remove NA 
survival.data <- na.omit(survival.data)

#get median survival in each group 
survival.data %>% 
  group_by(k_means_3_cluster) %>% 
  summarise(median=median(time), 
            Q1=quantile(time, 0.25), 
            Q3=quantile(time, 0.75))

#create survival object
library(survival)
library(survminer)
library(survMisc)

surv <- Surv(survival.data$time, survival.data$status)

#create survival df
sfit <- survfit(Surv(time, status)~k_means_3_cluster, data=survival.data)
sfit

survival.plot <- ggsurvplot(sfit, data = survival.data,
                            fun = "pct", ggtheme = theme_pubr(), 
                            conf.int=FALSE, pval=TRUE, risk.table=TRUE, 
                            legend.labs=c("Cluster 1", "Cluster 2", "Cluster 3"), legend.title="K means Cluster",  
                            #palette=c("brown2", "burlywood2", "chocolate2"), 
                            size=1,
                            title="Kaplan-Meier Curve for two year mortaltiy  According to K means Cluster",
                            ylab="Overall mortality Probability",
                            xlab="Time(Days)",
                            # tables.col = "strata", #color of numbers inside risk table
                            tables.height = 0.15,
                            fontsize = 4,
                            risk.table.y.text.col = T,
                            cumcensor = TRUE, 
                            tables.theme = theme_cleantable(),
                            # risk.table.y.text = FALSE,
                            risk.table.title = "Number at risk (cumulative number of recurrence)",
                            cumcensor.title = "Cumulative number of censored subjects"
)

#save 
survival.plot


#compare beween clusters stats 
# cluster 1 vs 2 p =0.1
subset_data <- survival.data %>% filter(k_means_3_cluster %in% c(1, 2))
surv_obj <- Surv(subset_data$time, subset_data$status)
surv_diff <- survdiff(surv_obj ~ k_means_3_cluster, data = subset_data)
p_value <- 1 - pchisq(surv_diff$chisq, df = 1) # df = 1 for two groups
p_value

# cluster 1 vs 2 p =0.06
subset_data <- survival.data %>% filter(k_means_3_cluster %in% c(1, 3))
surv_obj <- Surv(subset_data$time, subset_data$status)
surv_diff <- survdiff(surv_obj ~ k_means_3_cluster, data = subset_data)
p_value <- 1 - pchisq(surv_diff$chisq, df = 1) # df = 1 for two groups
p_value

# cluster 1 vs 2 p =0.2389217 
subset_data <- survival.data %>% filter(k_means_3_cluster %in% c(2, 3))
surv_obj <- Surv(subset_data$time, subset_data$status)
surv_diff <- survdiff(surv_obj ~ k_means_3_cluster, data = subset_data)
p_value <- 1 - pchisq(surv_diff$chisq, df = 1) # df = 1 for two groups
p_value


#### beta diversity according to sample type with 3 clusters 


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
centroids <- aggregate(cbind(PC1,PC2)~ Sample_Type,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Sample_Type",suffixes=c("",".centroid"))

newResults$k_means_3_cluster <- factor(newResults$k_means_3_cluster)
centroids$Sample_Type <- factor(centroids$Sample_Type, levels = c("BKG", "Supraglottic", "BAL"))
newResults$Sample_Type <- factor(newResults$Sample_Type, levels = c("BKG", "Supraglottic", "BAL"))

#stats
x <- adonis2(vegdist ~ newResults$k_means_3_cluster)

pdf(file = "metatranscriptomics/Figures/Beta.Diversity.Bray.k_means_3_clusters_and_Sample_type_al_samples_unpruned.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= k_means_3_cluster)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("#F8766D", "#00ba38", "#619cff", "grey","goldenrod", "blue")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Sample_Type), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Sample_Type)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("BAL", "BKG", "UA"), color=Sample_Type), size=20) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

##### between clusters p value 

#1 vs 2 (p=0.001)
ps.re.1 <- subset_samples(ps.rel, k_means_3_cluster %in% c("1", "2"))
vegdist   = vegdist(t(otu_table(ps.re.1)), method = "bray")
adonis2(vegdist~sample_data(ps.re.1)$k_means_3_cluster)

#2 vs 3 (p=0.001)
ps.re.1 <- subset_samples(ps.rel, k_means_3_cluster %in% c("1", "3"))
vegdist   = vegdist(t(otu_table(ps.re.1)), method = "bray")
adonis2(vegdist~sample_data(ps.re.1)$k_means_3_cluster)

#1 vs 3 (p=0.001)
ps.re.1 <- subset_samples(ps.rel, k_means_3_cluster %in% c("3", "2"))
vegdist   = vegdist(t(otu_table(ps.re.1)), method = "bray")
adonis2(vegdist~sample_data(ps.re.1)$k_means_3_cluster)





####### alpha diversity 
#calculate alpha diversity 
alpha.measures <- data.frame("Shannon" = estimate_richness(ps, measures = "Shannon"), 
                             "Observed" = estimate_richness(ps, measures = "Observed"), 
                             "Simpson" = estimate_richness(ps, measures = "Simpson"),
                             "k_means_3_cluster" = sample_data(ps)$k_means_3_cluster, 
                             "Sample_Type"=sample_data(ps)$Sample_Type)
alpha.measures <- alpha.measures %>% 
  mutate(k_means_3_cluster= factor(k_means_3_cluster))

# plot it
pdf(file = "metatranscriptomics/Figures/Shannon_NAT_k_means_3_cluster_all_samples.pdf", width=6, height=9)
ggplot(alpha.measures, aes(x=k_means_3_cluster, y=Shannon, fill=k_means_3_cluster))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  #scale_fill_manual(values = c("#F8766D","#00bfc4"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c("Cluster 1", "Cluster 2", "Cluster 3"))+
  stat_compare_means(comparisons = list(c("1", "2"), c("2", "3"), c("1", "3")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()


#facet by sample type 
pdf(file = "metatranscriptomics/Figures/Shannon_NAT_k_means_3_cluster_all_samples_facet.pdf", width=11, height=9)
ggplot(alpha.measures, aes(x=k_means_3_cluster, y=Shannon, fill=k_means_3_cluster))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  #scale_fill_manual(values = c("#F8766D","#00bfc4"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c("Cluster 1", "Cluster 2", "Cluster 3"))+
  stat_compare_means(comparisons = list(c("1", "2"), c("2", "3"), c("1", "3")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  facet_wrap(~Sample_Type)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

#alpha for BAL samples only 
alpha.measures_bal <- alpha.measures %>% filter(Sample_Type=="BAL")
pdf(file = "metatranscriptomics/Figures/Shannon_NAT_k_means_3_cluster_BAL_samples_only.pdf", width=6, height=9)
ggplot(alpha.measures_bal, aes(x=k_means_3_cluster, y=Shannon, fill=k_means_3_cluster))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  #scale_fill_manual(values = c("#F8766D","#00bfc4"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c("Cluster 1", "Cluster 2", "Cluster 3"))+
  stat_compare_means(comparisons = list(c("1", "2"), c("2", "3"), c("1", "3")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()





########## lefse analysi for 3 clusters 

# set thresholds for LefSe 
pvalue.thres <- 0.05
lda.thres <- 2


ps <- BAL.physeq.pruned
ps.rel <- BAL.physeq.rel.pruned

#set columns names according to the package requirements 
taxa_table_original <- data.frame(tax_table(ps))
taxa_table_original$Kingdom <- taxa_table_original$kingdom
taxa_table_original$Phylum <- "<NA>"
taxa_table_original$Class <- "<NA>"
taxa_table_original$Order <- "<NA>"
taxa_table_original$Family <- "<NA>"
taxa_table_original$Genus <- taxa_table_original$genus
taxa_table_original$Species <- taxa_table_original$taxa
taxa_table_lefse <- taxa_table_original %>%dplyr::select(-c(taxa, taxid, genus, kingdom, superkingdom, clade, family))

#check, should be: colnames(tax_table(ps)) = c("Kingdom", "Phylum","Class","Order","Family", "Genus", "Species")
colnames(taxa_table_lefse)
tax_matrix <- as.matrix(taxa_table_lefse)
tax_table(ps) <- tax_matrix
tax_table(ps.rel) <- tax_matrix

#set columns names according to the package requirements 
colnames(tax_table(ps)) = c("Kingdom", "Phylum","Class","Order","Family", "Genus", "Species")
colnames(tax_table(ps.rel))=c("Kingdom", "Phylum","Class","Order","Family", "Genus", "Species")

ps@sam_data$k_means_3_cluster
#Run LefSe using microbiomeMarker package
lefse.res <- run_lefse(ps, 
                       group = "k_means_3_cluster",  #the column choosed to differntiate groups
                       subgroup = NULL, 
                       taxa_rank = "none", 
                       transform = "identity", 
                       norm = "CPM", 
                       kw_cutoff = 0.05,              
                       wilcoxon_cutoff = 0.05, 
                       lda_cutoff = 2, 
                       bootstrap_n = 30, 
                       bootstrap_fraction = 2/3,
                       multigrp_strat = TRUE)

#inspect results 
lefse.res
lefse.res@marker_table

#make marker table as a data.frame and store the results
res1 <- data.frame(marker_table(lefse.res))
#add column of marker no. and row no. 
res1 <- res1 %>% 
  mutate(marker = rownames(.)) %>% 
  mutate(number = c(1:nrow(.)))

res1$taxa <- res1$feature

#prepeare contam list 
contamlist 

#now merge the table with your results according to asv column 
merged <- inner_join(res1, contamlist, by="taxa")

merged <- merged %>% 
  mutate(contam_color = ifelse(contaminant=="TRUE", "red", "black")) %>% 
  mutate(FDR = p.adjust(merged$pvalue, method = "BH")) %>%
  mutate(enrich_group = factor(enrich_group, levels = unique(enrich_group))) %>% 
  mutate(feature = factor(feature, levels = unique(merged$feature))) %>% 
  arrange(enrich_group, desc(ef_lda)) %>% 
  mutate(ord= order(enrich_group), ef_lda)

#add color of potential contaminants 
merged$names.colored <- paste0("<span style=\"color: ", merged$contam_color, "\">", merged$feature, "</span>")

write.csv(merged, file = "metatranscriptomics/Results/lefse_results_k_means_3_cluster.csv")
#if needed minimze number to top 10 features of each group 

#1 leave only significant ones 
plot_dat <- merged %>% filter(FDR <=0.2)

sig_1 <- plot_dat %>% filter(enrich_group=="1") %>% dplyr::slice(1:20)
sig_2 <- plot_dat %>% filter(enrich_group=="2") %>% dplyr::slice(1:20)
sig_3 <- plot_dat %>% filter(enrich_group=="3") %>% dplyr::slice(1:20)

#now bind again 
plot.df <- bind_rows(sig_1, sig_2, sig_3)
plot.df <- plot.df %>% mutate(ord= order(enrich_group), ef_lda)

#barplot 
p1 <- ggplot(plot.df, aes(x=ef_lda, y=fct_reorder(names.colored,ord, .desc = T), fill=enrich_group), alpha=0.8)+
  geom_bar(stat= "identity", position = position_dodge(0.9))+
  scale_fill_manual(name = "Enriched Groups",
                    values=c("1"="#F8766D","2"="#00BA38", "3"="#619CFF"),
                    labels=c("Cluster 1","Cluster 2", "Cluster 3"))+
  guides(fill= guide_legend(reverse=TRUE))+
  theme_pubr()+
  theme(
    plot.title = element_text(face="bold"),
    axis.text.y = element_markdown(angle=0,hjust=1, face = "bold.italic", size = 16),
    axis.text.x = element_text(face = "bold", size = 16), 
    axis.title.x = element_text(face = "bold", size = 20))+
  xlab("LDA (log10)")+
  xlim(0,ceiling(max(merged$ef_lda)))+
  ylab("")

#save it individually 
pdf(file = "metatranscriptomics/Figures/lefse_barplot_k_means_3_clusters_top_20.pdf", height = 18, width = 12)
p1
dev.off()

#abundance figure 
#get sample data 
sample.data <- data.frame(sample_data(ps.rel)) %>% 
  dplyr::select(k_means_3_cluster) %>% 
  mutate(k_means_3_cluster = factor(k_means_3_cluster))

#get abundance data 
abund <- data.frame(otu_table(ps.rel))
abund.select <- subset(abund, row.names(abund)%in%plot.df$feature)

#transpose 
abund.select.transposed <- data.frame(t(abund.select))

#now merge the metadata with the abundance table 
abund.select.transposed$k_means_3_cluster <- sample.data$k_means_3_cluster

#organize test1 
test1 <- data.frame(t(abund.select.transposed)) %>% 
  filter(!row_number() %in% 1)
names(test1) <- test1[nrow(test1),]
test1 <- test1[-nrow(test1),]
test1 <- data.frame(sapply(test1[,1:ncol(test1)],as.numeric))
row.names(test1) <- colnames(abund.select.transposed[,2:(ncol(abund.select.transposed)-1)])
#create asv column in test1 and merge with your resutls
test1 <- test1 %>% 
  mutate(tax_name = rownames(test1)) 

#fix taxa names in test1 so they match plot.df$feature
tax_name <- plot.df$feature
tax_name
tax_name <- gsub(" ", ".", tax_name)
tax_name <- gsub("-", ".", tax_name)
plot.df$tax_name <- tax_name

tax_name
rownames(test1)

test1 <- test1%>% 
  inner_join(., plot.df, by="tax_name") %>% 
  mutate(enrich_group = factor(enrich_group)) %>% 
  arrange(enrich_group, desc(ef_lda)) %>% 
  dplyr::select(-enrich_group, -ef_lda, -pvalue, -padj, -marker, -FDR, -number, -taxa, -tax_name, -contaminant,
                -contam_color, -ord,taxa, -names.colored, -taxa)
#transpose DF 
test1 <- data.frame(t(test1))
#set names accoridng to taxa names 
names(test1) <- test1[nrow(test1),]
#remove names row 
test1 <- test1[-nrow(test1),]
#asign categories 
test1$category <- rownames(test1)

# Extract numbers before the dot (if any) and convert to integer
test1$category <- as.integer(sub("X([0-9]+).*", "\\1", test1$category))
test1 <- test1 %>% 
  mutate(category=case_when(
    category== "1"~"1", 
    category== "2"~"2", 
    category== "3"~"3"
  ))

test1$category <- factor(test1$category, levels=c("1","2", "3"))

##### mean abundance plot dataframe
library(dplyr)
a_genus <- pivot_longer(test1, cols = -category , names_to = "id", values_to = "value") %>% 
  as_tibble(.) %>% 
  group_by(id, category) %>% 
  mutate(value = as.numeric(value)) %>% 
  mutate(id = factor(id, levels = unique(id))) %>% 
  dplyr::mutate(count = n()) %>% 
  dplyr::mutate(mean = mean(value)) %>% 
  dplyr::mutate(sd = sd(value)) %>%
  dplyr::mutate(se = sd/sqrt(n())) %>%
  dplyr::mutate(ic=se * qt((1-0.05)/2 + .5, count-1 )) %>%
  mutate(meanpercent = mean * 100) %>%
  ungroup() %>%
  #distinct(.,id,category,count,mean,sd,se,ic, .keep_all= TRUE) %>% 
  arrange(category) %>% 
  mutate(ord = order(category))
a_genus



#relative abundance plot 
p2 <- a_genus %>%
  ggplot(., aes(x = meanpercent, y = fct_reorder(id, -ord), fill=category)) + 
  geom_bar(stat="identity", position=position_dodge(0.9), width = 0.9)+
  geom_bar(stat="identity", position=position_dodge(0.9), width = 0.9, show.legend  =FALSE, size = 0.2)+
  geom_errorbar( aes(xmin=(mean-se)*100, xmax=(mean+se)*100, group = category), width=0.4, colour="black", alpha=0.9, size=.25, position = position_dodge(0.9))+
  scale_fill_manual(name = "Enriched Groups",
                    values=c("1"="#F8766D","2"="#00BA38", "3"="#619CFF"),
                    labels=c("Cluster 1","Cluster 2", "Cluster 3"))+
  guides(fill= guide_legend(reverse=TRUE))+
  theme_classic()+
  theme(
    plot.title = element_text(face="bold"),
    axis.text.y = element_markdown(angle=0,hjust=1, face = "bold.italic", size = 16),
    axis.text.x = element_text(face = "bold", size = 16), 
    axis.title.x = element_text(face = "bold", size = 20))+
  xlab("Relative Abundance (%)")+
  ylab("")+
  ggtitle("")
#coord_cartesian(xlim = c(1, 100))  # Set limits from 1 to 100

p2

#combine two plots 
p3<- ggarrange(p1 + geom_vline(xintercept=lda.thres,linetype="longdash",color="black",size =0.5),
               p2 + 
                 theme(axis.text.y=element_blank(),
                       axis.title=element_text(face="bold"),
                       axis.text.x=element_text(colour="black", face="bold"),
                       axis.ticks=element_line(colour="black"),
                       legend.background = element_rect(color=NA),
                       legend.position = "top"),
               ncol=2,
               nrow=1,
               align = "h",
               labels =c("A","B"),
               common.legend=TRUE,
               legend = "bottom",
               widths = c(1,1)
)

p3

pdf(file = "metatranscriptomics/Figures/lefse_full_k_means_3_clusters_top_20.pdf", height = 22, width = 22)
p3
dev.off()

##################### produce another figure with top 10 taxa only for publication
#1 leave only significant ones 
plot_dat <- merged %>% filter(FDR <=0.2)

plot_dat <- plot_dat%>% dplyr::filter(feature!="synthetic construct")
plot_dat <- plot_dat%>% dplyr::filter(feature!="Equine infectious anemia virus")


sig_1 <- plot_dat %>% filter(enrich_group=="1") %>% dplyr::slice(1:10)
sig_2 <- plot_dat %>% filter(enrich_group=="2") %>% dplyr::slice(1:10)
sig_3 <- plot_dat %>% filter(enrich_group=="3") %>% dplyr::slice(1:10)

#now bind again 
plot.df <- bind_rows(sig_1, sig_2, sig_3)
plot.df <- plot.df %>% mutate(ord= order(enrich_group), ef_lda)

#barplot 
p1 <- ggplot(plot.df, aes(x=ef_lda, y=fct_reorder(names.colored,ord, .desc = T), fill=enrich_group), alpha=0.8)+
  geom_bar(stat= "identity", position = position_dodge(0.9))+
  scale_fill_manual(name = "Enriched Groups",
                    values=c("1"="#F8766D","2"="#00BA38", "3"="#619CFF"),
                    labels=c("Cluster 1","Cluster 2", "Cluster 3"))+
  guides(fill= guide_legend(reverse=TRUE))+
  theme_pubr()+
  theme(
    plot.title = element_text(face="bold"),
    axis.text.y = element_markdown(angle=0,hjust=1, face = "bold.italic", size = 16),
    axis.text.x = element_text(face = "bold", size = 16), 
    axis.title.x = element_text(face = "bold", size = 20))+
  xlab("LDA (log10)")+
  xlim(0,ceiling(max(merged$ef_lda)))+
  ylab("")

#save it individually 
pdf(file = "metatranscriptomics/Figures/lefse_barplot_k_means_3_clusters_top_10_for_paper.pdf", height = 18, width = 12)
p1
dev.off()

#abundance figure 
#get sample data 
sample.data <- data.frame(sample_data(ps.rel)) %>% 
  dplyr::select(k_means_3_cluster) %>% 
  mutate(k_means_3_cluster = factor(k_means_3_cluster))

#get abundance data 
abund <- data.frame(otu_table(ps.rel))
abund.select <- subset(abund, row.names(abund)%in%plot.df$feature)

#transpose 
abund.select.transposed <- data.frame(t(abund.select))

#now merge the metadata with the abundance table 
abund.select.transposed$k_means_3_cluster <- sample.data$k_means_3_cluster

#organize test1 
test1 <- data.frame(t(abund.select.transposed)) %>% 
  filter(!row_number() %in% 1)
names(test1) <- test1[nrow(test1),]
test1 <- test1[-nrow(test1),]
test1 <- data.frame(sapply(test1[,1:ncol(test1)],as.numeric))
row.names(test1) <- colnames(abund.select.transposed[,2:(ncol(abund.select.transposed)-1)])
#create asv column in test1 and merge with your resutls
test1 <- test1 %>% 
  mutate(tax_name = rownames(test1)) 

#fix taxa names in test1 so they match plot.df$feature
tax_name <- plot.df$feature
tax_name
tax_name <- gsub(" ", ".", tax_name)
tax_name <- gsub("-", ".", tax_name)
plot.df$tax_name <- tax_name

tax_name
rownames(test1)

test1 <- test1%>% 
  inner_join(., plot.df, by="tax_name") %>% 
  mutate(enrich_group = factor(enrich_group)) %>% 
  arrange(enrich_group, desc(ef_lda)) %>% 
  dplyr::select(-enrich_group, -ef_lda, -pvalue, -padj, -marker, -FDR, -number, -taxa, -tax_name, -contaminant,
                -contam_color, -ord,taxa, -names.colored, -taxa)
#transpose DF 
test1 <- data.frame(t(test1))
#set names accoridng to taxa names 
names(test1) <- test1[nrow(test1),]
#remove names row 
test1 <- test1[-nrow(test1),]
#asign categories 
test1$category <- rownames(test1)

# Extract numbers before the dot (if any) and convert to integer
test1$category <- as.integer(sub("X([0-9]+).*", "\\1", test1$category))
test1 <- test1 %>% 
  mutate(category=case_when(
    category== "1"~"1", 
    category== "2"~"2", 
    category== "3"~"3"
  ))

test1$category <- factor(test1$category, levels=c("1","2", "3"))

##### mean abundance plot dataframe
library(dplyr)
a_genus <- pivot_longer(test1, cols = -category , names_to = "id", values_to = "value") %>% 
  as_tibble(.) %>% 
  group_by(id, category) %>% 
  mutate(value = as.numeric(value)) %>% 
  mutate(id = factor(id, levels = unique(id))) %>% 
  dplyr::mutate(count = n()) %>% 
  dplyr::mutate(mean = mean(value)) %>% 
  dplyr::mutate(sd = sd(value)) %>%
  dplyr::mutate(se = sd/sqrt(n())) %>%
  dplyr::mutate(ic=se * qt((1-0.05)/2 + .5, count-1 )) %>%
  mutate(meanpercent = mean * 100) %>%
  ungroup() %>%
  #distinct(.,id,category,count,mean,sd,se,ic, .keep_all= TRUE) %>% 
  arrange(category) %>% 
  mutate(ord = order(category))
a_genus



#relative abundance plot 
p2 <- a_genus %>%
  ggplot(., aes(x = meanpercent, y = fct_reorder(id, -ord), fill=category)) + 
  geom_bar(stat="identity", position=position_dodge(0.9), width = 0.9)+
  geom_bar(stat="identity", position=position_dodge(0.9), width = 0.9, show.legend  =FALSE, size = 0.2)+
  geom_errorbar( aes(xmin=(mean-se)*100, xmax=(mean+se)*100, group = category), width=0.4, colour="black", alpha=0.9, size=.25, position = position_dodge(0.9))+
  scale_fill_manual(name = "Enriched Groups",
                    values=c("1"="#F8766D","2"="#00BA38", "3"="#619CFF"),
                    labels=c("Cluster 1","Cluster 2", "Cluster 3"))+
  guides(fill= guide_legend(reverse=TRUE))+
  theme_classic()+
  theme(
    plot.title = element_text(face="bold"),
    axis.text.y = element_markdown(angle=0,hjust=1, face = "bold.italic", size = 16),
    axis.text.x = element_text(face = "bold", size = 16), 
    axis.title.x = element_text(face = "bold", size = 20))+
  xlab("Relative Abundance (%)")+
  ylab("")+
  ggtitle("")
#coord_cartesian(xlim = c(1, 100))  # Set limits from 1 to 100

p2

#combine two plots 
p3<- ggarrange(p1 + geom_vline(xintercept=lda.thres,linetype="longdash",color="black",size =0.5),
               p2 + 
                 theme(axis.text.y=element_blank(),
                       axis.title=element_text(face="bold"),
                       axis.text.x=element_text(colour="black", face="bold"),
                       axis.ticks=element_line(colour="black"),
                       legend.background = element_rect(color=NA),
                       legend.position = "top"),
               ncol=2,
               nrow=1,
               align = "h",
               labels =c("A","B"),
               common.legend=TRUE,
               legend = "bottom",
               widths = c(1,0.5)
)

p3

pdf(file = "metatranscriptomics/Figures/lefse_full_k_means_3_clusters_top_10_for_paper.pdf", height = 12, width = 10)
p3
dev.off()








##### add ddPCR / bacterial load 

ddpcr_dat <- data.frame(ps@sam_data)

ddpcr_dat <- ddpcr_dat %>% dplyr::mutate(ddPCR=as.numeric(ddPCR))
#calculate stats for ddpcr
my.comparisons <- compare_means(ddPCR~k_means_3_cluster, data = ddpcr_dat)
#save it
#write.csv(my.comparisons, file = "Results/metagenome/ddPCR.stats.Sample_Type.csv")

ddpcr_dat$k_means_3_cluster <- factor(ddpcr_dat$k_means_3_cluster, levels = c("1", "2", "3"))

# plot it
pdf(file = "metatranscriptomics/Figures/ddPCR_k_means_3_clusters.pdf", width=6, height=9)
ggplot(ddpcr_dat, aes(x=k_means_3_cluster, y=ddPCR, fill=k_means_3_cluster))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF"))+
  geom_jitter(color="black", size=1, alpha=0.7)+
  xlab("")+ylab("16S rRNA copies/uL")+
  scale_x_discrete(labels = c("Cluster 1", "Cluster 2", "Cluster 3"))+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(comparisons = list(c("1", "2"), c("2", "3"), c("1", "3")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()


#get thee clusters to compare RNA host 
temp <- data.frame(sample_data(BAL.physeq.pruned))

write.csv(temp, file = "metatranscriptomics/k_means_3_clusters_MT.csv")

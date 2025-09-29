# 0. README -----------
# raw files folder: ~/rawfiles
# PKC file = Mm_R_NGS_WTA_v1.0.pkc


# 1. Set working directory, load libraries, create functions ----------
setwd("setwd")

library(tidyverse)
library(NanoStringNCTools)
library(GeomxTools)
library(SpatialDecon)
library(GeoDiff)
library(GeoMxWorkflows)
library(knitr)
library(dplyr)
library(ggforce)
library(ggplot2)
library(scales)
library(reshape2)
library(cowplot)
library(patchwork)
library(umap)
library(Rtsne)
library(pheatmap)
library(ggrepel)
library(viridis)
library(ggsci)
library(gginnards)

'%!in%' <- function(x,y)!('%in%'(x,y))

# 2. Load data, build object----------

# Reference the main folder 'file.path' containing the sub-folders
datadir <- file.path("~/rawfiles")

# automatically list files in each directory for use
DCCFiles <- dir(file.path(datadir, "DCC-20230503"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)

PKCFiles <- dir(file.path(datadir, "pkcs"), pattern = ".pkc$",
                                full.names = TRUE, recursive = TRUE)
SampleAnnotationFile <-
  dir(file.path(datadir, "annotation"), pattern = ".xlsx$",
      full.names = TRUE, recursive = TRUE)

annotation_file <- readxl::read_excel(SampleAnnotationFile)

demoData <-
  readNanoStringGeoMxSet(dccFiles = DCCFiles,
                         pkcFiles = PKCFiles,
                         phenoDataFile = SampleAnnotationFile,
                         phenoDataSheet = "Sheet1",
                         phenoDataDccColName = "Sample_ID",
                         protocolDataColNames = c("aoi", "roi", "region","exposure","tissueID", "segment_BK", "aoi_BK", "region_exposure_segment"),
                         experimentDataColNames = c("panel"))

# 3. Study design------
## 3.1 modules used-------
pkcs <- annotation(demoData)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))

## 3.2 sample overview -------
count_mat <- dplyr::count(sData(demoData), `slide name`, exposure, region, segment_BK)
test_gr <- gather_set_data(count_mat, 1:4)
test_gr$x <- factor(test_gr$x, levels = c(1:4))
                    # levels = c("exposure", "slide name", "region", "segment_BK"))
# plot Sankey
ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = region), alpha = 0.5, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.2) +
  geom_parallel_sets_labels(color = "white", size = 4) +
  theme_classic(base_size = 17) + 
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank()) +
  scale_y_continuous(expand = expansion(0)) + 
  scale_x_discrete(expand = expansion(0.1), labels = c("slide","exposure","region","segment")) +
  labs(x = "", y = "")
ggsave("sankey_diagram.pdf", plot = last_plot(), width = 5, height = 5, units = "in")


# 4 QC & Pre-processing-------
# Shift counts to one to enable downstream transformations
demoData <- shiftCountsOne(demoData, useDALogic = TRUE)

## 4.1 Segment QC-----------
### 4.1.1 Select Segment QC -------
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10)
       maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 20,         # Minimum # of nuclei estimated (100)
       minArea = 1000)         # Minimum segment area (5000)
demoData <-
  setSegmentQCFlags(demoData, 
                    qcCutoffs = QC_params)  

# Collate QC Results
QCResults <- protocolData(demoData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

QC_Summary

### 4.1.2 Visualize Segment QC------
col_by <- "segment"
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

pdf("QC_histograms_colby_segment.pdf", height = 5, width = 6)
    QC_histogram(sData(demoData), "Trimmed (%)", col_by, 80)
    QC_histogram(sData(demoData), "Stitched (%)", col_by, 80)
    QC_histogram(sData(demoData), "Aligned (%)", col_by, 75)
    
    # % Sequencing saturation ([1-deduplicated reads/aligned reads]%):
    QC_histogram(sData(demoData), "Saturated (%)", col_by, 50) +
      labs(title = "Sequencing Saturation (%)",
           x = "Sequencing Saturation (%)")
    
    # Area:
    QC_histogram(sData(demoData), "area", col_by, 1000, scale_trans = "log10")
    
    
    # Negative Count:
    negativeGeoMeans <- 
      esBy(negativeControlSubset(demoData), 
           GROUP = "Module", 
           FUN = function(x) { 
             assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
           }) 
    protocolData(demoData)[["NegGeoMean"]] <- negativeGeoMeans
    
    # explicitly copy the Negative geoMeans from sData to pData
    negCols <- paste0("NegGeoMean_", modules)
    pData(demoData)[, negCols] <- sData(demoData)[["NegGeoMean"]]
    pData(demoData) %>% names()
    for(ann in negCols) {
      plt <- QC_histogram(pData(demoData), ann, col_by, 2, scale_trans = "log10")
      print(plt)
    }
    dev.off()
    
    # doing histograms for colby = "segment_BK"
    col_by <- "segment_BK"
    pdf("QC_histograms_colby_segment_BK.pdf", height = 5, width = 6)
    QC_histogram(sData(demoData), "Trimmed (%)", col_by, 80)
    QC_histogram(sData(demoData), "Stitched (%)", col_by, 80)
    QC_histogram(sData(demoData), "Aligned (%)", col_by, 75)
    QC_histogram(sData(demoData), "Saturated (%)", col_by, 50) +
      labs(title = "Sequencing Saturation (%)",
           x = "Sequencing Saturation (%)")
    QC_histogram(sData(demoData), "area", col_by, 1000, scale_trans = "log10")
    # QC_histogram(sData(demoData), "nuclei", col_by, 20)
    for(ann in negCols) {
      plt <- QC_histogram(sData(demoData), ann, col_by, 2, scale_trans = "log10")
      print(plt)
    }
    pData(demoData) <- pData(demoData)[, !colnames(pData(demoData)) %in% negCols]
dev.off()

col_by <- "segment"

# No Template Control (NTC) count:
pData(demoData) <- pData(demoData)[, !colnames(pData(demoData)) %in% negCols]

kable(table(NTC_Count = sData(demoData)$NTC),
      col.names = c("NTC Count", "# of Segments"))

# plot all of the QC Summary information in a table
kable(QC_Summary, caption = "QC Summary Table for each Segment")


### 4.1.3 Remove flagged segments-------
demoData <- demoData[, QCResults$QCStatus == "PASS"]

## 4.2 Probe QC-------
### 4.2.1 Set Probe QC Flags -----
demoData <- setBioProbeQCFlags(demoData, 
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = TRUE)

ProbeQCResults <- fData(demoData)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))

### 4.2.2 Exclude Outlier Probes-------
ProbeQCPassed <- 
  subset(demoData, 
         fData(demoData)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(demoData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)

demoData <- ProbeQCPassed 

## 4.3 Create Gene-level Count Data------
length(unique(featureData(demoData)[["TargetName"]]))
target_demoData <- aggregateCounts(demoData)

## 4.4 Limit of Quantification --------
cutoff <- 2
minLOQ <- 2

LOQ <- data.frame(row.names = colnames(target_demoData))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_demoData)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_demoData)[, vars[1]] * 
             pData(target_demoData)[, vars[2]] ^ cutoff)
  }
}
colnames(LOQ) <- "LOQ"
pData(target_demoData)$LOQ <- LOQ

## 4.5 Filtering------
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_demoData)$Module == module
  Mat_i <- t(esApply(target_demoData[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, "LOQ"]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
LOQ_Mat <- LOQ_Mat[fData(target_demoData)$TargetName, ]

### 4.5.1 Segment Gene Detection------
pData(target_demoData)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_demoData)$GeneDetectionRate <-
  pData(target_demoData)$GenesDetected / nrow(target_demoData)

pData(target_demoData)$DetectionThreshold <- 
  cut(pData(target_demoData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 1), 
      labels = c("<1%", "1-2%", "2-3%", "3-4%", "4-5%", "5-10%", "10-15%", ">15%"))

target_demoData <-
  target_demoData[, pData(target_demoData)$GeneDetectionRate >= 0.03] 

### 4.5.2 Gene Detection Rate--------
# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_demoData)]
fData(target_demoData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_demoData)$DetectionRate <-
  fData(target_demoData)$DetectedSegments / nrow(pData(target_demoData))

# Gene of interest detection table
goi <- c("Pdcd1", "Cd274", "Ifng", "Cd8a", "Cd68", "Epcam", "Krt18", "Nphs1", "Nphs2", "Calb1", "Cldn8")

goi_df <- data.frame(
  Gene = goi,
  Number = fData(target_demoData)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(target_demoData)[goi, "DetectionRate"]))

### 4.5.3 Gene Filtering-------
# Subset to target genes detected in at least 10% of the samples.
negativeProbefData <- subset(fData(target_demoData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_demoData <- 
  target_demoData[fData(target_demoData)$DetectionRate >= 0.1 |
                    fData(target_demoData)$TargetName %in% neg_probes, ]

# retain only detected genes of interest
goi <- goi[goi %in% rownames(target_demoData)]



# 5 Normalization--------
library(reshape2)  # for melt
library(cowplot)   # for plot_grid

# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
target_demoData <- normalize(target_demoData ,
                             norm_method = "quant", 
                             desiredQuantile = .75,
                             toElt = "q_norm")

# Background normalization for WTA/CTA without custom spike-in
target_demoData <- normalize(target_demoData ,
                             norm_method = "neg",
                             fromElt = "exprs",
                             toElt = "neg_norm")



# 7 Differential Expression ------
# convert test variables to factors
pData(target_demoData)$testRegion <- factor(sData(target_demoData)$region, c("outside_tumor", "solid_tumor"))
pData(target_demoData)[["slide"]] <- factor(sData(target_demoData)[["slide name"]])
assayDataElement(object = target_demoData, elt = "log_q") <- assayDataApply(target_demoData, 2, FUN = log, base = 2, elt = "q_norm")


# using standR package ------
exprs(target_demoData)  # raw counts
assayDataElement(target_demoData, elt = "q_norm") # q3 norm
assayDataElement(target_demoData, elt = "log_q") # log2 of q3 norm

getwd()
library(standR)
library(SpatialExperiment)
library(ggin)

## convert NanoStringGeoMxSet to SpatialExperiment, with normalized data assay ------
demospe <- as.SpatialExperiment(target_demoData,
                                normData = "q_norm",
                                forceRaw = F)

assay(demospe, "counts") <- assayDataElement(target_demoData, elt = "exprs")
assay(demospe, "logcounts") <- assayDataElement(target_demoData, elt = "q_norm")
assay(demospe, "log_q") <- assayDataElement(target_demoData, elt = "log_q")

## Dimensionality reduction----------
### PCA -------
set.seed(100)
demospe <- scater::runPCA(demospe, exprs_values = "log_q")

### MDS -------
standR::plotMDS(demospe, assay = 2, color = slide)

### UMAP ------
set.seed(42)
demospe <- scater::runUMAP(demospe, dimred = "PCA", exprs_values = "log_q")

### TSNE ------
set.seed(42)
demospe <- scater::runTSNE(demospe, dimred = "PCA", exprs_values = "log_q")

## batch correction using RUV4 (Remove Unwanted Variation 4)------
demospe <- findNCGs(demospe, n_assay = 3, batch_name = "slide", top_n = 300) # this is done on the log2(q3norm) = log_q counts...

NCGs_log_q <- metadata(demospe)$NCGs

spe_ruv <- geomxBatchCorrection(demospe, factors = c("exposure", "region","segment"), 
                                NCGs = metadata(demospe)$NCGs, k = 1, method = "RUV4", isLog = T,
                                n_assay = 3)
set.seed(100)
spe_ruv <- scater::runPCA(spe_ruv)#, exprs_values = "log_q")
pca_results_ruv <- reducedDim(spe_ruv, "PCA")
plotPairPCA(spe_ruv, precomputed = pca_results_ruv, color = exposure, title = "RUV4, k = 1", n_dimension = 4)

## DEG=======
colData(spe_ruv)[,seq(ncol(colData(spe_ruv))-1, ncol(colData(spe_ruv)))] |>
  head()

### Establishing a design matrix and contrast -----------------
## 4c. exposure_region: MOCsolid vs PBSsolid --------
library(edgeR)
library(limma)

spe_ruv$exposure_region <- paste0(spe_ruv$exposure,"_",spe_ruv$region)
spe_ruv_solid <- spe_ruv[, colData(spe_ruv)$region == "solid_tumor"]

dge <- SE2DGEList(spe_ruv_solid)

spe_ruv_solid$exposure_region <- factor(spe_ruv_solid$exposure_region, levels = c("PBS_solid_tumor","MOC_solid_tumor"))

design <- model.matrix(~0 + exposure_region + ruv_W1, data = colData(spe_ruv_solid) )
colnames(design)

colnames(design) <- gsub("^exposure_region","",colnames(design))

contr.matrix <- makeContrasts(
  MOC_solid__v__PBS_solid = MOC_solid_tumor - PBS_solid_tumor,
  levels = colnames(design))

contr.matrix

keep <- filterByExpr(dge, design)
table(keep)
rownames(dge)[!keep]
dge_all <- dge[keep, ]

### BCV check ----------
dge_all <- estimateDisp(dge_all, design = design, robust = TRUE)

plotBCV(dge_all, legend.position = "topleft", ylim = c(0, 1.3))

### Differential expression---------
v <- voom(dge_all, design, plot = TRUE)

fit <- lmFit(v)

fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix)

efit <- eBayes(fit_contrast, robust = TRUE)

results_efit<- decideTests(efit, p.value = 0.05)
summary_efit <- summary(results_efit)

summary_efit

### Visualisation -------
# We can obtain the DE results by using the TopTable function.
library(ggrepel)
library(tidyverse)

de_results_MOC_solid__v__PBS_solid <- topTable(efit, coef = 1, sort.by = "P", n = Inf)

de_genes_toptable_MOC_solid__v__PBS_solid <- topTable(efit, coef = 1, sort.by = "P", n = Inf, p.value = 0.05)

# Categorize Results based on P-value & FDR for plotting
de_results_MOC_solid__v__PBS_solid$Color <- "NS or log2FC < 0.5"
de_results_MOC_solid__v__PBS_solid$Color[de_results_MOC_solid__v__PBS_solid$adj.P.Val < 0.05] <- "FDR < 0.05"
de_results_MOC_solid__v__PBS_solid$Color[de_results_MOC_solid__v__PBS_solid$adj.P.Val < 0.001] <- "FDR < 0.001"
de_results_MOC_solid__v__PBS_solid$Color[abs(de_results_MOC_solid__v__PBS_solid$logFC) < 0.5] <- "NS or log2FC < 0.5"
de_results_MOC_solid__v__PBS_solid$Color <- factor(de_results_MOC_solid__v__PBS_solid$Color,
                                                   levels = c("NS or log2FC < 0.5", "FDR < 0.05", "FDR < 0.001"))

de_results_MOC_solid__v__PBS_solid$invert_P <- (-log10(de_results_MOC_solid__v__PBS_solid$adj.P.Val)) * sign(de_results_MOC_solid__v__PBS_solid$logFC)
top_g <- c()

de_results_MOC_solid__v__PBS_solid$Subset <- 
  ifelse(de_results_MOC_solid__v__PBS_solid$logFC > 0, "MOC_solid_tumor",
         ifelse(de_results_MOC_solid__v__PBS_solid$logFC < 0, "PBS_solid_tumor","NA"))

for(cond in c("MOC_solid_tumor","PBS_solid_tumor")) {
  ind <- de_results_MOC_solid__v__PBS_solid$Subset == cond
  top_g <- c(top_g,
             de_results_MOC_solid__v__PBS_solid[ind, 'TargetName'][
               order(de_results_MOC_solid__v__PBS_solid[ind, 'invert_P'], decreasing = TRUE)[1:15]],
             de_results_MOC_solid__v__PBS_solid[ind, 'TargetName'][
               order(de_results_MOC_solid__v__PBS_solid[ind, 'invert_P'], decreasing = FALSE)[1:15]])
}

top_g <- unique(top_g)
top_g <- top_g[top_g %!in% c("mCherry/tdTomato")]
top_g

de_results_MOC_solid__v__PBS_solid <- de_results_MOC_solid__v__PBS_solid %>% dplyr::select(-invert_P)

ggplot(de_results_MOC_solid__v__PBS_solid,
       aes(x = logFC, y = -log10(adj.P.Val),
           color = Color, label = TargetName)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "Enriched in PBS solid_tumor <- log2(FC) -> Enriched in MOC solid_tumor",
       y = "Significance, -log10(p_adj)",
       color = "Significance") +
  ggtitle("MOC outside_tumor vs MOC solid_tumor")+
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS or log2FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = de_results_MOC_solid__v__PBS_solid %>% dplyr::filter(TargetName %in% top_g) %>% head(20),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

## 4d. exposure_region: MOCoutside vs PBSoutside --------
# To incorporate the limma-voom pipeline, we recommend using the DGElist infrastructure. Our SpatialExperiment can be easily transformed into a DGElist object by using the SE2DGEList function from the edgeR package. For more information about DGEList see ?DGEList.

library(edgeR)
library(limma)

dge <- SE2DGEList(spe_ruv)
spe_ruv$exposure_region <- paste0(spe_ruv$exposure,"_",spe_ruv$region)
spe_ruv_outside <- spe_ruv[, colData(spe_ruv)$region == "outside_tumor"]

dge <- SE2DGEList(spe_ruv_outside)

spe_ruv_outside$exposure_region <- factor(spe_ruv_outside$exposure_region, levels = c("PBS_outside_tumor","MOC_outside_tumor"))

design <- model.matrix(~0 + exposure_region + ruv_W1, data = colData(spe_ruv_outside) )

colnames(design) <- gsub("^exposure_region","",colnames(design))
contr.matrix <- makeContrasts(
  MOC_outside__v__PBS_outside = MOC_outside_tumor - PBS_outside_tumor,
  levels = colnames(design))

keep <- filterByExpr(dge, design)
table(keep)
rownames(dge)[!keep]
dge_all <- dge[keep, ]

### BCV check ----------
dge_all <- estimateDisp(dge_all, design = design, robust = TRUE)

plotBCV(dge_all, legend.position = "topleft", ylim = c(0, 1.3))

### Differential expression---------
v <- voom(dge_all, design, plot = TRUE)

fit <- lmFit(v)

fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix)

efit <- eBayes(fit_contrast, robust = TRUE)

results_efit<- decideTests(efit, p.value = 0.05)
summary_efit <- summary(results_efit)

summary_efit

### Visualisation -------
# We can obtain the DE results by using the TopTable function.
library(ggrepel)
library(tidyverse)

de_results_MOC_outside__v__PBS_outside <- topTable(efit, coef = 1, sort.by = "P", n = Inf)
de_results_MOC_outside__v__PBS_outside %>% head()

de_genes_toptable_MOC_outside__v__PBS_outside <- topTable(efit, coef = 1, sort.by = "P", n = Inf, p.value = 0.05)

# Categorize Results based on P-value & FDR for plotting
de_results_MOC_outside__v__PBS_outside$Color <- "NS or log2FC < 0.5"
de_results_MOC_outside__v__PBS_outside$Color[de_results_MOC_outside__v__PBS_outside$adj.P.Val < 0.05] <- "FDR < 0.05"
de_results_MOC_outside__v__PBS_outside$Color[de_results_MOC_outside__v__PBS_outside$adj.P.Val < 0.001] <- "FDR < 0.001"
de_results_MOC_outside__v__PBS_outside$Color[abs(de_results_MOC_outside__v__PBS_outside$logFC) < 0.5] <- "NS or log2FC < 0.5"
de_results_MOC_outside__v__PBS_outside$Color <- factor(de_results_MOC_outside__v__PBS_outside$Color,
                                                       levels = c("NS or log2FC < 0.5", "FDR < 0.05", "FDR < 0.001"))

de_results_MOC_outside__v__PBS_outside$invert_P <- (-log10(de_results_MOC_outside__v__PBS_outside$adj.P.Val)) * sign(de_results_MOC_outside__v__PBS_outside$logFC)
top_g <- c()

de_results_MOC_outside__v__PBS_outside$Subset <- 
  ifelse(de_results_MOC_outside__v__PBS_outside$logFC > 0, "MOC_outside_tumor",
         ifelse(de_results_MOC_outside__v__PBS_outside$logFC < 0, "PBS_outside_tumor","NA"))

for(cond in c("MOC_outside_tumor","PBS_outside_tumor")) {
  ind <- de_results_MOC_outside__v__PBS_outside$Subset == cond
  top_g <- c(top_g,
             de_results_MOC_outside__v__PBS_outside[ind, 'TargetName'][
               order(de_results_MOC_outside__v__PBS_outside[ind, 'invert_P'], decreasing = TRUE)[1:15]],
             de_results_MOC_outside__v__PBS_outside[ind, 'TargetName'][
               order(de_results_MOC_outside__v__PBS_outside[ind, 'invert_P'], decreasing = FALSE)[1:15]])
}

top_g <- unique(top_g)
top_g <- top_g[top_g %!in% c("mCherry/tdTomato")]
top_g

de_results_MOC_outside__v__PBS_outside <- de_results_MOC_outside__v__PBS_outside %>% dplyr::select(-invert_P)
de_results_MOC_outside__v__PBS_outside
ggplot(de_results_MOC_outside__v__PBS_outside,
       aes(x = logFC, y = -log10(adj.P.Val),
           color = Color, label = TargetName)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "Enriched in PBS outside_tumor <- log2(FC) -> Enriched in MOC outside_tumor",
       y = "Significance, -log10(p_adj)",
       color = "Significance") +
  ggtitle("MOC outside_tumor vs MOC outside_tumor")+
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS or log2FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = de_results_MOC_outside__v__PBS_outside %>% dplyr::filter(TargetName %in% top_g) %>% head(20),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")


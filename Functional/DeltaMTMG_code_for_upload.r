######## Delta MT-MG analysis of lower airway microbiome 
######## paper: immunotherapy lung cancer paper 
######## Fares Darawshy MD
######## Sep 2025 



##/////////////////// note: the code includes data and figures for comparison according to clinical outcomes, which are not shown in the paper ########

#Load Packages
library(DESeq2)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(pathfindR)
library(scales)
library(data.table)
library(fBasics)
library(forcats)
library(omu)
library(dplyr)
library(maptools)
library(phyloseq)
library(SpiecEasi)
library(vegan)
library(tibble)
library(themetagenomics)
library(grid)
library(dplyr)
library(tibble)
library(formattable)
library("htmltools")
library("webshot")    
library(splitstackshape)
library(decontam)
library(cowplot)
library(wesanderson)
library(colorspace)
library(factoextra)
library(Rtsne)
library(cluster)
library(umap)
library(magrittr)
library(tidyr)

#Set Theme for Figures
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    strip.background=element_blank(),axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),
    plot.margin=unit(c(1,1,1,1),"line"),
    legend.key=element_rect(fill='white'))

theme2<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    strip.background=element_blank(),axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),
    plot.margin=unit(c(1,1,1,1),"line"),legend.position = "none",
    legend.key=element_rect(fill='white'))

theme3<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    strip.background=element_blank(),axis.text.x=element_text(colour=a),
    axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),
    plot.margin=unit(c(1,1,1,1),"line"),legend.position = "none",
    legend.key=element_rect(fill='white'))

theme4<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
    panel.grid.major = element_line(color="grey"),panel.grid.minor = element_blank(),
    strip.background=element_blank(),axis.text.x=element_text(colour="black", size = 16, face = "bold"),
    axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),
    strip.text.y = element_blank(),
    plot.margin=unit(c(1,1,1,1),"line"),
    legend.key=element_rect(fill='white'))


theme5<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
    panel.grid.major = element_line(color="grey"),panel.grid.minor = element_blank(),
    strip.background=element_blank(),axis.text.x=element_text(colour="black", size = 16, face = "bold"),
    strip.text.y = element_blank(),
    axis.text.y=element_blank(),axis.ticks=element_line(colour="black"),
    plot.margin=unit(c(1,1,1,1),"line"),
    legend.key=element_rect(fill='white'))

theme6<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
    panel.grid.major = element_line(color="grey"),panel.grid.minor = element_blank(),
    strip.background=element_rect(color="black", fill="#FC4E07", size=1.5, linetype="solid"),
    strip.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x=element_text(colour="black", size = 16, face = "bold"),
    axis.text.y=element_blank(),axis.ticks=element_line(colour="black"),
    plot.margin=unit(c(1,1,1,1),"line"),
    legend.key=element_rect(fill='white'))





############## PD as an outcome ##################

#//////////////////////////////////
##read and prepare METAGENOME DataSet
#//////////////////////////////////
# load metadata 
coldata <- read.csv(file = "MG_metabolic_metadata.csv")
rownames(coldata) <- coldata$X #Metagenome
coldata <- coldata[,-1]


# load metagenome dataset 
mycounts <-read.delim2("MG_metabolic_functional_counts.txt", sep="\t", row.names=1)
#mycounts <-read.delim2("counts.metagenome.ko.module.fmap.txt", sep="\t", row.names=1)
colnames(mycounts) <- gsub("X","",colnames(mycounts))

#Create Relative Abundance Table
rel <-as.data.frame(mycounts) #make counts dataframe
relcounts <-
  rel %>% 
  tibble::rownames_to_column('gs') %>%
  group_by(gs) %>% 
  summarise_all(funs(sum)) %>%
  mutate_if(is.numeric, funs(./sum(.))) %>%
  tibble::column_to_rownames('gs')

#Convert any NAs to 0
mycounts[is.na(mycounts)] <- 0
relcounts[is.na(relcounts)] <- 0
#Copy of Count Table
mycounts2 <- mycounts
relcounts2 <- relcounts
#Convert Count Table into a Numeic Data Frame
d1 = data.frame(lapply(mycounts2, function(x) as.numeric(as.character(x))),
                check.names=F, row.names = rownames(mycounts2))
d2 = data.frame(lapply(relcounts2, function(x) as.numeric(as.character(x))),
                check.names=F, row.names = rownames(relcounts2))
#Remove 0 counts in relative Table
d2 <- d2[rowSums(d2[, -1] > 0) != 0, ]
d2 <- d2 %>% select(which(!colSums(d2, na.rm=TRUE) %in% 0))

#get the columns to match between count table and relative table 
wanted<-which(colnames(d1) %in% colnames(d2))
wanted2<-which(rownames(coldata) %in% colnames(d2))
d1<-d1[,wanted]
coldata<-coldata[wanted2,]

#df <- relcounts3
#Subset group 1 
metag_group_1_coldata <- coldata[coldata$ultimate_response_2_lev=="PD",]
# now subset thier relative data and counts data 
metg_group_1_counts <- d1[,rownames(metag_group_1_coldata)]
metg_group_1_rel_counts <- d2[,rownames(metag_group_1_coldata)]

#subset group 2 
metag_group_2_coldata <- coldata[coldata$ultimate_response_2_lev=="non_PD",]
# now subset thier relative data and counts data 
metg_group_2_counts <- d1[,rownames(metag_group_2_coldata)]
metg_group_2_rel_counts <- d2[,rownames(metag_group_2_coldata)]

# all metagenome metadata 
metag_metadata <- coldata
# all metagenome counts 
metag_counts <- d1
# all metagenome relative counts 
metag_rel_counts <- d2
  

############# metatranscriptome dataset
# load metadata 
coldata <- read.csv(file = "MT_metabolic_metadata.csv")
rownames(coldata) <- coldata$X #Metagenome
coldata <- coldata[,-1]


# load  dataset 
mycounts <-read.delim2("MT_metabolic_functional_counts.txt", sep="\t", row.names=1)
#mycounts <-read.delim2("counts.metagenome.ko.module.fmap.txt", sep="\t", row.names=1)
colnames(mycounts) <- gsub("X","",colnames(mycounts))

#Create Relative Abundance Table
rel <-as.data.frame(mycounts) #make counts dataframe
relcounts <-
  rel %>% 
  tibble::rownames_to_column('gs') %>%
  group_by(gs) %>% 
  summarise_all(funs(sum)) %>%
  mutate_if(is.numeric, funs(./sum(.))) %>%
  tibble::column_to_rownames('gs')

#Convert any NAs to 0
mycounts[is.na(mycounts)] <- 0
relcounts[is.na(relcounts)] <- 0
#Copy of Count Table
mycounts2 <- mycounts
relcounts2 <- relcounts
#Convert Count Table into a Numeic Data Frame
d1 = data.frame(lapply(mycounts2, function(x) as.numeric(as.character(x))),
                check.names=F, row.names = rownames(mycounts2))
d2 = data.frame(lapply(relcounts2, function(x) as.numeric(as.character(x))),
                check.names=F, row.names = rownames(relcounts2))
#Remove 0 counts in relative Table
d2 <- d2[rowSums(d2[, -1] > 0) != 0, ]
d2 <- d2 %>% select(which(!colSums(d2, na.rm=TRUE) %in% 0))

#get the columns to match between count table and relative table 
wanted<-which(colnames(d1) %in% colnames(d2))
wanted2<-which(rownames(coldata) %in% colnames(d2))
d1<-d1[,wanted]
coldata<-coldata[wanted2,]

#Subset group 1 
metat_group_1_coldata <- coldata[coldata$ultimate_response_2_lev=="PD",]
# now subset thier relative data and counts data 
mett_group_1_counts <- d1[,rownames(metat_group_1_coldata)]
mett_group_1_rel_counts <- d2[,rownames(metat_group_1_coldata)]

#subset group 2 
metat_group_2_coldata <- coldata[coldata$ultimate_response_2_lev=="non_PD",]
# now subset thier relative data and counts data 
mett_group_2_counts <- d1[,rownames(metat_group_2_coldata)]
mett_group_2_rel_counts <- d2[,rownames(metat_group_2_coldata)]

# all metatranscriptome metadata 
metat_metadata <- coldata
# all metatrans riptome counts 
metat_counts <- d1
# all metatranscripotome relative counts 
metat_rel_counts <- d2




########### delta table ####### 
#match two tables of counts 
setdiff(rownames(metag_counts),rownames(metat_counts))

#set the same order 
#Order Tables
metag_counts <- metag_counts[,order(colnames(metag_counts))]
metag_counts <- metag_counts[order(rownames(metag_counts)),]
metat_counts <- metat_counts[,order(colnames(metat_counts))]
metat_counts <- metat_counts[order(rownames(metat_counts)),]

metag_rel_counts <- metag_rel_counts[,order(colnames(metag_rel_counts))]
metag_rel_counts <- metag_rel_counts[order(rownames(metag_rel_counts)),]
metat_rel_counts <- metat_rel_counts[,order(colnames(metat_rel_counts))]
metat_rel_counts <- metat_rel_counts[order(rownames(metat_rel_counts)),]

#Match Tables
rows1<-which(rownames(metag_counts)  %in%  rownames(metat_counts))
metag_counts_2<-metag_counts[rows1,]

rows2<-which(rownames(metat_counts)  %in%  rownames(metag_counts_2))
metat_counts_2<-metat_counts[rows2,]

cols1 <- which(colnames(metag_counts_2) %in% colnames(metat_counts_2))
metag_counts_3<-metag_counts_2[,cols1]
metat_counts_3<-metat_counts_2

#match relative tables 
rows1<-which(rownames(metag_rel_counts)  %in%  rownames(metat_rel_counts))
metag_rel_counts_2<-metag_rel_counts[rows1,]

rows2<-which(rownames(metat_rel_counts)  %in%  rownames(metag_rel_counts_2))
metat_rel_counts_2<-metat_rel_counts[rows2,]

cols1 <- which(colnames(metag_rel_counts_2) %in% colnames(metat_rel_counts_2))
metag_rel_counts_3<-metag_rel_counts_2[,cols1]
metat_rel_counts_3<-metat_rel_counts_2

#add one
metag_counts_4 <- metag_counts_3+1
metat_counts_4 <- metat_counts_3+1

metag_rel_counts_4 <- metag_rel_counts_3+1
metat_rel_counts_4 <- metat_rel_counts_3+1

#multiply by 100
metag_rel_counts_5 <- metag_rel_counts_3*100
metat_rel_counts_5 <- metat_rel_counts_3*100

#Create Ratio Table
ratio <- cbind(metat_counts_4[1],round(metat_counts_4[-1]/metag_counts_4[-1],1))

#Create Delta Table
delta <- metat_counts_4-metag_counts_4
#deltarel <- metatcountrel4-metagcountrel4
deltarel <- metat_rel_counts_5-metag_rel_counts_5

#Transpose
deltarel_trans <- as.data.frame(t(deltarel))
deltarel_trans$ID <- rownames(deltarel_trans)

#get Metadata of interest
deltarel_trans_meta <- coldata %>% select(ID,ultimate_response_2_lev)
setdiff(deltarel_trans_meta$ID,deltarel_trans$ID)
deltarel_trans_meta$ID <- gsub("_", ".", deltarel_trans_meta$ID)

#merge MetaData
deltarel_trans <- merge(deltarel_trans,deltarel_trans_meta,by="ID",all.x=TRUE)
deltarel_trans$ultimate_response_2_lev
# remove the ID number
df<-deltarel_trans%>%select(-ID)
deltarel_trans$ultimate_response_2_lev
#LOAD FILE
#load(file="COPD.TAXA.Compare.RData")
#BAL and UA
#df <- df[df$Sample.Type!="BKG",]

# get the means of all the variables
means<-df%>%
  group_by(ultimate_response_2_lev)%>%
  summarise_all(list(mean), na.rm=TRUE )%>%
  gather("Variable", "Mean", -ultimate_response_2_lev)%>%
  spread(ultimate_response_2_lev, Mean)%>%
  rename("Mean_PD"='PD', "Mean_non_PD"="non_PD")

# get the standard deviation of all the variables
sds<-df%>%
  group_by(ultimate_response_2_lev)%>%
  summarise_all(list(sd), na.rm=TRUE )%>%
  gather("Variable", "SD", -ultimate_response_2_lev)%>%
  spread(ultimate_response_2_lev, SD)%>%
  rename("SD_PD"='PD', "SD_non_PD"="non_PD")

# get the median of all the variables
medians<-df%>%
  group_by(ultimate_response_2_lev)%>%
  summarise_all(list(median), na.rm=TRUE )%>%
  gather("Variable", "Median", -ultimate_response_2_lev)%>%
  spread(ultimate_response_2_lev, Median)%>%
  rename("Median_PD"='PD', "Median_non_PD"="non_PD")


# join the tables 
summary_report<-means%>%
  inner_join(sds, "Variable")%>%
  inner_join(medians, "Variable")%>%
  mutate(DiffInMedians=Median_PD-Median_non_PD)

############ perform Mann Whitney test and export as a table 
# Select variables of interest
df2 <- df %>% select(-ultimate_response_2_lev)
variables <- colnames(df2)

# Initialize vectors to store results
pvals <- numeric(length(variables))
vars <- character(length(variables))

# Loop through variables
for (i in seq_along(variables)) {
  
  # Extract data for the Mann-Whitney U test
  data <- df %>% select(ultimate_response_2_lev, variables[i])
  
  # Extract samples for BAL and UA
  x1 <- data %>%
    filter(ultimate_response_2_lev == "PD") %>%
    pull(variables[i]) %>%
    na.omit()
  
  x2 <- data %>%
    filter(ultimate_response_2_lev == "non_PD") %>%
    pull(variables[i]) %>%
    na.omit()
  
  # Perform the Mann-Whitney U test
  wc <- wilcox.test(x1, x2)
  
  # Store results
  pvals[i] <- round(wc$p.value, 4)
  vars[i] <- variables[i]
}


wc_df<-data.frame(Variable=vars,pvalues=pvals)

wc_df$Variable<-as.character(wc_df$Variable)
#Create Summary Report
summary_report<-summary_report%>%inner_join(wc_df, by="Variable")
#Order Summary Report
summary_report <- summary_report[order(summary_report$DiffInMedians),]
#Subset Significant Only
summary_report_sig <- summary_report[summary_report$pvalues<0.05,]
#Order by Median Diff
summary_report_sig <- summary_report_sig[order(summary_report_sig$DiffInMedians),]
summary_report_sig$order <- 1:nrow(summary_report_sig)



########### data preparing for figures ###########
#If not using contaminants
summary_report_sig_c <- summary_report_sig

#different Ordering by group 1
summary_report_sig_c <- summary_report_sig_c[order(-summary_report_sig_c$Median_PD),]
summary_report_sig_c$order <- 1:nrow(summary_report_sig_c)

#different Ordering by group 2
summary_report_sig_c <- summary_report_sig_c[order(-summary_report_sig_c$Median_non_PD),]
summary_report_sig_c$order2 <- 1:nrow(summary_report_sig_c)

deltaorder <- summary_report_sig_c %>% select(Variable,order,order2)
deltaorder <- as.data.frame(deltaorder)
deltaorder <- deltaorder[!is.na(deltaorder$Variable),]
rownames(deltaorder) <- deltaorder$Variable
deltaorder$taxa <- rownames(deltaorder)

#Set new order
deltaorder$neworder <- 1:nrow(deltaorder)
#Get the rel abundance tables
rows4 <- which(rownames(deltarel) %in% rownames(deltaorder))
deltarel_sig <-deltarel[rows4,]
#get Metadata of interest
deltameta <- coldata %>% select(ID,ultimate_response_2_lev)
deltameta$ID <- gsub("_", ".", deltameta$ID)
#Transpose Data
deltarel_sig$taxa <- rownames(deltarel_sig)
ncol(deltarel_sig)
#select only columns of samples 
deltarel_sig2 <- deltarel_sig %>% tidyr::gather(ID, abundance,1:76)
#merge MetaData
deltarel_sig3 <- merge(deltarel_sig2,deltameta,by="ID",all.x=TRUE)
#merge order
deltaorder2 <- deltaorder %>% select(taxa,neworder)
deltarel_sig3 <- merge(deltarel_sig3,deltaorder2, by="taxa")

#set Order of Figure
deltarel_sig3$ultimate_response_2_lev <- factor(deltarel_sig3$ultimate_response_2_lev, levels = c("non_PD", "PD"))
deltarel_sig3 <- deltarel_sig3 %>% arrange(neworder)

#Set division of Figure by selecting to p3 pathways. in this case we will select top 4 pathways for PD
deltarel_sig3$abundance
#deltarel_sig3$top <- ifelse(deltarel_sig3$neworder>25,2,1)
deltarel_sig3$top <- ifelse(deltarel_sig3$neworder>5,"Top_non_PD","Top_PD")
#set order of top to adapt for the figure 
deltarel_sig3$top <- factor(deltarel_sig3$top, levels = c("Top_PD","Top_non_PD"))

#Get relabundance tables for each group 
library(tidyverse)
#for PD gorup 
rel.otuC1 <-               
  metat_counts %>% dplyr::select(rownames(metat_group_1_coldata)) %>% 
  rownames_to_column('gs') %>%
  group_by(gs) %>% 
  summarise_all(funs(sum)) %>%
  mutate_if(is.numeric, funs(./sum(.))) %>%
  column_to_rownames('gs')
#for non PD group 
rel.otuC2 <-
  metat_counts %>% dplyr::select(rownames(metat_group_2_coldata)) %>% 
  rownames_to_column('gs') %>%
  group_by(gs) %>% 
  summarise_all(funs(sum)) %>%
  mutate_if(is.numeric, funs(./sum(.))) %>%
  column_to_rownames('gs')
rel.otuC1 <- rel.otuC1*100
rel.otuC2 <- rel.otuC2*100

#decide what otu to save 
deltarel_sig3_non_PD <- deltarel_sig3[deltarel_sig3$ultimate_response_2_lev=="non_PD",]
deltarel_sig3_PD <- deltarel_sig3[deltarel_sig3$ultimate_response_2_lev=="PD",]
otu.to.save1 <-as.character(deltarel_sig3_non_PD$taxa)
otu.to.save2 <-as.character(deltarel_sig3_PD$taxa)
#convert to dataframe
df.1.df <- data.frame(rel.otuC1)
df.2.df <- data.frame(rel.otuC2)
#from relative table we should get the mean across the row of the otu table
df.1.meanRA <- rowMeans(df.1.df)
df.2.meanRA <- rowMeans(df.2.df)
#need to subset AND reorder just the otus that we have 
df.1.meanRA.save <- df.1.meanRA[otu.to.save1]
df.2.meanRA.save <- df.2.meanRA[otu.to.save2]
#add the abundnace data for the res dataframe
deltarel_sig3_non_PD$relabundance <- df.1.meanRA.save
deltarel_sig3_PD$relabundance <- df.2.meanRA.save

#Keep abundant ones

#Barplot data
deltarel_sig3_non_PD_bar <- deltarel_sig3_non_PD %>% distinct(taxa, .keep_all=TRUE)
deltarel_sig3_PD_bar <- deltarel_sig3_PD %>% distinct(taxa, .keep_all=TRUE)
barplot  <- rbind(deltarel_sig3_non_PD_bar,deltarel_sig3_PD_bar)


########### figures pltos ############
p1 <- ggplot(barplot, aes(x=reorder(taxa,neworder), y=relabundance, fill=ultimate_response_2_lev)) + 
  #stat_boxplot(geom ='errorbar', width=0.1)+
  geom_bar(stat="identity", aes(fill=ultimate_response_2_lev), width=0.5,position = position_dodge(width = 0.8))+
  #xmin=as.numeric(id)-.4,xmax=as.numeric(id)+.4, x=id, ymin=10E-7, ymax=ymean, fill=var)
  #geom_rect(position=position_dodge(.8), aes(xmin=as.numeric(reorder(taxa,+order))-.4,xmax=as.numeric(reorder(taxa,+order))+.4, ymin=0, ymax=relabundance,fill=Sample.Type))+
  #geom_rect(aes(xmin=as.numeric(reorder(taxa, +order))-0.4,xmax=as.numeric(reorder(taxa, +order))+.4, ymax=relabundance, ymin=0, fill=Sample.Type),color="black",position=position_dodge(0.9)) +
  #geom_boxplot(aes(color=Sample.Type),outlier.shape = NA, width=0.5)+
  #geom_boxplot(aes(ymin=..lower.., ymax=..upper..))+
  #geom_point(aes(color=Sample.Type),position = position_jitter(width = 0.25, seed = 123), alpha=0.7)+
  #geom_line(aes(group = SampleID),position = position_jitter(width = 0.25, seed = 123),linetype = "dashed",color="grey")+
  #geom_jitter(aes(color=Subject_Type_2),shape=1, position=position_jitter(0.2))+
  #geom_line(aes(group = SampleID), position=position_jitter(0.2),linetype = "dashed")+
  #facet_wrap(top ~ . , nrow=2,ncol=1)+
  #facet_grid(top ~., scales = "free_y", space = "free_y") 
facet_grid(top ~ . , scales = "free", space = "free")+
  #scale_color_manual(values=c("Asymptomatic_SC"="#4DAC26","Symptomatic_SC"="#FFD479","COPD"="#D01C8B")) +
  scale_fill_manual(values=c("PD"="red","non_PD"="blue"),guide=FALSE) +
  #scale_x_discrete(labels = c('BKG','BPT','SPT','Sup'))+ 
  #geom_text_repel(aes(label=ifelse(res$sig > 3 , as.character(res$o),'')),size=3,force=25) +
  #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
  #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
  #scale_y_continuous(name="Log10 Rel Ab",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
  #scale_shape_manual(values=c("MG"=17,"MT"=16), guide=guide_legend(title="Sequencing",size=c(3, 3)))+
  #ylab("[MT-MG]") + 
  xlab("")+
  coord_flip()+
  #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
  #geom_point(color=cols) +
  #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
  theme4

p1

p2 <-    ggplot(deltarel_sig3_non_PD, aes(x=reorder(taxa,neworder), y=abundance)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(aes(color=ultimate_response_2_lev),outlier.shape = NA, width=0.5)+
  #geom_boxplot(aes(ymin=..lower.., ymax=..upper..))+
  geom_point(aes(color=ultimate_response_2_lev),position = position_jitter(width = 0.25, seed = 123), alpha=0.7)+
  coord_flip()+
  #geom_line(aes(group = SampleID),position = position_jitter(width = 0.25, seed = 123),linetype = "dashed",color="grey")+
  #geom_jitter(aes(color=Subject_Type_2),shape=1, position=position_jitter(0.2))+
  #geom_line(aes(group = SampleID), position=position_jitter(0.2),linetype = "dashed")+
  #facet_wrap(~ Sample.Type, ncol=3)+
  facet_grid(top ~ ., scales = "free_y", space = "free")+
  #scale_color_manual(values=c("Asymptomatic_SC"="#4DAC26","Symptomatic_SC"="#FFD479","COPD"="#D01C8B")) +
  scale_color_manual(values=c("PD"="red","non_PD"="blue"),guide=FALSE) +
  #scale_x_discrete(labels = c('BKG','BPT','SPT','Sup'))+ 
  #geom_text_repel(aes(label=ifelse(res$sig > 3 , as.character(res$o),'')),size=3,force=25) +
  #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
  #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
  #scale_y_continuous(name="Delta MT-MG",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
  #scale_shape_manual(values=c("MG"=17,"MT"=16), guide=guide_legend(title="Sequencing",size=c(3, 3)))+
  ylab("[MT-MG]") + 
  xlab("")+
  scale_y_continuous(limits = c(-10, 10))+
  #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
  #geom_point(color=cols) +
  #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
  theme6

p3 <-    ggplot(deltarel_sig3_PD, aes(x=reorder(taxa,neworder), y=abundance)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(aes(color=ultimate_response_2_lev),outlier.shape = NA, width=0.5)+
  #geom_boxplot(aes(ymin=..lower.., ymax=..upper..))+
  geom_point(aes(color=ultimate_response_2_lev),position = position_jitter(width = 0.25, seed = 123), alpha=0.7)+
  coord_flip()+
  #geom_line(aes(group = SampleID),position = position_jitter(width = 0.25, seed = 123),linetype = "dashed",color="grey")+
  #geom_jitter(aes(color=Subject_Type_2),shape=1, position=position_jitter(0.2))+
  #geom_line(aes(group = SampleID), position=position_jitter(0.2),linetype = "dashed")+
  #facet_wrap(~ Sample.Type, ncol=3)+
  facet_grid(top ~ ., scales = "free_y", space = "free")+
  #scale_color_manual(values=c("Asymptomatic_SC"="#4DAC26","Symptomatic_SC"="#FFD479","COPD"="#D01C8B")) +
  scale_color_manual(values=c("PD"="red","non_PD"="blue"),guide=FALSE) +
  #scale_x_discrete(labels = c('BKG','BPT','SPT','Sup'))+ 
  #geom_text_repel(aes(label=ifelse(res$sig > 3 , as.character(res$o),'')),size=3,force=25) +
  #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
  #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
  #scale_y_continuous(name="Delta MT-MG",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
  #scale_shape_manual(values=c("MG"=17,"MT"=16), guide=guide_legend(title="Sequencing",size=c(3, 3)))+
  ylab("[MT-MG]") + 
  xlab("")+
  scale_y_continuous(limits = c(-10, 10))+
  #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
  #geom_point(color=cols) +
  #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
  theme5
dummy <- ggplot(data = deltarel_sig3_non_PD, aes(x=reorder(taxa,-neworder), y=abundance))+
  facet_grid(top ~ ., scales = "free_y", space = "free")+
  geom_rect(aes(fill=size), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme_minimal()

library(cowplot)

pdf("PD_vs_no_PD_.DELTA_sig_abundance_Functional_Module.pdf", height = 6, width = 16)
#plot_grid(p1,p3,p2, labels="", rel_widths=c(1.1,1.5,1.5),ncol=3)
plot_grid(p1,p3,p2, labels="", ncol=3)
dev.off()














#########two years mortality outcome #########

#//////////////////////////////////
#------------------->METAGENOME DataSet
#//////////////////////////////////
# load metadata 
coldata <- read.csv(file = "MG_metabolic_metadata.csv")
rownames(coldata) <- coldata$X #Metagenome
coldata <- coldata[,-1]


# load metagenome dataset 
mycounts <-read.delim2("MG_metabolic_functional_counts.txt", sep="\t", row.names=1)
#mycounts <-read.delim2("counts.metagenome.ko.module.fmap.txt", sep="\t", row.names=1)
colnames(mycounts) <- gsub("X","",colnames(mycounts))

#Create Relative Abundance Table
rel <-as.data.frame(mycounts) #make counts dataframe
relcounts <-
  rel %>% 
  tibble::rownames_to_column('gs') %>%
  group_by(gs) %>% 
  summarise_all(funs(sum)) %>%
  mutate_if(is.numeric, funs(./sum(.))) %>%
  tibble::column_to_rownames('gs')

#Convert any NAs to 0
mycounts[is.na(mycounts)] <- 0
relcounts[is.na(relcounts)] <- 0
#Copy of Count Table
mycounts2 <- mycounts
relcounts2 <- relcounts
#Convert Count Table into a Numeic Data Frame
d1 = data.frame(lapply(mycounts2, function(x) as.numeric(as.character(x))),
                check.names=F, row.names = rownames(mycounts2))
d2 = data.frame(lapply(relcounts2, function(x) as.numeric(as.character(x))),
                check.names=F, row.names = rownames(relcounts2))
#Remove 0 counts in relative Table
d2 <- d2[rowSums(d2[, -1] > 0) != 0, ]
d2 <- d2 %>% select(which(!colSums(d2, na.rm=TRUE) %in% 0))

#get the columns to match between count table and relative table 
wanted<-which(colnames(d1) %in% colnames(d2))
wanted2<-which(rownames(coldata) %in% colnames(d2))
d1<-d1[,wanted]
coldata<-coldata[wanted2,]

#df <- relcounts3
#Subset group 1 
metag_group_1_coldata <- coldata %>% filter(two_y_mort=="less_24")
# now subset thier relative data and counts data 
metg_group_1_counts <- d1[,rownames(metag_group_1_coldata)]
metg_group_1_rel_counts <- d2[,rownames(metag_group_1_coldata)]

#subset group 2 
metag_group_2_coldata <- coldata %>% filter(two_y_mort=="greater_24")
# now subset thier relative data and counts data 
metg_group_2_counts <- d1[,rownames(metag_group_2_coldata)]
metg_group_2_rel_counts <- d2[,rownames(metag_group_2_coldata)]

# all metagenome metadata 
metag_metadata <- coldata
# all metagenome counts 
metag_counts <- d1
# all metagenome relative counts 
metag_rel_counts <- d2


############# metatranscriptome dataset
# load metadata 
coldata <- read.csv(file = "MT_metabolic_metadata.csv")
rownames(coldata) <- coldata$X #Metagenome
coldata <- coldata[,-1]


# load  dataset 
mycounts <-read.delim2("MT_metabolic_functional_counts.txt", sep="\t", row.names=1)
#mycounts <-read.delim2("counts.metagenome.ko.module.fmap.txt", sep="\t", row.names=1)
colnames(mycounts) <- gsub("X","",colnames(mycounts))

#Create Relative Abundance Table
rel <-as.data.frame(mycounts) #make counts dataframe
relcounts <-
  rel %>% 
  tibble::rownames_to_column('gs') %>%
  group_by(gs) %>% 
  summarise_all(funs(sum)) %>%
  mutate_if(is.numeric, funs(./sum(.))) %>%
  tibble::column_to_rownames('gs')

#Convert any NAs to 0
mycounts[is.na(mycounts)] <- 0
relcounts[is.na(relcounts)] <- 0
#Copy of Count Table
mycounts2 <- mycounts
relcounts2 <- relcounts
#Convert Count Table into a Numeic Data Frame
d1 = data.frame(lapply(mycounts2, function(x) as.numeric(as.character(x))),
                check.names=F, row.names = rownames(mycounts2))
d2 = data.frame(lapply(relcounts2, function(x) as.numeric(as.character(x))),
                check.names=F, row.names = rownames(relcounts2))
#Remove 0 counts in relative Table
d2 <- d2[rowSums(d2[, -1] > 0) != 0, ]
d2 <- d2 %>% select(which(!colSums(d2, na.rm=TRUE) %in% 0))

#get the columns to match between count table and relative table 
wanted<-which(colnames(d1) %in% colnames(d2))
wanted2<-which(rownames(coldata) %in% colnames(d2))
d1<-d1[,wanted]
coldata<-coldata[wanted2,]

#Subset group 1 
metat_group_1_coldata <- coldata %>% filter(two_y_mort=="less_24")
# now subset thier relative data and counts data 
mett_group_1_counts <- d1[,rownames(metat_group_1_coldata)]
mett_group_1_rel_counts <- d2[,rownames(metat_group_1_coldata)]

#subset group 2 
metat_group_2_coldata <- coldata %>% filter(two_y_mort=="greater_24")
# now subset thier relative data and counts data 
mett_group_2_counts <- d1[,rownames(metat_group_2_coldata)]
mett_group_2_rel_counts <- d2[,rownames(metat_group_2_coldata)]

# all metatranscriptome metadata 
metat_metadata <- coldata
# all metatrans riptome counts 
metat_counts <- d1
# all metatranscripotome relative counts 
metat_rel_counts <- d2




########### delta table ####### 
#match two tables of counts 
setdiff(rownames(metag_counts),rownames(metat_counts))

#set the same order 
#Order Tables
metag_counts <- metag_counts[,order(colnames(metag_counts))]
metag_counts <- metag_counts[order(rownames(metag_counts)),]
metat_counts <- metat_counts[,order(colnames(metat_counts))]
metat_counts <- metat_counts[order(rownames(metat_counts)),]

metag_rel_counts <- metag_rel_counts[,order(colnames(metag_rel_counts))]
metag_rel_counts <- metag_rel_counts[order(rownames(metag_rel_counts)),]
metat_rel_counts <- metat_rel_counts[,order(colnames(metat_rel_counts))]
metat_rel_counts <- metat_rel_counts[order(rownames(metat_rel_counts)),]

#Match Tables
rows1<-which(rownames(metag_counts)  %in%  rownames(metat_counts))
metag_counts_2<-metag_counts[rows1,]

rows2<-which(rownames(metat_counts)  %in%  rownames(metag_counts_2))
metat_counts_2<-metat_counts[rows2,]

cols1 <- which(colnames(metag_counts_2) %in% colnames(metat_counts_2))
metag_counts_3<-metag_counts_2[,cols1]
metat_counts_3<-metat_counts_2

#match relative tables 
rows1<-which(rownames(metag_rel_counts)  %in%  rownames(metat_rel_counts))
metag_rel_counts_2<-metag_rel_counts[rows1,]

rows2<-which(rownames(metat_rel_counts)  %in%  rownames(metag_rel_counts_2))
metat_rel_counts_2<-metat_rel_counts[rows2,]

cols1 <- which(colnames(metag_rel_counts_2) %in% colnames(metat_rel_counts_2))
metag_rel_counts_3<-metag_rel_counts_2[,cols1]
metat_rel_counts_3<-metat_rel_counts_2

#add one
metag_counts_4 <- metag_counts_3+1
metat_counts_4 <- metat_counts_3+1

metag_rel_counts_4 <- metag_rel_counts_3+1
metat_rel_counts_4 <- metat_rel_counts_3+1

#multiply by 100
metag_rel_counts_5 <- metag_rel_counts_3*100
metat_rel_counts_5 <- metat_rel_counts_3*100

#Create Ratio Table
ratio <- cbind(metat_counts_4[1],round(metat_counts_4[-1]/metag_counts_4[-1],1))

#Create Delta Table
delta <- metat_counts_4-metag_counts_4
#deltarel <- metatcountrel4-metagcountrel4
deltarel <- metat_rel_counts_5-metag_rel_counts_5

#Transpose
deltarel_trans <- as.data.frame(t(deltarel))
deltarel_trans$ID <- rownames(deltarel_trans)

#get Metadata of interest
deltarel_trans_meta <- coldata %>% select(ID,two_y_mort)
setdiff(deltarel_trans_meta$ID,deltarel_trans$ID)
deltarel_trans_meta$ID <- gsub("_", ".", deltarel_trans_meta$ID)

#merge MetaData
deltarel_trans <- merge(deltarel_trans,deltarel_trans_meta,by="ID",all.x=TRUE)
deltarel_trans$two_y_mort
# remove the ID number
df<-deltarel_trans%>%select(-ID)
df <- df %>% filter(two_y_mort!="NA")
deltarel_trans$two_y_mort
df$two_y_mort
#LOAD FILE
#load(file="COPD.TAXA.Compare.RData")
#BAL and UA
#df <- df[df$Sample.Type!="BKG",]

# get the means of all the variables
means<-df%>%
  group_by(two_y_mort)%>%
  summarise_all(list(mean), na.rm=TRUE )%>%
  gather("Variable", "Mean", -two_y_mort)%>%
  spread(two_y_mort, Mean)%>%
  rename("Mean_less_24"='less_24', "Mean_greater_24"="greater_24")

# get the standard deviation of all the variables
sds<-df%>%
  group_by(two_y_mort)%>%
  summarise_all(list(sd), na.rm=TRUE )%>%
  gather("Variable", "SD", -two_y_mort)%>%
  spread(two_y_mort, SD)%>%
  rename("SD_less_24"='less_24', "SD_greater_24"="greater_24")

# get the median of all the variables
medians<-df%>%
  group_by(two_y_mort)%>%
  summarise_all(list(median), na.rm=TRUE )%>%
  gather("Variable", "Median", -two_y_mort)%>%
  spread(two_y_mort, Median)%>%
  rename("Median_less_24"='less_24', "Median_greater_24"="greater_24")


# join the tables 
summary_report<-means%>%
  inner_join(sds, "Variable")%>%
  inner_join(medians, "Variable")%>%
  mutate(DiffInMedians=Median_less_24-Median_greater_24)

############ perform Mann Whitney test 
# Select variables of interest
df2 <- df %>% select(-two_y_mort)
variables <- colnames(df2)

# Initialize vectors to store results
pvals <- numeric(length(variables))
vars <- character(length(variables))

# Loop through variables
for (i in seq_along(variables)) {
  
  # Extract data for the Mann-Whitney U test
  data <- df %>% select(two_y_mort, variables[i])
  
  # Extract samples for BAL and UA
  x1 <- data %>%
    filter(two_y_mort == "less_24") %>%
    pull(variables[i]) %>%
    na.omit()
  
  x2 <- data %>%
    filter(two_y_mort == "greater_24") %>%
    pull(variables[i]) %>%
    na.omit()
  
  # Perform the Mann-Whitney U test
  wc <- wilcox.test(x1, x2)
  
  # Store results
  pvals[i] <- round(wc$p.value, 4)
  vars[i] <- variables[i]
}


wc_df<-data.frame(Variable=vars,pvalues=pvals)

wc_df$Variable<-as.character(wc_df$Variable)
#Create Summary Report
summary_report<-summary_report%>%inner_join(wc_df, by="Variable")
#Order Summary Report
summary_report <- summary_report[order(summary_report$DiffInMedians),]
#Subset Significant Only
summary_report_sig <- summary_report[summary_report$pvalues<0.05,]
#Order by Median Diff
summary_report_sig <- summary_report_sig[order(summary_report_sig$DiffInMedians),]
summary_report_sig$order <- 1:nrow(summary_report_sig)



####### prepare data for figures #######
summary_report_sig_c <- summary_report_sig

#different Ordering by group 1
summary_report_sig_c <- summary_report_sig_c[order(-summary_report_sig_c$Median_less_24),]
summary_report_sig_c$order <- 1:nrow(summary_report_sig_c)

#different Ordering by group 2
summary_report_sig_c <- summary_report_sig_c[order(-summary_report_sig_c$Median_greater_24),]
summary_report_sig_c$order2 <- 1:nrow(summary_report_sig_c)

deltaorder <- summary_report_sig_c %>% select(Variable,order,order2)
deltaorder <- as.data.frame(deltaorder)
deltaorder <- deltaorder[!is.na(deltaorder$Variable),]
rownames(deltaorder) <- deltaorder$Variable
deltaorder$taxa <- rownames(deltaorder)

#Select bottom 25
#deltaorderbottom <- deltaorder[deltaorder$order>25,]
#deltaorderbottom <- deltaorderbottom[order(deltaorder$order2),]
#deltaorderbottom <- deltaorderbottom %>% slice_head(n=25)
#deltaorderbottom$order3 <- 1:nrow(deltaorderbottom)
#deltaorderbottom <- deltaorderbottom[deltaorderbottom$order3<10,]
#deltaorderbottom <- deltaorderbottom %>% select(!order3)
#Select top 25
#deltaordertop <- deltaorder[order(deltaorder$order),]
#deltaordertop <- deltaordertop %>% slice_head(n=25)
#merge again
#deltaorder <- rbind(deltaordertop,deltaorderbottom)

#Set new order by the leading group if few pathways 
#deltaorder$neworder <- 1:nrow(deltaorder)
deltaorder$neworder <- deltaorder$order
deltaorder <- deltaorder %>% arrange(neworder)

#Get the rel abundance tables
rows4 <- which(rownames(deltarel) %in% rownames(deltaorder))
deltarel_sig <-deltarel[rows4,]
#get Metadata of interest
deltameta <- coldata %>% select(ID,two_y_mort)
deltameta$ID <- gsub("_", ".", deltameta$ID)
#Transpose Data
deltarel_sig$taxa <- rownames(deltarel_sig)
ncol(deltarel_sig)
#select only columns of samples 
deltarel_sig2 <- deltarel_sig %>% tidyr::gather(ID, abundance,1:76)
#merge MetaData
deltarel_sig3 <- merge(deltarel_sig2,deltameta,by="ID",all.x=TRUE)
#merge order
deltaorder2 <- deltaorder %>% select(taxa,neworder)
deltarel_sig3 <- merge(deltarel_sig3,deltaorder2, by="taxa")
#Get contam Data
#contam <- read.delim2("rank_taxa.combined.txt", sep="\t")
##Merge with DESEQ data
#deltarel_sig3 <- merge(deltarel_sig3,contam,by="taxa",all.x=TRUE)
##Remove Contams
#deltarel_sig3 <- deltarel_sig3[deltarel_sig3$probable_metatranscriptome_contaminant=="FALSE",]
#deltarel_sig3 <- deltarel_sig3[deltarel_sig3$probable_metagenome_contaminant=="FALSE",]
#Remove BKG
#deltarel_sig3 <- deltarel_sig3[deltarel_sig3$two_y_mort %in% c("PD","non_PD"),]

#set Order of Figure
deltarel_sig3 <- deltarel_sig3 %>% filter(!is.na(two_y_mort))
#deltarel_sig3$two_y_mort <- factor(deltarel_sig3$two_y_mort, levels = c("non_PD", "PD"))
deltarel_sig3 <- deltarel_sig3 %>% arrange(neworder)

#Set division of Figure by selecting to p3 pathways. in this case we will select top 4 pathways for PD
deltarel_sig3$abundance
#deltarel_sig3$top <- ifelse(deltarel_sig3$neworder>25,2,1)
deltarel_sig3$top <- ifelse(deltarel_sig3$neworder>5,"Top_greater_24","Top_less_24")
#set order of top to adapt for the figure 
deltarel_sig3$top <- factor(deltarel_sig3$top, levels = c("Top_less_24","Top_greater_24"))

#Need to remove some of the modules with very low abundance
#deltarel_sss1 <- deltarel_sig3[deltarel_sig3$neworder<7,]
#deltarel_sss2 <- deltarel_sig3[deltarel_sig3$neworder>25,]
#deltarel_sss2 <- deltarel_sss2[deltarel_sss2$neworder<30,]
#create a spare in case
#spare <- deltarel_sig3
#Bind them back together
#deltarel_sig3 <- rbind(deltarel_sss1,deltarel_sss2)

#Get relabundance tables for each group 
library(tidyverse)
#for PD gorup 
rel.otuC1 <-               
  metat_counts %>% dplyr::select(rownames(metat_group_1_coldata)) %>% 
  rownames_to_column('gs') %>%
  group_by(gs) %>% 
  summarise_all(funs(sum)) %>%
  mutate_if(is.numeric, funs(./sum(.))) %>%
  column_to_rownames('gs')
#for non PD group 
rel.otuC2 <-
  metat_counts %>% dplyr::select(rownames(metat_group_2_coldata)) %>% 
  rownames_to_column('gs') %>%
  group_by(gs) %>% 
  summarise_all(funs(sum)) %>%
  mutate_if(is.numeric, funs(./sum(.))) %>%
  column_to_rownames('gs')
rel.otuC1 <- rel.otuC1*100
rel.otuC2 <- rel.otuC2*100

#decide what otu to save 
deltarel_sig3_non_PD <- deltarel_sig3[deltarel_sig3$two_y_mort=="greater_24",]
deltarel_sig3_PD <- deltarel_sig3[deltarel_sig3$two_y_mort=="less_24",]
otu.to.save1 <-as.character(deltarel_sig3_non_PD$taxa)
otu.to.save2 <-as.character(deltarel_sig3_PD$taxa)
#convert to dataframe
df.1.df <- data.frame(rel.otuC1)
df.2.df <- data.frame(rel.otuC2)
#from relative table we should get the mean across the row of the otu table
df.1.meanRA <- rowMeans(df.1.df)
df.2.meanRA <- rowMeans(df.2.df)
#need to subset AND reorder just the otus that we have 
df.1.meanRA.save <- df.1.meanRA[otu.to.save1]
df.2.meanRA.save <- df.2.meanRA[otu.to.save2]
#add the abundnace data for the res dataframe
deltarel_sig3_non_PD$relabundance <- df.1.meanRA.save
deltarel_sig3_PD$relabundance <- df.2.meanRA.save

#Keep abundant ones
#deltarel_sig3_ua <- deltarel_sig3_ua[deltarel_sig3_ua$relabundance>0.2,]
#deltarel_sig3_bal <- deltarel_sig3_bal[deltarel_sig3_ua$relabundance>0.2,]
#Barplot data
deltarel_sig3_non_PD_bar <- deltarel_sig3_non_PD %>% distinct(taxa, .keep_all=TRUE)
deltarel_sig3_PD_bar <- deltarel_sig3_PD %>% distinct(taxa, .keep_all=TRUE)
barplot  <- rbind(deltarel_sig3_non_PD_bar,deltarel_sig3_PD_bar)



p1 <- ggplot(barplot, aes(x=reorder(taxa,-neworder), y=relabundance, fill=two_y_mort)) + 
  #stat_boxplot(geom ='errorbar', width=0.1)+
  geom_bar(stat="identity", aes(fill=two_y_mort), width=0.5,position = position_dodge(width = 0.8))+
  #xmin=as.numeric(id)-.4,xmax=as.numeric(id)+.4, x=id, ymin=10E-7, ymax=ymean, fill=var)
  #geom_rect(position=position_dodge(.8), aes(xmin=as.numeric(reorder(taxa,+order))-.4,xmax=as.numeric(reorder(taxa,+order))+.4, ymin=0, ymax=relabundance,fill=Sample.Type))+
  #geom_rect(aes(xmin=as.numeric(reorder(taxa, +order))-0.4,xmax=as.numeric(reorder(taxa, +order))+.4, ymax=relabundance, ymin=0, fill=Sample.Type),color="black",position=position_dodge(0.9)) +
  #geom_boxplot(aes(color=Sample.Type),outlier.shape = NA, width=0.5)+
  #geom_boxplot(aes(ymin=..lower.., ymax=..upper..))+
  #geom_point(aes(color=Sample.Type),position = position_jitter(width = 0.25, seed = 123), alpha=0.7)+
  #geom_line(aes(group = SampleID),position = position_jitter(width = 0.25, seed = 123),linetype = "dashed",color="grey")+
  #geom_jitter(aes(color=Subject_Type_2),shape=1, position=position_jitter(0.2))+
  #geom_line(aes(group = SampleID), position=position_jitter(0.2),linetype = "dashed")+
  #facet_wrap(top ~ . , nrow=2,ncol=1)+
  #facet_grid(top ~., scales = "free_y", space = "free_y") 
facet_grid(top ~ . , scales = "free", space = "free")+
  #scale_color_manual(values=c("Asymptomatic_SC"="#4DAC26","Symptomatic_SC"="#FFD479","COPD"="#D01C8B")) +
  scale_fill_manual(values=c("less_24"="red3","greater_24"="orange3"),guide=FALSE) +
  #scale_x_discrete(labels = c('BKG','BPT','SPT','Sup'))+ 
  #geom_text_repel(aes(label=ifelse(res$sig > 3 , as.character(res$o),'')),size=3,force=25) +
  #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
  #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
  #scale_y_continuous(name="Log10 Rel Ab",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
  #scale_shape_manual(values=c("MG"=17,"MT"=16), guide=guide_legend(title="Sequencing",size=c(3, 3)))+
  #ylab("[MT-MG]") + 
  xlab("")+
  coord_flip()+
  #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
  #geom_point(color=cols) +
  #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
  theme4

p1

p2 <-    ggplot(deltarel_sig3_non_PD, aes(x=reorder(taxa,-neworder), y=abundance)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(aes(color=two_y_mort),outlier.shape = NA, width=0.5)+
  #geom_boxplot(aes(ymin=..lower.., ymax=..upper..))+
  geom_point(aes(color=two_y_mort),position = position_jitter(width = 0.25, seed = 123), alpha=0.7)+
  coord_flip()+
  #geom_line(aes(group = SampleID),position = position_jitter(width = 0.25, seed = 123),linetype = "dashed",color="grey")+
  #geom_jitter(aes(color=Subject_Type_2),shape=1, position=position_jitter(0.2))+
  #geom_line(aes(group = SampleID), position=position_jitter(0.2),linetype = "dashed")+
  #facet_wrap(~ Sample.Type, ncol=3)+
  facet_grid(top ~ ., scales = "free_y", space = "free")+
  #scale_color_manual(values=c("Asymptomatic_SC"="#4DAC26","Symptomatic_SC"="#FFD479","COPD"="#D01C8B")) +
  scale_color_manual(values=c("less_24"="red3","greater_24"="orange3"),guide=FALSE) +
  #scale_x_discrete(labels = c('BKG','BPT','SPT','Sup'))+ 
  #geom_text_repel(aes(label=ifelse(res$sig > 3 , as.character(res$o),'')),size=3,force=25) +
  #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
  #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
  #scale_y_continuous(name="Delta MT-MG",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
  #scale_shape_manual(values=c("MG"=17,"MT"=16), guide=guide_legend(title="Sequencing",size=c(3, 3)))+
  ylab("[MT-MG]") + 
  xlab("")+
  scale_y_continuous(limits = c(-10, 10))+
  #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
  #geom_point(color=cols) +
  #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
  theme6

p3 <-    ggplot(deltarel_sig3_PD, aes(x=reorder(taxa,-neworder), y=abundance)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(aes(color=two_y_mort),outlier.shape = NA, width=0.5)+
  #geom_boxplot(aes(ymin=..lower.., ymax=..upper..))+
  geom_point(aes(color=two_y_mort),position = position_jitter(width = 0.25, seed = 123), alpha=0.7)+
  coord_flip()+
  #geom_line(aes(group = SampleID),position = position_jitter(width = 0.25, seed = 123),linetype = "dashed",color="grey")+
  #geom_jitter(aes(color=Subject_Type_2),shape=1, position=position_jitter(0.2))+
  #geom_line(aes(group = SampleID), position=position_jitter(0.2),linetype = "dashed")+
  #facet_wrap(~ Sample.Type, ncol=3)+
  facet_grid(top ~ ., scales = "free_y", space = "free")+
  #scale_color_manual(values=c("Asymptomatic_SC"="#4DAC26","Symptomatic_SC"="#FFD479","COPD"="#D01C8B")) +
  scale_color_manual(values=c("less_24"="red3","greater_24"="orange3"),guide=FALSE) +
  #scale_x_discrete(labels = c('BKG','BPT','SPT','Sup'))+ 
  #geom_text_repel(aes(label=ifelse(res$sig > 3 , as.character(res$o),'')),size=3,force=25) +
  #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
  #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
  #scale_y_continuous(name="Delta MT-MG",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
  #scale_shape_manual(values=c("MG"=17,"MT"=16), guide=guide_legend(title="Sequencing",size=c(3, 3)))+
  ylab("[MT-MG]") + 
  xlab("")+
  scale_y_continuous(limits = c(-10, 10))+
  #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
  #geom_point(color=cols) +
  #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
  theme5
dummy <- ggplot(data = deltarel_sig3_non_PD, aes(x=reorder(taxa,-neworder), y=abundance))+
  facet_grid(top ~ ., scales = "free_y", space = "free")+
  geom_rect(aes(fill=size), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme_minimal()

library(cowplot)

pdf("two_y_mort_.DELTA_sig_abundance_Functional_Module.pdf", height = 6, width = 16)
#plot_grid(p1,p3,p2, labels="", rel_widths=c(1.1,1.5,1.5),ncol=3)
plot_grid(p1,p3,p2, labels="", ncol=3)
dev.off()


################ k means MT clusters 3 clusters ##########

#//////////////////////////////////
#------------------->METAGENOME DataSet
#//////////////////////////////////
# load metadata 
coldata <- read.csv(file = "MG_metabolic_metadata.csv")
rownames(coldata) <- coldata$X #Metagenome
coldata <- coldata[,-1]
#add k means clusters 
k_means_data <- read.csv(file = "k_means_3_clusters_MT.csv")
rownames(k_means_data) <- k_means_data$X
k_means_data <- k_means_data %>% dplyr::select(k_means_3_cluster)
coldata <- merge(coldata, k_means_data, by="row.names", all=TRUE)
rownames(coldata) <- coldata$Row.names
coldata <- coldata[,-1]

# load metagenome dataset 
mycounts <-read.delim2("MG_metabolic_functional_counts.txt", sep="\t", row.names=1)
#mycounts <-read.delim2("counts.metagenome.ko.module.fmap.txt", sep="\t", row.names=1)
colnames(mycounts) <- gsub("X","",colnames(mycounts))

#Create Relative Abundance Table
rel <-as.data.frame(mycounts) #make counts dataframe
relcounts <-
  rel %>% 
  tibble::rownames_to_column('gs') %>%
  group_by(gs) %>% 
  summarise_all(funs(sum)) %>%
  mutate_if(is.numeric, funs(./sum(.))) %>%
  tibble::column_to_rownames('gs')

#Convert any NAs to 0
mycounts[is.na(mycounts)] <- 0
relcounts[is.na(relcounts)] <- 0
#Copy of Count Table
mycounts2 <- mycounts
relcounts2 <- relcounts
#Convert Count Table into a Numeic Data Frame
d1 = data.frame(lapply(mycounts2, function(x) as.numeric(as.character(x))),
                check.names=F, row.names = rownames(mycounts2))
d2 = data.frame(lapply(relcounts2, function(x) as.numeric(as.character(x))),
                check.names=F, row.names = rownames(relcounts2))
#Remove 0 counts in relative Table
d2 <- d2[rowSums(d2[, -1] > 0) != 0, ]
d2 <- d2 %>% select(which(!colSums(d2, na.rm=TRUE) %in% 0))

#get the columns to match between count table and relative table 
wanted<-which(colnames(d1) %in% colnames(d2))
wanted2<-which(rownames(coldata) %in% colnames(d2))
d1<-d1[,wanted]
coldata<-coldata[wanted2,]

#df <- relcounts3
#Subset group 1 
metag_group_1_coldata <- coldata %>% filter(k_means_3_cluster=="1")
# now subset thier relative data and counts data 
metg_group_1_counts <- d1[,rownames(metag_group_1_coldata)]
metg_group_1_rel_counts <- d2[,rownames(metag_group_1_coldata)]

#subset group 2 
metag_group_2_coldata <- coldata %>% filter(k_means_3_cluster=="2")
# now subset thier relative data and counts data 
metg_group_2_counts <- d1[,rownames(metag_group_2_coldata)]
metg_group_2_rel_counts <- d2[,rownames(metag_group_2_coldata)]

#subset group 3
metag_group_3_coldata <- coldata %>% filter(k_means_3_cluster=="3")
# now subset thier relative data and counts data 
metg_group_3_counts <- d1[,rownames(metag_group_3_coldata)]
metg_group_3_rel_counts <- d2[,rownames(metag_group_3_coldata)]

# all metagenome metadata 
metag_metadata <- coldata
# all metagenome counts 
metag_counts <- d1
# all metagenome relative counts 
metag_rel_counts <- d2


############# metatranscriptome dataset
# load metadata 
coldata <- read.csv(file = "MT_metabolic_metadata.csv")
rownames(coldata) <- coldata$X 
coldata <- coldata[,-1]

#add k means clusters 
k_means_data <- read.csv(file = "k_means_3_clusters_MT.csv")
rownames(k_means_data) <- k_means_data$X
k_means_data <- k_means_data %>% dplyr::select(k_means_3_cluster)
coldata <- merge(coldata, k_means_data, by="row.names", all=TRUE)
rownames(coldata) <- coldata$Row.names
coldata <- coldata[,-1]


# load metagenome dataset 
mycounts <-read.delim2("MT_metabolic_functional_counts.txt", sep="\t", row.names=1)
#mycounts <-read.delim2("counts.metagenome.ko.module.fmap.txt", sep="\t", row.names=1)
colnames(mycounts) <- gsub("X","",colnames(mycounts))

#Create Relative Abundance Table
rel <-as.data.frame(mycounts) #make counts dataframe
relcounts <-
  rel %>% 
  tibble::rownames_to_column('gs') %>%
  group_by(gs) %>% 
  summarise_all(funs(sum)) %>%
  mutate_if(is.numeric, funs(./sum(.))) %>%
  tibble::column_to_rownames('gs')

#Convert any NAs to 0
mycounts[is.na(mycounts)] <- 0
relcounts[is.na(relcounts)] <- 0
#Copy of Count Table
mycounts2 <- mycounts
relcounts2 <- relcounts
#Convert Count Table into a Numeic Data Frame
d1 = data.frame(lapply(mycounts2, function(x) as.numeric(as.character(x))),
                check.names=F, row.names = rownames(mycounts2))
d2 = data.frame(lapply(relcounts2, function(x) as.numeric(as.character(x))),
                check.names=F, row.names = rownames(relcounts2))
#Remove 0 counts in relative Table
d2 <- d2[rowSums(d2[, -1] > 0) != 0, ]
d2 <- d2 %>% select(which(!colSums(d2, na.rm=TRUE) %in% 0))

#get the columns to match between count table and relative table 
wanted<-which(colnames(d1) %in% colnames(d2))
wanted2<-which(rownames(coldata) %in% colnames(d2))
d1<-d1[,wanted]
coldata<-coldata[wanted2,]

#Subset group 1 
metat_group_1_coldata <- coldata %>% filter(k_means_3_cluster=="1")
# now subset thier relative data and counts data 
mett_group_1_counts <- d1[,rownames(metat_group_1_coldata)]
mett_group_1_rel_counts <- d2[,rownames(metat_group_1_coldata)]

#subset group 2 
metat_group_2_coldata <- coldata %>% filter(k_means_3_cluster=="2")
# now subset thier relative data and counts data 
mett_group_2_counts <- d1[,rownames(metat_group_2_coldata)]
mett_group_2_rel_counts <- d2[,rownames(metat_group_2_coldata)]

metat_group_3_coldata <- coldata %>% filter(k_means_3_cluster=="3")
# now subset thier relative data and counts data 
mett_group_3_counts <- d1[,rownames(metat_group_3_coldata)]
mett_group_3_rel_counts <- d2[,rownames(metat_group_3_coldata)]

# all metatranscriptome metadata 
metat_metadata <- coldata
# all metatrans riptome counts 
metat_counts <- d1
# all metatranscripotome relative counts 
metat_rel_counts <- d2




########### delta table ####### 
#match two tables of counts 
setdiff(rownames(metag_counts),rownames(metat_counts))

#set the same order 
#Order Tables
metag_counts <- metag_counts[,order(colnames(metag_counts))]
metag_counts <- metag_counts[order(rownames(metag_counts)),]
metat_counts <- metat_counts[,order(colnames(metat_counts))]
metat_counts <- metat_counts[order(rownames(metat_counts)),]

metag_rel_counts <- metag_rel_counts[,order(colnames(metag_rel_counts))]
metag_rel_counts <- metag_rel_counts[order(rownames(metag_rel_counts)),]
metat_rel_counts <- metat_rel_counts[,order(colnames(metat_rel_counts))]
metat_rel_counts <- metat_rel_counts[order(rownames(metat_rel_counts)),]

#Match Tables
rows1<-which(rownames(metag_counts)  %in%  rownames(metat_counts))
metag_counts_2<-metag_counts[rows1,]

rows2<-which(rownames(metat_counts)  %in%  rownames(metag_counts_2))
metat_counts_2<-metat_counts[rows2,]

cols1 <- which(colnames(metag_counts_2) %in% colnames(metat_counts_2))
metag_counts_3<-metag_counts_2[,cols1]
metat_counts_3<-metat_counts_2

#match relative tables 
rows1<-which(rownames(metag_rel_counts)  %in%  rownames(metat_rel_counts))
metag_rel_counts_2<-metag_rel_counts[rows1,]

rows2<-which(rownames(metat_rel_counts)  %in%  rownames(metag_rel_counts_2))
metat_rel_counts_2<-metat_rel_counts[rows2,]

cols1 <- which(colnames(metag_rel_counts_2) %in% colnames(metat_rel_counts_2))
metag_rel_counts_3<-metag_rel_counts_2[,cols1]
metat_rel_counts_3<-metat_rel_counts_2

#add one
metag_counts_4 <- metag_counts_3+1
metat_counts_4 <- metat_counts_3+1

metag_rel_counts_4 <- metag_rel_counts_3+1
metat_rel_counts_4 <- metat_rel_counts_3+1

#multiply by 100
metag_rel_counts_5 <- metag_rel_counts_3*100
metat_rel_counts_5 <- metat_rel_counts_3*100

#Create Ratio Table
ratio <- cbind(metat_counts_4[1],round(metat_counts_4[-1]/metag_counts_4[-1],1))

#Create Delta Table
delta <- metat_counts_4-metag_counts_4
#deltarel <- metatcountrel4-metagcountrel4
deltarel <- metat_rel_counts_5-metag_rel_counts_5

#Transpose
deltarel_trans <- as.data.frame(t(deltarel))
deltarel_trans$ID <- rownames(deltarel_trans)

#get Metadata of interest
deltarel_trans_meta <- coldata %>% select(ID,k_means_3_cluster)
setdiff(deltarel_trans_meta$ID,deltarel_trans$ID)
deltarel_trans_meta$ID <- gsub("_", ".", deltarel_trans_meta$ID)

#merge MetaData
deltarel_trans <- merge(deltarel_trans,deltarel_trans_meta,by="ID",all.x=TRUE)
deltarel_trans$k_means_3_cluster
# remove the ID number
df<-deltarel_trans%>%select(-ID)
#df <- df %>% filter(two_y_mort!="NA")
deltarel_trans$k_means_3_cluster
df$k_means_3_cluster
#LOAD FILE
#load(file="COPD.TAXA.Compare.RData")
#BAL and UA
#df <- df[df$Sample.Type!="BKG",]
library(tidyverse)
# get the means of all the variables
means<-df%>%
  group_by(k_means_3_cluster)%>%
  summarise_all(list(mean), na.rm=TRUE )%>%
  gather("Variable", "Mean", -k_means_3_cluster)%>%
  spread(k_means_3_cluster, Mean)%>%
  rename("Mean_1"='1', "Mean_2"="2", "Mean_3"="3")

# get the standard deviation of all the variables
sds<-df%>%
  group_by(k_means_3_cluster)%>%
  summarise_all(list(sd), na.rm=TRUE )%>%
  gather("Variable", "SD", -k_means_3_cluster)%>%
  spread(k_means_3_cluster, SD)%>%
  rename("SD_1"='1', "SD_2"="2", "SD_3"="3")

# get the median of all the variables
medians<-df%>%
  group_by(k_means_3_cluster)%>%
  summarise_all(list(median), na.rm=TRUE )%>%
  gather("Variable", "Median", -k_means_3_cluster)%>%
  spread(k_means_3_cluster, Median)%>%
  rename("Median_1"='1', "Median_2"="2", "Median_3"="3")


# join the tables 
summary_report<-means%>%
  inner_join(sds, "Variable")%>%
  inner_join(medians, "Variable")%>%
  mutate(Diff_1_2 = `Median_1` - `Median_2`,
         Diff_1_3 = `Median_1` - `Median_3`,
         Diff_2_3 = `Median_2` - `Median_3`)

############ perform Mann Whitney test 
# Select variables of interest
df2 <- df %>% select(-k_means_3_cluster)
variables <- colnames(df2)

# Initialize a data frame to store results
results <- data.frame(
  Variable = character(0),
  Comparison = character(0),
  P_Value = numeric(0)
)

# Loop through variables and perform pairwise comparisons
for (var in variables) {
  
  # Extract data for the variable
  data <- df %>% select(k_means_3_cluster, all_of(var))
  
  # Pairwise comparisons
  for (pair in list(c(1, 2), c(2, 3), c(1, 3))) {
    group1 <- pair[1]
    group2 <- pair[2]
    
    # Extract values for the groups
    x1 <- data %>%
      filter(k_means_3_cluster == group1) %>%
      pull(var) %>%
      na.omit()
    
    x2 <- data %>%
      filter(k_means_3_cluster == group2) %>%
      pull(var) %>%
      na.omit()
    
    # Perform Mann-Whitney U test
    if (length(x1) > 1 && length(x2) > 1) { # Ensure there are enough data points
      wc <- wilcox.test(x1, x2)
      p_value <- round(wc$p.value, 4)
    } else {
      p_value <- NA  # Not enough data for testing
    }
    
    # Store the results
    results <- rbind(results, data.frame(
      Variable = var,
      Comparison = paste(group1, "vs", group2),
      P_Value = p_value
    ))
  }
}

# View the results
print(results)


wc_df<-results

wc_df$Variable<-as.character(wc_df$Variable)
#Create Summary Report
summary_report<-summary_report%>%inner_join(wc_df, by="Variable")

#export 
write.csv(summary_report, file = "MT_clusters_delta_results.csv")


########### prepare for plots ##########
#Subset Significant Only
summary_report_sig <- summary_report[summary_report$P_Value<0.05,]
summary_report_sig_c <- summary_report_sig

#different Ordering by group 1
summary_report_sig_c <- summary_report_sig_c[order(-summary_report_sig_c$Median_1),]
summary_report_sig_c$order_1 <- 1:nrow(summary_report_sig_c)

#different Ordering by group 2
summary_report_sig_c <- summary_report_sig_c[order(-summary_report_sig_c$Median_2),]
summary_report_sig_c$order_2 <- 1:nrow(summary_report_sig_c)

#different Ordering by group 3
summary_report_sig_c <- summary_report_sig_c[order(-summary_report_sig_c$Median_3),]
summary_report_sig_c$order_3 <- 1:nrow(summary_report_sig_c)

summary_report_sig_c


deltaorder <- summary_report_sig_c %>% select(Variable,order_1,order_2, order_3)
deltaorder <- as.data.frame(deltaorder)
deltaorder <- deltaorder[!is.na(deltaorder$Variable),]
#rownames(deltaorder) <- deltaorder$Variable
deltaorder$taxa <- deltaorder$Variable

deltaorder

#Select bottom and top 25 pathways according to each cluster 

# Select top and bottom pathways for each cluster
top_n <- 25  # Number of top pathways
bottom_n <- 25  # Number of bottom pathways

top_cluster_1 <- deltaorder %>%
  arrange(order_1) %>%
  slice_head(n = top_n)

bottom_cluster_1 <- deltaorder %>%
  arrange(desc(order_1)) %>%
  slice_head(n = bottom_n)

top_cluster_2 <- deltaorder %>%
  arrange(order_2) %>%
  slice_head(n = top_n)

bottom_cluster_2 <- deltaorder %>%
  arrange(desc(order_2)) %>%
  slice_head(n = bottom_n)

top_cluster_3 <- deltaorder %>%
  arrange(order_3) %>%
  slice_head(n = top_n)

bottom_cluster_3 <- deltaorder %>%
  arrange(desc(order_3)) %>%
  slice_head(n = bottom_n)

# Combine all selected pathways
deltaorder <- bind_rows(
  top_cluster_1, bottom_cluster_1,
  top_cluster_2, bottom_cluster_2,
  top_cluster_3, bottom_cluster_3
) %>%
  arrange(order_1) %>% 
  #distinct(Variable, .keep_all = TRUE) %>%
  mutate(neworder = row_number())

# now create top variable for each clsuter and will plot only top pathways for each cluster 
deltaorder <- deltaorder %>%
  mutate(top = case_when(
      order_1 <= 50 ~ "Top_Cluster_1",
      order_2 <= 50 ~ "Top_Cluster_2",
      order_3 <= 50 ~ "Top_Cluster_3"))
#remove those that are not top clusters 
deltaorder <- deltaorder %>% filter(!is.na(top)) 

#now set new order 
deltaorder <- deltaorder %>% mutate(neworder=row_number())


#Get the rel abundance tables
rows4 <- which(rownames(deltarel) %in% deltaorder_group_1$Variable)
deltarel_sig <-deltarel[rows4,]
#get Metadata of interest
deltameta <- coldata %>% select(ID,k_means_3_cluster)
deltameta$ID <- gsub("_", ".", deltameta$ID)
#Transpose Data
deltarel_sig$taxa <- rownames(deltarel_sig)
ncol(deltarel_sig)
#select only columns of samples 
deltarel_sig2 <- deltarel_sig %>% tidyr::gather(ID, abundance,1:76)
#merge MetaData
deltarel_sig3 <- merge(deltarel_sig2,deltameta,by="ID",all.x=TRUE)
#merge order
deltaorder2 <- deltaorder %>% select(taxa,neworder, top)
deltarel_sig3 <- merge(deltarel_sig3,deltaorder2, by="taxa")

#set Order of Figure
deltarel_sig3 <- deltarel_sig3 %>% arrange(neworder)
deltarel_sig3$abundance

#set order of top to adapt for the figure 
deltarel_sig3$top <- factor(deltarel_sig3$top, levels = c("Top_Cluster_1","Top_Cluster_2", "Top_Cluster_3"))


#Get relabundance tables for each group 
library(tidyverse)
#for cluster 1 gorup 
rel.otuC1 <-               
  metat_counts %>% dplyr::select(rownames(metat_group_1_coldata)) %>% 
  rownames_to_column('gs') %>%
  group_by(gs) %>% 
  summarise_all(funs(sum)) %>%
  mutate_if(is.numeric, funs(./sum(.))) %>%
  column_to_rownames('gs')
#for cluster 2 group 
rel.otuC2 <-
  metat_counts %>% dplyr::select(rownames(metat_group_2_coldata)) %>% 
  rownames_to_column('gs') %>%
  group_by(gs) %>% 
  summarise_all(funs(sum)) %>%
  mutate_if(is.numeric, funs(./sum(.))) %>%
  column_to_rownames('gs')
# for cluster 3 group
rel.otuC3 <-
  metat_counts %>% dplyr::select(rownames(metat_group_3_coldata)) %>% 
  rownames_to_column('gs') %>%
  group_by(gs) %>% 
  summarise_all(funs(sum)) %>%
  mutate_if(is.numeric, funs(./sum(.))) %>%
  column_to_rownames('gs')

rel.otuC1 <- rel.otuC1*100
rel.otuC2 <- rel.otuC2*100
rel.otuC3 <- rel.otuC3*100

#decide what taxa to save 
deltarel_sig3_1 <- deltarel_sig3[deltarel_sig3$k_means_3_cluster=="1",]
deltarel_sig3_2 <- deltarel_sig3[deltarel_sig3$k_means_3_cluster=="2",]
deltarel_sig3_3 <- deltarel_sig3[deltarel_sig3$k_means_3_cluster=="3",]

otu.to.save1 <-as.character(deltarel_sig3_1$taxa)
otu.to.save2 <-as.character(deltarel_sig3_2$taxa)
otu.to.save3 <-as.character(deltarel_sig3_3$taxa)

#convert to dataframe
df.1.df <- data.frame(rel.otuC1)
df.2.df <- data.frame(rel.otuC2)
df.3.df <- data.frame(rel.otuC3)

#from relative table we should get the mean across the row of the otu table
df.1.meanRA <- rowMeans(df.1.df)
df.2.meanRA <- rowMeans(df.2.df)
df.3.meanRA <- rowMeans(df.3.df)

#need to subset AND reorder just the otus that we have 
df.1.meanRA.save <- df.1.meanRA[otu.to.save1]
df.2.meanRA.save <- df.2.meanRA[otu.to.save2]
df.3.meanRA.save <- df.3.meanRA[otu.to.save3]

#add the abundnace data for the res dataframe
deltarel_sig3_1$relabundance <- df.1.meanRA.save
deltarel_sig3_2$relabundance <- df.2.meanRA.save
deltarel_sig3_3$relabundance <- df.3.meanRA.save

#Barplot data
deltarel_sig3_1_bar <- deltarel_sig3_1 %>% distinct(taxa, .keep_all=TRUE)
deltarel_sig3_2_bar <- deltarel_sig3_2 %>% distinct(taxa, .keep_all=TRUE)
deltarel_sig3_3_bar <- deltarel_sig3_3 %>% distinct(taxa, .keep_all=TRUE)

barplot  <- rbind(deltarel_sig3_1_bar,deltarel_sig3_2_bar, deltarel_sig3_3_bar)
barplot$k_means_3_cluster <- factor(barplot$k_means_3_cluster, levels = c("3", "2", "1"))
#remove map annotation 
barplot$taxa_new <- sub("^map\\d{5}:\\s*", "", barplot$taxa)

#keep main pahtwyas only 
barplot_filtered <- barplot %>%
  filter(!taxa_new %in% c("Tryptophan metabolism", "Carbon metabolism"))

p2_dataframe$taxa
p2_dataframe <- rbind(deltarel_sig3_1, deltarel_sig3_2, deltarel_sig3_3)
p2_dataframe_filtered <- p2_dataframe %>% filter(!taxa%in%c("map01200: Carbon metabolism", "map00380: Tryptophan metabolism"))

p1 <- ggplot(barplot_filtered, aes(x=reorder(taxa_new,-neworder), y=relabundance, fill=k_means_3_cluster)) + 
  #stat_boxplot(geom ='errorbar', width=0.1)+
  geom_bar(stat="identity", aes(fill=k_means_3_cluster), width=0.5,position = position_dodge(width = 0.8))+
  #xmin=as.numeric(id)-.4,xmax=as.numeric(id)+.4, x=id, ymin=10E-7, ymax=ymean, fill=var)
  #geom_rect(position=position_dodge(.8), aes(xmin=as.numeric(reorder(taxa,+order))-.4,xmax=as.numeric(reorder(taxa,+order))+.4, ymin=0, ymax=relabundance,fill=Sample.Type))+
  #geom_rect(aes(xmin=as.numeric(reorder(taxa, +order))-0.4,xmax=as.numeric(reorder(taxa, +order))+.4, ymax=relabundance, ymin=0, fill=Sample.Type),color="black",position=position_dodge(0.9)) +
  #geom_boxplot(aes(color=Sample.Type),outlier.shape = NA, width=0.5)+
  #geom_boxplot(aes(ymin=..lower.., ymax=..upper..))+
  #geom_point(aes(color=Sample.Type),position = position_jitter(width = 0.25, seed = 123), alpha=0.7)+
  #geom_line(aes(group = SampleID),position = position_jitter(width = 0.25, seed = 123),linetype = "dashed",color="grey")+
  #geom_jitter(aes(color=Subject_Type_2),shape=1, position=position_jitter(0.2))+
  #geom_line(aes(group = SampleID), position=position_jitter(0.2),linetype = "dashed")+
  #facet_wrap(top ~ . , nrow=2,ncol=1)+
  #facet_grid(top ~., scales = "free_y", space = "free_y") 
facet_grid(top ~ . , scales = "free", space = "free")+
  #scale_color_manual(values=c("Asymptomatic_SC"="#4DAC26","Symptomatic_SC"="#FFD479","COPD"="#D01C8B")) +
  scale_fill_manual(values=c("1"="#F8766D","2"="#00BA38", "3"="#619CFF"),guide=FALSE) +
  #scale_x_discrete(labels = c('BKG','BPT','SPT','Sup'))+ 
  #geom_text_repel(aes(label=ifelse(res$sig > 3 , as.character(res$o),'')),size=3,force=25) +
  #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
  #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
  #scale_y_continuous(name="Log10 Rel Ab",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
  #scale_shape_manual(values=c("MG"=17,"MT"=16), guide=guide_legend(title="Sequencing",size=c(3, 3)))+
  #ylab("[MT-MG]") + 
  xlab("")+
  coord_flip()+
  #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
  #geom_point(color=cols) +
  #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
  theme4

p1 <- p1+theme(axis.text.y=element_text(size = 40, face = "bold"), 
               axis.text.x=element_text(size = 30, face = "bold"))

deltarel_sig3_1$k_means_3_cluster <- factor(deltarel_sig3_1$k_means_3_cluster)
deltarel_sig3_2$k_means_3_cluster <- factor(deltarel_sig3_2$k_means_3_cluster)
deltarel_sig3_3$k_means_3_cluster <- factor(deltarel_sig3_3$k_means_3_cluster)

p2_dataframe_filtered$k_means_3_cluster <- factor(p2_dataframe_filtered$k_means_3_cluster, levels = c("3", "2", "1"))

p_extra <- ggplot(p2_dataframe_filtered, aes(x=reorder(taxa,-neworder), y=abundance)) + 
  geom_point(aes(color = k_means_3_cluster), 
             position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.6, seed = 123), 
             alpha = 0.7, 
             size = 1.5) +  # INCREASE POINT SIZE
  stat_boxplot(aes(color=k_means_3_cluster),  
               geom = 'errorbar', 
               width = 0.2,  # INCREASE ERROR BAR WIDTH
               position = position_dodge(width = 0.6)) + 
  geom_boxplot(aes(color=k_means_3_cluster), 
               outlier.shape = NA, 
               width = 0.7,  # INCREASE BOX WIDTH
               position = position_dodge(width = 0.6), 
               size = 1.2) +  # INCREASE BOX BORDER SIZE
  coord_flip() +
  facet_grid(top ~ ., scales = "free_y", space = "free") + 
  scale_color_manual(values=c("1"="#F8766D","2"="#00BA38", "3"="#619CFF"), guide=FALSE) + 
  ylab("[MT-MG]") + 
  xlab("") +
  scale_y_continuous(limits = c(-10, 20)) +
  theme6

p_extra
p_extra <- p_extra+theme(axis.text.x=element_text(size = 30, face = "bold"))

p2 <-    ggplot(deltarel_sig3_1, aes(x=reorder(taxa,-neworder), y=abundance)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(aes(color=k_means_3_cluster),outlier.shape = NA, width=0.5)+
  #geom_boxplot(aes(ymin=..lower.., ymax=..upper..))+
  geom_point(aes(color=k_means_3_cluster),position = position_jitter(width = 0.25, seed = 123), alpha=0.7)+
  coord_flip()+
  #geom_line(aes(group = SampleID),position = position_jitter(width = 0.25, seed = 123),linetype = "dashed",color="grey")+
  #geom_jitter(aes(color=Subject_Type_2),shape=1, position=position_jitter(0.2))+
  #geom_line(aes(group = SampleID), position=position_jitter(0.2),linetype = "dashed")+
  #facet_wrap(~ Sample.Type, ncol=3)+
  facet_grid(top ~ ., scales = "free_y", space = "free")+
  #scale_color_manual(values=c("Asymptomatic_SC"="#4DAC26","Symptomatic_SC"="#FFD479","COPD"="#D01C8B")) +
  scale_color_manual(values=c("1"="#F8766D","2"="#00BA38", "3"="#619CFF"),guide=FALSE) +
  #scale_x_discrete(labels = c('BKG','BPT','SPT','Sup'))+ 
  #geom_text_repel(aes(label=ifelse(res$sig > 3 , as.character(res$o),'')),size=3,force=25) +
  #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
  #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
  #scale_y_continuous(name="Delta MT-MG",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
  #scale_shape_manual(values=c("MG"=17,"MT"=16), guide=guide_legend(title="Sequencing",size=c(3, 3)))+
  ylab("[MT-MG]") + 
  xlab("")+
  scale_y_continuous(limits = c(-10, 10))+
  #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
  #geom_point(color=cols) +
  #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
  theme6

p2

p3 <-    ggplot(deltarel_sig3_2, aes(x=reorder(taxa,-neworder), y=abundance)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(aes(color=k_means_3_cluster),outlier.shape = NA, width=0.5)+
  #geom_boxplot(aes(ymin=..lower.., ymax=..upper..))+
  geom_point(aes(color=k_means_3_cluster),position = position_jitter(width = 0.25, seed = 123), alpha=0.7)+
  coord_flip()+
  #geom_line(aes(group = SampleID),position = position_jitter(width = 0.25, seed = 123),linetype = "dashed",color="grey")+
  #geom_jitter(aes(color=Subject_Type_2),shape=1, position=position_jitter(0.2))+
  #geom_line(aes(group = SampleID), position=position_jitter(0.2),linetype = "dashed")+
  #facet_wrap(~ Sample.Type, ncol=3)+
  facet_grid(top ~ ., scales = "free_y", space = "free")+
  #scale_color_manual(values=c("Asymptomatic_SC"="#4DAC26","Symptomatic_SC"="#FFD479","COPD"="#D01C8B")) +
  scale_color_manual(values=c("1"="#F8766D","2"="#00BA38", "3"="#619CFF"),guide=FALSE) +
  #scale_x_discrete(labels = c('BKG','BPT','SPT','Sup'))+ 
  #geom_text_repel(aes(label=ifelse(res$sig > 3 , as.character(res$o),'')),size=3,force=25) +
  #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
  #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
  #scale_y_continuous(name="Delta MT-MG",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
  #scale_shape_manual(values=c("MG"=17,"MT"=16), guide=guide_legend(title="Sequencing",size=c(3, 3)))+
  ylab("[MT-MG]") + 
  xlab("")+
  scale_y_continuous(limits = c(-10, 10))+
  #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
  #geom_point(color=cols) +
  #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
  theme5

p3

p4 <-    ggplot(deltarel_sig3_3, aes(x=reorder(taxa,-neworder), y=abundance)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(aes(color=k_means_3_cluster),outlier.shape = NA, width=0.5)+
  #geom_boxplot(aes(ymin=..lower.., ymax=..upper..))+
  geom_point(aes(color=k_means_3_cluster),position = position_jitter(width = 0.25, seed = 123), alpha=0.7)+
  coord_flip()+
  #geom_line(aes(group = SampleID),position = position_jitter(width = 0.25, seed = 123),linetype = "dashed",color="grey")+
  #geom_jitter(aes(color=Subject_Type_2),shape=1, position=position_jitter(0.2))+
  #geom_line(aes(group = SampleID), position=position_jitter(0.2),linetype = "dashed")+
  #facet_wrap(~ Sample.Type, ncol=3)+
  facet_grid(top ~ ., scales = "free_y", space = "free")+
  #scale_color_manual(values=c("Asymptomatic_SC"="#4DAC26","Symptomatic_SC"="#FFD479","COPD"="#D01C8B")) +
  scale_color_manual(values=c("1"="#F8766D","2"="#00BA38", "3"="#619CFF"),guide=FALSE) +
  #scale_x_discrete(labels = c('BKG','BPT','SPT','Sup'))+ 
  #geom_text_repel(aes(label=ifelse(res$sig > 3 , as.character(res$o),'')),size=3,force=25) +
  #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
  #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
  #scale_y_continuous(name="Delta MT-MG",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
  #scale_shape_manual(values=c("MG"=17,"MT"=16), guide=guide_legend(title="Sequencing",size=c(3, 3)))+
  ylab("[MT-MG]") + 
  xlab("")+
  scale_y_continuous(limits = c(-10, 10))+
  #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
  #geom_point(color=cols) +
  #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
  theme5

p4



library(cowplot)

pdf("MT_k_means_cluster_.DELTA_sig_abundance_Functional_Module_one_figure_for_publication.pdf", height = 26, width = 34)
#plot_grid(p1,p3,p2, labels="", rel_widths=c(1.1,1.5,1.5),ncol=3)
plot_grid(p1,p_extra, labels="", ncol=2, rel_widths = c(1.2, 0.5))
dev.off()


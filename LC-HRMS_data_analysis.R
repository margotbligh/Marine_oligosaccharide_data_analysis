#Script written by Margot Bligh
#RData objects on Github, they can be directly loaded
#The pre-processing steps will take a long time to run

#Please note that for GlycoAnnotateR to run you MUST have python installed
#on your computer

#1: Import packages -----
library(tidyverse)
library(xcms)
library(data.table)
library(MSnbase)
library(GlycoAnnotateR)
library(ggpubr)
library(scales)
library(BiocParallel)
library(SummarizedExperiment)
library(viridis)
library(magrittr)
library(stringr)
library(ggrepel)
library(CAMERA)
library(patchwork)
library(factoextra)
library(cluster)
library(NbClust)
library(pheatmap)


#2: Set working directory-----
#CHANGE THIS TO YOUR WORKING DIRECTORY
dir <- "/Users/mbligh/ownCloud/marglyco/Data/LC-MS_data/GCC-extracts/20221004"
setwd(dir)

#3: Plotting functions----
#make function for plotting PCAs
plotPCA <- function(pcDf, pcSummary, pal = viridis(n = 3, option = 'E')){
  ggplot(pcDf, aes(x = PC1, y = PC2, shape = cruise, 
                   fill = sample_type))+
    geom_text_repel(pcDf %>% filter(sample == "MSM107-DCM2"),
                    mapping = aes(label = sample),
                    size = 3,
                    force = 10, segment.size = 0.2) +
    geom_point(size = 4, stroke = 0.5, colour = "black") +
    scale_fill_manual(values = pal) +
    scale_shape_manual(values = shapes) +
    labs(x = paste0("PC1: ", format(pcSummary$importance[2, 1] * 100,
                                    digits = 3), " % variance"),
         y = paste0("PC2: ", format(pcSummary$importance[2, 2] * 100,
                                    digits = 3), " % variance")) +
    guides(fill = guide_legend(override.aes = list(shape = 22))) +
    theme_light() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          panel.grid = element_blank(),
          legend.position = "none",
          text = element_text(family = "Arial"),
          axis.line = element_line(),
          axis.ticks = element_line(color = "black"),
          panel.grid.major=element_line(color="#ECECEC"),
          panel.border = element_blank(),
          panel.grid.minor = element_blank())
}

#4. Data import ####
#get file paths to mzML files
fp <- dir(path = "mzML", all.files = FALSE, full.names = TRUE)

#create phenodata data.frame
#each sample must have a unique name!
pd <- data.frame(name = basename(fp) %>%
                   gsub("MS31_20221004_|.mzML|_500ng", "", .) %>% 
                   sub("_\\d{2}$", "", .),
                 sample_type = basename(fp) %>% 
                   sub(".*MSM107.*|.*EN683.*", "cruise extract", .) %>% 
                   sub(".*blank.*", "solvent blank", .) %>%
                   sub(".*spiked.*", "spiked QC", .) %>% 
                   sub(".*_QC.*", "QC", .) %>% 
                   sub(".*std.*", "standard mix", .) %>% 
                   sub(".*cell.*|.*lam.*|.*starch.*", "standard", .),
                 cruise = basename(fp) %>% 
                   sub(".*MSM.*", "MSM107", .) %>% 
                   sub(".*EN683.*", "EN683", .) %>% 
                   sub(".*MS31_.*", "", .),
                 extract_type = basename(fp) %>% 
                   sub(".*ASW.*", "ASW", .) %>% 
                   sub(".*DCM.*", "DCM", .) %>% 
                   sub(".*deep.*", "deep", .) %>% 
                   sub(".*MS31_.*", "", .),
                 stringsAsFactors = FALSE)

#read in ms data with phenodata 
data <- readMSData(files = fp, pdata = new("NAnnotatedDataFrame", pd), 
                   mode = "onDisk")

#5. Peak picking-----
#create CentWaveParam object
cwp <- CentWaveParam()
#assign IPO optimized parameters
cwp@ppm <- 15.4
cwp@peakwidth <- c(3, 130)
cwp@mzdiff <- -0.02465

#find peaks using CentWave algorithm
data.peaks <- findChromPeaks(data, param = cwp)

#NB! data.peaks object can be downloaded from Github
#To proceed with the object on Github you will need to re-assign
#The directory name of the data as below:

#load data.peaks object 
load("data.peaks.RData")
#fix directoy name to the directory where you have the data
dirname(data.peaks) <- rep(paste0(dir, "/mzML"), length(fp))


#6. Peak merging-----
#make parameter object
mnp <- MergeNeighboringPeaksParam()
#assign IPO optimized parameters
mnp@expandRt <- 5
mnp@minProp <- 0.5

#merge neighbouring peaks
data.peaks.ref <- refineChromPeaks(data.peaks, param = mnp)


#NB! data.peaks.ref object can be downloaded from Github
#To proceed with the object on Github you will need to re-assign
#The directory name of the data as below:
load("data.peaks.ref.RData")
dirname(data.peaks.ref) <- rep(paste0(dir, "/mzML"), length(fp))

#7. Retention time alignment (subset based)-----
#change 'standard mix' to standard
#all standards will be set as one group
pd$sample_type[pd$sample_type == "standard mix"] <- "standard"
pData(data.peaks.ref) <- pd

#create peak grouping parameter object
pdp <- PeakDensityParam(sampleGroups = pd$sample_type, #set sample_type as grouping variable
                        bw = 6, minFraction = 2/3)

#initial peak grouping
data.peaks.ref.gsub <- groupChromPeaks(data.peaks.ref, param = pdp)

#NB: can also load object after downloading from Github
load("data.peaks.ref.gsub.RData")
dirname(data.peaks.ref.gsub) <- rep(paste0(dir, "/mzML"), length(fp))

#make object subset-alignment
pdp_subs <- PeakGroupsParam(minFraction = 2/3,
                            subset = which(data.peaks.ref.gsub$sample_type == "QC"),
                            subsetAdjust = "average", span = 0.4)
#perform alignment
data.peaks.ref.gsub.rt <- adjustRtime(data.peaks.ref.gsub, param = pdp_subs)

#NB! can also load object after downloading from Github
load("data.peaks.ref.gsub.rt.RData")
dirname(data.peaks.ref.gsub.rt) <- rep(paste0(dir, "/mzML"), length(fp))

#8: Peak grouping (correspondence) -----
#redefine sample grouping so that different extract types
#or different sample types
pd$sample_type[pd$extract_type == "ASW"] <- "ASW"
pd$sample_type[pd$extract_type == "DCM"] <- "DCM"
pd$sample_type[pd$extract_type == "deep"] <- "deep"
#change the phenodata frame
pData(data.peaks.ref.gsub.rt) <- pd

#make the peak density param object
pdp <- PeakDensityParam(sampleGroups = data.peaks.ref.gsub.rt$sample_type,
                        minFraction = 1/6, bw = 5)

#perform grouping
data.peaks.ref.gsub.rt.grp <- groupChromPeaks(data.peaks.ref.gsub.rt, 
                                              param = pdp)

#NB! can also load object after downloading from Github
load("data.peaks.ref.gsub.rt.grp.RData")
dirname(data.peaks.ref.gsub.rt.grp) <- rep(paste0(dir, "/mzML"), length(fp))

#9. Gap filling----
#define parameter object, use defaults
fpp <- FillChromPeaksParam()
#gap fill
data.peaks.ref.gsub.rt.grp.fill <- fillChromPeaks(data.peaks.ref.gsub.rt.grp,
                                                  param = fpp)

#NB! can also load object after downloading from Github
load("data.peaks.ref.gsub.rt.grp.fill.RData")
dirname(data.peaks.ref.gsub.rt.grp.fill) <- rep(paste0(dir, "/mzML"), 
                                                length(fp))

#rename final object to res
res <- data.peaks.ref.gsub.rt.grp.fill

#10. Check overlapping features------
over.feat <- overlappingFeatures(res) 
length(over.feat) == 0
#TRUE

#11. Summarize results-----
#summarise experiment before filling
res.exp <- quantify(data.peaks.ref.gsub.rt.grp, value = "into")
#feature definitions = rowData(res.exp)
#feature abundances = assay(res.exp)

#add filled peaks 
assays(res.exp)$raw_filled <- featureValues(res, filled = TRUE)

#12. Make feature table for NeatMS-----
#This section DIRECTLY follows code from the NeatMS documentation+
#Code source: https://github.com/bihealth/NeatMS/blob/master/docs/first-steps/data-format.md

#You can also directly load the table from Github after deisotoping!

# Create dataframe containing adjusted retention times
df_chrom_peaks_fill <- as.data.frame(chromPeaks(res))
# Create dataframe with raw retention times
feature_dataframe <- as.data.frame(chromPeaks(dropAdjustedRtime(res)))
# Get the peaks that have been recovered by the gapfilling step
df_filled_peaks <- df_chrom_peaks_fill[!row.names(df_chrom_peaks_fill)
                                        %in% c(rownames(feature_dataframe)),]
# Add them to the raw retention time dataframe
#(filled peaks do not have raw retention times so we use the adjusted ones)
feature_dataframe <- bind_rows(feature_dataframe, df_filled_peaks)

# Rename the retention time columns of the adjusted rt dataframe
df_chrom_peaks_fill <- df_chrom_peaks_fill %>%
   dplyr::rename(rt_adjusted = rt, rtmin_adjusted = rtmin, rtmax_adjusted = rtmax)

# Add the adjusted rt columns to the dataframe containing raw rt
# To have a dataframe containing both raw and adjusted rt information
feature_dataframe <- left_join(
  rownames_to_column(feature_dataframe),
  rownames_to_column(
    df_chrom_peaks_fill[,c("rt_adjusted","rtmin_adjusted","rtmax_adjusted")]),
  by="rowname")

# Remove the rownames as we won't need them
feature_dataframe$rowname <- NULL

# Retrieve the sample names and store them as a dataframe
sample_names_df <- as.data.frame(basename(fileNames(res)))

 # Rename the unique column "sample_name"
colnames(sample_names_df) <- c("sample_name")

# Generate the correct sample ids for matching purposes
# XCMS sampleNames() function returns sample names ordered by their ids
sample_names_df$sample <- seq.int(nrow(sample_names_df))

# Attach the sample names to the main dataframe by matching ids (sample column)
feature_dataframe <- left_join(feature_dataframe,sample_names_df, by="sample")

### Feature information addition ###

# Here we will bring the feature alignment information
#stored in the XCMSnExp object to the dataframe that we have already created
featuresDef <- featureDefinitions(res) 
featuresDef_df = data.frame(featuresDef)

# Adjust variable
# Only keep the information we need (column named 'peakidx')
# Get the index of the peakidx column
column_index <- which(colnames(featuresDef_df)=="peakidx")
# Extract the peakidx column
features_df = data.frame(featuresDef_df[,column_index])
# Rename the column
peak_colummn_name <- colnames(features_df)
features_df = dplyr::rename(features_df, "peak_id"=all_of(peak_colummn_name))
features_df <- cbind(feature_id= row.names(features_df),features_df)

# # We'll use data.table for the next step
# # Get all the peak_id for each feature_id
features_df <- data.table(features_df)
features_df = features_df[, list(peak_id = unlist(peak_id)), by=feature_id]

# Bring the feature_id to the original peak dataframe
feature_dataframe = cbind(peak_id= row.names(feature_dataframe),feature_dataframe)
feature_dataframe$peak_id = as.character(feature_dataframe$peak_id)
features_df$peak_id = as.character(features_df$peak_id)
feature_dataframe = left_join(feature_dataframe, features_df, by="peak_id")

# Note: The dataframe contains an extra column called peak_id,
#but this won't affect NeatMS and will simply be ignored (as would any other
#column not present in the list above).

#13: Deisotoping-----
#creat xcmsSet object
xset <- as(res, "xcmsSet")
#set sample names
sampnames(xset) <- pData(res)$name
#set sample classes/groups
sampclass(xset) <- pData(res)$sample_type
#build xsAnnotate object
an <- xsAnnotate(xset)
#group peaks by retention time (full width half maximum)
an <- groupFWHM(an, perfwhm = 0.6)
#find 13C isotopes 
an <- findIsotopes(an, ppm = 3)
#get feature list
pl <-getPeaklist(an)

#remove 13C isotopes
pl_deisotope <- pl %>% 
  #set rownames to column to keep ids consistent
  rownames_to_column() %>% 
  #make new column with isotope 'state'
  mutate(isotope_state = gsub('.*\\[M|\\]\\+$', '', isotopes) %>% 
           sub('^\\+(\\d).*', '\\1', .) %>% 
           sub('^\\].*', '', .)) %>% 
  #remove any 13C isotopes
  filter(isotope_state == '')
#get feature IDs of features retained after deisotoping
featureNamesDeiso <- featuresDef@rownames[as.numeric(pl_deisotope$rowname)]
#filter feature dataframe to keep only feature IDs without isotopes
featuresDef_df_deiso <- featuresDef_df[rownames(featuresDef_df) %in% 
                                         featureNamesDeiso,]

#filter NeatMS dataframe to match deisotoping output
feature_dataframe_deiso <- feature_dataframe %>% 
  drop_na(feature_id) %>% 
  filter(feature_id %in% as.numeric(pl_deisotope$rowname))

#check dimensions match before proceeding
dim(feature_dataframe_deiso)
#[1] 270061     17

#write to file for neatms
file_path <- "aligned_feature_table_deiso.csv"
write.csv(feature_dataframe_deiso,file_path, row.names = FALSE)

#aligned_feature_table_deiso.csv was used as the input for NeatMS denoising
#denoising run with the model on Github

#denoising output is also on GitHub! can download table there

#14: Import NeatMS output-----
ntms <- fread("neatms_export_with_extra_properties.csv") %>% 
  #rename NeatMS style names to XCMS style names
  dplyr::rename(maxo=height, into=area,intb=area_bc) %>% 
  #round (required for matching, digits changed during NeatMS processing)
  mutate(maxo = round(maxo, digits = 2),
         into = round(into, digits = 2),
         intb = round(intb, digits = 2))

#NOTE: NeatMS filters based on a minimum scan number, peaks lost due to this

#add phenodata
ntms_pd <- ntms %>% 
  mutate(name = gsub("MS31_20221004_|_500ng", "", sample) %>% 
           sub("_\\d{2}$", "", .)) %>% 
  left_join(pd, by = "name")

#set levels of names as sorted
pd.sort <- pd[order(pd$sample_type),]
ntms_pd$name <- factor(ntms.pd$name, levels = pd.sort$name)

#histogram of labels
ggplot(ntms_pd, 
       aes( x= name, fill = label)) +
  geom_histogram(stat = "count") +
  scale_fill_manual(values = viridis(n = 3, option = "E"),
                    name = "Label") +
  labs(x = "Sample name", y = "Count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#produces supplementary figure 14

#15: Merge NeatMS output with input-----
#'maxo' (height), 'into' (area), 'intb' (baseline corrected area)
#round digits on input table for matching 
feature_dataframe_n <- feature_dataframe_deiso %>% 
  drop_na(feature_id) %>% 
  mutate(name = gsub("MS31_20221004_|_500ng|.mzML", "", sample_name) %>% 
           sub("_\\d{2}$", "", .),
         maxo = round(maxo, digits = 2),
         into = round(into, digits = 2),
         intb = round(intb, digits = 2))

#join input and output tables and filter by retention time
ntms_fdf <- ntms_pd %>% 
  #join
  left_join(., feature_dataframe_n, by = c("name", "maxo", "into", "intb")) %>% 
  #add feature ID column
  mutate(featID = paste0("FT",str_pad(feature_id, 5, pad = "0"))) %>% 
  #filter retention time between 1 and 35 minutes
  filter(between(rt, 1*60, 35*60))

#keep only high and low quality peaks and extract samples
ntms_HQLQ_ex <- ntms_fdf %>% 
  filter(label != 'Noise' & sample_type %in% c('ASW', 'DCM', 'deep'))


#16: Look at high and low quality features in extracts----
#get filtered feature definition table
featA_def <- featuresDef_df_deiso[rownames(featuresDef_df_deiso) %in%
                                    ntms_HQLQ_ex$featID,] %>% 
  rownames_to_column('featID')
#should have 4,786 features

#get filtered feature intensity matrix (filled values)
featA_matrix <- assay(res.exp, "raw_filled") %>%
  as.data.frame() %>% 
  rownames_to_column(., "featID") %>%  
  filter(featID %in% ntms_HQLQ_ex$featID) %>% 
  column_to_rownames(., "featID") %>% 
  as.matrix()
#filter to only have extract samples
featA_matrix <- featA_matrix[,grepl("EN683|MSM107", colnames(featA_matrix))]
featA_matrix <- featA_matrix[rowSums(is.na(featA_matrix)) != ncol(featA_matrix), ]

#log2 transform values
featA_matrix_log2 <- featA_matrix
featA_matrix_log2[is.na(featA_matrix_log2)] <- 1
featA_matrix_log2 <- log2(featA_matrix_log2)

#make vectors of variables for plotting 
short_samples <- colnames(featA_matrix_log2) %>% gsub('.*4_|_\\d{2}.mzML', '', .)
colnames(featA_matrix_log2) <- short_samples
cruises <- str_split_i(short_samples, "-|_", 1)
cruise_shapes <- c(21, 24)
names(cruise_shapes) <- c("MSM107", "EN683")
shapes <- cruise_shapes[cruises]
sample_types <- str_split_i(short_samples, "-|_", 2) %>% 
  sub('\\d', '', .)

#perform PCA on log2-transformed matrix
pc <- prcomp(t(featA_matrix_log2), center = TRUE)
#make dataframe with output
pc_x <- pc$x %>% as.data.frame() %>% select(c("PC1", "PC2")) %>% 
  mutate(sample = short_samples,
         cruise = cruises,
         sample_type = sample_types)
#get summary of PCA
pcSummary <- summary(pc)

#plot using plotting function
plotPCA(pc_x, pcSummary)

#this should generate Fig. 2B (but with different colours)



#17: Annotate putative oligosaccharides-----
#make paramater object
param <- GlycoAnnotateR::glycoPredictParam(dp = c(2,6),
                                           modifications = c("carboxylicacid", 
                                                             "deoxy", "nacetyl", 
                                                             "oacetyl", 
                                                             "omethyl", 
                                                             "amino", 
                                                             "anhydrobridge", 
                                                             "alditol",
                                                             "sulfate"),
                                           pent_option = T, 
                                           polarity = 'pos',
                                           ion_type = 'ESI')
#generate predicted oligosaccharides
#for the same results the tool should be run with the version from Feb. 2024
pred_table <- glycoPredict(param)

#perform annotation
#detailed descriptions for this can be found at
#https://margotbligh.github.io/GlycoAnnotateR/
featA_def_annot <- glycoAnnotate(featA_def, 
                                 mz_column = 'mzmed', #name of mz column in data
                                 pred_table =  pred_table, #parameter object
                                 collapse = T, #collapse annotations
                                 collapse_columns = c('IUPAC name', 'ion'),
                                 error = 3) #mass deviation allowed
#keep only features with annotations = 338 features
featA_def_annotonly <- featA_def_annot %>%  drop_na(annotations)

#get matrix of those features intensities
featA_matrix_annotonly <- featA_matrix[rownames(featA_matrix) %in% 
                                                  featA_def_annotonly$featID,]

#18: Filter out features abundant in processing blanks (ASW)----
#Remove features that have >60% of the maximum abundance 
#in one of the ASW processing blanks
#first set NA values to 0
featA_matrix_asw <- featA_matrix_annotonly
featA_matrix_asw[is.na(featA_matrix_asw)] <- 0

featA_matrix_asw <- featA_matrix_asw %>% 
  as.data.frame() %>% 
  rownames_to_column('featID') %>% 
  rowwise() %>% 
  mutate(maxval = max(across(contains("mzML")))) %>% 
  filter(!if_any(contains("ASW"), ~ .x > 0.6*maxval )) %>% 
  select(-maxval) %>% 
  column_to_rownames('featID') %>% 
  as.matrix()

#193 features retained

#19: Remove potential outliers----
#transpose matrix and add sample information
#then pivot to long format
featA_matrix_asw_t <- t(featA_matrix_asw) %>% 
  as.data.frame() %>% 
  set_rownames(res$name[grepl("MSM107|EN683", res$name)]) %>% 
  rownames_to_column("sample") %>% 
  mutate(sample_type = res$sample_type[grepl("MSM107|EN683", res$name)],
         cruise = res$cruise[grepl("MSM107|EN683", res$name)]) %>% 
  pivot_longer(cols = starts_with("F"), names_to = "featID", values_to = "intensity") %>% 
  filter(intensity > 0)

#plot feature intensities across samples with log10 scale
ggplot(featA_matrix_asw_t, 
       aes(x=sample, y=log10(intensity), fill = sample_type)) + 
  geom_boxplot() +
  geom_jitter(shape = 21, alpha = 0.5) +
  geom_hline(yintercept = 8, colour = 'darkred',
             linewidth = 1.5) +
  labs(x = "Sample", y = "Log10 feature intensities") +
  theme_light()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#remove features with intensities greater than 1E8 as potential outliers
featA_matrix_out <- featA_matrix_asw %>% 
  as.data.frame() %>% 
  rownames_to_column('featID') %>% 
  rowwise() %>% 
  filter(!if_any(contains("mzML"), ~ .x > 1e8)) %>% 
  column_to_rownames('featID') %>% 
  as.matrix()

#189 features remaining

#20: Make df for  plotting of chromatograms for manual review-----
#Extracted ion chromatograms for final 192 features were manually reviewed
ntms_annotonly_asw_out <- ntms_HQLQ_ex %>% 
  filter(featID %in% rownames(featA_matrix_out))

#21: Import and plot chromatograms for manual checking------

#Extract and plot chromatograms
#Make empty directory to write plots to
dir.create('chromatogram_checks')

#Loop through peaks
#Will be quite slow
for(i in 1:nrow(ntms_annotonly_asw_out)){
  #Get info on peak
  i_mz = ntms_annotonly_asw_out$mz[i]
  i_rtmin = ntms_annotonly_asw_out$rtmin_adjusted[i]
  i_rtmax = ntms_annotonly_asw_out$rtmax_adjusted[i]
  i_samplenum = ntms_annotonly_asw_out$sample.y[i]

  #Subset data
  i_data <- filterFile(data.peaks.ref.gsub.rt.grp.fill,
                       i_samplenum)
  
  #Extract EIC +/- 5 minutes and +/- 3 ppm
  #Make sure that the minimum time is not below zero
  i_rtmin_eic = i_rtmin - 300 
  if(i_rtmin_eic < 0){i_rtmin_eic = 0}
  i_rtmax_eic = i_rtmin + 300
  
  i_eic <- chromatogram(i_data, 
                        mz = c(i_mz - ppm_to_mz(i_mz, 3),
                               i_mz + ppm_to_mz(i_mz, 3)),
                        rt = c(i_rtmin_eic, i_rtmax_eic))
  
  #Make plot object
  plot(i_eic)
  
  #Make filename for plot
  filename <- paste0("chromatogram_checks/",i, ".png")
  
  #Plot chromatogram to file
  png(filename = filename,
      width = 9, height = 3, units = "in", res = 300)
  plot(i_eic)
  dev.off()
  
}

#22: Filter features after manual checking------
#Results of the manual check can be downloaded from GitHUB 
peaks_manual_check <- fread('peaks_manually_checked.csv')

#Filter feature definition dataframe
featD_def <- featA_def_annotonly %>% 
  filter(featID %in% peaks_manual_check$featID)

#Filter intensity matrix
featD_matrix <- featA_matrix_out[rownames(featA_matrix_out) %in%
                                   featD_def$featID,]

#Make log2 matrix
featD_matrix_log2 <- log2(featD_matrix+1)

#23: PCA of oligosaccharide features----

#Perform pca
pc <- prcomp(t(featD_matrix_log2), center = TRUE)
pc_x <- pc$x %>% as.data.frame() %>% select(c("PC1", "PC2")) %>% 
  mutate(sample = short_samples,
         cruise = cruises,
         sample_type = sample_types)
pcSummary <- summary(pc)

#Plot PCA
plotPCA(pc_x, pcSummary)

#should produce Figure 2C



#24: K-means clustering of features-----
#First - Decide on the number of clusters

#Elbow method
n_clusters <- 10
wss <- numeric(n_clusters)
set.seed(123)
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(featD_matrix_log2, centers = i, nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)

#Supports 4 clusters
ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  geom_hline(yintercept = wss, linetype = 'dashed', 
             col = c(rep('grey',3),'darkred', rep('grey', 6))) +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  labs(x = 'Number of clusters', y = 'Within sum of squares') +
  theme_light()

#Silhouette method
avg_sil_values <-  numeric(n_clusters)
for (i in 1:n_clusters) {
  if(i == 1){next}
  # Fit the model: km.out
  km.out <- kmeans(featD_matrix_log2, centers = i, nstart = 20)
  # Save the within cluster sum of squares
  ss <- silhouette(km.out$cluster, dist(featD_matrix_log2))
  avg_sil_values[i] <-  mean(ss[, 3])
}

sil_df <- data.frame(clusters = 1:10, ss = avg_sil_values)

#Supports 2 or 4 clusters
ggplot(sil_df, aes(x = clusters, y = ss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  geom_vline(xintercept = 4, linetype = 'dashed', col = 'darkred') +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  labs(x = 'Number of clusters', y = 'Average silhouettes') +
  theme_light()

#Ideal number = 4
k <- 4
#Perform clustering on features
set.seed(123)
km.out <- kmeans(featD_matrix_log2, centers = k, nstart = 20)
#Add cluster info to feature matrix
data_km <- featD_matrix_log2 %>% 
  as.data.frame() %>% 
  rownames_to_column("featID")
#Make factor
data_km$clusterID <- factor(km.out$cluster)
#Make palette
pal_k <- c("#613A3A", "#F5B841", "#475B63", "#B185A7")

#Perform PCA on features
pc <- prcomp(featD_matrix_log2, center = TRUE)
#Get dataframe
pc_x <- pc$x %>% 
  as.data.frame() %>% 
  select(c("PC1", "PC2", "PC3"))
pc_x$clusterID <- data_km$clusterID 
#Summarise
pcSummary <- summary(pc)
#Plot with clusters on first two dimensions
ggscatter(pc_x, x = "PC1", y = "PC2", palette = pal_k,
  color = "clusterID", ellipse = TRUE, ellipse.type = "convex",
  shape = 20, size = 2,  legend = "right", ggtheme = theme_light(),
  stroke = 0.5, 
  xlab = paste0("PC1 (", format(pcSummary$importance[2, 1] * 100,
                                digits = 3), "% )" ),
  ylab = paste0("PC2 (", format(pcSummary$importance[2, 2] * 100,
                                  digits = 3), "% )" ))
#Plot with clusters on first and third dimensions
ggscatter(pc_x, x = "PC1", y = "PC3", palette = pal_k,
  color = "clusterID", ellipse = TRUE, ellipse.type = "convex",
  shape = 20, size = 2,  legend = "right", ggtheme = theme_light(),
  stroke = 0.5, 
  xlab = paste0("PC1 (", format(pcSummary$importance[2, 1] * 100,
                                digits = 3), "% )" ),
  ylab = paste0("PC2 (", format(pcSummary$importance[2, 3] * 100,
                                digits = 3), "% )" ))

#look in 3d
#not necessary but interesting, shows separation more clearly in 3 dimensions
library(plotly)
plot_ly(x=pc_x$PC1, y=pc_x$PC2, z=pc_x$PC3, type="scatter3d", mode="markers", 
        color=pc_x$clusterID)

#25: Look at mz distribution by cluster----
#Add clusters to feature definition table
featD_def_k <- featD_def %>% 
  left_join(data_km %>% select(featID, clusterID))

#Perform kruskal test
kruskal.test(mzmed ~ clusterID, data = featD_def_k)

#Kruskal-Wallis rank sum test
#data:  mzmed by clusterID
#Kruskal-Wallis chi-squared = 17.782, df = 3, p-value = 0.0004878

#Points to significant difference between clusters

#Pairwise comparisons using Wilcoxon rank sum exact test 
pairwise.wilcox.test(featD_def_k$mzmed, 
                     featD_def_k$clusterID,
                     p.adjust.method = "BH")

#Two comparisons have significant differences


#Plot violin plot with significant differences
my_comparisons <- list( c(1, 2), c(1, 3), c(2, 3), c(1,4), c(2,4), c(3,4))
names(pal_k) <- as.character(1:4)
ggplot(featD_def_k, aes(x = clusterID, y = mzmed, fill= clusterID)) +
  geom_violin() +
  scale_fill_manual(values = pal_k) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",
                     vjust = 0.25) +
  labs(x = "Cluster", y = expression(italic("m/z")))+
  theme_light()

#NOTE THAT CLUSTERS MAY NOT BE NUMBERED THE SAME AS IN THE PAPER
#Result is the same, just potentially with the numbers assigned
#to clusters different

#26: K-means clusters of samples ----
#Elbow method
n_clusters <- 10
wss <- numeric(n_clusters)
set.seed(123)
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(t(featD_matrix_log2), centers = i, 
                   nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)

#4 looks ideal
ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  geom_hline(yintercept = wss, linetype = 'dashed', 
             col = c(rep('grey',3),'darkred', rep('grey', 6))) +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  labs(x = 'Number of clusters', y = 'Within sum of squares') +
  theme_light()

#Silhouette method
avg_sil_values <-  numeric(n_clusters)
for (i in 1:n_clusters) {
  if(i == 1){next}
  # Fit the model: km.out
  km.out <- kmeans(t(featD_matrix_log2), 
                   centers = i, nstart = 20)
  # Save the within cluster sum of squares
  ss <- silhouette(km.out$cluster, dist(t(featD_matrix_log2)))
  avg_sil_values[i] <-  mean(ss[, 3])
}

sil_df <- data.frame(clusters = 1:10, ss = avg_sil_values)

#4 clusters
ggplot(sil_df, aes(x = clusters, y = ss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  geom_vline(xintercept = 4, linetype = 'dashed', col = 'darkred') +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  labs(x = 'Number of clusters', y = 'Average silhouettes') +
  theme_light()

#ideal number = 4
k <- 4
#Perform clustering
set.seed(123)
samples_km_out <- kmeans(t(featD_matrix_log2), 
                         centers = k, nstart = 20)
samples_data_km <- t(featD_matrix_log2) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample")
samples_data_km$clusterID <- factor(samples_km_out$cluster)
#Perform PCA
pc <- prcomp(t(featD_matrix_log2), center = TRUE)
pc_x <- pc$x %>% 
  as.data.frame() %>% 
  select(c("PC1", "PC2")) %>% 
  rownames_to_column('sample')
pc_x$clusterID <- samples_data_km$clusterID 
pcSummary <- summary(pc)
ggscatter(pc_x, x = "PC1", y = "PC2", 
  color = "clusterID", ellipse = TRUE, ellipse.type = "convex",
  shape = "clusterID", size = 3,  legend = "right", ggtheme = theme_light(),
  label = "sample",
  xlab = paste0("PC1 (", format(pcSummary$importance[2, 1] * 100,
                                digits = 3), "% )" ),
  ylab = paste0("Dim 2 (", format(pcSummary$importance[2, 2] * 100,
                                  digits = 3), "% )" ))

#27: Heatmap of features with sample and feature clusters----
#Set up annotation dataframes for colours
sample_colour = data.frame(sample_type = sample_types)
rownames(sample_colour) <- colnames(featD_matrix_log2)
sample_colours <- list(sample_type = c(deep = "#1D223C",
                                       DCM = "#DE7C1A",
                                       ASW = "#D1D1FE"),
                       cluster = c("1" = "#613A3A",
                                   "2" = "#F5B841",
                                   "3" = "#475B63",
                                   "4" = "#B185A7"))
row_annot_df <- data.frame(row.names = rownames(featD_matrix_log2), 
                           cluster = data_km$clusterID)
#Order data by cluster ID
data_km_ordered <- data_km %>% 
  arrange(clusterID)
data_km_samples_ordered <- samples_data_km %>% 
  arrange(clusterID)

#Order matrix by cluster ID  
m_ordered <- featD_matrix_log2[data_km_ordered$featID,
                               data_km_samples_ordered$sample]

#Add m/z as rownames
m_ordered_mz <- round(featD_def$mzmed[match(rownames(m_ordered), 
                                            featD_def$featID)], 3)

#Make short version of column names 
colnames_short <- str_split_i(colnames(m_ordered),'_', 4)

#Plot heatmap
#NOTE again that clusters may be numbered differently
pheatmap(mat = m_ordered,
         cellheight = 4,
         cellwidth = 10,
         cluster_rows = F, 
         cluster_cols = F,
         scale = "none",
         color = viridis(n = 1000),
         fontsize = 5,
         labels_row = sprintf('%.3f', m_ordered_mz),
         labels_col = colnames_short,
         annotation_col = sample_colour,
         annotation_row = row_annot_df,
         annotation_colors = sample_colours,
         annotation_legend = T, legend = T,
         angle_col = 45,
         family = "Arial",
         annotation_names_row = F,
         annotation_names_col = F)

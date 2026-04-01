#==============================================================================#
#==============================================================================#
######                                                                    ######
######            Clone maps to demonstrate preinvasive seeding           ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Katherine Honan
# Date: 2026-01-19

setwd("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/")

####################################################
#### Source required functions & load libraries ####
####################################################

library(fst)
library(readr)
library(tibble)
library(tidyr)
library(stringr)
library(data.table) 
library(dplyr) 
library(ggplot2) 
library(cowplot)
library(RColorBrewer) 
library(patchwork)

#############################################
#### Make a folder for this analysis run ####
#############################################

date <- gsub("-","",Sys.Date())

analysis_name <- 'preinvasive_seeding_clonemap'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep='/')
if( !file.exists(out_dir_general) ) dir.create( out_dir_general )

out_dir_logs <- paste(out_dir_general, 'logs', sep='/')
if( !file.exists(out_dir_logs) ) dir.create( out_dir_logs )

outputs.folder <- paste0( out_dir_general, "/", date, "/" )

if( !file.exists(outputs.folder) ) dir.create( outputs.folder )


##############################################
#### Get Inputs required for all analyses ####
##############################################

# read in data
trees_path <- "/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/inputs/_RELEASE/_aggregate/conipher_trees/"
sample_info <- fread("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/outputs/cohort_metadata/20250929/full_cohort_sample_level_data.tsv")

# source functions
source("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/personalis_ctDNA_mets_analysis/scripts/mets_functions.R")
source("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/personalis_ctDNA_mets_analysis/scripts/ctdna_mets_functions.R")
source("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/personalis_ctDNA_mets_analysis/scripts/useful_functions.R")
source( '/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/personalis_ctDNA_mets_analysis/scripts/fishplots/cloneMap.R' )


#############################################
#### Concatenate data from multiple runs ####
#############################################

trees_files <- list.files(trees_path, full.names = TRUE)
trees_list <- lapply(trees_files, readRDS)
trees <- unlist(trees_list, recursive = FALSE)
names(trees) <- vapply(trees, function(x) x$parameters$sampleID, character(1))

##############################################################################
### Get example case and plot clonemap for preinvasive vs invasive disease ###
##############################################################################

patient <- "EPA00020_1"
pat_tree <- trees[[patient]]$graph_pyclone$default_tree
pat_clonality <- clonalities_updated[[patient]]
pat_seeding <- cloneInfo.list[[patient]]

# get enough colours for whole tree
clone_names <- unique(c(pat_tree[,1], pat_tree[,2]))
n_col <- length(clone_names)
default_pallete <- brewer.pal(12, "Paired") 
clone_cols <- colorRampPalette(default_pallete)(n_col)
names(clone_cols) <- clone_names

# create clone map for each individual sample and save
for (i in 1:length(colnames(pat_clonality))){
  sample <- colnames(pat_clonality)[i]
  if (grepl("[0-9a-f]{6,}", sample)){
    sample_nohash <- sub("_[^_]+$", "", sample)
    sample_type <- sample_info[sample_info$sample_name_hash == sample,]$sample_type
    plot_title <- paste0(sample_nohash, " (", sample_type, ")")
    out_name <- paste0(sample_nohash, "_", sample_type, ".png")
  } else {
    plot_title <- sample
    out_name <- paste0(sample, ".png")
  }

  sample_ccfs <- as.data.frame(trees[[patient]]$nested_pyclone$ccf_cluster_table[, sample , drop = FALSE])
  sample_ccfs$clones <- rownames(sample_ccfs)
  sample_ccfs <- sample_ccfs %>%
    rename(CCF = sample) %>%
    mutate(CCF = CCF/100) %>%
    select(clones, CCF)
  
  delete <- c("4", "8")
  sample_ccfs[sample_ccfs$clones %in% delete,]$CCF <- 0
    
  png(
    filename = paste0(outputs.folder, out_name),
    width = 1200,
    height = 1200,
    res = 150
  )
  
  cloneMap(pat_tree, sample_ccfs, clone.cols = clone_cols, inc_parents_ccf = T)
  mtext(plot_title, side = 3, line = 1, font = 2)
  
  dev.off()
}







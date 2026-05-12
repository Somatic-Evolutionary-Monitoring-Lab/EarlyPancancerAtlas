#==============================================================================#
#==============================================================================#
######                                                                    ######
######                      Preinvasive-met seeding?                      ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Katherine Honan
# Date: 2026-05-01

setwd("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas")

#============================================#
# Source required functions & load libraries #
#============================================#

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

#=====================================#
# Make a folder for this analysis run #
#=====================================#

date <- gsub("-","",Sys.Date())

analysis_name <- 'preinv_to_met_seeding'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep='/')
if( !file.exists(out_dir_general) ) dir.create( out_dir_general )

out_dir_logs <- paste(out_dir_general, 'logs', sep='/')
if( !file.exists(out_dir_logs) ) dir.create( out_dir_logs )

outputs.folder <- paste0( out_dir_general, "/", date, "/" )

if( !file.exists(outputs.folder) ) dir.create( outputs.folder )


#======================================#
# Get Inputs required for all analyses #
#======================================#

# read in data on RDS
trees <- readRDS("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/_RELEASE/release-20260306/_aggregate/conipher_trees/2026_03_06_cohort_conipher_trees.RDS")
inventory <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.tsv")
metadata <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/map_files/20250514_chen_2017_ESCC_metadata_EPA_mapping.txt")
clone_info <- readRDS("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/outputs/preinvasive_seeding/20260501/20260501_cloneInfo.list.rds")

# source both original (tissue) mets functions and ctdna adapted functions
source("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/scripts/00.useful_functions/mets_functions.R")
source("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/scripts/00.useful_functions/useful_functions.R")
source("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/scripts/00.useful_functions/cloneMap.R")

#====================================================#
# Check if any mets are seeded by preinvasive clones #
#====================================================#

met_cases <- inventory %>%
  filter(sample_type_simple == "cancer_met") %>%
  pull(patient)

names(trees) <- sapply(trees, function(x) x$parameters$sampleID)
met_case_tumours <- names(trees)[gsub("*_.", "", names(trees)) %in% met_cases]

# clonality tables
met_case_clonalities <- lapply(met_case_tumours, function(i){
  clonality <- trees[[i]]$clonality_out$clonality_table_corrected
})
names(met_case_clonalities) <- met_case_tumours

# investigate metastatic seeding
metastatic_seeding <- lapply(met_case_tumours, function(i){
  # get clonality table
  clonality_pat <- met_case_clonalities[[i]]
  # list all tree regions
  regions <- colnames(clonality_pat)
  # identify sample types
  cancer_regions <- as.character(
    inventory %>% 
      filter(
        sample_name_hash %in% regions, sample_type_simple == "cancer") %>% 
      pull(sample_name_hash))
  preinv_regions <- as.character(
    inventory %>% 
      filter(
        sample_name_hash %in% regions, sample_type_simple == "preinvasive") %>% 
      pull(sample_name_hash))
  met_regions <- as.character(
    inventory %>% 
      filter(
        sample_name_hash %in% regions, sample_type_simple == "cancer_met") %>% 
      pull(sample_name_hash))
  # dataframe of sample types for plotting
  n_samples_df <- data.frame(
    patient = i, 
    n_cancer_regions = length(cancer_regions), 
    n_preinv_regions = length(preinv_regions),
    n_met_regions = length(met_regions)
  )
  
  # Get sample type level clonality
  clonality_pat$primary_clonality <- get.tumLevel.clonality(clonality_pat, cancer_regions)
  clonality_pat$preinvasive_clonality <- get.tumLevel.clonality(clonality_pat, preinv_regions)
  clonality_pat$met_clonality <- get.tumLevel.clonality(clonality_pat, met_regions)
  
  # do any clones appear in prein and met but not primary?
  preinvasive_seeding_clones <- rownames(clonality_pat[
    clonality_pat$met_clonality != "absent" & clonality_pat$preinvasive_clonality != "absent" & clonality_pat$primary_clonality == "absent",
      ]
  )
  
  # tree
  tree_structure <- get.treeStructure(trees[[i]]$graph_pyclone$default_tree, trees[[i]]$graph_pyclone$trunk)
  
  return(list(sample_info = n_samples_df, preinvasive_seeding_clones = preinvasive_seeding_clones, tree_structure = tree_structure))
})

names(metastatic_seeding) <- met_case_tumours




#============================#
# Plot per sample clone maps #
#============================#
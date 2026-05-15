#==============================================================================#
#==============================================================================#
######                                                                    ######
######                         Mutation timing                            ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Katherine Honan
# Date: 2025-05-30

setwd("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas")

#============================================#
# Source required functions & load libraries #
#============================================#

library(readr)
library(tidyr)
library(dplyr) 
library(ggplot2) 
library(RColorBrewer) 
library(cowplot)
library(data.table)
library(fst)

#=====================================#
# Make a folder for this analysis run #
#=====================================#

date <- gsub("-","",Sys.Date())

analysis_name <- 'mutation_timing'
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

# read in full cohort data
trees <- readRDS("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/_RELEASE/release-20260306/_aggregate/conipher_trees/2026_03_06_cohort_conipher_trees.RDS")
mutations <- read.fst("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/_RELEASE/release-20260306/_aggregate/mutation_table/2026_03_06_cohort_concise_mutation_table.fst")
clone_info <- readRDS("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/outputs/preinvasive_seeding/20260501/20260501_cloneInfo.list.rds")

#===============#
# Preprocessing #
#===============#

names(trees) <- sapply(trees, function(x) x$parameters$sampleID)

#========================================#
# Classify clones into timing categories #
#========================================#

# classify clones into timing categories and add number of mutations per clone
clone_timing_categories <- lapply(names(clone_info), function(pat){
  # MRCA, shared, seeding, preinvasive specific and primary specific
  clone_timing <- clone_info[[pat]]$df.cloneInfo %>%
    mutate(
      timing_category = case_when(
        clonalClones & !seedingClones ~ "MRCA Non-initiating",
        seedingClones & clonalClones ~ "MRCA Initiating",
        seedingClones & !clonalClones ~ "Subclonal initiating",
        sharedClones & !seedingClones ~ "Shared",
        preinvasiveClones ~ "Preinvasive private",
        tumourClones ~ "Tumour private",
        TRUE ~ NA_character_
      )
    ) %>%
    select(
      patient_tumour = Patient,
      clone = clones,
      timing_category
    ) %>%
    filter(!is.na(timing_category)) # remove clones from met samples
  
  # get edge lengths
  edgelengths_pat <- trees[[pat]]$graph_pyclone$edgelength
  
  # add num mutations per clone
  clone_timing <- clone_timing %>%
    mutate(n_mutations = as.numeric(edgelengths_pat[as.character(clone)]),
           prop_muts = n_mutations/sum(n_mutations)) 
  
})
clone_timing_categories_df <- bind_rows(clone_timing_categories)

write_delim(clone_timing_categories_df, paste0(outputs.folder, date, "_mutation_timing_per_clone.tsv"), delim = "\t")

#======================================#
# Count mutations in timing categories #
#======================================#

clone_timing_categories_counts <- clone_timing_categories_df %>%
  group_by(patient_tumour, timing_category) %>%
  summarise(
    num_mutations = sum(n_mutations),
    prop_muts = sum(prop_muts),
    .groups = "drop"
)
  
write_delim(clone_timing_categories_counts, paste0(outputs.folder, date, "_mutation_count_per_timing.tsv"), delim = "\t")


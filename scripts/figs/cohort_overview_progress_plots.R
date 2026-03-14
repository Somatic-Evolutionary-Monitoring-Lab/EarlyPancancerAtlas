#!/usr/bin/env Rscript
# Author: Katherine Honan
# Date: 2026-03-12
# Description: EPA data acquisition progress plots

#================#
# Load libraries #
#================#

library(readr)
library(dplyr)
library(ggplot2)

options(scipen=999)

#=========================#
# Create output directory #
#=========================#

setwd("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/")

date <- gsub("-","",Sys.Date())

analysis_name <- 'cohort_overview_plots'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep='/')
if( !file.exists(out_dir_general) ) dir.create( out_dir_general )

out_dir_logs <- paste(out_dir_general, 'logs', sep='/')
if( !file.exists(out_dir_logs) ) dir.create( out_dir_logs )

outputs.folder <- paste0( out_dir_general, "/", date, "/" )

if( !file.exists(outputs.folder) ) dir.create( outputs.folder )


#===========#
# Load Data #
#===========#

# The inventory contains info needed to run alignment pipeline + sample type info, 
# this is iteratively added to as more datasets are being downloaded
inventory <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.tsv")

# There is one mapping file per study/ dataset, it exists to link back the 
# assigned EPA project IDs with the original metadata obtained from the repository 
# where the data was downloaded
map_pattern <- "*metadata_EPA_mapping.txt"



#=========================================#
# Generate plots to show dataset overview #
#=========================================#




#=====================================================#
# Generate cohort data quality plots for aligned data #
#=====================================================#











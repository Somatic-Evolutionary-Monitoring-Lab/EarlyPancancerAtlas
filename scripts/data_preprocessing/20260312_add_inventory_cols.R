#!/usr/bin/env Rscript
# Author: Katherine Honan
# Date: 2026-03-12
# Description: add sample type descriptor to inventory file for previous datasets

#================#
# Load libraries #
#================#

library(readr)
library(dplyr)
library(ggplot2)

#===========#
# Load Data #
#===========#

inventory_curr <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.tsv")
htan_lung_map <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/map_files/20260123_HTAN_Lung_metadata_EPA_mapping.tsv")
chen_escc_map <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/map_files/20250514_chen_2017_ESCC_metadata_EPA_mapping.txt")
black_eac_map <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/map_files/20250609_E_Black_EAC_WES_metadata_EPA_mapping.txt")

#==============================#
# add sample_type to inventory #
#==============================#

# Chen ESCC was the first data to be downloaded and for this dataset I renamed the
# raw fastqs but decided for later downloads I would just keep a mapping file per 
# study to keep track of fq -> EPA sample name has mapping. Handling this slightly 
# differently to others, which have fq mapping info in their map files.

inventory_curr$sample_type_detailed <- chen_escc_map$descriptor[
  match(inventory_curr$sample_name_hash, chen_escc_map$EPA)
]

# fill NAs from Black EAC
idx <- is.na(inventory_curr$sample_type_detailed)
# match on smaple_name_hash
inventory_curr$sample_type_detailed[idx] <- black_eac_map$specimen_phenotype[
  match(inventory_curr$sample_name_hash[idx], black_eac_map$sample_name_hash)
]

# fill remaining NAs from HTAN Lung
idx <- is.na(inventory_curr$sample_type_detailed)
# match on smaple_name_hash
inventory_curr$sample_type_detailed[idx] <- htan_lung_map$sample_type[
  match(inventory_curr$sample_name_hash[idx], htan_lung_map$sample_name_hash)
]

# create simplified sample_type column with consistant naming

sample_map <- c(
  germline = "germline",
  normal = "germline",
  normal_tissue = "germline",
  dysplasia_low = "preinvasive",
  dysplasia_high = "preinvasive",
  dysplasia = "preinvasive",
  barretts = "preinvasive",
  IM = "preinvasive",
  GM = "preinvasive",
  GM_IM = "preinvasive",
  IMC = "preinvasive",
  DYS_GM = "preinvasive",
  DYS_IMC = "preinvasive",
  preinvasive = "preinvasive",
  carcinoma = "cancer",
  cancer = "cancer",
  oac = "cancer",
  metastasis_lymphnode = "cancer_met"
)

inventory_curr$sample_type_simple <- sample_map[inventory_curr$sample_type_detailed]

# Write out modified inventory file
write_tsv(inventory_curr, "/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.tsv")





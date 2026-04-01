#!/usr/bin/env Rscript
# Author: Katherine Honan
# Date: 2026-03-13
# Description: Annotate phs001424 GBM project data
# Project annotation number: 6

# We want to create two files with this script:
# 1. map file: biospecimen/sample level info  with one row per biospecimen which has study specific information
# along with two columns for the fastq files (r1, r2) and the sample type
# 2. dataframe to append to the central inventory file used for alignment. this
# should include one row per fastq file (sometimes biosamples split over multiple lanes etc)

#================#
# Load libraries # 
#================#

library(dplyr)
library(readr)
library(stringr)
library(openssl)

options(scipen=999)
set.seed(123)

#=========================#
# Create output directory #
#=========================#

setwd("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/")

date <- gsub("-","",Sys.Date())

analysis_name <- 'mapping_and_inventory_files/phs001424_Brain_Brastianos'
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

biospecimen_metadata <- read_csv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/study_metadata/phs001424_Brain_SRA_Run_metadata.csv")
latest_inventory <- read_tsv('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.tsv')
# note: flow cell and lane could not be reconstructed from fastq files so will need to amend
fastq_seq_dat <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/study_metadata/phs001424_Brain_Brastianos_seq_platform_info.tsv")
fq_fail_list <- read_delim("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/fq_fail_list.txt")

#=================#
# Create map file #
#=================#

map_file <- biospecimen_metadata %>%
  mutate(germline = ifelse(Is_Tumor == "No", TRUE, FALSE),
         sample_type_simple = ifelse(germline, "normal", "cancer"),
         sample_type_detailed = case_when(
           grepl("-P", `Sample Name`) ~ "pre-treatment gbm",
           grepl("-MT", `Sample Name`) ~ "post-treatment gbm autopsy",
           TRUE ~ "normal"
         )
  ) %>%
  group_by(submitted_subject_id) %>%
  mutate(
    # assign sequential identifiers for germline samples
    GL_counter = as.integer(dense_rank(ifelse(germline, `Sample Name`, NA))),
    # assign independent counters for tumour samples
    S_counter = as.integer(dense_rank(ifelse(!germline, `Sample Name`, NA)))
  ) %>%
  # add sample ids
  mutate(sample_num = ifelse(germline, paste0("GL", GL_counter), paste0("S", S_counter))) %>%
  ungroup() %>%
  select(-GL_counter, -S_counter)

# generate unique EPA IDs per biosample
#curr_max_epa <- max(as.integer(sub("^EPA(\\d+).*", "\\1", latest_inventory$sample_name)))
curr_max_epa <- 228

map_file <- map_file %>%
  mutate(patient = sprintf("EPA%05d", dense_rank(submitted_subject_id) + curr_max_epa)) %>%
  # get random unique sample_hash for each sample
  mutate(sample_hash = sapply(1:n(), function(x) substr(as.character(md5(as.character(x))), 1, 12))) %>%
  # combine EPA patient number, sample_num, and sample_hash into the final ID
  mutate(sample_name_hash = paste0(patient, "_", sample_num, "_", sample_hash),
         sample_name = paste0(patient, "_", sample_num))

# check if there are any samples with no germline
participants_no_germline <- map_file %>%
  group_by(submitted_subject_id) %>%
  summarise(has_germline = any(germline == TRUE, na.rm = TRUE)) %>%
  filter(!has_germline) # none

# broadcast germline sample info to all samples
map_file <- map_file %>% 
  mutate(germline_use = germline & grepl("GL1", sample_name_hash)) %>%
  group_by(submitted_subject_id) %>%
  mutate(
    gl_sample_name_hash = if_else(
      !germline,
      sample_name_hash[germline_use][1],
      NA_character_
    )
  ) %>%
  ungroup()

# Add clinical info and some other metadata needed for inventory
map_file_clinical <- map_file %>%
  mutate(
    sample_class = ifelse(germline, "gl", "ffpe"),
    run_dir = "phs001424_Brain_Brastianos",
    sequencer = "HiSeq2500",
    ffpe = ifelse(sample_class == "ffpe", TRUE, FALSE),
    seq_id = `Sample Name`,
    seq_type = "WES",
    date_sequenced = format(as.Date("2017-10-12"), "%d/%m/%Y"),
    clinical_sex = tolower(sex),
    surgical_sample = ifelse(body_site == "Blood", FALSE, TRUE)
  ) 


# write out map file
write_delim(map_file_clinical, paste0(outputs.folder, date, "_phs001424_Brain_Brastianos_metadata_epa_mapping.tsv"), delim = "\t")


#=======================#
# Add data to inventory #
#=======================#

# first, write out the inventory unchanged but with .bak so that the last safe
# copy is preserved in case something needs correcting
write_delim(latest_inventory, paste0(outputs.folder, "epa-project.inventory.latest.pass.bak.tsv"), delim = "\t")

# remove excluded fastqs and column cleaning
fastq_seq_dat_alt <- fastq_seq_dat %>%
  mutate(Run = str_split_i(sample, "_", 1),
         lane = "NA") %>%
  left_join(map_file_clinical, by = "Run") %>%
  mutate(fq1 = basename(r1),
         fq2 = basename(r2),
         passing_lanes = 1,
         flowcell = gsub("\\\\,\\s*", "_", `flowcell_barcode (run)`)
         )

# check that the number of fastqs corresponds 1:1 with the number of unique 
# flowcell/lane combinations
lane_count_sanity_check <- fastq_seq_dat_alt %>%
  group_by(sample_name_hash) %>%
  summarise(
    rows = n(),
    unique_lane_fc = n_distinct(flowcell, lane),
    flowcells = n_distinct(flowcell),
    .groups = "drop"
  ) %>%
  arrange(desc(unique_lane_fc))

# any unexpected?
sample <- lane_count_sanity_check[!lane_count_sanity_check$rows == lane_count_sanity_check$unique_lane_fc,]$sample_name_hash

# Generate new inventory
run_loc <- '/rds/project/rds-LH0AvU65IRI/datasets/EPA/Brain/fastqs/'
inventory_new <- fastq_seq_dat_alt %>%
  # add any remaining columns needed
  mutate(
    fq_pass = ifelse(fq1 %in% fq_fail_list$failed_fastq | fq2 %in% fq_fail_list$failed_fastq, FALSE, TRUE),
    uid = paste0(run_dir, "_", seq_id, "_", flowcell, "_", lane),
    fq1 = paste0(run_loc, fq1),
    fq2 = paste0(run_loc, fq2)
  ) %>%
  select(colnames(latest_inventory))

# write out inventory for this study
write_delim(inventory_new, paste0(outputs.folder, date, "_epa-project.inventory.latest.pass.tsv"), delim = "\t")

# append to previous inventory
inventory_updated <- rbind(latest_inventory, inventory_new)

# write out updated inventory
write_delim(inventory_updated, "/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.tsv", delim = "\t")

#=======#
# Notes #
#=======#

## Note 1
#' This data was downloaded in bam format and then converted to fastqs using samtools
#' based on the read group information stored in the original bams, e.g., 
#' samtools view -H SRR6069480.bam | grep '^@RG'
#' @RG	ID:C5E28.1	PL:illumina	PU:C5E28ANXX150504.1.AGGTGCGA-GACATTAA	LB:Pond-400424	PI:0	DT:2015-05-05T00:00:00-0400	SM:GS-006-P	CN:BI
#' @RG	ID:C5KNA.3	PL:illumina	PU:C5KNAANXX150504.3.AGGTGCGA-GACATTAA	LB:Pond-400424	PI:0	DT:2015-05-05T00:00:00-0400	SM:GS-006-P	CN:BI
#' @RG	ID:C5KNA.4	PL:illumina	PU:C5KNAANXX150504.4.AGGTGCGA-GACATTAA	LB:Pond-400424	PI:0	DT:2015-05-05T00:00:00-0400	SM:GS-006-P	CN:BI
#' @RG	ID:C5KNA.5	PL:illumina	PU:C5KNAANXX150504.5.AGGTGCGA-GACATTAA	LB:Pond-400424	PI:0	DT:2015-05-05T00:00:00-0400	SM:GS-006-P	CN:BI
#' The samples were sequenced over many lanes and flowcells, but since we have converted
#' back to fastqs, we only have one per sample. Therefore have set passing_lanes = 1






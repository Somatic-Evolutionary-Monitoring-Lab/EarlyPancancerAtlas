#!/usr/bin/env Rscript
# Author: Katherine Honan
# Date: 2026-03-13
# Description: Annotate HTAN Stanford Colon Project Samples

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

analysis_name <- 'mapping_and_inventory_files/Stanford_Colon'
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

biosamples <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/data_downloads/HTAN/Stanford_colon/HTAN_Stanford_colon_Precancer_Biosample_Metadata.tsv")
fastq_metadata <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/data_downloads/HTAN/Stanford_colon/HTAN_Stanford_colon_Precancer_fq_Metadata.tsv")
fastq_seq_dat <- read_tsv('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/data_downloads/HTAN/Stanford_colon/HTAN_Stanford_colon_Precancer_seq_platform_info.tsv', col_names = T)
cases <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/data_downloads/HTAN/Stanford_colon/HTAN_Stanford_colon_Precancer_case_level_data.tsv")
latest_inventory <- read_tsv('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.tsv')
map_file_example <- read_tsv('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/map_files/20260123_HTAN_Lung_metadata_EPA_mapping.tsv')

## Process metadata files
# Keep relevent columns
cols_to_keep <- c("Biospecimen", "HTAN Parent ID", "Timepoint Label", "Biospecimen Type",
                  "Acquisition Method Type", "Participant ID", "Tumor Tissue Type")

clinical_cols_to_keep <- c("HTAN Participant ID", "Primary Diagnosis", "Ethnicity", 
                           "Gender", "Race")              

biosample_metadata_alt <- biosamples %>%
  select_if(~!all(is.na(.))) %>%
  rename(Biospecimen = `HTAN Biospecimen ID`) %>%
  select(all_of(cols_to_keep))

clinical_dat_alt <- cases %>%
  select_if(~!all(is.na(.))) %>%
  select(all_of(clinical_cols_to_keep)) %>%
  rename(`Participant ID` = `HTAN Participant ID`)


#=================#
# Create map file #
#=================#

# Notes:
# This data has a mix of WES and WGS. Some biosamples have one or the other, others
# have both. Annoyingly, for the germline and carcinoma samples there is only WGS.
# As a result, I am going to first assign each biosample a unique sample number/
# EPA ID, but then when it comes to making the inventory file I will exclude the 
# wgs fastqs in cases where a biosample has both wgs and wes.
# If the germline/ wgs only samples have decent coverage over the exome I will
# subset them to exome and then use these in the pipeline. If not, the dataset
# might not be useable - kh723 20260313

map_file <- biosample_metadata_alt %>%
  mutate(germline = ifelse(`Acquisition Method Type` == "Blood draw", TRUE, FALSE),
         sample_type = case_when(
           germline ~ "normal",
           `Tumor Tissue Type` == "Primary" ~ "cancer",
           `Tumor Tissue Type` == "Premalignant" ~ "preinvasive",
           TRUE ~ NA_character_
         )
  )  %>%
  group_by(`Participant ID`) %>%
  mutate(
    # assign sequential identifiers for germline samples
    GL_counter = as.integer(dense_rank(ifelse(germline, `Biospecimen`, NA))),
    # assign independent counters for tumour samples
    S_counter = as.integer(dense_rank(ifelse(!germline, `Biospecimen`, NA)))
  ) %>%
  # add sample ids
  mutate(sample_num = ifelse(germline, paste0("GL", GL_counter), paste0("S", S_counter))) %>%
  ungroup() %>%
  select(-GL_counter, -S_counter)

# generate unique EPA IDs per biosample
curr_max_epa <- max(as.integer(sub("^EPA(\\d+).*", "\\1", latest_inventory$sample_name)))

map_file <- map_file %>%
  mutate(patient = sprintf("EPA%05d", dense_rank(`Participant ID`) + curr_max_epa)) %>%
  # get random unique sample_hash for each sample
  mutate(sample_hash = sapply(1:n(), function(x) substr(as.character(md5(as.character(x))), 1, 12))) %>%
  # combine EPA patient number, sample_num, and sample_hash into the final ID
  mutate(sample_name_hash = paste0(patient, "_", sample_num, "_", sample_hash),
         sample_name = paste0(patient, "_", sample_num))

# check if there are any samples with no germline
participants_no_germline <- map_file %>%
  group_by(`Participant ID`) %>%
  summarise(has_germline = any(germline == TRUE, na.rm = TRUE)) %>%
  filter(!has_germline) # none

# broadcast germline sample info to all samples
map_file <- map_file %>% 
  mutate(germline_use = germline & grepl("GL1", sample_name_hash)) %>%
  group_by(`Participant ID`) %>%
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
  left_join(clinical_dat_alt, by = "Participant ID") %>%
  mutate(
    sample_class = ifelse(germline, "gl", "ffpe"),
    run_dir = "HTAN_Stanford_Colon",
    sequencer = "Novaseq",
    ffpe = ifelse(sample_class == "ffpe", TRUE, FALSE),
    seq_id = Biospecimen,
    date_sequenced = format(as.Date("2022-01-01"), "%d/%m/%Y"),
    clinical_sex = tolower(Gender),
    surgical_sample = ifelse(germline, FALSE, TRUE)
  ) %>%
  # remove redundant gender col
  select(-c("Gender", "Biospecimen"))

# write out map file
write_delim(map_file_clinical, paste0(outputs.folder, date, "_HTAN_Stanford_Colon_metadata_epa_mapping.tsv"), delim = "\t")

#=======================#
# Add data to inventory #
#=======================#

# first, write out the inventory unchanged but with .bak so that the last safe
# copy is preserved in case something needs correcting
write_delim(latest_inventory, "/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.bak.tsv", delim = "\t")

# start with the fastq seq info data, add the fastq metadata info required, filter
# out those we dont want, then join sample info from the map_file_clinical df
fastq_metadata_alt <- fastq_metadata %>%
  select_if(~!all(is.na(.))) %>%
  mutate(
    sample = basename(Filename),
    seq_type = ifelse(grepl("wgs", Filename), "WGS", "WES")
    ) %>%
  group_by(Biospecimen) %>%
  # drop wgs rows if biospecimen also has wes
  filter(!(any(seq_type == "WES") & seq_type == "WGS")) %>%
  ungroup()

# add info to the sequencing info file 
fastq_seq_dat_alt <- fastq_seq_dat %>%
  left_join(fastq_metadata_alt, by = "sample")

# get a list of wgs files that wont be aligned
wgs_excluded <- fastq_seq_dat_alt %>%
  filter(is.na(Biospecimen)) %>%
  select_if(~!all(is.na(.)))

#write_delim(
#  wgs_excluded, 
#  "/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/data_downloads/HTAN/Stanford_colon/wgs_exclude_list.tsv",
#  delim = "\t"
#  )

# remove excluded fastqs and column cleaning
fastq_seq_dat_alt <- fastq_seq_dat_alt %>%
  filter(!is.na(Biospecimen)) %>%
  rename("seq_id" = "Biospecimen") %>%
  mutate(fq1 = basename(r1),
         fq2 = basename(r2)) %>%
  left_join(map_file_clinical, by = "seq_id") %>%
  group_by(sample_name_hash) %>%
  mutate(passing_lanes = n_distinct(flowcell, lane)) %>%
  ungroup()

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

# only one which I think is just a fastq file that has been split, will need 
# pipeline to re-merge this after alignment so modify passing_lanes to 2
fastq_seq_dat_alt[fastq_seq_dat_alt$sample_name_hash == sample, ]$passing_lanes <- 2
# also change the lane number so uid is correctly unique 
fastq_seq_dat_alt[fastq_seq_dat_alt$sample_name_hash == sample, ]$lane <- c(1, 2)

run_loc <- '/rds/project/rds-LH0AvU65IRI/datasets/EPA/Colon/fastqs/'
inventory_new <- fastq_seq_dat_alt %>%
  # add any remaining columns needed
  mutate(
    fq_pass = TRUE,
    uid = paste0(run_dir, "_", seq_id, "_", flowcell, "_", lane),
    sample_type_detailed = sample_type,
    sample_type_simple = sample_type,
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






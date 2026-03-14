#!/usr/bin/env Rscript
# Author: Katherine Honan
# Date: 2026-03-14
# Description: Annotate HTAN Pilot Preinvasive Pancreatic dataset

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

analysis_name <- 'mapping_and_inventory_files/HTAN_Pilot_Pancreas'
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

biospecimen_metadata <- read_csv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/study_metadata/HTAN_Pilot_Precancer_Pancreatic_SRA_metadata.csv")
fastq_seq_dat <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/study_metadata/HTAN_Pilot_Precancer_Pancreatic_seq_platform_info.tsv")
latest_inventory <- read_tsv('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.tsv')

#=================#
# Create map file #
#=================#

# These samples were described as 'normal pancreas' on the SRA website 
normal_tissue <- c(
  "SRR25512619",
  "SRR25512627",
  "SRR25512690",
  "SRR25512744",
  "SRR25512789",
  "SRR25512814",
  "SRR25512844",
  "SRR25512859",
  "SRR25512874",
  "SRR25512912",
  "SRR25512905",
  "SRR25512914",
  "SRR25512939",
  "SRR25512938",
  "SRR25512991",
  "SRR25512964",
  "SRR25512995"
)

# These samples were described as PDAC samples on the SRA website 
pdac_samples <- c(
  "SRR25512801",
  "SRR25512822")

# This sample was described as 'whole slide ffpe' on the SRA website- assume premalignant due to metadata labelling
whole_slide <- c(
  "SRR25512675"
  )

# Described only as "pancreas acinar tissue"- assume premalignant due to metadata labelling
acinar_tissue <- c(
  "SRR25512628", 
  "SRR25512832", 
  "SRR25512840",
  "SRR25512857", 
  "SRR25512873", 
  "SRR25512913",
  "SRR25512906", 
  "SRR25512924", 
  "SRR25512959", 
  "SRR25512987", 
  "SRR25512965", 
  "SRR25513008",
  "SRR25513018"
  )

# These samples were described as germline controls on the SRA website 
germline_controls <- c(
  "SRR25512642", 
  "SRR25512661", 
  "SRR25512753", 
  "SRR25512779", 
  "SRR25512792", 
  "SRR25512811",
  "SRR25512828", 
  "SRR25512831", # v small fastq size, think will fail on coverage
  "SRR25512838", 
  "SRR25512883", 
  "SRR25512894",
  "SRR25512926", 
  "SRR25512941", 
  "SRR25512952", 
  "SRR25512970", 
  "SRR25512982", # descirbed as 'targeted seq' in description but WXS in metadata, patient has other gl samples 
  "SRR25513012",
  "SRR25513025", 
  "SRR25513043", 
  "SRR25513055"
  )

# all other samples in the SRA repo aside from those specified above were described
# explicitly as premalignant, 

# assign germline samples
biospecimen_metadata$germline <- ifelse(biospecimen_metadata$Run %in% c(germline_controls, normal_tissue), TRUE, FALSE)

# Not all patients have at least one germline sample
no_germline_patients <- biospecimen_metadata %>%
  group_by(submitted_subject_id) %>%
  summarise(has_germline = any(germline == TRUE, na.rm = TRUE), .groups = "drop") %>%
  filter(!has_germline) %>%
  pull(submitted_subject_id)

# list of samples to remove
samples_to_exclude <- biospecimen_metadata %>%
  filter(submitted_subject_id %in% no_germline_patients) %>%
  pull(Run)

# filter out all the samples for patients without any germline
biospecimen_metadata_alt <- biospecimen_metadata %>%
  filter(!Run %in% samples_to_exclude)

# Also eremove from fastq data
fastq_seq_dat_alt <- fastq_seq_dat %>%
  mutate(
    fq1 = basename(r1),
    fq2 = basename(r2),
    Run = str_split_i(sample, "_", 1)
    ) %>%
  filter(!Run %in% samples_to_exclude)

# write out list of SRR Run IDs for the fastqs to move to cold storage
# writeLines(samples_to_exclude, "/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/datasets/EPA/Pancreas/fastqs/samples_to_exclude.txt", sep = "\n")

# Now create initial map file
map_file <- biospecimen_metadata_alt %>%
  mutate(sample_type = case_when(
           germline ~ "normal",
           Run %in% pdac_samples ~ "cancer",
           Run %in% c(whole_slide, acinar_tissue) ~ "preinvasive",
           TRUE ~ "preinvasive"
         )
  )

# Assign which germline to use 
# in cases where there are multiple, first prioritise those that are in 'germline_controls'
# over those in 'normal_samples' and take the largest
map_file <- map_file %>%
  group_by(submitted_subject_id) %>%
  mutate(
    germline_use =
      germline &
      row_number() == which.max(
        ifelse(
          germline,
          Bytes + ifelse(Run %in% germline_controls, 1e15, 0),
          -Inf
        )
      )
  ) %>%
  ungroup()

# Assign sequential sample IDs within each participant
map_file <- map_file %>%
  # make sure the GL_counter assigns the gl_use as GL1
  arrange(submitted_subject_id, desc(germline_use)) %>%
  group_by(submitted_subject_id) %>%
  mutate(
    # assign sequential identifiers for germline samples
    GL_counter = as.integer(row_number(ifelse(germline, submitted_subject_id, NA))),
    # assign independent counters for tumour samples
    S_counter = as.integer(row_number(ifelse(!germline, submitted_subject_id, NA)))
  ) %>%
  # add sample ids
  mutate(sample_num = ifelse(germline, paste0("GL", GL_counter), paste0("S", S_counter))) %>%
  ungroup() %>%
  select(-GL_counter, -S_counter)

# generate unique EPA IDs per biosample
curr_max_epa <- max(as.integer(sub("^EPA(\\d+).*", "\\1", latest_inventory$sample_name)))

# generate EPA IDs
map_file <- map_file %>%
  mutate(patient = sprintf("EPA%05d", dense_rank(submitted_subject_id) + curr_max_epa)) %>%
  # get random unique sample_hash for each sample
  mutate(sample_hash = sapply(1:n(), function(x) substr(as.character(md5(as.character(x))), 1, 12))) %>%
  # combine EPA patient number, sample_num, and sample_hash into the final ID
  mutate(sample_name_hash = paste0(patient, "_", sample_num, "_", sample_hash),
         sample_name = paste0(patient, "_", sample_num))

# broadcast germline sample info to all samples
map_file <- map_file %>% 
  group_by(submitted_subject_id) %>%
  mutate(
    gl_sample_name_hash = if_else(
      !germline,
      sample_name_hash[germline_use][1],
      NA_character_
    )
  ) %>%
  ungroup()

# add addition info, reformat and select relevant columns
map_file <- map_file %>%
  mutate(
    sample_category = case_when(
      Run %in% normal_tissue ~ "normal_tissue",
      Run %in% pdac_samples ~ "pdac_samples",
      Run %in% whole_slide ~ "whole_slide",
      Run %in% acinar_tissue ~ "acinar_tissue",
      Run %in% germline_controls ~ "germline_controls",
      TRUE ~ "preinvasive_panin"
    )
  )

cols_to_keep <- c(
  "sample_name",
  "sample_name_hash",
  "gl_sample_name_hash",
  "seq_id",
  "sequencer",
  "patient",
  "run_dir",
  "sample_type",
  "clinical_sex",
  "ffpe",
  "germline",
  "germline_use",
  "seq_type",
  "sample_class",
  "sample_type_detailed",
  "sample_type_simple",
  "date_sequenced",
  "surgical_sample",
  "Run",
  "BioProject",
  "BioSample",
  "biospecimen_repository",
  "Experiment"
)
map_file_clinical <- map_file %>%
  rename(seq_id = biospecimen_repository_sample_id) %>%
  mutate(
    seq_type = "WES",
    sample_class = ifelse(germline, "gl", "ffpe"),
    run_dir = "HTAN_Pilot_Pancreas",
    sequencer = "Nextseq",
    ffpe = ifelse(sample_class == "ffpe", TRUE, FALSE),
    date_sequenced = format(as.Date(ReleaseDate), "%d/%m/%Y"),
    clinical_sex = sex,
    surgical_sample = ifelse(Run %in% germline_controls, FALSE, TRUE),
    sample_type_simple = sample_type,
    sample_type_detailed = sample_category,
  ) %>%
  select(cols_to_keep)

# write out map file
write_delim(map_file_clinical, paste0(outputs.folder, date, "_HTAN_Pilot_Precancer_Pancreatic_metadata_epa_mapping.tsv"), delim = "\t")

#=======================#
# Add data to inventory #
#=======================#

# first, write out the inventory unchanged but with .bak so that the last safe
# copy is preserved in case something needs correcting
write_delim(latest_inventory, "/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.bak.tsv", delim = "\t")

# start with fastq seq info and join with map file columns
fastq_metadata_alt <- fastq_seq_dat_alt %>%
  rename(fastq_sample = sample) %>%
  left_join(map_file_clinical, by = "Run")

# check each biosample is only associated with one flow cell and lane
lane_count_sanity_check <- fastq_metadata_alt %>%
  group_by(sample_name_hash) %>%
  summarise(
    rows = n(),
    unique_lane_fc = n_distinct(flowcell, lane),
    flowcells = n_distinct(flowcell),
    .groups = "drop"
  ) %>%
  arrange(desc(unique_lane_fc))

# add a couple of extra columns and select the inventory ones
inventory_new <- fastq_metadata_alt %>%
  mutate(
    # above check showed that all samples are only sequenced over one lane
    passing_lanes = 1,
    uid = paste0(run_dir, "_", seq_id, "_", flowcell, "_", lane),
    fq_pass = TRUE
  ) %>%
  select(colnames(latest_inventory))

# write out inventory for this study
write_delim(inventory_new, paste0(outputs.folder, date, "_epa-project.inventory.latest.pass.tsv"), delim = "\t")

# append to previous inventory
inventory_updated <- rbind(latest_inventory, inventory_new)

# write out updated inventory
write_delim(inventory_updated, "/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.tsv", delim = "\t")











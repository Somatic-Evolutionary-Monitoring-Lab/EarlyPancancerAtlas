#!/usr/bin/env Rscript

# Preprocessing metadata for HTAN Pilot Preinvasive Breast dataset
# In this data it is unclear if most germlines are missing or if there is mislabelling
# going to perform an initial run of ASCAT before full pipeline to check for 
# more potential germlines, so for now going to save to a dummy/tmp inventory file


# Load libraries
library(dplyr)
library(readr)
library(stringr)
library(openssl)

setwd('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/data_downloads/HTAN_Precancer_Pilot/Breast/fastqs')

# Load data
sra_metadata <- read_csv('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/data_downloads/HTAN_Precancer_Pilot/Breast/fastqs/HTAN_Pilot_Breast_WES_Metadata.csv')
fastq_info <- read_tsv('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/data_downloads/HTAN_Precancer_Pilot/Breast/fastqs/HTAN_Precancer_Pilot_Breast_seq_info.tsv')
latest_inventory <- read_tsv('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.tsv')
chen_example <- read_tsv('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/data_downloads/HTAN_Precancer_Pilot/Breast/20250514_chen_2017_ESCC_metadata_EPA_mapping.txt')


# Default values
process_dir <- '/rds/project/rds-LH0AvU65IRI/datasets/EPA/Breast/fastqs/'

## Process metadata files

# some SRA data does not come with flowcell or lane info in fastq header for some reason
fastq_info$flowcell <- NA
fastq_info$lane <- NA
fastq_info$sra <- str_split_i(fastq_info$sample, "_", 1)
fastq_info$fq1 <- paste0(process_dir, basename(fastq_info$r1))
fastq_info$fq2 <- paste0(process_dir, basename(fastq_info$r2))


## Right now only 5 samples are labelled as germline. Need to do a minimal ASCAT
# run on this data to work out if any other samples are germline so will create
# EPA IDs but these may have to change

# format dates
dates <- as.Date(sub(" UTC", "", sra_metadata$ReleaseDate))
formatted_dates <- format(dates, "%-m/%-d/%y")

# get columns from SRA metadata
sample_info <- data.frame(
  sra = sra_metadata$Run,
  seq_id = sra_metadata$biospecimen_repository_sample_id,
  sequencer = sra_metadata$Instrument,
  run_dir = "HTAN_Precancer_Pilot_Breast",
  clinical_sex = sra_metadata$sex,
  passing_lanes = 1,
  seq_type = "WES",
  fq_pass = TRUE,
  date_sequenced = formatted_dates,
  uid = paste(sra_metadata$biospecimen_repository_sample_id, "_HTAN_Precancer_Breast_WES"),
  histological_type = sra_metadata$histological_type,
  study_participant = sra_metadata$submitted_subject_id
)

# join with fastq level info
sample_info <- sample_info %>%
  left_join(fastq_info %>% select(sra, fq1, fq2, flowcell, lane), by = "sra")

# annotate sample types and assign EPA IDs
sample_info <- sample_info %>%
  mutate(germline = case_when(
    histological_type == "Normal" ~ TRUE,
    TRUE ~ FALSE
  ),
  sample_type = case_when(
    histological_type == "Normal" ~ "normal",
    # for now, assume everything else is DCIS but will need to check with ASCAT
    TRUE ~ "preinvasive"
  )
  ) %>%
  group_by(study_participant) %>%
  mutate(
    # assign sequential identifiers for germline samples
    GL_counter = as.integer(dense_rank(ifelse(germline, seq_id, NA))),
    # assign independent counters for tumour samples
    S_counter = as.integer(dense_rank(ifelse(!germline, seq_id, NA)))
  ) %>%
  # add sample ids
  mutate(sample_id = ifelse(germline, paste0("GL", GL_counter), paste0("S", S_counter))) %>%
  ungroup() %>%
  select(-GL_counter, -S_counter)

# generate unique pre check EPA IDs
curr_max_epa <- max(as.integer(sub("^EPA(\\d+).*", "\\1", latest_inventory$sample_name)))

sample_info <- sample_info %>%
  mutate(epa_number = sprintf("EPA%05d", dense_rank(study_participant) + curr_max_epa)) %>%
  # get random unique hash for each sample
  mutate(hash = sapply(1:n(), function(x) substr(as.character(md5(as.character(x))), 1, 12))) %>%
  # combine EPA number, sample_id, and hash into the final ID
  mutate(EPA = paste0(epa_number, "_", sample_id, "_", hash))

# get list of patients with no assigned germline to check with ASCAT
participants_no_germline <- sample_info %>%
  group_by(study_participant) %>%
  summarise(has_germline = any(germline == TRUE, na.rm = TRUE)) %>%
  filter(!has_germline)

samples_to_check_ascat <- sample_info %>%
  filter(study_participant %in% participants_no_germline$study_participant) 

write_tsv(samples_to_check_ascat, "samples_no_germline.tsv")

# broadcast germline sample info to all samples
sample_info <- sample_info %>% 
  mutate(germline_use = germline & grepl("GL1", EPA)) %>%
  group_by(study_participant) %>%
  mutate(
    gl_sample_name_hash = if_else(
      !germline,
      EPA[germline_use][1],
      NA_character_
    )
  ) %>%
  ungroup()

# fill in remaining inventory columns and select only required cols
inventory <- sample_info %>%
  mutate(sample_name = paste(epa_number, sample_id, sep = "_"),
         sample_name_hash = EPA,
         patient = epa_number,
         sample_class = ifelse(germline, "gl", "ffpe"),
         ffpe = ifelse(sample_class == "ffpe", T, F)) %>%
  select(colnames(latest_inventory))
  
# save temporary inventory to the run directory
write_tsv(inventory, "/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/run_dir_align/20260303_HTAN_breast/temp_inventory_HTAN_Precancer_Breast.tsv")


# for now dont add this data to the main inventory because may need to rename if
# we discover some are germlines







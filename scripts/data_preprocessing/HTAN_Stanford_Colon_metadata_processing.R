#!/usr/bin/env Rscript

## Annotate HTAN Stanford Colon Project Samples

# Load libraries
library(dplyr)
library(readr)
library(stringr)
library(openssl)

# Load metadata files and inventory
biosamples <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/data_downloads/HTAN/Stanford_colon/HTAN_Stanford_colon_Precancer_Biosample_Metadata.tsv")
fastq_metadata <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/data_downloads/HTAN/Stanford_colon/HTAN_Stanford_colon_Precancer_fq_Metadata.tsv")
fastq_seq_dat <- read_tsv('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/data_downloads/HTAN/Stanford_colon/HTAN_Stanford_colon_Precancer_seq_platform_info.tsv', col_names = T)
cases <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/data_downloads/HTAN/Stanford_colon/HTAN_Stanford_colon_Precancer_case_level_data.tsv")
latest_inventory <- read_tsv('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.tsv')


## Process metadata files

# Keep relevent columns
cols_to_keep <- c("Biospecimen", "HTAN Parent ID", "Timepoint Label", "Biospecimen Type",
                  "Acquisition Method Type", "Participant ID", "Tumor Tissue Type")

clinical_cols_to_keep <- c("HTAN Participant ID", "Age at Diagnosis (years)", "Primary Diagnosis",
                           "Ethnicity", "Gender", "Race")              

biosample_metadata_alt <- biosamples %>%
  select_if(~!all(is.na(.))) %>%
  rename(Biospecimen = `HTAN Biospecimen ID`) %>%
  select(all_of(cols_to_keep))

clinical_dat_alt <- cases %>%
  select_if(~!all(is.na(.))) %>%
  select(all_of(clinical_cols_to_keep)) %>%
  rename(`Participant ID` = `HTAN Participant ID`)

fastq_metadata_alt <- fastq_metadata %>%
  select_if(~!all(is.na(.)))


# Join biosample info with fastq metadata to annotate
fastq_biosample_dat <- fastq_metadata_alt %>%
  left_join(biosample_metadata_alt, by = "Biospecimen") %>%
  left_join(clinical_dat_alt, by = "Participant ID")


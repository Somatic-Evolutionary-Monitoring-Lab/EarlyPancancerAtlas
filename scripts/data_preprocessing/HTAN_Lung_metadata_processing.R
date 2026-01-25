#!/usr/bin/env Rscript

##### this script contains errors in the way I handled the data annotation for
# this project, have corrected these downstream but keeping this script here 
# with additional comments inserted as a record. K Honan 20260123. ##### 

library(dplyr)
library(readr)
library(stringr)
library(openssl)

rcs_path <- "/rcs/project/zz485/rcs-zz485-sem-lab-cold/pancancer_kh/Lung/HTAN/WES/"

fastq_metadata <- read_tsv(paste0(rcs_path, 'HTAN_Lung_PreCancer_fq_Metadata.tsv'))
fastq_seq_dat <- read_tsv(paste0(rcs_path, 'HTAN_Lung_PreCancer_seq_platform_info.tsv'), col_names = F)
biosample_metadata <- read_tsv(paste0(rcs_path, 'HTAN_Lung_PreCancer_Biosample_Metadata.tsv'))
clinical_dat <- read_tsv(paste0(rcs_path, 'HTAN_Lung_PreCancer_case_level_data.tsv'))
latest_inventory <- read_tsv('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.tsv')


## Process metadata files

# Keep relevent columns
cols_to_keep <- c("Biospecimen", "HTAN Parent ID", "Timepoint Label", "Biospecimen Type",
                  "Acquisition Method Type", "Participant ID", "Adjacent Biospecimen IDs",
                  "Preinvasive Morphology", "Degree of Dysplasia", "Tumor Tissue Type")

clinical_cols_to_keep <- c("HTAN Participant ID", "Age at Diagnosis (years)", "Primary Diagnosis",
                           "Ethnicity", "Gender", "Race", "Prior Malignancy", "Prior Treatment",
                           "AJCC Clinical Stage", "AJCC Clinical T", "AJCC Pathologic M",
                           "AJCC Pathologic N", "AJCC Pathologic Stage", "AJCC Pathologic T",
                           "AJCC Staging System Edition")                

biosample_metadata_alt <- biosample_metadata %>%
  select_if(~!all(is.na(.))) %>%
  rename(Biospecimen = `HTAN Biospecimen ID`) %>%
  select(all_of(cols_to_keep))

fastq_metadata_alt <- fastq_metadata %>%
  select_if(~!all(is.na(.)))

clinical_dat_alt <- clinical_dat %>%
    select_if(~!all(is.na(.))) %>%
    select(all_of(clinical_cols_to_keep)) %>%
    rename(`Participant ID` = `HTAN Participant ID`)

# Join biosample info with fastq metadata to annotate
fastq_biosample_dat <- fastq_metadata_alt %>%
    left_join(biosample_metadata_alt, by = "Biospecimen") %>%
    left_join(clinical_dat_alt, by = "Participant ID")


## Create EPA ID mapping and number the germline, preinvasive and tumour samples per patient

# define germline conditions for ambigous cases
# add column to indicate germline status
fastq_biosample_dat <- fastq_biosample_dat %>% 
  mutate(
    germline = case_when(
      # Explicit normal calls
      `Tumor Tissue Type` %in% c("Normal adjacent", "Normal") ~ TRUE,
      `Biospecimen Type` == "Blood Biospecimen Type" ~ TRUE,
      # Rescue nos cases if other cols suggest normal tissue -- delete this
      `Tumor Tissue Type` == "Not Otherwise Specified" &
        `Preinvasive Morphology` == "Normal WDA" &
        str_detect(
          str_to_lower(`Degree of Dysplasia`),
          "normal"
        ) ~ TRUE,
      # All other NOS remain NA
      `Tumor Tissue Type` == "Not Otherwise Specified" ~ NA,
      # Everything else is premalignant/tumour
      TRUE ~ FALSE
    )
  ) %>%
  # filter out case where we don't have enough info to determine sample type
  filter(!is.na(germline)) %>%
  mutate(
    fastq_file = basename(Filename),
    sample_type = case_when(
      # Explicit primary tumour calls
      `Tumor Tissue Type` %in% c("Primary", "Additional Primary") ~ "cancer",

      # Explicit preinvasive calls
      `Tumor Tissue Type` %in% c(
        "Atypia - hyperplasia",
        "Premalignant - in situ",
        "Premalignant"
      ) ~ "preinvasive",

      # Normal calls
      germline == TRUE ~ "normal",

      # Everything else
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(`Participant ID`) %>%
  mutate(
    # assign sequential identifiers for germline samples
    GL_counter = as.integer(dense_rank(ifelse(germline, `Biospecimen`, NA))),
    # assign independent counters for tumour samples
    S_counter = as.integer(dense_rank(ifelse(!germline, `Biospecimen`, NA)))
  ) %>%
  # add sample ids
  mutate(sample_id = ifelse(germline, paste0("GL", GL_counter), paste0("S", S_counter))) %>%
  ungroup() %>%
  select(-GL_counter, -S_counter)

# **** 23-01-26: decided HTA3_50044_1001588 was too ambiguous to call germline  
# so removed this and other samples associated with HTA3_50044 from inventory
# manually after running this script ****

# **** 23-01-26: Made an error here while processing this dataset. Should not
# have assigned EPA IDs until after I matched up the fq1 and fq2 files to be on
# the same line. This basically duplicated the number of fq pairs per sample,
# see script inventory_fix_20260123.R for details/ code to fix this ****

# generate unique EPA IDs
curr_max_epa <- max(as.integer(sub("^EPA(\\d+).*", "\\1", latest_inventory$sample_name)))

fastq_biosample_dat <- fastq_biosample_dat %>%
  mutate(epa_number = sprintf("EPA%05d", dense_rank(`Participant ID`) + curr_max_epa)) %>%
  # get random unique hash for each sample
  mutate(hash = sapply(1:n(), function(x) substr(as.character(md5(as.character(x))), 1, 12))) %>%
  # combine EPA number, sample_id, and hash into the final ID
  mutate(EPA = paste0(epa_number, "_", sample_id, "_", hash))

# check if there are any samples with no germline
participants_no_germline <- fastq_biosample_dat %>%
  group_by(`Participant ID`) %>%
  summarise(has_germline = any(germline == TRUE, na.rm = TRUE)) %>%
  filter(!has_germline)

# exclude all these cases
exclude_list <- unique(participants_no_germline$`Participant ID`)

fastq_biosample_dat <- fastq_biosample_dat %>%
  filter(!`Participant ID` %in% exclude_list)

# broadcast germline sample info to all samples
fastq_biosample_dat <- fastq_biosample_dat %>% 
  mutate(germline_use = germline & grepl("GL1", EPA)) %>%
  group_by(`Participant ID`) %>%
  mutate(
    gl_sample_name_hash = if_else(
      !germline,
      EPA[germline_use][1],
      NA_character_
    )
  ) %>%
  ungroup()

# mark repeats
fastq_biosample_dat <-  fastq_biosample_dat %>%
  group_by(Biospecimen) %>%
  mutate(repeat_sample = row_number() > 2) %>%
  ungroup() %>%
  filter(!repeat_sample)


# Write out file that contains fastq file level info on samples
# *** 20260123 The below file has now been regenerated to reflect the corrections  
# made to the inventory for this project i.e., redundant EPA IDs have been removed ***
write_delim(fastq_biosample_dat, paste0(rcs_path, "20260105_HTAN_Lung_metadata_epa_mapping.tsv"), delim = '\t')

## Inventory
# append study information to inventory file

fastq_seq_dat <- fastq_seq_dat %>%
  rename(sample = X1,
  fq1_path = X2,
  fq2_path = X3,
  fq1_size = X4,
  fq2_size = X5,
  flowcell = X6,
  lane = X7) %>%
  mutate(fastq_file = basename(fq1_path)) %>%
  left_join(fastq_biosample_dat %>% select(fastq_file, Biospecimen), by = "fastq_file")

fastq_biosample_dat <- fastq_biosample_dat %>%
  left_join(fastq_seq_dat %>% select(Biospecimen, fq1_path, fq2_path, flowcell, lane), by = "Biospecimen") %>%
  mutate(Gender = ifelse(Gender == "Not Reported", "unknown", Gender))
  
study_inventory <- data.frame(sample_name = gsub("_[^_]+$", "", fastq_biosample_dat$EPA),
                              sample_name_hash = fastq_biosample_dat$EPA,
                              gl_sample_name_hash = fastq_biosample_dat$gl_sample_name_hash,
                              seq_id = fastq_biosample_dat$Biospecimen,
                              sequencer = "Novaseq",
                              patient = gsub("_.*", "", fastq_biosample_dat$epa_number),
                              run_dir = "HTAN_Lung",
                              clinical_sex = tolower(fastq_biosample_dat$Gender),
                              germline = fastq_biosample_dat$germline,
                              germline_use = fastq_biosample_dat$germline_use,
                              passing_lanes = 1,
                              seq_type = "WES",
                              fq_pass = TRUE,
                              sample_class = ifelse(fastq_biosample_dat$germline, "gl", "ffpe"),
                              fq1 = fastq_biosample_dat$fq1_path,
                              fq2 = fastq_biosample_dat$fq2_path,
                              lane = as.character(fastq_biosample_dat$lane),
                              flowcell = fastq_biosample_dat$flowcell,
                              date_sequenced = as.Date("2000-01-01"),
                              ffpe = ifelse(!fastq_biosample_dat$germline, TRUE, FALSE),
                              uid = paste("HTAN_Lung",
                              fastq_biosample_dat$Biospecimen,
                              fastq_biosample_dat$flowcell,
                              fastq_biosample_dat$lane, sep = "_")
                              )

latest_inventory <- bind_rows(latest_inventory, study_inventory)

# 20260123 the below file was updated after mistakes in this script were discovered.
# write_delim(latest_inventory, "/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.tsv", "\t")



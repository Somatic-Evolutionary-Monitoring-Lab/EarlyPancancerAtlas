#==============================================================================#
#==============================================================================#
######                                                                    ######
######               Creating a cohort level metadata df                  ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Katherine Honan
# Date: 2025-06-12

setwd("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/")

####################################################
#### Source required functions & load libraries ####
####################################################

library(readr)
library(tidyr)
library(dplyr) 
library(stringr)

#############################################
#### Make a folder for this analysis run ####
#############################################

date <- gsub("-","",Sys.Date())

analysis_name <- 'cohort_metadata'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep='/')
if( !file.exists(out_dir_general) ) dir.create( out_dir_general )

out_dir_logs <- paste(out_dir_general, 'logs', sep='/')
if( !file.exists(out_dir_logs) ) dir.create( out_dir_logs )

outputs.folder <- paste0( out_dir_general, "/", date, "/" )

if( !file.exists(outputs.folder) ) dir.create( outputs.folder )


##############################################
#### Get Inputs required for all analyses ####
##############################################

# read in individual study sample metadata/ clinical data

# ESCC
chen_2017_ESCC_sample_metadata <- read_tsv("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/inputs/20250514_chen_2017_ESCC_metadata_EPA_mapping.txt")
chen_2017_ESCC_clinical <- read_tsv("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/inputs/Chen_2017_ESCC_clinical_info.txt")

# EAC
E_black_clinical <- read_tsv("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/inputs/20250610_clinical_data.tsv")
E_black_metadata <- read_tsv("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/inputs/20250609_E_Black_EAC_WES_metadata_mapping.txt")


##############################################
####  Function for counting sample types  ####
##############################################

count_samples_by_type <- function(metadata_df, cancer_label = "cancer", preinv_label = "preinv") {
  metadata_df %>%
    filter(sample_type %in% c(cancer_label, preinv_label)) %>%
    distinct(patient, uid, sample_type) %>%
    group_by(patient, sample_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = sample_type,
      values_from = n,
      values_fill = 0,
      names_prefix = "num_"
    ) %>%
    rename_with(
      ~ ifelse(. == paste0("num_", cancer_label), "num_tumour_samples",
               ifelse(. == paste0("num_", preinv_label), "num_preinv_samples", .))
    )
}


################################################
####  Define project specific sample types  ####
################################################

normal <- c("normal", "germline", "normal_tissue")

preinv <- list(escc = c("dysplasia", "dysplasia_high", "dysplasia_low", "dysplasia"),
               eac = c("IM","GM","GM_IM","IMC","DYS_GM","barretts","DYS_IMC"))

cohort_df_cols <- c("study_patient", "patient", "age", "plt_gender", "plt_smoking",
                    "plt_stage", "plt_cohort_cancer", "num_tumour_samples", "num_preinv_samples")

##################################################
####      Preprocess data for plotting        ####
##################################################

##################################################
################      ESCC        ################
##################################################

# Add EPA ID to metadata
chen_2017_ESCC_sample_metadata <- chen_2017_ESCC_sample_metadata %>%
  rename(sample_name_hash = EPA) %>%
  mutate(EPA = gsub("_.*", "", sample_name_hash)) %>%
  mutate(sample_type = case_when(
    descriptor == "carcinoma" ~ "cancer",
    descriptor == "metastasis_lymphnode" ~ "cancer",
    descriptor %in% normal ~ "normal",
    descriptor %in% preinv$escc ~ "preinv",
    TRUE ~ NA_character_))

# get sample type counts
chen_2017_ESCC_sample_metadata$patient = chen_2017_ESCC_sample_metadata$EPA
chen_2017_sample_counts <- count_samples_by_type(chen_2017_ESCC_sample_metadata)

# add EPA IDs to clinical data
chen_2017_combined <- chen_2017_ESCC_clinical %>%
  rename(study_patient = Patient, age = Age) %>%
  left_join(chen_2017_ESCC_sample_metadata %>% 
              select(study_patient, EPA) %>% 
              group_by(EPA) %>%
              slice_head(n = 1), 
            by = "study_patient") %>%
  rename(patient = EPA) %>%
  #left_join(chen_2017_sample_counts, by = c("EPA" = "patient")) %>%
  left_join(chen_2017_sample_counts, by = "patient") %>%
  mutate(Stage = ifelse(Stage == "T1bN1M0/_B", "IIB", Stage), # fix one 
         plt_stage = gsub(".*/", "", Stage),
         plt_smoking = case_when(is.na(Tobacco) ~ "Unknown", grepl("Y", Tobacco) ~ "Yes", TRUE ~ "No"),
         plt_gender = case_when(Gender == "M" ~ "Male", Gender == "F" ~ "Female"),
         plt_cohort_cancer = "ESCC") %>%
  select(all_of(cohort_df_cols))
  

##################################################
################      EAC        #################
##################################################

E_black_metadata <- E_black_metadata %>%
  # normalise sample type columns
  mutate(sample_type = case_when(
    specimen_phenotype == "oac" ~ "cancer",
    specimen_phenotype %in% normal ~ "normal",
    specimen_phenotype %in% preinv$eac ~ "preinv",
    TRUE ~ NA_character_))

# count number of cancer and preinvasive samples for each patient
E_black_sample_counts <- count_samples_by_type(E_black_metadata)

# add EPA IDs to clinical data
E_black_combined <- E_black_clinical %>%
  rename(study_patient = OCCAMS_ID, age = DI_ageAtDiagnosis) %>%
  # add EPA patient number
  left_join(E_black_metadata %>%
              select(case_id, patient) %>%
              group_by(patient) %>%
              slice_head(n = 1),
            by = c("study_patient" = "case_id")) %>%
  # 3 duplicate cases, remove for now
  distinct(patient, .keep_all = TRUE) %>%
  left_join(E_black_sample_counts, by = "patient") %>%
  # remove samples with no tumour
  filter(!is.na(num_tumour_samples)) %>%
  mutate(plt_gender = str_to_sentence(DI_PatientGender),
         plt_smoking = case_when(EX_IsSmoker %in% c("current", "former") ~ "Yes",
                                 EX_IsSmoker == "never" ~ "No",
                                 EX_IsSmoker %in% c("unknown", NA) ~ "Unknown"),
         # Normliase TNM columns
         T_clean = case_when(
           str_detect(PS_TStage_PrimaryTumour_FinalPretreatmentStaging, "^T[0-9a-b]+$") ~ PS_TStage_PrimaryTumour_FinalPretreatmentStaging,
           TRUE ~ NA_character_),
         N_clean = case_when(
           str_detect(PS_NStage_PrimaryTumour_FinalPretreatmentStaging_TNM7, "^N[0-3]$") ~ PS_NStage_PrimaryTumour_FinalPretreatmentStaging_TNM7,
           TRUE ~ NA_character_),
         M_clean = case_when(
           PS_MStage_PrimaryTumour_FinalPretreatmentStaging == "M0" ~ "M0",
           PS_MStage_PrimaryTumour_FinalPretreatmentStaging == "Mx" ~ "M1",
           TRUE ~ NA_character_),
         # Map TNM staging to roman numeral formal
         plt_stage = case_when(
           is.na(T_clean) | is.na(N_clean) | is.na(M_clean) ~ "Not recorded",
           M_clean == "M1" ~ "IV",
           T_clean %in% c("T1", "T1a", "T1b") & N_clean == "N0" ~ "IA",
           T_clean == "T2" & N_clean == "N0" ~ "IB",
           T_clean %in% c("T3") & N_clean == "N0" ~ "IIA",
           T_clean %in% c("T2") & N_clean == "N1" ~ "IIB",
           T_clean == "T3" & N_clean == "N1" ~ "IIIA",
           T_clean == "T4" & N_clean == "N1" ~ "IIIB",
           N_clean %in% c("N2", "N3") ~ "IIIC",
           TRUE ~ "Not recorded"),
         # add cohort cancer type
         plt_cohort_cancer = "EAC") %>%
  select(all_of(cohort_df_cols))



################################################################
################      Cohort Level Data        #################
################################################################

cohort_data <- bind_rows(chen_2017_combined, E_black_combined)

write_tsv(cohort_data, paste0(outputs.folder, "full_cohort_data.tsv"))





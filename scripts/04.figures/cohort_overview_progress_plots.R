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
library(tidytext)

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

map_dir <- "/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/map_files/"
suffix <- "*metadata_EPA_mapping.tsv"

inventory <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.tsv")
htan_lung_map <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/map_files/20260123_HTAN_Lung_metadata_EPA_mapping.tsv")
chen_escc_map <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/map_files/20250514_chen_2017_ESCC_metadata_EPA_mapping.txt")
black_eac_map <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/map_files/20250609_E_Black_EAC_WES_metadata_EPA_mapping.txt")
htan_stanford_colon_map <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/map_files/20260313_HTAN_Stanford_Colon_metadata_epa_mapping.tsv")
htan_pancreas_map <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/map_files/20260314_HTAN_Pilot_Precancer_Pancreatic_metadata_epa_mapping.tsv")
phs001424_brain_map <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/map_files/20260327_phs001424_Brain_Brastianos_metadata_epa_mapping.tsv")
phs001460_brain_map <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/map_files/20260331_phs001460_Brain_metadata_epa_mapping.tsv")
phs002612_brain_map <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/map_files/20260401_phs002612_Brain_metadata_epa_mapping.tsv")

# source colour palettes
source("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/scripts/04.figures/EPA_colour_palletes.R")

#====================#
# Data preprocessing #
#====================#

cohort_samples <- inventory %>%
  filter(!germline) %>%
  select(
    sample_name_hash,
    patient,
    run_dir,
    sample_type_simple
    ) %>%
  distinct() %>%
  mutate(
    organ = case_when(
      grepl("ESCC", run_dir) ~ "Esophageal (ESCC)",
      grepl("EAC", run_dir) ~ "Esophageal (EAC)",
      grepl("Lung", run_dir) ~ "Lung",
      grepl("Colon", run_dir) ~ "Colon",
      grepl("Pancreas", run_dir) ~ "Pancreas",
      grepl("Brain", run_dir) ~ "Brain",
      TRUE ~ NA
    )
  ) %>%
  filter(!is.na(sample_type_simple),
         !sample_type_simple == "cancer_met")

sample_type_counts <- cohort_samples %>%
  count(patient, organ, sample_type_simple) %>%
  group_by(organ, patient) %>%
  mutate(total_samples = sum(n)) %>%
  ungroup() %>%
  mutate(organ = factor(organ, levels = sort(unique(organ)))) %>%
  group_by(organ) %>%
  mutate(patient = reorder_within(patient, total_samples, organ))

ggplot(sample_type_counts, aes(x = patient, y = n, fill = sample_type_simple)) +
  geom_bar(stat = "identity", colour = "black", width = 0.8) +
  scale_fill_manual(values = c(preinvasive = unname(lesion_type_pallete[2]), cancer = unname(lesion_type_pallete[3])))+
  facet_grid(~ organ, scales = "free_x", space = "free_x") +
  scale_x_reordered() +
  labs(
    x = "Patient",
    y = "Number of samples",
    fill = "Sample type"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )

ggsave(paste0(outputs.folder, date, "_cohort_sample_overview.png"), width=40, height=10)

#=====================================================#
# Generate cohort data quality plots for aligned data #
#=====================================================#











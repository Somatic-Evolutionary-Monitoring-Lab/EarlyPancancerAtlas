#!/usr/bin/env Rscript
# Author: Katherine Honan
# Date: 2026-03-06
# Description: Analysis of T Cell Fraction in EPA

#================#
# Load libraries #
#================#

library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lme4)
library(lmerTest)
library(stringr)
library(ggeffects)

options(scipen=999)

#=========================#
# Create output directory #
#=========================#

setwd("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/")

date <- gsub("-","",Sys.Date())

analysis_name <- 'T_cell_fraction'
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

cohort_tcra <- read_tsv('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/_RELEASE/release-20260306/_aggregate/tcell_extract/cohort_TCellExTRECT.tsv')
chen_ESCC_metadata <- read_tsv('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/map_files/20250514_chen_2017_ESCC_metadata_EPA_mapping.txt')
clone_info <- readRDS('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/outputs/preinvasive_seeding/20260306/20260306_cloneInfo.list.rds')

source('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/scripts/EPA_colour_palletes.R')

#=================#
# Data processing #
#=================#

cohort_tcra_df <- cohort_tcra %>%
  left_join(chen_ESCC_metadata %>% select(EPA, descriptor), by = c("sample" = "EPA")) %>%
  # exclude mets
  filter(!descriptor == "metastasis_lymphnode")

# make all normals called germline
cohort_tcra_df[cohort_tcra_df$descriptor == "normal_tissue", ]$descriptor <- "germline"

# for now just classify any dysplasia as preinvasive
any_dys <- unique(cohort_tcra_df$descriptor)[grepl("dysplasia", unique(cohort_tcra_df$descriptor))]
cohort_tcra_df[cohort_tcra_df$descriptor %in% any_dys, ]$descriptor <- "preinvasive"


#=========================================================#
# Analysis 1: Overview plot of T cell fractions in cohort #
#=========================================================#

# T cell fractions per tumour - need to incorperate info about which samples 
# from which tumours (where multiple clonally unrelated)
overview_plot_df <- cohort_tcra_df %>%
  filter(!descriptor == "germline") %>%
  mutate(tumour = str_split_i(sample, "_", 1)) %>%
  group_by(tumour) %>%
  mutate(tumour_mean = mean(TCRA.tcell.fraction)) %>%
  ungroup() %>%
  mutate(tumour = reorder(tumour, tumour_mean))

# look at distributions grouped by sample type
ggplot(cohort_tcra_df, aes(x = descriptor, y = TCRA.tcell.fraction, fill = descriptor)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  labs(
    x = "",
    y = "T cell fraction",
    title = "T cell infiltration in different sample types"
  ) +
  scale_fill_manual(values = c(
    "germline" = unname(lesion_type_pallete["germline"]),
    "preinvasive" = unname(lesion_type_pallete["preinvasive"]),
    "carcinoma" = unname(lesion_type_pallete["tumour"])
  )) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(paste0(outputs.folder, date, "_Tcell_frac_by_sample_type.png"), width = 6, height = 6)

#==============================================================#
# Analysis 2: T cell fraction in preinvasive vs tumour samples #
#==============================================================#

# Is there a difference in T cell fraction between tumour and preinvasive samples?
tumour_vs_preinv_df <- cohort_tcra_df %>%
  filter(descriptor %in% c("carcinoma", "preinvasive")) %>%
  mutate(descriptor = factor(descriptor, levels = c("preinvasive", "carcinoma")))

# stats test between tumour and preinvasive
wilcox.test(TCRA.tcell.fraction ~ descriptor, data = tumour_vs_preinv_df)

#=====================================================================#
# Analysis 2: T cell fraction in initiating vs non-initiating samples #
#=====================================================================#

# get initiating clones from preinvasive seeding script output
initiating <- lapply(names(clone_info), function(pat){

  init_region_info <- clone_info[[pat]]$seedingRegionInfo
  
  if (is.null(init_region_info) || !is.data.frame(init_region_info)) {
    return(NA)
  }
  
  init_regions <- init_region_info %>%
    filter(Initiating) %>%
    pull(PreinvasiveRegion)
  return(init_regions)
})

initiating_regions <- as.character(na.omit(unlist(initiating)))

# get non initiating regions
cohort_tcra_init <- cohort_tcra_df %>%
  filter(descriptor == "preinvasive") %>%
  mutate(initiating = ifelse(sample %in% initiating_regions, TRUE, FALSE)) %>%
  select(sample, TCRA.tcell.fraction, initiating)

# compare 
initiating_vs_not_df <- cohort_tcra_init %>%
  filter(!is.na(initiating)) %>%
  mutate(
    TCRA.tcell.fraction = ifelse(TCRA.tcell.fraction <= 0, 1e-6, TCRA.tcell.fraction),
    initiating = factor(initiating, levels = c(FALSE, TRUE),
                        labels = c("Non-initiating", "Initiating"))
  )

wilcox.test(TCRA.tcell.fraction ~ initiating, data = initiating_vs_not_df)

ggplot(initiating_vs_not_df, aes(x = initiating, y = TCRA.tcell.fraction, fill = initiating)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  labs(
    x = "",
    y = "T cell fraction",
    title = "T cell infiltration in initiating vs \nnon-initiating regions"
  ) +
  stat_compare_means(method = "wilcox.test") +
  theme_classic() +
  theme(legend.position = "none")


ggsave(paste0(outputs.folder, date, "_Tcell_frac_by_initiation_status.png"), width = 5, height = 6)

# Are initiating regions depleted for detectable T cell infiltration?
initiating_vs_not_df <- initiating_vs_not_df %>%
  mutate(tcell_present = TCRA.tcell.fraction > 1e-5)

f_test <- fisher.test(table(initiating_vs_not_df$initiating, initiating_vs_not_df$tcell_present))
p_label <- paste0("Fisher p = ", signif(f_test$p.value, 2))

# Stacked bar plot of Tcell +ve and Tcell -ve split by initiating or not
tcell_pos_neg_plot_df <- as.data.frame(table(initiating_vs_not_df$initiating, initiating_vs_not_df$tcell_present))
colnames(tcell_pos_neg_plot_df) <- c("initiating", "tcell_status", "count")

tcell_pos_neg_plot_df$tcell_status <- ifelse(
  tcell_pos_neg_plot_df$tcell_status == "TRUE", "T cell Positive", "T cell Negative"
  ) 

tcell_pos_neg_plot_df <- tcell_pos_neg_plot_df %>%
  group_by(initiating) %>%
  mutate(
    total = sum(count),
    prop = count / total,
    label = paste0(count)
  )

ggplot(tcell_pos_neg_plot_df, aes(x = initiating, y = count, fill = tcell_status)) +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label = count),
            position = position_fill(vjust = 0.5),
            size = 4) +
  annotate("text", x = 1.5, y = 1.05, label = p_label, size = 5) +
  labs(
    x = "",
    y = "Proportion of samples",
    fill = "T cell status"
  ) +
  scale_fill_manual(values = c(
    "T cell Positive" = unname(t_cell_infiltration_pallete["present"]),
    "T cell Negative" = unname(t_cell_infiltration_pallete["absent"])
    )) +
  theme_classic()

ggsave(paste0(outputs.folder, date, "_Tcell_frac_by_initiation_stat_binary.png"), width = 5, height = 6)

#=======================================#
# Analysis 3: Mixed logistic regression #
#=======================================#

model <- glmer(tcell_present ~ initiating + (1|tumour),
               data = initiating_vs_not_df,
               family = binomial)

summary(model)


#==========================================================#
# Analysis 4: Investigate patterns across different organs #
#==========================================================#


#==============================================================================#
#==============================================================================#
######                                                                    ######
######          Early Pancancer Atlas Pilot Overview Figures              ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Katherine Honan
# Date: 2025-05-30

setwd("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/")

####################################################
#### Source required functions & load libraries ####
####################################################

library(readr)
library(tidyr)
library(dplyr) 
library(ggplot2) 
library(RColorBrewer) 
library(cowplot)
library(data.table)
library(fst)

#############################################
#### Make a folder for this analysis run ####
#############################################

date <- gsub("-","",Sys.Date())

analysis_name <- 'cohort_overview'
out_name <- 'outputs/plots'
out_dir_general <- paste(out_name, analysis_name, sep='/')
if( !file.exists(out_dir_general) ) dir.create( out_dir_general )

out_dir_logs <- paste(out_dir_general, 'logs', sep='/')
if( !file.exists(out_dir_logs) ) dir.create( out_dir_logs )

outputs.folder <- paste0( out_dir_general, "/", date, "/" )

if( !file.exists(outputs.folder) ) dir.create( outputs.folder )


##############################################
#### Get Inputs required for all analyses ####
##############################################

# read in full cohort data
cohort <- read_tsv("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/outputs/cohort_metadata/20250612/full_cohort_data.tsv")
metadata <- fread("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/inputs/20250514_chen_2017_ESCC_metadata_EPA_mapping.txt")
trees <- readRDS("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/inputs/_RELEASE/_aggregate/conipher_trees/2025_06_28_cohort_conipher_trees.RDS")
mutations <- read.fst("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/inputs/_RELEASE/_aggregate/mutation_table/2025_06_28_cohort_concise_mutation_table.fst")
scna_ith <- read_tsv("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/inputs/_RELEASE/_aggregate/refphase/2025_06_26_cohort_scna_ith_metrics.tsv")
wgd <- read_tsv("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/inputs/_RELEASE/_aggregate/pgddetect/20250628_cohort_GDs_per_tumour.tsv")
seeding <- readRDS("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/outputs/preinveasive_seeding/20250702/20250702_seeding_clonality.rds")

##################################
#### Format data for plotting ####
##################################

patient_tumours <- mutations %>%
  filter(tree_cluster) %>%
  distinct(patient_id, patient_tumour, .keep_all = FALSE)

seeding <- seeding %>%
  mutate(seeding_pattern = paste0(seeding_clonality, " ", seeding_phyletic),
         patient_tumour = sub("_(\\d+)", "_tumour\\1", tumour_id))

# filter and number cohort
cohort <- cohort %>%
  filter(patient %in% patient_tumours$patient_id) %>%
  mutate(plt_stage = ifelse(is.na(plt_stage), "Unknown", plt_stage)) %>%
  left_join(patient_tumours, by = c("patient" = "patient_id")) %>%
  group_by(plt_cohort_cancer) %>%
  arrange(plt_cohort_cancer, plt_stage) %>%
  mutate(patient_order = row_number()) %>%
  ungroup() %>%
  left_join(seeding %>% select(patient_tumour, seeding_pattern), by = "patient_tumour")

##################################
######## Annotate mutations ######
##################################

# get all samples that contributed to tree building
tree_samples <- unique(unlist(lapply(trees, function(x) colnames(x$clonality_table))))
names(trees) <- sapply(trees, function(x) x$parameters$sampleID)

# annotate the patient_tumour that corresponds to each sample

# filter out those that didn't get included in a tree

metadata <- metadata %>%
  mutate(
    sample_type = case_when(
      grepl("dysplasia", descriptor, ignore.case = TRUE) ~ "preinvasive",
      grepl("carcinoma", descriptor, ignore.case = TRUE) ~ "cancer",
      grepl("metastasis_lymphnode", descriptor, ignore.case = TRUE) ~ "cancer",
      descriptor %in% c("germline", "normal_tissue") ~ "normal",
      TRUE ~ NA_character_),
    sample = sub("^(([^_]+)_([^_]+)).*$", "\\1", EPA)
  ) %>%
  rename(sample_name_hash = EPA)


# label mutations as follows:
# 1. pre-malignant private (unique to a single pre-malignant sample)
# 2. pre-malignant shared (common across multiple pre-malignant samples but absent from all tumour samples)
# 3. shared clonal in all pre-malignant and tumour
# 4. shared subclonal in all pre-malignant and tumour
# 5. shared tumour (shared across multiple tumour samples but absent from all pre-malignant samples)
# 6. tumour private (unique to a single tumour sample)

mutations_labelled <- mutations %>% 
  left_join(metadata %>% select(sample_name_hash, sample_type), by = "sample_name_hash") %>%
  filter(!is.na(is_trunk), !is.na(sample_type), mutations$tree_cluster) %>%
  group_by(patient_id, patient_tumour, mutation_id, mutation_cluster) %>%
  summarise(
    #n_pre_samples = n_distinct(sample_name_hash[sample_type == "preinvasive"]),
    #n_cancer_samples = n_distinct(sample_name_hash[sample_type == "cancer"]),
    n_pre_present = sum(is_present & sample_type == "preinvasive"),
    n_cancer_present = sum(is_present & sample_type == "cancer"),
    pre_sample_count = sum(sample_type == "preinvasive"),
    cancer_sample_count = sum(sample_type == "cancer"),
    
    # Clonality info
    all_pre_clonal = all(cluster_region_clonality[sample_type == "preinvasive"] == "clonal", na.rm = TRUE),
    all_cancer_clonal = all(cluster_region_clonality[sample_type == "cancer"] == "clonal", na.rm = TRUE),
    #all_pre_clonal = all(cluster_region_clonality[sample_type == "preinvasive" & is_present] == "clonal", na.rm = TRUE), 
    #all_cancer_clonal = all(cluster_region_clonality[sample_type == "cancer" & is_present] == "clonal", na.rm = TRUE),
    any_pre_present = any(is_present & sample_type == "preinvasive"),
    any_cancer_present = any(is_present & sample_type == "cancer"),
    .groups = "drop"
    ) %>%  
  mutate(
    category = case_when(
    n_pre_present == 1 & n_cancer_present == 0 ~ "premalignant_private",
    n_pre_present > 1 & n_cancer_present == 0 ~ "premalignant_shared",
    any_pre_present & any_cancer_present & all_pre_clonal & all_cancer_clonal ~ "shared_clonal",
    any_pre_present & any_cancer_present & (!all_pre_clonal | !all_cancer_clonal) ~ "shared_subclonal",
    n_pre_present == 0 & n_cancer_present > 1 ~ "tumour_shared",
    n_pre_present == 0 & n_cancer_present == 1 ~ "tumour_private",
    TRUE ~ NA_character_)
  ) #%>%
  #filter(!is.na(category))


mutation_proportions <- mutations_labelled %>%
  count(patient_id, patient_tumour, category) %>%
  group_by(patient_id, patient_tumour) %>%
  mutate(total = sum(n), proportion = n / total) %>%
  select(-n) %>%
  pivot_wider(names_from = category, values_from = proportion, values_fill = 0)


####################################
######## Get copy number info ######
####################################

scna_metrics <- scna_ith %>%
  filter(patient_tumour %in% patient_tumours$patient_tumour) %>%
  select(patient_tumour, frac_abberant_genom_subcl) %>%
  mutate(frac_abberant_genom_clonal = 1 - frac_abberant_genom_subcl)

####################################
########### Get WGD info ###########
####################################

wgd_metrics <- wgd %>%
  filter(tumour_id %in% patient_tumours$patient_tumour) %>%
  rename(patient_tumour = tumour_id) %>%
  select(patient_tumour, num_clonal_gds, num_subclonal_gds)
  

#####################################
####### Create colour mapping #######
#####################################

sample_type_colors <- c("Tumour" = "#009E73", "Pre-invasive" = "#D55E00")

stage_colours <- c("IA" = "#f5f5f5", 
                   "IB" = "#c7eae5", 
                   "IIA" = "#80cdc1", 
                   "IIB" = "#35978f", 
                   "IIIA" = "#01665e", 
                   "IIIB" = "#003c30",
                   "IIIC" = "#002620",
                   "Unknown" = "grey30")

mutation_category_colors <- c(
  "premalignant_private" = "#6A3D9A",
  "premalignant_shared"  = "#CAB2D6",
  "shared_clonal"        = "#1F78B4",
  "shared_subclonal"     = "#A6CEE3",
  "tumour_shared"        = "#B15928",
  "tumour_private"       = "#E31A1C"
)

clonality_colours <- c("Clonal" = "#5F9EA0", "Subclonal" =  "#E6AB02")

gender_colours <- c("Male" = "#377eb8", "Female" = "#CC79A7")

smoking_colours <- c("Yes" = "black", "No" = "#187FC3", "Unknown" = "grey80")

seeding_pattern_colours <- c("monoclonal monophyletic" = "#00A087",
                             "polyclonal monophyletic" = "#3C5488",
                             "polyclonal polyphyletic" = "#F39B7F")

blank_theme <- theme(
  axis.text.x = element_blank(),
  axis.text.y = element_text(size = 14),
  axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5, , size = 14),
  axis.title = element_blank(),
  panel.grid = element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  strip.text = element_blank()
  #panel.border = element_rect(color = "black", fill = NA, size = 0.7)
)

bar_theme <- theme(axis.text.x = element_blank(),
                   axis.ticks = element_blank(),
                   strip.text = element_blank(),
                   axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5, , size = 14),
                   legend.position = "right")


# Mutation bar plot labels
mutation_category_labels <- c(
  "premalignant_private" = "Premalignant private",
  "premalignant_shared"  = "Premalignant shared",
  "shared_clonal"        = "Shared clonal",
  "shared_subclonal"     = "Shared subclonal",
  "tumour_shared"        = "Tumour shared",
  "tumour_private"       = "Tumour private"
)

#####################################
#### Create different plot types ####
#####################################

## Get data for bar plots

# pivot num samples -- rework this code to have post copy number QC numbers (refer to trees)
#num_samples_long <- cohort %>%
#  select(patient_tumour, patient_order, plt_cohort_cancer, num_tumour_samples, num_preinv_samples) %>%
#  pivot_longer(cols = starts_with("num_"),
#               names_to = "sample_type",
#               values_to = "count") %>%
#  mutate(sample_type = recode(sample_type,
#                              num_tumour_samples = "Tumour",
#                              num_preinv_samples = "Pre-invasive"))

# pivot num samples
num_samples_long <- mutations_labelled %>%
  select(patient_id, patient_tumour, pre_sample_count, cancer_sample_count) %>%
  group_by(patient_tumour) %>%
  slice(1) %>%
  left_join(cohort %>% select(patient_tumour, plt_cohort_cancer, patient_order), by = "patient_tumour") %>%
  pivot_longer(cols = ends_with("_count"),
               names_to = "sample_type",
               values_to = "count") %>%
  mutate(sample_type = recode(sample_type,
                              cancer_sample_count = "Tumour",
                              pre_sample_count = "Pre-invasive"))

# pivot mutation category data
mutation_long <- mutation_proportions %>%
  left_join(cohort %>% select(patient_tumour, patient_order, plt_cohort_cancer), by = "patient_tumour") %>%
  pivot_longer(cols = -c(patient_id, patient_tumour, patient_order, plt_cohort_cancer, total),
               names_to = "mutation_category",
               values_to = "proportion") %>%
  filter(!mutation_category %in% c("total"))

# pivot copy number data 
copy_number_long <- scna_metrics %>%
  pivot_longer(cols = c(frac_abberant_genom_subcl, frac_abberant_genom_clonal), names_to = "cna_category", values_to = "proportion") %>%
  mutate(cna_category = recode(cna_category, "frac_abberant_genom_subcl" = "Subclonal", "frac_abberant_genom_clonal" = "Clonal")) %>%
  left_join(cohort %>% select(patient_tumour, patient_order, plt_cohort_cancer), by = "patient_tumour")

# pivot wgd data
wgd_long <- wgd_metrics %>%
  pivot_longer(cols = c(num_clonal_gds, num_subclonal_gds), names_to = "category", values_to = "count") %>%
  mutate(category = recode(category, "num_clonal_gds" = "Clonal", "num_subclonal_gds" = "Subclonal")) %>%
  left_join(cohort %>% select(patient_tumour, patient_order, plt_cohort_cancer), by = "patient_tumour") %>%
  mutate(Patient = sub("_.*", "", patient_tumour))


## Get data for tile plots

tile_data <- cohort %>%
  select(patient_order, plt_cohort_cancer, plt_stage, plt_smoking, plt_gender, seeding_pattern) %>%
  pivot_longer(cols = c(plt_stage, plt_smoking, plt_gender, seeding_pattern),
               names_to = "attribute",
               values_to = "value") %>%
  mutate(attribute = recode(attribute,
                            plt_stage = "Stage",
                            plt_gender = "Gender",
                            plt_smoking = "Smoking Status",
                            seeding_pattern = "Seeding Pattern"))


## Make tile plots

tile_stage <- ggplot(tile_data %>% filter(attribute == "Stage")) +
  geom_tile(aes(x = factor(patient_order), y = "", fill = value), color = "black") +
  labs(y = "Stage") +
  scale_fill_manual(name = "Stage", values = stage_colours) +
  facet_wrap(~plt_cohort_cancer, scales = "free_x", nrow = 1) +
  theme_minimal() +
  blank_theme

tile_smoking <- ggplot(tile_data %>% filter(attribute == "Smoking Status")) +
  geom_tile(aes(x = factor(patient_order), y = "", fill = value), color = "black") +
  scale_fill_manual(name = "Smoking Status", values = smoking_colours) +
  labs(y = "Smoking Status") +
  facet_wrap(~plt_cohort_cancer, scales = "free_x", nrow = 1) +
  theme_minimal() +
  blank_theme

tile_gender <- ggplot(tile_data %>% filter(attribute == "Gender")) +
  geom_tile(aes(x = factor(patient_order), y = "", fill = value), color = "black") +
  scale_fill_manual(name = "Gender", values = gender_colours) +
  labs(y = "Gender") +
  facet_wrap(~plt_cohort_cancer, scales = "free_x", nrow = 1) +
  theme_minimal() +
  blank_theme

tile_seeding <- ggplot(tile_data %>% filter(attribute == "Seeding Pattern")) +
  geom_tile(aes(x = factor(patient_order), y = "", fill = value), color = "black") +
  scale_fill_manual(name = "Initiation Pattern", values = seeding_pattern_colours) +
  labs(y = "Initiation Pattern") +
  facet_wrap(~plt_cohort_cancer, scales = "free_x", nrow = 1) +
  theme_minimal() +
  blank_theme


## Make bar plots

# number of samples bar plot
num_sample_bars <- ggplot(num_samples_long, aes(x = factor(patient_order), y = count, fill = sample_type)) +
  geom_bar(stat = "identity", width = 0.9) +
  facet_wrap(~plt_cohort_cancer, scales = "free_x", nrow = 1) +
  theme_minimal()+
  scale_fill_manual(values = sample_type_colors, name = "Sample Type") +
  labs(y = "Num Regions", x = NULL) +
  bar_theme

# mutation bar plot
mutation_bar_plot <- ggplot(mutation_long, aes(x = factor(patient_order), y = proportion, fill = mutation_category)) +
  geom_bar(stat = "identity", width = 0.9) +
  facet_wrap(~plt_cohort_cancer, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = mutation_category_colors, labels = mutation_category_labels, name = "Mutation Category") +
  labs(y = "Proportion of Muts", x = NULL) +
  theme_minimal() +
  bar_theme

# copy number barplot
copy_number_bar <- ggplot(copy_number_long, aes(x = factor(patient_order), y = proportion, fill = cna_category)) +
  geom_bar(stat = "identity", width = 0.9) +
  facet_wrap(~plt_cohort_cancer, scales = "free_x", nrow = 1) +
  scale_fill_manual(name = 'Clonality', values = clonality_colours) +
  labs(y = "Proportion Abberant\nGenome SCNA", x = NULL) +
  theme_minimal() +
  bar_theme

# wgd bar plot
patient_labels <- setNames(wgd_long$Patient, wgd_long$patient_order)

gds_bar <- ggplot(wgd_long, aes(x = factor(patient_order), y = count, fill = category)) +
  geom_bar(stat = "identity", width = 0.9) +
  facet_wrap(~plt_cohort_cancer, scales = "free_x", nrow = 1) +
  scale_x_discrete(labels = patient_labels) +
  scale_fill_manual(values = clonality_colours, name = "Clonality") +
  labs(y = "Num WGD", x = NULL) +
  theme_minimal() +
  bar_theme 
  
# Extract legends
legend_num_samples <- get_legend(num_sample_bars)
legend_mutation <- get_legend(mutation_bar_plot)
legend_copy_number <- get_legend(copy_number_bar)
legend_stage <- get_legend(tile_stage)
legend_smoking <- get_legend(tile_smoking)
legend_gender <- get_legend(tile_gender)
legend_seeding <- get_legend(tile_seeding)


# Final plot
final_plot <- plot_grid(
  num_sample_bars + theme(legend.position = "none"),
  tile_stage + theme(legend.position = "none"),
  tile_smoking + theme(legend.position = "none"),
  tile_gender + theme(legend.position = "none"),
  tile_seeding + theme(legend.position = "none"),
  mutation_bar_plot + theme(legend.position = "none"),
  copy_number_bar + theme(legend.position = "none"),
  gds_bar + theme(legend.position = "none", axis.text.x = element_text(size = 10, angle = 90)),
  ncol = 1,
  align = "v",
  axis = "tb",
  rel_heights = c(2, 0.3, 0.3, 0.3, 0.3, 1, 1, 1)  
)

# Legend plot
legends <- plot_grid(legend_num_samples,
                     legend_stage,
                     legend_smoking,
                     legend_gender,
                     legend_seeding,
                     legend_mutation,
                     legend_copy_number,
                     ncol = 7,
                     align = "h",
                     axis = "t",
                     rel_widths = c(0.9, 0.8, 0.8, 0.8, 1.2, 1.2, 1))

legends <- plot_grid(legend_num_samples,
                     legend_stage,
                     legend_smoking,
                     legend_gender,
                     legend_seeding,
                     legend_mutation,
                     legend_copy_number,
                     ncol = 1,
                     align = "v",
                     rel_widths = c(0.7, 1.3, 1, 0.1, 1.2, 1.4, 0.5))
ggsave(paste0(outputs.folder, date, "_overview_legend_vertical.png"), width = 2, height = 12, dpi = 300)


plot_and_leegend <- plot_grid(final_plot, legends, ncol = 1, align = "v", rel_heights = c(2, 0.6))

ggsave(paste0(outputs.folder, date, "_ESCC_cohort_overview.png"), width = 9, height = 12, dpi = 300)


######################################
#### Look at drivers for cohort #####
######################################

drivers <- lapply(names(trees), function(tumour){
  tumour_name <- sub("_(\\d+)$", "_tumour\\1", tumour)
  drivers <- mutations %>% 
    filter(patient_tumour == tumour_name, is_driver_mut, is_present, tree_cluster) %>%
    select(patient_tumour, sample_name_hash, mutation_id, mutation_cluster, gene, aachange_refgene, is_trunk)
  return(drivers)
})

names(drivers) <- names(trees)


# how many clonal driver mutations
trunk_counts <- sapply(drivers, function(df) {
  sum(!is.na(df$is_trunk) & df$is_trunk == TRUE & !duplicated(df$mutation_id))
})

nontrunk_counts <- sapply(drivers, function(df) {
  sum(!is.na(df$is_trunk) & df$is_trunk == FALSE & !duplicated(df$mutation_id))
})

# count average number of drivers
patient_mut_counts <- lapply(names(drivers), function(pat) {
  df <- drivers[[pat]]
  if (nrow(df) > 0){
  aggregate(mutation_id ~ patient_tumour, data = df, FUN = function(x) length(unique(x)))
  } else {
    data.frame(patient_tumour = pat, mutation_id = 0)
  }
})

# combine all the results into one data frame
all_counts <- do.call(rbind, patient_mut_counts)

average_unique_drivers <- mean(all_counts$mutation_id)


# most frequent driver genes
gene_mutations <- lapply(drivers, function(df) {
  unique(df[, c("mutation_id", "gene")])
})

all_gene_mutations <- do.call(rbind, gene_mutations)
gene_counts <- aggregate(mutation_id ~ gene, data = all_gene_mutations, FUN = length)
gene_counts_ordered <- gene_counts[order(gene_counts$mutation_id, decreasing = TRUE), ]

print(gene_counts_ordered)


######################################
####### Deep dive on EPA00020 ########
######################################

EPA00020_muts <- read.fst("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/inputs/_RELEASE/EPA00020/mutation_table/tumour_1/EPA00020_tumour1_mutation_table_concise.fst")













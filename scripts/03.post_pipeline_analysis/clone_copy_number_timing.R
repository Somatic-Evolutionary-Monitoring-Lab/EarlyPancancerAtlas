#==============================================================================#
#==============================================================================#
######                                                                    ######
######           Investigating timing of copy number acquisition          ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Katherine Honan
# Date: 2026-05-11

setwd("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas")

#============================================#
# Source required functions & load libraries #
#============================================#

library(fst)
library(readr)
library(tibble)
library(tidyr)
library(stringr)
library(data.table) 
library(dplyr) 
library(ggplot2) 
library(cowplot)
library(RColorBrewer) 
library(patchwork)

#=====================================#
# Make a folder for this analysis run #
#=====================================#

date <- gsub("-","",Sys.Date())

analysis_name <- 'copy_number_timing'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep='/')
if( !file.exists(out_dir_general) ) dir.create( out_dir_general )

out_dir_logs <- paste(out_dir_general, 'logs', sep='/')
if( !file.exists(out_dir_logs) ) dir.create( out_dir_logs )

outputs.folder <- paste0( out_dir_general, "/", date, "/" )

if( !file.exists(outputs.folder) ) dir.create( outputs.folder )


#======================================#
# Get Inputs required for all analyses #
#======================================#

# temporarily read in the inputs from the extra directory I accidentally made
trees <- readRDS("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/_RELEASE/rds/project/rds-LH0AvU65IRI/_RELEASE/release-20260306/_aggregate/conipher_trees/2026_05_02_cohort_conipher_trees.RDS")
inventory <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.tsv")
cn_change_to_ancestor_paths <- "/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/_RELEASE/rds/project/rds-LH0AvU65IRI/_RELEASE/release-20260306/*/EPA*/alpaca/results/cn_change_to_ancestor.csv"
seeding <- readRDS("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/outputs/preinvasive_seeding/20260501/20260501_seeding_clonality.rds")
clone_info <- readRDS("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/outputs/preinvasive_seeding/20260501/20260501_cloneInfo.list.rds")
# get files
files <- Sys.glob(cn_change_to_ancestor_paths)

# extract patient ids
patient_ids <- sub(".*?/(EPA[^/]+)/alpaca/results/cn_change_to_ancestor\\.csv$", "\\1", files)

# store each patients csv as list element
cn_change_to_ancestor_list <- setNames(
  lapply(files, read.csv),
  patient_ids
)

# read in data on RDS
#trees <- readRDS("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/_RELEASE/release-20260306/_aggregate/conipher_trees/2026_03_06_cohort_conipher_trees.RDS")
#inventory <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.tsv")
#metadata <- read_tsv("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/map_files/20250514_chen_2017_ESCC_metadata_EPA_mapping.txt")

# source both original (tissue) mets functions and ctdna adapted functions
source("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/scripts/00.useful_functions/mets_functions.R")
source("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/scripts/00.useful_functions/useful_functions.R")
source("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/scripts/04.figures/EPA_colour_palettes.R")

#===============================#
# EPA ids & data pre-processing #
#===============================#

cohort_overview <- inventory %>%
  filter(!sample_type_simple %in% c("unknown_histology", "cancer_met"))

names(trees) <- sapply(trees, function(x) x$parameters$sampleID)

#===========#
# Functions #
#===========#

#### Functions from https://github.com/McGranahanLab/ALPACA-paper/blob/master/bin/ANALYSIS/analysis_functions.R#L703

### Function to compute interval events between a pair of clones based on
# their cross genome copy number profile
get_interval_events <- function(parent_cn, child_cn, change_threshold_cap = NA) {
  cn_change <- child_cn - parent_cn
  
  if (!is.na(change_threshold_cap)) {
    cn_change[cn_change > change_threshold_cap] <- change_threshold_cap
    cn_change[cn_change < -change_threshold_cap] <- -change_threshold_cap
  }
  
  # count interval events
  interval_events <- c(cn_change[1])
  
  if (length(cn_change) > 1) {
    for (seg_change_id in 2:length(cn_change)) {
      seg_change <- cn_change[seg_change_id]
      if (interval_events[length(interval_events)] != seg_change) {
        interval_events <- c(interval_events, seg_change)
      }
    }
  }
  
  gains <- interval_events[interval_events > 0]
  losses <- interval_events[interval_events < 0] * -1
  
  # add binary events:
  gains_bin <- as.numeric(gains > 0)
  losses_bin <- as.numeric(losses > 0)
  
  interval_events_df <- data.frame(interval_gains = sum(gains),
                                   interval_losses = sum(losses),
                                   interval_gains_bin = sum(gains_bin),
                                   interval_losses_bin = sum(losses_bin))
  
  interval_events_df$interval_events = interval_events_df$interval_gains + interval_events_df$interval_losses
  interval_events_df$interval_events_bin = interval_events_df$interval_gains_bin + interval_events_df$interval_losses_bin
  
  return(interval_events_df)
}

### Function to compute interval events
compute_interval_events_alpaca <- function(edge_df, alpaca_tum_out) {  
  alpaca_clone_cn <- unique(alpaca_tum_out[, .(segment, clone = gsub('clone', '', clone), pred_CN_A, pred_CN_B)])
  alpaca_clone_cn[, chr := mapply(function(x) as.numeric(strsplit(x, "_")[[1]][1]), segment)]
  alpaca_clone_cn[, start := mapply(function(x) as.numeric(strsplit(x, "_")[[1]][2]), segment)]
  
  setorder(alpaca_clone_cn, chr, start)
  all_clones <- unique(edge_df[, clone])
  chromosomes <- sort(unique(as.numeric(alpaca_clone_cn$chr)))
  chr_events <- rbindlist(lapply(chromosomes, function(chromosome) {
    int_events <- rbindlist(lapply(all_clones, function(child) {
      parent <- edge_df[clone == child, parent_clone]
      child_cn_A <- alpaca_clone_cn[chr == chromosome & clone == child, pred_CN_A]
      child_cn_B <- alpaca_clone_cn[chr == chromosome & clone == child, pred_CN_B]
      parent_cn_A <- alpaca_clone_cn[chr == chromosome & clone == parent, pred_CN_A]
      parent_cn_B <- alpaca_clone_cn[chr == chromosome & clone == parent, pred_CN_B]
      interval_events_A <- get_interval_events(parent_cn = parent_cn_A, child_cn = child_cn_A)
      interval_events_B <- get_interval_events(parent_cn = parent_cn_B, child_cn = child_cn_B)
      interval_events_A$allele <- "A"
      interval_events_B$allele <- "B"
      interval_events_df <- rbind(interval_events_A, interval_events_B)
      # sum interval events across alleles:
      cols_to_sum <- colnames(interval_events_df)[colnames(interval_events_df) != "allele"]
      total_events <- paste0("total_", cols_to_sum)
      interval_events_df <- as.data.table(interval_events_df)
      interval_events_df[, (total_events) := lapply(.SD, sum), .SDcols = cols_to_sum]
      interval_events_df[, clone := child]
      interval_events_df[, chr := chromosome]
      return(interval_events_df)
    }))
    return(int_events)
  }))
  # Now sum across chromosomes
  cols_to_keep <- c("allele", "clone", "chr")
  cols_to_sum <- colnames(chr_events)[!colnames(chr_events) %in% cols_to_keep]
  fullgenome_events <- copy(chr_events)
  for (c in cols_to_sum) {
    fullgenome_events[, (c) := sum(get(c)), by = .(clone, allele)]
  }
  final_cols <- c(cols_to_sum, "clone", "allele")
  fullgenome_events <- unique(fullgenome_events[, ..final_cols])
  return(fullgenome_events)
}


#==========================#
# Define timing categories #
#==========================#

# classify clones into timing categories
clone_timing_categories <- lapply(names(clone_info), function(pat){
  # MRCA, shared, seeding, primary specific and metastasis specific
  clone_timing <- clone_info[[pat]]$df.cloneInfo %>%
    mutate(
      timing_category = case_when(
        clonalClones & !seedingClones ~ "MRCA Non-initiating",
        seedingClones & clonalClones ~ "MRCA Initiating",
        seedingClones ~ "Initiating",
        sharedClones & !seedingClones ~ "Shared",
        preinvasiveClones ~ "Preinvasive private",
        tumourClones ~ "Tumour private",
        TRUE ~ NA_character_
      )
    ) %>%
    select(
      patient_tumour = Patient,
      clone = clones,
      timing_category
    ) %>%
    filter(!is.na(timing_category)) # remove clones from met samples
    
})
clone_timing_categories_df <- bind_rows(clone_timing_categories)

# want a dataframe like:
# patient, clone, timing category, num copy gain, num copy loss
# clone timing categories:
# clonal, preinvasive private, prinvasive shared, seeding, tumour shared, tumour private
# do we want to compare these sets: [preinvasive private + shared] vs [seeding]?
# preinvasive unique and shared together covers non-seeding preinvasive clones


# To investigate timing of copy number acquisition across the premalignant
# switch, quantify allele specific copy number events within timing bins

# per-segment change and total number of interval events
# Example patient:
pat = "EPA00015_2"
tree <- trees[[pat]]$graph_pyclone$default_tree
trunk <- trees[[pat]]$graph_pyclone$trunk
seeding <- clone_info[[pat]]$seedingClones
edgelength <- trees[[pat]]$graph_pyclone$edgelength
tree_clones <- unique(c(tree[,1], tree[,2]))

# color assignment
preinvasive_unique <- clone_info[[pat]]$preinvasiveClones[clone_info[[pat]]$preinvasiveClones %in% tree_clones]
cancer_unique <- clone_info[[pat]]$tumourClones[clone_info[[pat]]$tumourClones %in% tree_clones]
shared_nodes <- clone_info[[pat]]$sharedClones[clone_info[[pat]]$sharedClones %in% tree_clones]

# set node colours
node_colors <- c(setNames(rep("#D55E00", length(preinvasive_unique)), preinvasive_unique),
                 setNames(rep("#009E73", length(cancer_unique)), cancer_unique),
                 setNames(rep("#D0DDE2", length(shared_nodes)), shared_nodes))


# plot tree
plottingTxTree(
  tree, 
  trunk, 
  edgelength,
  seedingCluster = seeding,
  color.mapping = "perNode",
  col.palette = node_colors,
  tumorID = pat,
  add.title = TRUE,
  node.size = 15,
  node.label.type = "nodeNumber"
)

# Compute number of SCNAs on tree edges for each patient
# Described in https://www.nature.com/articles/s41586-025-09398-w#Sec10
# under: "Computing number of SCNAs on a tree edge"

events_along_tree_edge <- lapply(names(cn_change_to_ancestor_list), function(pat){
  print(pat)
  cn_change <- data.table(cn_change_to_ancestor_list[[pat]])
  # Explicitly add in the diploid baseline to cn table
  diploid_cn <- unique(cn_change[, .(tumour_id, segment)])
  diploid_cn[, clone := "diploid"]
  diploid_cn[, pred_CN_A := 1L]
  diploid_cn[, pred_CN_B := 1L]
  # columns required but not meaningful for diploid root
  diploid_cn[, parent := NA_character_]
  diploid_cn[, parent_pred_cpnA := NA_integer_]
  diploid_cn[, parent_pred_cpnB := NA_integer_]
  diploid_cn[, cn_dist_to_parent_A := NA_integer_]
  diploid_cn[, cn_dist_to_parent_B := NA_integer_]
  # reorder columns to match original table
  setcolorder(diploid_cn, colnames(cn_change))
  # append
  cn_change <- rbind(cn_change, diploid_cn)

  tree_structure <- trees[[pat]]$graph_pyclone$default_tree
  trunk <- trees[[pat]]$graph_pyclone$trunk
  tree_structure <- rbind(tree_structure, c("diploid", trunk))
  tree_df <- data.table(tree_structure)
  colnames(tree_df) <- c("parent_clone", "clone")

  # 1. Compute number of SCNAs on a tree edge, i.e. actual number of copy-number 
  # events that occurred in the copy-number evolutionary tree inferred by ALPACA 
  # while taking into account neighboring SCNAs
  edge_events <- compute_interval_events_alpaca(
    edge_df = tree_df, 
    alpaca_tum_out = cn_change
    ) 
  edge_events[, patient_tumour := pat]
  
  # pivot allele-specific metrics wider
  edge_events_wide <- dcast(
    edge_events,
    patient_tumour + clone ~ allele,
    value.var = c(
      "interval_gains",
      "interval_losses",
      "interval_gains_bin",
      "interval_losses_bin",
      "interval_events",
      "interval_events_bin"
    )
  )
  # get total metrics per clone
  edge_totals <- unique(
    edge_events[, .(
      patient_tumour,
      clone,
      total_interval_gains,
      total_interval_losses,
      total_interval_gains_bin,
      total_interval_losses_bin,
      total_interval_events,
      total_interval_events_bin
    )]
  )
  # add back
  edge_events_wide <- merge(
    edge_events_wide,
    edge_totals,
    by = c("patient_tumour", "clone")
  )
})
events_along_tree_edge_df <- bind_rows(events_along_tree_edge)

# join with clone timing df
clone_timing_events <- events_along_tree_edge_df %>%
    left_join(clone_timing_categories_df, by = c("patient_tumour", "clone"))
    

# Figure 1. Box plot comparing the number of SCNAs in seeding versus non-initiating primary clones
# y axis = num events, 
# x axis = sample type (initiating clone vs on initiating clone)
# 3 panels, one for all events, one for gains, one for losses
# wilcoxon test

# if multiple seeding clones along a branch -> count only top level ancestral clone?

# Collapse timing categories
plot_1_df <- clone_timing_events %>%
  filter(timing_category %in% c("Preinvasive private", "Shared", "Initiating")) %>%
  mutate(clone_class = ifelse(timing_category == "Initiating", "Initiating", "Non-initiating"))

# Create long-format table
plot_1_long <- bind_rows(
  plot_1_df %>%
    transmute(patient_tumour, clone, clone_class, event_type = "All events", value = total_interval_events_bin),
  plot_1_df %>%
    transmute(patient_tumour, clone, clone_class, event_type = "Gains", value = total_interval_gains_bin),
  plot_1_df %>% transmute(patient_tumour, clone, clone_class, event_type = "Losses", value = total_interval_losses_bin)
  ) %>%
  mutate(event_type = factor(event_type, levels = c("All events", "Gains", "Losses")),
         clone_class = factor(clone_class, levels = c("Non-initiating", "Initiating")))

# Wilcoxon p-values
pvals <- plot_1_long %>%
  group_by(event_type) %>%
  summarise(
    p = wilcox.test(
      value[clone_class == "Non-initiating"],
      value[clone_class == "Initiating"]
    )$p.value
  ) %>%
  mutate(
    label = paste0("P = ", signif(p, 2))
  )

# Plot
ggplot(plot_1_long, aes(x = clone_class, y = value, fill = clone_class)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  facet_wrap(~event_type, scales = "free_y") +
  geom_text(
    data = pvals, 
    aes(x = 1.5, y = Inf, label = label), 
    inherit.aes = FALSE, 
    vjust = 1.5, 
    size = 4
    ) +
  scale_fill_manual(values = initiating_clone_colours) +
  labs(x = "Clone class", y = "No. SCNAs") +
  theme_classic() +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )
ggsave(paste0(outputs.folder, date, "_init_vs_non_init_clone_events.png"), width = 9, height = 7)

# Figure 2. Fraction of seeding (purple) and non-seeding (green) clones with 
# gain (top) or loss (bottom) at each genomic locus ()
# https://www.nature.com/articles/s41586-025-09398-w/figures/13


# 2. Compute proportion of genome altered relative to parent in each clone 
# weighted by segment size 


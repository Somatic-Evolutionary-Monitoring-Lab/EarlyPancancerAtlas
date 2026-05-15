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
library(rstatix)
library(ggpubr)
library(GenomicRanges)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg19)
library(RColorBrewer) 
library(scales)
library(purrr)
library(igraph)

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
clone_wgd_ratio_paths <- "/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/inputs/_RELEASE/rds/project/rds-LH0AvU65IRI/_RELEASE/release-20260306/*/EPA*/alpaca/results/wgd_ratio_scores.csv"
seeding <- readRDS("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/outputs/preinvasive_seeding/20260501/20260501_seeding_clonality.rds")
clone_info <- readRDS("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/outputs/preinvasive_seeding/20260501/20260501_cloneInfo.list.rds")

# get cn change to ancestor files
cn_change_to_ancestor_files <- Sys.glob(cn_change_to_ancestor_paths)

# extract patient ids
patient_ids <- sub(".*?/(EPA[^/]+)/alpaca/results/cn_change_to_ancestor\\.csv$", "\\1", cn_change_to_ancestor_files)

# store each patients csv as list element
cn_change_to_ancestor_list <- setNames(
  lapply(cn_change_to_ancestor_files, read.csv),
  patient_ids
)

# get cn change to ancestor files
clone_wgd_ratio_files <- Sys.glob(clone_wgd_ratio_paths)

# extract patient ids
patient_ids <- sub(".*?/(EPA[^/]+)/alpaca/results/wgd_ratio_scores\\.csv$", "\\1", clone_wgd_ratio_files)

# store each patients csv as list element
clone_wgd_ratio_list <- setNames(
  lapply(clone_wgd_ratio_files, read.csv),
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
  # MRCA, shared, seeding, preinvasive specific and primary specific
  clone_timing <- clone_info[[pat]]$df.cloneInfo %>%
    mutate(
      timing_category = case_when(
        clonalClones & !seedingClones ~ "MRCA Non-initiating",
        seedingClones & clonalClones ~ "MRCA Initiating",
        seedingClones & !clonalClones ~ "Subclonal initiating",
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


# per-segment change and total number of interval events
# Example patient:
pat = "EPA00020_1"
tree <- trees[[pat]]$graph_pyclone$default_tree
trunk <- trees[[pat]]$graph_pyclone$trunk
seeding_clones <- clone_info[[pat]]$seedingClones
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
  seedingCluster = seeding_clones,
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
  # add patient tumour id
  edge_events[, patient_tumour := pat]
  # add back parent
  edge_events[tree_df, parent := i.parent_clone, on = .(clone)]
  
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

clone_timing_events <- clone_timing_events %>%
  group_by(patient_tumour) %>%
  mutate(
    prop_scna_all = total_interval_events_bin/sum(total_interval_events_bin),
    prop_gains = total_interval_gains_bin/sum(total_interval_gains_bin),
    prop_losses = total_interval_losses_bin/sum(total_interval_losses_bin),
  ) %>%
  ungroup()

write_delim(clone_timing_events, paste0(outputs.folder, date, "_scna_timing_per_clone.tsv"), delim = "\t")

# collapse across categories
clone_timing_categories_counts <- clone_timing_events %>%
  group_by(patient_tumour, timing_category) %>%
  summarise(
    total_num_scna = sum(total_interval_events_bin),
    num_gain = sum(total_interval_gains_bin),
    num_loss = sum(total_interval_losses_bin),
    prop_num_scna = sum(prop_scna_all),
    prop_gain = sum(prop_gains),
    prop_loss = sum(prop_losses),
    .groups = "drop"
  )

write_delim(clone_timing_categories_counts, paste0(outputs.folder, date, "_scna_count_per_timing.tsv"), delim = "\t")

# Figure 1a. Box plot comparing the number of SCNAs in seeding versus non-initiating primary clones

# if multiple seeding clones along a branch -> count only top level ancestral clone?

# Collapse timing categories
plot_1_df <- clone_timing_events %>%
  filter(timing_category %in% c("Preinvasive private", "Subclonal initiating")) %>%
  mutate(clone_category = ifelse(timing_category == "Subclonal initiating", "Subclonal initiating", "Non-initiating"))

# Create long-format table
plot_1_long <- bind_rows(
  plot_1_df %>%
    transmute(patient_tumour, clone, clone_category, event_type = "All events", value = total_interval_events_bin),
  plot_1_df %>%
    transmute(patient_tumour, clone, clone_category, event_type = "Gains", value = total_interval_gains_bin),
  plot_1_df %>% transmute(patient_tumour, clone, clone_category, event_type = "Losses", value = total_interval_losses_bin)
  ) %>%
  mutate(event_type = factor(event_type, levels = c("All events", "Gains", "Losses")),
         clone_category = factor(clone_category, levels = c("Non-initiating", "Subclonal initiating")))

# Wilcoxon p-values
pvals <- plot_1_long %>%
  group_by(event_type) %>%
  summarise(
    p = wilcox.test(
      value[clone_category == "Non-initiating"],
      value[clone_category == "Subclonal initiating"]
    )$p.value
  ) %>%
  mutate(
    label = paste0("P = ", signif(p, 2))
  )

# Plot
ggplot(plot_1_long, aes(x = clone_category, y = value, fill = clone_category)) +
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
  scale_fill_manual(values = clone_category_colours) +
  labs(x = "Clone Category", y = "Num. SCNAs") +
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



# Figure 1b. Box plot comparing the number of SCNAs all timing categories

# Set timing category order
timing_levels <- unique(na.omit(clone_timing_events$timing_category))
timing_levels <- c(
  "MRCA Non-initiating",
  "MRCA Initiating",
  "Shared",
  "Subclonal initiating",
  "Preinvasive private",
  "Tumour private"
  )

# Prepare plotting df
plot_2_df <- clone_timing_events %>%
  filter(timing_category %in% timing_levels) %>%
  mutate(
    timing_category = factor(
      timing_category,
      levels = timing_levels
    )
  )

# Long format
plot_2_long <- bind_rows(
  plot_2_df %>%
    transmute(patient_tumour, clone, timing_category, event_type = "All events", value = total_interval_events_bin ),
  plot_2_df %>%
    transmute(patient_tumour, clone, timing_category, event_type = "Gains", value = total_interval_gains_bin),
  plot_2_df %>%
    transmute(patient_tumour, clone, timing_category, event_type = "Losses", value = total_interval_losses_bin)
  ) %>%
  mutate(event_type = factor(event_type,levels = c("All events", "Gains", "Losses")))

pairwise_tests <- plot_2_long %>%
  group_by(event_type) %>%
  pairwise_wilcox_test(value ~ timing_category, p.adjust.method = "BH")

ggplot(plot_2_long, aes(x = timing_category, y = value, fill = timing_category)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  facet_wrap(~event_type, scales = "free_y") +
  labs(x = "Clone timing category", y = "Num. SCNAs") +
  scale_fill_manual(values = clone_category_colours) +
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

ggsave(paste0(outputs.folder, date, "_all_categories_clone_events.png"), width = 12, height = 7)

#=========================#
# SCNA vs initiation type #
#=========================#

## plot total number of SCNA events across all tree edges split by mono vs polyclonal

total_scna_per_pat <- events_along_tree_edge_df %>%
  group_by(patient_tumour) %>%
  summarise(
    total_events = sum(total_interval_events_bin),
    total_gains = sum(total_interval_gains_bin),
    total_losses = sum(total_interval_losses_bin)
    )

total_scna_vs_clonality <- total_scna_per_pat %>%
  left_join(seeding, by = c("patient_tumour" = "tumour_id")) %>%
  filter(!is.na(seeding_clonality)) %>%
  mutate(
    initiation_class = 
      ifelse(
        seeding_clonality == "monoclonal", "monoclonal", paste0(seeding_clonality, " ", seeding_phyletic)
        )
    )

plot_df_long <- total_scna_vs_clonality %>%
  pivot_longer(
    cols = c(total_events, total_gains, total_losses),
    names_to = "event_type",
    values_to = "num_events"
  )


ggplot(plot_df_long, aes(x = seeding_clonality, y = num_events, fill = seeding_clonality)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("monoclonal", "polyclonal")),
    label = "p.format"
  ) +
  facet_wrap(~ event_type, scales = "free_y") +
  labs(x = "Initiation class", y = "Avg. SCNA across tree edges") +
  scale_fill_manual(values = c("monoclonal" = "#00A087", "polyclonal" = "#3C5488")) +
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

ggsave(paste0(outputs.folder, date, "_mean_scna_per_edge_vs_seeding.png"), width = 7, height = 5)


## plot number of events across tree edge at each depth vs mono/polyclonal

## work out clone tree depths

tree_depth_df <- imap_dfr(
  trees, function(tree_obj, tumour_id) {
    edge_df <- as.data.frame(tree_obj$graph_pyclone$default_tree, stringsAsFactors = FALSE)
    colnames(edge_df) <- c("parent", "child")
    trunk <- tree_obj$graph_pyclone$trunk
    
    # build directed graph
    g <- graph_from_data_frame(edge_df, directed = TRUE)
    
    # shortest path from trunk to every node
    depths <- distances(g, v = trunk, to = V(g), mode = "out")
    
    tibble(
      patient_tumour = tumour_id,
      clone = names(depths[1, ]),
      tree_depth = as.numeric(depths[1, ]) + 1
    )
  }
)

tree_depth_events <- events_along_tree_edge_df %>%
  left_join(
    tree_depth_df,
    by = c("patient_tumour", "clone")
  ) %>%
  group_by(patient_tumour, tree_depth) %>%
  summarise(
    total_events = sum(total_interval_events_bin),
    total_gains = sum(total_interval_gains_bin),
    total_losses = sum(total_interval_losses_bin),
    .groups = "drop"
  )

total_scna_vs_clonality <- tree_depth_events %>%
  left_join(seeding, by = c("patient_tumour" = "tumour_id")) %>%
  filter(!is.na(seeding_clonality)) %>%
  mutate(
    initiation_class = 
      ifelse(
        seeding_clonality == "monoclonal", "monoclonal", paste0(seeding_clonality, " ", seeding_phyletic)
      )
  )

total_scna_vs_clonality <- tree_depth_events %>%
  left_join(
    seeding,
    by = c("patient_tumour" = "tumour_id")
  ) %>%
  filter(!is.na(seeding_clonality)) %>%
  mutate(
    initiation_class = ifelse(
      seeding_clonality == "monoclonal",
      "monoclonal",
      paste0(seeding_clonality, " ", seeding_phyletic)
    )
  )
plot_df <- total_scna_vs_clonality %>%
  mutate(
    tree_depth_group = case_when(
      tree_depth >= 5 ~ "5+",
      TRUE ~ as.character(tree_depth)
    ),
    tree_depth_group = factor(
      tree_depth_group,
      levels = c("1", "2", "3", "4", "5+")
    )
  )

# Wilcoxon test at each tree depth
wilcox_df <- plot_df %>%
  group_by(tree_depth_group) %>%
  summarise(
    p_value = wilcox.test(
      total_events ~ seeding_clonality
    )$p.value,
    
    # position annotation slightly above max value
    y_position = max(total_events, na.rm = TRUE) * 1.08,
    
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("p = ", signif(p_value, 2))
  )

ggplot(plot_df, aes(x = tree_depth_group, y = total_events, fill = seeding_clonality)) +
  geom_jitter(aes(colour = seeding_clonality), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8 )) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.shape = NA, alpha = 0.7) +
  geom_text(data = wilcox_df, aes(x = tree_depth_group, y = y_position, label = label), inherit.aes = FALSE, size = 3.5) +
  scale_fill_manual(
    values = c(
      "monoclonal" = "#00A087",
      "polyclonal" = "#3C5488"
    )
  ) +
  scale_colour_manual(
    values = c(
      "monoclonal" = "#00A087",
      "polyclonal" = "#3C5488"
    )
  ) +
  labs(
    x = "Tree depth",
    y = "Total SCNA events"
  ) +
  theme_classic() +
  theme(
    legend.position = "right"
  )

ggsave(paste0(outputs.folder, date, "_scna_by_tree_depth_vs_seeding.png"), width = 7, height = 5)

#======================#
# Genome-wide analysis #
#======================#

# Compute fraction of initiating and non-initiating clones with gain or loss at 
# each genomic locus
# Plot example: https://www.nature.com/articles/s41586-025-09398-w/figures/13

cn_change_to_ancestor_df <- rbindlist(cn_change_to_ancestor_list)


# Common segmentation approach
cn_change_pat_seg_clone <- cn_change_to_ancestor_df %>%
  separate(
    segment,
    into = c("chr", "start", "end"),
    sep = "_",
    convert = TRUE
  ) %>%
  select(tumour_id, chr, start, end, clone, cn_dist_to_parent_A, cn_dist_to_parent_B) %>%
  mutate(
    gain_event =
      cn_dist_to_parent_A > 0 |
      cn_dist_to_parent_B > 0,
    loss_event =
      cn_dist_to_parent_A < 0 |
      cn_dist_to_parent_B < 0,
    patient_tumour = str_replace(tumour_id, "-", "_")
  ) %>%
  left_join(
    clone_timing_categories_df %>% 
      mutate(clone = paste0("clone", clone)), 
    by = c("patient_tumour", "clone")
  )

# convert segments to granges
cn_gr <- GRanges(
  seqnames = cn_change_pat_seg_clone$chr,
  ranges = IRanges(
    cn_change_pat_seg_clone$start,
    cn_change_pat_seg_clone$end
  )
)

# get disjoint segmentation
common_gr <- disjoin(cn_gr)

# get table for clone X has gain/loss event in common interval Y
ov <- findOverlaps(common_gr, cn_gr)

mapped_df <- data.frame(
  segment_id = queryHits(ov),
  chr =
    as.character(seqnames(common_gr))[queryHits(ov)],
  start =
    start(common_gr)[queryHits(ov)],
  end =
    end(common_gr)[queryHits(ov)],
  patient_tumour =
    cn_change_pat_seg_clone$patient_tumour[subjectHits(ov)],
  clone =
    cn_change_pat_seg_clone$clone[subjectHits(ov)],
  timing_category =
    cn_change_pat_seg_clone$timing_category[subjectHits(ov)],
  gain_event =
    cn_change_pat_seg_clone$gain_event[subjectHits(ov)],
  loss_event =
    cn_change_pat_seg_clone$loss_event[subjectHits(ov)]
)

# now ever clone needs to be evaluated at every interval
all_clones <- cn_change_pat_seg_clone %>%
  distinct(
    patient_tumour,
    clone,
    timing_category
  )

all_segments <- data.frame(
  segment_id = seq_along(common_gr),
  chr = as.character(seqnames(common_gr)),
  start = start(common_gr),
  end = end(common_gr)
)

full_matrix <- crossing(
  all_clones,
  all_segments
)

full_matrix <- full_matrix %>%
  left_join(
    mapped_df,
    by = c(
      "patient_tumour",
      "clone",
      "timing_category",
      "segment_id",
      "chr",
      "start",
      "end"
    )
  ) %>%
  mutate(
    gain_event =
      ifelse(is.na(gain_event), FALSE, gain_event),
    loss_event =
      ifelse(is.na(loss_event), FALSE, loss_event),
    chr = paste0("chr", chr)
  )

full_matrix <- full_matrix %>%
  filter(timing_category %in% c("Preinvasive private", "Subclonal initiating")) %>%
  mutate(
    seeding_class =
      ifelse(
        timing_category == "Subclonal initiating",
        "Subclonal initiating",
        "Non-initiating"
      )
  )

freq_df <- full_matrix %>%
  group_by(
    segment_id,
    chr,
    start,
    end,
    seeding_class
  ) %>%
  summarise(
    pct_gain =
      mean(gain_event) * 100,
    pct_loss =
      mean(loss_event) * 100,
    .groups = "drop"
  )

# hg19
chr_lengths <- seqlengths(Hsapiens)

# cumulative offsets
chr_offsets <- cumsum(c(0, chr_lengths[-length(chr_lengths)]))
names(chr_offsets) <- names(chr_lengths)

freq_df <- freq_df %>%
  mutate(
    genomic_position =
      start +
      chr_offsets[as.character(chr)]
  )

# gain plot
ggplot(
  freq_df,
  aes(
    genomic_position,
    pct_gain,
    colour = seeding_class
  )
) +
  scale_colour_manual(values = clone_category_colours) +
  geom_line()+
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



# Binning approach

# prepare copy number events df: binary event labels, focus on preinv clones
cn_events <- cn_change_to_ancestor_df %>%
  separate(
    segment,
    into = c("chr", "start", "end"),
    sep = "_",
    convert = TRUE
  ) %>%
  mutate(
    patient_tumour = str_replace(tumour_id, "-", "_"),
    
    gain_event =
      cn_dist_to_parent_A > 0 |
      cn_dist_to_parent_B > 0,
    
    loss_event =
      cn_dist_to_parent_A < 0 |
      cn_dist_to_parent_B < 0,
    
    chr = paste0("chr", chr)
  ) %>%
  left_join(
    clone_timing_categories_df %>%
      mutate(clone = paste0("clone", clone)),
    by = c("patient_tumour", "clone")
  ) %>%
  filter(
    timing_category %in% c(
      "Subclonal initiating",
      "Preinvasive private"
    )
  ) %>%
  mutate(
    initiation_status = ifelse(
      timing_category == "Subclonal initiating",
      "Subclonal initiating",
      "Non-initiating"
    )
  )

# greate genomic bins
hg19_lengths <- seqlengths(Hsapiens)

hg19_lengths <- hg19_lengths[
  paste0("chr", c(1:22))
]

bins <- tileGenome(
  seqlengths = hg19_lengths,
  tilewidth = 5000,
  cut.last.tile.in.chrom = TRUE
)

bins_df <- data.frame(
  bin_id = seq_along(bins),
  chr = as.character(seqnames(bins)),
  start = start(bins),
  end = end(bins)
)

# make genomic ranges from copy numbver event df
gain_df <- cn_events %>%
  filter(gain_event)

gain_gr <- GRanges(
  seqnames = gain_df$chr,
  ranges = IRanges(
    start = gain_df$start,
    end = gain_df$end
  ),
  patient_tumour = gain_df$patient_tumour,
  clone = gain_df$clone,
  initiation_status = gain_df$initiation_status
)

loss_df <- cn_events %>%
  filter(loss_event)

loss_gr <- GRanges(
  seqnames = loss_df$chr,
  ranges = IRanges(
    start = loss_df$start,
    end = loss_df$end
  ),
  patient_tumour = loss_df$patient_tumour,
  clone = loss_df$clone,
  initiation_status = loss_df$initiation_status
)


# find overlaps with bins
gain_ov <- findOverlaps(bins, gain_gr)

gain_hits <- data.frame(
  bin_id = queryHits(gain_ov),
  patient_tumour = mcols(gain_gr)$patient_tumour[subjectHits(gain_ov)],
  clone = mcols(gain_gr)$clone[subjectHits(gain_ov)],
  initiation_status = mcols(gain_gr)$initiation_status[subjectHits(gain_ov)],
  gain = TRUE
) %>%
  distinct()

loss_ov <- findOverlaps(bins, loss_gr)

loss_hits <- data.frame(
  bin_id = queryHits(loss_ov),
  patient_tumour = mcols(loss_gr)$patient_tumour[subjectHits(loss_ov)],
  clone = mcols(loss_gr)$clone[subjectHits(loss_ov)],
  initiation_status = mcols(loss_gr)$initiation_status[subjectHits(loss_ov)],
  loss = TRUE
) %>%
  distinct()


# expand all clones/ bins into matrix
all_clones <- cn_events %>%
  distinct(
    patient_tumour,
    clone,
    initiation_status
  )

full_matrix <- crossing(
  all_clones,
  bin_id = bins_df$bin_id
)


# join with copy number info
full_matrix <- full_matrix %>%
  left_join(
    gain_hits,
    by = c(
      "patient_tumour",
      "clone",
      "initiation_status",
      "bin_id"
    )
  ) %>%
  left_join(
    loss_hits,
    by = c(
      "patient_tumour",
      "clone",
      "initiation_status",
      "bin_id"
    )
  ) %>%
  mutate(
    gain = ifelse(is.na(gain), FALSE, gain),
    loss = ifelse(is.na(loss), FALSE, loss)
  ) %>%
  left_join(
    bins_df,
    by = "bin_id"
  )

# calculate frequency of gains/losses in each bin
freq_df <- full_matrix %>%
  group_by(
    bin_id,
    chr,
    start,
    end,
    initiation_status
  ) %>%
  summarise(
    pct_gain = mean(gain) * 100,
    pct_loss = mean(loss) * 100,
    .groups = "drop"
  )

# losses plotted negative
freq_df <- freq_df %>%
  mutate(
    pct_loss = -pct_loss
  )

# create cumulative genomic coordinates
chr_df <- data.frame(
  chr = names(hg19_lengths),
  chr_length = as.numeric(hg19_lengths)
)

chr_df <- chr_df %>%
  mutate(
    offset = cumsum(lag(chr_length, default = 0))
  )

freq_df <- freq_df %>%
  left_join(chr_df, by = "chr") %>%
  mutate(
    genomic_position = start + offset
  )

# chromosome centers for labels
chr_centers <- chr_df %>%
  mutate(
    center = offset + (chr_length / 2)
  )


# plot
ggplot(freq_df) +
  geom_line(aes(genomic_position, pct_gain, colour = initiation_status), linewidth = 0.5) +
  geom_line(aes(genomic_position, pct_loss, colour = initiation_status), linewidth = 0.5) +
  geom_vline(xintercept = chr_df$offset, colour = "grey80", linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "black" ) +
  scale_colour_manual(values = clone_category_colours) +
  scale_x_continuous(breaks = chr_centers$center, labels = gsub("chr", "", chr_centers$chr)) +
  labs(x = "Chromosome", y = "% clones with event", colour = "Clone class") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 15), 
    panel.border = element_blank(),
    legend.position = "top"
    )

ggsave(paste0(outputs.folder, date, "_5kb_bin_pct_clones_altered_across_genome.png"), width = 16, height = 5)


#=============================#
# Next #
#=============================#

# Figure 2 Classification of tumours based on their initiation clone type

# 1. MRCA seeding
# 2. MRCA and subclones seed
# 2. only subclone seeding, no non-seeding paths
# 3. only subclone seeding, some non-seeding paths




#=============================#
# Clone level genome doubling #
#=============================#

#If the mean ratio, weighted by segment length, between the child clone and the 
#parent clone was above 1.5, we classified the child clone as genome-doubled. If 
#this ratio was below 1.5, but the ratio between the child and the grandparent was 
#above 1.5 and the parent clone was not already classified as genome-doubled, we 
#also classified the child clone as genome-doubled.

clone_wgd_status <- lapply(names(clone_wgd_ratio_list), function(pat){
  wgd_pat <- clone_wgd_ratio_list[[pat]]

  # initialise
  wgd_pat$wgd <- FALSE
  wgd_pat$cumulative_ratio <- NA_real_
  
  # root clones relative to diploid
  root_idx <- wgd_pat$parent == "diploid"
  
  wgd_pat$cumulative_ratio[root_idx] <- wgd_pat$ratio_A[root_idx]
  wgd_pat$wgd[root_idx] <- wgd_pat$ratio_A[root_idx] > 1.5
  
  # iterate through descendants
  repeat {
    
    updated <- FALSE
    
    for(i in seq_len(nrow(wgd_pat))) {
      
      if(!is.na(wgd_pat$cumulative_ratio[i]))
        next
      
      parent_row <- match(wgd_pat$parent[i], wgd_pat$clone)
      
      if(is.na(parent_row))
        next
      
      if(is.na(wgd_pat$cumulative_ratio[parent_row]))
        next
      
      # cumulative ratio to root
      wgd_pat$cumulative_ratio[i] <-
        wgd_pat$cumulative_ratio[parent_row] *
        wgd_pat$ratio_A[i]
      
      # direct parent test
      direct_wgd <- wgd_pat$ratio_A[i] > 1.5
      
      # grandparent rescue
      rescue_wgd <-
        !wgd_pat$wgd[parent_row] &&
        wgd_pat$cumulative_ratio[i] > 1.5
      
      wgd_pat$wgd[i] <- direct_wgd || rescue_wgd
      
      updated <- TRUE
    }
    
    if(!updated)
      break
  }
  return(wgd_pat)
})
clone_wgd_status_df <- rbindlist(clone_wgd_status)


# plot number of clones with genome doubling relative to parent split by timing category
clone_wgd_status_timing_df <- clone_wgd_status_df %>%
  mutate(patient_tumour = str_replace(tumour_id, "-", "_")) %>%
  left_join(
    clone_timing_categories_df %>% 
      mutate(clone = paste0("clone", clone)),
    by = c("patient_tumour", "clone")
    )

write_delim(clone_wgd_status_timing_df, paste0(outputs.folder, date, "_clone_wgd_status_timing.tsv"), delim = "\t")

wgd_plot_df <- clone_wgd_status_timing_df %>%
  filter(timing_category %in% timing_levels) %>%
  group_by(timing_category, wgd) %>%
  summarise(
    n_clones = n(),
    .groups = "drop"
  ) %>%
  group_by(timing_category) %>%
  mutate(
    proportion = n_clones / sum(n_clones),
    wgd_status = ifelse(
      wgd,
      "Genome doubled",
      "Non-genome doubled"
    )
  ) %>%
  ungroup() %>%
  mutate(
    timing_category = factor(
      timing_category,
      levels = timing_levels
    )
  )

# change plot, remove met clone, add black bar and coloutrs and counts

ggplot(wgd_plot_df, aes(x = timing_category, y = proportion, fill = wgd_status)) +
  geom_col(width = 0.7, color = "black") +
  scale_fill_manual(
    values = c(
      "Genome doubled" = "#765396",
      "Non-genome doubled" = "grey80"
    )
  ) +
  
  scale_y_continuous(
    labels = percent_format(),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  labs(
    x = NULL,
    y = "Proportion of clones",
    fill = NULL
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "top"
  )

ggsave(paste0(outputs.folder, date, "_prop_clones_wgd.png"), width = 7, height = 5)




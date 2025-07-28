#==============================================================================#
#==============================================================================#
######                                                                    ######
######           Seeding patterns from pre-invasive to carcinoma          ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Katherine Honan
# Date: 2025-06-30

setwd("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/")

####################################################
#### Source required functions & load libraries ####
####################################################

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

#############################################
#### Make a folder for this analysis run ####
#############################################

date <- gsub("-","",Sys.Date())

analysis_name <- 'preinveasive_seeding'
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

# read in data
trees <- readRDS("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/inputs/_RELEASE/_aggregate/conipher_trees/2025_06_28_cohort_conipher_trees.RDS")
metadata <- fread("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/EarlyPancancerAtlas/inputs/20250514_chen_2017_ESCC_metadata_EPA_mapping.txt")

# source both original (tissue) mets functions and ctdna adapted functions
source("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/personalis_ctDNA_mets_analysis/scripts/mets_functions.R")
source("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/personalis_ctDNA_mets_analysis/scripts/ctdna_mets_functions.R")
source("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/Somatic-Evolutionary-Monitoring-Lab/personalis_ctDNA_mets_analysis/scripts/useful_functions.R")

##########################
######## Defaults ########
##########################

use.original.clonality <- TRUE
ccf.buffer.toUse <- 10
save.output <- TRUE

##################################################
######### EPA ids & data pre-processing ##########
##################################################

cohort_overview <- metadata %>%
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

names(trees) <- sapply(trees, function(x) x$parameters$sampleID)

##############################################################
############## Annotate tumour level clonality ###############
##############################################################

clonalities <- lapply(trees, function(trees) {
  if (use.original.clonality) {
    clonality <- trees$clonality_out$clonality_table_original
  } else {
    clonality <- trees$clonality_out$clonality_table_corrected
  }
})


clonalities_updated <- lapply(names(clonalities), function(name) {
  clonality_pat <- clonalities[[name]]
  regions <- colnames(clonality_pat)
  # check if there are tumour regions in the tree
  cancer_regions <- as.character(cohort_overview %>% filter(sample_name_hash %in% regions, sample_type == "cancer") %>% pull(sample_name_hash))
  if (length(cancer_regions) == 0){
    print(paste0(name, " has no tumour regions in the tree"))
  } else {
    # get tumour level clonality
    clonality_pat$primary_clonality <- get.tumLevel.clonality(clonality_pat, cancer_regions)
  }
  # check if there are preinvasive regions in the tree
  preinvasive_regions <- as.character(cohort_overview %>% filter(sample_name_hash %in% regions, sample_type == "preinvasive") %>% pull(sample_name_hash))
  if (length(preinvasive_regions) == 0){
    print(paste0(name, " has no pre-invasive regions in the tree"))
  } else {
    # get pre-invasive clonality
    clonality_pat$preinvasive_clonality <- get.tumLevel.clonality(clonality_pat, preinvasive_regions)
  }
  return(clonality_pat)
})

names(clonalities_updated) <- names(clonalities)


############################################################
############## Determining seeding clonality ###############
############################################################

seedingClonality <- lapply(trees, function(fullTreeOutput) {
  
  print(fullTreeOutput$parameters$sampleID)
  
  if (use.original.clonality) {
    clonality <- fullTreeOutput$clonality_out$clonality_table_original
  } else {
    clonality <- fullTreeOutput$clonality_out$clonality_table_corrected
  }
  
  regions <- colnames(clonality)
  preinvasive.regions <- as.character(cohort_overview %>% filter(sample_name_hash %in% regions, sample_type == "preinvasive") %>% pull(sample_name_hash))
  cancer.regions <- as.character(cohort_overview %>% filter(sample_name_hash %in% regions, sample_type == "cancer") %>% pull(sample_name_hash))
  
  if (length(preinvasive.regions) == 0 | length(cancer.regions) == 0){
    print(paste0("Either no preinvaseive or cancer samples availible for clonality calculations"))
    return(NULL)
  } else {
    tree <- fullTreeOutput$graph_pyclone$default_tree
    
    # get seeding clonality for each met region
    seedingPattern.region.cancer <- sapply(cancer.regions, get.seedingPattern.region, tree = tree, clonality = clonality, primary.regions = preinvasive.regions)
  
    # assign overall, case level seeding clonality using all mets
    if (any(seedingPattern.region.cancer == "polyclonal")) {
      overallSeeding.cancer <- "polyclonal"
    } else {
      if (length(seedingPattern.region.cancer) == 1) {
        overallSeeding.cancer <- "monoclonal"
      } else {
        # may need to use get.overallClonality function
        overallSeeding.cancer <- get.overallClonality(tree, clonality, cancer.regions, preinvasive.regions)
      }
    }
    return(overallSeeding.cancer)
  }
})

seedingPhyletic <- lapply(trees, function(fullTreeOutput) {
  
  if (use.original.clonality) {
    clonality <- fullTreeOutput$clonality_out$clonality_table_original
  } else {
    clonality <- fullTreeOutput$clonality_out$clonality_table_corrected
  }
  
  pat <- fullTreeOutput$parameters$sampleID
  print(pat)
  
  regions <- colnames(clonality)
  preinvasive.regions <- as.character(cohort_overview %>% filter(sample_name_hash %in% regions, sample_type == "preinvasive") %>% pull(sample_name_hash))
  cancer.regions <- as.character(cohort_overview %>% filter(sample_name_hash %in% regions, sample_type == "cancer") %>% pull(sample_name_hash))
  
  if (length(preinvasive.regions) == 0 | length(cancer.regions) == 0){
    print(paste0("Either no preinvaseive or cancer samples availible for clonality calculations"))
    return(NULL)
  } else {

    tree            <- fullTreeOutput$graph_pyclone$default_tree
    trunk           <- fullTreeOutput$graph_pyclone$trunk
    
    seedingPattern.region.cancer <- sapply(cancer.regions, get.seedingPattern.region, tree = tree, clonality = clonality, primary.regions = preinvasive.regions)
  
    if (any(seedingPattern.region.cancer == "polyclonal")) {
      overallSeeding.cancer <- "polyclonal"
    } else {
      if (length(seedingPattern.region.cancer) == 1) {
        overallSeeding.cancer <- "monoclonal"
      } else {
        overallSeeding.cancer <- get.overallClonality(tree, clonality, cancer.regions, preinvasive.regions)
      }
    }
    
    if (overallSeeding.cancer == "monoclonal") {
      overallPhyletic.cancer <- "monophyletic"
    } else {
      overallPhyletic.cancer <- get.overallPhyletic(tree, trunk, clonality, cancer.regions, preinvasive.regions)
    }
    return(overallPhyletic.cancer)
  
  }
})

names(seedingPhyletic) <- names(trees)

# Seeding clonality dataframe
seeding_table <- data.frame(tumour_id = names(unlist(seedingClonality)), 
                            seeding_clonality = unlist(seedingClonality), 
                            seeding_phyletic = unlist(seedingPhyletic))

if (save.output) {
  saveRDS(seeding_table, file = paste0(outputs.folder, date, "_seeding_clonality.rds"))
}

##########################################################################################
############## Work out which preinvasive clones seeded the primary tumour ###############
##########################################################################################

cloneInfo.list <- lapply(names(trees), function(pat) {
  print(pat)
  
  treeStructure <- get.treeStructure(trees[[pat]]$graph_pyclone$default_tree, trees[[pat]]$graph_pyclone$trunk)
  
  if (use.original.clonality) {
    clonality <- trees[[pat]]$clonality_out$clonality_table_original
  } else {
    clonality <- trees[[pat]]$clonality_out$clonality_table_corrected
  }
  
  allPatientClones <- rownames(clonality)
  regions <- colnames(clonality)
  preinvasive.regions <- as.character(cohort_overview %>% filter(sample_name_hash %in% regions, sample_type == "preinvasive") %>% pull(sample_name_hash))
  cancer.regions <- as.character(cohort_overview %>% filter(sample_name_hash %in% regions, sample_type == "cancer") %>% pull(sample_name_hash))
  
  # Cases with no preinvasive regions
  if (length(preinvasive.regions) == 0) {
    print(paste0(pat, ": No preinvasive clones in the tree, skipping seeding clone estimation."))
    clonality$TumourClonality <- get.tumLevel.clonality(clonality, cancer.regions)
    tmp.return <- list(seedingClones = NA,
                       allTreeClones = treeStructure$allTreeClones,
                       clonalClones  = rownames(clonality[clonality$TumourClonality == "clonal", , drop = FALSE]), 
                       sharedClones  = NA,
                       preinvasiveClones = NA,
                       tumourClones     = rownames(clonality))
    df <- data.frame(Patient = pat, clones = allPatientClones, stringsAsFactors = FALSE)
    df$treeClones    <- df$clones %in% tmp.return$allTreeClones
    df$seedingClones <- df$clones %in% tmp.return$seedingClones
    df$clonalClones  <- df$clones %in% tmp.return$clonalClones
    df$sharedClones  <- df$clones %in% tmp.return$sharedClones
    df$preinvasiveClones <- df$clones %in% tmp.return$preinvasiveClones
    df$tumourClones     <- df$clones %in% tmp.return$tumourClones
    return(c(tmp.return, list(df.cloneInfo = df, individualSeedingClones = NA, seedingRegions = NA, seedingRegionInfo = NA)))
    }
  
  # Cases with no cancer regions
  if (length(cancer.regions) == 0){
    print(paste0(pat, ": No tumour clones in the tree, skipping seeding clone estimation."))
    clonality$PreinvasiveClonality <- get.tumLevel.clonality(clonality, preinvasive.regions)
    tmp.return <- list(seedingClones = NA,
                       allTreeClones = treeStructure$allTreeClones,
                       clonalClones  = rownames(clonality[clonality$PreinvasiveClonality == "clonal", , drop = FALSE]), 
                       sharedClones  = NA,
                       preinvasiveClones = rownames(clonality),
                       tumourClones     = NA)
    df <- data.frame(Patient = pat, clones = allPatientClones, stringsAsFactors = FALSE)
    df$treeClones    <- df$clones %in% tmp.return$allTreeClones
    df$seedingClones <- df$clones %in% tmp.return$seedingClones
    df$clonalClones  <- df$clones %in% tmp.return$clonalClones
    df$sharedClones  <- df$clones %in% tmp.return$sharedClones
    df$preinvasiveClones <- df$clones %in% tmp.return$preinvasiveClones
    df$tumourClones     <- df$clones %in% tmp.return$tumourClones
    return(c(tmp.return, list(df.cloneInfo = df, individualSeedingClones = NA, seedingRegions = NA, seedingRegionInfo = NA)))
    
  }
  
  # Cases with at least one preinvasive and tumour sample, work out seeding clonlity
  clonality$PreinvasiveClonality <- get.tumLevel.clonality(clonality, preinvasive.regions)
  clonality$TumourClonality <- get.tumLevel.clonality(clonality, cancer.regions)
  clonality.full <- clonality
  ### keep all clones - even non-tree to be able to classify more mutations 
  # clonality <- subset(clonality, rownames(clonality) %in% treeStructure$allTreeClones)
  if (nrow(clonality) == 0) {
    print(paste0(pat, ": No clones left after removing clones not in tree."))
    return(NULL)
  }
  
  allSharedClones <- intersect(rownames(subset(clonality, clonality$PreinvasiveClonality != "absent")),
                               rownames(subset(clonality, clonality$TumourClonality != "absent")))
  
  # Cases with no shared clones
  if (length(allSharedClones) == 0) {
    print(paste0(pat, ": No shared clones."))
    return(NULL)
  }
  
  # Cases with one shared clone
  if (length(allSharedClones) == 1) {
    if (clonality[which(rownames(clonality) == allSharedClones), "PreinvasiveClonality"] == "clonal" & clonality[which(rownames(clonality) == allSharedClones), "TumourClonality"] == "clonal") {
      print(paste0(pat, ": Only clonal cluster left."))
      individualSeedingClones <- setNames(rep(allSharedClones, length(cancer.regions)), cancer.regions)
      seedingRegions <- list(preinvasive.regions)
      names(seedingRegions) <- allSharedClones
      
      seedingRegionsDF <- data.frame(Patient = pat, PreinvasiveRegion = preinvasive.regions, Initiating = ifelse(cancer.regions %in% unique(unlist(seedingRegions)), TRUE, FALSE), Carcinoma = "", stringsAsFactors = FALSE)
      seedingRegionsDF$Carcinoma <- sapply(seedingRegionsDF$PreinvasiveRegion, function(reg) {
        tmp.indx <- grep(reg, seedingRegions)
        if (length(tmp.indx) == 0) return("")
        
        tmp.seeding <- names(seedingRegions[tmp.indx])
        
        tmp.out <- sapply(tmp.seeding, function(x) {
          names(individualSeedingClones[grep(x, individualSeedingClones)])
        })
        tmp.out <- paste0(unique(unlist(tmp.out)), collapse = ";")
        return(tmp.out)
      })
      
      tmp.return <- list(seedingClones = allSharedClones,
                         allTreeClones = treeStructure$allTreeClones,
                         clonalClones  = rownames(subset(clonality.full, clonality.full$PreinvasiveClonality == "clonal" & clonality.full$TumourClonality == "clonal")), 
                         sharedClones  = rownames(subset(clonality.full, clonality.full$PreinvasiveClonality %in% c("subclonal", "clonal") & clonality.full$TumourClonality %in% c("subclonal", "clonal"))), 
                         preinvasiveClones = rownames(subset(clonality.full, clonality.full$PreinvasiveClonality %in% c("subclonal", "clonal") & clonality.full$TumourClonality == "absent")), 
                         tumourClones     = rownames(subset(clonality.full, clonality.full$PreinvasiveClonality == "absent" & clonality.full$TumourClonality %in% c("subclonal", "clonal"))))
      df <- data.frame(Patient = pat, clones = allPatientClones, stringsAsFactors = FALSE)
      df$treeClones    <- df$clones %in% tmp.return$allTreeClones
      df$seedingClones <- df$clones %in% tmp.return$seedingClones
      df$clonalClones  <- df$clones %in% tmp.return$clonalClones
      df$sharedClones  <- df$clones %in% tmp.return$sharedClones
      df$preinvasiveClones <- df$clones %in% tmp.return$preinvasiveClones
      df$tumourClones     <- df$clones %in% tmp.return$tumourClones
      return(c(tmp.return, list(df.cloneInfo = df, individualSeedingClones = individualSeedingClones, seedingRegions = seedingRegions, seedingRegionInfo = seedingRegionsDF)))
    } else {
      print(paste0(pat, ": Only clone left is not clonal in both primary and met."))
      return(NULL)
    }
  }
  
  treeBranches <- strsplit(treeStructure$structure, split = ":")
  
  individualSeedingClones <- lapply(cancer.regions, function(met) {
    sharedClones  <- subset(allSharedClones, allSharedClones %in% rownames(subset(clonality, clonality[, met] != "absent")))
    seedingClones <- unique(get.seedingClones(clonality, sharedClones, treeBranches, prim = "PreinvasiveClonality", met = met))
    
    if (length(seedingClones) == 0) {
      print(paste0(pat, ": No seeding clones found for region"))
      return(NULL)
    }
    
    if (length(seedingClones) == 1) {
      if (clonality[which(rownames(clonality) %in% seedingClones), met] == "clonal") {
        return(seedingClones)
      } else {
        tmp.sharedClones <- subset(sharedClones, !sharedClones %in% seedingClones)
        while (length(tmp.sharedClones) != 0) {
          tmp.seedingClones <- get.seedingClones(clonality, tmp.sharedClones, treeBranches, met = met)
          seedingClones <- c(seedingClones, tmp.seedingClones)
          if (clonality[which(rownames(clonality) %in% tmp.seedingClones), met] == "clonal") {
            break
          } else {
            tmp.sharedClones <- subset(tmp.sharedClones, !tmp.sharedClones %in% tmp.seedingClones)
          }
        }
        return(seedingClones)
      }
    }
    
    if (length(seedingClones) > 1) {
      clonality.corrected <- trees[[pat]]$clonality_out$clonality_table_corrected
      if (all(clonality[which(rownames(clonality) %in% seedingClones), met] == "clonal") & all(clonality.corrected[which(rownames(clonality.corrected) %in% seedingClones), met] == "clonal")) {
        print(paste0(pat, ": All seeding clones defined as clonal even after correction of clonality. Check again."))
        return(NULL)
      } else {
        # get pairwise comparison of seeding clones in this region
        cloneComparisons <- combn(seedingClones, 2, simplify = FALSE)
        # check if these pairs are on same branch of tree
        tmp <- sapply(cloneComparisons, function(clones) {
          f.sameBranch(clones[1], clones[2], treeStructure$structure)
        })
        # start with all seeding clones at this region 
        seedingClones.toCount <- seedingClones
        # if any clones are on the same branch, remove the child of the pair from seedingClones.toCount
        if (any(tmp)) {
          for (x in cloneComparisons) {
            if (f.sameBranch(x[1], x[2], treeStructure$structure)) {
              seedingClones.toCount <- seedingClones.toCount[-which(seedingClones.toCount %in% x[which(x != f.findParent(x[1], x[2], treeStructure$structure))])]
            }
          }
        } 
        ccf.values <- trees[[pat]]$nested_pyclone$ccf_cluster_table[, met]
        ccf.buffer <- trees[[pat]]$ccf_buffer
        if (is.null(ccf.buffer)) {
          ccf.buffer <- ccf.buffer.toUse
        }
        # make sure the sum of the remaining seeing clones (parents across branches) sum to > 90 
        # or that any of the seeding clones are defined as clonal in the met region
        if (sum(ccf.values[seedingClones.toCount]) > 100 - ccf.buffer | any(clonality[which(rownames(clonality) %in% seedingClones), met] == "clonal")) {
          return(seedingClones)
        } else {
          # if not, identify additional potential seeding clones from shared clones
          tmp.sharedClones <- subset(sharedClones, !sharedClones %in% seedingClones)
          while (length(tmp.sharedClones) != 0) {
            # go through remaining shared clones and add to the set of seeding clones
            # until either:
            # total ccf of selected seeding clones reaches 90
            # some of the newly identified clones are clonal in the met
            # no more shared clones remain to be tested
            
            # get new candidate seeding clones from shared clones that are not yet classified as seeding clones
            tmp.seedingClones <- get.seedingClones(clonality, tmp.sharedClones, treeBranches, met = met)
            
            # for each newly identified clone (x)
            for (x in tmp.seedingClones) {
              # and each previously selected seeding clone (y)
              for (y in seedingClones.toCount) {
                # if same branch, keep the one with the higher ccf in the met region and remove the other
                if (f.sameBranch(x, y, treeStructure$structure)) {
                  if (ccf.values[x] - ccf.values[y] > 0) {
                    seedingClones <- c(seedingClones, x)
                    seedingClones.toCount <- c(seedingClones.toCount, x)
                    seedingClones.toCount <- seedingClones.toCount[-which(seedingClones.toCount == y)]    
                  } else {
                    next
                  }
                } else {
                  next
                }
              }
            }
            # remove duplicates
            seedingClones <- unique(seedingClones)
            seedingClones.toCount <- unique(seedingClones.toCount)
            # break out of the loop if the total ccf of the selected clones exceeds 90
            # or if one of the new clones is clonal at this region
            if (sum(ccf.values[seedingClones.toCount]) > 100 - ccf.buffer | any(clonality[which(rownames(clonality) %in% tmp.seedingClones), met] == "clonal")) {
              break
            } else {
              # remove processed clones from shared clones and continue testing others
              tmp.sharedClones <- subset(tmp.sharedClones, !tmp.sharedClones %in% tmp.seedingClones)
            }
          }
          return(seedingClones)
        }
      }
    }
  })
  
  if (any(sapply(individualSeedingClones, is.null))) {
    print(paste0(pat, ": Something is off. Please check"))
    return(NULL)
  }
  
  names(individualSeedingClones) <- cancer.regions
  seedingClones <- unique(unlist(individualSeedingClones))
  
  seedingRegions <- lapply(seedingClones, function(x) {
    preinvasive.regions[which(clonality[which(rownames(clonality) == x), preinvasive.regions] != "absent")]
  })
  
  names(seedingRegions) <- seedingClones
  
  seedingRegionsDF <- data.frame(Patient = pat, PreinvasiveRegion = preinvasive.regions, Initiating = ifelse(preinvasive.regions %in% unique(unlist(seedingRegions)), TRUE, FALSE), Carcinoma = "", stringsAsFactors = FALSE)
  
  seedingRegionsDF$Carcinoma <- sapply(seedingRegionsDF$PreinvasiveRegion, function(reg) {
    tmp.indx <- grep(reg, seedingRegions)
    if (length(tmp.indx) == 0) return("")
    
    tmp.seeding <- names(seedingRegions[tmp.indx])
    
    tmp.out <- sapply(tmp.seeding, function(x) {
      names(individualSeedingClones[grep(x, individualSeedingClones)])
    })
    
    tmp.out <- paste0(unique(unlist(tmp.out)), collapse = ";")
    return(tmp.out)
  })
  
  tmp.return <- list(seedingClones = seedingClones, 
                     allTreeClones = treeStructure$allTreeClones,
                     clonalClones  = rownames(subset(clonality.full, clonality.full$PreinvasiveClonality == "clonal" & clonality.full$TumourClonality == "clonal")), 
                     sharedClones  = rownames(subset(clonality.full, clonality.full$PreinvasiveClonality %in% c("subclonal", "clonal") & clonality.full$TumourClonality %in% c("subclonal", "clonal"))), 
                     preinvasiveClones = rownames(subset(clonality.full, clonality.full$PreinvasiveClonality %in% c("subclonal", "clonal") & clonality.full$TumourClonality == "absent")), 
                     tumourClones     = rownames(subset(clonality.full, clonality.full$PreinvasiveClonality == "absent" & clonality.full$TumourClonality %in% c("subclonal", "clonal"))))
  
  df <- data.frame(Patient = pat, clones = allPatientClones, stringsAsFactors = FALSE)
  df$treeClones    <- df$clones %in% tmp.return$allTreeClones
  df$seedingClones <- df$clones %in% tmp.return$seedingClones
  df$clonalClones  <- df$clones %in% tmp.return$clonalClones
  df$sharedClones  <- df$clones %in% tmp.return$sharedClones
  df$preinvasiveClones <- df$clones %in% tmp.return$preinvasiveClones
  df$tumourClones     <- df$clones %in% tmp.return$tumourClones
  
  return(c(tmp.return, list(df.cloneInfo = df, individualSeedingClones = individualSeedingClones, seedingRegions = seedingRegions, seedingRegionInfo = seedingRegionsDF)))
})

names(cloneInfo.list) <- names(trees)

if (save.output) saveRDS(cloneInfo.list, file = paste0(outputs.folder, date, "_cloneInfo.list.rds"))


#########################################
############## Plot trees ###############
#########################################

tree_plots <- list()
# Loop through all tumour IDs in 'trees'
for (tumour_id in names(trees)) {
  print(tumour_id)
  
  # Plot for those where pre-invasive samples weren't included on tree
  if (is.null(cloneInfo.list[[tumour_id]])) {
    tree <- trees[[tumour_id]]$graph_pyclone$default_tree
    tree_clones <- unique(c(tree[, 1], tree[, 2]))
    trunk <- trees[[tumour_id]]$graph_pyclone$trunk
    edgelength <- trees[[tumour_id]]$graph_pyclone$edgelength

    # Create the plot
    p <- plottingTxTree(
      tree, trunk, edgelength,
      seedingCluster = NA,
      color.mapping = "single", 
      col.palette = c("#009E73"),
      tumorID = tumour_id,
      add.title = TRUE,
      node.size = 7
    )
    
    # Append to list
    tree_plots[[tumour_id]] <- p
    
    } else {
      
      tree <- trees[[tumour_id]]$graph_pyclone$default_tree
      tree_clones <- unique(c(tree[, 1], tree[, 2]))
      trunk <- trees[[tumour_id]]$graph_pyclone$trunk
      edgelength <- trees[[tumour_id]]$graph_pyclone$edgelength
      seedingCluster <- cloneInfo.list[[tumour_id]]$seedingClones
    
      # Color assignment
      preinvasive_unique <- cloneInfo.list[[tumour_id]]$preinvasiveClones[cloneInfo.list[[tumour_id]]$preinvasiveClones %in% tree_clones]
      cancer_unique <- cloneInfo.list[[tumour_id]]$tumourClones[cloneInfo.list[[tumour_id]]$tumourClones %in% tree_clones]
      shared_nodes <- cloneInfo.list[[tumour_id]]$sharedClones[cloneInfo.list[[tumour_id]]$sharedClones %in% tree_clones]
    
      node_colors <- c(setNames(rep("#D55E00", length(preinvasive_unique)), preinvasive_unique),
                       setNames(rep("#009E73", length(cancer_unique)), cancer_unique),
                       setNames(rep("#D0DDE2", length(shared_nodes)), shared_nodes))
    
      # Create the plot
      p <- plottingTxTree(
        tree, trunk, edgelength,
        seedingCluster = seedingCluster,
        color.mapping = "perNode",
        col.palette = node_colors,
        tumorID = tumour_id,
        add.title = TRUE,
        node.size = 7
        )
    
      # Append to list
      tree_plots[[tumour_id]] <- p
  }
}


tree_plots <- lapply(tree_plots, function(p) {
  p + theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
            plot.title = element_text(size = 20, hjust = 0.5, face = "bold") )
})

wrap_plots(tree_plots, ncol = 6) 

ggsave(paste0(outputs.folder, date, "_forest.png"), width = 19, height = 20, dpi = 300)

# create legend with dummy data
legend_data <- data.frame(
  x = 1:5,
  y = 1,
  type = c("Preinvasive unique", "Cancer unique", "Shared", "Seeding clone", "Non-seeding clone"),
  fill = c("#D55E00", "#009E73", "#D0DDE2", "white", "white"),
  stroke = c(1, 1, 1, 2, 1),  # outline width
  shape = 21
)

# Build legend as a custom ggplot
legend_plot <- ggplot(legend_data, aes(x = x, y = y)) +
  geom_point(aes(fill = fill, shape = shape, stroke = stroke), size = 8, color = "black") +
  geom_text(aes(label = type), vjust = -2, size = 6) +
  scale_fill_identity() +
  scale_shape_identity() +
  scale_size_identity() +
  theme_void() +
  xlim(0.5, 5.5) +
  ylim(0.8, 1.5) +
  theme(legend.position = "none")

ggsave(paste0(outputs.folder, date, "_legend.png"), width = 9, height = 1, dpi = 300)



#!/usr/bin/env Rscript
# Author: Katherine Honan
# Date: 2026-03-26
# Description: Flag failed fastqs in inventory

#================#
# Load libraries # 
#================#

library(dplyr)
library(readr)

#=========================#
# Create output directory #
#=========================#

setwd("/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/EarlyPancancerAtlas/")

date <- gsub("-","",Sys.Date())

analysis_name <- 'inventory_update/flag_fastqs'
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

latest_inventory <- read_tsv('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/epa-project.inventory.latest.pass.tsv')
fastq_blacklist <- read_delim('/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/epa-project/inventory/fq_fail_list.txt', delim = '\t')

#========#
# script #
#========#

# before editing the inventory, create a backup of most recent version in case of problem
write_delim(latest_inventory, paste0(outputs.folder, "epa-project.inventory.latest.pass.bak.tsv"), delim = "\t")

# set fq_pass to FALSE for all files in fastq_blacklist if either R1 or R2 fail
failed_fastqs <- fastq_blacklist$failed_fastq
inventory_updated[inventory_updated$fq1 %in% failed_fastqs,]$fq_pass <- FALSE
inventory_updated[inventory_updated$fq2 %in% failed_fastqs,]$fq_pass <- FALSE

# write out edited inventory
write_delim(inventory_updated, paste0(outputs.folder, "epa-project.inventory.latest.pass.tsv"), delim = "\t")



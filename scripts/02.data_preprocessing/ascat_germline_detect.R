#!/usr/bin/env Rscript
# Author: Katie Honan
# Description: Perform minimal ASCAT run in order to determine germline samples with flat copy number profiles

#================#
# Load libraries #
#================#

suppressPackageStartupMessages({
  library(ASCAT)
  library(optparse)
})

#=================#
# Parse arguments #
#=================#

option_list <- list(
  make_option("--bam",
              type = "character",
              help = "Tumour BAM/CRAM file",
              metavar = "character"),

  make_option("--sample_id",
              type = "character",
              help = "Sample name hash",
              metavar = "character"),

  make_option("--outdir",
              type = "character",
              help = "Output directory",
              metavar = "character")
)

parser <- OptionParser(option_list = option_list,
                       description = "Run minimal ASCAT")

args <- parse_args(parser)

required_args <- c("bam", "sample_id", "outdir")
missing_args <- required_args[
  !required_args %in% names(args) |
  sapply(required_args, function(x) is.null(args[[x]]))
]

if (length(missing_args) > 0) {
  stop(paste("Missing required arguments:", paste(missing_args, collapse = ", ")))
}

bam       <- normalizePath(args$bam)
sample_id <- args$sample_id
outdir    <- args$outdir

# Create output directory if needed
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

setwd(outdir)

#==========#
# Defaults #
#==========#

exome_bed <- "/rds/project/rds-LH0AvU65IRI/_PIPELINE/TRACERx-assets/v3/capture_targets/SureSelectV5/S04380110_hg38.padded50.reduced.bed"
loci.prefix = "/rds/project/rds-LH0AvU65IRI/_PIPELINE/TRACERx-assets/v3/ascat_files/hg38/G1000_lociAll_hg38/G1000_loci_hg38_chr"
alleles.prefix = "/rds/project/rds-LH0AvU65IRI/_PIPELINE/TRACERx-assets/v3/ascat_files/hg38/G1000_allelesAll_hg38/G1000_alleles_hg38_chr"
GCcontentfile <- "/rds/project/rds-LH0AvU65IRI/_PIPELINE/TRACERx-assets/v3/ascat_files/hg38/GC_G1000_hg38.txt"
replictimingfile <- "/rds/project/rds-LH0AvU65IRI/_PIPELINE/TRACERx-assets/v3/ascat_files/hg38/RT_G1000_hg38.txt"

tumour_logR <- file.path(outdir, paste0(sample_id, "_Tumor_LogR.txt"))
tumour_BAF  <- file.path(outdir, paste0(sample_id, "_Tumor_BAF.txt"))

#=============#
# STEP 1      #
#=============#
print("Starting ascat.prepareHTS...")
ascat.prepareHTS(
  tumourseqfile                   = bam,
  tumourname                      = sample_id,
  allelecounter_exe               = "/home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/_ENV_v2/conda/ascat_hts/bin/alleleCounter",
  alleles.prefix                  = alleles.prefix,
  loci.prefix                     = loci.prefix,
  genomeVersion                   = "hg38",
  gender                          = "XX",
  chrom_names                     = c(1:22, "X"),
  nthreads                        = 4,
  minCounts                       = 3,
  tumourLogR_file                 = tumour_logR,
  tumourBAF_file                  = tumour_BAF,
  BED_file                        = exome_bed,
  additional_allelecounter_flags  = "-r /rds/project/rds-LH0AvU65IRI/_PIPELINE/TRACERx-assets/v3/reference/hg38/hg38_u2af1/GRCh38.d1.vd1.chr21fix.fa"
)

print("ascat.prepareHTS complete!")

#=============#
# STEP 2      #
#=============#

ascat.bc <- ascat.loadData(
  Tumor_LogR_file = tumour_logR,
  Tumor_BAF_file  = tumour_BAF,
  genomeVersion   = "hg38"
)

#=============#
# STEP 3      #
#=============#

print("Starting ascat.plotRawData (Before)...")
ascat.plotRawData(ascat.bc, img.prefix = paste0(sample_id, "_Before_"))
print("ascat.plotRawData (Before) complete!")

print("Starting ascat.plotRawData (log corrected)...")
ascat.bc <- ascat.correctLogR(ascat.bc,
                              GCcontentfile = GCcontentfile,
                              replictimingfile = replictimingfile)

ascat.plotRawData(ascat.bc, img.prefix = paste0(sample_id, "_After_"))
print("ascat.plotRawData (log corrected) complete!")

#=============#
# STEP 4      #
#=============#

ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc)

ascat.bc <- ascat.aspcf(ascat.bc, penalty = 50, ascat.gg = ascat.gg)

ascat.plotSegmentedData(
  ascat.bc,
  logr.y_values = c(-2, 2),
  img.prefix = paste0(sample_id, "_ASPCF_")
)

#=============#
# STEP 5      #
#=============#

tryCatch({
  ascat.output_initial <- ascat.runAscat(
    ascat.bc,
    write_segments = TRUE,
    pdfPlot = TRUE,
    gamma = 1
  )}, error = function(e) {
  message("ASCAT failed for ", sample_id, ": ", e$message)
  return(NULL)
})

cat("Initial ASCAT purity: ", ascat.output_initial$aberrantcellfraction, "\n")
cat("Initial ASCAT ploidy: ", ascat.output_initial$ploidy, "\n")

QC <- ascat.metrics(ascat.bc, ascat.output_initial)

save(ascat.bc, ascat.output_initial, QC,
  file = paste0(outdir, "/", "ascat_step5_objects_", sample_id, ".RData"))


solution_grid <- ascat.output_initial$purityploidy

if (is.null(solution_grid) || nrow(as.data.frame(solution_grid)) == 0) {

  solution_grid <- data.frame(
    purity = ascat.output_initial$aberrantcellfraction,
    ploidy = ascat.output_initial$ploidy,
    goodnessOfFit = ascat.output_initial$goodnessOfFit
  )

} else {

  solution_grid <- as.data.frame(solution_grid)

}

solution_grid$sample <- sample_id
solution_grid <- solution_grid[order(-solution_grid$goodnessOfFit), ]

write.table(solution_grid,
            file = file.path(outdir, paste0(sample_id, "_purity_ploidy_grid.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

#=============#
# STEP 6      #
#=============#

top_n <- 10
available_n <- nrow(solution_grid)
top_solutions <- head(solution_grid, n = min(top_n, available_n))

solutions_base_dir <- file.path(outdir, "ascat_solutions")
dir.create(solutions_base_dir, showWarnings = FALSE)

for (i in seq_len(nrow(top_solutions))) {

  rho <- top_solutions$purity[i]
  psi <- top_solutions$ploidy[i]
  gof <- top_solutions$goodnessOfFit[i]

  sol_dir <- file.path(
    solutions_base_dir,
    paste0(sample_id, "_rho", sprintf("%.2f", rho), "_psi", sprintf("%.2f", psi))
  )

  dir.create(sol_dir, recursive = TRUE, showWarnings = FALSE)

  ascat.bc$aberrantcellfraction <- NULL
  ascat.bc$ploidy <- NULL

  alt_output <- ascat.runAscat(
    ascat.bc,
    write_segments = TRUE,
    pdfPlot = TRUE,
    gamma = 1,
    rho_manual = rho,
    psi_manual = psi
  )

  saveRDS(alt_output, file = file.path(sol_dir, "ASCAT_output.rds"))

  write.table(alt_output$segments,
              file = file.path(sol_dir, "segments.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  write.table(alt_output$segments_raw,
              file = file.path(sol_dir, "segments_raw.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

}


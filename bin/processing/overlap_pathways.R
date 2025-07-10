### Loading libraries
required_libraries <- c(
    "data.table",    
    "dplyr",
    "optparse",
    "rlang",
    "purrr",
    "stringr")

for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
}

### Options
options(stringsAsFactors = FALSE)

### Command line options
option_list <- list(

    optparse::make_option(
        c("-i", "--pathways"),
        type = "character",
      #  action = "append",
        default = NULL,
        help = "Path to the differential genesets file.",
        metavar = "character"),
    # optparse::make_option(
    #     c("-g", "--geneset"),
    #     type = "character",
    #     default = NULL,
    #     help = "Path to the gene set file (gmt).",
        # metavar = "character"),
    # optparse::make_option(
    #     c("-t", "--type"),
    #     type = "character",
    #     default = NULL,
    #     help = "The wildcard thing",
    #     metavar = "character")
    optparse::make_option(
        c("-o", "--output_dir"),
        type = "character",
        default = NULL,
        help = "Path to the output directory.",
        metavar = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser, positional_arguments = TRUE)

## Initialize variable
PATHWAY_FILES <- opt$args
#GENESET_FILE <- opt$geneset
#DATA_TYPE <- opt$type
OUTPUT_DIR <- opt$options$output_dir

## Debug
#PATHWAY_FILES <- list("/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results/differential_indegrees/cosgrove/differential_genesets.RData", 
#                        "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results/differential_indegrees/aurora/differential_genesets.RData")
# OUTPUT_DIR <- "test"

#PATHWAY_FILES <- as.list(PATHWAY_FILES)

## Functions
source("bin/processing/overlap_pathways_fn.R")

## Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

PATHWAY_FILES <- strsplit(PATHWAY_FILES, split = ",")
dss <- list("COS", "AUR")
column_to_merge <- "pathway"
all_pathways <- merge_all_pathways(PATHWAY_FILES, dss, column_to_merge)



fwrite(
    all_pathways,
    file = file.path(OUTPUT_DIR, "overlapping_pathways_all.tsv"),
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE
)
### Loading libraries
required_libraries <- c(
    "data.table",     ## Many (gg)plots in the same frame
    "tidyverse",
    "optparse",
    "limma",
    "rlang")

for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
}

### Options
options(stringsAsFactors = FALSE)

### Command line options
option_list <- list(

    optparse::make_option(
        c("-i", "--indegrees"),
        type = "character",
        default = NULL,
        help = "Path to the indegree file.",
        metavar = "character"),
    optparse::make_option(
        c("-m", "--metadata"),
        type = "character",
        default = NULL,
        help = "Path to the metadata file.",
        metavar = "character"),
    optparse::make_option(
        c("-o", "--output_dir"),
        type = "character",
        default = NULL,
        help = "Path to the output directory.",
        metavar = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## Initialize variable
METADATA_FILE <- opt$metadata
INDEGREES_FILE <- opt$indegrees
OUTPUT_DIR <- opt$output_dir

## Debug
# INDEGREES_FILE <- "data/geo_brca/filtered_cos_indegree.csv"
# METADATA_FILE <- "data/geo_brca/metadata.csv"
# OUTPUT_DIR <- "test"


## Functions
source("bin/preprocessing/compute_diff_indegrees_fn.R")

## Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Computing
indegree_df <- read_indegree(INDEGREES_FILE)
metadata_df <- read_and_match_metadata(
    METADATA_FILE,
    indegree_df)

diff_indegree_df <- create_toptable_paired(indegree_df, metadata_df, "patient", "sample_type")

## Exporting
save(
    diff_indegree_df,
    file.path(OUTPUT_DIR, "differential_indegrees.RData"))
data.table::fwrite(
    diff_indegree_df,
    file = file.path(OUTPUT_DIR, "differential_indegrees.tsv"),
    sep = "\t",
    col.names = TRUE)
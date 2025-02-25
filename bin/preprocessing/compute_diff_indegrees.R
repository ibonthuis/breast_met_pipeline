### Loading libraries
required_libraries <- c(
    "data.table",    
    "dplyr",
    "optparse",
    "limma",
    "rlang",
    "purrr")

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
opt <- optparse::parse_args(opt_parser)

## Initialize variable
METADATA_FILE <- opt$metadata
INDEGREES_FILE <- opt$indegrees
#DATA_TYPE <- opt$type
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
metadata_df <- data.table::fread(METADATA_FILE)
metadata_df <- match_metadata(
    metadata_df,
    indegree_df)

nr_of_prim <- sum(metadata_df$sample_type == "Primary")
nr_of_met <- sum(metadata_df$sample_type == "Metastasis")

if (nr_of_prim == nr_of_met) {
    paired_toptables <- create_toptable_paired(indegree_df, metadata_df, "patient", "sample_type")
    paired_toptables$gene <- rownames(paired_toptables)
    ranked_paired_toptables <- paired_toptables[, c("gene", "t")]
} else {
    list_of_pairs <- list()
    names_for_list <- list("first", "second", "third", "fourth", "fifth")
   for(seednr in 1:5){
        select_pair <- filter_pairs(indegree_df, metadata_df, seednr)
        select_pair <- list(select_pair)
        names(select_pair) <- names_for_list[seednr]
        list_of_pairs <- append(list_of_pairs, select_pair) # go from here
   }
    paired_toptables <- purrr::map(list_of_pairs, ~ create_toptable_paired(.x, metadata_df, "patient", "sample_type"))
    ranked_paired_toptables <- purrr::map(paired_toptables, ~ {
        df <- as.data.frame(.x)
        df <- df[, c("t"), drop = FALSE]
        df$gene <- rownames(df)
        df <- df[, c("gene", "t")]
        df
        })
}

## Exporting
save(
    paired_toptables,
    file = file.path(OUTPUT_DIR, "differential_indegrees.RData"))
# data.table::fwrite(
#     diff_indegree_df,
#     file = file.path(OUTPUT_DIR, "differential_indegrees.tsv"),
#     sep = "\t",
#     col.names = TRUE,
#     row.names = TRUE)

save(
    ranked_paired_toptables,
    file = file.path(OUTPUT_DIR, "differential_indegrees_rank_file.RData")
)
# data.table::fwrite(
#     ranked_indegree_df,
#     file = file.path(OUTPUT_DIR, "differential_indegrees_ranked_files.RData"),
#     sep = "\t",
#     col.names = FALSE,
#     row.names = FALSE)
### Loading libraries
required_libraries <- c(
    "data.table",    
    "dplyr",
    "optparse",
    "rlang",
    "ggplot2",
    "purrr")

for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
}

### Options
options(stringsAsFactors = FALSE)

### Command line options
option_list <- list(

    optparse::make_option(
        c("-g", "--gene_set_file"),
        type = "character",
        default = NULL,
        help = "Path to the file containing gene sets (.gmt)",
        metavar = "character"),
    optparse::make_option(
        c("-i", "--ranks"),
        type = "character",
        default = NULL,
        help = "Path to the ranked file. But snakemake will know this because it was 
        built upstream",
        metavar = "character") #,
    # optparse::make_option(
    #     c("-v", "--visualization_var"),
    #     type = "character",
    #     default = NULL,
    #     help = "Variable to visualize, one of the column names of the metadata",
    #     metavar = "character"),
    # optparse::make_option(
    #     c("-o", "--output_dir"),
    #     type = "character",
    #     default = NULL,
    #     help = "Path to the output directory.",
    #     metavar = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## Initialize variable
GENE_SET <- opt$gene_set_file
RANK_FILE <- opt$ranks

# INDEGREES_FILE <- opt$indegrees
# OUTPUT_DIR <- opt$output_dir
# VARIABLE <- opt$visualization_var

## Debug
# INDEGREES_FILE <- "data/geo_brca/filtered_cos_indegree.csv"
# METADATA_FILE <- "data/geo_brca/metadata.csv"
# OUTPUT_DIR <- "test"


## Functions
source("bin/preprocessing/compute_gse_fn.R")

## Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Gene set enrichment analysis
load(ranks)
rank_list <- ranked_paired_toptables

if (length(rank_list) > 1) {
   list_of_gsea <- purrr::map(rank_list, ~ perform_gsea(.x, GENE_SET))
} else {
    list_of_gsea <- perform_gsea(rank_list, GENE_SET)
}


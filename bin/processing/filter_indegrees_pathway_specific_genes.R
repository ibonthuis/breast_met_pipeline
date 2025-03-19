
### Loading libraries
required_libraries <- c(
    "data.table",    
    "dplyr",
    "optparse",
    "rlang",
    "ggplot2",
    "purrr",
    "tibble",
    "stringr",
    "limma",
    "tidyr")

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
        c("-g", "--gene_set_file"),
        type = "character",
        default = NULL,
        help = "Path to the geneset file.",
        metavar = "character"),
    optparse::make_option(
        c("-e", "--edge_file"),
        type = "character",
        default = NULL,
        help = "Path to the edge file.",
        metavar = "character"),
    optparse::make_option(
        c("-m", "--metadata"),
        type = "character",
        default = NULL,
        help = "Path to the metadata file.",
        metavar = "character"),
    optparse::make_option(
        c("-p", "--pathway_of_interest"),
        type = "character",
        default = NULL,
        help = "Character string with pathway of interest that was probably significant",
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
INDEGREES_FILE <- opt$indegrees
METADATA_FILE <- opt$metadata
EDGE_FILE <- opt$edge_file
PATHWAY <- opt$pathway_of_interest
GENE_SET_FILE <- opt$gene_set_file
OUTPUT_DIR <- opt$output_dir

## Debug
# INDEGREES_FILE <- "data/input_data/aurora/filtered_indegree.csv"
# PERMUTATION_INDEGREES_FILE <- "subselection_diff_indegrees/selected_metprim_pairs_indegrees.RData"
# METADATA_FILE <- "data/input_data/aurora/metadata.csv"
# PATHWAY <- "REACTOME_INTERFERON_SIGNALING"
# GENE_SET_FILE <- "data/gene_sets/c2.cp.reactome.v5.0.symbols.gmt"
# #EDGE_FILE <- "/storage/kuijjerarea/ine/breast_met/breast_other/tissue_batch_corr_allAUR/output/lioness_filtered_for_genes.csv"
# EDGE_FILE <- "/storage/kuijjerarea/ine/breast_met/breast_other/tissue_batch_corr_allAUR/edge_table_slightly_adjusted.tsv"
# OUTPUT_DIR <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/subselection_diff_indegrees"

## Functions
source("bin/preprocessing/compute_diff_indegrees_fn.R")
source("bin/processing/filter_analyze_pathway_specific_genes_fn.R")

# To do:
# filter expression data (PREP a file first that has the corrct column names matching to metadata)
metadata <- fread(METADATA_FILE)
indegrees_paired <- read_indegree_if(INDEGREES_FILE)

indegrees_filtered_by_pathway <- purrr::map(indegrees_paired, ~ filter_indegree_df(GENE_SET_FILE, PATHWAY, .x))

filtered_diff_indegrees <- purrr::map(indegrees_filtered_by_pathway, ~ create_toptable_paired(.x, metadata, "patient", "sample_type"))


save(
    filtered_diff_indegrees,
    file = file.path(OUTPUT_DIR, "diff_indegrees_pathway_specific_genes.RData")
)

binary <- grep(".RData", INDEGREES_FILE)
class(binary)

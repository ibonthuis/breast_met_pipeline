
### Loading libraries
required_libraries <- c(
    "data.table",    
    "dplyr",
    "optparse",
    "rlang",
    "ggplot2")

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
        help = "Pathway of interest that was probably significant",
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

## Functions
source("bin/preprocessing/compute_diff_indegrees_fn.R")
source("bin/processing/filter_analyze_pathway_specific_genes.R")


# reading_gene_set_file
# select_genes_from_pathway_of_interest
# filter indegrees by filter_df_by_pathway_of_interest
# 
# filter edges by filtering the target column filter_df_by_pathway_of_interest
# 

#
# I KIND OF DID THIS ALREADY:
# source compute_diff_indegrees:
# create_toptable_paired for subset of indegrees
# create_toptable_paired for subset of edges data
# 

# To do:
# filter expression data (PREP a file first that has the corrct column names matching to metadata)

indegrees_filtered_by_pathway <- filter_indegree_file(GENE_SET_FILE, PATHWAY, INDEGREES_FILE)

filtered_diff_indegrees <- create_toptable_paired_from_file(INDEGREES_FILE, METADATA_FILE, "patient", "sample_type")

filtered_diff_edges <- filter_edges(EDGE_FILE, METADATA_FILE, PATHWAY, 1.5)
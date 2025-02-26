### Loading libraries
required_libraries <- c(
    "data.table",    
    "dplyr",
    "optparse",
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
        c("-i", "--pathways"),
        type = "character",
        default = NULL,
        help = "Path to the differential genesets file.",
        metavar = "character"),
    optparse::make_option(
        c("-g", "--geneset"),
        type = "character",
        default = NULL,
        help = "Path to the gene set file (gmt).",
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
PATHWAY_FILE <- opt$pathways
GENESET_FILE <- opt$geneset
#DATA_TYPE <- opt$type
OUTPUT_DIR <- opt$output_dir

## Debug
PATHWAY_FILE <- list("snakemake_results/differential_indegrees/cosgrove/differential_genesets.RData", "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results/differential_indegrees/aurora/differential_genesets.RData")
# METADATA_FILE <- "data/geo_brca/metadata.csv"
# OUTPUT_DIR <- "test"


## Functions
source("bin/processing/overlap_pathways_fn.R")

## Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

enriched_pathways1 <- PATHWAY_FILE[[1]]
enriched_pathways2 <- PATHWAY_FILE[[2]]
ds1 <- "COS"
ds2 <- "AUR"

all_pathways <- merge_all_pathways(enriched_pathways1, enriched_pathways2, ds1, ds2, "pathway")
head(enriched_pathways2[[1]])
class(enriched_pathways2[[1]])

test <- merge_gsea_results(enriched_pathways2, "AUR", "pathway")

enr <- get(load(enriched_pathways2))
head(enr)

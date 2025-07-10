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
    # optparse::make_option(
    #     c("-e", "--edge_file"),
    #     type = "character",
    #     default = NULL,
    #     help = "Path to the edge file.",
        # metavar = "character"),
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
INDEGREES_FILE <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/data/input_data/cosgrove/filtered_indegree.csv"
METADATA_FILE <- "data/input_data/cosgrove/metadata.csv"
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
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


metadata <- fread(METADATA_FILE)
indegrees_paired <- read_indegree_if(INDEGREES_FILE)



diff_indegrees_cov <- purrr::map(indegrees_paired, ~ run_limma(.x, metadata = metadata, covariates = "tumor_purity", type_col = "sample_type"))

indegrees_filtered_by_pathway <- purrr::map(indegrees_paired, ~ filter_indegree_df(GENE_SET_FILE, PATHWAY, .x))

# indegrees_paired <- indegrees_paired[[1]]
# diff_indegrees_cov <- run_limma(indegrees_paired, metadata = metadata, covariates = "tumor_purity", type_col = "sample_type")
diff_indegrees_cov <- diff_indegrees_cov[[1]]
diff_indegrees_cov$genes <- rownames(diff_indegrees_cov)
diff_indegrees_cov_rnk <- diff_indegrees_cov[, c("genes", "t")]
# head(diff_indegrees_cov)

filtered_diff_indegrees <- purrr::map(indegrees_filtered_by_pathway, ~ create_toptable_paired(.x, metadata, "patient", "sample_type"))

# for (i in 1:length(filtered_diff_indegrees)) {
#    fwrite(filtered_diff_indegrees[[i]], paste0("/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/interferon_signaling_diff_aurora_indegrees/diff_indegrees_aur", i, ".tsv"), row.names = TRUE, sep = "\t")
# }


save(
    filtered_diff_indegrees,
    file = file.path(OUTPUT_DIR, "diff_indegrees_pathway_specific_genes.RData")
)

# filtered_diff_indegrees <- filtered_diff_indegrees[[1]]
# data.table::fwrite(
#     filtered_diff_indegrees,
#     file = "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/interferon_signaling_diff_cosgrove_indegrees/diff_indegrees_pathway_specific_genes.tsv",
#     sep = "\t",
#     col.names = TRUE,
#     row.names = TRUE)


# binary <- grep(".RData", INDEGREES_FILE)
# class(binary)


# dir.create("confounding_factor_diff_indegrees")
# fwrite(diff_indegrees_cov, "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/confounding_factor_diff_indegrees/diff_limma_analysis_with_confounding_factor.tsv", row.names = TRUE, sep = "\t")
# fwrite(diff_indegrees_cov_rnk, "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/confounding_factor_diff_indegrees/diff_limma_analysis_with_confounding_factor.rnk", col.names = FALSE, sep = "\t")

### Loading libraries
required_libraries <- c(
    "data.table",    
    "dplyr",
    "optparse",
    "limma",
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
        c("-c", "--covariate"),
        type = "character",
        default = NULL,
        help = "The covariate to correct for",
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
#DATA_TYPE <- opt$type
OUTPUT_DIR <- opt$output_dir
COVARIATE <- opt$covariate

## Debug
# INDEGREES_FILE <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/data/input_data/cosgrove/filtered_indegree.csv"
# METADATA_FILE <- "/storage/kuijjerarea/ine/breast_met/breast_other/R/data/metadata_cosgrove_purity_discrete.csv"

# INDEGREES_FILE <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/data/input_data/aurora/filtered_indegree.csv"
# METADATA_FILE <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/data/input_data/aurora/metadata.csv"

# OUTPUT_DIR <- "cov_analysis"
# COVARIATE <- "tumor_purity"

## Functions
source("bin/preprocessing/compute_diff_indegrees_fn.R") #/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/

## Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Computing
indegree_df <- read_indegree(INDEGREES_FILE)
metadata_df <- data.table::fread(METADATA_FILE)
metadata_df <- match_metadata(
    metadata_df,
    indegree_df)



metadata_df <- as.data.frame(metadata_df)
metadata_df <- metadata_df[!is.na(metadata_df[, str_which(colnames(metadata_df), COVARIATE)]), ]
indegree_df <- indegree_df[, colnames(indegree_df) %in% metadata_df$sample_name]


corrected_toptable <- run_limma(indegree_df, metadata_df, COVARIATE, "sample_type")
corrected_toptable$gene <- rownames(corrected_toptable)
ranked_corrected_toptable <- corrected_toptable[, c("gene", "t")]
#head(ranked_corrected_toptable, n = 10)



# which number of the column equals that of the covariate?




# check if covariate (tumor_purity) value is actually present. subselect those samples that have a tumor_purity value


save(
    ranked_corrected_toptable,
    file = file.path(OUTPUT_DIR, "ranked_corrected_tumor_purity_differential_indegrees.RData")
)



data.table::fwrite(
    corrected_toptable,
    file = file.path(OUTPUT_DIR, "corrected_tumor_purity_differential_indegrees.tsv"),
    sep = "\t",
    col.names = TRUE,
    row.names = TRUE)


#  if(covariate_correction == TRUE){
#         corrected_toptable <- run_limma(indegree_df, metadata_df, "tumor_purity", "sample_type")
#         corrected_toptable$gene <- rownames(corrected_toptable)
#         ranked_corrected_paired_toptable <- corrected_toptable[, c("gene", "t")] 
#         }
            # corrected_toptable <- purrr::map(list_of_pairs, ~ run_limma(.x, metadata_df, "tumor_purity", "sample_type"))
    # ranked_corrected_paired_toptable <- purrr::map(corrected_toptable, ~ {
    #     df <- as.data.frame(.x)
    #     df <- df[, c("t"), drop = FALSE]
    #     df$gene <- rownames(df)
    #     df <- df[, c("gene", "t")]
    #     df
    #     })
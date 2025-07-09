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

# INDEGREES_FILE <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/data/input_data/cosgrove/filtered_indegree.csv"
# METADATA_FILE <- "/storage/kuijjerarea/ine/breast_met/breast_other/R/data/metadata_cosgrove_purity_discrete.csv"
# OUTPUT_DIR <- "test2"

# INDEGREES_FILE <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/data/input_data/aurora/filtered_indegree.csv"
# METADATA_FILE <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/data/input_data/aurora/metadata.csv"

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

# arrange(colnames(indegree_df))

# head(metadata_df)

# metadata_ordered <- as.data.frame(metadata_df)
  
# metadata_ordered <- metadata_ordered %>%
#     filter(sample_name %in% colnames(indegree_df))
#   #nrow(metadata_ordered)
# patient_col <- "patient"
# type_col <- "sample_type"
# tumor_purity_col <- "tumor_purity_discrete"
#   patient_col <- rlang::ensym(patient_col)
#   type_col <- rlang::ensym(type_col)
#   tumor_purity_col <- rlang::ensym(tumor_purity_col)
#   pat <- factor(metadata_ordered[[as.character(patient_col)]])
#   tumor_purity <- factor(metadata_ordered[[as.character(tumor_purity_col)]])
#   type <- factor(metadata_ordered[[as.character(type_col)]])
#    tp_effect <- factor(paste(metadata_ordered$sample_type,metadata_ordered$tumor_purity_discrete,sep="."))
#   design <- model.matrix(~0 + tp_effect)
#   colnames(design)
#   corfit <- duplicateCorrelation(indegree_df, design, block = metadata_ordered$patient) #where patient is a categorical variable denoting which patient is which#

#   fit_adj <- lmFit(indegree_df, design, block=pat, correlation = corfit$consensius) 


if (nr_of_prim == nr_of_met) {
    paired_toptables <- create_toptable_paired(indegree_df, metadata_df, "patient", "sample_type")
    paired_toptables$gene <- rownames(paired_toptables)
    ranked_paired_toptables <- paired_toptables[, c("gene", "t")]

} else {
    #metadata_df <- metadata_df[!is.na(metadata_df$tumor_purity), ]
    #indegree_df <- indegree_df[, colnames(indegree_df) %in% metadata_df$sample_name]
    #nrow(metadata_df)
   
    list_of_pairs <- list()
    names_for_list <- list("first", "second", "third", "fourth", "fifth")
    for(seednr in 1:5){ # 5 permutations
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


save(
    ranked_paired_toptables,
    file = file.path(OUTPUT_DIR, "differential_indegrees_rank_file.RData")
)

# save(
#     corrected_toptable,
#     file = file.path(OUTPUT_DIR, "covariate_corr_differential_indegrees.RData"))

# save(
#     list_of_pairs,
#     file = file.path(OUTPUT_DIR, "selected_metprim_pairs.RData"))

# data.table::fwrite(
#     list_of_pairs,
#     file = file.path(OUTPUT_DIR, "selected_indegree_pairs.tsv"),
#     sep = "\t",
#     col.names = TRUE,
#     row.names = TRUE)



# save(
#     ranked_corrected_paired_toptable,
#     file = file.path(OUTPUT_DIR, "covariate_corr_differential_indegrees_rank_file.RData")

# )
# data.table::fwrite(
#     ranked_indegree_df,
#     file = file.path(OUTPUT_DIR, "differential_indegrees_ranked_files.RData"),
#     sep = "\t",
#     col.names = FALSE,
#     row.names = FALSE)


# save(
#     corrected_toptable,
#     file = file.path(OUTPUT_DIR, "corrected_tumor_purity_differential_indegrees.tsv")
# )

#OUTPUT_DIR <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/corr_1_aur"


# data.table::fwrite(
#     toptable,
#     file = file.path(OUTPUT_DIR, "corrected_tumor_purity_differential_indegrees_aurora.tsv"),
#     sep = "\t",
#     col.names = TRUE,
#     row.names = TRUE)

# data.table::fwrite(
#     indegree_df,
#     file = file.path(OUTPUT_DIR, "indegrees_with_tumor_purity_level_aurora.tsv"),
#     sep = "\t",
#     col.names = TRUE,
#     row.names = TRUE)

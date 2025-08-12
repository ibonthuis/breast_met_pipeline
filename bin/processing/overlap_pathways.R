### Loading libraries
required_libraries <- c(
    "data.table",    
    "dplyr",
    "optparse",
    "rlang",
    "purrr",
    "stringr",
    "ggvenn")

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
     #   action = "append",
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
opt <- optparse::parse_args(opt_parser)

## Initialize variable
PATHWAY_FILES <- opt$pathways
#GENESET_FILE <- opt$geneset
#DATA_TYPE <- opt$type
OUTPUT_DIR <- opt$options$output_dir

## Debug
# PATHWAY_FILES <- list("/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results/differential_indegrees/cosgrove/differential_genesets.RData", 
#                        "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results/differential_indegrees/aurora/differential_genesets.RData")
# OUTPUT_DIR <- "test"

# PATHWAY_FILES <- as.list(PATHWAY_FILES)


#   bin/processing/overlap_pathways.R 
# -i {snakemake_results_code_review/differential_indegrees/aurora/differential_genesets.RData, 
# snakemake_results_code_review/differential_indegrees/cosgrove/differential_genesets.RData }
# it should look like this eventually:
# -i {},{}

# -o snakemake_results_code_review/overlapping_pathways

## Functions
source("bin/processing/overlap_pathways_fn.R")

## Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
print(PATHWAY_FILES)
PATHWAY_FILES <- unlist(strsplit(PATHWAY_FILES, split = ","))
print(PATHWAY_FILES)

dss <- list("COS", "AUR")
column_to_merge <- "pathway"
all_pathways <- merge_all_pathways(PATHWAY_FILES, dss, column_to_merge)
# er gaat al wat mis voor merge_gsea_results want head all_pathways, geeft nog maar alleen voor AUR1 de resultaten.
#head(all_pathways)

#sig_pathways_aur <- filter_by_sig_level(all_pathways, 0.05, "AUR")
#head(sig_pathways_aur)
# # # head(sig_pathways_aur)
# # # nrow(sig_pathways_aur)
# # # # sig_pathways_aur[[2]]
# # # pathway_df <- all_pathways
# # # nrow(all_pathways)
# # # nrow(pathway_df)
# # # # sig_level <- 0.05
# # # # dataset <- "AUR"
# # # head(pathway_df)
# # # # class(colnames_sig[[1]])

#sig_pathways_cos <- filter_by_sig_level(all_pathways, 0.05, "COS")
# # #packageVersion("rlang")

# gene_sets <- list(
#   Discovery = sig_pathways_cos$pathway,
#   Validation = sig_pathways_aur$pathway
# )



# venn <- venn_plot(gene_sets, c("#984EA3", "#4DAF4A"))

fwrite(
    all_pathways,
    file = file.path(OUTPUT_DIR, "overlapping_pathways_all.tsv"),
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE
)

# fwrite(
#     sig_pathways_aur,
#     file = file.path(OUTPUT_DIR, "overlapping_pathways_sig_aurora.tsv"),
#     sep = "\t",
#     col.names = TRUE,
#     row.names = FALSE
# )
# pdf(file.path(OUTPUT_DIR, "overlapping_pathways_metastasis.pdf"), width = 10)
# venn
# dev.off()

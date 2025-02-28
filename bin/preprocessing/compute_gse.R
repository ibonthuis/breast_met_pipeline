### Loading libraries
required_libraries <- c(
    "data.table",    
    "dplyr",
    "optparse",
    "rlang",
    "ggplot2",
    "purrr",
    "fgsea",
    "RColorBrewer")

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
        metavar = "character"),
    optparse::make_option(
        c("-p", "--p_threshold"),
        type = "character",
        default = NULL,
        help = "Threshold value for significance in adjusted p-values for the enriched gene sets",
        metavar = "numeric"),
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
GENE_SET <- opt$gene_set_file
RANK_FILE <- opt$ranks
P_THRESH <- opt$p_threshold
OUTPUT_DIR <- opt$output_dir


## Debug
# GENE_SET <- "data/gene_sets/c2.cp.reactome.v5.0.symbols.gmt"
# RANK_FILE <- "test_again/differential_indegrees.rnk"
# P_THRESH <- 0.05
# OUTPUT_DIR <- "test"


## Functions
source("bin/preprocessing/compute_gse_fn.R")

## Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Gene set enrichment analysis
rank_list <- get(load(RANK_FILE))
column_to_merge <- "pathway" # possible variable to change when the headers are changing
dataset <- "AUR" # possible variable to change when the headers are changing
blues <- RColorBrewer::brewer.pal(n = 8, name = "Blues")
pathways <- gmtPathways(GENE_SET)

if (class(rank_list) == "list") {
   list_of_gsea <- purrr::map(rank_list, ~ perform_gsea(.x, pathways))
   selected_gsea <- purrr::map(list_of_gsea, ~ {
    df <- as.data.frame(.x)
    df <- df[order(df$padj, na.last = TRUE, decreasing = FALSE ),]
    df <- df[1:20, ]
    df
   })
   list_of_bubble <- purrr::map(selected_gsea, ~ plot_bubble_plot(.x, blues))
   list_for_tsv <- purrr::map(list_of_gsea, ~ {
    df <- as.data.frame(.x)
    df <- df[df$padj<P_THRESH, ]
    df
   })
  df_for_tsv <- merge_enriched_pathways_list_into_df(list_for_tsv, column_to_merge, dataset)

} else {
    list_of_gsea <- perform_gsea(rank_list, pathways)
    selected_gsea <- list_of_gsea[order(list_of_gsea$padj, na.last = TRUE, decreasing = FALSE),]
    selected_gsea <- selected_gsea[1:20,]
    list_of_bubble <- plot_bubble_plot(selected_gsea, blues)
    df_for_tsv <- list_of_gsea[list_of_gsea$padj<P_THRESH, ]
    list_of_gsea <- list(list_of_gsea)
}


save(
    list_of_gsea,
    file = file.path(OUTPUT_DIR, paste0("differential_genesets.RData")) # , P_THRESH,
)

data.table::fwrite(
    df_for_tsv,
    file = file.path(OUTPUT_DIR, "differential_genesets.tsv"),
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE
)

if (class(rank_list) == "list") {
    pdf(
        file.path(OUTPUT_DIR, paste0("enrichment_bubble_plot_.pdf")),
        #file.path(OUTPUT_DIR, paste0("enrichment_bubble_plot_", P_THRESH, ".pdf")),
        width = 18,
        height = 10
    )
    for (i in 1:length(list_of_bubble)) {
        show(list_of_bubble[[i]])
    }
    dev.off()
} else {
    pdf(
        file.path(OUTPUT_DIR, paste0("enrichment_bubble_plot_.pdf")),
        #file.path(OUTPUT_DIR, paste0("enrichment_bubble_plot_", P_THRESH, ".pdf")),
        width = 18,
        height = 10)
   show(list_of_bubble)
   dev.off()
}
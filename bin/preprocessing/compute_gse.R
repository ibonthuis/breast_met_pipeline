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
load(RANK_FILE)
rank_list <- ranked_paired_toptables
blues <- RColorBrewer::brewer.pal(n = 8, name = "Blues")
pathways <- gmtPathways(GENE_SET)

if (class(rank_list) == "list") {
   list_of_gsea <- purrr::map(rank_list, ~ perform_gsea(.x, pathways, P_THRESH))
   list_of_bubble <- purrr::map(list_of_gsea, ~ plot_bubble_plot(.x, blues))
} else {
    list_of_gsea <- perform_gsea(rank_list, GENE_SET, P_THRESH)
    list_of_bubble <- plot_bubble_plot(list_of_gsea, blues)
}


save(
    list_of_gsea,
    file = file.path(OUTPUT_DIR, paste0("differential_genesets", P_THRESH, ".RData"))
)


if (class(list_of_bubble) == "list") {
    pdf(
        file.path(OUTPUT_DIR, paste0("bubble_plot_padj", P_THRESH, ".pdf"),
        width = 16,
        height = 16)
    )
    for (i in 1:length(list_of_bubble)) {
       show(list_of_bubble[i])
    }
    dev.off()
} else {
    pdf(
        file.path(OUTPUT_DIR, paste0("bubble_plot_padj", P_THRESH, ".pdf")),
        width = 16,
        height = 16
    )
   show(list_of_bubble)
   dev.off()
}


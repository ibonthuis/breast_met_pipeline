### Loading libraries
required_libraries <- c(
    "data.table",    
    "dplyr",
    "optparse",
    "rlang",
    "ggplot2",
    "ggrepel")

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
        c("-v", "--visualization_var"),
        type = "character",
        default = NULL,
        help = "Variable to visualize, one of the column names of the metadata",
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
OUTPUT_DIR <- opt$output_dir
VARIABLE <- opt$visualization_var

## Debug
# INDEGREES_FILE <- "data/geo_brca/filtered_cos_indegree.csv"
# METADATA_FILE <- "data/geo_brca/metadata.csv"
# OUTPUT_DIR <- "test"


## Functions
source("bin/preprocessing/perform_dimensionality_reduction_fn.R")

## Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Dimensionality reduction
indegree_df <- read_indegree(INDEGREES_FILE)
metadata_df <- fread(METADATA_FILE)

if (ncol(indegree_df) < 100) {
    print(paste0("indegree_df has less than 100 columns"))
    # metadata_df <- metadata_df %>%
    #     filter(sample_type == "Metastasis")
}

# Comment out the following for pseudobulks.
metadata_df <- match_metadata_order_with_df(metadata_df, indegree_df)

## Plotting
standard_colors <- RColorBrewer::brewer.pal(8, "Set1")
the_theme <- ggplot2::theme_set(
    theme_bw() +
        theme(
            plot.title = element_text(size = 28),
            text = element_text(size = 26, colour = "black"),
            axis.title = element_text(size = 28, colour = "black"),
            axis.text = element_text(size = 18, colour = "black"),
            legend.title = element_text(size = 28, colour = "black"),
            legend.text = element_text(size = 26, colour = "black"),
            
        )
)
print(nrow(indegree_df))
indegree_df <- indegree_df[!rowSums(indegree_df) == 0,]
print(nrow(indegree_df))

#head(metadata_df)
# var <- ensym(VARIABLE)
# TOPRINT <- metadata_df$var
# print(class(TOPRINT))
print(VARIABLE)
#print(class(VARIABLE))

values <- metadata_df[[VARIABLE]]
print(class(values))
pca <- do_pca(indegree_df)


if(is.integer(values) == TRUE) {
   pcaplot <- pca_visualize_continuous(pca, metadata_df, VARIABLE)
} else {
   pcaplot <- plot_12_discrete(pca, metadata_df, VARIABLE)
}


## Exporting
pdf( 
    file.path(OUTPUT_DIR, 
    paste0("PCA_", VARIABLE, "_pc12.pdf"))
    )
pcaplot
dev.off()
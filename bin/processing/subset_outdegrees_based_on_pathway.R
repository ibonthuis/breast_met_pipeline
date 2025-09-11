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
    "tidyr",
    "ggpubr")

for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
}

### Options
options(stringsAsFactors = FALSE)

### Command line options
option_list <- list(

    optparse::make_option(
        c("-i", "--outdegrees"),
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
        c("-d", "--diff_outdegrees"),
        type = "character",
        default = NULL,
        help = "Path to the differential_outdegrees file.",
        metavar = "character"),
    # optparse::make_option(
    #     c("-e", "--edge_file"),
    #     type = "character",
    #     default = NULL,
    #     help = "Path to the edge file.",
    #     metavar = "character"),
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
METADATA_FILE <- opt$metadata
OUTDEGREES_FILE <- opt$outdegrees
DIFF_OUT_FILE <- opt$diff_outdegrees
GENE_SET_FILE <- opt$gene_set_file
PATHWAY <- opt$pathway_of_interest
OUTPUT_DIR <- opt$output_dir

## Debug
# METADATA_FILE <- "data/input_data/cosgrove/metadata.csv"
# OUTDEGREES_FILE <- "data/input_data/cosgrove/filtered_outdegrees.csv"
# DIFF_OUT_FILE <- "snakemake_results_with_outdegree/differential_outdegrees/cosgrove/differential_outdegrees.RData"
# GENE_SET_FILE <- "data/gene_sets/c2.cp.reactome.v5.0.symbols.gmt"
# PATHWAY <- "REACTOME_INTERFERON_SIGNALING"
# OUTPUT_DIR <- "testetst"
# dir.create(OUTPUT_DIR)
# For aurora
# OUTDEGREES_FILE <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/data/input_data/aurora/filtered_outdegrees.csv"
# METADATA_FILE <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/data/input_data/aurora/metadata.csv"

## Functions 
source("bin/processing/filter_analyze_pathway_specific_genes_fn.R")
source("bin/preprocessing/compute_diff_indegrees_fn.R")


## Computing
metadata <- fread(METADATA_FILE)
gene_set <- reading_gene_set_file(GENE_SET_FILE)
#genes <- select_genes_from_pathway_of_interest(gene_set, PATHWAY)

genes_of_pathway_of_interest <- gene_set[PATHWAY]
head(genes_of_pathway_of_interest)
outdegrees <- read_indegree(OUTDEGREES_FILE)
nrow(outdegrees)

#pathway_oi_genes <- select_genes_from_pathway_of_interest(gene_set, PATHWAY)
outdegrees_pathway_oi <- filter_indegrees_by_pathway_of_interest(outdegrees, genes_of_pathway_of_interest)
nrow(outdegrees_pathway_oi)

# head(rownames(outdegrees))

# outdegrees_pathway_oi <- outdegrees[rownames(outdegrees) %in% genes_of_pathway_of_interest, ]
# nrow(outdegrees_pathway_oi)
# which(genes_of_pathway_of_interest %in% rownames(outdegrees))

# intersect(genes_of_pathway_of_interest, rownames(outdegrees))

# I could create paired plots now
# Maybe good to know how many TFs we're talking about = 12 TFs in interferon signal
# head(metadata)
#differential_outdegrees <- create_toptable_paired(outdegrees_pathway_oi, metadata, "patient", "sample_type") 
# the above only works for cosgrove not for aurora
#head(differential_outdegrees)




# diff_outs <- get(load(DIFF_OUT_FILE))
# head(diff_outs)
# diff_outdegrees_pathway_oi <- filter_indegrees_by_pathway_of_interest(diff_outs, genes_of_pathway_of_interest)

fwrite(outdegrees_pathway_oi,
        file = file.path(OUTPUT_DIR, "outdegrees_pathway_oi.tsv"),
        row.names = TRUE,
        sep = "\t"
)

# pairedplots
# also this one doesn't really work for aurora data...........
# make unpaired plot?
pdf(file.path(OUTPUT_DIR, "outdegrees_pathway_paired_plot.pdf"))
for(i in 1:nrow(outdegrees_pathway_oi)) {
  tf <- rownames(outdegrees_pathway_oi)[i]
  idx <- which(row.names(outdegrees_pathway_oi)==tf)
  stat1_ind <- outdegrees_pathway_oi[rownames(outdegrees_pathway_oi) == tf,]
stat1_ind <- make_longer(stat1_ind, "outdegree")
  stat1_ind <- as.data.frame(stat1_ind)
  class(metadata)
  metadata <- as.data.frame(metadata)
  stat1_ <- merge(stat1_ind, metadata, by.x = 1, by.y = 2)
  test <- stat1_[, c("patient", "sample_type", "outdegree")]
  colnames(test) <- c("pat","type","outdegree")
  #test[which(test[,2]=="met"),2] <- "zmet"
  test <- as.data.frame(test)
  test[,1] <- as.factor(test[,1])
  test[,2] <- factor(test[,2], levels = c("Primary", "Metastasis"))
  test[,3] <- as.numeric(test[,3])
  pairplot <- ggpaired(test, x = "type", y = "outdegree",
           color = "type", line.color = "gray", line.size = 0.4,
           palette = c("lightsteelblue2","#ff69b4"), 
           id = "pat", 
           main=tf)+
           ylab("Outdegree")+
           theme(axis.title.y = element_text(size = 22), axis.title.x = element_blank(), axis.text.x = element_text(size = 22))+
    stat_compare_means(paired = TRUE)
    
  show(pairplot)    
}
dev.off()

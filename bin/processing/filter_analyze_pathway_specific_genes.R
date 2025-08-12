
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
INDEGREES_FILE <- "data/input_data/aurora/filtered_indegree.csv"
PERMUTATION_INDEGREES_FILE <- "subselection_diff_indegrees/selected_metprim_pairs_indegrees.RData"
METADATA_FILE <- "data/input_data/aurora/metadata.csv"
PATHWAY <- "REACTOME_INTERFERON_SIGNALING"
GENE_SET_FILE <- "data/gene_sets/c2.cp.reactome.v5.0.symbols.gmt"
#EDGE_FILE <- "/storage/kuijjerarea/ine/breast_met/breast_other/tissue_batch_corr_allAUR/output/lioness_filtered_for_genes.csv"
#EDGE_FILE <- "/storage/kuijjerarea/ine/breast_met/breast_other/tissue_batch_corr_allAUR/edge_table_slightly_adjusted.tsv"
EDGE_FILE <- "/storage/kuijjerarea/ine/breast_met/breast_other/tissue_batch_corr_allAUR/output_v2_samplenames_update/lioness_output_filtered_for_genes.csv"
EDGE_FILE <- "/storage/kuijjerarea/ine/breast_met/breast_other/R/data/edges_filtered_for_interferon_signallingpathway.tsv"
OUTPUT_DIR <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/subselection_diff_indegrees"

## Functions
source("bin/preprocessing/compute_diff_indegrees_fn.R")
source("bin/processing/filter_analyze_pathway_specific_genes_fn.R")


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
metadata <- fread(METADATA_FILE)
gene_set <- reading_gene_set_file(GENE_SET_FILE)
genes <- select_genes_from_pathway_of_interest(gene_set, PATHWAY)


indegrees_filtered_by_pathway <- purrr::map(indegrees_paired, ~ filter_indegree_df(GENE_SET_FILE, PATHWAY, .x))
nrow(indegrees_filtered_by_pathway[[1]])

filtered_diff_indegrees <- purrr::map(indegrees_filtered_by_pathway, ~ create_toptable_paired(.x, metadata, "patient", "sample_type"))
head(filtered_diff_indegrees[[1]])

#filtered_diff_edges <- filter_edges(EDGE_FILE, METADATA_FILE, PATHWAY, 1.5)

# in case of permutations:

# read indegree file into df (or list of indegrees)
indegrees <- read_indegree(INDEGREES_FILE)
colnames(indegrees) %in% colnames(edges_df)
colnames(indegrees)[51]
colnames(indegrees)[52]

colnames(edges_df) %in% colnames(indegrees) 

indegrees_paired <- read_indegrees_from_permutation(PERMUTATION_INDEGREES_FILE)
head(indegrees_paired[[1]])


selected_columns_for_paired_edge_analysis <- purrr::map(indegrees_paired, ~ select_column_names(.x))
edges_df <- fread(EDGE_FILE)
colnames(edges_df)

edges_df[1:5, 1:5]
nrow(edges_df)

# for aurora edges
filtered_edges <- filter_edges_by_pathway_of_interest(edges_df, GENE_SET_FILE, PATHWAY)
nrow(edges_df)
head(filtered_edges)
numeric_values <- filtered_edges %>%
  select(-TF, -Target) %>%  # Exclude TF and Target columns
  unlist()  


# for cosgrove edges
numeric_values <- filtered_edges %>%
  select(-reg, -tar) %>%  # Exclude TF and Target columns
  unlist()  
head(numeric_values)
standard_colors <- RColorBrewer::brewer.pal(8, "Set1")

colors_more_than_8 <- viridis::plasma(14)
colors_12 <- colors_more_than_8[1:12]

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
    
  ggplot(numeric_values, aes(x = logFC)) +
  geom_histogram(aes(y=after_stat(density)), binwidth = 0.1, fill = "grey") +
  labs(
    title = "Validation dataset",
    x = "logFC edge weights",
    y = "Density"
  ) 

ggplot(data.frame(value = numeric_values), aes(x = value)) +
  geom_histogram(aes(y=after_stat(density)), binwidth = 0.1, fill = "grey") +
  labs(
    title = "Interferon pathway edges discovery dataset",
    x = "Edge weights",
    y = "Density"
  ) +
  geom_density(alpha = 0.2, fill = "#FF6666") +
  theme(title = element_text(size = 14)) #+
  #xlim(-10, 10)


ggsave("/storage/kuijjerarea/ine/breast_met/breast_other/R/results/cosgrove/edges_plots/histogram_ifn_filtered_edges_cosgrove_densitytoo.png",
    dpi = 300)

pdf("/storage/kuijjerarea/ine/breast_met/breast_other/R/results/cosgrove/edges_plots/histogram_ifn_filtered_edges_cosgrove_densitytoo.pdf")
ggplot(data.frame(value = numeric_values), aes(x = value)) +
  geom_histogram(aes(y=after_stat(density)), binwidth = 0.1, fill = "grey") +
  labs(
    title = "Interferon pathway edges discovery dataset",
    x = "Edge weights",
    y = "Density"
  ) +
  geom_density(alpha = 0.2, fill = "#FF6666") +
  theme(title = element_text(size = 14)) #+
  #xlim(-10, 10)
dev.off()

prepped_edge_df <- prep_edge_df_for_differential_analysis(filtered_edges)
selected_edge_columns <- purrr::map(selected_columns_for_paired_edge_analysis, ~ select_edges_for_permutations(prepped_edge_df, .x))


#fwrite(filtered_edges_df, "/storage/kuijjerarea/ine/breast_met/breast_other/R/data/interferon_signaling_edges_aurora_data_target_filtered.tsv", sep = "\t")


## NEED TO RUN FILTERED_EDGES AGAIN FROM LINE 119


# function from compute_diff_indegrees
paired_edge_toptables <- purrr::map(selected_edge_columns, ~ create_toptable_paired(.x, metadata, "patient", "sample_type"))
#paired_edge_toptables <- purrr::map(paired_edge_toptables, ~ split_edge_df_rownames(.x))

nrow(paired_edge_toptables[[1]])

save(
    paired_edge_toptables,
    file = file.path(OUTPUT_DIR, "differential_edges_aurora.RData"))

save(
    filtered_diff_indegrees,
    file = file.path(OUTPUT_DIR, "diff_indegrees_pathway_specific_genes.RData")
)

for (i in 1:length(paired_edge_toptables)) {
   fwrite(paired_edge_toptables[[i]], paste0("/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/interferon_signaling_diff_aurora_edges/diff_edges_aur", i, ".tsv"), row.names = TRUE, sep = "\t")
}


edges_aur1 <- fread("/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/interferon_signaling_diff_aurora_edges/diff_edges_aur1.tsv")
top_to_select <- fread("/storage/kuijjerarea/ine/breast_met/breast_other/R/results/cosgrove/limma_analysis/interferon_signalling/differential_top150_edges_logfc1.5.tsv")

head(top_to_select)

top_to_select <- prep_edge_df_for_differential_analysis(top_to_select)

edges_aur1_fil <- edges_aur1 %>%
    filter(V1 %in% rownames(top_to_select))
nrow(edges_aur1_fil)

tail(edges_aur1_fil, n = 100)
edges_aur1_fil <- as.data.frame(edges_aur1_fil)
rownames(edges_aur1_fil) <- edges_aur1_fil$V1

edges_aur1_filt <- split_edge_df_rownames(edges_aur1_fil)
head(edges_aur1_filt)
edges_aur1_filt$V1 <- NULL


fwrite(edges_aur1_filt, "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/subselection_diff_indegrees/differential_top150_edges_logfc1.5_fromcosgrove_selected_aurora.tsv", sep = "\t")

### Loading libraries
required_libraries <- c(
    "data.table",    
    "dplyr",
    "optparse",
    "rlang",
    "ggplot2",
    "purrr",
    "fgsea",
    "RColorBrewer",
    "stringr",
    "ComplexHeatmap")

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
        help = "Path to the ranked file. But snakemake will know this because it was built upstream",
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
GENE_SET <- "data/gene_sets/c2.cp.reactome.v5.0.symbols.gmt"
RANK_FILE <- "snakemake_results/differential_indegrees/cosgrove/differential_indegrees_rank_file.RData"
P_THRESH <- 0.05
OUTPUT_DIR <- "test_jaccard_index_v1"

# For covariate corrected ranked genes of cosgrove
RANK_FILE <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/corr_1/ranked_corrected_tumor_purity_differential_indegrees.RData"
OUTPUT_DIR <- "corr_1"
setwd("/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline")

## Functions
source("bin/preprocessing/compute_gse_fn.R")
source("bin/preprocessing/test_JI_stuff.R")

## Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Gene set enrichment analysis
rank_list <- get(load(RANK_FILE))
#rank_list <- fread(RANK_FILE)
column_to_merge <- "pathway" # possible variable to change when the headers are changing
dataset <- "AUR" # possible variable to change when the headers are changing
blues <- RColorBrewer::brewer.pal(n = 8, name = "Blues")
pathways <- gmtPathways(GENE_SET)

head(rank_list)
head(pathways)
filtered_genes_in_pathways <- purrr::map(pathways, ~ .x[.x %in% rank_list$gene])
list_of_gsea <- perform_gsea(rank_list, filtered_genes_in_pathways)

#filtered_genes_in_pathways_df <- as.data.frame(filtered_genes_in_pathways)

# The following lines of code are only to be run if the output doesn't contain a column with pathway genes
# there is really no use of this code since I discovered that you have to have leading edge genes
# which have a certain order that determines which kind of help in prioritizing pathways somehow
pathways_inter <- purrr::map(filtered_genes_in_pathways, ~ stringr::str_flatten(.x, ";"))
pathways_df <- map_df(names(pathways_inter), ~ data.frame(
  pathway_name = .x,
  gene = pathways_inter[[.x]]
))
ready_for_ji <- merge(list_of_gsea, pathways_df, by.x = 1, by.y = 1)
head(ready_for_ji)
ji_cluster <- cluster_gsea_enrichment(ready_for_ji, "gene", "pathway")
sig_for_ji <- ready_for_ji %>%
    filter(pval < 0.05) %>%
    filter(size > 10)
nrow(sig_for_ji)
head(pathways_df)

# some self checks
head(list_of_gsea)

list_of_gsea[list_of_gsea$pathway == "REACTOME_INTERFERON_SIGNALING",]
length(list_of_gsea$pathway)

class(list_of_gsea)
list_of_gsea <- as.data.frame(list_of_gsea)

head(ji_cluster)


sig_for_ji <- list_of_gsea %>%
    filter(pval < 0.05) %>%
    filter(size > 10)
nrow(sig_for_ji) #180 if pval < 0.05, 110 with padj < 0.05
head(sig_for_ji) # 234 for covariate corrected gsea

sig_for_ji <- sig_for_ji %>%
    filter(padj < 0.05)
nrow(sig_for_ji)

# try collapsePathways() and follow JI
head(rank_list)
class(rank_list)
colnames(rank_list) <- c("geneName", "stat")
ranks <- tibble::deframe(rank_list)
class(ranks)
head(ranks)
str(ranks)
collapsed_gsea <- collapsePathways(sig_for_ji, pathways = filtered_genes_in_pathways, stats = ranks)
nrow(collapsed_gsea)

collapsed_gsea$mainPathways
collapsed_gsea$parentPathways

sig_for_ji_filt <- sig_for_ji[sig_for_ji$pathway %in% collapsed_gsea$mainPathways, ]

sig_for_ji_filt <- sig_for_ji_filt %>% 
    dplyr::arrange(padj, )


#hmap_outpath <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/test_jaccard_index_v1/heatmap_leadingEdge_padj0.05.png"
hmap_outpath <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/test_jaccard_index_v1/heatmap_collapse_jaccard_fdr.png"

# change sig_for_ji
ji_cluster_sig <- cluster_gsea_enrichment(sig_for_ji_filt, "leadingEdge", "pathway", lead_genes_split = ',',
                    jaccard_hclust_cuts = c(0.5, 1, 1.2, 1.5), hmap_outpath = hmap_outpath, 
                    hmap_title = "metastasis_vs_primary")


head(ji_cluster_sig)
class(ji_cluster_sig)
head(sig_cl_hmap)
length(unique(ji_cluster_sig$path_cluster_cut_18))
nrow(ji_cluster_sig)

summarizing_pathway_12 <- best_pathway_clust(ji_cluster_sig, "path_cluster_cut_12")
head(summarizing_pathway_12)
nrow(summarizing_pathway_12)

summarizing_pathway_15 <- best_pathway_clust(ji_cluster_sig, "path_cluster_cut_15")
head(summarizing_pathway_15)
nrow(summarizing_pathway_15)

summarizing_pathway_18 <- best_pathway_clust(ji_cluster_sig, "path_cluster_cut_18")
head(summarizing_pathway_18$pathway)
nrow(summarizing_pathway_18)

summarizing_pathway_12[8, 1] <- "RESPIRATORY ELECTRON TRANSPORT ATP SYNTHESIS \nBY CHEMIOSMOTIC COUPLING AND HEAT PRODUCTION BY UNCOUPLING PROTEINS"
summarizing_pathway_12$pathway <- gsub("_", " ", summarizing_pathway_12$pathway)
summarizing_pathway_12$pathway <- gsub("REACTOME", "", summarizing_pathway_12$pathway)

#   summarizing_pathway_12 <- summarizing_pathway_12 %>%
#                   mutate(log_padj = log(padj)) %>%
#                   dplyr::arrange(log_padj, )
  #gsea_result <- gsea_result[1:nr_of_pathways, ]


  gsea_result <- summarizing_pathway_12 %>%
                  mutate(log_padj = -log(padj)) %>%
                  dplyr::arrange(desc(log_padj), )

list_of_bubble <- plot_bubble_plot(summarizing_pathway_12, blues, 10)

summarizing_pathway_15$pathway <- gsub("_", " ", summarizing_pathway_15$pathway)
summarizing_pathway_15$pathway <- gsub("REACTOME", "", summarizing_pathway_15$pathway)


summarizing_pathway_15[10, 1] <- "RESPIRATORY ELECTRON TRANSPORT ATP SYNTHESIS \nBY CHEMIOSMOTIC COUPLING AND HEAT PRODUCTION BY UNCOUPLING PROTEINS"
list_of_bubble <- plot_bubble_plot(summarizing_pathway_15, blues, 10)

# For covariate adjusted:
summarizing_pathway_18$pathway <- gsub("_", " ", summarizing_pathway_18$pathway)
summarizing_pathway_18$pathway <- gsub("REACTOME", "", summarizing_pathway_18$pathway)

summarizing_pathway_18[12, 1] <- "RESPIRATORY ELECTRON TRANSPORT ATP SYNTHESIS \nBY CHEMIOSMOTIC COUPLING AND HEAT PRODUCTION BY UNCOUPLING PROTEINS"
list_of_bubble <- plot_bubble_plot(summarizing_pathway_18, blues, 10)


clusters <- ji_cluster_sig %>%
    arrange(path_cluster_cut_12)


ncol(clusters)
clusters <- clusters[, c(1, 8, 11)]
head(clusters)

clusters$pathway <- gsub("_", " ", clusters$pathway)
clusters$pathway <- gsub("REACTOME", "", clusters$pathway)

cluster <- as.data.frame(clusters[, 'pathway'])

fwrite(cluster, file.path(OUTPUT_DIR, "table_clusters_12_just_pathway_names.tsv"), sep = "\t")



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

class(rank_list)

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
        file.path(OUTPUT_DIR, paste0("enrichment_collapse_JI_bubble_plot_ji12_top10.pdf")),
       # file.path(OUTPUT_DIR, paste0("enrichment_bubble_plot_", P_THRESH, ".pdf")),
        width = 12,
        height = 8)
   show(list_of_bubble)
   dev.off()
}

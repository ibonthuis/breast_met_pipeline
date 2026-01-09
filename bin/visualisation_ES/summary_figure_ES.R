required_libraries <- c(
    "data.table",    
    "dplyr",
    #"optparse",
    "rlang",
    "purrr",
    "stringr",
    "limma",
    #"ggvenn",
    #"ggpubr",
    "ggplot2",
    "ggrepel")

for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
}

vignette(package = "limma")


## Functions
#source("bin/processing/overlap_pathways_fn.R")
## Theme
standard_colors <- RColorBrewer::brewer.pal(8, "Set1")
color_scheme <- RColorBrewer::brewer.pal(9, "Set1")
the_theme <- ggplot2::theme_set(
    theme_bw() +
        theme(
            plot.title = element_text(size = 28),
            text = element_text(size = 30, colour = "black"),
            axis.title = element_text(size = 30, colour = "black"),
            axis.text = element_text(size = 18, colour = "black"),
            legend.title = element_text(size = 28, colour = "black"),
            legend.text = element_text(size = 26, colour = "black"),
            
        )
)


OUTPUT_DIR <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/results_summary_20251219"

GSEA_COSGROVE <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results_05112025/differential_indegrees/cosgrove/differential_genesets.tsv"
GSEA_AURORA <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results_05112025/differential_indegrees/aurora/differential_genesets.tsv"
GSEA_COSGROVE_EXP <- "snakemake_results_expression_20251203/differential_indegrees/cosgrove/differential_genesets.tsv"
GSEA_AURORA_EXP <- "snakemake_results_expression_20251203/differential_indegrees/aurora/differential_genesets.tsv"
GSEA_EPITHELIAL <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results_Epithelial_Myeloid_IFNpathway/differential_indegrees/Epithelial/differential_genesets.tsv"
# GSEA_EPITHELIAL_EXP
GSEA_T <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results_T_Epithelial_IFNpathway/differential_indegrees/T/differential_genesets.tsv"
GSEA_MYELOID <- "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results_Epithelial_Myeloid_IFNpathway/differential_indegrees/Myeloid/differential_genesets.tsv"
# GSEA_T_EXP
# GSEA_MYELOID
# GSEA_MYELOID_EXP

# Define the file paths and their corresponding variable names
file_paths <- list(
  gsea_cosgrove = "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results_05112025/differential_indegrees/cosgrove/differential_genesets.tsv",
  gsea_aurora = "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results_05112025/differential_indegrees/aurora/differential_genesets.tsv",
  gsea_cosgrove_exp = "snakemake_results_expression_20251203/differential_indegrees/cosgrove/differential_genesets.tsv",
  gsea_aurora_exp = "snakemake_results_expression_20251203/differential_indegrees/aurora/differential_genesets.tsv",
  gsea_epithelial = "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results_Epithelial_Myeloid_IFNpathway/differential_indegrees/Epithelial/differential_genesets.tsv",
  gsea_t = "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results_T_Epithelial_IFNpathway/differential_indegrees/T/differential_genesets.tsv",
  gsea_myeloid = "/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results_Epithelial_Myeloid_IFNpathway/differential_indegrees/Myeloid/differential_genesets.tsv"
)

# Read all files into data.tables and assign them to variables dynamically
for (name in names(file_paths)) {
  assign(name, fread(file_paths[[name]]))
}

# Check the loaded data.tables
ls(pattern = "^gsea_")



PSEUDOBULK_NORM <- "/storage/kuijjerarea/ine/projects/BRCA_SINGLECELL_TO_PSEUDOBULK/celltype_specific_pseudobulks/snakemake_results_myeloid/myeloid/pseudobulk_normalized.tsv"


pseudo_norm <- fread(PSEUDOBULK_NORM)
head(pseudo_norm)
pseudo_norm <- as.data.frame(pseudo_norm)
rownames(pseudo_norm) <- pseudo_norm$genes
pseudo_norm$genes <- NULL

histnorm <- ggplot(pseudo_norm, aes(patient10_Primary_Myeloid)) +
    geom_histogram(bins = 50)+
    xlim(NA, 20)
histnorm
ncol(pseudo_norm)

PSEUDOBULK <- "/storage/kuijjerarea/ine/projects/BRCA_SINGLECELL_TO_PSEUDOBULK/celltype_specific_pseudobulks/snakemake_results_50cells/Myeloid/pseudobulk_normalized.tsv"
pseudobulk <- fread(PSEUDOBULK)
head(pseudobulk)

pseudobulk <- as.data.frame(pseudobulk)
rownames(pseudobulk) <- pseudobulk$genes
pseudobulk$genes <- NULL
ncol(pseudobulk) # 31

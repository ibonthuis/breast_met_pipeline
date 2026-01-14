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

#vignette(package = "limma")


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




# Read all files into data.tables and assign them to variables dynamically
for (name in names(file_paths)) {
  assign(name, fread(file_paths[[name]]))
}

# Check the loaded data.tables
ls(pattern = "^gsea_")

length(ls(pattern = "^gsea_"))

rowname <- paste0(ls(pattern = "^gsea_")[1])
df <- get(rowname)
gsea_aurora_exp_rowname <- paste0(ls(pattern = "^gsea_")[2])
gsea_aurora_exp <- gsea_aurora_exp[gsea_aurora_exp$pathway == "REACTOME_INTERFERON_SIGNALING", ]
df <- df[df$pathway == "REACTOME_INTERFERON_SIGNALING", ]


head(gsea_t)


head(df)

all_gsea <- data.frame(dataset_analysis = c(rowname, gsea_aurora_exp_rowname), 
                        padj = c(df$PADJ_AUR1, gsea_aurora_exp$PADJ_AUR1),
                        es = c(df$ES_AUR1, gsea_aurora_exp$ES_AUR1))
head(all_gsea)

for (i in 4:length(ls(pattern = "^gsea_"))) {
  print(i)
  rowname <- print(paste0(ls(pattern = "^gsea_")[i]))
  df <- get(rowname)
  if(nrow(df) > 1) {
    df <- df[df$pathway == "REACTOME_INTERFERON_SIGNALING", ]
  }
  df$pathway <- rowname
  row <- df[, c("pathway", "padj", "ES")]
  colnames(row) <- c("dataset_analysis", "padj", "es")
  all_gsea <- rbind(all_gsea, row)
 # print(head(all_gsea))
}

head(all_gsea)
all_gsea
#regular_expression <- paste0("(FDR|PADJ).*", dataset, "|", dataset, ".*(FDR|PADJ)")

all_gsea <- all_gsea %>%
  mutate(sig = case_when(padj > 0.1 ~ "ns",
                            TRUE ~ "significant"))

all_gsea$dataset_analysis <- c("Validation", 
                              "Validation expression", 
                              "Discovery", 
                              "Discovery expression",
                              "Epithelial cells",
                              "Epithelial cells expression",
                              "Myeloid cells",
                              "Myeloid cells expression",
                              "T cells",
                              "T cells expression")
class(all_gsea$dataset_analysis)

all_gse_plot <- ggplot(all_gsea, aes(x = dataset_analysis, y = es, fill = sig)) +
                  geom_col()+
                  coord_flip()+
                  ylab("ES")+
                  scale_fill_manual(values = c("grey60", "#377EB8"))+
#                  labs(sig = "significant")+
                  theme(axis.title.y = element_blank())
all_gse_plot              

standard_colors

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

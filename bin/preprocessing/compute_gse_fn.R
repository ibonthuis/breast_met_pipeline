#' Perform Gene Set Enrichment Analysis (GSEA)
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) on differential expression results.
#'
#' @param diff_results A data frame containing differential expression results. 
#'        It should include columns for gene identifiers and their associated statistics (e.g., log fold change, p-values).
#' @param gene_set A list or data frame containing the gene sets to be tested. 
#'        Each gene set should be a vector of gene identifiers.
#' @param p_threshold A numeric value specifying the p-value threshold for significance.
#'
#' @return A data frame with the results of the GSEA, including enrichment scores and p-values for each gene set.
#'
#' @examples
#' # Example usage:
#' diff_results <- data.frame(gene = c("gene1", "gene2", "gene3"), logFC = c(2.3, -1.5, 0.8), p_value = c(0.01, 0.05, 0.2))
#' gene_set <- list(set1 = c("gene1", "gene3"), set2 = c("gene2"))
#' p_threshold <- 0.05
#' gsea_results <- perform_gsea(diff_results, gene_set, p_threshold)
#'
#' @export
perform_gsea <- function(diff_results, gene_set) {
  colnames(diff_results) <- c("geneName", "stat")
  ranks <- tibble::deframe(diff_results)
  set.seed(1)
  fgseaRes <- fgseaMultilevel(pathways=gene_set, stats=ranks, maxSize=1000)
 # fgseaRes <- fgseaRes[, c(1:7)]
 # fgseaRes <- fgseaRes[fgseaRes$padj<p_threshold,] 
  return(fgseaRes)
}


#' Plot Bubble Plot for GSEA Results
#'
#' This function generates a bubble plot to visualize the results of Gene Set Enrichment Analysis (GSEA).
#'
#' @param gsea_result A data frame containing the results of GSEA. It should include columns for 
#' gene set name (NAME), enrichment scores (ES), p-values (padj), and other relevant metrics.
#' @param blues A color palette consisting of multiple colours, preferably getting increasingly 
#' darker. 
#'
#' @return A ggplot2 object representing the bubble plot of the GSEA results.
#'
#' @examples
#' \dontrun{
#'   # Assuming gsea_result is a data frame with the required structure
#'   bubble_plot <- plot_bubble_plot(gsea_result)
#'   print(bubble_plot)
#' }
#'
#' @export
plot_bubble_plot <- function(gsea_result, blues, nr_of_pathways) {
  gsea_result <- as.data.frame(gsea_result)
  # gsea_result <- gsea_result %>%
  #                   mutate(log_padj = log(padj))
  gsea_result <- gsea_result %>%
                  mutate(log_padj = -log(padj)) %>%
                  dplyr::arrange(desc(log_padj), )
  gsea_result <- gsea_result[1:nr_of_pathways, ]
    gsea_result <- gsea_result %>%
                  dplyr::arrange(log_padj, )
  gsea_result$pathway <- factor(gsea_result$pathway, levels = gsea_result$pathway)
  #gsea_result <- gsea_result[1:nr_of_pathways, ]

  g <- ggplot(gsea_result, aes(x=log_padj,y=pathway,size=size, fill = ES))+
          geom_point(shape=21, fg="grey")+ 
          scale_size_area(max_size=10)+
          scale_fill_gradient2(high = "red", mid = "white", low = blues[length(blues)])+
          theme(text = element_text(size = 20), axis.text.y = element_text(size = 14), axis.title.y = element_blank()
          #axis.title.x = element_text(size = 16)
          )+
          xlab("-log(FDR)")
  return(g)
}


merge_enriched_pathways_list_into_df <- function(gsea_results_list, column_to_merge, ds) {  
  gsea_results_list <- purrr::map2(gsea_results_list, 1:length(gsea_results_list), function(df, i) {
    df <- as.data.frame(df)
    df <- df[, c(1, 3, 5, 7)]
    colnames(df) <- c(column_to_merge, paste0("PADJ_", ds, i), paste0("ES_", ds, i), paste0("SIZE_", ds, i) )
    return(df)
  })

  merged_gsea_results <- purrr::reduce(gsea_results_list, function(df1, df2) {
    merge(df1, df2, by = column_to_merge)
  })
  
  return(merged_gsea_results)
}
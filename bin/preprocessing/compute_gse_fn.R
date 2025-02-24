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
perform_gsea <- function(diff_results, gene_set, p_threshold) {
  colnames(diff_results) <- c("geneName", "stat")
  ranks <- tibble::deframe(diff_results)
  set.seed(1)
  fgseaRes <- fgseaMultilevel(pathways=gene_set, stats=ranks, maxSize=1000)
  fgseaRes <- fgseaRes[, c(1:7)]
  fgseaRes <- fgseaRes[fgseaRes$padj<p_threshold,] 
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
plot_bubble_plot <- function(gsea_result, blues) {
  gsea_result <- as.data.frame(gsea_result)
  g <- ggplot(gsea_result, aes(x=padj,y=pathway,size=size, fill = ES))+
          geom_point(shape=21, fg="grey")+ 
          scale_size_area(max_size=10)+
          scale_fill_gradient2(high = "red", mid = "white", low = blues[length(blues)])+
          theme(axis.text.y = element_text(size = 14), axis.title.y = element_blank(), )
  return(g)
}

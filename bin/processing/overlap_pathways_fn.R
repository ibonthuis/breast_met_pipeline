#' Merge GSEA Results
#'
#' This function merges multiple GSEA (Gene Set Enrichment Analysis) result files into a single data frame.
#'
#' @param gsea_results_multi A list of file paths to GSEA result files.
#' @param ds A string representing the dataset identifier.
#' @param column_to_merge A string representing a column to merge the data by.
#'
#' @return A data frame containing the merged GSEA results.
#'
#' @importFrom purrr map map2 reduce
#' @importFrom data.table fread
#'
#' @examples
#' gsea_results_multi <- list("path/to/result1.txt", "path/to/result2.txt")
#' ds <- "dataset1"
#' merged_results <- merge_gsea_results(gsea_results_multi, ds, "pathway")
#' 
#' @export
merge_gsea_results <- function(gsea_results_multi, ds, column_to_merge) {
  gsea_results_list <- map(gsea_results_multi[1:length(gsea_results_multi)], function(ls) {
    get(load(ls))
  })
  gsea_results_list <- gsea_results_list[[1]]
  gsea_results_list <- map2(gsea_results_list, 1:length(gsea_results_list), function(df, i) {
    df <- as.data.frame(df)
    df <- df[, c(1, 3, 5, 7)]
    colnames(df) <- c(column_to_merge, paste0("PADJ_", ds, i), paste0("ES_", ds, i), paste0("SIZE_", ds, i) )
    return(df)
  })
  
  if (length(gsea_results_list) < 2) {
    merged_gsea_results_df <- gsea_results_list[[1]]
  } else {
    merged_gsea_results_df <- reduce(gsea_results_list, function(df1, df2) {
      merge(df1, df2, by = column_to_merge)
    })
  }
  return(merged_gsea_results_df)
}



merge_all_pathways <- function(pathway_files, list_of_datasets, column_to_merge) {
  gsea_results_list <- purrr::map2(pathway_files, list_of_datasets, function(file, dataset) {
    merge_gsea_results(file, dataset, column_to_merge)
  })

  merged_gsea_results_df <- purrr::reduce(gsea_results_list, function(df1, df2) {
    merge(df1, df2, by = column_to_merge)
  })
  
  return(merged_gsea_results_df)
}


# #' Merge All Pathways
# #'
# #' This function merges pathway results from two different sources. It first checks if the inputs are lists and, if so, processes them using the `merge_gsea_results` function. Then, it merges the processed results across datasets using the `merge_pathway_results_across_datasets` function.
# #'
# #' @param pathway_results1 A list or data frame containing pathway results from the first source.
# #' @param pathway_results2 A list or data frame containing pathway results from the second source.
# #' @param ds1 A character string specifying the first dataset.
# #' @param ds2 A character string specifying the second dataset.
# #'
# #' @return A data frame containing the merged pathway results from both sources.
# #'
# #' @examples
# #' \dontrun{
# #' pathway_results1 <- list() # Replace with actual data
# #' pathway_results2 <- list() # Replace with actual data
# #' merged_pathways <- merge_all_pathways(pathway_results1, pathway_results2)
# #' }
# #' 
# #' @export
# merge_all_pathways <- function(pathway_results1, pathway_results2, ds1, ds2, column_to_merge) {
#     # Check if pathway_results1 and pathway_results2 are lists
#     if (length(pathway_results1) > 1) {
#         pathway_results1 <- merge_gsea_results(pathway_results1, ds1)
#     }

#     if (length(pathway_results2) > 1) {
#         pathway_results2 <- merge_gsea_results(pathway_results2, ds2)
#     }
    
#     # Merge the results across datasets
#     all_pathways <- merge_pathway_results_across_datasets(pathway_results1, pathway_results2, column_to_merge)
    
#     return(all_pathways)
# }

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
    merged_gsea_results_df <- purrr::reduce(gsea_results_list, function(df1, df2) {
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


filter_by_sig_level_aur <- function(pathway_df, sig_level, dataset){
  # pathway_df contains all pathways (so nrow 640)
  #colnames_sig <- colnames(pathway_df)[grep("FDR|PADJ", colnames(pathway_df))]
  # dataset either "AUR" or "COS"

  regular_expression <- paste0("(FDR|PADJ).*", dataset, "|", dataset, ".*(FDR|PADJ)")
  colnames_sig <- colnames(pathway_df)[grep(regular_expression, colnames(pathway_df))]
# colnames_sig <- rlang::syms(colnames_sig)
  col1 <- data_sym(colnames_sig[[1]])
  col2 <- data_sym(colnames_sig[[2]])
  col3 <- data_sym(colnames_sig[[3]])
  col4 <- data_sym(colnames_sig[[4]]) 
  col5 <- data_sym(colnames_sig[[5]])
  pathway_df <- pathway_df %>%
   dplyr::filter(col1 <= sig_level &  col2 <= sig_level & col3 <= sig_level & col4 <= sig_level & col5 <= sig_level) %>%
   #dplyr::filter(!!data_sym(colnames_sig[[1]]) <= sig_level)
    dplyr::filter(ES_AUR1 < 0)
    return(pathway_df)
}


filter_by_sig_level_cos <- function(pathway_df, sig_level, dataset) {
  regular_expression <- paste0("(FDR|PADJ).*", dataset, "|", dataset, ".*(FDR|PADJ)")
  colnames_sig <- colnames(pathway_df)[grep(regular_expression, colnames(pathway_df))]
  pathway_df <- pathway_df %>%
    filter(!!rlang::sym(colnames_sig[1]) <= sig_level) %>%
    filter(ES_COS1 < 0)
  return(pathway_df)
}

filter_by_sig_level <- function(pathway_df, sig_level, dataset){
  if(dataset == "AUR"){
    filtered_df <- filter_by_sig_level_aur(pathway_df, sig_level, dataset)
  } else {
    filtered_df <- filter_by_sig_level_cos(pathway_df, sig_level, dataset)
  }
  return(filtered_df)
}

venn_plot <- function(gene_sets, colors) {
 venn <- ggvenn(gene_sets, c(names(gene_sets)[1], names(gene_sets)[2]), 
        fill_color = c(colors[1], colors[2]),
        set_name_size = 8,
        text_size = 12,
        show_percentage = FALSE,
        auto_scale = TRUE
        )
  return(venn)
}

# !!rlang::sym(colnames_sig[1]) <= sig_level & !!rlang::sym(colnames_sig[2]) <= sig_level & !!rlang::sym(colnames_sig[3]) <= sig_level & !!rlang::sym(colnames_sig[4]) <= sig_level & !!rlang::sym(colnames_sig[5]) <= sig_level & !!rlang::sym(colnames_sig[6]) <= sig_level ~ "all",
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



#' Generate Vector of Column Names that Match a Specified Input String
#'
#' This function takes a data frame of merged pathways and a string pattern to match,
#' and returns a list of column names that match the given string pattern.
#'
#' @param merged_pathways A data frame containing merged pathways.
#' @param string_match A string pattern to match column names with theagainst. For example, PADJ or FDR.
#'
#' @return A character vector of column names that match the given string pattern.
#'
#' @examples
#' merged_pathways <- data.frame(
#'   pathway1 = c(1, 2, 3),
#'   pathway2_sig = c(4, 5, 6),
#'   pathway3 = c(7, 8, 9),
#'   pathway4_sig = c(10, 11, 12)
#' )
#' string_match <- "_sig"
#' generate_list_significance_colnames(merged_pathways, string_match)
#' # Returns: "pathway2_sig" "pathway4_sig"
#'
#' @importFrom stringr str_subset
#' @export
get_significance_colnames <- function(merged_pathways, string_match) {
  all_cols <- colnames(merged_pathways)
  colnames_sig <- stringr::str_subset(all_cols, string_match)
  return(colnames_sig)
}



#' Categorize Significance of Pathways
#'
#' This function categorizes the significance of pathways based on specified significance levels.
#' It takes a data frame of merged pathways and a vector of column names representing significance values,
#' and assigns a category to each pathway based on the provided significance level.
#'
#' @param merged_pathways A data frame containing merged pathways with significance values.
#' @param colnames_sig A character vector of column names in `merged_pathways` representing significance values.
#' @param sig_level A numeric value representing the significance threshold.
#'
#' @return A data frame with an additional column `enriched` categorizing the significance of each pathway.
#' The categories are:
#' \itemize{
#'   \item "all": All significance values are below or equal to the threshold.
#'   \item "validation dataset only": Only the first significance value is above the threshold, while others are below or equal to it.
#'   \item "discovery dataset only": Only the first significance value is below or equal to the threshold, while others are above it.
#'   \item "none": All significance values are above the threshold.
#'   \item "none": Default category if none of the above conditions are met.
#' }
#'
#' @examples
#' \dontrun{
#' merged_pathways <- data.frame(
#'   pathway = c("pathway1", "pathway2"),
#'   sig1 = c(0.01, 0.05),
#'   sig2 = c(0.02, 0.06),
#'   sig3 = c(0.03, 0.07),
#'   sig4 = c(0.04, 0.08),
#'   sig5 = c(0.05, 0.09),
#'   sig6 = c(0.06, 0.10)
#' )
#' colnames_sig <- c("sig1", "sig2", "sig3", "sig4", "sig5", "sig6")
#' sig_level <- 0.05
#' categorized_pathways <- categorize_significance(merged_pathways, colnames_sig, sig_level)
#' }
#' @export
categorize_significance <- function(merged_pathways, colnames_sig, sig_level) {
  merged_pathways <- merged_pathways %>%
    mutate(enriched = case_when(
      !!rlang::sym(colnames_sig[1]) <= sig_level & !!rlang::sym(colnames_sig[2]) <= sig_level & !!rlang::sym(colnames_sig[3]) <= sig_level & !!rlang::sym(colnames_sig[4]) <= sig_level & !!rlang::sym(colnames_sig[5]) <= sig_level & !!rlang::sym(colnames_sig[6]) <= sig_level ~ "all",
      !!rlang::sym(colnames_sig[1]) >= sig_level & !!rlang::sym(colnames_sig[2]) <= sig_level & !!rlang::sym(colnames_sig[3]) <= sig_level & !!rlang::sym(colnames_sig[4]) <= sig_level & !!rlang::sym(colnames_sig[5]) <= sig_level & !!rlang::sym(colnames_sig[6]) <= sig_level ~ "validation dataset only",
      !!rlang::sym(colnames_sig[1]) <= sig_level & !!rlang::sym(colnames_sig[2]) > sig_level & !!rlang::sym(colnames_sig[3]) > sig_level & !!rlang::sym(colnames_sig[4]) > sig_level & !!rlang::sym(colnames_sig[5]) > sig_level & !!rlang::sym(colnames_sig[6]) > sig_level ~ "discovery dataset only",
      !!rlang::sym(colnames_sig[1]) > sig_level & !!rlang::sym(colnames_sig[2]) > sig_level & !!rlang::sym(colnames_sig[3]) > sig_level & !!rlang::sym(colnames_sig[4]) > sig_level & !!rlang::sym(colnames_sig[5]) > sig_level & !!rlang::sym(colnames_sig[6]) > sig_level ~ "none",
      TRUE ~ "none"
        )
      )
  return(merged_pathways)
}



create_color_scheme <- function(nr_of_colors) {
  color_scheme <- RColorBrewer::brewer.pal(8, "Set1")
  color_scheme <- c(color_scheme, "#f2f2f2")
  return(color_scheme)
}


# Make a scatterplot from the merged_pathways resulting from categorize_significance.
# But maybe break up in more different functions for layout of the plot.
plot_scatterplot <- function(merged_pathways_categories, colname_x, colname_y, color_scheme) {
  corpl <- ggplot(merged_pathways_categories, aes(x = !!sym(colname_x), y = !!sym(colname_y), color = enriched)) +#, label = my_label))+#, label = NAME))+
            geom_point(alpha = 0.6, size = 3) +
            xlab("ES discovery dataset")+
            ylab("ES validation dataset")+
            scale_color_manual("Significant in",
                values = c(color_scheme[5], color_scheme[2], color_scheme[9], color_scheme[4]),
              #  breaks = c("1", "2", "3", "4", "5"),
                labels = c("all", "discovery dataset only", "none", "validation dataset only")
                )
  return(corpl)
}

# ES_level is sensitive to minus sign
label_most_significant <- function(merged_pathways_categories, colnames_sig, colnames_ES, sig_level, ES_level) {
  merged_pathways_categories <- merged_pathways_categories %>%
    mutate(to_label = ifelse(
      enriched == "all" & !!sym(colnames_sig[1]) < sig_level & !!sym(colnames_sig[2]) < sig_level & !!sym(colnames_ES[1]) < ES_level & !!sym(colnames_ES[2]) < ES_level
      ,
      NAME, ""
    )) %>%
    ungroup()
    merged_pathways_categories$to_label <- stringr::str_sub(merged_pathways_categories$to_label, start = 10)
    merged_pathways_categories$to_label <- gsub("_", " ", merged_pathways_categories$to_label)
    return(merged_pathways_categories)
}


plot_scatterplot_labeled <- function(merged_pathways_categories, colname_x, colname_y, color_scheme) {
  corpl <- ggplot(merged_pathways_categories, aes(x = !!sym(colname_x), y = !!sym(colname_y), color = enriched, label = to_label)) +
            geom_point(alpha = 0.6, size = 3) +
            geom_text_repel(min.segment.length = 0,
                            colour = "black",
                            #verbose = TRUE,
                            seed = 123,
                            box.padding = 2,
                           # max.time = 3,
                           # max.iter = Inf,
                            size = 5,
                            max.overlaps = Inf,
                            nudge_y = 0.2,
                            nudge_x = 0.2
                            )+
            xlab("ES discovery dataset")+
            ylab("ES validation dataset")+
            scale_color_manual("Significant in",
                values = c(color_scheme[5], color_scheme[2], color_scheme[9], color_scheme[4]),
              #  breaks = c("1", "2", "3", "4", "5"),
                labels = c("all", "discovery dataset only", "none", "validation dataset only")
                )+
            xlim(-1, 1)+
            ylim(-1, 1)
  return(corpl)
}

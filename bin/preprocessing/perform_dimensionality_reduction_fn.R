#' Read Indegree Data
#'
#' This function reads indegree data from a specified file.
#'
#' @param indegree_file A string representing the path to the file containing the indegree data.
#' 
#' @return A data frame containing the indegree data.
#' 
#' @examples
#' \dontrun{
#'   indegree_data <- read_indegree("path/to/indegree_file.csv")
#' }
#' 
#' @export
read_indegree <- function(indegree_file) {
  # Reading in in- or outdegree file and making gene names the rownames. Returns the degree dataframe
  degree <- data.table::fread(indegree_file)
  degree <- as.data.frame(degree)
  rownames(degree) <- degree[, 1]
  degree[, 1] <- NULL
  return(degree)
}

#' Perform Principal Component Analysis (PCA)
#'
#' This function performs Principal Component Analysis (PCA) on the given data frame.
#' It transposes the input data, centers it, and scales it before applying PCA.
#' The function sets a seed for reproducibility.
#'
#' @param data_df A data frame containing the input data for PCA.
#' @return A list containing the results of the PCA.
#' @examples
#' # Example usage:
#' result <- do_pca(data_df)
#' summary(result)
#' @export

do_pca <- function(data_df){
    set.seed(1)
    PCA <- stats::prcomp(t(data_df), center = TRUE, scale. = TRUE)
    return(PCA)
}

match_metadata_order_with_df <- function(metadata_df, indegree_df) {
    metadata_ordered <- metadata_df[match(colnames(indegree_df), metadata_df$sample_name), ]
    return(metadata_ordered)
}

#' Merge PCA Results with Metadata
#'
#' This function merges the results of a Principal Component Analysis (PCA) with a metadata dataframe.
#'
#' @param pca A PCA object, typically the result of a call to `prcomp`, containing the principal components.
#' @param metadata_df A dataframe containing metadata, with a column named `sample_name` that matches the row names of the PCA results.
#'
#' @return A dataframe that combines the PCA results and the metadata, with rows corresponding to the samples.
#'
#' @details The function converts the PCA results and metadata to dataframes, then merges them based on the sample names. The merge is performed such that all rows from the PCA results are retained (`all.x = TRUE`).
#'
#' @examples
#' # Assuming `pca_result` is the output of prcomp and `metadata` is a dataframe with sample metadata
#' merged_data <- merge_pca_met(pca_result, metadata)
#'
#' @export
merge_pca_met <- function(pca, metadata_df) {
    pca <- as.data.frame(pca$x)
    metadata_df <- as.data.frame(metadata_df)
    merged <- merge(pca, metadata_df, by.x = "row.names", by.y = "sample_name", all.x = TRUE)
    # rownames(merged) <- merged$Row.names
    # merged$Row.names <- NULL
    return(merged)
}

#' Plot PCA results with discrete color mapping
#'
#' This function takes PCA results and metadata, merges them, and creates a scatter plot of the first two principal components (PC1 and PC2) with points colored based on a specified metadata column.
#'
#' @param pca A PCA object containing the principal component analysis results.
#' @param metadata_df A data frame containing metadata to be merged with the PCA results.
#' @param to_visualize A character string specifying the column name in the metadata to be used for coloring the points.
#' @param nr_of_colors An integer specifying the number of colors to use in the plot.
#'
#' @return A ggplot object representing the PCA plot with discrete color mapping.
#'
#' @import ggplot2
#' @import rlang
#'
#' @examples
#' # Assuming `pca_results` is a PCA object and `metadata` is a data frame with metadata
#' plot <- plot_12_discrete(pca_results, metadata, "sample_type", 5)
#' print(plot)
plot_12_discrete <- function(pca, metadata_df, to_visualize) { # To do: incorporate, first_pc, second_pc with rlang package ...., then I don't need plot_12 separate from plot_13
    merged <- merge_pca_met(pca, metadata_df)
    summ <- summary(pca)
    summ <- as.data.frame(summ$importance)
    pc1_label <- paste("PC1", sprintf("%0.1f%%", summ[2, 1]*100), sep = " ")
    pc2_label <- paste("PC2", sprintf("%0.1f%%", summ[2, 2]*100), sep = " ")
    which_col <- merged[, c(to_visualize)]
    nr_of_colors <- length(unique(which_col))
    to_visualize <- rlang::sym(to_visualize)
    pc_plot <- ggplot(merged, aes(x = PC1, y = PC2)) +
                    geom_point(aes(color = !!to_visualize), size = 3)+
                    xlab(pc1_label)+
                    ylab(pc2_label)+
                    scale_color_manual(values = standard_colors[1:nr_of_colors])
    return(pc_plot)
}


#' Visualize PCA with Continuous Variable
#'
#' This function creates a scatter plot of the first two principal components (PC1 and PC2) from a PCA object, 
#' with points colored by a continuous variable from the metadata.
#'
#' @param pca A PCA object obtained from a PCA analysis.
#' @param to_visualize A string representing the name of the column in the metadata to visualize as a continuous variable.
#' 
#' @return A ggplot2 object representing the scatter plot of PC1 vs PC2 with points colored by the specified continuous variable.
#'
#' @details The function merges the PCA results with the metadata, extracts the importance of the principal components, 
#' and constructs axis labels that include the percentage of variance explained by PC1 and PC2. 
#' It then creates a scatter plot using ggplot2, with points colored according to the specified continuous variable.
#'
#' @examples
#' # Assuming `pca_result` is a PCA object and `metadata_df` is a data frame with metadata:
#' pca_plot <- pca_visualize_continuous(pca_result, "Age")
#' print(pca_plot)
#'
#' @import ggplot2
#' @import rlang
pca_visualize_continuous <- function(pca, metadata_df, to_visualize){
    # to_visualize must be the name of a column to visualize
    merged <- merge_pca_met(pca, metadata_df)
    print(colnames(merged))
    to_visualize <- rlang::sym(to_visualize)
    merged <- merged %>%
        mutate(cellnr_label = case_when(
            !!to_visualize < 130 ~ Var1,
            TRUE ~ ""
    ))
    summ <- summary(pca)
    summ <- as.data.frame(summ$importance)
    pc1_label <- paste("PC1", sprintf("%0.1f%%", summ[2, 1]*100), sep = " ")
    pc2_label <- paste("PC2", sprintf("%0.1f%%", summ[2, 2]*100), sep = " ")
    to_visualize <- rlang::sym(to_visualize)
    pc_plot <- ggplot(merged, aes(x = PC1, y = PC2, label = cellnr_label)) +
        geom_point(aes(color = !!to_visualize), size = 3)+
        xlab(pc1_label)+
        ylab(pc2_label)+
       # scale_color_gradient(to_visualize)+
        # een scale waarin je lager dan 100 beter kan zien:
        scale_colour_gradient(limits = c(0, 500), oob = scales::squish)+
        geom_text_repel()
    return(pc_plot)
}
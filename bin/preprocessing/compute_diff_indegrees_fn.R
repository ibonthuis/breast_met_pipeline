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

#' Read and Match Metadata
#'
#' This function reads metadata from a specified file and matches it with an indegree dataframe.
#'
#' @param metadata_file A string specifying the path to the metadata file.
#' @param indegree_df A dataframe containing indegree information.
#' 
#' @return A dataframe with matched metadata and indegree information.
#' 
#' @examples
#' metadata <- read_and_match_metadata("path/to/metadata_file.csv", indegree_df)
#' 
#' @export
match_metadata <- function(metadata_df, indegree_df) {
    if (!all(colnames(indegree_df) %in% metadata_df$sample_name)) {
        stop("column names of indegree_df don't match sample_name")
    }
    metadata_ordered <- metadata_df[match(colnames(indegree_df), metadata_df$sample_name), ]
    return(metadata_ordered)
}


#' Create Top Table for Paired Indegree Data
#'
#' This function generates a top table for paired indegree data based on the provided metadata.
#'
#' @param indegree_df A data frame containing the indegree data.
#' @param meta A data frame containing the metadata.
#' @param patient_col A string specifying the column name in the metadata that contains patient identifiers.
#' @param type_col A string specifying the column name in the metadata that contains the type information.
#'
#' @return A data frame representing the top table for paired indegree data.
#' @export
#'
#' @examples
#' \dontrun{
#' indegree_df <- data.frame(patient_id = c("P1", "P2"), indegree = c(5, 10))
#' meta <- data.frame(patient_id = c("P1", "P2"), type = c("A", "B"))
#' create_toptable_paired(indegree_df, meta, "patient_id", "type")
#' }
create_toptable_paired <- function(indegree_df, meta, patient_col, type_col){
  # indegree_df is a df containing your samples in columns and indegrees(target genes) as rows
  # meta is a df containing the metadata of your indegree df, with rows in the same order as your column containing sample names (here sample_name)
  # patient_col should be the column name of interest, for example "patient"
  # type_col should be the column name of interest, for example "sample_type" for the design matrix of the limma analysis
  #metadata_ordered <- meta[match(colnames(indegree_df), meta$sample_name), ]
  
  metadata_ordered <- as.data.frame(meta)
  
  metadata_ordered <- metadata_ordered %>%
    filter(sample_name %in% colnames(indegree_df))
  #nrow(metadata_ordered)
 
  patient_col <- rlang::ensym(patient_col)
  type_col <- rlang::ensym(type_col)
  pat <- factor(metadata_ordered[[as.character(patient_col)]])
  type <- factor(metadata_ordered[[as.character(type_col)]])
  design <- model.matrix(~pat+type)
  fit <- lmFit(indegree_df, design); fit <- eBayes(fit); topTable(fit, coef="typePrimary")
  toptable <- topTable(fit, coef="typePrimary", number=nrow(indegree_df))

  return(toptable)
}


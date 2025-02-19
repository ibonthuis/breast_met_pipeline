#' Create a List of Gene Ranked Files
#'
#' This function sets the working directory to the specified directory name and 
#' retrieves a list of files with the ".rnk" extension.
#'
#' @param dir_name A character string specifying the directory path where the 
#' gene ranked files are located. When running the script using snakemake, the 
#' directory path should be upstream of where this script is living.
#' 
#' @return A character vector containing the names of the files with the ".rnk" 
#' extension in the specified directory.
#'
#' @examples
#' \dontrun{
#'   gene_ranked_files <- create_list_of_gene_ranked_files("/path/to/directory")
#'   print(gene_ranked_files)
#' }
create_list_of_gene_ranked_files <- function(dir_name){
    setwd(dir_name)
    analyses <- list.files(pattern = ".rnk")
    return(analyses)
}

run_gsea <- function(list_of_analyses){

}
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


# run_gsea <- function(ranked_df){
#     ranked_df_name <- names(ranked_df)
#     assign(paste("gsea", ranked_df_name, sep=""), 
#         paste(
#         "java -Xmx4096m -cp gsea2-2.0.13.jar xtools.gsea.GseaPreranked -gmx ",
#         signaturename,
#         " -collapse false -mode Max_probe -norm meandiv -nperm ", 
#         nriter,
#         " -rnk ./", 
#         ranked_df_name
#         ,
#         " -scoring_scheme weighted -rpt_label ", 
#         substr(ranked_df_name, 1, nchar(ranked_df_name)-4),
#         " -include_only_symbols true -make_sets true -plot_top_x 0 -rnd_seed timestamp -set_max 250 -set_min 1 -zip_report false -out ./ -gui false", sep=""
#         ) )
#       set.seed(1)
#       system(eval(parse(text=paste("gsea", ranked_df_name, sep=""))))

#     # read in results and save results and significant results
#       gseadir <- list.files(pattern=paste(substr(ranked_df, 1, nchar(ranked_df)-4),".GseaPreranked", sep=""))
#       setwd(gseadir) 
#       #getwd()
#       gseafiles <- list.files(pattern="gsea_report") # gsea outputs .html and .xls files for signatures with positive and negative enrichment scores
#       neg <- read.delim(gseafiles[2]) # signatures with negative enrichment scores
#       pos <- read.delim(gseafiles[4]) # signatures with positive enrichment scores
#       head(pos, n = 50)
#       head(pos[order(pos$FDR.q.val),])
#       head(neg[order(neg$FDR.q.val),])
# }

#' @name perform_gsea
#'
#' @description This function performs Gene Set Enrichment Analysis (GSEA) on
#' results of differential analysis based on JDR factors. See
#' \code{\link{differential_analysis}} for details.
#'
#' @inheritParams differential_analysis
#' @param diff_results A data frame containing either a limma toptable or the
#' output of \code{link{differential_analysis}}.
#' @param gene_set A .gmt file containing gene sets of interest or a list of
#' gene sets to test.
#' @param differential Logiacal. Whether differential analysis has been
#' performed on the data or this is just being run on weights. I will probably
#' eventually remove the differential options entirely, but for now.
#' @param limma Logical. Whether limma was used for differential. Only needed
#' if differential = TRUE. If false, it assumes wilcoxon was used.
#' @param ... Any other parameters accepted by fgsea.
#' @seealso \code{\link{fgsea::fgsea}}
#'
#' @returns Results of GSEA
#' @export

perform_gsea <- function(diff_results, differential = FALSE, limma = TRUE,
                         gene_set, save_file = TRUE,
                         file_name = NULL, ...) {
  # get gene set
  if (is.character(gene_set)) {
    gene_set <- fgsea::gmtPathways(gene_set)
  }

  # sanity checks
  # check that omic and gene_set use the same annotation
  gs_names <- unlist(unique(gene_set))
  omic_names <- rownames(diff_results)
  .check_names(omic_names, gs_names, partial = TRUE,
               err_msg = "the gene sets and the omic data use the same annotation") # nolint

  # create rank
  rnk <- .create_rank(diff_results, differential, limma)

  gsea_res <- fgsea::fgsea(pathways = gene_set, rnk, ...)


  if (save_file) {
    if (is.null(file_name)) {
      file_name <- "gsea_analysis_results.RData"
    }
    save(rnk, gsea_res, file = file_name)
  }
  return(gsea_res)
}
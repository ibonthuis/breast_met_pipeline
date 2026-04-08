#' Cluster GSEA Enrichment
#'
#' This function clusters pathways based on the Jaccard index of their leading gene sets and optionally generates a heatmap of the clustered pathways.
#'
#' @param gsea_sign A data frame containing GSEA results, including pathway names, leading genes, and NES scores.
#' @param lead_genes_colname A string specifying the column name in `gsea_sign` that contains the leading genes for each pathway.
#' @param path_colname A string specifying the column name in `gsea_sign` that contains the pathway names.
#' @param jaccard_hclust_cuts A numeric vector specifying the cut heights for hierarchical clustering of pathways based on the Jaccard index. Default is `c(0.5, 1, 1.2, 1.5)`.
#' @param lead_genes_split A string specifying the delimiter used to split the leading genes in the `lead_genes_colname`. Default is `';'`.
#' @param nes_colname A string specifying the column name in `gsea_sign` that contains the normalized enrichment scores (NES). Default is `'NES'`.
#' @param hmap_outpath A string specifying the file path to save the clustered pathway heatmap. If `NULL`, no heatmap is generated. Default is `NULL`.
#' @param hmap_title A string specifying the title of the heatmap. Default is `NULL`.
#'
#' @return A modified version of the `gsea_sign` data frame with additional columns indicating pathway clusters for each cut height in `jaccard_hclust_cuts`.
#'
#' @examples
#' # Example usage:
#' clustered_gsea <- cluster_gsea_enrichment(
#'   gsea_sign = gsea_results,
#'   lead_genes_colname = "leadingEdge",
#'   path_colname = "pathway",
#'   jaccard_hclust_cuts = c(0.5, 1),
#'   lead_genes_split = ";",
#'   nes_colname = "NES",
#'   hmap_outpath = "clustered_heatmap.pdf",
#'   hmap_title = "Clustered Pathway Heatmap"
#' )


# clustering gsea enrichment results by jaccard idx (+ heatmap)
cluster_gsea_enrichment <- function(gsea_sign, lead_genes_colname, path_colname, jaccard_hclust_cuts = c(0.5, 1, 1.2, 1.5), lead_genes_split = ';',
                                    nes_colname = 'NES', hmap_outpath = NULL, hmap_title = NULL){
                                      
  # cluster pathways based on jaccard idx -----------------------------------
    paths_genes_list <- lapply(gsea_sign[[lead_genes_colname]], function(x){
    genelist <- unlist(strsplit(x, split=lead_genes_split, fixed=T)) # TODO wtf it looks like any split works
  })

 names(paths_genes_list) <- gsea_sign[[path_colname]]
  
  # calculate jaccard score between each pathway leading gene set
  path_jaccard <- lapply(paths_genes_list, function(x){
    p1 <- lapply(paths_genes_list, function(y){
      jacc_idx <- as.numeric(round(length(intersect(x, y)) / length(union(x,y)), digits = 4))
    })
    return(unlist(p1))
  })
  
 path_jaccard_mtx <- do.call('cbind', path_jaccard)
  # heatmap(path_jaccard_mtx)
  # clustering with hclust
  path_hclust <- hclust(dist(path_jaccard_mtx), method = "average")
  #plot(path_hclust, hang = -1, cex = 0.4)
  # making clustered pathway heatmap
  if(!is.null(hmap_outpath)){
    make_clustered_gsea_hmap(gsea_sign, path_jaccard_mtx, path_colname, nes_colname, hmap_outpath, hmap_title)
  }

  for(cutnr in jaccard_hclust_cuts){
    # cut the hclust tree at given point
    path_hclust_cut <- cutree(path_hclust, h = cutnr)
    # merge with gsea result
    if(identical(gsea_sign[[path_colname]], names(path_hclust_cut))){
      gsea_sign[[paste0('path_cluster_cut_', gsub('\\.', '', as.character(cutnr)))]] <- path_hclust_cut
    }
  }
  return(gsea_sign)
 # return(path_jaccard_mtx)
}

make_clustered_gsea_hmap <- function(gsea_sign, path_jaccard_mtx, path_colname, nes_colname, hmap_outpath, hmap_title){
  # make clustered heatmap
  nes_anno <- sapply(colnames(path_jaccard_mtx), function(path){
    nes <- ifelse(sign(gsea_sign[[nes_colname]][gsea_sign[[path_colname]] == path]) == 1, 'pos', 'neg')
  })
  ha = HeatmapAnnotation(
    NES = anno_simple(nes_anno, col = c("pos" = "green", "neg" = "blue")),
    annotation_name_side = "left")
  options(bitmapType='cairo')
  png(filename=hmap_outpath, width=12, height=12,units="in",res=1000)
  # pdf(hmap_outpath, w=8, h=6)
  condition_heat <- Heatmap(as.matrix(path_jaccard_mtx), row_names_max_width = max_text_width(
                            rownames(as.matrix(path_jaccard_mtx)), 
                            gp = gpar(fontsize = 14)), 
                            border="white",
                            rect_gp = gpar(col = "white", lwd = 2), column_title = hmap_title,
                            cluster_columns = T, cluster_rows= T, col = brewer.pal(5, "YlOrRd"),
                            show_heatmap_legend = F, top_annotation = ha,
                            row_names_gp = gpar(fontsize = 6),
                            column_names_gp = gpar(fontsize = 6))
  
    show(condition_heat)
  dev.off()
}

################################################
# choose best pathway from jaccard clustering
#' Best pathway clustering
#'
#' This function finds the best pathway clusters based on the Jaccard index of their leading gene sets.
#'
#' @param gsea_clust A data frame containing the Jaccard index clustering results from the cluster_gsea_enrichment function, including pathway names, leading genes, and NES scores, and some columns named path_cluster_cut_..
#' @param clust_colname A column name that is present in `gsea_clust` that contains the clusters that a pathway would end up in, e.g. path_cluster_cut_12
#'
#' @return A dataframe containing the pathways that summarize best the clustered pathways.

best_pathway_clust <- function(gsea_clust, clust_colname, nes_colname = 'NES'){
  gsea_main_all <- lapply(c('pos', 'neg'), function(sign){
    if(sign == 'pos'){
      gsea_main <- gsea_clust[sign(gsea_clust[[nes_colname]]) == 1, ]
      gsea_main <- arrange(gsea_main, desc(get(nes_colname)))
    } else{
      gsea_main <- gsea_clust[sign(gsea_clust[[nes_colname]]) == -1, ]
      gsea_main <- arrange(gsea_main, get(nes_colname))
    }
    # keep 1st pathway from cluster
    gsea_main <- gsea_main[!(duplicated(gsea_main[[clust_colname]])), ]
    return(gsea_main)
  })
  gsea_main_all <- do.call(rbind, gsea_main_all)
  if(nrow(gsea_main_all) > 0){
    return(gsea_main_all)
  }
}

# gsea_res <- fread(file.path(gsea_dir_path, 'gseapy.gene_set.prerank.report.csv'))
#   # filter to significant results
#   gsea_sign <- gsea_res[gsea_res$`FDR q-val` <= fdr_thr, ]
#   gsea_sign <- gsea_sign[gsea_sign$`FWER p-val` <= fwer_thr]
#   if(nrow(gsea_sign) > 0){
#     # cluster pathways by jaccard idx and make heatmap
#     hmap_outpath <- file.path(gsea_dir_path, paste0('hmap_fdr', as.character(fdr_thr), '_fwer',as.character(fwer_thr), '.png'))
#     gsea_clust <- cluster_gsea_enrichment(gsea_sign, 'Lead_genes', 'Term', hmap_outpath = hmap_outpath, hmap_title = 'stroma vs tumor')
#   }
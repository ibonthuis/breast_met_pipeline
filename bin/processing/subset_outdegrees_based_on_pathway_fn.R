reading_gene_set_file <- function(gene_set_file) {
    geneset_list <- fgsea::gmtPathways(gene_set_file)
    return(geneset_list)
}

# pathway_of_interest being a character string matching exactly one of the pathway names
# in the geneset_list
select_genes_from_pathway_of_interest <- function(geneset_list, pathway_of_interest) {
    genes_of_pathway_of_interest <- geneset_list[pathway_of_interest]
    return(genes_of_pathway_of_interest)
}

filter_indegrees_by_pathway_of_interest <- function(indegree_df, genes_of_pathway_of_interest) {
    filtered_indegrees <- indegree_df[rownames(indegree_df) %in% unlist(genes_of_pathway_of_interest), ]
    return(filtered_indegrees)
}

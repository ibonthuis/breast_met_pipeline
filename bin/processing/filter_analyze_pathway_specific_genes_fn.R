reading_gene_set_file <- function(gene_set_file) {
    geneset_list <- fgsea::gmtPathways(gene_set_file)
    return(geneset_list)
}

# pathway_of_interest being a character string matching exactly one of the pathway names
# in the geneset_list
select_genes_from_pathway_of_interest <- function(geneset_list, pathway_of_interest) {
    genes_of_pathway_of_interest <- geneset_list[pathway_of_interest]
    genes_of_pathway_of_interest <- unlist(genes_of_pathway_of_interest)
    return(genes_of_pathway_of_interest)
}




# filter any df, such as edges, indegree and expression matrix with the genes in the rownames 
filter_df_by_pathway_of_interest <- function(df, genes_of_pathway_of_interest) {
    filtered_indegrees <- df[rownames(df) %in% genes_of_pathway_of_interest, ]
    return(filtered_indegrees)
}

filter_indegree_file <- function(gene_set_file, pathway_of_interest, indegree_file) {
    indegree_df <- read_indegree(indegree_file)
    gene_set <- reading_gene_set_file(gene_set_file)
    genes <- select_genes_from_pathway_of_interest(gene_set, pathway_of_interest)
    df_filtered <- filter_df_by_pathway_of_interest(indegree_df, genes)
    return(df_filtered)
}

# filter an edge df Target genes by genes_of_pathway_of_interest and an integer for logfc_filter

read_edge_file <- function(edge_file) {
    edge_df <- fread(edge_file)
    return(edge_df)
}

filter_edges_by_pathway_of_interest <- function(df_edges, pathway_of_interest) {
    gene_set <- reading_gene_set_file(gene_set_file)
    genes <- select_genes_from_pathway_of_interest(gene_set, pathway_of_interest)
    second_column <- colnames(df_edges)[2]
    filtered_df_edges <- df_edges %>%
        filter(!!rlang::sym(second_column) %in% genes)
    return(filtered_df_edges)
}


prep_edge_df_for_differential_analysis <- function(filtered_df_edges) {
    filtered_df_edges <- as.data.frame(filtered_df_edges)
    second_column <- colnames(filtered_df_edges)[2]
    first_column <- colnames(filtered_df_edges)[1]

    rownames(filtered_df_edges) <- paste(filtered_df_edges[[first_column]], filtered_df_edges[[second_column]], sep = "_")
    filtered_df_edges[[first_column]] <- NULL
    filtered_df_edges[[second_column]] <- NULL
    return(filtered_df_edges)
}

# A function for after toptables with limma which will always result in logFC
filter_diff_by_lfc <- function(differential_df, lfc_filter) {
    differential_df <- differential_df %>%
        filter(abs(logFC) >= lfc_filter)
}

split_edge_df_rownames <- function(diff_edges) {
    #toptable_edge_filt$edge <- rownames(toptable_edge_filt)
    diff_edges <- diff_edges %>%
        rownames_to_column(var = "edge")

    limma_edges <- diff_edges %>%
        separate(edge, into = c("TF", "Target"), sep = "_")

    return(limma_edges)
}



filter_edges <- function(edge_file, metadata_file, pathway_of_interest, lfc_filter){
    df_edges <- read_edge_file(edge_file)
    edges_filtered <- filter_edges_by_pathway_of_interest(df_edges, pathway_of_interest)
    prepped_edges <- prep_edge_df_for_differential_analysis(edges_filtered)
    metadata <- fread(metadata_file)
    diff_edges <- create_toptable_paired(edges_filtered, metadata, "patient", "sample_type")
    filtered_diff_edges <- filter_diff_by_lfc(diff_edges, lfc_filter)
    filtered_diff_edges_for_cytoscape <- split_edge_df_rownames(filtered_diff_edges)
    return(filtered_diff_edges_for_cytoscape)
}

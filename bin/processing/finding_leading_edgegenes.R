corr <- fread("/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results_expression/limma_covariate_correction/cosgrove/differential_genesets.tsv")

norm <- fread("/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results_expression/differential_indegrees/cosgrove/differential_genesets.tsv")

head(corr)
corr[corr$V1 == "REACTOME_IMMUNE_SYSTEM",]

corr_le <- corr[corr$V1 == "REACTOME_IMMUNE_SYSTEM", ]
corr_le_genes <- corr_le$V8
corr_le_genes <- strsplit(corr_le_genes, "|", fixed = TRUE)
length(corr_le_genes[[1]])

norm_le <- norm[norm$V1 == "REACTOME_IMMUNE_SYSTEM",]
norm_le_genes <- norm_le$V8
norm_le_genes <- strsplit(norm_le_genes, "|", fixed = TRUE)
length(norm_le_genes[[1]])

dif <- setdiff(norm_le_genes, corr_le_genes) #190 (norm_le_genes has 243 genes)
str(dif)
str(norm_le_genes)
setdiff(corr_le_genes, norm_le_genes) #229 (corr_le_genes has 282)

dif <- strsplit(dif, "|", fixed = TRUE)

for (i in 1:10) {
    norm_le_genes <- unlist(norm_le_genes)
   print(as.name((norm_le_genes[i])))
}

for (i in 1:10) {
    corr_le_genes <- unlist(corr_le_genes)
   print(as.name((corr_le_genes[i])))
}





diff_ind_leading <- fread("/storage/kuijjerarea/ine/projects/BRCA_MET/breast_met_pipeline/snakemake_results/differential_indegrees/cosgrove/differential_genesets.tsv")

diff_ind_le <- diff_ind_leading[diff_ind_leading$V1 == "REACTOME_INTERFERON_SIGNALING",]
diff_ind_le <- diff_ind_le$V8
diff_ind_le <- strsplit(diff_ind_le, "|", fixed = TRUE)
for (i in 1:10) {
    diff_ind_le <- unlist(diff_ind_le)
   print(as.name((diff_ind_le[i])))
}

print(unlist(diff_ind_le ))

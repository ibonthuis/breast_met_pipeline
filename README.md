# breast_met_pipeline
Pipeline to analyse breast cancer metastasis using networks


## Folder structure 
The folder structure right now is roughly like this. The bin folder contains all the scripts, the data folder all the input data. This one is a bit polluted with different inputs, but they can all be used when setting different options in the config file.
The output of snakemake will always be like the folder "snakemake_results_codereview_10102025".

```
breast_met_pipeline
├── bin
│    ├── preprocessing
│    │   ├── compute_covariate_diff_indegrees.R
│    │   ├── compute_diff_indegrees_fn.R
│    │   ├── compute_diff_indegrees.R
│    │   ├── compute_diff_outdegrees.R
│    │   ├── compute_gse_fn.R
│    │   ├── compute_gse.R
│    │   ├── perform_dimensionality_reduction_fn.R
│    │   └── perform_dimensionality_reduction.R
│    ├── processing
│    │   ├── filter_analyze_pathway_specific_genes_fn.R
│    │   ├── filter_analyze_pathway_specific_genes.R
│    │   ├── filter_indegrees_pathway_specific_genes.R
│    │   ├── overlap_indegrees.R
│    │   ├── overlap_outdegrees_fn.R
│    │   ├── overlap_outdegrees.R
│    │   ├── overlap_pathways_fn.R
│    │   ├── overlap_pathways.R
│    │   ├── subset_outdegrees_based_on_pathway_fn.R
│    │   └── subset_outdegrees_based_on_pathway.R
│    └── visualisation_in_cytoscape
│    data
│    ├── gene_sets
│    │   ├── c2.cp.reactome.v5.0.symbols.gmt
│    │   └── gene_set.gmt
│    ├── github_data
│    │   ├── aurora2
│    │   │   ├── filtered_indegree.csv
│    │   │   └── metadata.csv
│    │   └── cosgrove2
│    │       ├── filtered_indegree.csv
│    │       └── metadata.csv
│    ├── input_data
│    │   ├── aurora
│    │   │   ├── edges_filtered_for_interferon_signallingpathway.tsv
│    │   │   ├── filtered_expression_nobatch.tsv
│    │   │   ├── filtered_expression_only_metastases.tsv
│    │   │   ├── filtered_expression_raw_counts_metastases.tsv
│    │   │   ├── filtered_expression_raw_counts.tsv
│    │   │   ├── filtered_expression.tsv
│    │   │   ├── filtered_indegree.csv
│    │   │   ├── filtered_outdegrees.csv
│    │   │   └── metadata.csv
│    │   └── cosgrove
│    │       ├── differential_edges_filtered_for_interferon_signallingpathway_filtered_by_lfc1.5 copy.tsv
│    │       ├── edges_filtered_for_interferon_signallingpathway.tsv
│    │       ├── filtered_expression.tsv
│    │       ├── filtered_indegree.csv
│    │       ├── filtered_outdegrees.csv
│    │       ├── interferon_pathway_specific_edges.tsv
│    │       └── metadata.csv
│    snakemake_results_codereview_10102025
│    ├── differential_indegrees
│    │   ├── aurora
│    │   │   ├── differential_genesets.RData
│    │   │   ├── differential_genesets.tsv
│    │   │   ├── differential_indegrees_rank_file.RData
│    │   │   ├── differential_indegrees.RData
│    │   │   └── enrichment_bubble_plot_.pdf
│    │   └── cosgrove
│    │       ├── differential_genesets.RData
│    │       ├── differential_genesets.tsv
│    │       ├── differential_indegrees_rank_file.RData
│    │       ├── differential_indegrees.RData
│    │       └── enrichment_bubble_plot_.pdf
│    ├── dimensionality_reduction
│    │   ├── aurora
│    │   │   └── PCA_sample_type_pc12.pdf
│    │   └── cosgrove
│    │       └── PCA_sample_type_pc12.pdf
│    └── overlapping_pathways
│        └── overlapping_pathways_all.tsv
├── config.yaml
├── dag.pdf
├── .gitignore
├── LICENSE
├── README.md
└── Snakefile
```

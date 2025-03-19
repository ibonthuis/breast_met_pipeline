## How to run this?
## snakemake --cores 1 -np ### For dry run

## Libraries
import os 
import sys
import glob
from pathlib import Path
import time

## Config

global CONFIG_PATH
CONFIG_PATH = "config.yaml"
configfile: CONFIG_PATH

## Container
# container: config["container"]
# print("\n; Available container:\n %s" % config["container"])

## Directories
DATA_DIR = config["data_dir"]
DATASET_NAMES = config["types"]
BASE_OUTPUT_DIR = config["output_dir"]
DIFFERENTIAL_INDEGREE_OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, "differential_indegrees")
DIMENSIONALITY_REDUCTION_OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, "dimensionality_reduction")
OVERLAPPING_PATHWAYS_DIR = os.path.join(BASE_OUTPUT_DIR, "overlapping_pathways")
PATHWAY_SPECIFIC_SUBSETS_DIR = os.path.join(BASE_OUTPUT_DIR, "pathway_specific_subsets")

## Input files ##
INPUT_METADATA = os.path.join(DATA_DIR, "{dataset_type}", "metadata.csv")
INPUT_INDEGREES = os.path.join(DATA_DIR, "{dataset_type}", "filtered_indegree.csv")
GENE_SET_FILE = config["gene_set_file"]

## Other inputs ##
VIS_VAR = config["visualisation_var"]
P_THRESH = config["p_threshold"]
PATHWAY = config["pathway_of_interest"]

## Output files ##
# DIFFERENTIAL_INDEGREES_TAB = os.path.join(DIFFERENTIAL_INDEGREE_OUTPUT_DIR, "{dataset_type}", "differential_indegrees.tsv")
DIFFERENTIAL_INDEGREES_RDATA = os.path.join(DIFFERENTIAL_INDEGREE_OUTPUT_DIR, "{dataset_type}", "differential_indegrees.RData")
DIFFERENTIAL_INDEGREES_RANKED_RDATA = os.path.join(DIFFERENTIAL_INDEGREE_OUTPUT_DIR, "{dataset_type}", "differential_indegrees_rank_file.RData")
ENRICHMENT_RESULTS_RDATA = os.path.join(DIFFERENTIAL_INDEGREE_OUTPUT_DIR, "{dataset_type}", "differential_genesets.RData")
ENRICHMENT_RESULTS_TSV = os.path.join(DIFFERENTIAL_INDEGREE_OUTPUT_DIR, "{dataset_type}", "differential_genesets.tsv")
OVERLAPPING_RESULTS_TSV = os.path.join(OVERLAPPING_PATHWAYS_DIR, "overlapping_pathways_all.tsv")
PATHWAY_SPECIFIC_INDEGREES = os.path.join(PATHWAY_SPECIFIC_SUBSETS_DIR, "{dataset_type}", "diff_indegrees_pathway_specific_genes.RData")

## Output plots ##
PCA_PLOT_PDF = os.path.join(DIMENSIONALITY_REDUCTION_OUTPUT_DIR, "{dataset_type}", "PCA_{visualisation_var}_pc12.pdf") # In the R script it's written as follows:  pdf(pcaplot, file.path(OUTPUT_DIR, paste0("PCA", VARIABLE, "pc12.pdf")))
ENRICHMENT_RESULTS_PDF = os.path.join(DIFFERENTIAL_INDEGREE_OUTPUT_DIR, "{dataset_type}", "enrichment_bubble_plot_.pdf")

## Rule ALL ##
rule all:
    input:
        expand(DIFFERENTIAL_INDEGREES_RDATA, dataset_type = DATASET_NAMES), \
        expand(PCA_PLOT_PDF, dataset_type = DATASET_NAMES, visualisation_var = VIS_VAR), \
        expand(ENRICHMENT_RESULTS_RDATA, dataset_type = DATASET_NAMES),\
        expand(ENRICHMENT_RESULTS_TSV, dataset_type = DATASET_NAMES), \
        expand(ENRICHMENT_RESULTS_PDF, dataset_type = DATASET_NAMES), \
        OVERLAPPING_RESULTS_TSV

        

## Rules ##

rule compute_differential_indegrees:
    """
    This rule takes original input and computes differential indegrees.

    Inputs
    ------
    INPUT_METADATA:
        Path to the metadata file containing sample information.
    INPUT_INDEGREES:
        Path to the file containing indegree data for analysis.
    Outputs
    -------
    DIFFERENTIAL_INDEGREES_TAB:
        blablabla in tab
    DIFFERENTIAL_INDEGREES_RDATA:
        blablabla in RData
    """
    input:
        metadata = INPUT_METADATA, \
        indegrees = INPUT_INDEGREES
    output:
        # DIFFERENTIAL_INDEGREES_TAB, \
        DIFFERENTIAL_INDEGREES_RDATA, \
        DIFFERENTIAL_INDEGREES_RANKED_RDATA
    message:
        "; Running differential indegree computation on {input}."
    params:
        bin = os.path.join(config["bin"], "preprocessing"), \
        output_dir = os.path.join(BASE_OUTPUT_DIR, "differential_indegrees", "{dataset_type}")
    shell:
        """
        Rscript {params.bin}/compute_diff_indegrees.R \
            -i {input.indegrees} \
            -m {input.metadata} \
            -o {params.output_dir}
        """

rule perform_dimensionality_reduction:
    """
    This rule takes original input and reduces dimensionality of the indegree data and plots this.

    Inputs
    ------
    INPUT_METADATA:
        Path to the metadata file containing sample information.
    INPUT_INDEGREES:
        Path to the file containing indegree data for analysis.
    VIS_VAR:
        Variable to visualize.
    Outputs
    -------
    PCA_PLOT_PDF:
        blablabla in pdf

    """
    input:
        metadata = INPUT_METADATA, \
        indegrees = INPUT_INDEGREES
    output:
        PCA_PLOT_PDF
    message:
        "; Running dimensionality reduction of indegrees on {input}."
    container:
        config["container"]
    params:
        bin = os.path.join(config["bin"], "preprocessing"), \
        output_dir = os.path.join(BASE_OUTPUT_DIR, "dimensionality_reduction", "{dataset_type}"), \
        visualization_var = VIS_VAR, \
      #  dataset_type = "{dataset_type}"
    shell:
        """
        Rscript {params.bin}/perform_dimensionality_reduction.R \
            -i {input.indegrees} \
            -m {input.metadata} \
            -v {params.visualization_var} \
            -o {params.output_dir}
        """

rule run_gsea_on_ranks:
    """
    This rule takes ranked genes as input and computes enriched gene sets.

    Inputs
    ------
    DIFFERENTIAL_INDEGREES_RANKED_RDATA:
        blalbla
    GENE_SET_FILE:
        blablabla
    Outputs
    -------
    ENRICHMENT_RESULTS_RDATA:
        Table or list of tables containing the gene sets passing the p-value threshold 
    ENRICHMENT_RESULTS_PDF:
        A bubble plot (or bubble plots) in one pdf file.

    """
    input:
        ranks = DIFFERENTIAL_INDEGREES_RANKED_RDATA
    output:
        # rdata = expand(ENRICHMENT_RESULTS_RDATA, dataset_type="{dataset_type}"), \
        # pdf = expand(ENRICHMENT_RESULTS_PDF, dataset_type="{dataset_type}")
        ENRICHMENT_RESULTS_RDATA, \
        ENRICHMENT_RESULTS_TSV, \
        ENRICHMENT_RESULTS_PDF
    message:
        "; Running GSEA."
    params:
        bin = os.path.join(config["bin"], "preprocessing"), \
        output_dir = os.path.join(BASE_OUTPUT_DIR, "differential_indegrees", "{dataset_type}"), \
        p_threshold = P_THRESH,
        genes = GENE_SET_FILE
    shell:
        """
        echo "; I love snakemake" ;
        Rscript {params.bin}/compute_gse.R \
            -i {input.ranks} \
            -g {params.genes} \
            -p {params.p_threshold} \
            -o {params.output_dir}
        """
    
rule find_overlapping_pathways:
    """
    This rule takes the significant pathways output from both datasets as computed by run_gsea_on_ranks rule. 
    
    Inputs
    ------
    pathways:
        The RData files from run_gsea_on_ranks
    Outputs
    -------
    ENRICHMENT_RESULTS_RDATA:
        Table or list of tables containing the gene sets passing the p-value threshold 
    ENRICHMENT_RESULTS_PDF:
        A bubble plot (or bubble plots) in one pdf file.
    
    Params
    ------
    GENE_SET_FILE:
        gene set file is not dependent on wildcards, input and output should be something produced by snakemake


    """
    input:
    output:
        OVERLAPPING_RESULTS_TSV #, \
        # OVERLAPPING_PATHWAYS_HEATMAP_PDF, \
        # SCATTERPLOT_PATHWAYS_PDF
    message:
        "; Overlapping the pathways"
    params:
        bin = os.path.join(config["bin"], "processing"), \
        output_dir = os.path.join(BASE_OUTPUT_DIR, "overlapping_pathways"), \
        pathways = expand(ENRICHMENT_RESULTS_RDATA, dataset_type=DATASET_NAMES)
    shell:
        """
        echo "; I love snakemake more" ;
        echo {params.pathways}
        Rscript {params.bin}/overlap_pathways.R \
            -i {params.pathways} \
            -o {params.output_dir}
         """

rule filter_indegrees_on_pathway:
    """
    This rule takes indegrees and one specific pathway as input and computes differential indegrees for that specific pathway.

    Inputs
    ------
    INPUT_METADATA:
        Path to the metadata file containing sample information.
    INPUT_INDEGREES:
        Path to the file containing indegree data for analysis.
    Outputs
    -------
    PATHWAY_SPECIFIC_INDEGREES:
        Table or list of tables containing the indegrees with differential values between primary and metastasis. 
    """
    input:
        metadata = INPUT_METADATA, \
        indegrees = INPUT_INDEGREES
    output:
        PATHWAY_SPECIFIC_INDEGREES
    message:
        "; Running pathway specific indegree filtering."
    params:
        bin = os.path.join(config["bin"], "processing"), \
        output_dir = os.path.join(BASE_OUTPUT_DIR, "pathway_specific_subsets", "{dataset_type}"), \
        pathway = PATWHAY,
        genes = GENE_SET_FILE
    shell:
        """
        echo "; I love indegrees in snakemake" ;
        Rscript {params.bin}/filter_indegrees_pathway_specific_genes.R \
            -i {input.indegrees} \
            -m {params.metadata} \
            -p {params.pathway} \
            -G {params.genes} \
            -o {params.output_dir}
        """
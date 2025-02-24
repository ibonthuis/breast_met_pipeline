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

## Input files ##
INPUT_METADATA = os.path.join(DATA_DIR, "{dataset_type}", "metadata.csv")
INPUT_INDEGREES = os.path.join(DATA_DIR, "{dataset_type}", "filtered_indegree.csv")
VIS_VAR = config["visualisation_var"]
GENE_SET_FILE = config["gene_set_file"]

## Output files ##
# DIFFERENTIAL_INDEGREES_TAB = os.path.join(DIFFERENTIAL_INDEGREE_OUTPUT_DIR, "{dataset_type}", "differential_indegrees.tsv")
DIFFERENTIAL_INDEGREES_RDATA = os.path.join(DIFFERENTIAL_INDEGREE_OUTPUT_DIR, "{dataset_type}", "differential_indegrees.RData")
DIFFERENTIAL_INDEGREES_RANKED_RDATA = os.path.join(DIFFERENTIAL_INDEGREE_OUTPUT_DIR, "{dataset_type}", "differential_indegrees_rank_file.RData")
ENRICHMENT_RESULTS = os.path.join(DIFFERENTIAL_INDEGREE_OUTPUT_DIR, "{dataset_type}", "enrichment.RData")

## Output plots ##
PCA_PLOT_PDF = os.path.join(DIMENSIONALITY_REDUCTION_OUTPUT_DIR, "{dataset_type}_pc12.pdf") # In the R script it's written as follows:  pdf(pcaplot, file.path(OUTPUT_DIR, paste0("PCA", VARIABLE, "pc12.pdf")))


## Rule ALL ##
rule all:
    input:
        expand(PCA_PLOT_PDF, dataset_type = DATASET_NAMES), \
        expand(ENRICHMENT_RESULTS, dataset_type = DATASET_NAMES)
        

## Rules ##

rule compute_differential_indegrees:
    """
    This rule takes original input and computes differential indegrees.

    Inputs
    ------
    INPUT_METADATA:
        blalbla
    INPUT_INDEGREES:
        blablabla
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
        blalbla
    INPUT_INDEGREES:
        blablabla
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
        output_dir = os.path.join(BASE_OUTPUT_DIR, "dimensionality_reduction"), \
        visualization_var = VIS_VAR
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
    VIS_VAR:
        Variable to visualize.
    Outputs
    -------
    ENRICHMENT_RESULTS_RDATA:
        Table or list of tables containing the gene sets passing the p-value threshold 
    ENRICHMENT_RESULTS_PDF:
        A bubble plot (or bubble plots) in one pdf file.

    """
    input:
        ranks = DIFFERENTIAL_INDEGREES_RANKED_RDATA, \
        genes = GENE_SET_FILE
    output:
        ENRICHMENT_RESULTS_RDATA
        ENRICHMENT_RESULTS_PDF
    message:
        "; Running GSEA."
    params:
        bin = os.path.join(config["bin"], "gsea"), \
        output_dir = os.path.join(BASE_OUTPUT_DIR, "differential_indegrees")
    shell:
        """
        echo "; I love snakemake" ;
        # Rscript {params.bin}/compute_gse.R \
        #     -i {input.ranks} \
        #     -g {input.genes} \
        #     -o {params.output_dir}
        """
    
    
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

## Directories
DATA_DIR = config["data_dir"]
BASE_OUTPUT_DIR = config["output_dir"]
DIFFERENTIAL_INDEGREE_OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, "differential_indegrees")
DIMENSIONALITY_REDUCTION_OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, "dimensionality_reduction")

## Input files ##
INPUT_METADATA = config["metadata_file"]
INPUT_INDEGREES = config["indegree_file"]
VIS_VAR = config["visualisation_var"]

# ## Output files ##
# DIFFERENTIAL_INDEGREES_TAB = os.path.join(DIFFERENTIAL_INDEGREE_OUTPUT_DIR, "differential_indegrees.tsv")
# DIFFERENTIAL_INDEGREES_RDATA = os.path.join(DIFFERENTIAL_INDEGREE_OUTPUT_DIR, "differential_indegrees.RData")

## Output plots ##
PCA_PLOT_PDF = os.path.join(DIMENSIONALITY_REDUCTION_OUTPUT_DIR, "pc12.pdf") # In the R script it's written as follows:  pdf(pcaplot, file.path(OUTPUT_DIR, paste0("PCA", VARIABLE, "pc12.pdf")))


## Rule ALL ##
rule all:
    input:
        DIFFERENTIAL_INDEGREES_TAB, \
        DIFFERENTIAL_INDEGREES_RDATA, \
        DIFFERENTIAL_INDEGREES_RANKED, \
        PCA_PLOT_PDF

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
        DIFFERENTIAL_INDEGREES_TAB, \
        DIFFERENTIAL_INDEGREES_RDATA \
        DIFFERENTIAL_INDEGREES_RANKED
    message:
        "; Running differential indegree computation on {input}."
    params:
        bin = os.path.join(config["bin"], "preprocessing"), \
        output_dir = os.path.join(BASE_OUTPUT_DIR, "differential_indegrees")
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
        indegrees = INPUT_INDEGREES, \
        visualization_var = VIS_VAR
    output:
        PCA_PLOT_PDF
    message:
        "; Running dimensionality reduction of indegrees on {input}."
    params:
        bin = os.path.join(config["bin"], "preprocessing"), \
        output_dir = os.path.join(BASE_OUTPUT_DIR, "dimensionality_reduction")
    shell:
        """
        Rscript {params.bin}/compute_diff_indegrees.R \
            -i {input.indegrees} \
            -m {input.metadata} \
            -v {input.visualization_var} \
            -o {params.output_dir}
        """
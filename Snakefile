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

## Input files ##
INPUT_METADATA = config["metadata_file"]
INPUT_INDEGREES = config["indegree_file"]

## Output files ##
DIFFERENTIAL_INDEGREES_TAB = os.path.join(DIFFERENTIAL_INDEGREE_OUTPUT_DIR, "differential_indegrees.tsv")
DIFFERENTIAL_INDEGREES_RDATA = os.path.join(DIFFERENTIAL_INDEGREE_OUTPUT_DIR, "differential_indegrees.RData")


## Rule ALL ##
rule all:
    input:
        DIFFERENTIAL_INDEGREES_TAB, \
        DIFFERENTIAL_INDEGREES_RDATA

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
        DIFFERENTIAL_INDEGREES_RDATA
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
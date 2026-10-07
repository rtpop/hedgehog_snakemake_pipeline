## ------------------------------------------------------------------------------------------- ##
## HOW TO RUN                                                                                  ##
## Run from directory containing Snakefile                                                     ##
## ------------------------------------------------------------------------------------------- ##
## For dry run                                                                                 ##
## snakemake --cores 1 -np                                                                     ##
## ------------------------------------------------------------------------------------------- ##
## For local run                                                                               ##
## snakemake --cores 1                                                                         ##
## ------------------------------------------------------------------------------------------- ##
## For running with singularity container                                                      ##
## snakemake --cores 1 --use-singularity --singularity-args '\-e' --cores 1                    ##
## ------------------------------------------------------------------------------------------- ## ------------------------------ ##
## For running in the background                                                                                                 ##
## nohup snakemake --cores 1 --use-singularity --singularity-args '\-e' 2>&1 > logs/snakemake_$(date +'%Y-%m-%d_%H-%M-%S').log & ##
## ----------------------------------------------------------------------------------------------------------------------------- ##

##-----------##
## Libraries ##
##-----------##

import os 
import sys
import glob
from pathlib import Path
import time

## ----------------- ##
## Global parameters ##
## ----------------- ##

configfile: "config.yaml"

NCORES = config["ncores"]

# Containers
PYTHON_CONTAINER = config["python_container"]
R_CONTAINER = config["r_container"]

# Directories
DATA_DIR = config["data_dir"]
OUTPUT_DIR = config["output_dir"]
SRC = config["src_dir"]
HEDGEHOG_DIR = os.path.join(OUTPUT_DIR, "{tissue_type}", "hedgehog_bug_fix_q_score")
BIHIDEF_RUN_DIR = os.path.join(HEDGEHOG_DIR, "bihidef_run")
GO_DIR = os.path.join(OUTPUT_DIR, "{tissue_type}", "go_enrichment")

# Other params
DELIMITER = config["delimiter"]
TISSUE = config["tissue"] # wildcard
TAR_TAG = config["target_tag"]

## ------------------------------ ##
## Download and process GTEx data ##
## ------------------------------ ##
GTEX_DATA_FILE = os.path.join(DATA_DIR, "download", config["gtex_data_file"])
PROCESSING_LOG = os.path.join(DATA_DIR, config["processing_log"])
MOTIF_PRIOR = os.path.join(DATA_DIR, config["motif_prior"])

##-------------------------------------##
## Filtering PANDA network for BiHiDef ##
##-------------------------------------##
PANDA_NET = os.path.join(DATA_DIR, "{tissue_type}", config["panda_net_file"])
FILTERING_METHOD = config["filtering_method"]

## ------------------- ##
## Filtering benchmark ##
## ------------------- ##
BENCHMARK = config["benchmark"]
BENCH_RESOLUTION = config["bench_resolution"]
UNFILTERED = config["plot_unfiltered"]

# set benchmarking params if benchmarking is enabled
if BENCHMARK:
    BENCHMARK_DIR = os.path.join(HEDGEHOG_DIR, "benchmark")
    MAX_COMMUNITIES = config["max_communities"]
    BENCH_RESOLUTION = config["bench_resolution"]
    UNFILTERED = config["plot_unfiltered"]
    FILTERING_BENCH = os.path.join(BENCHMARK_DIR, "filtering_benchmark_res_{bench_resolution}.txt")
    FILTERING_METHOD = "both"  # always run both for benchmarking
    FILTERING_BENCH_DF = os.path.join(BENCHMARK_DIR, "filtering_benchmark_df.txt")
    FILTERING_BENCH_CONSOLIDATED = os.path.join(OUTPUT_DIR, "filtering_benchmark_consolidated.txt")
    BENCHMARK_PLOT = os.path.join(BENCHMARK_DIR, "filtering_benchmark_plot.pdf")
    FILTERING_HEATMAP = os.path.join(OUTPUT_DIR, "filtering_benchmark_heatmaps")

# set filtering params
if FILTERING_METHOD == "both":
    PANDA_NET_FILTERED = [
        os.path.join(HEDGEHOG_DIR, "panda_network_filtered_prior.txt"),
        os.path.join(HEDGEHOG_DIR, "panda_network_filtered_hedgehog.txt")
    ]
elif FILTERING_METHOD == "prior":
    PANDA_NET_FILTERED = [os.path.join(HEDGEHOG_DIR, "panda_network_filtered_prior.txt")]
elif FILTERING_METHOD == "hedgehog":
    PANDA_NET_FILTERED = [os.path.join(HEDGEHOG_DIR, "panda_network_filtered_hedgehog.txt")]
else:
    raise ValueError("Unknown filtering method: {}".format(FILTERING_METHOD))

## ----------- ##
## Run BiHiDef ##
## ----------- ##
GENE_COMMUNITIES = os.path.join(BIHIDEF_RUN_DIR, TAR_TAG + ".nodes")
MAX_COMMUNITIES = config["max_communities"]
MAX_RESOLUTION = config["max_resolution"]
FILTERING_BIHIDEF = config["filtering_bihidef"]
REG_TAG = config["regulator_tag"]
TAR_TAG = config["target_tag"]

## ------------------ ##
## Select communities ##
## ------------------ ##
SELECTED_COMMUNITIES = os.path.join(BIHIDEF_RUN_DIR, TAR_TAG + "_selected_communities.gmt")
COMMUNITY_STATS = os.path.join(BIHIDEF_RUN_DIR, TAR_TAG + "_community_stats.txt")
MAX_GENES = config["max_genes"]
MIN_GENES = config["min_genes"]



# set filtering params
if FILTERING_METHOD == "both":
    PANDA_NET_FILTERED = [
        os.path.join(HEDGEHOG_DIR, "panda_network_filtered_prior.txt"),
        os.path.join(HEDGEHOG_DIR, "panda_network_filtered_hedgehog.txt")
    ]
elif FILTERING_METHOD == "prior":
    PANDA_NET_FILTERED = [os.path.join(HEDGEHOG_DIR, "panda_network_filtered_prior.txt")]
elif FILTERING_METHOD == "hedgehog":
    PANDA_NET_FILTERED = [os.path.join(HEDGEHOG_DIR, "panda_network_filtered_hedgehog.txt")]
else:
    raise ValueError("Unknown filtering method: {}".format(FILTERING_METHOD))

##-------##
## RULES ##
##-------##

rule all:
    input:
        expand(PANDA_NET_FILTERED, tissue_type=TISSUE),
        expand(GENE_COMMUNITIES, tissue_type=TISSUE),
        expand(SELECTED_COMMUNITIES, tissue_type=TISSUE),
        expand(COMMUNITY_STATS, tissue_type=TISSUE),
        expand(FILTERING_BENCH, tissue_type=TISSUE, bench_resolution=BENCH_RESOLUTION) if BENCHMARK else [],
        expand(FILTERING_BENCH_DF, tissue_type=TISSUE) if BENCHMARK else [],
        FILTERING_BENCH_CONSOLIDATED if BENCHMARK else []

## ---------------------------- ##
## Download & process GTEX data ##
## ---------------------------- ##

rule downaload_gtex:
    output:
        GTEX_DATA_FILE
    params:
        out_dir = os.path.join(DATA_DIR, "download"), \
        log_file = os.path.join(DATA_DIR, "download", "download_gtex.log")
    container:
        PYTHON_CONTAINER
    message:
        "; Downloading GTEx data."
    shell:
        """
        mkdir -p {params.out_dir}
        curl --output {output} https://zenodo.org/records/838734/files/GTEx_PANDA_tissues.RData?download=1 > {params.log_file} 2>&1
        """

rule process_gtex:
    input:
        gtex_data = GTEX_DATA_FILE
    output:
        output_log = PROCESSING_LOG,
        prior = MOTIF_PRIOR
    params:
        out_dir = DATA_DIR,
        extract_edges = True
    container:
        R_CONTAINER
    message:
        "; Processing GTEx data."
    script:
        os.path.join(SRC, "process_gtex.R")

##-------------------------------------##
## Filtering PANDA network for BiHiDef ##
##-------------------------------------##
rule process_and_filter_panda:
    input:
        panda = PANDA_NET,
        prior = MOTIF_PRIOR
    output:
        filtered_panda = PANDA_NET_FILTERED
    params:
        out_dir = os.path.join(HEDGEHOG_DIR),
        delimiter = DELIMITER,
        filtering_method = FILTERING_METHOD
    container:
        PYTHON_CONTAINER
    script:
        os.path.join(SRC, "filter_panda.py")

## ------------------- ##
## Benchmark filtering ##
## ------------------- ##

rule panda_filtering_benchmark:
    """
    This rule benchmarks filtering methods for the PANDA network.

    Inputs
    ------
    PANDA_NET:
        A TXT file with the PANDA network.
    MOTIF_PRIOR:
        A TXT file with the motif prior.
    ------
    Outputs
    -------
    BENCHMARK_FILTERED:
        A TXT file with the benchmark data filtered.
    """
    input:
        panda = PANDA_NET, \
        panda_filtered = PANDA_NET_FILTERED, \
        prior = MOTIF_PRIOR
    output:
        filtering_bench = FILTERING_BENCH
    params:
        script = os.path.join(SRC, "filter_benchmark.py"), \
        out_dir = BENCHMARK_DIR, \
        delimiter = DELIMITER, \
        resolution = '{bench_resolution}', \
        max_communities = MAX_COMMUNITIES
    container:
        PYTHON_CONTAINER
    script:
        params.script

rule consolidate_benchmark_resolutions:
    """
    This rule consolidates the different resolutions into one file per tissue.

    Inputs
    ------
    FILTERING_BENCH_LIST:
        A list of TXT files with the benchmark data for all resolutions.
    ------
    Outputs
    -------
    FILTERING_BENCH_DF:
        A TXT file with the consolidated benchmark data.
    """
    input:
        filtering_bench_list = lambda wildcards: expand(
            FILTERING_BENCH,
            tissue_type=wildcards.tissue_type,
            bench_resolution=BENCH_RESOLUTION
        )
    output:
        filtering_bench_df=FILTERING_BENCH_DF
    params:
        function_file=os.path.join(SRC, "consolidate_benchmark_fn.R"),
        tissue_type="{tissue_type}"
    container:
        R_CONTAINER
    message:
        "; Consolidating benchmark data for {wildcards.tissue_type}"
    script:
        os.path.join(SRC, "consolidate_resolutions.R")

rule consolidate_benchmark_all:
    """
    This rule consolidates all benchmark data into a single file.

    Inputs
    ------
    FILTERING_BENCH_DF:
        A list of TXT files with the benchmark data for all tissues and resolutions.
    ------
    Outputs
    -------
    FILTERING_BENCH_CONSOLIDATED:
        A TXT file with the consolidated benchmark data for all tissues.
    """
    input:
        filtering_bench_dfs = expand(FILTERING_BENCH_DF, tissue_type=TISSUE)
    output:
        filtering_benchmark_consolidated=FILTERING_BENCH_CONSOLIDATED
    params:
        function_file=os.path.join(SRC, "consolidate_benchmark_fn.R")
    container:
        R_CONTAINER
    message:
        "; Consolidating all benchmark data"
    script:
        os.path.join(SRC, "consolidate_benchmark.R")

## --------------- ##
## Running BiHiDeF ##
## --------------- ##

rule run_bihidef:
    """
    This rule runs the BiHiDeF algorithm.

    BiHiDeF is available at
    """
    input:
        net = PANDA_NET_FILTERED
    output:
        gene_communities = GENE_COMMUNITIES
    params:
        run_script = os.path.join(SRC, "run_bihidef.py"), \
        max_communities = MAX_COMMUNITIES, \
        max_resolution = MAX_RESOLUTION, \
        output_prefix_reg = REG_TAG, \
        output_prefix_tar = TAR_TAG, \
        outdir = BIHIDEF_RUN_DIR, \
        resource_log = os.path.join(BIHIDEF_RUN_DIR, "run_resources.log"), \
        filtering_method = FILTERING_BIHIDEF,
        ncores = NCORES
    container:
        PYTHON_CONTAINER
    script:
        params.run_script

## --------------------- ##
## Selecting communities ##
## --------------------- ##

rule select_communities:
    """
    This rule selects the communities from the BiHiDeF output and formats them as a GMT file.

    Inputs
    ------
    GENE_COMMUNITIES:
        A TXT file with the communities from BiHiDeF.
    ------
    Outputs
    -------
    SELECTED_COMMUNITIES:
        A TXT file with the selected communities.
    COMMUNITY_STATS:
        A TXT file with statistics about the communities.
    """
    input:
        communities = GENE_COMMUNITIES
    output:
        selected_communities = SELECTED_COMMUNITIES, \
        stats = COMMUNITY_STATS
    params:
        script = os.path.join(SRC, "select_communities.py"), \
        max_genes = MAX_GENES, \
        min_genes = MIN_GENES
    container:
        PYTHON_CONTAINER
    script:
        params.script
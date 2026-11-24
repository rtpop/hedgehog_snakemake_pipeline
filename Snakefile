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

# Existing CONDOR and Leiden scripts, using the same prior-filtered PANDA inputs.
LEIDEN_BENCHMARK_DIR = os.path.join(OUTPUT_DIR, "leiden_benchmark")
LEIDEN_SUMMARY = os.path.join(LEIDEN_BENCHMARK_DIR, "modularity_condor_leiden_gene_projection.tsv")
LEIDEN_TISSUE_SUMMARY = os.path.join(
    LEIDEN_BENCHMARK_DIR, "{tissue_type}", "modularity_condor_leiden_gene_projection.tsv"
)
LEIDEN_MEMBERSHIP = os.path.join(LEIDEN_BENCHMARK_DIR, "{tissue_type}", "leiden_output",
                                 "genes_memb_{leiden_resolution}.txt")
LEIDEN_RESOLUTIONS = config["leiden_resolutions"]
CONDOR_RESOLUTIONS = config["leiden_condor_resolutions"]
CONDOR_BENCHMARK_DIR = os.path.join(LEIDEN_BENCHMARK_DIR, config["leiden_condor_membership_dir"])
CONDOR_MEMBERSHIP = os.path.join(CONDOR_BENCHMARK_DIR, "genes_memb_{condor_resolution}.txt")
CONDOR_REG_MEMBERSHIP = os.path.join(CONDOR_BENCHMARK_DIR, "tf_memb_{condor_resolution}.txt")
CONDOR_SELECTED_COMMUNITIES = os.path.join(
    CONDOR_BENCHMARK_DIR, "gene_selected_communities_condor_{condor_resolution}.gmt"
)

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

JACCARD_HEATMAP = os.path.join(OUTPUT_DIR, "jaccard_heatmaps")
JACCARD_BOXPLOT = os.path.join(OUTPUT_DIR, config["jaccard_boxplot"])
JACCARD_BOXPLOT_TISSUES = config.get("jaccard_boxplot_tissues", TISSUE)
if not set(JACCARD_BOXPLOT_TISSUES).issubset(TISSUE):
    raise ValueError("jaccard_boxplot_tissues must be included in tissue")

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
COMMUNITY_STATS_CONSOLIDATED = os.path.join(OUTPUT_DIR, config["community_stats_consolidated"])
COMMUNITY_TISSUE_SUMMARY = os.path.join(OUTPUT_DIR, config["community_tissue_summary"])
COMMUNITY_PLOT = os.path.join(OUTPUT_DIR, config["community_plot"])

## --------------------------- ##
## GO enrichment and qscores ##
## --------------------------- ##
AUTO_BG = config["auto_bg"]
GENE_BACKGROUND = os.path.join(GO_DIR, "network_background.tsv")
GO_ENRICHMENT = os.path.join(GO_DIR, config["go_results"])
GO_DETAILS = os.path.join(GO_DIR, "community_terms")
GO_ENRICHMENT_PLOT = os.path.join(OUTPUT_DIR, config["go_enrichment_plot"])
QSCORE_TAR = os.path.join(BIHIDEF_RUN_DIR, "qscore_tar.csv")
QSCORE_REG = os.path.join(BIHIDEF_RUN_DIR, "qscore_reg.csv")
BIHIDEF_RESOURCE_LOG = os.path.join(BIHIDEF_RUN_DIR, "run_resources.log")
SELECTED_QSCORE_PLOT = os.path.join(BIHIDEF_RUN_DIR, "selected_community_qscores.pdf")
SELECTED_QSCORE_MATCHES = os.path.join(BIHIDEF_RUN_DIR, "selected_community_qscores.csv")
SELECTED_QSCORE_COVERAGE = os.path.join(BIHIDEF_RUN_DIR, "selected_community_qscore_coverage.csv")



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
        COMMUNITY_STATS_CONSOLIDATED,
        COMMUNITY_TISSUE_SUMMARY,
        COMMUNITY_PLOT,
        expand(GO_ENRICHMENT, tissue_type=TISSUE),
        GO_ENRICHMENT_PLOT,
        expand(SELECTED_QSCORE_PLOT, tissue_type=TISSUE),
        expand(SELECTED_QSCORE_MATCHES, tissue_type=TISSUE),
        expand(SELECTED_QSCORE_COVERAGE, tissue_type=TISSUE),
        LEIDEN_SUMMARY if BENCHMARK else [],
        expand(FILTERING_BENCH, tissue_type=TISSUE, bench_resolution=BENCH_RESOLUTION) if BENCHMARK else [],
        expand(FILTERING_BENCH_DF, tissue_type=TISSUE) if BENCHMARK else [],
        FILTERING_BENCH_CONSOLIDATED if BENCHMARK else [],
        FILTERING_HEATMAP if BENCHMARK else [],
        JACCARD_HEATMAP,
        JACCARD_BOXPLOT

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

rule run_condor:
    input:
        network = os.path.join(HEDGEHOG_DIR, "panda_network_filtered_prior.txt")
    output:
        target_membership = CONDOR_MEMBERSHIP,
        regulator_membership = CONDOR_REG_MEMBERSHIP,
        selected_communities = CONDOR_SELECTED_COMMUNITIES
    params:
        resolution = "{condor_resolution}",
        delimiter = DELIMITER,
        min_genes = MIN_GENES,
        max_genes = MAX_GENES
    threads: 1
    container:
        PYTHON_CONTAINER
    message:
        "Running CONDOR for {wildcards.tissue_type} at resolution {wildcards.condor_resolution}"
    script:
        os.path.join("src", "old", "hedgehog", "run_condor.py")


rule run_leiden:
    input:
        summaries = expand(LEIDEN_TISSUE_SUMMARY, tissue_type=TISSUE)
    output:
        summary = LEIDEN_SUMMARY
    threads: 1
    container:
        PYTHON_CONTAINER
    message:
        "Combining CONDOR-Leiden comparison summaries across tissues"
    script:
        os.path.join("src", "new", "consolidate_leiden_summaries.py")


rule run_leiden_tissue:
    input:
        networks = lambda wildcards: expand(
            os.path.join(HEDGEHOG_DIR, "panda_network_filtered_prior.txt"),
            tissue_type=[wildcards.tissue_type]
        ),
        condor_memberships = lambda wildcards: expand(
            CONDOR_MEMBERSHIP, tissue_type=[wildcards.tissue_type],
            condor_resolution=CONDOR_RESOLUTIONS
        )
    output:
        summary = LEIDEN_TISSUE_SUMMARY,
        memberships = expand(LEIDEN_MEMBERSHIP, tissue_type="{tissue_type}",
                             leiden_resolution=LEIDEN_RESOLUTIONS)
    params:
        tissues = lambda wildcards: [wildcards.tissue_type],
        condor_resolutions = CONDOR_RESOLUTIONS,
        leiden_resolutions = LEIDEN_RESOLUTIONS,
        delimiter = DELIMITER
    threads: 1
    container:
        config.get("leiden_container", PYTHON_CONTAINER)
    message:
        "Running the existing CONDOR-Leiden comparison for {wildcards.tissue_type}"
    script:
        os.path.join("src", "run_leiden.py")


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

rule plot_bench_heatmap:
    """
    This rule plots the benchmark heatmap.

    Inputs
    ------
    FILTERING_BENCH_CONSOLIDATED:
        A TXT file with the consolidated benchmark data.
    ------
    Outputs
    -------
    BENCHMARK_HEATMAP:
        A PDF file with the benchmark heatmap.
    """
    input:
        filtering_bench_df = FILTERING_BENCH_CONSOLIDATED
    output:
        benchmark_heatmap = directory(FILTERING_HEATMAP)
    params:
        function_file = os.path.join("src", "old", "analysis", "plot_filtering_bench_fn.R"), \
        metric = "Modularity",  \
        filtering_method = "all", \
        separate_files = True
    container:
        R_CONTAINER
    message:
        "; Plotting benchmark heatmap"
    script:
        os.path.join(SRC, "plot_filtering_bench_heatmap.R")

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
        gene_communities = GENE_COMMUNITIES,
        qscore_tar = QSCORE_TAR,
        qscore_reg = QSCORE_REG
    params:
        run_script = os.path.join(SRC, "run_bihidef.py"), \
        max_communities = MAX_COMMUNITIES, \
        max_resolution = MAX_RESOLUTION, \
        output_prefix_reg = REG_TAG, \
        output_prefix_tar = TAR_TAG, \
        outdir = BIHIDEF_RUN_DIR, \
        resource_log = BIHIDEF_RESOURCE_LOG,
        filtering_method = FILTERING_BIHIDEF,
        ncores = NCORES
    container:
        PYTHON_CONTAINER
    script:
        os.path.join(SRC, "run_bihidef.py")

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

## ---------------------------------- ##
## Compare communities across tissues ##
## ---------------------------------- ##
rule plot_jaccard_heatmaps:
    """
    This rule plots Jaccard heatmaps for all configured tissue pairs.

    Inputs
    ------
    SELECTED_COMMUNITIES:
        A list of BiHiDef community files for all configured tissues.
    ------
    Outputs
    -------
    JACCARD_HEATMAP:
        A directory containing one heatmap for each tissue pair.
    """
    input:
        selected_communities = expand(SELECTED_COMMUNITIES, tissue_type=TISSUE)
    output:
        jaccard_heatmaps = directory(JACCARD_HEATMAP)
    params:
        tissues = TISSUE
    container:
        R_CONTAINER
    message:
        "; Plotting Jaccard heatmaps for all configured tissue pairs"
    script:
        os.path.join(SRC, "plot_jaccard_heatmaps.R")

rule plot_jaccard_boxplot:
    input:
        selected_communities = expand(SELECTED_COMMUNITIES, tissue_type=JACCARD_BOXPLOT_TISSUES)
    output:
        jaccard_boxplot = JACCARD_BOXPLOT
    params:
        tissues = JACCARD_BOXPLOT_TISSUES,
        min_jaccard = config["jaccard_boxplot_min_jaccard"]
    threads: 1
    container:
        R_CONTAINER
    message:
        "Plotting selected-community Jaccard distributions for tissue pairs"
    script:
        os.path.join(SRC, "plot_jaccard_boxplot.R")


## ---------------------------- ##
## GO enrichment of communities ##
## ---------------------------- ##

rule derive_go_background:
    input:
        network = PANDA_NET
    output:
        background = GENE_BACKGROUND
    params:
        delimiter = DELIMITER
    container:
        PYTHON_CONTAINER
    message:
        "Deriving GO background from all PANDA targets for {wildcards.tissue_type}"
    script:
        os.path.join(SRC, "derive_go_background.py")


rule go_enrichment:
    input:
        gmt = SELECTED_COMMUNITIES,
        bg = [] if AUTO_BG else GENE_BACKGROUND
    output:
        go_enrichment = GO_ENRICHMENT,
        details = directory(GO_DETAILS)
    params:
        function_file = os.path.join(SRC, "GO_enrichment_fn.R"),
        auto_bg = AUTO_BG,
        save_all = config["save_all"],
        sig_thresh = config["sig_thresh"],
        statistic = config["statistic"],
        algorithm = config["algorithm"]
    container:
        R_CONTAINER
    message:
        "GO enrichment of selected communities for {wildcards.tissue_type}"
    script:
        os.path.join(SRC, "GO_enrichment.R")


rule plot_n_enriched_communities:
    input:
        tissue_summary = COMMUNITY_TISSUE_SUMMARY
    output:
        go_enrichment_plot = GO_ENRICHMENT_PLOT
    container:
        R_CONTAINER
    script:
        os.path.join(SRC, "go_enrichment_plot.R")


rule plot_selected_community_qscores:
    input:
        communities = SELECTED_COMMUNITIES,
        qscores = QSCORE_TAR
    output:
        plot = SELECTED_QSCORE_PLOT,
        matched_scores = SELECTED_QSCORE_MATCHES,
        coverage = SELECTED_QSCORE_COVERAGE
    params:
        function_file = os.path.join(SRC, "selected_qscores_fn.R")
    container:
        R_CONTAINER
    message:
        "Plotting qscores for exact selected-community matches in {wildcards.tissue_type}"
    script:
        os.path.join(SRC, "plot_selected_qscores.R")

## --------------- ##
## Community stats ##
## --------------- ##
rule consolidate_community_stats:
    input:
        community_stats = expand(COMMUNITY_STATS, tissue_type=TISSUE),
        selected_communities = expand(SELECTED_COMMUNITIES, tissue_type=TISSUE),
        go_summary = expand(GO_ENRICHMENT, tissue_type=TISSUE)
    output:
        consolidated_community_stats = COMMUNITY_STATS_CONSOLIDATED,
        tissue_summary = COMMUNITY_TISSUE_SUMMARY
    params:
        tissues = TISSUE
    threads: 1
    container:
        PYTHON_CONTAINER
    message:
        "Consolidating selected community statistics and GO enrichment across tissues"
    script:
        os.path.join(SRC, "consolidate_community_stats.py")


rule plot_n_communities:
    input:
        tissue_summary = COMMUNITY_TISSUE_SUMMARY
    output:
        community_plot = COMMUNITY_PLOT
    threads: 1
    container:
        R_CONTAINER
    message:
        "Plotting total and selected community counts across tissues"
    script:
        os.path.join(SRC, "community_plot.R")

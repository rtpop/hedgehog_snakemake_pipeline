"""Run CONDOR on prior-filtered PANDA networks and export selected communities."""
import os

import condor
import pandas as pd


## ------------------------------- ##
## Detection and community export ##
## ------------------------------- ##
def run_condor(network, target_output, regulator_output, selected_output,
               resolution, delimiter=",", min_genes=10, max_genes=200):
    for path in (target_output, regulator_output, selected_output):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    # Activate the detection call from the original script. Filtered PANDA files
    # have no header: header=None retains the first TF-gene edge as data.
    condor.run_condor(
        network_file=network,
        tar_output=target_output,
        reg_output=regulator_output,
        index_col=None,
        header=None,
        sep=delimiter,
        return_output=True,
        resolution=resolution,
    )

    # Preserve the original target-name cleanup, strict size bounds and export
    # format. The full target membership table is used by run_leiden.py; the
    # selected-community file is an additional output, not its scoring input.
    targets = pd.read_csv(target_output, sep=",", header=0)
    targets["tar"] = targets["tar"].str.replace("tar_", "", regex=False)
    selected = targets.groupby("community").filter(
        lambda group: min_genes < len(group) < max_genes
    )
    with open(selected_output, "w") as output:
        for community_id in sorted(selected["community"].unique()):
            genes = selected.loc[selected["community"] == community_id, "tar"].tolist()
            output.write(f"community_{community_id}: {' '.join(genes)}\n")


## --------------------- ##
## Snakemake entry point ##
## --------------------- ##
if __name__ == "__main__":
    if "snakemake" in globals():
        run_condor(
            network=str(snakemake.input.network),
            target_output=str(snakemake.output.target_membership),
            regulator_output=str(snakemake.output.regulator_membership),
            selected_output=str(snakemake.output.selected_communities),
            resolution=float(snakemake.params.resolution),
            delimiter=snakemake.params.delimiter,
            min_genes=int(snakemake.params.min_genes),
            max_genes=int(snakemake.params.max_genes),
        )
    else:
        # Original standalone tissue list, paths and resolutions.
        tissues = [
            "Adipose_subcutaneous", "Adipose_visceral", "Adrenal_gland", "Artery_aorta",
            "Artery_coronary", "Artery_tibial", "Brain_basal_ganglia", "Brain_cerebellum", "Brain_other",
            "Breast", "Colon_sigmoid", "Colon_transverse", "Esophagus_mucosa", "Esophagus_muscularis",
            "Heart_atrial_appendage", "Heart_left_ventricle", "Intestine_terminal_ileum",
            "Kidney_cortex", "Liver", "Lung", "Minor_salivary_gland", "Ovary", "Pancreas", "Pituitary",
            "Prostate", "Skin", "Spleen", "Stomach", "Testis", "Thyroid", "Tibial_nerve", "Uterus",
            "Vagina", "Whole_blood",
        ]
        for resolution in (1.25, 1.5, 1.75):
            for tissue in tissues:
                directory = os.path.join("output", tissue, "condor_output")
                run_condor(
                    network=os.path.join("output", tissue, "hedgehog", "panda_network_filtered_prior.txt"),
                    target_output=os.path.join(directory, f"genes_memb_{resolution}.txt"),
                    regulator_output=os.path.join(directory, f"tf_memb_{resolution}.txt"),
                    selected_output=os.path.join(directory, f"gene_selected_communities_condor_{resolution}.gmt"),
                    resolution=resolution,
                )

from netZooPy import condor
import pandas as pd
import os
tissues = ["Adipose_subcutaneous", "Adipose_visceral", "Adrenal_gland", "Artery_aorta",
    "Artery_coronary", "Artery_tibial", "Brain_basal_ganglia", "Brain_cerebellum", "Brain_other",
    "Breast", "Colon_sigmoid", "Colon_transverse", "Esophagus_mucosa", "Esophagus_muscularis",
    "Heart_atrial_appendage","Heart_left_ventricle", "Intestine_terminal_ileum",
    "Kidney_cortex", "Liver", "Lung", "Minor_salivary_gland", "Ovary", "Pancreas", "Pituitary",
    "Prostate", "Skin", "Spleen","Stomach", "Testis", "Thyroid", "Tibial_nerve","Uterus", "Vagina", "Whole_blood"]
data_dir = "output"
# run condor on prior filtered panda nets

resolution = [1.25, 1.5, 1.75]

for res in resolution:
    # for tissue in tissues:
    #     net_path = os.path.join(data_dir, tissue, "hedgehog", "panda_network_filtered_prior.txt")
    #     output_dir = os.path.join(data_dir, tissue, "condor_output")
    #     tar_output = os.path.join(output_dir, f"genes_memb_{res}.txt")
    #     reg_output = os.path.join(output_dir, f"tf_memb_{res}.txt")
    #     reg_comm_output = os.path.join(output_dir, f"tf_comm_{res}.gmt")
    #     gene_comm_output = os.path.join(output_dir, f"gene_comm_{res}.gmt")
    #     if not os.path.exists(output_dir):
    #         os.makedirs(output_dir)
    #     #run condor
    #     print(f"Running condor for tissue: {tissue}")
    #     condor_object = condor.run_condor(network_file=net_path, tar_output=tar_output, reg_output=reg_output,index_col=None,header = 0, sep = ",",return_output = True, resolution = res)
    #     print(f"Saving condor outputs for tissue: {tissue}")
    for tissue in tissues:
        # Read the file
        directory = os.path.join(data_dir, tissue, "condor_output")
        tar_file = os.path.join(directory, f"genes_memb_{res}.txt")
        tar_comm = pd.read_csv(tar_file, sep=",", header=0)
        tar_file_to_save = os.path.join(directory, f"gene_selected_communities_condor_{res}.gmt")
        # process gene names to remove tar_ prefix
        tar_comm['tar'] = tar_comm['tar'].str.replace('tar_', '', regex=False)
        # select only communities with more than 10 and less than 200 genes
        filtered_communities = tar_comm.groupby('community').filter(lambda x: 10 < len(x) < 200)
        # Write communities to file
        with open(tar_file_to_save, 'w') as f:
            for community_id in sorted(filtered_communities['community'].unique()):
                genes = filtered_communities[filtered_communities['community'] == community_id]['tar'].tolist()
                f.write(f"community_{community_id}: {' '.join(genes)}\n")

import numpy as np
import pandas as pd
import scipy.stats as stats
from sklearn import metrics
import os
from itertools import combinations

tissues=["Adipose_subcutaneous","Adipose_visceral","Adrenal_gland","Artery_aorta","Artery_coronary","Artery_tibial","Brain_basal_ganglia","Brain_cerebellum","Brain_other","Breast","Colon_sigmoid","Colon_transverse","Esophagus_mucosa","Esophagus_muscularis","Heart_atrial_appendage","Heart_left_ventricle","Intestine_terminal_ileum","Kidney_cortex","Liver","Lung","Minor_salivary_gland","Ovary","Pancreas","Pituitary","Prostate","Skin","Spleen","Stomach","Testis","Thyroid","Tibial_nerve","Uterus","Vagina","Whole_blood"]
data_dir="output"
resolutions=[1,1.25, 1.5, 1.75, 2,5]
heatmap_dir=os.path.join(data_dir,"filtering_benchmark_heatmaps")
os.makedirs(heatmap_dir,exist_ok=True)
for res in resolutions:
    print(f"=== Processing resolution {res} ===")
    pooled_path=os.path.join(data_dir,f"condor_pooled_communities_res{res}.tsv")
    if not os.path.exists(pooled_path):
        pooled_data=[]
        missing_tissues=[]
        for tissue in tissues:
            gmt_path=os.path.join(data_dir,tissue,"condor_output",f"gene_selected_communities_condor_{res}.gmt")
            if not os.path.exists(gmt_path):
                missing_tissues.append(tissue)
                continue
            data=[]
            with open(gmt_path,"r") as f:
                for line in f:
                    line=line.strip()
                    if not line:
                        continue
                    parts=line.split()
                    community_id=parts[0].rstrip(":")
                    genes=parts[1:]
                    if not genes:
                        continue
                    data.append({"community":community_id,"genes":" ".join(genes),"tissue":tissue})
            if data:
                df=pd.DataFrame(data)
                pooled_data.append(df)
        if missing_tissues:
            print(f"missing GMT for res {res} in tissues: {','.join(missing_tissues)}")
        if not pooled_data:
            print(f"no data found for resolution {res}; skipping")
            continue
        pooled_df=pd.concat(pooled_data,ignore_index=True)
        pooled_df.to_csv(pooled_path,sep="\t",index=False)
    else:
        pooled_df=pd.read_csv(pooled_path,sep="\t")
    all_genes=set()
    for genes_str in pooled_df["genes"]:
        all_genes.update(str(genes_str).split())
    total_genes=len(all_genes)
    if total_genes==0:
        print(f"no genes found for resolution {res}; skipping")
        continue
    results=[]
    for tissue1,tissue2 in combinations(tissues,2):
        df1=pooled_df[pooled_df["tissue"]==tissue1]
        df2=pooled_df[pooled_df["tissue"]==tissue2]
        if df1.empty or df2.empty:
            continue
        for _,row1 in df1.iterrows():
            genes1=set(str(row1["genes"]).split())
            if not genes1:
                continue
            for _,row2 in df2.iterrows():
                genes2=set(str(row2["genes"]).split())
                if not genes2:
                    continue
                intersection=genes1.intersection(genes2)
                union=genes1.union(genes2)
                if not union:
                    continue
                jaccard_index=len(intersection)/len(union)
                contingency_table=np.array([[len(intersection),len(genes1)-len(intersection)],[len(genes2)-len(intersection),total_genes-len(union)]])
                try:
                    _,p_value=stats.fisher_exact(contingency_table,alternative="greater")
                except Exception:
                    p_value=1.0
                results.append({"resolution":res,"tissue1":tissue1,"community1":row1["community"],"tissue2":tissue2,"community2":row2["community"],"jaccard_index":jaccard_index,"fisher_p_value":p_value})
    if not results:
        print(f"no comparisons for resolution {res}; skipping outputs")
        continue
    results_df=pd.DataFrame(results)
    p_value_adj=stats.false_discovery_control(results_df["fisher_p_value"],method="bh")[0]
    results_df["fisher_p_value_adj"]=p_value_adj
    all_out_path=os.path.join(data_dir,f"community_comparisons_condor_res{res}.tsv")
    results_df.to_csv(all_out_path,sep="\t",index=False)
    significant_results_df=results_df[results_df["jaccard_index"]>0.7]
    sig_out_path=os.path.join(data_dir,f"significant_community_comparisons_fisher_condor_res{res}_jaccard.tsv")
    significant_results_df.to_csv(sig_out_path,sep="\t",index=False)
    print(f"wrote {len(results_df)} comparisons to {all_out_path}")
    print(f"wrote {len(significant_results_df)} significant comparisons to {sig_out_path}")
    if significant_results_df.empty:
        print(f"no significant overlaps for resolution {res}; skipping heatmap")
        continue

# import numpy as np
# import pandas as pd
# import scipy.stats as stats
# from sklearn import metrics
# import os
# from itertools import combinations

# # tissues to compare
# tissues = ["Adipose_subcutaneous", "Adipose_visceral", "Adrenal_gland", "Artery_aorta",
#     "Artery_coronary", "Artery_tibial", "Brain_basal_ganglia", "Brain_cerebellum", "Brain_other",
#     "Breast", "Colon_sigmoid", "Colon_transverse", "Esophagus_mucosa", "Esophagus_muscularis",
#     "Heart_atrial_appendage","Heart_left_ventricle", "Intestine_terminal_ileum",
#     "Kidney_cortex", "Liver", "Lung", "Minor_salivary_gland", "Ovary", "Pancreas", "Pituitary",
#     "Prostate", "Skin", "Spleen","Stomach", "Testis", "Thyroid", "Tibial_nerve","Uterus", "Vagina", "Whole_blood"]
# data_dir = "output"

# # create one dataframe with all communities
# if not os.path.exists(os.path.join(data_dir, "condor_pooled_communities.tsv")):
#     pooled_data = []
#     for tissue in tissues:
#         with open(os.path.join(data_dir, tissue, "condor_output/gene_selected_communities_condor.gmt"), 'r') as f:
#             data = []
#             for line in f:
#                 parts = line.strip().split(' ')
#                 genes = parts[1:]
#                 data.append({'community': parts[0], 'genes': genes, "tissue": tissue})
#         df = pd.DataFrame(data)
#         pooled_data.append(df)
#     pooled_df = pd.concat(pooled_data, ignore_index=True)
#     pooled_df.to_csv(os.path.join(data_dir, "condor_pooled_communities.tsv"), sep="\t", index=False)

# # compare communities across tissues
# pooled_df = pd.read_csv(os.path.join(data_dir, "condor_pooled_communities.tsv"), sep="\t")
# # get all genes across all tissues
# all_genes = set()
# for genes_str in pooled_df['genes']:
#     all_genes.update(genes_str.split())

# total_genes = len(all_genes)

# results = []
# for (tissue_pair) in combinations(tissues, 2):
#     tissue1, tissue2 = tissue_pair
#     df1 = pooled_df[pooled_df['tissue'] == tissue1]
#     df2 = pooled_df[pooled_df['tissue'] == tissue2]
#     for _, row1 in df1.iterrows():
#         genes1 = set(row1['genes'].split())
#         for _, row2 in df2.iterrows():
#             genes2 = set(row2['genes'].split())
#             intersection = genes1.intersection(genes2)
#             union = genes1.union(genes2)
#             jaccard_index = len(intersection) / len(union) if len(union) > 0 else 0
#             # check if overlap is significant more than expected by chance with fisher's exact test
#             contingency_table = np.array([[len(intersection), len(genes1) - len(intersection)],
#                                            [len(genes2) - len(intersection), total_genes - len(union)]])
#             _, p_value = stats.fisher_exact(contingency_table, alternative='greater')
#             results.append({
#                 'tissue1': tissue1,
#                 'community1': row1['community'],
#                 'tissue2': tissue2,
#                 'community2': row2['community'],
#                 'jaccard_index': jaccard_index,
#                 'fisher_p_value': p_value
#                 })

# results_df = pd.DataFrame(results)
# # adjust p-value for multiple testing using BH FDR
# p_value_adj = stats.false_discovery_control(results_df['fisher_p_value'], method='bh')[0]

# results_df['fisher_p_value_adj'] = p_value_adj
# results_df.to_csv(os.path.join(data_dir, "community_comparisons_condor.tsv"), sep="\t", index=False)

# # sort to only keep significant overlaps
# significant_results_df = results_df[(results_df['fisher_p_value_adj'] <= 0.05)]
# significant_results_df.to_csv(os.path.join(data_dir, "significant_community_comparisons_fisher_condor.tsv"), sep="\t", index=False)

# # plot heatmap of jaccard indices clustered by tissue

# import seaborn as sns
# import matplotlib.pyplot as plt

# heatmap_data = significant_results_df.pivot_table(index=['tissue1', 'community1'],
#                                                   columns=['tissue2', 'community2'],
#                                                   values='jaccard_index',
#                                                   fill_value=0)
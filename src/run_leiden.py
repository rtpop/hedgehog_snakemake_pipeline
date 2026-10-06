#!/usr/bin/env python
import os
import pandas as pd
import igraph as ig
import leidenalg as la
from scipy.sparse import coo_matrix
tissues=["Adipose_subcutaneous","Adipose_visceral","Adrenal_gland","Artery_aorta","Artery_coronary","Artery_tibial","Brain_basal_ganglia","Brain_cerebellum","Brain_other","Breast","Colon_sigmoid","Colon_transverse","Esophagus_mucosa","Esophagus_muscularis","Heart_atrial_appendage","Heart_left_ventricle","Intestine_terminal_ileum","Kidney_cortex","Liver","Lung","Minor_salivary_gland","Ovary","Pancreas","Pituitary","Prostate","Skin","Spleen","Stomach","Testis","Thyroid","Tibial_nerve","Uterus","Vagina","Whole_blood"]
data_dir="output"
condor_resolutions=[1.25,1.5,1.75]
leiden_resolutions=[1.25,1.5,1.75]
summary_out=os.path.join(data_dir,"modularity_condor_leiden_gene_projection.tsv")

def build_gene_projection(panda_path):
    df=pd.read_csv(panda_path,header=None,names=["tf","gene","weight"])
    df["tf"]=df["tf"].astype(str)
    df["gene"]=df["gene"].astype(str)
    tfs=sorted(df["tf"].unique())
    genes=sorted(df["gene"].unique())
    tf_to_idx={t:i for i,t in enumerate(tfs)}
    gene_to_idx={g:i for i,g in enumerate(genes)}
    row=df["tf"].map(tf_to_idx).to_numpy()
    col=df["gene"].map(gene_to_idx).to_numpy()
    data=df["weight"].to_numpy()
    B=coo_matrix((data,(row,col)),shape=(len(tfs),len(genes)))
    W=(B.T@B).tocsr()
    W.setdiag(0)
    W.eliminate_zeros()
    coo=W.tocoo()
    mask=coo.row<coo.col
    rows=coo.row[mask]
    cols=coo.col[mask]
    weights=coo.data[mask]
    g=ig.Graph()
    g.add_vertices(len(genes))
    g.vs["name"]=genes
    edges=list(zip(rows.tolist(),cols.tolist()))
    if len(edges)>0:
        g.add_edges(edges)
        g.es["weight"]=weights.tolist()
    else:
        g.es["weight"]=[]
    return g

def load_condor_gene_membership(tissue,res):
    condor_dir=os.path.join(data_dir,tissue,"condor_output")
    memb_path=os.path.join(condor_dir,f"genes_memb_{res}.txt")
    if not os.path.exists(memb_path):
        return None
    df=pd.read_csv(memb_path,sep=",",header=0)
    df["gene"]=df["tar"].str.replace("tar_","",regex=False)
    df=df[["gene","community"]].drop_duplicates(subset="gene")
    membership=dict(zip(df["gene"],df["community"].astype(int)))
    return membership

def modularity_from_membership(graph,membership_dict):
    names=graph.vs["name"]
    max_comm=max(membership_dict.values()) if membership_dict else 0
    singleton_id=max_comm+1
    membership=[]
    for name in names:
        if name in membership_dict:
            membership.append(membership_dict[name])
        else:
            membership.append(singleton_id)
            singleton_id+=1
    return graph.modularity(membership)

def run_leiden_and_modularity(graph,resolution):
    partition=la.find_partition(graph,la.CPMVertexPartition,resolution_parameter=resolution,weights=graph.es["weight"] if "weight" in graph.es.attribute_names() else None)
    membership=partition.membership
    mod = graph.modularity(membership)
    return membership,mod

def main():
    rows=[]
    for tissue in tissues:
        print(f"Processing tissue: {tissue}")
        panda_path=os.path.join(data_dir,tissue,"hedgehog","panda_network_filtered_prior.txt")
        if not os.path.exists(panda_path):
            print(f"  PANDA file not found, skipping: {panda_path}")
            continue
        print("  Building gene projection graph...")
        g=build_gene_projection(panda_path)
        for res in condor_resolutions:
            condor_memb=load_condor_gene_membership(tissue,res)
            if condor_memb is None:
                print(f"  Condor membership not found for res={res}, skipping condor modularity.")
                condor_mod=None
            else:
                condor_mod=modularity_from_membership(g,condor_memb)
            rows.append({"tissue":tissue,"method":"condor","resolution":res,"modularity":condor_mod})
        for res in leiden_resolutions:
            print(f"  Running Leiden for resolution={res}...")
            leiden_memb, leiden_mod = run_leiden_and_modularity(g, res)

            # save Leiden partition (one file per tissue x resolution)
            leiden_dir = os.path.join(data_dir, tissue, "leiden_output")
            os.makedirs(leiden_dir, exist_ok=True)
            leiden_path = os.path.join(leiden_dir, f"genes_memb_{res}.txt")

            leiden_df = pd.DataFrame({
                "gene": g.vs["name"],
                "community": leiden_memb
            })
            leiden_df.to_csv(leiden_path, sep=",", index=False)

            # keep modularity summary as before
            rows.append({
                "tissue": tissue,
                "method": "leiden",
                "resolution": res,
                "modularity": leiden_mod
            })
    if rows:
        out_df=pd.DataFrame(rows)
        out_df.to_csv(summary_out,sep="\t",index=False)
        print(f"Saved modularity summary to: {summary_out}")
    else:
        print("No results to save (no tissues processed).")
        
if __name__=="__main__":
    main()

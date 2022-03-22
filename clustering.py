#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function to perform hierarchical clustering analysis on each dataset using in turn the 
three signatures

input parameters:
data           = dataframe of the dataset to cluster
data_name      = name of the dataset to cluster
signature      = list of genes to use as signature
signature name = name of the signatures
subtype        = list of subtype labels of the samples in data
z_score        = True/False
plot           = True/False

variable needed:
df_with_subtypes = dictionary of the dataset whose subtype information is available, where keys are the data name and 
                   values are tuples (df, subtypes), e.g.: {"data X": (discovery_dataset_used_for_signature_X, subtypes_related_to_signature_X),...}
"""

def clustering(data, data_name, signature, signature_name, subtype, z_score, plot):

    import pandas as pd
    from scipy.cluster.hierarchy import fcluster
    from scipy.stats import zscore
    import seaborn as sns
    import matplotlib.pyplot as plt        
            
    # if i want z-score
    if z_score: 
        data = data.apply(zscore)
    else:
        data = data.copy()
        
    print("Cluster:  ",data_name+" data - ",signature_name," signatures")
    
    # becomes True in case we know the real subtypes of the current dataset and we want them compared with the clusters
    compareWreal = False    
    if data_name in df_with_subtypes.keys():
        real_subtype = df_with_subtypes[data_name][1]
        compareWreal = True 
        
    nclust = len(set(subtype))
    print("Find ->", nclust, " clusters")
    
    # filter genes of the current classifier     
    data = data[list(set(data.columns).intersection(signature))]

    # print the amount of genes in common between classifier and dataset
    overl_size = data.shape[1]
    print("\nOverlapping genes: ", overl_size, "/", len(signature))

    # find clusters
    g = sns.clustermap(data, metric="correlation", z_score=1); plt.close();
    clust_lab = fcluster(g.dendrogram_row.linkage, t=nclust, criterion='maxclust')

    if plot:
        # assign colors to clusters
        data["Clusters"] = ["plum" if i==1 else "cadetblue" if i==2 else "gold" if i==3 else "skyblue" for i in clust_lab] 
    
        if compareWreal:
            data["Real"]  = list(real_subtype)
    
            labels_real = {"Basal":"orange", "Classical":"dodgerblue", "Classical PDA":"mediumvioletred", "Exocrine-like PDA":"dodgerblue", "QM-PDA":"lightgrey",
                           "ADEX":"mediumvioletred", "Pancreatic Progenitor":"dodgerblue", "Squamous":"lightgrey", "Immunogenic":"orange",
                           "KRT81-pos.":"dodgerblue", "HNF1A-pos.":"mediumvioletred", "Double negative":"lightgrey"}
    
            labels_real_name = list(data["Real"])       
            data["Real"] = data.Real.map(labels_real)
                        
            data = data.sort_values("Clusters")
            g=sns.clustermap(data = data.drop(["Clusters", "Real"],axis=1), cmap="Spectral", figsize=(7,7),
                           cbar_kws = {"shrink": 0.5}, xticklabels=False, yticklabels=False, row_cluster=False, col_cluster=True, z_score=1,
                           row_colors=pd.DataFrame([data["Clusters"],data["Real"]]).T, colors_ratio=.03, metric="correlation")
        
            for label in set(labels_real_name):
                g.ax_row_dendrogram.bar(0, 0, color=labels_real[label], label=label, linewidth=0)
            g.ax_row_dendrogram.legend(loc="best", ncol=1, facecolor="white", fontsize=10)
            
        else:
            data = data.sort_values("Clusters")
            g=sns.clustermap(data = data.drop(["Clusters"],axis=1), cmap="Spectral", figsize=(7,7),
                           cbar_kws = {"shrink": 0.3}, xticklabels=False, yticklabels=False, row_cluster=False, col_cluster=True, z_score=1,
                           row_colors=data["Clusters"], colors_ratio=.03, metric="correlation")
            
        g.fig.suptitle(signature_name+" signatures on "+data_name+" dataset", y=1.04, x=.58, fontsize=14)
        ax = g.ax_heatmap
        ax.set_xlabel("Signatures", fontsize=12)
        ax.set_ylabel("Samples", fontsize=12)
        plt.show()    
        return(clust_lab)
        
    else:
        return(clust_lab)
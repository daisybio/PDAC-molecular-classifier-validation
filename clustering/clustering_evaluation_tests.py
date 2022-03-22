#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Variables needed:
signature_dict = dictionary with list of genes of the signatures and subtypes,
                 e.g. {"Signature X": (list_genes_signature_X, subtypes_related_to_signature_X),...}
data_dict      = dictionary containing all the dataframes of the datasets we want to make predictions on 
                  (e.g.: {"data X":data_X_df, "data Y":data_y_df, ...})
z_score = True/False
"""

def test_clusters(data, data_name, signature, signature_name, subtype, z_score):
    
    from scipy.cluster.hierarchy import fcluster
    from scipy.stats import zscore
    import numpy as np
    
    # if i want z-score
    if z_score: 
        data = data.apply(zscore)
    else:
        data = data.copy()
        
    print("Cluster:  ",data_name+" data - ",signature_name," signatures")
        
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
    
    data.index  = list(map(lambda x: str(x), clust_lab))
    all_pval_list = []
    
    if nclust==2 and min([list(data.index).count(x) for x in set(data.index)]) >1:
        pval = []
        for i in data.columns:
            F, p = stats.ranksums(data.loc["1", i], data.loc["2", i])
            pval.append(p)
        a, correct_pval, b, c = multipletests(pvals=pval, alpha=0.05, method="bonferroni")
        return(-np.log10(correct_pval))
    
    elif nclust==3 and min([list(data.index).count(x) for x in set(data.index)]) >1:
        pval = []
        for i in data.columns:
            F, p = stats.kruskal(data.loc["1", i], data.loc["2", i], data.loc["3", i])
            pval.append(p)
        a, correct_pval, b, c = multipletests(pvals=pval, alpha=0.05, method="bonferroni")
        return(-np.log10(correct_pval))
    
    elif nclust==4 and min([list(data.index).count(x) for x in set(data.index)]) >1:
        pval = []
        for i in data.columns:
            F, p = stats.kruskal(data.loc["1", i], data.loc["2", i], data.loc["3", i], data.loc["4", i])
            pval.append(p)
        a, correct_pval, b, c = multipletests(pvals=pval, alpha=0.05, method="bonferroni")
        return(-np.log10(correct_pval))
    
    elif min([list(data.index).count(x) for x in set(data.index)]) == 1:
        cluster_1_sample = min([(x,list(data.index).count(x)) for x in set(data.index)], key = lambda t: t[1])[0]
        data = data.drop([cluster_1_sample], axis=0)
        data_clusters_left = list(set(data.index))
        pval = []
        for i in data.columns:
            F, p = stats.kruskal(data.loc[data_clusters_left[0], i], data.loc[data_clusters_left[1], i], data.loc[data_clusters_left[2], i])
            pval.append(p)
        a, correct_pval, b, c = multipletests(pvals=pval, alpha=0.05, method="bonferroni")
        return(-np.log10(correct_pval))

# =============================================================================
#  Run function and generate boxplot to attach below the clustering heatmaps
# =============================================================================

pval_distr_dict = {k:[] for k in data_dict.keys()}
for data_name, data in data_dict.items():
    for signature_name, signature in signature_dict.items():

        pval_distr = test_clusters(data  = data,
                      data_name       = data_name,
                      signature       = signature[0], 
                      signature_name  = signature_name, 
                      subtype         = signature[1], 
                      z_score         = True)
        
        pval_distr_dict[data_name].append(pval_distr)
 
import matplotlib.pyplot as plt
       
for k in pval_distr_dict.keys():
    for i in list(enumerate(["Moffitt", "Collisson", "Bailey"])):
        plt.figure(figsize=(5,1.2))
        plt.boxplot(pval_distr_dict[k][i[0]], vert=False); plt.xlabel("-log10(adjusted pvalue)", fontsize=16); plt.xticks(fontsize=14);
        plt.axvline(-np.log10(0.05), color="darkred", linestyle='--', lw=1.2); plt.title(str(i[1]))
        
    

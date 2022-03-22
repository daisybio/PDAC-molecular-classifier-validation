#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Variables needed:
signature_dict = dictionary with list of genes of the signatures and subtypes,
                 e.g. {"Signature X": (list_genes_signature_X, subtypes_related_to_signature_X),...}
data_dict      = dictionary containing all the dataframes of the datasets we want to make predictions on 
                  (e.g.: {"data X":data_X_df, "data Y":data_y_df, ...})
z_score = True/False
n_times = number of times clustering is repeated
measure = type of measure used to compare clusters
"""

def clusters_comparison(signature_dict, data_dict, z_score, n_times, measure):
    
    import itertools
    import zscore
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import fcluster
    import pandas as pd
    
    compare_dict = {k:{k:[] for k in data_dict.keys()} for k in signature_dict.keys()}
    compare_dict_random = {k:{k:[] for k in data_dict.keys()} for k in signature_dict.keys()}
    
    # if i want z-score
    if z_score: 
        data_dict = {k: v.apply(zscore) for k,v in data_dict.items()}
                
    clusters_dict = {k:{k:[] for k in data_dict.keys()} for k in signature_dict.keys()}
    # compute clusters
    for signature_name, signature_ in signature_dict.items():
        for data_name, data in data_dict.items():
            signature = signature_[0] 
            subtype   = signature_[1]

            data_original = data.copy()

            nclust = len(set(subtype))

            # filter genes of the current classifier     
            in_common = list(set(data.columns).intersection(signature))
            data = data[in_common]

            # predict subtypes with real signatures
            g = sns.clustermap(data, metric="correlation"); plt.close();
            clust_lab = fcluster(g.dendrogram_row.linkage, t=nclust, criterion='maxclust')
            clusters_dict[signature_name][data_name] = clust_lab
            
    for n in range(n_times):
        for signature_name, signature_ in signature_dict.items():
            for data_name, data in data_dict.items():
                signature = signature_[0] 
                subtype   = signature_[1]

                data_original = data.copy()

                nclust = len(set(subtype))

                # filter genes of the current classifier     
                in_common = list(set(data.columns).intersection(signature))
                data = data[in_common]

                # predict subtypes with random signatures
                data_random_genes = data_original[random.sample(list(set(data_original.columns)), len(in_common))]          # -set(signature)                      
                g_r = sns.clustermap(data_random_genes, metric="correlation"); plt.close();
                random_clust_lab = fcluster(g_r.dendrogram_row.linkage, t=nclust, criterion='maxclust')            

                compare_dict[signature_name][data_name].append(measure(clusters_dict[signature_name][data_name], random_clust_lab))
                compare_dict_random[signature_name][data_name].append(random_clust_lab)
                
    moff_res = pd.DataFrame(compare_dict["Moffitt et al."])
    coll_res = pd.DataFrame(compare_dict["Collisson et al."])
    bailey_res = pd.DataFrame(compare_dict["Bailey et al."])
    
    # =============================================================================
    # add pairwised random comparison
    # =============================================================================
    
    # compute pairwise score
    compute_scores = {k:{k:[] for k in data_dict.keys()} for k in signature_dict.keys()}
    
    for sign_key in compute_scores.keys():
        for data_key in compute_scores[sign_key].keys():
            to_df = pd.DataFrame(compare_dict_random[sign_key][data_key]).T
            pair_list = list(itertools.combinations(to_df.columns, 2))
            for pair in pair_list:
                a = to_df[list(pair)[0]]
                b = to_df[list(pair)[1]]
                
                pairwise_score = measure(a,b)
                compute_scores[sign_key][data_key].append(pairwise_score) #list of pairwised score between clusters obtained using random genes
                
    moff_res_rand = pd.DataFrame(compute_scores["Moffitt et al."])
    coll_res_rand = pd.DataFrame(compute_scores["Collisson et al."])
    bailey_res_rand = pd.DataFrame(compute_scores["Bailey et al."])
    
    data_moff   = []
    data_coll   = []
    data_bailey = []
    for i in data_dict.keys():
        data_moff.append(moff_res[i]); data_moff.append(moff_res_rand[i]) #list of scores: Real vs Random - Random vs Random
        data_coll.append(coll_res[i]); data_coll.append(coll_res_rand[i])
        data_bailey.append(bailey_res[i]); data_bailey.append(bailey_res_rand[i])

    return moff_res, coll_res, bailey_res, moff_res_rand, coll_res_rand, bailey_res_rand, data_moff, data_coll, data_bailey

# =============================================================================
#  Run function
# =============================================================================

from sklearn.metrics.cluster import adjusted_rand_score
moff_res, coll_res, bailey_res, moff_res_rand, coll_res_rand, bailey_res_rand, data_moff, data_coll, data_bailey = clusters_comparison(signature_dict, data_dict, 
                                                                                                                    z_score=True, n_times=1000) 
# =============================================================================
#  Figure S3
# =============================================================================
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D 
import numpy as np
                                                                             
for s in [("Moffitt",data_moff),("Collisson",data_coll),("Bailey",data_bailey)]:

    title_bp = s[0]
    data_bp  = s[1]

    fig, ax1 = plt.subplots(figsize=(13, 6))
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

    bp = ax1.boxplot(data_bp, notch=0, sym='o', vert=1, whis=1.5)
    plt.setp(bp['boxes'], color='grey', linewidth=1)
    plt.setp(bp['medians'], color='grey', linewidth=1)
    plt.setp(bp['whiskers'], color='grey', linewidth=1)
    plt.setp(bp['fliers'], color='grey', marker='o', markersize=4, alpha=0.5)

    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

    ax1.set_axisbelow(True)
    ax1.set_title(title_bp + " signatures - Adjusted Rand Index between clusters", fontsize=16, y=1.13)
    ax1.set_ylabel('Adjusted Rand Index', fontsize=14)
    ax1.set_xticks(list(range(1,17)))
    ax1.set_ylim(-0.25, 1.035)
    box_colors = ['orange', 'purple']
    num_boxes = len(data_bp)
    medians = np.empty(num_boxes)
    for i in range(num_boxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        box_coords = np.column_stack([boxX, boxY])
        ax1.add_patch(Polygon(box_coords, facecolor=box_colors[i % 2], alpha=0.6))

    ax1.set_xticklabels("")

    for i in list(zip([1,3,5,7,9,11,13,15], data_dict.keys())):
        ax1.text(i[0]-0.1, -0.33, i[1], fontsize=13)  

    for i in [2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5]:            
        ax1.axvline(i, color="grey", linestyle='--', lw=1)

    purple = Line2D([], [], color='orange', marker='s', markersize=8, label='Real vs Random', alpha=0.6)
    orang = Line2D([], [], color='purple', marker='s', markersize=8, label='Random vs Random', alpha=0.6)
    ax1.legend(handles=[purple, orang], loc=(0.32,1.03), fontsize=13, ncol=2)

    ax1.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax1.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1], fontsize=13)
    ax1.grid(color='grey', axis='y', linestyle='--', linewidth=0.5, alpha=0.7)
    plt.show()
    
# =============================================================================
#  Figure 6
# =============================================================================
import pandas as pd

# compute strictly standardized mean difference
def compute_ssmd(a,b): return (a.mean()-b.mean())/np.sqrt(a.var()+b.var())


ssmd_dict_moff = {k:[] for k in data_dict.keys()}
for i in data_moff: ssmd_dict_moff[i.name].append(i)

ssmd_dict_coll = {k:[] for k in data_dict.keys()}
for i in data_coll: ssmd_dict_coll[i.name].append(i)

ssmd_dict_bailey = {k:[] for k in data_dict.keys()}
for i in data_bailey: ssmd_dict_bailey[i.name].append(i)


ssmd_df = pd.DataFrame(index=["Moffitt", "Collisson", "Bailey"], columns=data_dict.keys(), dtype="float")

for k,v in ssmd_dict_moff.items():
    ssmd_df.loc["Moffitt"][k] = compute_ssmd(v[0],v[1])

for k,v in ssmd_dict_coll.items():
    ssmd_df.loc["Collisson"][k] = compute_ssmd(v[0],v[1])
    
for k,v in ssmd_dict_bailey.items():
    ssmd_df.loc["Bailey"][k] = compute_ssmd(v[0],v[1])
    
plt.figure(figsize=(5,2))
sns.heatmap(-np.array(ssmd_df), annot=True, xticklabels=ssmd_df.columns, yticklabels=ssmd_df.index, linewidths=1); 
plt.xlabel("Datasets"); plt.ylabel("Signatures"); plt.title("Strictly standardized mean difference")
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
variables needed:
data_dict = dictionary containing all the dataframes of the datasets we want to make predictions on 
            (e.g.: {"data X":data_X_df, "data Y":data_y_df, ...})

signature_dict = dictionary with list of genes of the signatures and subtypes,
                 e.g. {"Signature X": (list_genes_signature_X, subtypes_related_to_signature_X),...}
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon 
import matplotlib.lines as mlines
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy import stats

###############################################################################
# Figure 1 A
###############################################################################
boxplot_distributions = dict()

for signature_name, signature_tup in signature_dict.items():
    
    distr_boxplot = []
    for dataset_name in signature_dict.keys():
        print("Data:", dataset_name, "- Genes:", signature_name)
        
        dataset_df = data_dict[dataset_name].copy()
        study_df_filt = dataset_df[set(dataset_df.columns).intersection(signature_tup[0])]
        print(study_df_filt.shape)
        distr = np.std(study_df_filt)
        distr_boxplot.append(list(distr))
        
    boxplot_distributions[signature_name] = distr_boxplot

fig, ax1 = plt.subplots(figsize=(15, 6))
fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

bp = ax1.boxplot(boxplot_distributions["Moffitt et al. signatures"]+boxplot_distributions["Collisson et al. signatures"]+boxplot_distributions["Bailey et al. signatures"], notch=0, sym='+', vert=1, whis=1.5)

plt.setp(bp['boxes'], color='grey', linewidth=1.1)
plt.setp(bp['whiskers'], color='grey', linewidth=1.1)
plt.setp(bp['medians'], color='grey', linewidth=1.1)
plt.setp(bp['fliers'], color='gainsboro', marker='o', linewidth=1.1)
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

ax1.set_axisbelow(True)
ax1.set_title("Standard deviation of the three signatures in the three study cohorts", fontsize=16, y=1.03)
ax1.set_xlabel("")
ax1.set_ylabel('Standard deviation distribution', fontsize=14)
ax1.set_xticks(list(range(1,19)))
ax1.set_yticks(list(range(0,5)))

box_colors = ['darkorange', 'royalblue']
num_boxes = 18
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

ax1.set_xticklabels([])

ax1.text(2, -1.12, "Moffitt et al. signatures", fontsize=14)  
ax1.text(8, -1.12, "Collisson et al. signatures", fontsize=14) 
ax1.text(14, -1.12, "Bailey et al. signatures", fontsize=14) 
            
ax1.text(0.8, -0.6, "Moffitt et al.\n    data", fontsize=13)  
ax1.text(2.8, -0.6, "Collisson et al.\n    data", fontsize=13)  
ax1.text(4.8, -0.6, "Bailey et al.\n    data", fontsize=13)  
ax1.text(6.8, -0.6, "Moffitt et al.\n    data", fontsize=13)  
ax1.text(8.8, -0.6, "Collisson et al.\n    data", fontsize=13)  
ax1.text(10.8, -0.6, "Bailey et al.\n    data", fontsize=13)  
ax1.text(12.8, -0.6, "Moffitt et al.\n    data", fontsize=13)  
ax1.text(14.8, -0.6, "Collisson et al.\n    data", fontsize=13)  
ax1.text(16.8, -0.6, "Bailey et al.\n    data", fontsize=13)  

ax1.set_yticks([0, 1, 2, 3, 4])
ax1.set_yticklabels([0, 1, 2, 3, 4], fontsize=12)

for i in [6.5, 12.5]:            
    ax1.axvline(i, color="grey", lw=1.2)
    
for i in [2.5, 4.5, 8.5, 10.5, 14.5, 16.5]:            
    ax1.axvline(i, color="grey", lw=0.4, linestyle="--")

orange = mlines.Line2D([], [], color='darkorange', marker='s', markersize=11, label='Signature', alpha=0.6)
blue = mlines.Line2D([], [], color='royalblue', marker='s', markersize=11, label='All genes', alpha=0.6)
ax1.legend(handles=[orange, blue], title="Genes used", loc=2, fontsize=12, title_fontsize=12)
plt.show()

###############################################################################
# Figure 1 B
###############################################################################

all_pval_list = []
for s_genes in signature_dict.keys(): 
    # first take signatures 
    g = signature_dict[s_genes][0]
    print("\n signatures: ", s_genes)
    
    for s_data in signature_dict.keys():

        print("data: ", s_data)
        df = data_dict[s_data].copy()
        sub = signature_dict[s_data][1]
        
        df = df[set(df.columns).intersection(g)]     
        df.index  = sub
        
        if s_data=="Moffitt et al.":
            pval = []
            for i in df.columns:
                F, p = stats.ranksums(df.loc["Basal", i], df.loc["Classical", i])
                pval.append(p)

        elif s_data=="Collisson et al.":
            pval = []
            for i in df.columns:
                F, p = stats.kruskal(df.loc["Exocrine-like PDA", i], df.loc["QM-PDA", i], df.loc["Classical PDA", i])
                pval.append(p)

        elif s_data=="Bailey et al.":
            pval = []
            for i in df.columns:
                F, p = stats.kruskal(df.loc["ADEX", i], df.loc["Pancreatic Progenitor", i], df.loc["Immunogenic", i], df.loc["Squamous", i])
                pval.append(p)
        
        a, correct_pval, b, c = multipletests(pvals=pval, alpha=0.05, method="bonferroni")
        all_pval_list.append(-np.log10(correct_pval))   
         

all_boxes = all_pval_list
fig, ax1 = plt.subplots(figsize=(10, 5))
fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

bp = ax1.boxplot(all_boxes, notch=0, sym='+', vert=1, whis=1.5)

plt.setp(bp['boxes'], color='grey', linewidth=1)
plt.setp(bp['whiskers'], color='grey', linewidth=1)
plt.setp(bp['medians'], color='grey', linewidth=1)
plt.setp(bp['fliers'], color='gainsboro', marker='o', linewidth=1, markersize=4, alpha=0.5)
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

ax1.set_axisbelow(True)
ax1.set_title("Comparison of signatures expression level across subtypes", fontsize=16)
ax1.set_xlabel("")
ax1.set_ylabel('-log10(adjusted pvalue)', fontsize=14)
ax1.set_xticks(list(range(1,10)))
box_colors = ['gold', 'darkcyan', 'deeppink']
num_boxes = 9
medians = np.empty(num_boxes)
for i in range(num_boxes):
    box = bp['boxes'][i]
    boxX = []
    boxY = []
    for j in range(5):
        boxX.append(box.get_xdata()[j])
        boxY.append(box.get_ydata()[j])
    box_coords = np.column_stack([boxX, boxY])
    ax1.add_patch(Polygon(box_coords, facecolor=box_colors[i % 3], alpha=0.5))

ax1.set_xticklabels([])

gold = mlines.Line2D([], [], color='gold', marker='s', markersize=11, label='Moffitt et al.', alpha=0.5)
cyan = mlines.Line2D([], [], color='darkcyan', marker='s', markersize=11, label='Collisson et al.', alpha=0.5)
pink = mlines.Line2D([], [], color='deeppink', marker='s', markersize=11, label='Bailey et al.', alpha=0.5)
ax1.legend(handles=[gold, cyan, pink], title="Subtypes", loc=1, fontsize=12, title_fontsize=12)

ax1.text(1.6, -1.3, "Moffitt et al.", fontsize=13)  
ax1.text(4.6, -1.3, "Collisson et al.", fontsize=13) 
ax1.text(7.6, -1.3, "Bailey et al.", fontsize=13) 

ax1.axhline(-np.log10(0.05), color="darkred", linestyle='--', lw=1.2)
ax1.text(4.6, -2.1, "Signatures", fontsize=14) 
                     
for i in [3.5, 6.5]:            
    ax1.axvline(i, color="grey", linestyle='--', lw=0.5)

ax1.set_yticks([0, 2, 4, 6, 8])
ax1.set_yticklabels([0, 2, 4, 6, 8], fontsize=12)
ax1.grid(color='grey', axis='y', linestyle='--', linewidth=0.3, alpha=0.7)
plt.show()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
imput parameters:
n_rep            = number of repetitions
k                = number of k-fold cross-validation
RF_dict_accuracy = nested dictionary to store the cross-validation results; keys are the 9 models, 
                 values are dictionary containing 3 lists to collect cross-validation accuracy values computed using real 
                 or random signatures,
                 or shuffled subtype labels (e.g. {"Moffitt data - Moffitt genes":{"Real":[],"Random":[],"Shuffled":[]}, ...})
train_data_dict  = dictionary containing the information needed for training the model, in the format:
                 {"signature X":(discovery_dataset_used_for_signature_X, list_genes_in_signature_X, subtypes_related_to_signature_X),...}
signature_dict   = dictionary with list of genes of the signatures and subtypes,
                 e.g. {"Signature X": (list_genes_signature_X, subtypes_related_to_signature_X),...}
z_score          = True/False
"""

# =============================================================================
#  Run Random Forest k-fold cross-validation
# =============================================================================
def RF_CV_compare(n_rep, k, RF_dict_accuracy, train_data_dict, signature_dict , z_score):
                  
    import numpy as np
    from scipy.stats import zscore
    from sklearn.ensemble import RandomForestClassifier
    import random
    from sklearn.model_selection import cross_val_score
    
    for r in range(n_rep):
        for model_name in RF_dict_accuracy.keys():

            train_data_name = model_name.split("-")[0].split(" ")[0]

            #list of genes
            train_genes = signature_dict[model_name.split("-")[1].split(" ")[1]][0]
    
            df = train_data_dict[train_data_name][0]
            y  = train_data_dict[train_data_name][2]

            # z-score?
            if z_score: df = df.apply(zscore)

            # ACCURACY -----------------------------------------------------------------
            # train model using REAL signatures
            df_real = df[set(df.columns).intersection(train_genes)]
            rf_real = RandomForestClassifier(n_estimators=100)
            rf_real.fit(df_real, y) 
    
            cv  = cross_val_score(rf_real, df_real, y, cv=k, scoring="accuracy") 
            RF_dict_accuracy[model_name]["Real"].append(np.mean(cv))

            # we want to built the models only using the signatures from Moffitt, Collisson and Bailey
            # but want to keep Badea et al. for performance evaluation since it provides real subtype labels
            if train_data_name != "Badea et al.":
                # train model using Random signatures
                df_rand = df[random.sample(list(df.columns), df_real.shape[1])]

                rf_cl_rand = RandomForestClassifier(n_estimators=100)
                rf_cl_rand.fit(df_rand, y)
                cv  = cross_val_score(rf_cl_rand, df_rand, y, cv=k, scoring="accuracy") 
                RF_dict_accuracy[model_name]["Random"].append(np.mean(cv))

                # train model using Shuffled subtypes
                df_real = df[set(df.columns).intersection(train_genes)]
                shuffle_subtype = y.copy()
                random.shuffle(shuffle_subtype)            

                rf_cl_shuff = RandomForestClassifier(n_estimators=100)
                rf_cl_shuff.fit(df_real, shuffle_subtype)
                cv  = cross_val_score(rf_cl_shuff, df_real, shuffle_subtype, cv=k, scoring="accuracy") 
                RF_dict_accuracy[model_name]["Shuffled"].append(np.mean(cv))
            
    return(RF_dict_accuracy)

# =============================================================================
#  Figure 4
# =============================================================================
res = RF_CV_compare(n_rep=1000, k=5, RF_dict_accuracy, train_data_dict, signature_dict , z_score)

import matplotlib.pyplot as plt  
from matplotlib.patches import Polygon 
import matplotlib.lines as mlines

fig, ax1 = plt.subplots(figsize=(13, 6))
fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

bp = ax1.boxplot([res['Moffitt data - Moffitt genes']["Real"], res['Collisson data - Moffitt genes']["Real"], 
                  res['Bailey data - Moffitt genes']["Real"],  res['Badea data - Moffitt genes']["Real"]]
                +[res['Moffitt data - Collisson genes']["Real"], res['Collisson data - Collisson genes']["Real"], 
                  res['Bailey data - Collisson genes']["Real"],  res['Badea data - Collisson genes']["Real"]]
                +[res['Moffitt data - Bailey genes']["Real"], res['Collisson data - Bailey genes']["Real"], 
                  res['Bailey data - Bailey genes']["Real"],  res['Badea data - Bailey genes']["Real"]],
                 notch=0, sym='+', vert=1, whis=1.5)

plt.setp(bp['boxes'], color='grey', linewidth=1)
plt.setp(bp['whiskers'], color='grey', linewidth=1)
plt.setp(bp['medians'], color='grey', linewidth=1)
plt.setp(bp['fliers'], color='gainsboro', marker='o', linewidth=1)
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

ax1.set_axisbelow(True)
ax1.set_title("1000 repeated K-fold cross validation (mean accuracy)", fontsize=16)
ax1.set_xlabel("")
ax1.set_xticks(list(range(1,13)))
ax1.set_yticks(list(range(0,1)))
box_colors = ['darkorange', 'royalblue', 'mediumseagreen', 'purple']
num_boxes = 12
medians = np.empty(num_boxes)
for i in range(num_boxes):
    box = bp['boxes'][i]
    boxX = []
    boxY = []
    for j in range(5):
        boxX.append(box.get_xdata()[j])
        boxY.append(box.get_ydata()[j])
    box_coords = np.column_stack([boxX, boxY])
    ax1.add_patch(Polygon(box_coords, facecolor=box_colors[i % 4], alpha=0.6))

ax1.set_xticklabels([])

ax1.text(1.5, -0.07, "Moffitt et al. signatures", fontsize=14)  
ax1.text(5.5, -0.07, "Collisson et al. signatures", fontsize=14) 
ax1.text(9.5, -0.07, "Bailey et al. signatures", fontsize=14) 
                     
for i in [4.5, 8.5]:            
    ax1.axvline(i, color="grey", linestyle='--', lw=1)

import matplotlib.lines as mlines
orange = mlines.Line2D([], [], color='darkorange', marker='s', markersize=11, label='Moffitt et al.', alpha=0.6)
blue = mlines.Line2D([], [], color='royalblue', marker='s', markersize=11, label='Collisson et al.', alpha=0.6)
green = mlines.Line2D([], [], color='mediumseagreen', marker='s', markersize=11, label='Bailey et al.', alpha=0.6)
violet = mlines.Line2D([], [], color='purple', marker='s', markersize=11, label='Badea et al.', alpha=0.6)
ax1.legend(handles=[orange, blue, green, violet], title="Dataset", loc=3, fontsize=12, title_fontsize=12)

ax1.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
ax1.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8], fontsize=12)
ax1.grid(color='grey', axis='y', linestyle='--', linewidth=0.4, alpha=0.6)
plt.show()

# =============================================================================
#  Figure 5
# =============================================================================

moff_merge_box_moff_genes = [res['Moffitt data - '+"Moffitt genes"]["Real"], res['Moffitt data - '+"Moffitt genes"]["Random"], res['Moffitt data - '+"Moffitt genes"]["Shuffled"]]
moff_merge_box_coll_genes = [res['Moffitt data - '+"Collisson genes"]["Real"], res['Moffitt data - '+"Collisson genes"]["Random"], res['Moffitt data - '+"Collisson genes"]["Shuffled"]]
moff_merge_box_bail_genes = [res['Moffitt data - '+"Bailey genes"]["Real"], res['Moffitt data - '+"Bailey genes"]["Random"], res['Moffitt data - '+"Bailey genes"]["Shuffled"]]
merged_boxes_moffitt = moff_merge_box_moff_genes + moff_merge_box_coll_genes + moff_merge_box_bail_genes

coll_merge_box_moff_genes = [res['Collisson data - '+"Moffitt genes"]["Real"],   res['Collisson data - '+"Moffitt genes"]["Random"],   res['Collisson data - '+"Moffitt genes"]["Shuffled"]]
coll_merge_box_coll_genes = [res['Collisson data - '+"Collisson genes"]["Real"], res['Collisson data - '+"Collisson genes"]["Random"], res['Collisson data - '+"Collisson genes"]["Shuffled"]]
coll_merge_box_bail_genes = [res['Collisson data - '+"Bailey genes"]["Real"],    res['Collisson data - '+"Bailey genes"]["Random"],    res['Collisson data - '+"Bailey genes"]["Shuffled"]]
merged_boxes_collisson = coll_merge_box_moff_genes + coll_merge_box_coll_genes + coll_merge_box_bail_genes

bail_merge_box_moff_genes = [res['Bailey data - '+"Moffitt genes"]["Real"],   res['Bailey data - '+"Moffitt genes"]["Random"],   res['Bailey data - '+"Moffitt genes"]["Shuffled"]]
bail_merge_box_coll_genes = [res['Bailey data - '+"Collisson genes"]["Real"], res['Bailey data - '+"Collisson genes"]["Random"], res['Bailey data - '+"Collisson genes"]["Shuffled"]]
bail_merge_box_bail_genes = [res['Bailey data - '+"Bailey genes"]["Real"],    res['Bailey data - '+"Bailey genes"]["Random"],    res['Bailey data - '+"Bailey genes"]["Shuffled"]]
merged_boxes_bailey = bail_merge_box_moff_genes + bail_merge_box_coll_genes + bail_merge_box_bail_genes 

all_boxes = moff_merge_box_moff_genes + moff_merge_box_coll_genes + moff_merge_box_bail_genes + coll_merge_box_moff_genes + coll_merge_box_coll_genes + coll_merge_box_bail_genes + bail_merge_box_moff_genes + bail_merge_box_coll_genes + bail_merge_box_bail_genes 
fig, ax1 = plt.subplots(figsize=(15, 6))
fig.canvas.set_window_title('A Boxplot Example')
fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

bp = ax1.boxplot(all_boxes, notch=0, sym='+', vert=1, whis=1.5)

plt.setp(bp['boxes'], color='grey', linewidth=1)
plt.setp(bp['whiskers'], color='grey', linewidth=1)
plt.setp(bp['medians'], color='grey', linewidth=1)
plt.setp(bp['fliers'], color='gainsboro', marker='o', linewidth=1, markersize=4, alpha=0.5)
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

ax1.set_axisbelow(True)
ax1.set_title("Mean K-fold cross validation accuracy (repeated 1000 times)", fontsize=16, y=1.18)
ax1.set_xlabel("")
ax1.set_ylabel('K-fold CV (mean accuracy)', fontsize=14)
ax1.set_xticks(list(range(1,28)))
ax1.set_yticks(list(range(0,1)))
box_colors = ['orchid', 'mediumaquamarine', 'sandybrown']
num_boxes = 27
medians = np.empty(num_boxes)
for i in range(num_boxes):
    box = bp['boxes'][i]
    boxX = []
    boxY = []
    for j in range(5):
        boxX.append(box.get_xdata()[j])
        boxY.append(box.get_ydata()[j])
    box_coords = np.column_stack([boxX, boxY])
    ax1.add_patch(Polygon(box_coords, facecolor=box_colors[i % 3], alpha=0.6))
ax1.set_xticklabels([])

ax1.text(4, -0.07, "Moffitt et al.",    fontsize=14)  
ax1.text(13, -0.07, "Collisson et al.", fontsize=14) 
ax1.text(22.2, -0.07, "Bailey et al.",  fontsize=14) 

ax1.text(1,  1.04, "Moffitt et al.",   fontsize=14)  
ax1.text(4,  1.04, "Collisson et al.", fontsize=14) 
ax1.text(7,  1.04, "Bailey et al.",    fontsize=14) 
ax1.text(10, 1.04, "Moffitt et al.",   fontsize=14)  
ax1.text(13, 1.04, "Collisson et al.", fontsize=14) 
ax1.text(16, 1.04, "Bailey et al.",    fontsize=14) 
ax1.text(19, 1.04, "Moffitt et al.",   fontsize=14)  
ax1.text(22, 1.04, "Collisson et al.", fontsize=14) 
ax1.text(25, 1.04, "Bailey et al.",    fontsize=14) 

ax1.text(13.3, -0.15, "Subtypes", fontsize=15) 
ax1.text(13, 1.12, "Signatures:", fontsize=15) 
                     
for i in [3.5, 6.5, 12.5, 15.5, 18.5, 21.5, 24.5]:            
    ax1.axvline(i, color="grey", linestyle='--', lw=0.5)
ax1.axvline(9.5, color="grey", lw=1.4)
ax1.axvline(18.5, color="grey", lw=1.4)

import matplotlib.lines as mlines
orange = mlines.Line2D([], [], color='orchid', marker='s', markersize=8, label='Real signatures', alpha=0.6)
blue = mlines.Line2D([], [], color='mediumaquamarine', marker='s', markersize=8, label='Random signatures', alpha=0.6)
green = mlines.Line2D([], [], color='sandybrown', marker='s', markersize=8, label='Shuffled labels', alpha=0.6)
ax1.legend(handles=[orange, blue, green], title="Criteria", loc=3, fontsize=12, title_fontsize=12)

ax1.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
ax1.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1], fontsize=12)
ax1.grid(color='grey', axis='y', linestyle='--', linewidth=0.5, alpha=0.7)
plt.show()
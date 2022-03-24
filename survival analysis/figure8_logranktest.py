#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
variables needed:
surv = dictionary with names of data with survival info as keys and survival data as values (df with time and event column)
data = dictionary with names of data with survival info as keys and samples name as values
"""

# =============================================================================
#  Compute log-rank test
# =============================================================================

import json
with open('.../predicted_labels.json') as json_data:
    pred = json.load(json_data)

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test, pairwise_logrank_test
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

surv_pval = {k:{k:[] for k in surv.keys()} for k in pred.keys()}
surv_pval_moff = {'Moffitt signature':[], 'Collisson signature':[], 'Bailey signature':[]}
surv_pval_coll = {'Moffitt signature':[], 'Collisson signature':[], 'Bailey signature':[]}
surv_pval_bail = {'Moffitt signature':[], 'Collisson signature':[], 'Bailey signature':[]}
    
for model in pred.keys():
    for surv_name, surv_df in surv.items():
        pred_list = pred[model][surv_name] # for that dataset i get the list of prediction for a specific model
        idx = data[surv_name]
        pred_df = pd.DataFrame(pred_list, index = idx, columns=["predicted"])
        
        index_common = list(set(idx).intersection(surv_df.index))  # i extract the samples index for which we have survival data
        pred_df = pred_df[pred_df.index.isin(index_common)]
        surv_df = surv_df[surv_df.index.isin(index_common)]

        combined = surv_df.join(pred_df)
        combined = combined.astype({'time': 'float'})
        combined = combined.astype({'event': 'int'})
        subt = list(set(combined["predicted"])) # list of the unique names of the subtypes
        
        if len(subt)==2:
            g0 = combined[combined["predicted"] == subt[0]]
            g1 = combined[combined["predicted"] == subt[1]]
            pval = logrank_test(g0.iloc[:,0], g1.iloc[:,0], g0.iloc[:,1], g1.iloc[:,1]).p_value
            
            # save for heatmap:
            surv_pval_moff[model.split(" - ")[1].split(" ")[0]+" signature"].append(pval)

        elif len(subt)==3:

            # Pairwise log-rank test
            pair_test = pairwise_logrank_test(combined["time"], combined['predicted'], combined['event']).summary
            pair_test_coll = list(map(lambda x: x[0]+" - "+x[1], pair_test.index))*5
            
            collect_pval = []
            for i in pair_test.index:
                # save for heatmap:
                surv_pval_coll[model.split(" - ")[1].split(" ")[0]+" signature"].append(pair_test.loc[i,"p"])
        
            collect_pval_df = pd.DataFrame(collect_pval)
            collect_pval_df.index = [i.split(":")[0] for i in collect_pval_df[0]]
            collect_pval_df["p"] = [": "+i.split(":")[1][1:] for i in collect_pval_df[0]]
            surv_pval[model][surv_name] = collect_pval_df['p']
            
        elif len(subt)==4:

            # Pairwise log-rank test
            pair_test = pairwise_logrank_test(combined["time"], combined['predicted'], combined['event']).summary
            pair_test_bail = list(map(lambda x: x[0]+" - "+x[1], pair_test.index))*5
            
            collect_pval = []
            for i in pair_test.index:
                
                # save for heatmap:
                surv_pval_bail[model.split(" - ")[1].split(" ")[0]+" signature"].append(pair_test.loc[i,"p"])
                    
            collect_pval_df = pd.DataFrame(collect_pval)
            collect_pval_df.index = [i.split(":")[0] for i in collect_pval_df[0]]
            collect_pval_df["p"] = [": "+i.split(":")[1][1:] for i in collect_pval_df[0]]
            surv_pval[model][surv_name] = collect_pval_df['p']


# =============================================================================
#  Generate heatmap
# =============================================================================

surv_pval_moff = pd.DataFrame(surv_pval_moff, index = surv.keys())
surv_pval_coll = pd.DataFrame(surv_pval_coll, index = pair_test_coll)
surv_pval_coll.index = [i.replace(" PDA","").replace("-PDA","") for i in surv_pval_coll.index]
surv_pval_bail = pd.DataFrame(surv_pval_bail, index = pair_test_bail)
surv_pval_bail.index = [i.replace("Pancreatic Progenitor","PP") for i in surv_pval_bail.index]
import matplotlib.lines as mlines

# heatmap 1
c = sns.clustermap(abs(-np.log10(surv_pval_moff.astype(float))), col_colors = ["orange", "peru", "brown"], row_colors = ["dodgerblue","violet","lightgrey","tomato","pink"], 
                   row_cluster=False, col_cluster=False, xticklabels=False, cmap='summer', figsize=(4,2), annot=True, annot_kws={"size": 11})
c.ax_heatmap.set_yticklabels(c.ax_heatmap.get_ymajorticklabels(), fontsize = 13)
c1 = mlines.Line2D([], [], color='orange', marker='s', markersize=11, label='Moffitt et al.', alpha=0.9)
c2 = mlines.Line2D([], [], color='peru', marker='s', markersize=11, label='Collisson et al.', alpha=0.9)
c3 = mlines.Line2D([], [], color='brown', marker='s', markersize=11, label='Bailey et al.', alpha=0.9)
plt.legend(handles=[c1, c2, c3], title="Signature", loc=2)

# heatmap 2
c = sns.clustermap(abs(-np.log10(surv_pval_bail.astype(float))), row_cluster=False, col_cluster=False, col_colors=["orange", "peru", "brown"], xticklabels=False, 
                   row_colors = ['dodgerblue']*6+['violet']*6+['lightgrey']*6+['tomato']*6+['pink']*6, cmap='summer', figsize=(7,10), annot=True, annot_kws={"size": 11})
c.ax_heatmap.set_yticklabels(c.ax_heatmap.get_ymajorticklabels(), fontsize = 13)
c4 = mlines.Line2D([], [], color='dodgerblue', marker='s', markersize=11, label='Moffitt et al.', alpha=0.9)
c5 = mlines.Line2D([], [], color='violet', marker='s', markersize=11, label='Collisson et al.', alpha=0.9)
c6 = mlines.Line2D([], [], color='lightgrey', marker='s', markersize=11, label='Bailey et al.', alpha=0.9)
c7 = mlines.Line2D([], [], color='tomato', marker='s', markersize=11, label='TCGA-PAAD', alpha=0.9)
c8 = mlines.Line2D([], [], color='pink', marker='s', markersize=11, label='Yang et al.', alpha=0.9)
plt.legend(handles=[c4,c5,c6,c7,c8], title="Cohort", loc=2)

# heatmap 3
c = sns.clustermap(abs(-np.log10(surv_pval_coll.astype(float))), row_cluster=False, col_cluster=False, col_colors =["orange", "peru", "brown"], xticklabels=False, 
                   row_colors = ['dodgerblue']*3+['orange']*3+['lightgrey']*3+['tomato']*3+['pink']*3, cmap='summer', figsize=(5,7), annot=True, annot_kws={"size": 11})
c.ax_heatmap.set_yticklabels(c.ax_heatmap.get_ymajorticklabels(), fontsize = 13)


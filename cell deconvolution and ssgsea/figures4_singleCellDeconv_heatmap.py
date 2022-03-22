#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
deconv_data = dictionary with keys as dataset name and values as sincle-cell deconvolution output matrices
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json

with open('.../predicted_labels.json') as json_data:
    pred = json.load(json_data)
    
for df_test_name in deconv_data.keys():

    df_toplot = deconv_data[df_test_name].copy()
    classes1 = pred["Moffitt data - Moffitt genes"][df_test_name]
    classes2 = pred["Collisson data - Moffitt genes"][df_test_name]
    classes3 = pred["Bailey data - Moffitt genes"][df_test_name]  
    classes4 = pred["Moffitt data - Collisson genes"][df_test_name]
    classes5 = pred["Collisson data - Collisson genes"][df_test_name]
    classes6 = pred["Bailey data - Collisson genes"][df_test_name]  
    classes7 = pred["Moffitt data - Bailey genes"][df_test_name]
    classes8 = pred["Collisson data - Bailey genes"][df_test_name]
    classes9 = pred["Bailey data - Bailey genes"][df_test_name]  

    colors = {"Basal":"orange", "Classical":"lightgrey", "Classical PDA":"mediumvioletred", "Exocrine-like PDA":"pink", "QM-PDA":"dodgerblue",
                   "ADEX":"tomato", "Pancreatic Progenitor":"skyblue", "Squamous":"gold", "Immunogenic":"grey"}

    title = "Cell deconvolution of "+df_test_name

    df_toplot = df_toplot.T.reset_index().sort_values("index").set_index("index").T
    df_toplot2 = df_toplot.T.copy()

    df_toplot2["Moffitt subtypes - Moffitt signatures"] = classes1;   df_toplot2["Moffitt subtypes - Moffitt signatures"] = df_toplot2["Moffitt subtypes - Moffitt signatures"].map(colors)
    df_toplot2["Collisson subtypes - Moffitt signatures"] = classes2; df_toplot2["Collisson subtypes - Moffitt signatures"] = df_toplot2["Collisson subtypes - Moffitt signatures"].map(colors)
    df_toplot2["Bailey subtypes - Moffitt signatures"] = classes3;    df_toplot2["Bailey subtypes - Moffitt signatures"] = df_toplot2["Bailey subtypes - Moffitt signatures"].map(colors)
    df_toplot2[" "] = ["white"]*df_toplot2.shape[0]
    df_toplot2["Moffitt subtypes - Collisson signatures"] = classes4;   df_toplot2["Moffitt subtypes - Collisson signatures"] = df_toplot2["Moffitt subtypes - Collisson signatures"].map(colors)
    df_toplot2["Collisson subtypes - Collisson signatures"] = classes5; df_toplot2["Collisson subtypes - Collisson signatures"] = df_toplot2["Collisson subtypes - Collisson signatures"].map(colors)
    df_toplot2["Bailey subtypes - Collisson signatures"] = classes6;    df_toplot2["Bailey subtypes - Collisson signatures"] = df_toplot2["Bailey subtypes - Collisson signatures"].map(colors)
    df_toplot2["  "] = ["white"]*df_toplot2.shape[0]
    df_toplot2["Moffitt subtypes - Bailey signatures"] = classes7;   df_toplot2["Moffitt subtypes - Bailey signatures"] = df_toplot2["Moffitt subtypes - Bailey signatures"].map(colors)
    df_toplot2["Collisson subtypes - Bailey signatures"] = classes8; df_toplot2["Collisson subtypes - Bailey signatures"] = df_toplot2["Collisson subtypes - Bailey signatures"].map(colors)
    df_toplot2["Bailey subtypes - Bailey signatures"] = classes9;    df_toplot2["Bailey subtypes - Bailey signatures"] = df_toplot2["Bailey subtypes - Bailey signatures"].map(colors)
    df_toplot2["   "] = ["white"]*df_toplot2.shape[0]

    df_toplot2 = df_toplot2.T

    df_toplot2.columns.name="Samples"
    df_toplot3 = pd.DataFrame(np.matrix(df_toplot2.iloc[:-12,:]).astype('float'), index=df_toplot2.index[:-12], columns=df_toplot2.columns)
    
    col_colors_df = pd.DataFrame(np.matrix(df_toplot2.T[df_toplot2.iloc[-12:,:].index]).T, 
                                 index=df_toplot2.iloc[-12:,:].index,columns=df_toplot2.columns)

    g=sns.clustermap(df_toplot3, cmap="Spectral", figsize=(9,6), row_cluster=True, col_cluster=True, col_colors=col_colors_df.T, vmin=0, vmax=1, yticklabels=True, xticklabels=False)

    sns.set(font_scale=0.7)
    g.fig.suptitle(title, y=1.07, x=.5, fontsize=11)
    plt.xticks(fontsize=12)
    box = g.ax_heatmap.get_position()
    # add legend
    for label in set(classes1):
        g.ax_row_dendrogram.bar(0, 0, color=colors[label], label=label, linewidth=0)
    g.ax_row_dendrogram.legend(loc="center left", ncol=1, facecolor="white", fontsize=7)
    for label in set(classes2):
        g.ax_row_dendrogram.bar(0, 0, color=colors[label], label=label, linewidth=0)
    g.ax_row_dendrogram.legend(loc="center left", ncol=1, facecolor="white", fontsize=7)
    for label in set(classes3):
        g.ax_row_dendrogram.bar(0, 0, color=colors[label], label=label, linewidth=0)
    g.ax_row_dendrogram.legend(loc="center left", ncol=1, facecolor="white", fontsize=7, bbox_to_anchor=(-0.3,0.5))

    plt.show()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gseapy as gp
# apply ssGSEA to every dataset
# gene_sets used: KEGG_2016, Reactome_2016, GO_Cellular_Component_2018, GO_Molecular_Function_2018, GO_Biological_Process_2018
ss_df_1_... = gp.ssgsea(data=df_1, gene_sets=...).res2d 

"""
Create summary of the enrichment terms:
collect all ssGSEA output matrices of every dataset and for every gene_set in a dictionary:
data_ssgsea = {"kegg":{data1:ss_df_1_kegg, data_2:ss_df_2_kegg,...}, "reactome":{data1:ss_df_1_reac, data2:ss_df_2_reac,...}}
"""

import json
with open('.../predicted_labels.json') as json_data:
    pred = json.load(json_data)


def summary_EA():

    import pandas as pd
    import scipy.stats import zscore
    from statsmodels.sandbox.stats.multicomp import multipletests
    import scipy.stats as stats
    
    summary_dict = {k:{} for k in data_ssgsea.keys()}

    for enrich_name in data_ssgsea.keys():

        df_ea_moffitt = pd.DataFrame(0, index=list(data_ssgsea["KEGG"].keys()), columns=[])
        df_ea_collisson = pd.DataFrame(0, index=list(data_ssgsea["KEGG"].keys()), columns=[])
        df_ea_bailey = pd.DataFrame(0, index=list(data_ssgsea["KEGG"].keys()), columns=[])

        for data_name, ea_table in data_ssgsea[enrich_name].items():
            #e.g. "TCGA", TCGA ssgsea matrix        
            if enrich_name == "MF": ea_table = ea_table.groupby(ea_table.index).first()

            ea_table = ea_table.drop_duplicates()
            
            for model in pred.keys():
                
                    pred_list = pred[model][data_name]
                    # if there's any subtype with only one sample -> skip
                    if min([pred_list.count(x) for x in pred_list]) == 1: 
                        pass
                    else:                            
                        ea_table.columns = pred_list               
    
                        # sort columns by predicted class
                        df_toplot = ea_table.T.reset_index().sort_values("index").set_index("index").T.copy()
                                
                        # compute t-test and ANOVA between groups of subtypes
                        if model.split(" ")[0]=="Moffitt":
                            pval = []
                            for i in df_toplot.index:
                                F, p = stats.ttest_ind(df_toplot.loc[i,"Basal"], df_toplot.loc[i,"Classical"])
                                pval.append(p)
    
                        elif model.split(" ")[0]=="Collisson":
                            pval = []
                            for i in df_toplot.index:
                                F, p = stats.f_oneway(df_toplot.loc[i,"Exocrine-like PDA"], df_toplot.loc[i,"QM-PDA"], df_toplot.loc[i,"Classical PDA"])
                                pval.append(p)
    
                        elif model.split(" ")[0]=="Bailey":
                            pval = []
                            for i in df_toplot.index:
                                F, p = stats.f_oneway(df_toplot.loc[i,"ADEX"], df_toplot.loc[i,"Pancreatic Progenitor"], df_toplot.loc[i,"Immunogenic"], df_toplot.loc[i,"Squamous"])
                                pval.append(p)
    
                        df_toplot["removenan"] = pval
                        df_toplot = df_toplot.dropna()
                        newpval = list(df_toplot.removenan)
    
                        df_toplot = df_toplot.drop(["removenan"], axis=1)
    
                        a, corr_pval, b, c = multipletests(pvals=newpval, alpha=0.05, method="bonferroni") 
    
                        df_toplot["test"] = corr_pval
                        df_toplot = df_toplot.sort_values("test")
                        df_toplot = df_toplot[df_toplot["test"]<0.05]
                        df_toplot = df_toplot.drop(["test"], axis=1)
    
                        # for each term, see where is up or downregulated
                        df_toplot = df_toplot.T.apply(zscore).T
                        # now the dataset has subtypes as column names and the terms are filtered by adjusted pvalue obtained from t-test/anova
                        for i in df_toplot.index:
                            if model.split(" ")[0]=="Moffitt":
    
                                # for each term i look at the mean z-score in every subtype if it's >0 or <0
                                s1 = ("Basal", df_toplot.loc[i,"Basal"].mean(), "+" if df_toplot.loc[i,"Basal"].mean()>0 else "-")
                                s2 = ("Classical", df_toplot.loc[i,"Classical"].mean(), "+" if df_toplot.loc[i,"Classical"].mean()>0 else "-")
    
                                # i save only the subtypes with the highest mean zscore for this term
                                s3 = i+ " _ " + max([s1,s2], key=lambda item:abs(item[1]))[0]+ " ("+max([s1,s2], key=lambda item:abs(item[1]))[2]+")"
                                if s3 in df_ea_moffitt.columns: 
                                    df_ea_moffitt.loc[data_name, s3] += 1
                                else: 
                                    df_ea_moffitt[s3] = 0
                                    df_ea_moffitt.loc[data_name, s3] += 1
                                    
    
                            elif model.split(" ")[0]=="Collisson":
                                s1 = ("Exocrine-like PDA", df_toplot.loc[i,"Exocrine-like PDA"].mean(), "+" if df_toplot.loc[i,"Exocrine-like PDA"].mean()>0 else "-")
                                s2 = ("QM-PDA", df_toplot.loc[i,"QM-PDA"].mean(), "+" if df_toplot.loc[i,"QM-PDA"].mean()>0 else "-")
                                s3 = ("Classical PDA", df_toplot.loc[i,"Classical PDA"].mean(), "+" if df_toplot.loc[i,"Classical PDA"].mean()>0 else "-")
                                
                                s4 = i+ " _ " + max([s1,s2,s3], key=lambda item:abs(item[1]))[0]+ " ("+max([s1,s2,s3], key=lambda item:abs(item[1]))[2]+")"
                                if s4 in df_ea_collisson.columns: 
                                    df_ea_collisson.loc[data_name, s4] += 1
                                else: 
                                    df_ea_collisson[s4] = 0
                                    df_ea_collisson.loc[data_name, s4] += 1
    
                            elif model.split(" ")[0]=="Bailey":
                                s1 = ("ADEX", df_toplot.loc[i,"ADEX"].mean(), "+" if df_toplot.loc[i,"ADEX"].mean()>0 else "-")
                                s2 = ("Pancreatic Progenitor", df_toplot.loc[i,"Pancreatic Progenitor"].mean(), "+" if df_toplot.loc[i,"Pancreatic Progenitor"].mean()>0 else "-")
                                s3 = ("Immunogenic", df_toplot.loc[i,"Immunogenic"].mean(), "+" if df_toplot.loc[i,"Immunogenic"].mean()>0 else "-")
                                s4 = ("Squamous", df_toplot.loc[i,"Squamous"].mean(), "+" if df_toplot.loc[i,"Squamous"].mean()>0 else "-")
    
                                s5 = i+ " _ " + max([s1,s2,s3,s4], key=lambda item:abs(item[1]))[0]+ " ("+max([s1,s2,s3,s4], key=lambda item:abs(item[1]))[2]+")"
                                if s5 in df_ea_bailey.columns: 
                                    df_ea_bailey.loc[data_name, s5] += 1
                                else: 
                                    df_ea_bailey[s5] = 0
                                    df_ea_bailey.loc[data_name, s5] += 1

        df_ea_moffitt = df_ea_moffitt.T
        df_ea_collisson = df_ea_collisson.T
        df_ea_bailey = df_ea_bailey.T
        # i want to sort the terms according to the amount of datasets in which they are enriched (if==0 means not enriched, so i count the nonzero elements)
        df_ea_moffitt["nonzero"] = df_ea_moffitt[df_ea_moffitt!=0].count(axis=1)
        df_ea_collisson["nonzero"] = df_ea_collisson[df_ea_collisson!=0].count(axis=1)
        df_ea_bailey["nonzero"] = df_ea_bailey[df_ea_bailey!=0].count(axis=1)
        
        summary_dict[enrich_name]["Moffitt"] = df_ea_moffitt.sort_values(by="nonzero", ascending=False).iloc[:,:-1]
        summary_dict[enrich_name]["Collisson"] = df_ea_collisson.sort_values(by="nonzero", ascending=False).iloc[:,:-1]
        summary_dict[enrich_name]["Bailey"] = df_ea_bailey.sort_values(by="nonzero", ascending=False).iloc[:,:-1]
    return(summary_dict)
    
# =============================================================================
# Run fuction and generate plot
# =============================================================================

summary_EA_res = summary_EA()

def plot_summary(ea_type):
    
    import matplotlib.pyplot as plt
    
    plt.figure()
    summary_count_mat = pd.concat([summary_EA_res[ea_type]["Moffitt"], summary_EA_res[ea_type]["Collisson"], summary_EA_res[ea_type]["Bailey"]])
    # to group same terms
    summary_count_mat["term"] = [i.split(" _")[0] for i in summary_count_mat.index]
    summary_count_mat["type"] = [i.split(" _")[1] for i in summary_count_mat.index]
    summary_count_mat = summary_count_mat.groupby(['term'],as_index=False).agg({'type': '~'.join, 'Moffitt': 'sum', 'Collisson':'sum', 'Bailey':'sum', 
                                             'TCGA-PAAD':'sum', 'Badea et al.':'sum','Yang et al':'sum', 'ICGC-Array':'sum', 'Sandhu et al':'sum'})
    summary_count_mat.index = summary_count_mat["term"]+" _ "+summary_count_mat["type"]
    summary_count_mat = summary_count_mat.drop(["term", "type"], axis=1)

    # sort for amount of enriched times
    summary_count_mat["sum"] = summary_count_mat.sum(axis=1)
    summary_count_mat = summary_count_mat.sort_values(by="sum", ascending=False).iloc[:,:-1]

    mat_toplot = summary_count_mat.copy()
    mat_toplot[["Basal", "Classical", "Classical PDA", "QM-PDA", "Exocrine-like PDA", "Squamous", "Pancreatic Progenitor", "ADEX", "Immunogenic"]] = 0
    mat_toplot.index = [i.split(" _")[0] for i in mat_toplot.index]

    terms_to_drop = []
    for term in summary_count_mat.index:
        
        # i want to see every subtype if it has different "direction" in different model, in the same term.. if yes -> drop
        counting = []
        for j in term.split("_ ")[1].split("~"):
            subtype_ = j[1:].split("(")[0][:-1]
            if j[-2] == "-": direction = "cornflowerblue"
            elif j[-2] == "+": direction = "darkorange"
            mat_toplot.loc[term.split(" _")[0], subtype_] = direction
            counting.append(subtype_)
        for i in counting:
            if counting.count(i)>1: terms_to_drop.append(term.split("_ ")[0][:-1])
    terms_to_drop = list(set(terms_to_drop))
    print(len(terms_to_drop))

    col_df=pd.DataFrame([mat_toplot["Basal"], mat_toplot["Classical"], mat_toplot["Classical PDA"], mat_toplot["QM-PDA"], mat_toplot["Exocrine-like PDA"], 
                                    mat_toplot["Squamous"], mat_toplot["Pancreatic Progenitor"], mat_toplot["ADEX"], mat_toplot["Immunogenic"]]).T
    col_df["Classical"] = ["cornflowerblue" if i=="darkorange" else "darkorange" if i=="cornflowerblue" else i for i in mat_toplot["Basal"]]
            
    col_df = col_df.drop(terms_to_drop, axis=0)
    mat_toplot = mat_toplot.drop(terms_to_drop, axis=0)
    
    c=sns.clustermap(mat_toplot.iloc[:30,:8], row_colors=col_df.replace(0,"white"), linewidths=0.6,
                   cmap="summer", row_cluster=False, col_cluster=False, yticklabels=True, annot=True, figsize=(8,7), annot_kws={"size": 6})
    # row_colors 
    ax_row_colors = c.ax_row_colors; box = ax_row_colors.get_position()
    box_heatmap = c.ax_heatmap.get_position()
    # adjust coordinates according to the gene_set
    ax_row_colors.set_position([box_heatmap.x0-0.13, box.y0, box.width*0.65, box.height])
    c.ax_row_colors.tick_params(pad=1, labelsize=6.5)
    c.ax_row_colors.set_xticklabels(["Basal", "Classical", "Classical PDA", "QM-PDA", "Exocrine-like PDA", "Squamous", "Pancreatic Progenitor", "ADEX", "Immunogenic"])

    # heatmap, ticks and labels
    c.ax_heatmap.tick_params(labelsize=7, length=2, width=0.8)
    c.ax_heatmap.set_yticklabels(c.ax_heatmap.get_ymajorticklabels(), fontsize = 6.9) 
    c.ax_heatmap.set_xticklabels(['Moffitt et al.','Collisson et al.','Bailey et al.','TCGA-PAAD','Badea et al.','Yang et al.','ICGC-Array','Sandhu et al.'], 
                                 fontsize = 7, rotation=90)
    
    # title
    if ea_type == 'CC': ea_title = "Cellular components"
    elif ea_type == "MF": ea_title = "Molecular functions"
    elif ea_type == "BP": ea_title = "Biological process"
    elif ea_type in ["KEGG", "Reactome"]: ea_title = ea_type + " pathways"
    
    # adjust coordinates according to the gene_set
    c.fig.suptitle("ssGSEA - "+ea_title, x=0.38, y=0.91, fontsize=9)
    
    # adjust colorbar according to the gene_set
    c.cax.set_position([.08, .4, .01, .3]) 
    c.cax.tick_params(labelsize=6.5, length=2, width=0.8)
    
    # legend
    label_col1 = {"down": "cornflowerblue", "up": "darkorange"}
    labels1 = ["down" if i[-2]=="-" else "up" for i in summary_count_mat.iloc[:30,:].index]
    for label in set(labels1):
        c.ax_col_dendrogram.bar(0, 0, color=label_col1[label], label=label, linewidth=0)
    c.ax_col_dendrogram.legend(loc="lower center", ncol=2, facecolor="white", fontsize=8)   
    
    plt.show()
    
    return mat_toplot.iloc[:,:8]

mat = plot_summary("KEGG")

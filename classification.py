#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Function to perform Random Forest classification on each dataset using in turn the nine models built on
the three signatures using their discovery dataset

input parameters:
train_model    = name of the model you want to built, which includes the signature and the dataset name to use for training
test_data_name = name of the dataset whose sample's subtypes we want to predict
z_score        = True/False

variables needed:
data_dict       = dictionary containing all the dataframes of the datasets we want to make predictions on 
                  (e.g.: {"data X":data_X_df, "data Y":data_y_df, ...})
train_data_dict = dictionary containing the information needed for training the model, in the format:
                  {"signature X":(discovery_dataset_used_for_signature_X, list_genes_in_signature_X, subtypes_related_to_signature_X),...}
"""

def classification(data_dict, train_data_dict, train_model, test_data_name, z_score):

    from scipy.stats import zscore
    from sklearn.ensemble import RandomForestClassifier

    train_model_name = train_model.split("-")[0].split(" ")[0]

    X_train = train_data_dict[train_model_name][0].copy()
    y_train = train_data_dict[train_model_name][2].copy()
    
    genes   = train_data_dict[train_model.split("-")[1].split(" ")[1]][1]
    X_test  = data_dict[test_data_name].copy()

    # if i want z-score normalization
    if z_score: 
        X_train = X_train.apply(zscore)
        X_test  = X_test.apply(zscore)
        
    print("\nTrain on:  ",train_model+" - Test on ", test_data_name," data")
    
    # filter genes of the current classifier     
    common_genes = list(set.intersection(set(X_train.columns), set(X_test.columns), set(genes)))
    X_train      = X_train.loc[:,common_genes]
    X_test       = X_test.loc[:,common_genes]

    # train RF  ----------------------------------------------------------------------------------   
    clf=RandomForestClassifier(random_state=0)
    clf.fit(X_train,y_train)
    y_pred=clf.predict(X_test)

    return (y_pred.tolist())

#######################################################################################
# Run the function using the nine models
RF_models = ["Moffitt data - Moffitt genes", "Moffitt data - Collisson genes", "Moffitt data - Bailey genes", 
             "Collisson data - Moffitt genes", "Collisson data - Collisson genes", "Collisson data - Bailey genes",
             "Bailey data - Moffitt genes", "Bailey data - Collisson genes", "Bailey data - Bailey genes"]

res_dict = {model:{} for model in RF_models}
"""
Dictionary with predictions: {"classification_train_model_1":{"test_data_1":predicted_labels, "test_data_2":predicted_labels,...}, ...}
"""
for model in RF_models:
    for test_data_name in data_dict.keys():       
        
        pred = classification(model, test_data_name, z_score=True, batchCorrected=True, heatmap=True, pca=False)
        res_dict[model][test_data_name] = pred
 
# =============================================================================
#  Save predictions     
# =============================================================================
import json
with open('.../predicted_labels.json', 'w') as fp: json.dump(res_dict, fp)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
variable needed:
df_with_subtypes = dictionary of the dataset whose subtype information is available, where keys are the data name and 
                   values are tuples (df, subtypes), e.g.: {"data X": (discovery_dataset_used_for_signature_X, subtypes_related_to_signature_X),...}   
"""

# =============================================================================
#  Import predictions
# =============================================================================
import json
with open('.../predicted_labels.json') as json_data:
    pred = json.load(json_data)
    
"""
Dictionary with predictions: {"classification_train_model_1":{"test_data_1":predicted_labels, "test_data_2":predicted_labels,...}, ...}
"""

import matplotlib.pyplot as plt 
import pandas as pd

for model in pred.keys():
    for cohort in pred[model].keys():
        if cohort in ["Moffitt", "Collisson", "Bailey"]: 
            if cohort+" et al." in df_with_subtypes.keys():
                real_pred = pd.DataFrame({'Real': list(df_with_subtypes[cohort+" et al."][1]), 'Predicted':pred[model][cohort]})
                real_pred.groupby(['Predicted','Real']).size().unstack().plot.bar(stacked=True)
                plt.title("Model: "+model+"\n Cohort: "+cohort+" et al.", fontsize=16)
                plt.xticks(rotation=45, fontsize=14)
                plt.legend(fontsize=12)
                plt.yticks(fontsize=13)
                plt.xlabel("Predicted", fontsize=14)
                plt.show()
        else:
            if cohort in df_with_subtypes.keys():
                real_pred = pd.DataFrame({'Real': list(df_with_subtypes[cohort][1]), 'Predicted':pred[model][cohort]}).sort_values(by='Predicted')
                real_pred.groupby(['Predicted','Real']).size().unstack().plot.bar(stacked=True)
                plt.legend(fontsize=12)
                plt.title("Model: "+model+"\n Cohort: "+cohort, fontsize=16)
                plt.xticks(rotation=45, fontsize=14)
                plt.yticks(fontsize=13)
                plt.xlabel("Predicted", fontsize=14)
                plt.show()

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib

" ------------------------- Load Data -------------------------- "

files_list = os.listdir(os.path.join('..', '..', 'Output', 'Transfer_Learning_Diseases', 'Prediction_Slim'))
print('files_list', len(files_list), files_list)

predictions_list = []

for file_name in files_list:

    short_name = file_name.replace('_Specific_Probabilities.csv', '')
    disease_name = file_name.split('_')[0]
    ps_num = file_name.split('_')[1]
    genes = file_name.split('_')[2]
    print("@", disease_name)
    path = os.path.join('..', '..', 'Output', 'Transfer_Learning_Diseases', 'Prediction_Slim', file_name)
    Disease_Cell_Probabilities = pd.read_csv(path)
    Cell_Median_Probabilities = pd.DataFrame(Disease_Cell_Probabilities.groupby('Cell_Type', as_index=True)['scRNA_ML_Probability'].median())#scRNA_ML_Probability, preferential_expression, expression

    Cell_Median_Probabilities = Cell_Median_Probabilities.T
    Cell_Median_Probabilities['Disease'] = disease_name

    predictions_list.append(Cell_Median_Probabilities)

All_Disease_Cell_Predictions_Median = pd.concat(predictions_list)
print(All_Disease_Cell_Predictions_Median)
All_Disease_Cell_Predictions_Median.set_index('Disease', inplace=True)
relevant_cols = [c for c in list(All_Disease_Cell_Predictions_Median) if c != 'index']
print(All_Disease_Cell_Predictions_Median[relevant_cols])

matplotlib.rcParams.update({'font.size': 0.1})
plt.figure(0, figsize=(1000, 500))
ax = sns.clustermap(All_Disease_Cell_Predictions_Median[relevant_cols], z_score=0, yticklabels=True, cmap="coolwarm", vmin=-4, vmax=4)#Blues"
locs, labels = plt.yticks()            # Get locations and labels
print(locs, labels )
path = os.path.join('..', '..', 'Output', 'Transfer_Learning_Diseases', 'Summery', 'All_Disease_Cell_Heatmap_Prediction_Blue_Slim.pdf')
plt.savefig(path)

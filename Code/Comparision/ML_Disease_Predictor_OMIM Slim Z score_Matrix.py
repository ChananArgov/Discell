import os

import numpy as np
import pandas as pd
import networkx as nx

import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import classification_report
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
import matplotlib.pyplot as plt
import sklearn as skl
from sklearn.metrics import average_precision_score
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from xgboost import XGBClassifier
from statsmodels.stats.multitest import multipletests
from sklearn import metrics



# from scipy import interp

pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 10)

" ------------------------- Load Data -------------------------- "

path= os.path.join('..', '..', 'Output', 'Full_scRNA_Datasets', 'Transferred_scRNA_Fisher_Full_Residuals_OMIM.csv')
Transferred_Fisher_Residuals = pd.read_csv(path, index_col=0)
print(Transferred_Fisher_Residuals)

Transferred_Fisher_Residuals['preferential_expression_BH'] = multipletests(Transferred_Fisher_Residuals['preferential_expression'].values, method='fdr_bh')[1]

prediction_cols = [c for c in list(Transferred_Fisher_Residuals) if c not in ['PS_ID', 'Cell_Type', 'Phenotype Name', 'Disease_Phenotype', 'Standardized_Residuals', 'preferential_expression_BH']]
print('prediction_cols', len(prediction_cols), prediction_cols)

Transferred_Fisher_Residuals[prediction_cols] = -1 * np.log10(Transferred_Fisher_Residuals[prediction_cols] )
print(Transferred_Fisher_Residuals)


critical_cols = prediction_cols + ['Standardized_Residuals']
Transferred_Fisher_Residuals.dropna(subset=critical_cols, inplace=True)

"-------------------------------- Create Synthetic Dataset --------------------------------"


threshold = 3.45
Transferred_Fisher_Residuals['Disease_Label'] = False
Transferred_Fisher_Residuals.loc[Transferred_Fisher_Residuals['Standardized_Residuals'] > threshold, 'Disease_Label'] = True

Positives = Transferred_Fisher_Residuals[(Transferred_Fisher_Residuals['Disease_Label'] == True)]
Negatives = Transferred_Fisher_Residuals[Transferred_Fisher_Residuals['Disease_Label'] == False].sample(n=len(Positives)*9, random_state=1234)
pos_num = len(Positives)
cell_pos_num = len(Positives['Cell_Type'].unique())
diseases_num = len(Positives['PS_ID'].unique())
print('pos_num', pos_num)
print('cell_pos_num', cell_pos_num)
print('diseases_num', diseases_num)


All_Synthetic = pd.concat([Positives, Negatives])
print(All_Synthetic)

"----------------------------- Add Multi-cell gene model for comparison ----------------"

path = os.path.join('..', '..', 'Output', 'Transfer_Learning_Diseases', 'Summery', 'All_Diseases_Prediction_Long_Slim_Z_score_Genes.csv')
All_Diseases_Prediction_Long = pd.read_csv(path)
print(All_Diseases_Prediction_Long)
minimum_genes = 3
All_Diseases_Prediction_Long = All_Diseases_Prediction_Long[All_Diseases_Prediction_Long['Num Disease Genes'] >= minimum_genes]
print('Number of PS with more than 5 genes: ', len(All_Diseases_Prediction_Long['Disease PS ID'].unique()))
All_Synthetic = pd.merge(All_Synthetic, All_Diseases_Prediction_Long,  how='left', left_on=['PS_ID','Cell_Type'], right_on = ['Disease PS ID','Cell_Type'])
All_Synthetic.dropna(subset=['Z_Score_Probability'], inplace=True)


print(All_Synthetic)
print(All_Synthetic['Disease_Label'].sum())

z_threshold = 0.9
All_Synthetic['Z_Prediction'] = False
All_Synthetic.loc[All_Synthetic['Z_Score_Probability'] >= z_threshold, 'Z_Prediction'] = True

actual, predicted = All_Synthetic['Disease_Label'], All_Synthetic['Z_Prediction']

confusion_matrix = metrics.confusion_matrix(actual, predicted, normalize='true')
print('@ Multi-cell confusion_matrix')
print(confusion_matrix)
tn, fp, fn, tp = metrics.confusion_matrix(actual, predicted).ravel()
print('tn', tn, 'fp', fp, 'fn', fn, 'tp', tp)

Accuracy = metrics.accuracy_score(actual, predicted)
Precision = metrics.precision_score(actual, predicted)
Sensitivity_recall = metrics.recall_score(actual, predicted)
F1_score = metrics.f1_score(actual, predicted)

results_multi = {"Accuracy":Accuracy,"Precision":Precision,"Sensitivity_recall":Sensitivity_recall,"F1_score":F1_score}
print(results_multi)

cm_display = metrics.ConfusionMatrixDisplay(confusion_matrix = confusion_matrix, display_labels = [False, True])
cm_display.plot(cmap=plt.cm.Blues)
plt.show()
"----------------------------- Add Pref for comparison ----------------"

All_Synthetic['preferential_expression_Prediction'] = False
All_Synthetic.loc[All_Synthetic['preferential_expression_BH'] <= 0.1, 'preferential_expression_Prediction'] = True


actual, predicted = All_Synthetic['Disease_Label'], All_Synthetic['preferential_expression_Prediction']
confusion_matrix = metrics.confusion_matrix(actual, predicted, normalize='true')
print('@ Preferential method confusion_matrix')
print(confusion_matrix)
tn, fp, fn, tp = metrics.confusion_matrix(actual, predicted).ravel()
print('tn', tn, 'fp', fp, 'fn', fn, 'tp', tp)

Accuracy = metrics.accuracy_score(actual, predicted)
Precision = metrics.precision_score(actual, predicted)
Sensitivity_recall = metrics.recall_score(actual, predicted)
F1_score = metrics.f1_score(actual, predicted)
results_pref = {"Accuracy":Accuracy,"Precision":Precision,"Sensitivity_recall":Sensitivity_recall,"F1_score":F1_score}
print(results_pref)

cm_display = metrics.ConfusionMatrixDisplay(confusion_matrix = confusion_matrix, display_labels = [False, True])
cm_display.plot(cmap=plt.cm.Blues)
plt.show()
"------------------------------- Save Files ---------------------------"
"""
path = os.path.join('..', '..', 'Output', 'Transfer_Learning_Diseases', 'Summery', 'All_Synthetic_Data_Genes.csv')
All_Synthetic.to_csv(path)
# # plt.show()
"""

Results1 = pd.DataFrame([results_multi]).T
Results1.reset_index(inplace=True)
Results1.rename(columns={'index': 'Measurement', 0: 'Value'}, inplace=True)
Results1['Method'] = 'Multi-cellular model'
print(Results1)

Results2 = pd.DataFrame([results_pref]).T
Results2.reset_index(inplace=True)
Results2.rename(columns={'index': 'Measurement', 0: 'Value'}, inplace=True)
Results2['Method'] = 'Preferential based method'
print(Results2)
All_Results = pd.concat([Results1, Results2], ignore_index=True)
print(All_Results)
sns.factorplot(x='Measurement', y='Value', hue='Method', data=All_Results, kind='bar', hue_order=['Preferential based method', 'Multi-cellular model'], palette={'Multi-cellular model':'indianred', 'Preferential based method':'lightgrey'})
# sns.factorplot("Method", "Value", col="Measurement", data=All_Results, kind="bar")

path = os.path.join('..', '..', 'Output', 'Transfer_Learning_Diseases', 'Summery', 'All_Comparison_Results.pdf')
plt.savefig(path)

path = os.path.join('..', '..', 'Output', 'Transfer_Learning_Diseases', 'Summery', 'All_Comparison_Results.csv')
All_Results.to_csv(path, index=False)

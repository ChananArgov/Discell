import os
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import classification_report
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
import matplotlib.pyplot as plt
import sklearn as skl
from sklearn.metrics import average_precision_score
import sys
" -------------------------------------- Load Data ------------------------------"
file_name = sys.argv[1]
path= os.path.join('..', 'Input',  file_name)
Query_Genes = pd.read_csv(path, header=None)
print(Query_Genes)
query_genes = Query_Genes[0].tolist()
print(query_genes)

path= os.path.join('..', 'Data', 'scRNA_Features_Fraction.csv')
scRNA_Dataset = pd.read_csv(path)
print(scRNA_Dataset)

scRNA_Dataset['Query_Gene'] = False
scRNA_Dataset.loc[scRNA_Dataset['Gene_ID'].isin(query_genes), 'Query_Gene'] = True
print(scRNA_Dataset)


" ------------------------------------- Create Dataset ------------------------------"

folds = 9

Disease_Genes_Data = scRNA_Dataset[scRNA_Dataset['Query_Gene'] == True]
data_list = [Disease_Genes_Data]
counts = len(Disease_Genes_Data)

Non_Query_Data = scRNA_Dataset[scRNA_Dataset['Query_Gene'] == False].sample(n=counts * folds, axis='index', random_state=1234)
data_list.append(Non_Query_Data)


Query_Synthetic_Dataset = pd.concat(data_list)
print(Query_Synthetic_Dataset)
path= os.path.join('..', 'Output', file_name.replace('.csv', '') + '_Query_Synthetic_Dataset.csv')
Query_Synthetic_Dataset.to_csv(path, index=None)

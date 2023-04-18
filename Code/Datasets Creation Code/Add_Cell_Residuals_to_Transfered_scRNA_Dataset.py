import os
import pandas as pd

path = os.path.join('../..', 'Output', 'Text_Analysis_Output', 'Genes_Cells_PS_Standardized_Residuals_Max_Long.csv')
Cell_Residuals = pd.read_csv(path)

print(Cell_Residuals)
print(list(Cell_Residuals))
Cell_Residuals.rename(columns={'Ensembl Gene ID (Ensembl)':'Gene_ID'}, inplace=True)

path= os.path.join('..', '..', 'Output', 'Full_scRNA_Datasets', 'Transferred_scRNA_Dataset_Fraction.csv')
Transferred_scRNA = pd.read_csv(path)

print(Transferred_scRNA)
Transferred_scRNA_Residuals = Transferred_scRNA.merge(Cell_Residuals, on=['Gene_ID', 'Cell_Type'], how='left')
print(Transferred_scRNA_Residuals)

path= os.path.join('..', '..', 'Output', 'Full_scRNA_Datasets', 'Transfer_scRNA_Dataset_Fraction_Residuals.csv')
Transferred_scRNA_Residuals.to_csv(path, index=False)

import os
import pandas as pd

"-------------------------------- Load Data ---------------------------------"

path= os.path.join('../..', 'Data', 'GSE156793_S6_gene_expression_celltype.txt')
scRNAseq_Expression = pd.read_csv(path)
scRNAseq_Expression['Gene_ID'] = [x[0] for x in scRNAseq_Expression['RowID'].str.split('.')]
scRNAseq_Expression.set_index('Gene_ID', inplace=True)

"----------------------- Calculate Q1, Q3, IQR --------------------"

relevant_cols = [c for c in list(scRNAseq_Expression) if c not in ['RowID']]
scRNAseq_Expression = scRNAseq_Expression[relevant_cols]
print(scRNAseq_Expression)
Q1 = scRNAseq_Expression.quantile(0.25, axis=1)
Q3 = scRNAseq_Expression.quantile(0.75, axis=1)
IQR = Q3 - Q1
print('IQR', type(IQR), IQR)
IQR_fixed = IQR.copy(deep=True)
IQR_fixed.values[IQR_fixed.values < 1] = 1.0
print('IQR_fixed', type(IQR_fixed), IQR_fixed)

median_col = scRNAseq_Expression.median(axis=1)
print(median_col)
Temp_data2 = scRNAseq_Expression.copy(deep=True)
Temp_data2['Median_Expression'] = median_col
Temp_data2['IQR'] = IQR
Temp_data2['IQR_fixed'] = IQR_fixed
Temp_data2['Q1'] = Q1
Temp_data2['Q3'] = Q3


"----------------------- Calculate Preferential Expression --------------------"

PrePref_All = scRNAseq_Expression.sub(median_col, axis='index')
PrePref_All = PrePref_All.div(IQR_fixed, axis='index')
print(PrePref_All)

path = os.path.join('../..', 'Output', 'General_Preferential_Expression.csv')
PrePref_All.to_csv(path)

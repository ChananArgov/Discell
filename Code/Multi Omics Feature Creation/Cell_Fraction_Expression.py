import os
import pandas as pd

path= os.path.join('../..', 'Data', 'Gene_Fraction_Cell_Type_Filtered.csv')
Fraction_Cell = pd.read_csv(path)
print(Fraction_Cell)
relevant_genes = Fraction_Cell['Gene_ID'].tolist()

path= os.path.join('../..', 'Data', 'GSE156793_S6_gene_expression_celltype.csv')
Gene_scExpression = pd.read_csv(path)

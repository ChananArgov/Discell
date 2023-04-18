import os
import pandas as pd
import numpy as np

path = os.path.join('../..', 'Output', 'General_Preferential_PPI_Interactome.csv')
PPI_Net = pd.read_csv(path)

print(list(PPI_Net))
relevant_cols = [c for c in list(PPI_Net) if c not in ['Unnamed: 0', 'Weight']]
Relevant_PPI_Net = PPI_Net[relevant_cols]

Relevant_PPI_Net.loc[Relevant_PPI_Net['Protein_1'] > Relevant_PPI_Net['Protein_2'], 'Interaction'] = Relevant_PPI_Net['Protein_1'] + '_' + Relevant_PPI_Net['Protein_2']
Relevant_PPI_Net.loc[Relevant_PPI_Net['Protein_1'] < Relevant_PPI_Net['Protein_2'], 'Interaction'] = Relevant_PPI_Net['Protein_2'] + '_' + Relevant_PPI_Net['Protein_1']
Relevant_PPI_Net = Relevant_PPI_Net.drop_duplicates(subset='Interaction', keep="first")
print(Relevant_PPI_Net)

cell_types = list(Relevant_PPI_Net)
cell_types = [c for c in cell_types if c not  in ['Unnamed: 0', 'Protein_1', 'Protein_2', 'Weight', 'Interaction']]

path = os.path.join('../..', 'Data', 'GSE156793_S6_gene_expression_celltype.csv')
scRNA_Expression = pd.read_csv(path)
print(scRNA_Expression)

print('cell_types', len(cell_types), cell_types)
Relevant_PPI_Net = Relevant_PPI_Net[Relevant_PPI_Net['Interaction'].notna()]
c = 0
for cell in cell_types:
    c += 1
    print(cell,  ' number: ', c)
    non_expressed_genes = scRNA_Expression['Gene_ID'][scRNA_Expression[cell] < 1].tolist()
    non_expressed_genes = [x for x in non_expressed_genes if pd.isna(x) == False]
    print('non_expressed_genes', len(non_expressed_genes))
    Relevant_PPI_Net.loc[Relevant_PPI_Net['Interaction'].str.contains("|".join(non_expressed_genes)), cell] = None



path = os.path.join('../..', 'Output', 'Relevant_General_Preferential_PPI_Interactome.csv')
Relevant_PPI_Net.to_csv(path)
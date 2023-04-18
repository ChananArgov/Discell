'------------------------------- Imports -------------------------------------'

import os
import pandas as pd

pd.set_option('display.max_columns', None)

'------------------------------- Load Data ----------------------------------------------'

path = os.path.join('../..', '..', 'Data', 'Paralogs_data', 'Paralogs_data_biomart_Mars2021.txt')
Paralogs_Identity = pd.read_csv(path, sep='\t')
print(Paralogs_Identity)

path = os.path.join('../..', '..', 'Data', 'Protein_coding_genes.txt')
Protein_Coding = pd.read_csv(path)
print(Protein_Coding)

protein_coding = Protein_Coding['Gene stable ID'].unique().tolist()
print('protein_coding', len(protein_coding))

path = os.path.join('../..', '..', 'Data', 'Scince_03_11_2020', 'GSE156793_S6_gene_expression_celltype.csv')
scRNA_Exp = pd.read_csv(path)
print(scRNA_Exp)

Relevant_Paralogs = Paralogs_Identity[(Paralogs_Identity['Gene stable ID'].isin(protein_coding)) & (Paralogs_Identity['Human paralogue gene stable ID'].isin(protein_coding)) & (Paralogs_Identity['Paralogue %id. target Human gene identical to query gene'] >= 40) & (Paralogs_Identity['Paralogue %id. query gene identical to target Human gene'] >= 40)]
print(Relevant_Paralogs)
Relevant_Paralogs['Mean_Identity'] = (Paralogs_Identity['Paralogue %id. target Human gene identical to query gene'] + Paralogs_Identity['Paralogue %id. query gene identical to target Human gene'])/2#Relevant_Paralogs.mean(numeric_only=True)
print(Relevant_Paralogs)

cell_types = [c for c in list(scRNA_Exp) if c not in ['RowID', 'Gene_ID']]
print('cell_types', len(cell_types), cell_types)

'------------------------------- Main Code ----------------------------------------------'

cell_paralogs_ratio_list = []
relevant_features_list = []
n = 0
for cell in cell_types:
    n += 1
    Cell_Paralogs = pd.merge(Relevant_Paralogs, scRNA_Exp[['Gene_ID', cell]], left_on='Gene stable ID', right_on='Gene_ID', how='left')
    Cell_Paralogs.rename(columns={cell:'Gene_Expression'}, inplace=True)

    Cell_Paralogs = pd.merge(Cell_Paralogs, scRNA_Exp[['Gene_ID', cell]], left_on='Human paralogue gene stable ID', right_on='Gene_ID', how='left')
    Cell_Paralogs.rename(columns={cell:'Paralog_Expression'}, inplace=True)
    Cell_Paralogs.dropna(subset=['Gene_ID_x', 'Gene_ID_y'], inplace=True)

    Cell_Paralogs['Cell_Type'] = cell
    Cell_Paralogs['Gene_Paralog_Exp_Ratio'] = 0
    Cell_Paralogs.loc[(Cell_Paralogs['Paralog_Expression'] >= 1)&(Cell_Paralogs['Gene_Expression'] >= 1), 'Gene_Paralog_Exp_Ratio'] = Cell_Paralogs['Gene_Expression']/Cell_Paralogs['Paralog_Expression']
    Cell_Paralogs.loc[(Cell_Paralogs['Paralog_Expression'] < 1)&(Cell_Paralogs['Gene_Expression'] >= 1), 'Gene_Paralog_Exp_Ratio'] = Cell_Paralogs['Gene_Expression']/1.0

    Sum_Ratios = Cell_Paralogs.groupby('Gene stable ID', as_index=False)['Paralog_Expression'].sum()
    Sum_Ratios = pd.merge(Sum_Ratios, scRNA_Exp[['Gene_ID', cell]], left_on='Gene stable ID', right_on='Gene_ID', how='left')
    Sum_Ratios = Sum_Ratios[['Gene_ID', cell, 'Paralog_Expression']]
    Sum_Ratios.rename(columns={cell:'Gene_Expression', 'Paralog_Expression':'Sum_Paralog_Expression'}, inplace=True)
    Sum_Ratios['Cell_Type'] = cell
    Sum_Ratios['Max_Identity_Paralog_Expression'] = Cell_Paralogs['Paralog_Expression'][Cell_Paralogs.groupby('Gene stable ID', as_index=False)['Mean_Identity'].idxmax()].tolist()
    Sum_Ratios['_paralogs_ratio_highest_identity'] = Cell_Paralogs['Gene_Paralog_Exp_Ratio'][Cell_Paralogs.groupby('Gene stable ID', as_index=False)['Mean_Identity'].idxmax()].tolist()
    Sum_Ratios['_paralogs_ratio_all'] = 0
    Sum_Ratios.loc[(Sum_Ratios['Sum_Paralog_Expression'] >= 1)&(Sum_Ratios['Gene_Expression'] >= 1), '_paralogs_ratio_all'] = Sum_Ratios['Gene_Expression']/Sum_Ratios['Sum_Paralog_Expression']
    Sum_Ratios.loc[(Sum_Ratios['Sum_Paralog_Expression'] < 1)&(Sum_Ratios['Gene_Expression'] >= 1), '_paralogs_ratio_all'] = Sum_Ratios['Gene_Expression']/1.0
    print('Sum_Ratios')
    print(Sum_Ratios)
    cell_paralogs_ratio_list.append(Sum_Ratios)

    Relevant_Features = Sum_Ratios[['Gene_ID', '_paralogs_ratio_highest_identity', '_paralogs_ratio_all']].copy(deep=True)
    Relevant_Features.rename(columns={'_paralogs_ratio_highest_identity': cell+ '_paralogs_ratio_highest_identity','_paralogs_ratio_all': cell+ '_paralogs_ratio_all'}, inplace=True)
    Relevant_Features.set_index('Gene_ID', inplace=True)

    relevant_features_list.append(Relevant_Features)

'------------------------------- Sava Output ----------------------------------------------'

print(Relevant_Features)
All_relevant_paralogs = pd.concat(relevant_features_list, axis=1)
print(All_relevant_paralogs)
path = os.path.join('../..', '..', 'Data', 'Paralogs_data', 'Paralogs_Ratio_All_Devo_scRNA_cells_Features.csv')
All_relevant_paralogs.to_csv(path, index=True)

All_Ratios = pd.concat(cell_paralogs_ratio_list, axis=0)
print(All_Ratios)
path = os.path.join('../..', '..', 'Data', 'Paralogs_data', 'Paralogs_Ratio_All_Devo_scRNA_cells_Final.csv')
All_Ratios.to_csv(path, index=False)
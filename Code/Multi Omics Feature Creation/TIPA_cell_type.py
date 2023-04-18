import os
import pandas as pd
import numpy as np
from scipy import stats

"---------------------------- Load Data ----------------------------"
path= os.path.join('../..', 'Output', 'Tissue_Specific_Preferential', 'General_Preferential_Expression_acros_All_Cell_Types.csv')
scRNAseq_Pref = pd.read_csv(path)
print(scRNAseq_Pref)
cell_types_list = list(scRNAseq_Pref)
cell_types_list = [c for c in cell_types_list if c != 'Gene_ID']
print('cell_types_list', len(cell_types_list), cell_types_list)

path= os.path.join('../..', '..', 'data', 'GO_data', 'GO2Genes.txt')
GO2Gene = pd.read_csv(path, sep='\t')
print(GO2Gene)

go_list = list(GO2Gene['GO term accession'].unique())
print('go_list', len(go_list))

path= os.path.join('../..', '..', 'data', 'GO_data', 'GO_description_name.txt')
GO_description = pd.read_csv(path, sep='\t')
print(GO_description)


" ----------------------- Main ---------------------------------------"


def caculate_tipa(preferential_values):
    sorted_values = sorted(preferential_values)
    tipa_score = stats.trim_mean(sorted_values, 0.1)
    return tipa_score


go_tipa_cell_dict = {'GO_ID': [], 'GO_Name': [],'Cell_Type':[], 'TIPA':[]}
n = 1
for go in go_list:
    go_genes = list(GO2Gene['Gene stable ID'][GO2Gene['GO term accession'] == go].unique())
    try:
        go_name = GO_description['GO term name'][GO_description['GO term accession'] == go].values[0]
    except:
        go_name = None

    print('@', go, n)
    n += 1
    # print('go_name', go_name)
    # print('go_genes', go_genes)

    for cell in cell_types_list:
        go_cell_pref = scRNAseq_Pref[cell][scRNAseq_Pref['Gene_ID'].isin(go_genes)].tolist()
        tipa_score = caculate_tipa(go_cell_pref)

        go_tipa_cell_dict['GO_ID'].append(go)
        go_tipa_cell_dict['GO_Name'].append(go_name)
        go_tipa_cell_dict['Cell_Type'].append(cell)
        go_tipa_cell_dict['TIPA'].append(tipa_score)

        # print('go_cell_pref', go_cell_pref)
        # print('tipa_score', tipa_score)
        # print('\n')
    #     break
    # break

GO_TIPA_Cell = pd.DataFrame.from_dict(go_tipa_cell_dict, orient='columns')
print(GO_TIPA_Cell)
path= os.path.join('../..', 'Output', 'General_Pref_TIPA_Developmental_scRNA_Cell_Types.csv')
GO_TIPA_Cell.to_csv(path, index=False)
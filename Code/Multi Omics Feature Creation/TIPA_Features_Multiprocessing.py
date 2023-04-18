
'------------------------ Imports -------------------------------------------'
import os
import pandas as pd

import multiprocessing as mp
from time import gmtime, strftime

'------------------------- Load Data ----------------------------------------'

path = os.path.join('../..', 'Data', 'General_Pref_TIPA_Developmental_scRNA_Cell_Types_GO_Domain.csv')
GO_TIPA_Cell = pd.read_csv(path)
print(GO_TIPA_Cell)

path = os.path.join('../..', '..', 'Oncogene_Prediction', 'Data', 'Protein_coding_genes.csv')
Proteins = pd.read_csv(path)
print(Proteins)

path= os.path.join('../..', 'Data', 'GO2Genes.txt')
GO2Gene = pd.read_csv(path, sep='\t')
print(GO2Gene)

protein_coding_genes = Proteins['Gene stable ID'].unique().tolist()
print('protein_coding_genes', len(protein_coding_genes))

cell_types = GO_TIPA_Cell['Cell_Type'].unique().tolist()
print('cell_types', len(cell_types))

go_types = GO_TIPA_Cell['GO domain'].unique().tolist()
print('go_types', len(go_types), go_types)


def tipa_features(cell):
    print('@', cell)
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))

    Cell_data = GO_TIPA_Cell[GO_TIPA_Cell['Cell_Type'] == cell].copy(deep=True)
    cell_tipa_dict = {'Gene_ID': [], cell + '_tipa_min': [], cell + '_tipa_mean': [], cell + '_tipa_max': []}
    for gene in protein_coding_genes:
        gos = GO2Gene['GO term accession'][GO2Gene['Gene stable ID'] == gene].unique().tolist()
        TIPA_gene = Cell_data['TIPA'][(Cell_data['GO_ID'].isin(gos))&(Cell_data['GO domain'] == 'biological_process')]
        tipa_max = TIPA_gene.max()
        tipa_mean = TIPA_gene.mean()
        tipa_min = TIPA_gene.min()
        cell_tipa_dict['Gene_ID'].append(gene)
        cell_tipa_dict[cell + '_tipa_min'].append(tipa_min)
        cell_tipa_dict[cell + '_tipa_mean'].append(tipa_mean)
        cell_tipa_dict[cell + '_tipa_max'].append(tipa_max)
    Cell_Features = pd.DataFrame.from_dict(cell_tipa_dict, orient='columns')
    Cell_Features.set_index('Gene_ID', inplace=True)
    print(Cell_Features)
    return cell, Cell_Features

    
def driver_func_tipa():
    PROCESSES = 30
    df_list = []
    
    with mp.Pool(PROCESSES) as pool:
        cells = cell_types
        results = [pool.apply_async(tipa_features, (c,)) for c in cells]
        print('results', len(results), results)
    
        for r in results:
            
            results_tuple = r.get(timeout=None)
            print('#', results_tuple[0], ' finished')
            TIPA_cell = results_tuple[1]
            df_list.append(TIPA_cell)
            
        TIPA_All = pd.concat(df_list, axis=1, join="inner")
        path = os.path.join('../..', 'Data', 'TIPA_General_Preferential_Biological_Processes_Features.csv')
        TIPA_All.to_csv(path)
        
if __name__ == '__main__':
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    driver_func_tipa()    
    




'------------------------ Imports -------------------------------------------'

import os
import pandas as pd
import multiprocessing as mp
from time import gmtime, strftime

'------------------------- Load Data ----------------------------------------'

path = os.path.join('../..', 'Data', 'Relevant_General_Preferential_PPI_Interactome.csv')
Relevant_PPI_Diffnet = pd.read_csv(path)
print(Relevant_PPI_Diffnet)

path = os.path.join('../..', 'Data', 'General_Preferential_Expression.csv')
General_Pref = pd.read_csv(path)
print(General_Pref)

path = os.path.join('../..', '..', 'Oncogene_Prediction', 'Data', 'Protein_coding_genes.csv')
Protein_genes = pd.read_csv(path)  # low_memory=False,
print(Protein_genes)

protein_coding = Protein_genes['Gene stable ID'].unique().tolist()
print('protein_coding', len(protein_coding))

cell_types = list(Relevant_PPI_Diffnet)
cell_types = [c for c in cell_types if c not in ['Protein_1', 'Protein_2', 'Unnamed: 0', 'Interaction']]
print('cell_types', len(cell_types), cell_types)
Relevant_General_Pref = General_Pref[General_Pref['Gene_ID'].isin(protein_coding)]
print(Relevant_General_Pref)

'------------------------- Main Code ----------------------------------------'

def remove_non_pref(cell):
    print(cell)
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    threshold_95 = Relevant_General_Pref[cell].quantile(0.95)
    non_pref_genes = Relevant_General_Pref['Gene_ID'][Relevant_General_Pref[cell] < threshold_95].tolist()
    non_pref_genes = [x for x in non_pref_genes if pd.isna(x) == False]
    print('non_pref_genes', len(non_pref_genes))

    Cell_Data = Relevant_PPI_Diffnet[['Interaction', cell]]
    Cell_Data.loc[Cell_Data['Interaction'].str.contains("|".join(non_pref_genes)), cell] = None
    Cell_Data.set_index('Interaction', inplace=True)

    return cell, Cell_Data


def driver_func3():
    PROCESSES = 40
    top95_results_list = []
    results = []
    with mp.Pool(PROCESSES) as pool:
        cells = cell_types
        #         print(cells)
        results = [pool.apply_async(remove_non_pref, (c,)) for c in cells]
        #         results.append(pool.apply_async(create_net_features, c) for c in cells)
        print('results', len(results), results)

        for r in results:
            results_tuple = r.get(timeout=None)
            print('#', results_tuple[0], ' finished')
            Top95_cell = results_tuple[1]
            top95_results_list.append(Top95_cell)

        TOP95_Pref = pd.concat(top95_results_list, axis=1, join="inner")
        path = os.path.join('../..', 'Output', 'Only_Top95_General_Preferential_PPI_scRNA_Interactome.csv')
        TOP95_Pref.to_csv(path)


if __name__ == '__main__':
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    driver_func3()
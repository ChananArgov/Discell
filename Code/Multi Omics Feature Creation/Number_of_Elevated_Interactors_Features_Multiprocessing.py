
'------------------------ Imports -------------------------------------------'
import os
import pandas as pd

import multiprocessing as mp
from time import gmtime, strftime


'------------------------- Load Data ----------------------------------------'


path = os.path.join('../..', 'Data', 'Only_Top95_General_Preferential_PPI_scRNA_Interactome.csv')
TOP95_Pref = pd.read_csv(path)
print(TOP95_Pref)

path = os.path.join('../..', '..', 'Oncogene_Prediction', 'Data', 'Protein_coding_genes.csv')
Protein_genes = pd.read_csv(path)  # low_memory=False,
print(Protein_genes)

protein_coding = Protein_genes['Gene stable ID'].unique().tolist()
print('protein_coding', len(protein_coding))

cell_types = list(TOP95_Pref)
cell_types = [c for c in cell_types if c not in ['Protein_1', 'Protein_2', 'Unnamed: 0', 'Interaction']]
print('cell_types', len(cell_types), cell_types)

'------------------------- Main Code ----------------------------------------'

def create_PrefNet_features(cell):
    
    print('@', cell)
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))

    pref_int_results_dict = {'Gene_ID': [], cell + '_Num_Pref_Int': []}
    c = 0
    for protein in protein_coding:

        Protein_Data = TOP95_Pref[TOP95_Pref['Interaction'].str.contains(protein)]
        num_interactors = Protein_Data[cell].count()
        
        pref_int_results_dict['Gene_ID'].append(protein)
        pref_int_results_dict[cell + '_Num_Pref_Int'].append(num_interactors)


    PrefInt_results = pd.DataFrame.from_dict(pref_int_results_dict, orient='columns')
    PrefInt_results.set_index('Gene_ID', inplace=True)
    
    return cell, PrefInt_results
    
            
def driver_func2():
    PROCESSES = 40
    pref_int_results_list = []
    results = []
    with mp.Pool(PROCESSES) as pool:
        cells = cell_types
        results = [pool.apply_async(create_PrefNet_features, (c,)) for c in cells]
        print('results', len(results), results)
    
        for r in results:
            print(r)
            #https://stackoverflow.com/questions/40773925/where-is-documentation-for-multiprocessing-pool-applyresult
            results_tuple = r.get(timeout=None)
            print('#', results_tuple[0], 'finished')
            PrefInt_results = results_tuple[1]
            pref_int_results_list.append(PrefInt_results)
            
        Pref_Int_All = pd.concat(pref_int_results_list, axis=1, join="inner")
        print(Pref_Int_All)

        path = os.path.join('../..', 'Output', 'Number_of_Preferential_Int_Devo_scRNA.csv')
        Pref_Int_All.to_csv(path)#, index=False

    
if __name__ == '__main__':
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    driver_func2()





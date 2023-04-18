

'------------------------ Imports -------------------------------------------'

import os
import pandas as pd

'------------------------- Load Data ----------------------------------------'

path = os.path.join('../..', 'Data', 'Relevant_General_Preferential_PPI_Interactome.csv')
Relevant_PPI_Diffnet = pd.read_csv(path)
print(Relevant_PPI_Diffnet)

path = os.path.join('../..', '..', 'Oncogene_Prediction', 'Data', 'Protein_coding_genes.csv')
Protein_genes = pd.read_csv(path)#low_memory=False,
print(Protein_genes)

protein_coding = Protein_genes['Gene stable ID'].unique().tolist()
print('protein_coding', len(protein_coding))

cell_types = list(Relevant_PPI_Diffnet)
cell_types = [c for c in cell_types if c not in ['Protein_1', 'Protein_2', 'Unnamed: 0', 'Interaction']]
print('cell_types', len(cell_types), cell_types)



# Example by - https://stackoverflow.com/questions/50937362/multiprocessing-on-python-3-jupyter
# https://medium.com/swlh/5-step-guide-to-parallel-processing-in-python-ac0ecdfcea09#id_token=eyJhbGciOiJSUzI1NiIsImtpZCI6Ijc4M2VjMDMxYzU5ZTExZjI1N2QwZWMxNTcxNGVmNjA3Y2U2YTJhNmYiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJodHRwczovL2FjY291bnRzLmdvb2dsZS5jb20iLCJuYmYiOjE2MTExMzY3MzksImF1ZCI6IjIxNjI5NjAzNTgzNC1rMWs2cWUwNjBzMnRwMmEyamFtNGxqZGNtczAwc3R0Zy5hcHBzLmdvb2dsZXVzZXJjb250ZW50LmNvbSIsInN1YiI6IjEwOTAwMTEwMDgxNzE4MDc5MTM4MSIsImhkIjoicG9zdC5iZ3UuYWMuaWwiLCJlbWFpbCI6ImNoYW5hbmFAcG9zdC5iZ3UuYWMuaWwiLCJlbWFpbF92ZXJpZmllZCI6dHJ1ZSwiYXpwIjoiMjE2Mjk2MDM1ODM0LWsxazZxZTA2MHMydHAyYTJqYW00bGpkY21zMDBzdHRnLmFwcHMuZ29vZ2xldXNlcmNvbnRlbnQuY29tIiwibmFtZSI6IkNoYW5hbiBBcmdvdiIsInBpY3R1cmUiOiJodHRwczovL2xoMy5nb29nbGV1c2VyY29udGVudC5jb20vLU5jRnVXWlNvcDQ4L0FBQUFBQUFBQUFJL0FBQUFBQUFBQUFBL0FNWnV1Y25zUDJ1QlM3Q3Z1TUZhQ2twQ0pIci1mV1dmTmcvczk2LWMvcGhvdG8uanBnIiwiZ2l2ZW5fbmFtZSI6IkNoYW5hbiIsImZhbWlseV9uYW1lIjoiQXJnb3YiLCJpYXQiOjE2MTExMzcwMzksImV4cCI6MTYxMTE0MDYzOSwianRpIjoiMmJhODAyZDZkOGVhZmMxMmRmY2ExZWE5YjA0ZmM1ODk1NzI1OGFlZiJ9.kNekY55wVAJ7YZ391HXJcxMUJALdFQ5rBAe3hhQqk-RtY52qb4_gCQnqS9he6NbD_-Rs5vK8ZBMd15ZbCT1NmxennQHx6cMAESjxT_Wz5QA7UyD-km9AzIGNhiS9m2GkecKd_hJUks5i623kB5SxsQi4CIj5aOAb5k9krIr3jU3DNlF7DCssckeVz-I6z0bzdKBd2rqsgPgBemrSj_xazzNi8D_LB_Ph5TA_KL109VnEkNu-vj4zSfyXjYNBrYqorEJHr84efw6ZuHpSxbC0sgrEQtM17YbVzE0NlCsIH-Pt4Kq-V97TchjVyEWwCxiGVJcCVlTOtKByFnkNTj_geA
import multiprocessing as mp
from time import gmtime, strftime


def create_net_features(cell):
    
    print('@', cell)
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))

    diff_results_dict = {'Gene_ID':[],  cell + '_DiffNet_Mean':[], cell + '_DiffNet_Max':[],cell + '_DiffNet_Min':[],}
    int_results_dict = {'Gene_ID': [], cell + '_Num_Interactions': []}
    c = 0
    for protein in protein_coding:

        Protein_Data = Relevant_PPI_Diffnet[(Relevant_PPI_Diffnet['Protein_1'] == protein)|(Relevant_PPI_Diffnet['Protein_2'] == protein)]
        num_interactors = Protein_Data[cell].count()
        diff_mean = Protein_Data[cell].mean()
        diff_max = Protein_Data[cell].max()
        diff_min = Protein_Data[cell].min()

        int_results_dict['Gene_ID'].append(protein)
        int_results_dict[cell + '_Num_Interactions'].append(num_interactors)

        diff_results_dict['Gene_ID'].append(protein)
        diff_results_dict[cell + '_DiffNet_Mean'].append(diff_mean)
        diff_results_dict[cell + '_DiffNet_Max'].append(diff_max)
        diff_results_dict[cell + '_DiffNet_Min'].append(diff_min)
        
    Diff_results = pd.DataFrame.from_dict(diff_results_dict, orient='columns')
    Diff_results.set_index('Gene_ID', inplace=True)

    Int_results = pd.DataFrame.from_dict(int_results_dict, orient='columns')
    Int_results.set_index('Gene_ID', inplace=True)
    
    return Diff_results, Int_results

            
def driver_func2():

    PROCESSES = 40 # The number of CPUs to use
    diff_results_list = []
    int_results_list = []

    with mp.Pool(PROCESSES) as pool:
        cells = cell_types
        results = [pool.apply_async(create_net_features, (c,)) for c in cells]
        print('results', len(results), results)
    
        for r in results:
            print(r)
            #https://stackoverflow.com/questions/40773925/where-is-documentation-for-multiprocessing-pool-applyresult
            results_tuple = r.get(timeout=None)
            print('\t', results_tuple)
            Diff_results = results_tuple[0]
            Int_results = results_tuple[1]
            diff_results_list.append(Diff_results)
            int_results_list.append(Int_results)
            
        Int_Results = pd.concat(int_results_list, axis=1, join="inner")
        Diff_Results = pd.concat(diff_results_list, axis=1, join="inner")
        print(Int_Results)
        print(Diff_Results)

        path = os.path.join('../..', 'Number_of_Interactions_Devo_scRNA_2.csv')
        Int_Results.to_csv(path)#, index=False

        path = os.path.join('../..', 'Output', 'DiffNet_Features_Devo_scRNA_2.csv')
        Diff_Results.to_csv(path)#, index=False
    
if __name__ == '__main__':
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    driver_func2()
    






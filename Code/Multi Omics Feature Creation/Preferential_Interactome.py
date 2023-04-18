import os
import pandas as pd

"-------------------------------- Load Data --------------------------"
pd.set_option('display.max_column',None)

path= os.path.join('../..', 'Output', 'General_Preferential_Expression.csv')
General_Pref= pd.read_csv(path)

print(General_Pref)

path= os.path.join('../..', 'Data', 'GlobalInteractome.tsv')
PPI_Net = pd.read_csv(path, '\t',  header=None)

print(PPI_Net)
print(PPI_Net.describe())

cell_types = list(General_Pref)
cell_types = [c for c in cell_types if c != 'Gene_ID']
print('cell_types', len(cell_types), cell_types)

Relevant_PPI_Net = PPI_Net[[0,1,5]].copy(deep=True)
Relevant_PPI_Net.rename(columns={0: 'Protein_1', 1: 'Protein_2', 5: 'Weight'}, inplace=True)
print(Relevant_PPI_Net)

def pref_interaction(cell):
    Temp_interactions = Relevant_PPI_Net[['Protein_1', 'Protein_2']].copy(deep=True)
    Temp_interactions = pd.merge(Temp_interactions, General_Pref[['Gene_ID', cell]], how="left", left_on='Protein_1', right_on='Gene_ID')
    Temp_interactions.rename(columns={cell:'Pref_1'}, inplace=True)
    Temp_interactions = pd.merge(Temp_interactions, General_Pref[['Gene_ID', cell]], how="left", left_on='Protein_2', right_on='Gene_ID')
    Temp_interactions.rename(columns={cell:'Pref_2'}, inplace=True)
    Temp_interactions[cell] = Temp_interactions['Pref_1'] + Temp_interactions['Pref_2']
    return Temp_interactions[cell]

for cell in cell_types:
    print('@', cell)
    Relevant_PPI_Net[cell] = pref_interaction(cell)

print(Relevant_PPI_Net)
path = os.path.join('../..', 'Output', 'General_Preferential_PPI_Interactome.csv')
Relevant_PPI_Net.to_csv(path)
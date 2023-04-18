import pandas as pd
import os
import plotly.graph_objects as go
from pyvis.network import Network
import networkx as nx
import matplotlib as mpl
import sys

mpl.rcParams['axes.spines.left'] = True
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.bottom'] = True

"---------------------------- Load Data --------------------------------------------"

path= os.path.join('..',  'Data', 'Edited_OMIM_Morbidmap_Phenotypic_Series.csv')
PS_OMIM_Data = pd.read_csv(path, index_col=0)
print(PS_OMIM_Data)

path= os.path.join('..', 'Data', 'Only_Top95_General_Preferential_PPI_scRNA_Interactome (1).csv')
Top95_Interactom = pd.read_csv(path)
print(Top95_Interactom)

path= os.path.join('..',  'Data', 'Relevant_General_Preferential_PPI_Interactome_Cell_Development_expressed_only.csv')
Interactom = pd.read_csv(path, index_col=0)
print(Interactom)

path= os.path.join('..', 'Data', 'General_Pref_TIPA_Developmental_scRNA_Cell_Types_GO_Domain.csv')
Processes_data = pd.read_csv(path)
print(Processes_data)

path= os.path.join('..', 'Data', 'GO2Genes.txt')
GO2Gene = pd.read_csv(path, sep='\t')
print(GO2Gene)

path= os.path.join('..', 'Data', 'GO_description_name.txt')
GO_description = pd.read_csv(path, sep='\t')
print(GO_description)

path= os.path.join('..', 'Data', 'Gene2Name_May2022.txt')
Gene2Name = pd.read_csv(path, sep='\t')
print('Gene2Name')
print(Gene2Name)
"----------------------- Get Input -----------------------------"

input_path = sys.argv[1] #os.path.join('..', 'Input', 'Hypuerinsulimic_Example.txt')
Input_Example = pd.read_csv(input_path)
relevant_genes = Input_Example['Gene_ID'].tolist()
tissue = sys.argv[2] #'Pancreas'
cell = sys.argv[3] #'Islet endocrine cells'
cell = cell.replace('@', ' ')
affected_cell = tissue + '-' + cell
print('sys.argv[3]', sys.argv[3])
print('affected_cell', affected_cell)


" ------------------ Predict Disease Genes in Cellular context ------------- "

Affected_Cell_Interactom = Interactom[['Protein_1',	'Protein_2', affected_cell, 'Interaction']]
print(Affected_Cell_Interactom)
Disease_Genes_Interactom = Affected_Cell_Interactom[(Affected_Cell_Interactom['Interaction'].str.contains('|'.join(relevant_genes))) & (~Affected_Cell_Interactom[affected_cell].isnull())]
print(Disease_Genes_Interactom)
threshold = Affected_Cell_Interactom[affected_cell].quantile(0.90)
print('threshold', threshold)

Disease_Genes_Processes = GO2Gene[(GO2Gene['Gene stable ID'].isin(relevant_genes)) & (~GO2Gene['GO term accession'].isna())]
print(Disease_Genes_Processes)
Disease_Genes_Processes_Counts = pd.DataFrame(Disease_Genes_Processes['GO term accession'].value_counts())
Disease_Genes_Processes_Counts.reset_index(inplace=True)
Disease_Genes_Processes_Counts.rename(columns={'GO term accession':'Number of disease genes', 'index':'GO term accession'}, inplace=True)
print(Disease_Genes_Processes_Counts)
Disease_Genes_Processes_Counts = Disease_Genes_Processes_Counts[Disease_Genes_Processes_Counts['Number of disease genes'] >= 3]
Disease_Genes_Processes_Counts = pd.merge(Disease_Genes_Processes_Counts, Processes_data[(Processes_data['Cell_Type'] == affected_cell)&(Processes_data['GO domain'] == 'biological_process')], left_on='GO term accession', right_on='GO_ID', how='inner')
Disease_Genes_Processes_Counts.sort_values(by=['Number of disease genes', 'TIPA'], inplace=True, ascending=[False, False])
print(Disease_Genes_Processes_Counts)
relevant_processes = Disease_Genes_Processes_Counts['GO term accession_x'][Disease_Genes_Processes_Counts['TIPA'] >= 1.5].tolist()
print(relevant_processes)
all_processes_genes = GO2Gene['Gene stable ID'][GO2Gene['GO term accession'].isin(relevant_processes)].tolist()

processes_genes_dict = {}
for pro in relevant_processes:
    process_genes = GO2Gene['Gene stable ID'][GO2Gene['GO term accession'] == pro].tolist()
    processes_genes_dict[pro] = process_genes

all_interactors_list = Disease_Genes_Interactom['Protein_1'].tolist() + Disease_Genes_Interactom['Protein_2'].tolist()
all_interactors_list = set(all_interactors_list)
print('all_interactors_list', len(all_interactors_list), all_interactors_list)

Interactors_data = pd.DataFrame(all_interactors_list, columns=['Interactor'])
Interactors_data['color'] = 'lightsteelblue'
Interactors_data.loc[Interactors_data['Interactor'].isin(relevant_genes), 'color'] = 'crimson'
print(Interactors_data)

N = Network(directed=False, height='100%', width='100%')  # height='100%', width='100%', bgcolor='#222222', font_color='white',

for int in all_interactors_list:
    title = ''
    label = int

    if int in Gene2Name['Gene stable ID'].tolist():
        label = Gene2Name['Gene name'][Gene2Name['Gene stable ID'] == int].values[0]

    if int in relevant_genes:  # if the node is part of the sub-graph
        color = 'crimson'
    elif int in all_processes_genes:
        color = 'yellow'
    else:
        color = 'lightsteelblue'
    for process in processes_genes_dict:
        if int in processes_genes_dict[process]:
            tipa = Disease_Genes_Processes_Counts['TIPA'][Disease_Genes_Processes_Counts['GO term accession_x'] == process].values[0]
            title += GO_description['GO term name'][GO_description['GO term accession'] == process].values[0]
            title += ' = ' + str(round(tipa, 2))
            title += ', '
    N.add_node(int, label=label, color=color, title=title.rstrip(', '))

Disease_Genes_Interactom.rename(columns={affected_cell: 'weight'}, inplace=True)
Disease_Genes_Interactom['width'] = 2
Disease_Genes_Interactom.loc[Disease_Genes_Interactom['weight'] > 0, 'width'] = 4
Disease_Genes_Interactom.loc[Disease_Genes_Interactom['weight'] >= threshold, 'width'] = 6

Disease_Genes_Interactom['color'] = 'lightsteelblue'
Disease_Genes_Interactom.loc[Disease_Genes_Interactom['weight'] > 0, 'color'] = 'lightcoral'
Disease_Genes_Interactom.loc[Disease_Genes_Interactom['weight'] >= threshold, 'color'] = 'crimson'

G = nx.from_pandas_edgelist(Disease_Genes_Interactom, 'Protein_1', 'Protein_2', edge_attr=True)
N.from_nx(G)
path = os.path.join('..', 'Output', 'Example_Network_Processes.html')
N.write_html(path)  # write_html -  save a html file in current dir

path = os.path.join('..', 'Output',  'Example_Processes.csv')
Disease_Genes_Processes_Counts.to_csv(path, index=False)


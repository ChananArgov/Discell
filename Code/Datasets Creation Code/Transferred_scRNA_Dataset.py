import pandas as pd
import os

# path= os.path.join('..', '..', 'Output', 'Full_scRNA_Datasets', 'Fix_Full_Features_dataset_scRNA_Development_TIPA_processes.csv')
path= os.path.join('..', '..', 'Output', 'Full_scRNA_Datasets', 'scRNA_Features_Fraction.csv')

scRNA_Features_dataset = pd.read_csv(path)
scRNA_Features_dataset.rename(columns={'Unnamed: 0': 'Gene_ID'}, inplace=True)
scRNA_Features_dataset.set_index('Gene_ID', inplace=True)
print(scRNA_Features_dataset)

cell_list = [c.replace('_expression', '') for c in list(scRNA_Features_dataset) if ('_expression' in c)&('_preferential' not in c)]
print('cell_list', len(cell_list), cell_list)
cell_data_list = []
for cell in cell_list:
    print('@', cell)
    cell_features = [f for f in list(scRNA_Features_dataset) if cell in f]
    print('cell_features', len(cell_features), cell_features)
    Cell_data = scRNA_Features_dataset[cell_features]
    Cell_data.rename(columns={f:f.replace(cell + '_', '') for f in cell_features}, inplace=True)
    Cell_data.reset_index(inplace=True)
    # print(Cell_data)
    Cell_data['Cell_Type'] = cell
    cell_data_list.append(Cell_data)
    # break

Transferred_data = pd.concat(cell_data_list, axis=0)
print(Transferred_data)
path= os.path.join('..', '..', 'Output', 'Full_scRNA_Datasets', 'Transferred_scRNA_Dataset_Fraction.csv')
Transferred_data.to_csv(path, index=False)

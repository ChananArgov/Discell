import pandas as pd
import os
'------------------------------- Load Files ----------------------------------------------'

path= os.path.join('..', '..', 'Output', 'Cell_Type_Features', 'TIPA_General_Preferential_Biological_Processes_Features.csv')
Tipa_Features = pd.read_csv(path)
Tipa_Features.set_index('Gene_ID', inplace=True)
print(Tipa_Features)
print('Tipa_Features', len(list(Tipa_Features)))#, list(Tipa_Features))
print('Tipa', len(Tipa_Features))

path= os.path.join('..', '..', 'Data', 'Protein_coding_genes.txt')
Proteins= pd.read_csv(path)
print(Proteins)
protein_coding_genes = Proteins['Gene stable ID'].unique().tolist()
print('protein_coding_genes', len(protein_coding_genes))

path= os.path.join('..', '..', 'Data', 'Scince_03_11_2020', 'GSE156793_S6_gene_expression_celltype.csv')
Expression_Features = pd.read_csv(path)
relevant_features = [f for f in list(Expression_Features) if f not in ['RowID']]
Expression_Features = Expression_Features[relevant_features]
Expression_Features = Expression_Features[Expression_Features['Gene_ID'].isin(protein_coding_genes)]
Expression_Features.set_index('Gene_ID', inplace=True)
Expression_Features.columns += '_expression'
print(Expression_Features)
print('Expression_Features', len(list(Expression_Features)))#, list(Expression_Features))
print('Expression_Features', len(Expression_Features))

path= os.path.join('..', '..', 'Data', 'Scince_03_11_2020', 'GSE156793_S7_gene_fraction_celltype.csv')
Fraction_Features = pd.read_csv(path)
print(Fraction_Features)
relevant_features = [f for f in list(Fraction_Features) if f not in ['RowID']]
Fraction_Features = Fraction_Features[relevant_features]
Fraction_Features = Fraction_Features[Fraction_Features['Gene_ID'].isin(protein_coding_genes)]
Fraction_Features.set_index('Gene_ID', inplace=True)
Fraction_Features.columns += '_Fraction'
print(Fraction_Features)
print('Fraction_Features', len(list(Fraction_Features)))#, list(Expression_Features))
print('Fraction_Features', len(Fraction_Features))

path= os.path.join('..', '..', 'Output', 'Tissue_Specific_Preferential', 'General_Preferential_Expression_acros_All_Cell_Types.csv')
Preferential_Features = pd.read_csv(path)
Preferential_Features = Preferential_Features[Preferential_Features['Gene_ID'].isin(protein_coding_genes)]
Preferential_Features.set_index('Gene_ID', inplace=True)
Preferential_Features.columns += '_preferential_expression'
print(Preferential_Features)
print('Preferential_Features', len(list(Preferential_Features)))#, list(Preferential_Features))
print('Preferential_Features', len(Preferential_Features))


path= os.path.join('..', '..', 'Output', 'Network_Features')
network_features_files = os.listdir(path)
print(network_features_files)

suffixes = ['_num_interactors_dif_mean', '_interactors_dif_med', '_num_elevated_interactors_dif_mean', '_interactors_dif_med', '_Num_Specific_Int', '_num_specific_interactions_dif_mean', '_num_specific_interactions_dif_med', '_Num_Interactions', '_Num_Pref_Int']
replace_suffix = ['_num_interactors_dif_mean', '_interactors_dif_med', '_num_elevated_interactors_dif_mean', '_interactors_dif_median', '_num_specific_interactions', '_num_specific_interactions_dif_mean', '_num_specific_interactions_dif_median',  '_num_interactors', '_num_elevated_interactors']

diffnet_suffix = [ '_DiffNet_Mean', '_DiffNet_Max', '_DiffNet_Min']
replace_diffnet_suffix = ['_diff_net_mean', '_diff_net_max', '_diff_net_min']

net_features_list = []
c = 0
for file in network_features_files:
    path = os.path.join('..', '..', 'Output', 'Network_Features', file)
    file_name = file.replace('.csv', '')
    print('@', file_name)
    print('path', path)
    Net_Features = pd.read_csv(path)
    Net_Features.set_index('Gene_ID', inplace=True)
    # print(Net_Features)
    print('file_name', len(list(Net_Features)), list(Net_Features))
    if file_name != 'DiffNet_Features_Devo_scRNA_2':
        replace_cols_dict = {col:col.replace(suffixes[c], replace_suffix[c]) for col in list(Net_Features)}
        Net_Features.rename(columns=replace_cols_dict, inplace=True)
        c += 1
    elif file_name == 'DiffNet_Features_Devo_scRNA_2':
        c2 = 0
        for suf in diffnet_suffix:
            replace_cols_dict = {col: col.replace(diffnet_suffix[c2], replace_diffnet_suffix[c2]) for col in list(Net_Features)}
            Net_Features.rename(columns=replace_cols_dict, inplace=True)
            c2 += 1

    net_features_list.append(Net_Features)
    print('file_name', len(list(Net_Features)))#, list(Net_Features))
    print('file_name', len(Net_Features))


All_Net_Features = pd.concat(net_features_list, axis=1, join="inner")
print('All_Net_Features', len(list(All_Net_Features)))#, list(All_Net_Features))
print('All_Net_Features', len(All_Net_Features))
print(All_Net_Features)

path= os.path.join('..', '..', 'Data', 'Paralogs_data', 'Paralogs_Ratio_All_Devo_scRNA_cells_Features.csv')
Paralogs_Features = pd.read_csv(path)
Paralogs_Features = Paralogs_Features[Paralogs_Features['Gene_ID'].isin(protein_coding_genes)]
Paralogs_Features.set_index('Gene_ID', inplace=True)
print(Paralogs_Features)
print('Paralogs_Features', len(list(Paralogs_Features)))#, list(Paralogs_Features))
print('Paralogs_Features', len(Paralogs_Features))

Full_Features_dataset = pd.concat([Expression_Features, Fraction_Features, Preferential_Features,  Tipa_Features, All_Net_Features, Paralogs_Features], axis=1)
print(Full_Features_dataset)
print('Full_Features_dataset', len(list(Full_Features_dataset)))#, list(Full_Features_dataset))

path= os.path.join('..', '..', 'Output', 'Full_scRNA_Datasets', 'scRNA_Features_Fraction.csv')
Full_Features_dataset.to_csv(path)
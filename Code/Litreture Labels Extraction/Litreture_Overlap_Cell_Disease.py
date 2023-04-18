import os
import pandas as pd
from scipy.stats import fisher_exact

path = os.path.join('../..', 'Output', 'Text_Analysis_Output', 'Cells_Development_PubMed_Occurrence.tsv')
Cell_Occurrence = pd.read_csv(path, sep='\t')
print(Cell_Occurrence)

path = os.path.join('../..', 'Output', 'Text_Analysis_Output', 'Phenotypic_Series_PubMed_Occurrence.tsv')
Disease_Occurrence = pd.read_csv(path, sep='\t')
print(Disease_Occurrence)

cell_list = Cell_Occurrence['Term'].to_list()
disease_list = Disease_Occurrence['Term'].to_list()

print('cell_list', len(cell_list), cell_list)
print('disease_list', len(disease_list), disease_list)
overlap_occurrence_dict = {'Cell_Type':[], 'Disease_Phenotype':[], 'PubMed_Overlaps':[], 'PubMed_IDs':[]}
print(Disease_Occurrence['PubMed_IDs'][Disease_Occurrence['Term'] == 'Short-rib thoracic dysplasia'])
for cell in cell_list:
    print('@', cell)
    for disease in disease_list:
        print(' #', disease, cell)
        try:
            cell_pmd_list = Cell_Occurrence['PubMed_IDs'][Cell_Occurrence['Term'] == cell].values[0].split('_')
            disease_pmd_list = Disease_Occurrence['PubMed_IDs'][Disease_Occurrence['Term'] == disease].values[0].split('_')
            intersection = [i for i in cell_pmd_list if i in disease_pmd_list]
            num_overlap = len(intersection)
        except:
            intersection = None
            num_overlap = None
            print('*** Cell: ', cell, 'or disease:', disease, 'have problem')
        print(' intersection', num_overlap)

        overlap_occurrence_dict['Cell_Type'].append(cell)
        overlap_occurrence_dict['Disease_Phenotype'].append(disease)
        overlap_occurrence_dict['PubMed_Overlaps'].append(num_overlap)
        overlap_occurrence_dict['PubMed_IDs'].append(intersection)

PubMed_Overlaps = pd.DataFrame.from_dict(overlap_occurrence_dict, orient='columns')
path = os.path.join('../..', 'Output', 'Text_Analysis_Output', 'Cell_Type_Phenotypic_Series_PubMed_Overlaps.csv')
PubMed_Overlaps.to_csv(path, index=False)

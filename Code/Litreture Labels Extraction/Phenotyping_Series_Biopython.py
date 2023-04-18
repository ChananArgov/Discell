import os
import pandas as pd
from Bio import Entrez

" ------------------------------- Load Disease Phenotypes and Cell Names --------------------- "

path = os.path.join('../..', 'Data', 'Diseases', 'Phenotypic_Series_OMIM_Names_Oct2021.csv')
PS_Names = pd.read_csv(path)
print(PS_Names)

ps_name_list = PS_Names['Phenotype Name'].to_list()
comma_names = [ps for ps in ps_name_list if ',' in ps]

print('ps_name_list', len(ps_name_list), ps_name_list)
print('comma_names', len(comma_names), comma_names)

path = os.path.join('../..', 'Output', 'TIPA_Cell_Type', 'General_Pref_TIPA_for_Diseases_Developmental_scRNA_Cell_Types_Genes.csv')
TIPA_data = pd.read_csv(path)
print(TIPA_data)

cell_types = TIPA_data['Cell_Type'].unique().tolist()
print('cell_types', len(cell_types), cell_types)

"------------------ Biopython Functions -----------------------------"

from Bio import Entrez
#https://marcobonzanini.com/2015/01/12/searching-pubmed-with-python/

def search(query):
    Entrez.email = 'chanana@post.bgu.ac.il'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax='100000',
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle)
    return results

def fetch_details(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'chanana@post.bgu.ac.il'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results

"---------------- Main --------------------------------------"

term_occurrence_dict = {'Term':[], 'PubMed_Occurrence':[], 'PubMed_IDs':[]}

for term in cell_types:
    print('@', term)
    query = term.replace('_', ' ')
    query = query.replace('-', ' ')

    results = search(query)
    term_occurrence_dict['Term'].append(term)
    term_occurrence_dict['PubMed_Occurrence'].append(int(results['Count']))
    term_occurrence_dict['PubMed_IDs'].append('_'.join(results['IdList']))

PubMed_Search = pd.DataFrame.from_dict(term_occurrence_dict, orient='columns')
path = os.path.join('../..', 'Output', 'Text_Analysis_Output', 'Cells_Development_PubMed_Occurrence.tsv')
PubMed_Search.to_csv(path, index=False, sep='\t')

term_occurrence_dict = {'Term': [], 'PubMed_Occurrence': [], 'PubMed_IDs': []}

for term in ps_name_list:
    print('@', term)
    results = search(term)
    term_occurrence_dict['Term'].append(term)
    term_occurrence_dict['PubMed_Occurrence'].append(int(results['Count']))
    term_occurrence_dict['PubMed_IDs'].append('_'.join(results['IdList']))

PubMed_Search = pd.DataFrame.from_dict(term_occurrence_dict, orient='columns')
path = os.path.join('../..', 'Output', 'Text_Analysis_Output', 'Phenotypic_Series_PubMed_Occurrence.tsv')
PubMed_Search.to_csv(path, index=False, sep='\t')
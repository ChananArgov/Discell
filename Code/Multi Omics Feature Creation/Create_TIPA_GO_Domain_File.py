import os
import pandas as pd

path= os.path.join('../..', 'Output', 'TIPA_Cell_Type', 'General_Pref_TIPA_Developmental_scRNA_Cell_Types.csv')
GO_TIPA_Cell= pd.read_csv(path)
print(GO_TIPA_Cell)

path= os.path.join('../..', '..', 'data', 'GO_data', 'GO2GeneDomain.txt')
GO2GeneDomain = pd.read_csv(path, sep='\t')
print(GO2GeneDomain)
Relevant_GO2GeneDomain = GO2GeneDomain[['GO term accession', 'GO domain']][~GO2GeneDomain['GO term accession'].isnull()]
Relevant_GO2GeneDomain.drop_duplicates(inplace=True)
print(Relevant_GO2GeneDomain)

TIPA_Cell_Go_Domain = pd.merge(GO_TIPA_Cell, Relevant_GO2GeneDomain, left_on='GO_ID', right_on='GO term accession', how='left')
print(TIPA_Cell_Go_Domain)

path= os.path.join('../..', 'Output', 'TIPA_Cell_Type', 'General_Pref_TIPA_Developmental_scRNA_Cell_Types_GO_Domain.csv')
TIPA_Cell_Go_Domain.to_csv(path, index=False)

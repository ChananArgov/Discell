import os
import pandas as pd


path = os.path.join('../..', 'Output', 'Text_Analysis_Output', 'Cell_Type_Phenotypic_Series_PubMed_Overlaps.csv')
PubMed_Overlaps = pd.read_csv(path)
print(PubMed_Overlaps)

PubMed_Overlaps_Wide = PubMed_Overlaps.pivot(index='Disease_Phenotype', columns='Cell_Type', values='PubMed_Overlaps')
print(PubMed_Overlaps_Wide)
path = os.path.join('../..', 'Output', 'Text_Analysis_Output', 'Cell_Type_Phenotypic_Series_PubMed_Overlaps_Wide.csv')
PubMed_Overlaps_Wide.to_csv(path)
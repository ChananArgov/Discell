import os
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt

path = os.path.join('../..', 'Output', 'Text_Analysis_Output', 'Cell_Type_Phenotypic_Series_PubMed_Overlaps_Wide.csv')
PubMed_Overlaps_Wide = pd.read_csv(path, index_col=0)
print(PubMed_Overlaps_Wide)
PubMed_Overlaps_Wide = PubMed_Overlaps_Wide.replace(0, np.nan)
PubMed_Overlaps_Wide.dropna(axis='columns', how='all', inplace=True)
PubMed_Overlaps_Wide.dropna(axis='index', how='all', inplace=True)
PubMed_Overlaps_Wide.fillna(0, inplace=True)
print(PubMed_Overlaps_Wide)
# path = os.path.join('../..', 'Output', 'Text_Analysis_Output', 'Cell_Type_Phenotypic_Series_Cline.csv')
# PubMed_Overlaps_Wide.to_csv(path)

table = sm.stats.Table(PubMed_Overlaps_Wide)
print(table)

PubMed_Overlaps_Standardized_Resids = pd.DataFrame(data=table.standardized_resids, index=PubMed_Overlaps_Wide.index.tolist(), columns=list(PubMed_Overlaps_Wide))
print(PubMed_Overlaps_Standardized_Resids)

# path = os.path.join('../..', 'Output', 'Text_Analysis_Output', 'Cell_Type_Phenotypic_Series_Standardized_Resids.csv')
# PubMed_Overlaps_Standardized_Resids.to_csv(path)


PubMed_Overlaps_Standardized_Resids.plot.density()

plt.xlim(-5, 10)
plt.legend().set_visible(False)
# plt.show()
plt.xlabel('Standardized_Residuals')
path = os.path.join('../..', 'Output', 'Text_Analysis_Output', 'Occurrence_Standardized_Residuals_Plot.pdf')
# plt.savefig(path)
print(PubMed_Overlaps_Standardized_Resids.describe())
Percentage = PubMed_Overlaps_Standardized_Resids.quantile([.9, .95, 0.99])
print(Percentage)

PubMed_Overlaps_Standardized_Resids.reset_index(inplace=True)
PubMed_Overlaps_Standardized_Resids.rename(columns={'index':'Disease_Phenotype'}, inplace=True)
print(PubMed_Overlaps_Standardized_Resids)

Standardized_Resids_Long = PubMed_Overlaps_Standardized_Resids.melt(id_vars=['Disease_Phenotype'], value_vars=list(PubMed_Overlaps_Standardized_Resids), var_name='Cell_Type', value_name='Standardized_Residuals')
print(Standardized_Resids_Long)
path = os.path.join('../..', 'Output', 'Text_Analysis_Output', 'Cell_Type_Phenotypic_Series_Standardized_Residuals_Long.csv')
# Standardized_Resids_Long.to_csv(path, index=False)

print(Standardized_Resids_Long.describe())
Percentage = Standardized_Resids_Long.quantile([.9, .95, 0.96, 0.97, 0.975, 0.98, 0.99])
print(Percentage)

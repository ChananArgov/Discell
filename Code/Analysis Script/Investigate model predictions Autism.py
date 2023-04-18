import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import matplotlib as mpl

mpl.rcParams['axes.spines.left'] = True
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.bottom'] = True

path = os.path.join('..', '..', '..', 'Data', 'GWAS Analysis', 'SFARI-Gene_genes_11-07-2022release_11-28-2022export.csv')
SFARI = pd.read_csv(path)
print(SFARI)


path = os.path.join('..', '..', 'Results', 'SFARI_Autism_Lower_90_Prediction.csv')
Prediction = pd.read_csv(path)
print(Prediction)

SFARI_Predictions = pd.merge(SFARI, Prediction, left_on='ensembl-id', right_on='Gene_ID', how='inner')
print(SFARI_Predictions)

path = os.path.join('..', '..', 'Results', 'SFARI_Autism_Lower_90_Prediction_MetaData.csv')
SFARI_Predictions.to_csv(path)
correlation = SFARI_Predictions['number-of-reports'].corr(SFARI_Predictions['Discell Score'])

print(correlation)

f, ax = plt.subplots(figsize=(7, 8))

SFARI_Predictions['More than 5 publications'] = False
SFARI_Predictions.loc[SFARI_Predictions['number-of-reports']>5, 'More than 5 publications'] = True

sns.boxplot(x = 'More than 5 publications', y = 'Discell Score', data=SFARI_Predictions)

path = os.path.join('..', '..', 'Results', 'SFARI_Autism_Lower_10_Prediction.pdf')
plt.savefig(path)

low_publication_scores = SFARI_Predictions['Discell Score'][SFARI_Predictions['More than 5 publications'] == False].tolist()
high_publication_scores = SFARI_Predictions['Discell Score'][SFARI_Predictions['More than 5 publications'] == True].tolist()

U1, p = mannwhitneyu(high_publication_scores, low_publication_scores, alternative='greater')#greater

print('Mann Whitney Results: ', U1, p )

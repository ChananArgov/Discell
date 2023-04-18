import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

"------------------------- Load Data --------------------------------"

path = os.path.join('../..', 'Output', 'UMAP', 'Metadata_Sample.csv')
Metadata_Sample = pd.read_csv(path)
print(Metadata_Sample)

organs = list({c.split('-')[0] for c in Metadata_Sample['Organ_cell_lineage'].unique()})
organs.sort()
print('organs :', len(organs), organs)

# path = os.path.join('../..', 'Output', 'Transfer_Learning_Diseases', 'Prediction', 'Hyperinsulinemic hypoglycemia_PS256450_7_Cell_Specific_Probabilities.csv')
path = os.path.join('../..', 'Output', 'Transfer_Learning_Diseases', 'Prediction_Slim', 'Bleeding disorder_PS231200_32_Cell_Specific_Probabilities_Slim.csv')
# path = os.path.join('../..', 'Output', 'Transfer_Learning_Diseases', 'Prediction_Slim', 'Cardiomyopathy_PS604169_35_Cell_Specific_Probabilities_Slim.csv')
# path = os.path.join('../..', 'Output', 'Transfer_Learning_Diseases', 'Prediction_Slim', 'Myopathy_PS161800_25_Cell_Specific_Probabilities_Slim.csv')

Disease_Prediction = pd.read_csv(path)

print(Disease_Prediction)
Disease_Prediction_Mean = Disease_Prediction.groupby('Cell_Type')['scRNA_ML_Probability'].mean()
Disease_Prediction_Mean = pd.DataFrame(Disease_Prediction_Mean)
Disease_Prediction_Mean.reset_index(inplace=True)

print(Disease_Prediction_Mean)

prob_max = Disease_Prediction_Mean['scRNA_ML_Probability'].max()
prob_med = Disease_Prediction_Mean['scRNA_ML_Probability'].median()

Disease_Prediction_Mean['scRNA_ML_Probability'] = Disease_Prediction_Mean['scRNA_ML_Probability'] / prob_max
# Disease_Prediction_Mean['scRNA_ML_Probability'] = Disease_Prediction_Mean['scRNA_ML_Probability'] - prob_med

print(Disease_Prediction_Mean)

"----------------------------- Data Integration -----------------------"

Metadata_Sample_Prediction = pd.merge(Metadata_Sample, Disease_Prediction_Mean, left_on='Organ_cell_lineage', right_on='Cell_Type')
print(Metadata_Sample_Prediction)
Metadata_Sample_Prediction['Degree'] = 'Median'
Metadata_Sample_Prediction.loc[Metadata_Sample_Prediction['scRNA_ML_Probability'] < 0.2, 'Degree'] = 'Low'
Metadata_Sample_Prediction.loc[Metadata_Sample_Prediction['scRNA_ML_Probability'] > 0.6, 'Degree'] = 'High'
Metadata_Sample_Prediction.loc[Metadata_Sample_Prediction['scRNA_ML_Probability'] > 0.9, 'Degree'] = 'Very_high'


"------------------------- Load Data --------------------------------"


umap_types = ['Global_umap_', 'subcluster_umap_', 'Main_cluster_umap_']
fig, axs = plt.subplots(4, 4, figsize=(10, 8))
n = 0
for i in range(0, 4):
    for j in range(0, 4):
        organ = organs[n]
        ax = axs[i, j]
        sns.scatterplot(x='Main_cluster_umap_1', y='Main_cluster_umap_2',data=Metadata_Sample_Prediction[Metadata_Sample_Prediction['Organ_cell_lineage'].str.contains(organ)], hue="Degree", ax=ax, legend = False, s=1, palette=dict(Low="lightgrey", Median="lightgrey", High="#ffbfb8", Very_high="red"))#, c =Metadata_Sample_Prediction[Metadata_Sample_Prediction['Organ_cell_lineage'].str.contains(organ)]['scRNA_ML_Probability'])#hue="Main_cluster_name", palette='RdBu_r', High="#ffdad4"
        n += 1
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.title.set_text(organ)
        ax.set_frame_on(False)
        # ax.axhline(y=1, color='grey', linestyle='-')
        cell_types = Metadata_Sample_Prediction['Main_cluster_name'][Metadata_Sample_Prediction['Organ_cell_lineage'].str.contains(organ)].unique()
        """
        for ct in cell_types:
            ct_probability = Metadata_Sample_Prediction['scRNA_ML_Probability'][Metadata_Sample_Prediction['Main_cluster_name'] == ct].values[0]
            if ct_probability > 0.6:
                ct_xmax = Metadata_Sample_Prediction['Main_cluster_umap_1'][Metadata_Sample_Prediction['Main_cluster_name'] == ct].max()
                ct_xmin = Metadata_Sample_Prediction['Main_cluster_umap_1'][Metadata_Sample_Prediction['Main_cluster_name'] == ct].min()
                ct_ymax = Metadata_Sample_Prediction['Main_cluster_umap_2'][Metadata_Sample_Prediction['Main_cluster_name'] == ct].max()
                ct_ymin = Metadata_Sample_Prediction['Main_cluster_umap_2'][Metadata_Sample_Prediction['Main_cluster_name'] == ct].min()
                ct_ymean = Metadata_Sample_Prediction['Main_cluster_umap_2'][Metadata_Sample_Prediction['Main_cluster_name'] == ct].median()
                ct_xmean = Metadata_Sample_Prediction['Main_cluster_umap_1'][Metadata_Sample_Prediction['Main_cluster_name'] == ct].median()

                ct_x_center = (ct_xmax - ct_xmin) / 2 + ct_xmin
                ax.text(ct_xmean, ct_ymean, ct, fontsize=8)#, color='r'
        """
        if n == 15:
            break

ax = axs[3, 3]
ax.set_frame_on(False)
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

plt.tight_layout()
path = os.path.join('../..', 'Output', 'UMAP', 'Bleeding disorder_PS231200_32_UMAP_Slim.png')
# plt.savefig(path)
plt.show()